#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string reference;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> annotations_to_keep;
    std::vector<std::string> format_annotations_to_remove;
    std::string output;
    std::string manifest;
    std::vector<int> gq_bands{20, 100};
    // GATK declares this as a double (despite the short-name documentation
    // often showing an integer).  Keep fractional cutoffs observable rather
    // than truncating them at the command-line boundary.
    double rgq_threshold = 0.0;
    double tree_score_threshold = 0.0;
    bool floor_blocks = false;
    bool drop_low_quals = false;
    bool keep_all_alts = false;
    bool do_qual_approx = false;
    bool allow_missing_hom_ref_data = false;
    bool keep_filters = false;
    bool add_filters_to_genotype = false;
    bool create_index = true;
    bool explicit_gq_bands = false;
};

struct GenotypeKernelTelemetry {
    std::uint64_t pl_remap_calls = 0;
    double pl_remap_prepare_seconds = 0.0;
    double pl_remap_execute_seconds = 0.0;
    std::string pl_remap_execution_space;
    std::uint64_t allele_field_remap_calls = 0;
    double allele_field_remap_prepare_seconds = 0.0;
    double allele_field_remap_execute_seconds = 0.0;
    std::string allele_field_remap_execution_space;
    std::uint64_t raw_genotype_count_annotations = 0;
    std::uint64_t qual_approx_depth_annotations = 0;
    std::uint64_t reverse_allele_trimmed_variants = 0;
    std::uint64_t deletion_gap_ref_blocks = 0;
    std::uint64_t deletion_gap_ref_block_fallbacks = 0;
    std::uint64_t non_ref_ad_zeroed = 0;
    std::uint64_t non_ref_ad_kernel_calls = 0;
    double non_ref_ad_kernel_prepare_seconds = 0.0;
    double non_ref_ad_kernel_execute_seconds = 0.0;
    std::string non_ref_ad_kernel_execution_space;
    std::uint64_t info_annotations_removed = 0;
    std::uint64_t format_annotations_removed = 0;
    std::uint64_t converted_high_confidence_variants = 0;
    std::uint64_t overlapping_ref_block_trimmed = 0;
    std::uint64_t overlapping_ref_block_split = 0;
    std::uint64_t overlapping_ref_block_dropped = 0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* long_name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, long_name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == long_name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto character : text) {
        if (character == '"' || character == '\\') escaped << '\\';
        if (character == '\n') escaped << "\\n";
        else if (character == '\r') escaped << "\\r";
        else if (character == '\t') escaped << "\\t";
        else escaped << character;
    }
    return escaped.str();
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

#if FASTGATK_HAS_HTSLIB
bool has_header_field(const bcf_hdr_t* header, int header_type, const char* id) {
    return bcf_hdr_get_hrec(header, header_type, "ID", id, nullptr) != nullptr;
}

bool info_annotation_is_allele_specific(const bcf_hdr_t* header,
                                        const std::string& id) {
    const auto* hrec = bcf_hdr_get_hrec(header, BCF_HL_INFO, "ID", id.c_str(), nullptr);
    if (hrec == nullptr) return false;
    for (int key = 0; key < hrec->nkeys; ++key) {
        if (std::strcmp(hrec->keys[key], "Number") == 0 && hrec->vals[key] != nullptr)
            return std::strcmp(hrec->vals[key], "A") == 0;
    }
    return false;
}

std::string site_filter_string(const bcf_hdr_t* header, const bcf1_t* record) {
    if (record->d.n_flt <= 0 || record->d.flt == nullptr) return "PASS";
    std::string result;
    for (int index = 0; index < record->d.n_flt; ++index) {
        const char* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name == nullptr || *name == '\0' || std::strcmp(name, ".") == 0) continue;
        if (!result.empty()) result.push_back(';');
        result += name;
    }
    return result.empty() ? "PASS" : result;
}

const std::vector<std::string>& fixed_info_annotations_to_remove() {
    // GATK's ReblockGVCF removes annotations that are stale or not
    // informative for a reblocked GVCF.  These are the stable
    // IDs from ReblockGVCF.infoFieldAnnotationKeyNamesToRemove; unknown
    // caller-specific INFO fields remain available through
    // --annotations-to-keep.
    static const std::vector<std::string> ids{
        "GVCFBlock", "HaplotypeScore", "InbreedingCoeff", "MLEAC", "MLEAF",
        "ExcessHet", "AS_InbreedingCoeff", "Downsampled", "RAW_RMS_MAPPING_QUALITY"};
    return ids;
}

std::vector<std::string> stale_info_annotations_to_remove(const bcf_hdr_t* header) {
    std::vector<std::string> ids(fixed_info_annotations_to_remove().begin(),
                                 fixed_info_annotations_to_remove().end());
    // GATK emits one INFO header per GQ band (for example GVCFBlock10-20),
    // and the exact IDs depend on the configured bands.  Discover those
    // structured header IDs instead of assuming the default pair.
    for (int index = 0; index < header->nhrec; ++index) {
        const auto* hrec = header->hrec[index];
        if (hrec == nullptr || hrec->type != BCF_HL_INFO) continue;
        for (int key = 0; key < hrec->nkeys; ++key) {
            if (std::strcmp(hrec->keys[key], "ID") != 0 || hrec->vals[key] == nullptr) continue;
            const std::string id(hrec->vals[key]);
            if (id.rfind("GVCFBlock", 0) == 0 &&
                std::find(ids.begin(), ids.end(), id) == ids.end())
                ids.push_back(id);
        }
    }
    return ids;
}

bool annotation_is_kept(const Options& options, const std::string& id) {
    return std::find(options.annotations_to_keep.begin(),
                     options.annotations_to_keep.end(), id) !=
           options.annotations_to_keep.end();
}

bool record_has_info(const bcf_hdr_t* header, bcf1_t* record, const std::string& id) {
    if (!has_header_field(header, BCF_HL_INFO, id.c_str())) return false;
    bcf_unpack(record, BCF_UN_INFO);
    const int key = bcf_hdr_id2int(header, BCF_DT_ID, id.c_str());
    for (int index = 0; index < record->n_info; ++index)
        if (record->d.info[index].key == key && record->d.info[index].vptr != nullptr)
            return true;
    return false;
}

int info_int(const bcf_hdr_t* header, bcf1_t* record, const char* tag, int fallback);
double info_float(const bcf_hdr_t* header, bcf1_t* record, const char* tag, double fallback);

void remove_reblock_info_annotations(const bcf_hdr_t* header, bcf1_t* record,
                                     const Options& options,
                                     GenotypeKernelTelemetry& telemetry) {
    for (const auto& id : stale_info_annotations_to_remove(header)) {
        if (annotation_is_kept(options, id) || !record_has_info(header, record, id)) continue;
        if (bcf_update_info(header, record, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove INFO/" + id);
        ++telemetry.info_annotations_removed;
    }
}

// ReblockGVCF.updateMQAnnotations() materializes the modern RAW_MQandDP
// tuple for every high-quality variant.  Older callers may only provide the
// deprecated RAW_MQ value; in that case GATK preserves RAW_MQ, adds MQ_DP,
// and uses RAW_MQ (rather than MQ^2*DP) as the first tuple element.
void update_mq_annotations(const bcf_hdr_t* header, bcf1_t* record) {
    if (!has_header_field(header, BCF_HL_INFO, "RAW_MQandDP") ||
        record_has_info(header, record, "RAW_MQandDP"))
        return;
    const int depth = std::max(0, info_int(header, record, "DP", 0));
    int32_t raw_mq = 0;
    if (record_has_info(header, record, "RAW_MQ")) {
        raw_mq = static_cast<int32_t>(std::lround(info_float(header, record, "RAW_MQ", 0.0)));
        if (has_header_field(header, BCF_HL_INFO, "MQ_DP")) {
            const int32_t mq_dp = depth;
            if (bcf_update_info_int32(header, record, "MQ_DP", &mq_dp, 1) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write INFO/MQ_DP");
        }
    } else {
        const double mq = info_float(header, record, "MQ", 60.0);
        raw_mq = static_cast<int32_t>(std::lround(mq * mq * static_cast<double>(depth)));
    }
    const int32_t tuple[2] = {raw_mq, static_cast<int32_t>(depth)};
    if (bcf_update_info_int32(header, record, "RAW_MQandDP", tuple, 2) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write INFO/RAW_MQandDP");
}

void remove_reblock_format_annotations(const bcf_hdr_t* header, bcf1_t* record,
                                       const Options& options,
                                       GenotypeKernelTelemetry& telemetry) {
    for (const auto& id : options.format_annotations_to_remove) {
        if (!has_header_field(header, BCF_HL_FMT, id.c_str()))
            throw std::runtime_error("BAD_INPUT: FORMAT annotation not found: " + id);
        if (bcf_update_format(header, record, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove FORMAT/" + id);
        ++telemetry.format_annotations_removed;
    }
}

void remove_floor_block_arrays(const bcf_hdr_t* header, bcf1_t* record) {
    // ReblockingGVCFWriter documents --floor-blocks as a storage mode: the
    // block GQ is replaced by its band lower bound and PL/MIN_DP are dropped
    // from reference blocks.  Variant records retain PL and other likelihood
    // fields, so apply this only after the record has been classified.
    for (const char* tag : {"PL", "MIN_DP"}) {
        if (has_header_field(header, BCF_HL_FMT, tag) &&
            bcf_update_format(header, record, tag, nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot remove floor-block FORMAT/") + tag);
    }
}

void finalize_reblock_annotation_header(bcf_hdr_t* header, const Options& options) {
    for (const auto& id : stale_info_annotations_to_remove(header))
        if (!annotation_is_kept(options, id)) bcf_hdr_remove(header, BCF_HL_INFO, id.c_str());
    for (const auto& id : options.format_annotations_to_remove)
        bcf_hdr_remove(header, BCF_HL_FMT, id.c_str());
    if (bcf_hdr_sync(header) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync cleaned ReblockGVCF header");
}
#endif

int positive_int(const std::string& value, const char* option) {
    try {
        std::size_t consumed = 0;
        const auto parsed = std::stoi(value, &consumed);
        if (consumed != value.size() || parsed < 1 || parsed > 100) throw std::invalid_argument("range");
        return parsed;
    } catch (...) {
        throw std::invalid_argument(std::string("invalid GQ band for ") + option + ": " + value);
    }
}

double nonnegative_double(const std::string& value, const char* option) {
    try {
        std::size_t consumed = 0;
        const auto parsed = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(parsed) || parsed < 0.0 || parsed > 1000.0)
            throw std::invalid_argument("range");
        return parsed;
    } catch (...) {
        throw std::invalid_argument(std::string("invalid non-negative threshold for ") + option + ": " + value);
    }
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-reblock-gvcf (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE              input GVCF (repeatable shard)\n"
                         "  -R, --reference FILE            reference FASTA (accepted)\n"
                         "  -L, --intervals REGION          interval subset (repeatable)\n"
                         "      --interval-set-rule RULE  UNION (default) or INTERSECTION\n"
                         "  -GQB, --gvcf-gq-bands INTEGER   exclusive GQ band upper bound (repeatable)\n"
                         "  --floor-blocks                  output band lower bounds\n"
                         "  --drop-low-quals               convert low-quality sites to hom-ref blocks\n"
                         "  --annotations-to-keep ID       preserve a caller INFO annotation (repeatable)\n"
                         "  --format-annotations-to-remove ID\n"
                         "                                 remove FORMAT annotation (repeatable)\n"
                         "  -O, --output FILE               reblocked GVCF\n"
                         "      --output-manifest FILE     OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-V" || is_option(argument, "--variant"))
            options.inputs.push_back(require_value(index, argc, argv, argument, "--variant", "-V"));
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(index, argc, argv, argument,
                argument.rfind("--intervals", 0) == 0 ? "--intervals" :
                (argument.rfind("--interval", 0) == 0 ? "--interval" :
                (argument.rfind("--region", 0) == 0 ? "--region" : "--intervals")), "-L"));
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "-GQB" || is_option(argument, "--gvcf-gq-bands")) {
            if (!options.explicit_gq_bands) {
                options.gq_bands.clear();
                options.explicit_gq_bands = true;
            }
            options.gq_bands.push_back(positive_int(
                require_value(index, argc, argv, argument, "--gvcf-gq-bands", "-GQB"), "-GQB"));
        }
        else if (is_option(argument, "--rgq-threshold-to-no-call") || is_option(argument, "--rgq-threshold")) {
            // Both spellings are GATK aliases.  The longer name must be
            // passed to require_value explicitly; testing only the common
            // prefix (--rgq-threshold) otherwise rejects the separated form
            // ``--rgq-threshold-to-no-call 10`` as a missing value.
            const char* rgq_name = argument.rfind("--rgq-threshold-to-no-call", 0) == 0
                ? "--rgq-threshold-to-no-call" : "--rgq-threshold";
            options.rgq_threshold = nonnegative_double(
                require_value(index, argc, argv, argument, rgq_name), "--rgq-threshold");
        }
        else if (is_option(argument, "--tree-score-threshold-to-no-call"))
            options.tree_score_threshold = std::stod(require_value(index, argc, argv, argument,
                "--tree-score-threshold-to-no-call"));
        else if (argument == "--floor-blocks") options.floor_blocks = true;
        else if (argument == "--drop-low-quals") options.drop_low_quals = true;
        else if (argument == "--keep-all-alts") options.keep_all_alts = true;
        else if (argument == "--do-qual-score-approximation" ||
                 argument == "--do-qual-score-approx" || argument == "--do-qual-approx")
            options.do_qual_approx = true;
        else if (argument == "--allow-missing-hom-ref-data") options.allow_missing_hom_ref_data = true;
        else if (is_option(argument, "--annotations-to-keep"))
            options.annotations_to_keep.push_back(require_value(
                index, argc, argv, argument, "--annotations-to-keep"));
        else if (is_option(argument, "--format-annotations-to-remove") ||
                 is_option(argument, "--annotations-to-remove")) {
            const char* annotation_option =
                argument.rfind("--format-annotations-to-remove", 0) == 0
                    ? "--format-annotations-to-remove" : "--annotations-to-remove";
            options.format_annotations_to_remove.push_back(require_value(
                index, argc, argv, argument, annotation_option));
        }
        else if (argument == "--keep-site-filters" || argument == "--keep-filters") options.keep_filters = true;
        else if (argument == "--add-site-filters-to-genotype") options.add_filters_to_genotype = true;
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility flags; HTSlib still validates records against the header.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("--variant/-V is required");
    if (options.output.empty()) throw std::invalid_argument("--output/-O is required");
    std::sort(options.gq_bands.begin(), options.gq_bands.end());
    options.gq_bands.erase(std::unique(options.gq_bands.begin(), options.gq_bands.end()), options.gq_bands.end());
    if (options.gq_bands.empty() || options.gq_bands.back() != 100)
        options.gq_bands.push_back(100);
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct Region {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
};

void append_region_selector(const std::string& selector, const bcf_hdr_t* header,
                            std::vector<Region>& regions,
                            fastgatk::io::IntervalFileStats& stats) {
    std::vector<fastgatk::io::IndexedInterval> parsed;
    fastgatk::io::append_interval_selector(selector, header, parsed, stats);
    for (const auto& interval : parsed)
        regions.push_back(Region{interval.rid, interval.begin, interval.end});
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_reblock_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.end <= region.begin) continue;
        if (!merged.empty() && merged.back().rid == region.rid &&
            region.begin <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, region.end);
        } else {
            merged.push_back(region);
        }
    }
    regions.swap(merged);
}

// Apply GATK's set operation to repeated -L selectors before record-span
// filtering.  INTERSECTION is a true half-open interval intersection rather
// than a post-hoc filter on the selector union.
void append_reblock_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_reblock_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    intersection.reserve(std::min(regions.size(), incoming.size()));
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].rid < incoming[right].rid) {
            ++left;
            continue;
        }
        if (incoming[right].rid < regions[left].rid) {
            ++right;
            continue;
        }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end)
            intersection.push_back(Region{regions[left].rid, begin, end});
        if (regions[left].end < incoming[right].end) ++left;
        else ++right;
    }
    regions.swap(intersection);
}

bool overlaps_regions(const bcf1_t* record, const std::vector<Region>& regions,
                      int end_inclusive) {
    if (regions.empty()) return true;
    const int end_exclusive = std::max(static_cast<int>(record->pos + 1), end_inclusive);
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        return record->rid == region.rid && record->pos < region.end &&
               end_exclusive > region.begin;
    });
}

struct RegionRecord {
    bcf1_t* value = nullptr;
    int rid = -1;
    int pos = -1;
    int end = -1;
    int gq = 0;
    int min_dp = -1;
    bool reference_block = false;
    std::vector<int32_t> gqs;
    std::vector<int32_t> min_dps;
    // Keep the per-input DP samples so a merged GVCF block can reproduce
    // GVCFBlock's rounded median DP rather than simply retaining the first
    // record's value.
    std::vector<std::vector<int32_t>> dp_history;
};

void destroy_records(std::vector<RegionRecord>& records) {
    for (auto& record : records) bcf_destroy(record.value);
    records.clear();
}

int info_int(const bcf_hdr_t* header, bcf1_t* record, const char* tag, int fallback) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, tag, &values, &count);
    const auto result = length > 0 && count > 0 ? values[0] : fallback;
    free(values);
    return result;
}

int format_int(const bcf_hdr_t* header, bcf1_t* record, const char* tag, int fallback,
               int sample_index = 0) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, tag, &values, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    const int offset = width > 0 && sample_index >= 0 && sample_index < sample_count
        ? sample_index * width : 0;
    const auto result = width > 0 && values[offset] != bcf_int32_missing &&
        values[offset] != bcf_int32_vector_end ? values[offset] : fallback;
    free(values);
    return result;
}

std::vector<int32_t> format_values(const bcf_hdr_t* header, bcf1_t* record,
                                   const char* tag, int sample_count, int fallback) {
    sample_count = std::max(1, sample_count);
    std::vector<int32_t> result(static_cast<std::size_t>(sample_count), fallback);
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, tag, &values, &count);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    if (width > 0) {
        for (int sample = 0; sample < sample_count; ++sample) {
            const auto value = values[sample * width];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end)
                result[static_cast<std::size_t>(sample)] = value;
        }
    }
    free(values);
    return result;
}

void set_format_values(const bcf_hdr_t* header, bcf1_t* record, const char* tag,
                       const std::vector<int32_t>& values) {
    if (values.empty() || !has_header_field(header, BCF_HL_FMT, tag)) return;
    if (bcf_update_format_int32(header, record, tag, values.data(),
                                static_cast<int>(values.size())) != 0)
        throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot update reblocked ") + tag);
}

void set_format_string_values(const bcf_hdr_t* header, bcf1_t* record, const char* tag,
                              const std::vector<std::string>& values) {
    if (values.empty() || !has_header_field(header, BCF_HL_FMT, tag)) return;
    std::vector<const char*> pointers;
    pointers.reserve(values.size());
    for (const auto& value : values) pointers.push_back(value.c_str());
    if (bcf_update_format_string(header, record, tag, pointers.data(),
                                 static_cast<int>(pointers.size())) != 0)
        throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot update reblocked ") + tag);
}

struct VariantDepthResult {
    int variant_depth = 0;
    std::vector<int32_t> allele_depths;
    bool has_allele_depths = false;
};

// Mirrors the raw depth portion of GATK QualByDepth/AS_QualByDepth for the
// native ReblockGVCF path.  Only called variant genotypes
// contribute; usable AD rows are preferred, and rows with <=1 alternate read
// are excluded from the AD-restricted depth just like the Java annotation.
VariantDepthResult variant_depth(const bcf_hdr_t* header, bcf1_t* record,
                                 int fallback_depth) {
    VariantDepthResult result;
    int32_t* gt = nullptr;
    int gt_count = 0;
    const auto gt_length = bcf_get_genotypes(header, record, &gt, &gt_count);
    const int sample_count = header->n[BCF_DT_SAMPLE];
    const int ploidy = gt_length > 0 && sample_count > 0 && gt_count % sample_count == 0
        ? gt_count / sample_count : 0;
    int32_t* ad = nullptr;
    int ad_count = 0;
    const auto ad_length = bcf_get_format_int32(header, record, "AD", &ad, &ad_count);
    const int ad_width = ad_length > 0 && sample_count > 0 && ad_count % sample_count == 0
        ? ad_count / sample_count : 0;
    int32_t* dp = nullptr;
    int dp_count = 0;
    const auto dp_length = bcf_get_format_int32(header, record, "DP", &dp, &dp_count);
    const int dp_width = dp_length > 0 && sample_count > 0 && dp_count % sample_count == 0
        ? dp_count / sample_count : 0;
    int depth = 0;
    int ad_restricted_depth = 0;
    if (ad_width > 0) {
        result.has_allele_depths = true;
        result.allele_depths.assign(static_cast<std::size_t>(ad_width), 0);
    }
    for (int sample = 0; sample < sample_count; ++sample) {
        bool variant = false;
        if (ploidy > 0) {
            for (int copy = 0; copy < ploidy; ++copy) {
                const auto encoded = gt[sample * ploidy + copy];
                if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
                if (bcf_gt_allele(encoded) > 0) {
                    variant = true;
                    break;
                }
            }
        }
        if (!variant) continue;
        if (ad_width > 0) {
            const auto* row = ad + sample * ad_width;
            int total = 0;
            for (int allele = 0; allele < ad_width; ++allele)
                if (row[allele] != bcf_int32_missing && row[allele] != bcf_int32_vector_end && row[allele] >= 0)
                    total += row[allele];
            const int ref_depth = row[0] >= 0 ? row[0] : 0;
            if (total != 0) {
                depth += total;
                if (total - ref_depth > 1) {
                    ad_restricted_depth += total;
                    for (int allele = 0; allele < ad_width; ++allele)
                        if (row[allele] != bcf_int32_missing && row[allele] != bcf_int32_vector_end && row[allele] >= 0)
                            result.allele_depths[static_cast<std::size_t>(allele)] += row[allele];
                }
                continue;
            }
        }
        if (dp_width > 0 && dp[sample * dp_width] >= 0)
            depth += dp[sample * dp_width];
    }
    free(gt);
    free(ad);
    free(dp);
    result.variant_depth = ad_restricted_depth > 0 ? ad_restricted_depth : depth;
    if (result.variant_depth == 0 && fallback_depth > 0) result.variant_depth = fallback_depth;
    if (result.variant_depth == 0) result.variant_depth = 1;
    return result;
}

int format_pl0(const bcf_hdr_t* header, bcf1_t* record, int fallback,
               int sample_index = 0) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "PL", &values, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    const int offset = width > 0 && sample_index >= 0 && sample_index < sample_count
        ? sample_index * width : 0;
    const auto result = width > 0 && values[offset] != bcf_int32_missing &&
        values[offset] != bcf_int32_vector_end ? values[offset] : fallback;
    free(values);
    return result;
}

double info_float(const bcf_hdr_t* header, bcf1_t* record, const char* tag, double fallback) {
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, tag, &values, &count);
    const auto result = length > 0 && count > 0 && !bcf_float_is_missing(values[0]) &&
        !bcf_float_is_vector_end(values[0]) ? static_cast<double>(values[0]) : fallback;
    free(values);
    return result;
}

int genotype_gq(const bcf_hdr_t* header, bcf1_t* record, int sample_index = 0) {
    const auto explicit_gq = format_int(header, record, "GQ", -1, sample_index);
    if (explicit_gq >= 0) return explicit_gq;
    int32_t* pls = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "PL", &pls, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    if (width < 2 || sample_index < 0 || sample_index >= sample_count) {
        free(pls);
        return -1;
    }
    std::vector<int32_t> values;
    values.reserve(static_cast<std::size_t>(width));
    const int offset = sample_index * width;
    for (int index = 0; index < width; ++index) {
        if (pls[offset + index] == bcf_int32_missing || pls[offset + index] == bcf_int32_vector_end) continue;
        values.push_back(pls[offset + index]);
    }
    free(pls);
    if (values.size() < 2) return -1;
    std::nth_element(values.begin(), values.begin(), values.end());
    const auto best = values.front();
    std::nth_element(values.begin() + 1, values.begin() + 1, values.end());
    return std::min(99, std::max(0, static_cast<int>(values[1] - best)));
}

int pl_gq(const bcf_hdr_t* header, bcf1_t* record, int sample_index = 0) {
    int32_t* pls = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "PL", &pls, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    if (width < 2 || sample_index < 0 || sample_index >= sample_count) {
        free(pls);
        return -1;
    }
    std::vector<int32_t> values;
    values.reserve(static_cast<std::size_t>(width));
    const int offset = sample_index * width;
    for (int index = 0; index < width; ++index) {
        if (pls[offset + index] == bcf_int32_missing || pls[offset + index] == bcf_int32_vector_end) continue;
        values.push_back(pls[offset + index]);
    }
    free(pls);
    if (values.size() < 2) return -1;
    std::nth_element(values.begin(), values.begin(), values.end());
    const auto best = values.front();
    std::nth_element(values.begin() + 1, values.begin() + 1, values.end());
    return std::min(99, std::max(0, static_cast<int>(values[1] - best)));
}

bool pl_best_genotype_is_ref(const bcf_hdr_t* header, bcf1_t* record,
                             int sample_index = 0) {
    int32_t* pls = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "PL", &pls, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    if (width < 1 || sample_index < 0 || sample_index >= sample_count) {
        free(pls);
        return false;
    }
    const int offset = sample_index * width;
    const auto best = pls[offset];
    bool valid = best != bcf_int32_missing && best != bcf_int32_vector_end;
    for (int index = 1; valid && index < width; ++index) {
        const auto value = pls[offset + index];
        if (value != bcf_int32_missing && value != bcf_int32_vector_end && value < best)
            valid = false;
    }
    free(pls);
    return valid;
}

bool has_non_ref(const bcf1_t* record) {
    for (int allele = 1; allele < record->n_allele; ++allele)
        if (std::strcmp(record->d.allele[allele], "<NON_REF>") == 0) return true;
    return false;
}

bool has_concrete_alt(const bcf1_t* record) {
    for (int allele = 1; allele < record->n_allele; ++allele) {
        const auto* value = record->d.allele[allele];
        if (value != nullptr && value[0] != '<') return true;
    }
    return false;
}

std::size_t genotype_width(int allele_count, int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t result = 1;
    for (int step = 1; step <= ploidy; ++step) {
        const auto numerator = static_cast<std::size_t>(allele_count + step - 1);
        if (result > std::numeric_limits<std::size_t>::max() / numerator) return 0;
        result *= numerator;
        result /= static_cast<std::size_t>(step);
    }
    return result;
}

void update_gt(const bcf_hdr_t* input_header, const bcf_hdr_t* output_header,
               bcf1_t* record, const std::vector<int>& kept) {
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "GT") < 0) return;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(input_header, record, &values, &count);
    if (length <= 0 || count < 1) { free(values); return; }
    const auto old_width = count;
    for (int index = 0; index < old_width; ++index) {
        if (bcf_gt_is_missing(values[index]) || values[index] == bcf_int32_vector_end) {
            values[index] = bcf_gt_missing;
            continue;
        }
        const auto old_allele = bcf_gt_allele(values[index]);
        const auto found = std::find(kept.begin(), kept.end(), old_allele);
        const bool phased = bcf_gt_is_phased(values[index]);
        values[index] = found == kept.end() ? bcf_gt_missing :
            (phased ? bcf_gt_phased(static_cast<int>(found - kept.begin())) :
                      bcf_gt_unphased(static_cast<int>(found - kept.begin())));
    }
    if (bcf_update_genotypes(output_header, record, values, count) != 0) {
        free(values);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update reblocked GT");
    }
    free(values);
}

void update_ad_pl(const bcf_hdr_t* input_header, const bcf_hdr_t* output_header,
                  bcf1_t* record, const std::vector<int>& kept, int old_alleles,
                  GenotypeKernelTelemetry& telemetry) {
    const int sample_count = std::max(1, input_header->n[BCF_DT_SAMPLE]);
    int ploidy = 2;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    if (bcf_get_genotypes(input_header, record, &genotypes, &genotype_count) > 0 && genotype_count > 0)
        ploidy = genotype_count / sample_count;
    if (ploidy <= 0) ploidy = 2;
    free(genotypes);
    const auto old_pl_width = genotype_width(old_alleles, ploidy);
    const auto new_pl_width = genotype_width(static_cast<int>(kept.size()), ploidy);
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AD") >= 0) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(input_header, record, "AD", &values, &count);
        const int source_width = length > 0 && count >= sample_count && count % sample_count == 0
            ? count / sample_count : 0;
        if (source_width >= old_alleles) {
            std::vector<int32_t> source(static_cast<std::size_t>(sample_count * old_alleles), bcf_int32_missing);
            for (int sample = 0; sample < sample_count; ++sample)
                for (int allele = 0; allele < old_alleles; ++allele)
                    source[static_cast<std::size_t>(sample * old_alleles + allele)] =
                        values[sample * source_width + allele];
            std::vector<int32_t> target_to_source;
            target_to_source.reserve(kept.size());
            for (const auto allele : kept) target_to_source.push_back(allele);
            const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
                source, static_cast<std::size_t>(sample_count), old_alleles,
                static_cast<int>(kept.size()), target_to_source);
            if (bcf_update_format_int32(output_header, record, "AD", remapped.values.data(),
                                        static_cast<int>(remapped.values.size())) != 0) {
                free(values);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update reblocked AD");
            }
            ++telemetry.allele_field_remap_calls;
            telemetry.allele_field_remap_prepare_seconds += remapped.prepare_seconds;
            telemetry.allele_field_remap_execute_seconds += remapped.seconds;
            telemetry.allele_field_remap_execution_space = remapped.execution_space;
        }
        free(values);
    }
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "PL") >= 0 && old_pl_width > 0 && new_pl_width > 0) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(input_header, record, "PL", &values, &count);
        const int source_width = length > 0 && count >= sample_count && count % sample_count == 0
            ? count / sample_count : 0;
        if (source_width >= static_cast<int>(old_pl_width)) {
            std::vector<std::int32_t> source(
                static_cast<std::size_t>(sample_count) * old_pl_width, bcf_int32_missing);
            for (int sample = 0; sample < sample_count; ++sample)
                for (std::size_t cell = 0; cell < old_pl_width; ++cell)
                    source[static_cast<std::size_t>(sample) * old_pl_width + cell] =
                        values[static_cast<std::size_t>(sample) * source_width + cell];
            std::vector<std::int32_t> kept_indices;
            kept_indices.reserve(kept.size());
            for (const auto allele : kept) kept_indices.push_back(static_cast<std::int32_t>(allele));
            const auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
                source, static_cast<std::size_t>(sample_count), old_alleles,
                static_cast<int>(kept.size()), ploidy, kept_indices);
            if (bcf_update_format_int32(output_header, record, "PL", remapped.pl.data(),
                                        static_cast<int>(remapped.pl.size())) != 0) {
                free(values);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update reblocked PL");
            }
            ++telemetry.pl_remap_calls;
            telemetry.pl_remap_prepare_seconds += remapped.prepare_seconds;
            telemetry.pl_remap_execute_seconds += remapped.seconds;
            telemetry.pl_remap_execution_space = remapped.execution_space;
        }
        free(values);
    }
}

void set_scalar_format(const bcf_hdr_t* header, bcf1_t* record, const char* tag, int value);
void update_gq_from_pl(const bcf_hdr_t* header, bcf1_t* record);

// GATK's removeNonRefADs prevents reads assigned to <NON_REF> from being
// propagated to every concrete ALT during later GVCF merges.  Keep the
// mutation on the Kokkos execution path so this annotation cleanup remains
// backend-portable, then apply the scalar DP correction at the HTSlib boundary.
bool remove_non_ref_ad_kokkos(const bcf_hdr_t* header, bcf1_t* record,
                              GenotypeKernelTelemetry& telemetry) {
    int non_ref_index = -1;
    for (int allele = 1; allele < record->n_allele; ++allele) {
        if (record->d.allele[allele] != nullptr &&
            std::strcmp(record->d.allele[allele], "<NON_REF>") == 0) {
            non_ref_index = allele;
            break;
        }
    }
    if (non_ref_index < 0 || !has_header_field(header, BCF_HL_FMT, "AD")) return false;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int width = length > 0 && count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    if (width < record->n_allele) {
        free(values);
        return false;
    }
    std::vector<int32_t> non_ref_depths(static_cast<std::size_t>(sample_count), 0);
    bool any_non_ref = false;
    for (int sample = 0; sample < sample_count; ++sample) {
        const auto value = values[sample * width + non_ref_index];
        if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0) {
            non_ref_depths[static_cast<std::size_t>(sample)] = value;
            any_non_ref = true;
        }
    }
    if (!any_non_ref) {
        free(values);
        return false;
    }
    using execution_space = Kokkos::DefaultExecutionSpace;
    Kokkos::View<std::int32_t*, execution_space> device_values(
        "reblock_non_ref_ad", static_cast<std::size_t>(count));
    auto host_values = Kokkos::create_mirror_view(device_values);
    fastgatk::core::HostBatch host_batch("reblock-non-ref-ad-v1");
    host_batch.records = static_cast<std::size_t>(count);
    host_batch.bytes = static_cast<std::size_t>(count) * sizeof(std::int32_t);
    fastgatk::core::KernelPlan<execution_space> kernel_plan("reblock-non-ref-ad");
    kernel_plan.begin_prepare(host_batch);
    for (int index = 0; index < count; ++index) host_values(index) = values[index];
    Kokkos::deep_copy(device_values, host_values);
    fastgatk::core::DeviceBatch<execution_space> device_batch(static_cast<std::size_t>(count));
    device_batch.bind("values", device_values);
    execution_space().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    Kokkos::parallel_for("reblock_zero_non_ref_ad",
                         Kokkos::RangePolicy<execution_space>(0, count),
                         KOKKOS_LAMBDA(const int index) {
                             if (index % width == non_ref_index) device_values(index) = 0;
                         });
    execution_space().fence();
    kernel_plan.end_execute();
    ++telemetry.non_ref_ad_kernel_calls;
    telemetry.non_ref_ad_kernel_prepare_seconds += kernel_plan.telemetry().prepare_seconds;
    telemetry.non_ref_ad_kernel_execute_seconds += kernel_plan.telemetry().execute_seconds;
    telemetry.non_ref_ad_kernel_execution_space = execution_space::name();
    Kokkos::deep_copy(host_values, device_values);
    std::vector<std::int32_t> updated(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) updated[static_cast<std::size_t>(index)] = host_values(index);
    free(values);
    if (bcf_update_format_int32(header, record, "AD", updated.data(), count) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot zero reblocked NON_REF AD");
    int32_t* dp_values = nullptr;
    int dp_count = 0;
    const auto dp_length = bcf_get_format_int32(header, record, "DP", &dp_values, &dp_count);
    const int dp_width = dp_length > 0 && dp_count >= sample_count && dp_count % sample_count == 0
        ? dp_count / sample_count : 0;
    if (dp_width > 0) {
        std::vector<int32_t> updated_dp(dp_values, dp_values + dp_count);
        for (int sample = 0; sample < sample_count; ++sample) {
            auto& value = updated_dp[static_cast<std::size_t>(sample * dp_width)];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end &&
                non_ref_depths[static_cast<std::size_t>(sample)] > 0)
                value = std::max(0, value - non_ref_depths[static_cast<std::size_t>(sample)]);
        }
        if (bcf_update_format_int32(header, record, "DP", updated_dp.data(), dp_count) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update reblocked DP");
    }
    free(dp_values);
    ++telemetry.non_ref_ad_zeroed;
    return true;
}

void set_scalar_format(const bcf_hdr_t* header, bcf1_t* record, const char* tag, int value) {
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    set_format_values(header, record, tag,
                     std::vector<int32_t>(static_cast<std::size_t>(sample_count), value));
}

// ReblockGVCF emits this GATK bookkeeping annotation for high-quality
// variants. ExcessHet consumes rounded genotype counts in the order
// hom-ref/het/hom-var.  Count every called sample independently rather than
// collapsing a multi-sample record to the first GT row.
bool update_raw_genotype_count(const bcf_hdr_t* header, bcf1_t* record) {
    if (!has_header_field(header, BCF_HL_INFO, "RAW_GT_COUNT")) return false;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &values, &count);
    if (length <= 0 || count <= 0) {
        free(values);
        return false;
    }
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    const int ploidy = count >= sample_count && count % sample_count == 0
        ? count / sample_count : 0;
    int32_t raw_count[3] = {0, 0, 0};
    if (ploidy > 0) {
        for (int sample = 0; sample < sample_count; ++sample) {
            bool called = false;
            bool all_ref = true;
            bool all_same_alt = true;
            int first_alt = -1;
            for (int copy = 0; copy < ploidy; ++copy) {
                const auto encoded = values[sample * ploidy + copy];
                if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                    called = false;
                    all_ref = false;
                    all_same_alt = false;
                    break;
                }
                const int allele = bcf_gt_allele(encoded);
                called = true;
                if (allele != 0) all_ref = false;
                if (allele == 0) all_same_alt = false;
                else if (first_alt < 0) first_alt = allele;
                else if (first_alt != allele) all_same_alt = false;
            }
            if (!called) continue;
            if (all_ref) ++raw_count[0];
            else if (all_same_alt) ++raw_count[2];
            else ++raw_count[1];
        }
    }
    free(values);
    if (bcf_update_info_int32(header, record, "RAW_GT_COUNT", raw_count, 3) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write RAW_GT_COUNT");
    return true;
}

void update_gq_from_pl(const bcf_hdr_t* header, bcf1_t* record) {
    if (bcf_hdr_id2int(header, BCF_DT_ID, "GQ") < 0) return;
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    std::vector<int32_t> gqs(static_cast<std::size_t>(sample_count), 0);
    bool any = false;
    for (int sample = 0; sample < sample_count; ++sample) {
        const auto gq = pl_gq(header, record, sample);
        if (gq >= 0) {
            gqs[static_cast<std::size_t>(sample)] = gq;
            any = true;
        }
    }
    if (any) set_format_values(header, record, "GQ", gqs);
}

void convert_to_ref_block(const bcf_hdr_t* header, bcf1_t* record, int end,
                          const std::vector<int32_t>& source_min_dps,
                          int fallback_min_dp) {
    // GATK's changeCallToHomRefVersusNonRef trims a deletion REF allele to
    // its leading base when converting a low-quality call to a reference
    // block.  A spanning-deletion or no-call GT triggers the same trim even
    // when the REF itself is not longer than one base; END retains the full
    // covered span below.
    bool trim_reference = record->d.allele[0] != nullptr &&
        std::strlen(record->d.allele[0]) > 1;
    const int sample_count = std::max(1, header->n[BCF_DT_SAMPLE]);
    int ploidy = 2;
    int32_t* old_genotypes = nullptr;
    int genotype_count = 0;
    if (bcf_get_genotypes(header, record, &old_genotypes, &genotype_count) > 0 && genotype_count > 0)
        ploidy = genotype_count / sample_count;
    if (ploidy <= 0) ploidy = 2;
    for (int index = 0; index < genotype_count; ++index) {
        if (old_genotypes[index] == bcf_int32_vector_end) continue;
        if (bcf_gt_is_missing(old_genotypes[index])) {
            trim_reference = true;
            continue;
        }
        const int allele = bcf_gt_allele(old_genotypes[index]);
        if (allele >= 0 && allele < record->n_allele &&
            std::strcmp(record->d.allele[allele], "*") == 0)
            trim_reference = true;
    }
    free(old_genotypes);
    std::string reference = record->d.allele[0] == nullptr ? "N" : record->d.allele[0];
    if (trim_reference && reference.size() > 1) reference.resize(1);
    const std::string alleles = reference + ",<NON_REF>";
    if (bcf_update_alleles_str(header, record, alleles.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create reference block alleles");
    // ReblockGVCF's changeCallToHomRefVersusNonRef() uses noAD()/noAttributes()
    // when a low-quality call becomes a reference block.  AD describes the
    // discarded concrete alleles and must not survive this conversion as a
    // stale annotation (the later NON_REF cleanup is only for high-quality
    // variants).
    if (has_header_field(header, BCF_HL_FMT, "AD") &&
        bcf_update_format(header, record, "AD", nullptr, 0, BCF_HT_STR) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove AD from reference block");
    std::vector<int32_t> gt(static_cast<std::size_t>(sample_count * ploidy), bcf_gt_unphased(0));
    if (bcf_update_genotypes(header, record, gt.data(), static_cast<int>(gt.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create reference block GT");
    const auto pl_width = genotype_width(2, ploidy);
    std::vector<int32_t> pl(static_cast<std::size_t>(sample_count) * pl_width, 0);
    if (bcf_hdr_id2int(header, BCF_DT_ID, "PL") >= 0 &&
        bcf_update_format_int32(header, record, "PL", pl.data(), static_cast<int>(pl.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create reference block PL");
    set_scalar_format(header, record, "GQ", 0);
    const int sample_count_for_min_dp = std::max(1, header->n[BCF_DT_SAMPLE]);
    std::vector<int32_t> min_dps(static_cast<std::size_t>(sample_count_for_min_dp), fallback_min_dp);
    for (int sample = 0; sample < sample_count_for_min_dp; ++sample) {
        if (sample < static_cast<int>(source_min_dps.size()) && source_min_dps[static_cast<std::size_t>(sample)] >= 0)
            min_dps[static_cast<std::size_t>(sample)] = source_min_dps[static_cast<std::size_t>(sample)];
    }
    if (fallback_min_dp >= 0 || !source_min_dps.empty())
        set_format_values(header, record, "MIN_DP", min_dps);
    const int32_t end_value = end;
    if (bcf_update_info_int32(header, record, "END", &end_value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update reference block END");
}

void clear_info_except_end(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_INFO);
    std::vector<std::string> ids;
    ids.reserve(static_cast<std::size_t>(record->n_info));
    for (int index = 0; index < record->n_info; ++index) {
        const char* id = bcf_hdr_int2id(header, BCF_DT_ID, record->d.info[index].key);
        if (id != nullptr && std::strcmp(id, "END") != 0) ids.emplace_back(id);
    }
    for (const auto& id : ids) {
        if (bcf_update_info(header, record, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear reference-block INFO/" + id);
    }
}

void clear_format_attributes_except_core(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FMT);
    std::vector<std::string> ids;
    ids.reserve(static_cast<std::size_t>(record->n_fmt));
    for (int index = 0; index < record->n_fmt; ++index) {
        const char* id = bcf_hdr_int2id(header, BCF_DT_ID, record->d.fmt[index].id);
        if (id == nullptr) continue;
        if (std::strcmp(id, "GT") != 0 && std::strcmp(id, "DP") != 0 &&
            std::strcmp(id, "GQ") != 0 && std::strcmp(id, "MIN_DP") != 0 &&
            std::strcmp(id, "PL") != 0)
            ids.emplace_back(id);
    }
    for (const auto& id : ids) {
        if (bcf_update_format(header, record, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear reference-block FORMAT/" + id);
    }
}

void order_reference_block_format(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FMT);
    const auto rank = [&](const bcf_fmt_t& field) {
        const char* id = bcf_hdr_int2id(header, BCF_DT_ID, field.id);
        if (id == nullptr) return 100;
        if (std::strcmp(id, "GT") == 0) return 0;
        if (std::strcmp(id, "DP") == 0) return 1;
        if (std::strcmp(id, "GQ") == 0) return 2;
        if (std::strcmp(id, "MIN_DP") == 0) return 3;
        if (std::strcmp(id, "PL") == 0) return 4;
        return 100;
    };
    std::stable_sort(record->d.fmt, record->d.fmt + record->n_fmt,
                     [&](const bcf_fmt_t& left, const bcf_fmt_t& right) {
                         return rank(left) < rank(right);
                     });
    record->d.indiv_dirty = 1;
}

// GATK's high-confidence hom-ref path first re-genotypes the concrete variant,
// then subsets it to REF plus the most likely ALT (which is represented by
// <NON_REF> in the resulting reference block).  This is distinct from the
// low-confidence conversion above, whose PL/GQ are intentionally all zero.
// Keep the allele-subset operation on the shared Kokkos Number=G kernel so the
// conversion has the same backend behavior as ordinary ALT compaction.
bool convert_high_confidence_homref_to_ref_block(
    const bcf_hdr_t* input_header, const bcf_hdr_t* output_header, bcf1_t* record,
    int end, const std::vector<int32_t>& source_min_dps, int fallback_min_dp,
    GenotypeKernelTelemetry& telemetry) {
    const int sample_count = std::max(1, input_header->n[BCF_DT_SAMPLE]);
    const int old_alleles = record->n_allele;
    if (old_alleles < 3) return false;

    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto gt_length = bcf_get_genotypes(input_header, record, &genotypes, &genotype_count);
    const int ploidy = gt_length > 0 && genotype_count > 0 && genotype_count % sample_count == 0
        ? genotype_count / sample_count : 2;
    if (ploidy <= 0) { free(genotypes); return false; }

    int32_t* pl_values = nullptr;
    int pl_count = 0;
    const auto pl_length = bcf_get_format_int32(input_header, record, "PL", &pl_values, &pl_count);
    const auto old_pl_width = genotype_width(old_alleles, ploidy);
    const auto target_pl_width = genotype_width(2, ploidy);
    const int source_width = pl_length > 0 && pl_count >= sample_count &&
        pl_count % sample_count == 0 ? pl_count / sample_count : 0;
    if (source_width < static_cast<int>(old_pl_width) || old_pl_width == 0 || target_pl_width == 0) {
        free(genotypes); free(pl_values); return false;
    }
    std::vector<int32_t> source(static_cast<std::size_t>(sample_count) * old_pl_width,
                                bcf_int32_missing);
    for (int sample = 0; sample < sample_count; ++sample)
        for (std::size_t cell = 0; cell < old_pl_width; ++cell)
            source[static_cast<std::size_t>(sample) * old_pl_width + cell] =
                pl_values[static_cast<std::size_t>(sample) * source_width + cell];

    int best_alt = -1;
    long long best_score = std::numeric_limits<long long>::max();
    std::vector<int32_t> best_projected;
    for (int alt = 1; alt < old_alleles; ++alt) {
        const char* allele = record->d.allele[alt];
        if (allele == nullptr || std::strcmp(allele, "*") == 0) continue;
        const auto projected = fastgatk::kernels::remap_genotype_pl_kokkos(
            source, static_cast<std::size_t>(sample_count), old_alleles, 2, ploidy, {0, alt});
        if (projected.pl.size() < static_cast<std::size_t>(sample_count) * target_pl_width)
            continue;
        long long score = 0;
        bool valid = true;
        for (int sample = 0; sample < sample_count; ++sample) {
            const auto offset = static_cast<std::size_t>(sample) * target_pl_width;
            const auto value = projected.pl[offset + (target_pl_width > 1 ? 1 : 0)];
            if (value == bcf_int32_missing || value == bcf_int32_vector_end) {
                valid = false;
                break;
            }
            score += value;
        }
        if (!valid || (best_alt >= 0 && score >= best_score)) continue;
        best_alt = alt;
        best_score = score;
        best_projected = projected.pl;
        ++telemetry.pl_remap_calls;
        telemetry.pl_remap_prepare_seconds += projected.prepare_seconds;
        telemetry.pl_remap_execute_seconds += projected.seconds;
        telemetry.pl_remap_execution_space = projected.execution_space;
    }
    free(pl_values);
    if (best_alt < 0 || best_projected.empty()) { free(genotypes); return false; }

    bool trim_reference = record->d.allele[0] != nullptr &&
        std::strlen(record->d.allele[0]) > 1;
    for (int index = 0; index < genotype_count; ++index) {
        if (genotypes[index] == bcf_int32_vector_end) continue;
        if (bcf_gt_is_missing(genotypes[index])) { trim_reference = true; continue; }
        const int allele = bcf_gt_allele(genotypes[index]);
        if (allele >= 0 && allele < old_alleles &&
            std::strcmp(record->d.allele[allele], "*") == 0)
            trim_reference = true;
    }
    free(genotypes);
    std::string reference = record->d.allele[0] == nullptr ? "N" : record->d.allele[0];
    if (trim_reference && reference.size() > 1) reference.resize(1);
    const std::string alleles = reference + ",<NON_REF>";
    if (bcf_update_alleles_str(output_header, record, alleles.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create high-confidence reference block alleles");

    std::vector<int32_t> gt(static_cast<std::size_t>(sample_count * ploidy), bcf_gt_unphased(0));
    if (bcf_update_genotypes(output_header, record, gt.data(), static_cast<int>(gt.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set high-confidence reference GT");
    if (bcf_update_format_int32(output_header, record, "PL", best_projected.data(),
                                static_cast<int>(best_projected.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set high-confidence reference PL");
    update_gq_from_pl(output_header, record);
    clear_format_attributes_except_core(output_header, record);
    // Re-add MIN_DP before PL to match the stable GATK reference-block
    // FORMAT order (GT:DP:GQ:MIN_DP:PL), even when the source concrete
    // variant used GT:AD:DP:GQ:PL:SB.
    if (bcf_update_format(output_header, record, "PL", nullptr, 0, BCF_HT_STR) != 0 ||
        bcf_update_format(output_header, record, "MIN_DP", nullptr, 0, BCF_HT_STR) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot reorder reference-block FORMAT fields");
    if (fallback_min_dp >= 0 || !source_min_dps.empty()) {
        std::vector<int32_t> min_dps(static_cast<std::size_t>(sample_count), fallback_min_dp);
        for (int sample = 0; sample < sample_count && sample < static_cast<int>(source_min_dps.size()); ++sample)
            if (source_min_dps[static_cast<std::size_t>(sample)] >= 0)
                min_dps[static_cast<std::size_t>(sample)] = source_min_dps[static_cast<std::size_t>(sample)];
        set_format_values(output_header, record, "MIN_DP", min_dps);
    }
    if (bcf_update_format_int32(output_header, record, "PL", best_projected.data(),
                                static_cast<int>(best_projected.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot restore high-confidence reference PL");
    order_reference_block_format(output_header, record);
    clear_info_except_end(output_header, record);
    const int32_t end_value = end;
    if (bcf_update_info_int32(output_header, record, "END", &end_value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set high-confidence reference END");
    if (bcf_update_filter(output_header, record, nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear high-confidence reference FILTER");
    std::uint32_t missing_bits = bcf_float_missing;
    std::memcpy(&record->qual, &missing_bits, sizeof(record->qual));
    ++telemetry.converted_high_confidence_variants;
    return true;
}

// ReblockGVCF calls GATKVariantContextUtils.trimAlleles(builder.make(),
// false, true) after unused concrete ALTs have been removed.  The relevant
// operation here is reverse trimming: remove the longest common suffix while
// keeping at least one base in every non-symbolic allele.  Symbolic and
// spanning-deletion alleles are intentionally left untouched.  GT/AD/PL
// indices do not change because only allele strings and the reference span
// are normalized.
bool trim_common_suffix(const bcf_hdr_t* header, bcf1_t* record, int& output_end) {
    if (record->n_allele <= 1) return false;
    std::vector<int> concrete;
    concrete.reserve(static_cast<std::size_t>(record->n_allele));
    std::size_t shortest = std::numeric_limits<std::size_t>::max();
    for (int index = 0; index < record->n_allele; ++index) {
        const char* allele = record->d.allele[index];
        if (allele == nullptr || allele[0] == '<' || std::strcmp(allele, "*") == 0)
            continue;
        const auto length = std::strlen(allele);
        // GATK returns the input unchanged when any concrete allele is a
        // single base; this also prevents an empty allele after clipping.
        if (length <= 1) return false;
        concrete.push_back(index);
        shortest = std::min(shortest, length);
    }
    if (concrete.empty() || shortest <= 1) return false;

    std::size_t trim = 0;
    while (trim + 1 < shortest) {
        const char suffix = record->d.allele[concrete.front()][
            std::strlen(record->d.allele[concrete.front()]) - trim - 1];
        bool matches = true;
        for (const auto index : concrete) {
            const auto* allele = record->d.allele[index];
            if (allele[std::strlen(allele) - trim - 1] != suffix) {
                matches = false;
                break;
            }
        }
        if (!matches) break;
        ++trim;
    }
    if (trim == 0) return false;

    std::string allele_string;
    for (int index = 0; index < record->n_allele; ++index) {
        if (index != 0) allele_string.push_back(',');
        const auto* allele = record->d.allele[index];
        if (allele == nullptr || allele[0] == '<' || std::strcmp(allele, "*") == 0) {
            if (allele != nullptr) allele_string += allele;
            continue;
        }
        allele_string.append(allele, std::strlen(allele) - trim);
    }
    if (bcf_update_alleles_str(header, record, allele_string.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot reverse-trim reblocked alleles");
    bcf_unpack(record, BCF_UN_ALL);
    // A concrete variant without an explicit END spans the trimmed REF.  An
    // existing gVCF END is retained by the caller because it may encode a
    // caller-specific covered interval rather than the allele stop.
    output_end = record->pos + static_cast<int>(std::strlen(record->d.allele[0]));
    return true;
}

// When reverse trimming removes the tail of a deletion, GATK inserts a small
// reference block for the uncovered tail.  Build that block from the original
// genotype likelihoods projected to REF versus the shortest dropped concrete
// ALT.  This keeps the PL math on the shared Kokkos Number=G path; Host only
// performs FAIDX lookup and HTSlib record materialization.
bcf1_t* make_deletion_gap_ref_block(const bcf_hdr_t* input_header,
                                    const bcf_hdr_t* output_header,
                                    const bcf1_t* original,
                                    faidx_t* reference,
                                    int gap_start_one_based,
                                    int gap_end_one_based,
                                    int dropped_allele,
                                    GenotypeKernelTelemetry& telemetry) {
    if (reference == nullptr || gap_start_one_based > gap_end_one_based ||
        dropped_allele <= 0 || dropped_allele >= original->n_allele)
        return nullptr;
    const char* contig = bcf_hdr_int2id(output_header, BCF_DT_CTG, original->rid);
    if (contig == nullptr) return nullptr;
    int fetched_length = 0;
    char* fetched = faidx_fetch_seq(reference, contig, gap_start_one_based - 1,
                                    gap_start_one_based - 1, &fetched_length);
    if (fetched == nullptr || fetched_length < 1) {
        free(fetched);
        return nullptr;
    }
    const char reference_base = static_cast<char>(std::toupper(
        static_cast<unsigned char>(fetched[0])));
    free(fetched);

    auto* mutable_original = const_cast<bcf1_t*>(original);
    int32_t* genotype_values = nullptr;
    int genotype_count = 0;
    const auto genotype_length = bcf_get_genotypes(input_header, mutable_original,
                                                   &genotype_values, &genotype_count);
    const int sample_count = std::max(1, input_header->n[BCF_DT_SAMPLE]);
    const int ploidy = genotype_length > 0 && genotype_count > 0 && genotype_count % sample_count == 0
        ? genotype_count / sample_count : 2;
    free(genotype_values);
    const auto source_width = genotype_width(original->n_allele, ploidy);
    if (source_width == 0) return nullptr;

    int32_t* pl_values = nullptr;
    int pl_count = 0;
    const auto pl_length = bcf_get_format_int32(input_header, mutable_original, "PL",
                                                &pl_values, &pl_count);
    if (pl_length <= 0 || pl_count < static_cast<int>(sample_count * source_width)) {
        free(pl_values);
        return nullptr;
    }
    const std::vector<std::int32_t> source(
        pl_values, pl_values + static_cast<std::size_t>(sample_count) * source_width);
    free(pl_values);
    const auto projected = fastgatk::kernels::remap_genotype_pl_kokkos(
        source, static_cast<std::size_t>(sample_count), original->n_allele, 2, ploidy,
        {0, static_cast<std::int32_t>(dropped_allele)});
    ++telemetry.pl_remap_calls;
    telemetry.pl_remap_prepare_seconds += projected.prepare_seconds;
    telemetry.pl_remap_execute_seconds += projected.seconds;
    telemetry.pl_remap_execution_space = projected.execution_space;
    if (projected.pl.empty()) return nullptr;
    std::vector<std::int32_t> normalized = projected.pl;
    std::vector<int32_t> gqs(static_cast<std::size_t>(sample_count), 0);
    for (int sample = 0; sample < sample_count; ++sample) {
        auto begin = normalized.begin() + static_cast<std::size_t>(sample) * projected.pl.size() / sample_count;
        auto end = begin + projected.pl.size() / sample_count;
        std::vector<int32_t> finite;
        for (auto iter = begin; iter != end; ++iter)
            if (*iter >= 0 && *iter != bcf_int32_missing && *iter != bcf_int32_vector_end)
                finite.push_back(*iter);
        if (finite.empty()) {
            std::fill(begin, end, 0);
            continue;
        }
        const auto minimum = *std::min_element(finite.begin(), finite.end());
        for (auto iter = begin; iter != end; ++iter)
            if (*iter >= 0) *iter -= minimum;
        std::sort(finite.begin(), finite.end());
        gqs[static_cast<std::size_t>(sample)] = finite.size() > 1
            ? std::min(99, std::max(0, static_cast<int>(finite[1] - finite[0]))) : 0;
    }

    auto* block = bcf_init();
    if (block == nullptr) return nullptr;
    block->rid = original->rid;
    block->pos = gap_start_one_based - 1;
    const std::string alleles = std::string(1, reference_base) + ",<NON_REF>";
    if (bcf_update_alleles_str(output_header, block, alleles.c_str()) != 0) {
        bcf_destroy(block);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create deletion-gap reference block alleles");
    }
    std::vector<int32_t> hom_ref(static_cast<std::size_t>(sample_count * ploidy), bcf_gt_unphased(0));
    if (bcf_update_genotypes(output_header, block, hom_ref.data(), static_cast<int>(hom_ref.size())) != 0 ||
        (bcf_hdr_id2int(output_header, BCF_DT_ID, "PL") >= 0 &&
         bcf_update_format_int32(output_header, block, "PL", normalized.data(),
                                  static_cast<int>(normalized.size())) != 0)) {
        bcf_destroy(block);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize deletion-gap PL/GT");
    }
    set_format_values(output_header, block, "GQ", gqs);
    const int info_depth = info_int(input_header, mutable_original, "DP", -1);
    std::vector<int32_t> depths(static_cast<std::size_t>(sample_count), 0);
    for (int sample = 0; sample < sample_count; ++sample) {
        const int depth = info_depth >= 0 ? info_depth :
            format_int(input_header, mutable_original, "DP", 0, sample);
        depths[static_cast<std::size_t>(sample)] = std::max(0, depth);
    }
    set_format_values(output_header, block, "DP", depths);
    set_format_values(output_header, block, "MIN_DP", depths);
    const int end = gap_end_one_based;
    if (bcf_update_info_int32(output_header, block, "END", &end, 1) != 0) {
        bcf_destroy(block);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update deletion-gap END");
    }
    return block;
}

int band_for_gq(const Options& options, int gq) {
    for (std::size_t index = 0; index < options.gq_bands.size(); ++index)
        if (gq < options.gq_bands[index]) return static_cast<int>(index);
    return static_cast<int>(options.gq_bands.size() - 1);
}

int lower_bound_for_band(const Options& options, int band) {
    return band <= 0 ? 0 : options.gq_bands[static_cast<std::size_t>(band - 1)];
}

// GVCFBlock keeps an element-wise minimum PL across every hom-ref record in a
// band.  Do this at the HTSlib boundary, after the values have been prepared
// by the caller, and leave all genotype-array remapping on the shared Kokkos
// path for allele-changing operations.
bool merge_reference_block_pl(const bcf_hdr_t* header, bcf1_t* previous,
                              const bcf1_t* current, int sample_count,
                              std::vector<int32_t>& gqs) {
    if (bcf_hdr_id2int(header, BCF_DT_ID, "PL") < 0) return false;
    int32_t* previous_values = nullptr;
    int previous_count = 0;
    int32_t* current_values = nullptr;
    int current_count = 0;
    const auto previous_length = bcf_get_format_int32(
        header, previous, "PL", &previous_values, &previous_count);
    const auto current_length = bcf_get_format_int32(
        header, const_cast<bcf1_t*>(current), "PL", &current_values, &current_count);
    const int rows = std::max(1, sample_count);
    const int previous_width = previous_length > 0 && previous_count >= rows &&
        previous_count % rows == 0 ? previous_count / rows : 0;
    const int current_width = current_length > 0 && current_count >= rows &&
        current_count % rows == 0 ? current_count / rows : 0;
    if (previous_width <= 0 || previous_width != current_width) {
        free(previous_values);
        free(current_values);
        return false;
    }
    std::vector<int32_t> merged(static_cast<std::size_t>(rows * previous_width));
    for (int sample = 0; sample < rows; ++sample) {
        for (int cell = 0; cell < previous_width; ++cell) {
            const auto left = previous_values[sample * previous_width + cell];
            const auto right = current_values[sample * current_width + cell];
            if (left == bcf_int32_missing || left == bcf_int32_vector_end)
                merged[static_cast<std::size_t>(sample * previous_width + cell)] = right;
            else if (right == bcf_int32_missing || right == bcf_int32_vector_end)
                merged[static_cast<std::size_t>(sample * previous_width + cell)] = left;
            else
                merged[static_cast<std::size_t>(sample * previous_width + cell)] = std::min(left, right);
        }
    }
    free(previous_values);
    free(current_values);
    if (bcf_update_format_int32(header, previous, "PL", merged.data(),
                                static_cast<int>(merged.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge reference-block PL");
    gqs.assign(static_cast<std::size_t>(rows), 0);
    for (int sample = 0; sample < rows; ++sample) {
        const auto gq = pl_gq(header, previous, sample);
        if (gq >= 0) gqs[static_cast<std::size_t>(sample)] = gq;
    }
    set_format_values(header, previous, "GQ", gqs);
    return true;
}

int rounded_median(std::vector<int32_t> values) {
    if (values.empty()) return -1;
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    double median = static_cast<double>(values[middle]);
    if (values.size() % 2 == 0)
        median = (static_cast<double>(values[middle - 1]) +
                  static_cast<double>(values[middle])) / 2.0;
    return static_cast<int>(std::lround(median));
}

void append_dp_history(RegionRecord& previous, const RegionRecord& current,
                       const bcf_hdr_t* header) {
    if (previous.dp_history.size() != current.dp_history.size()) return;
    for (std::size_t sample = 0; sample < previous.dp_history.size(); ++sample)
        previous.dp_history[sample].insert(previous.dp_history[sample].end(),
                                           current.dp_history[sample].begin(),
                                           current.dp_history[sample].end());
    std::vector<int32_t> medians;
    medians.reserve(previous.dp_history.size());
    for (const auto& values : previous.dp_history)
        medians.push_back(rounded_median(values));
    if (!medians.empty() && std::any_of(medians.begin(), medians.end(), [](int32_t value) { return value >= 0; }))
        set_format_values(header, previous.value, "DP", medians);
}

// ReblockingGVCFBlockCombiner buffers reference blocks because a concrete
// variant (in particular a deletion) may overlap one or more input blocks.
// Before GVCFBlock can merge the blocks, Java trims the covered portion and
// splits a block into left/right pieces when the variant is internal.  The
// old native path sorted records and merged blocks directly, which could emit
// a reference block over a concrete call.  Keep the operation at the HTSlib
// record boundary, but obtain the new REF base from the indexed reference so
// the same record remains valid after a left trim/split.
void move_reference_block_start(const bcf_hdr_t* header, bcf1_t* record,
                                faidx_t* reference, int new_start_one_based) {
    if (reference == nullptr)
        throw std::runtime_error(
            "BAD_INPUT: reference FASTA is required to split an overlapping reference block");
    const char* contig = bcf_hdr_int2id(header, BCF_DT_CTG, record->rid);
    if (contig == nullptr)
        throw std::runtime_error("BAD_INPUT: reference block contig is absent from output header");
    int fetched_length = 0;
    char* fetched = faidx_fetch_seq(reference, contig, new_start_one_based - 1,
                                    new_start_one_based - 1, &fetched_length);
    if (fetched == nullptr || fetched_length < 1) {
        free(fetched);
        throw std::runtime_error("BAD_INPUT: cannot fetch REF base while splitting overlapping reference block");
    }
    const char ref_base = static_cast<char>(std::toupper(
        static_cast<unsigned char>(fetched[0])));
    free(fetched);
    const std::string alleles = std::string(1, ref_base) + ",<NON_REF>";
    if (bcf_update_alleles_str(header, record, alleles.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split reference block REF");
    record->pos = new_start_one_based - 1;
}

void update_reference_block_span(const bcf_hdr_t* header, RegionRecord& block,
                                 int new_start_one_based, int new_end_one_based,
                                 faidx_t* reference) {
    if (new_start_one_based > new_end_one_based)
        throw std::invalid_argument("invalid empty reference block after overlap trimming");
    const int old_start_one_based = static_cast<int>(block.value->pos) + 1;
    if (new_start_one_based != old_start_one_based)
        move_reference_block_start(header, block.value, reference, new_start_one_based);
    const int32_t end_value = new_end_one_based;
    if (bcf_update_info_int32(header, block.value, "END", &end_value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update trimmed reference block END");
    block.pos = new_start_one_based - 1;
    block.end = new_end_one_based;
}

struct ReblockVariantSpan {
    int rid = -1;
    int begin = 0;  // zero-based inclusive
    int end = 0;    // one-based inclusive, represented as a coordinate value
};

// Trim/split reference blocks against concrete output records before the
// normal GQ-band combiner.  The pass is deliberately independent of input
// order: shards and VCF readers can expose a deletion before or after a block
// at the same locus, while Java's buffered writer produces the same final
// coverage in either case.
void trim_overlapping_reference_blocks(
    std::vector<RegionRecord>& records, const bcf_hdr_t* header,
    faidx_t* reference, GenotypeKernelTelemetry& telemetry) {
    std::vector<ReblockVariantSpan> variants;
    variants.reserve(records.size());
    for (const auto& item : records) {
        if (item.reference_block) continue;
        variants.push_back(ReblockVariantSpan{
            item.rid, item.pos, std::max(item.pos + 1, item.end)});
    }
    if (variants.empty()) return;

    std::vector<RegionRecord> trimmed;
    trimmed.reserve(records.size());
    for (auto& item : records) {
        if (!item.reference_block) {
            trimmed.push_back(std::move(item));
            continue;
        }
        std::vector<RegionRecord> pieces;
        pieces.push_back(std::move(item));
        for (const auto& variant : variants) {
            std::vector<RegionRecord> next;
            next.reserve(pieces.size() + 1);
            for (auto& piece : pieces) {
                const int block_begin = piece.pos;
                const int block_end = piece.end;
                if (piece.rid != variant.rid || variant.begin >= block_end ||
                    variant.end <= block_begin) {
                    next.push_back(std::move(piece));
                    continue;
                }
                const int variant_begin = variant.begin + 1;
                const int variant_end = variant.end;
                if (variant_begin <= block_begin + 1 && variant_end >= block_end) {
                    ++telemetry.overlapping_ref_block_dropped;
                    bcf_destroy(piece.value);
                    continue;
                }
                if (variant_begin <= block_begin + 1) {
                    const int new_start = variant_end + 1;
                    if (new_start <= block_end) {
                        update_reference_block_span(header, piece, new_start, block_end, reference);
                        ++telemetry.overlapping_ref_block_trimmed;
                        next.push_back(std::move(piece));
                    } else {
                        ++telemetry.overlapping_ref_block_dropped;
                        bcf_destroy(piece.value);
                    }
                    continue;
                }
                if (variant_end >= block_end) {
                    update_reference_block_span(header, piece, block_begin + 1,
                                                variant_begin - 1, reference);
                    ++telemetry.overlapping_ref_block_trimmed;
                    next.push_back(std::move(piece));
                    continue;
                }

                // The variant is strictly internal: retain the left piece and
                // duplicate it for the right tail, which needs a new REF
                // base and position just like Java moveBuilderStart().
                RegionRecord right;
                right.value = bcf_dup(piece.value);
                if (right.value == nullptr)
                    throw std::runtime_error(
                        "RESOURCE_EXHAUSTED: cannot duplicate split reference block");
                right.rid = piece.rid;
                right.pos = piece.pos;
                right.end = piece.end;
                right.gq = piece.gq;
                right.min_dp = piece.min_dp;
                right.reference_block = true;
                right.gqs = piece.gqs;
                right.min_dps = piece.min_dps;
                right.dp_history = piece.dp_history;
                update_reference_block_span(header, piece, block_begin + 1,
                                            variant_begin - 1, reference);
                update_reference_block_span(header, right, variant_end + 1,
                                            block_end, reference);
                ++telemetry.overlapping_ref_block_split;
                next.push_back(std::move(piece));
                next.push_back(std::move(right));
            }
            pieces.swap(next);
            if (pieces.empty()) break;
        }
        for (auto& piece : pieces) trimmed.push_back(std::move(piece));
    }
    records.swap(trimmed);
    std::sort(records.begin(), records.end(), [](const RegionRecord& left,
                                                 const RegionRecord& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.pos != right.pos) return left.pos < right.pos;
        return left.end < right.end;
    });
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    bcf_hdr_t* output_header = nullptr;
    std::vector<RegionRecord> records;
    std::uint64_t input_records = 0;
    std::uint64_t converted_low_quality = 0;
    std::uint64_t dropped_low_quality_blocks = 0;
    std::uint64_t dropped_low_quality_variants = 0;
    std::uint64_t converted_tree_score = 0;
    std::uint64_t compacted_alleles = 0;
    GenotypeKernelTelemetry genotype_kernel_telemetry;
    std::uint64_t merged_blocks = 0;
    std::uint64_t interval_skipped = 0;
    fastgatk::io::IntervalFileStats interval_file_stats;
    int sample_count = -1;
    std::vector<std::string> sample_names;
    std::vector<Region> regions;
    faidx_t* reference_index = nullptr;
    try {
        if (!options.reference.empty()) {
            reference_index = fai_load(options.reference.c_str());
            if (reference_index == nullptr)
                throw std::runtime_error("BAD_INPUT: cannot load reference FASTA/FAI: " + options.reference);
        }
        for (const auto& input_path : options.inputs) {
            htsFile* input = bcf_open(input_path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open GVCF: " + input_path);
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) { bcf_close(input); throw std::runtime_error("BAD_INPUT: cannot read GVCF header: " + input_path); }
            if (options.tree_score_threshold > 0.0 &&
                bcf_hdr_id2int(header, BCF_DT_ID, "TREE_SCORE") < 0) {
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error("BAD_INPUT: --tree-score-threshold-to-no-call requires INFO/TREE_SCORE");
            }
            if (sample_count < 0) {
                sample_count = header->n[BCF_DT_SAMPLE];
                sample_names.reserve(static_cast<std::size_t>(std::max(0, sample_count)));
                for (int sample = 0; sample < sample_count; ++sample) {
                    const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                    sample_names.emplace_back(name == nullptr ? "" : name);
                }
            }
            bool sample_header_matches = header->n[BCF_DT_SAMPLE] == sample_count;
            for (int sample = 0; sample_header_matches && sample < sample_count; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (name == nullptr || sample >= static_cast<int>(sample_names.size()) ||
                    sample_names[static_cast<std::size_t>(sample)] != name)
                    sample_header_matches = false;
            }
            if (!sample_header_matches) {
                bcf_hdr_destroy(header); bcf_close(input);
                throw std::runtime_error("BAD_INPUT: ReblockGVCF shard sample headers differ");
            }
            if (!output_header) {
                output_header = bcf_hdr_dup(header);
                if (!output_header) { bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate GVCF header"); }
                for (const auto& annotation : options.annotations_to_keep) {
                    if (!has_header_field(output_header, BCF_HL_INFO, annotation.c_str())) {
                        bcf_hdr_destroy(output_header); output_header = nullptr;
                        bcf_hdr_destroy(header); bcf_close(input);
                        throw std::runtime_error("BAD_INPUT: INFO annotation not found: " + annotation);
                    }
                    // GATK ReblockGVCF only carries fixed-length caller
                    // annotations through allele compaction.  Number=A
                    // fields would describe dropped/reordered ALTs and are
                    // therefore rejected instead of being silently stale.
                    if (info_annotation_is_allele_specific(output_header, annotation)) {
                        bcf_hdr_destroy(output_header); output_header = nullptr;
                        bcf_hdr_destroy(header); bcf_close(input);
                        throw std::runtime_error(
                            "BAD_INPUT: allele-specific INFO annotation cannot be kept: " + annotation);
                    }
                }
                for (const auto& annotation : options.format_annotations_to_remove) {
                    if (!has_header_field(output_header, BCF_HL_FMT, annotation.c_str())) {
                        bcf_hdr_destroy(output_header); output_header = nullptr;
                        bcf_hdr_destroy(header); bcf_close(input);
                        throw std::runtime_error("BAD_INPUT: FORMAT annotation not found: " + annotation);
                    }
                }
                bcf_hdr_append(output_header, "##fastgatk_reblock_gvcf_status=prototype-gq-blocking");
                if (!has_header_field(output_header, BCF_HL_FMT, "MIN_DP"))
                    bcf_hdr_append(output_header, "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum DP observed within the GVCF block>");
                if (!has_header_field(output_header, BCF_HL_FMT, "RGQ"))
                    bcf_hdr_append(output_header, "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=Unconditional reference genotype quality>");
                if (options.add_filters_to_genotype && !has_header_field(output_header, BCF_HL_FMT, "FT"))
                    bcf_hdr_append(output_header, "##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype-level filter>");
                if (options.do_qual_approx && !has_header_field(output_header, BCF_HL_INFO, "QUALapprox"))
                    bcf_hdr_append(output_header, "##INFO=<ID=QUALapprox,Number=1,Type=Integer,Description=Sum of PL[0] values; used to approximate QUAL>");
                if (options.do_qual_approx && !has_header_field(output_header, BCF_HL_INFO, "AS_QUALapprox"))
                    bcf_hdr_append(output_header, "##INFO=<ID=AS_QUALapprox,Number=A,Type=String,Description=Allele-specific QUAL approximations>");
                if (options.do_qual_approx && !has_header_field(output_header, BCF_HL_INFO, "VarDP"))
                    bcf_hdr_append(output_header, "##INFO=<ID=VarDP,Number=1,Type=Integer,Description=(informative) depth over variant genotypes>");
                if (options.do_qual_approx && !has_header_field(output_header, BCF_HL_INFO, "AS_VarDP"))
                    bcf_hdr_append(output_header, "##INFO=<ID=AS_VarDP,Number=1,Type=String,Description=\"Allele-specific (informative) depth over variant genotypes -- including ref, RAW format\">");
                if (!has_header_field(output_header, BCF_HL_INFO, "RAW_GT_COUNT"))
                    bcf_hdr_append(output_header, "##INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=Raw genotype counts for ExcessHet>");
                // These headers are always emitted by GATK's annotation
                // engine, even when the input only has legacy MQ fields.
                if (!has_header_field(output_header, BCF_HL_INFO, "RAW_MQandDP"))
                    bcf_hdr_append(output_header, "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation>");
                if (!has_header_field(output_header, BCF_HL_INFO, "MQ_DP"))
                    bcf_hdr_append(output_header, "##INFO=<ID=MQ_DP,Number=1,Type=Integer,Description=Depth for RMS Mapping Quality calculation (deprecated)>");
                if (bcf_hdr_sync(output_header) != 0) { bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync GVCF header"); }
                bool first_selector = true;
                for (const auto& text : options.regions) {
                    append_reblock_region_selector_with_rule(
                        text, output_header, regions, interval_file_stats,
                        options.interval_set_rule, first_selector);
                    first_selector = false;
                }
                normalize_reblock_regions(regions);
            }
            bcf1_t* record = bcf_init();
            if (!record) { bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed"); }
            while (bcf_read(input, header, record) == 0) {
                ++input_records;
                bcf_unpack(record, BCF_UN_ALL);
                if (record->rid < 0 || record->n_allele < 2) continue;
                const int explicit_end = info_int(header, record, "END", -1);
                const int end = explicit_end >= 0 ? explicit_end :
                    record->pos + static_cast<int>(std::strlen(record->d.allele[0]));
                const int original_ref_length = static_cast<int>(std::strlen(record->d.allele[0]));
                if (!overlaps_regions(record, regions, end)) {
                    ++interval_skipped;
                    continue;
                }
                auto* copy = bcf_dup(record);
                if (!copy) { bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_dup failed"); }
                bcf_unpack(copy, BCF_UN_ALL);
                const int sample_rows = std::max(1, sample_count);
                const bool reference_block = has_non_ref(record) && !has_concrete_alt(record);
                std::vector<int32_t> gqs(static_cast<std::size_t>(sample_rows), -1);
                std::vector<int32_t> rgqs(static_cast<std::size_t>(sample_rows), -1);
                std::vector<int32_t> depths(static_cast<std::size_t>(sample_rows), 0);
                std::vector<int32_t> min_dps(static_cast<std::size_t>(sample_rows), -1);
                for (int sample = 0; sample < sample_rows; ++sample) {
                    gqs[static_cast<std::size_t>(sample)] = genotype_gq(header, record, sample);
                    rgqs[static_cast<std::size_t>(sample)] = format_pl0(header, record, -1, sample);
                    const int info_depth = info_int(header, record, "DP", -1);
                    depths[static_cast<std::size_t>(sample)] = format_int(
                        header, record, "DP", info_depth, sample);
                    // GATK's reblocking writer recomputes MIN_DP from the
                    // per-site DP values for an already-compressed hom-ref
                    // block; an input MIN_DP is only an annotation from the
                    // previous block encoding.  Newly converted variants
                    // still use the caller's MIN_DP when available.
                    min_dps[static_cast<std::size_t>(sample)] = reference_block
                        ? depths[static_cast<std::size_t>(sample)]
                        : format_int(header, record, "MIN_DP", info_int(header, record, "MIN_DP",
                          depths[static_cast<std::size_t>(sample)]), sample);
                }
                const int gq = gqs.empty() ? -1 : *std::min_element(gqs.begin(), gqs.end());
                const int depth = depths.empty() ? 0 : *std::min_element(depths.begin(), depths.end());
                const int min_dp = min_dps.empty() ? -1 : *std::min_element(min_dps.begin(), min_dps.end());
                const bool low_quality = !rgqs.empty() && std::all_of(
                    rgqs.begin(), rgqs.end(), [&](const int value) {
                        return value >= 0 && value < options.rgq_threshold;
                    });
                const bool tree_low_quality = !reference_block && options.tree_score_threshold > 0.0 &&
                    info_float(header, record, "TREE_SCORE", 0.0) < options.tree_score_threshold;
                // Existing hom-ref blocks already carry the caller's PL/GQ
                // summary and must flow into GVCFBlock unchanged.  Only a
                // variant that is newly converted because of a threshold is
                // rewritten to the synthetic [0,0,0] hom-ref call.
                // GATK applies the explicit RGQ/PL[0] threshold independently
                // of --drop-low-quals.  The latter controls dropping GQ0
                // records, whereas --rgq-threshold-to-no-call converts a
                // low-reference-confidence variant into a GQ0 hom-ref block.
                const bool convert_low_quality = (options.drop_low_quals ||
                    options.rgq_threshold > 0) && low_quality && !reference_block;
                if (convert_low_quality || tree_low_quality) {
                    if (convert_low_quality) ++converted_low_quality;
                    if (tree_low_quality) ++converted_tree_score;
                    convert_to_ref_block(output_header, copy, end, min_dps,
                                         min_dp >= 0 ? min_dp : depth);
                } else if (options.drop_low_quals &&
                           std::all_of(gqs.begin(), gqs.end(), [](const int value) { return value < 0; }) &&
                           !options.allow_missing_hom_ref_data) {
                    bcf_destroy(copy);
                    continue;
                }
                bcf_unpack(copy, BCF_UN_ALL);
                bool now_block = has_non_ref(copy) && !has_concrete_alt(copy);
                // ReblockingGVCF.regenotypeVC drops an existing hom-ref
                // block when --drop-low-quals is enabled and its GQ is zero
                // (or below the explicit RGQ threshold).  Do this before
                // global block merging; otherwise a low-quality span would
                // incorrectly bridge two retained blocks.
                const bool drop_existing_block = options.drop_low_quals && now_block &&
                    std::all_of(gqs.begin(), gqs.end(), [&](const int value) {
                        return value < 0 || value == 0 || value < options.rgq_threshold;
                    });
                if (drop_existing_block) {
                    bcf_destroy(copy);
                    ++dropped_low_quality_blocks;
                    continue;
                }
                if (options.drop_low_quals && !now_block && has_concrete_alt(copy)) {
                    // GATK re-genotypes concrete variants in drop mode with
                    // the standard confidence threshold.  A lightweight
                    // equivalent at this boundary uses PL[0] and the
                    // PL-derived GQ: low-confidence ALT calls are dropped,
                    // while a confidently hom-ref PL row is projected to a
                    // REF/<NON_REF> block below.  Explicit RGQ threshold
                    // conversions above still take precedence.
                    long long reference_confidence = 0;
                    bool any_non_ref_best = false;
                    bool all_ref_best = true;
                    bool all_ref_high_confidence = true;
                    for (int sample = 0; sample < sample_rows; ++sample) {
                        const int pl0 = format_pl0(header, record, -1, sample);
                        if (pl0 >= 0) reference_confidence += pl0;
                        const bool ref_best = pl_best_genotype_is_ref(header, record, sample);
                        if (!ref_best) {
                            any_non_ref_best = true;
                            all_ref_best = false;
                        } else if (pl_gq(header, record, sample) < 30) {
                            all_ref_high_confidence = false;
                        }
                    }
                    const bool low_confidence_variant =
                        (any_non_ref_best && reference_confidence < 30) ||
                        (all_ref_best && !all_ref_high_confidence);
                    if (low_confidence_variant) {
                        bcf_destroy(copy);
                        ++dropped_low_quality_variants;
                        continue;
                    }
                    if (all_ref_best && all_ref_high_confidence) {
                        // GATK's genotyping engine can turn a concrete site
                        // whose best call is confidently hom-ref into a
                        // reference block.  Preserve the projected PL/GQ
                        // (rather than the all-zero low-quality encoding).
                        convert_high_confidence_homref_to_ref_block(
                            header, output_header, copy, end, min_dps,
                            min_dp >= 0 ? min_dp : depth, genotype_kernel_telemetry);
                        // The conversion changes the allele set in-place;
                        // refresh the classification before annotation and
                        // block-merging logic below.
                        bcf_unpack(copy, BCF_UN_ALL);
                        now_block = has_non_ref(copy) && !has_concrete_alt(copy);
                    }
                }
                bool compacted_this_record = false;
                int shortest_dropped_allele = -1;
                int output_end = end;
                if (now_block) {
                    // GATK's low-quality hom-ref conversion is unfiltered;
                    // site FILTERs are not carried onto reference blocks.
                    if (bcf_update_filter(output_header, copy, nullptr, 0) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear ref-block FILTER");
                } else {
                    const auto filters = site_filter_string(output_header, copy);
                    if (options.add_filters_to_genotype) {
                        set_format_string_values(output_header, copy, "FT",
                            std::vector<std::string>(static_cast<std::size_t>(sample_rows), filters));
                    }
                    if (!options.keep_filters &&
                        bcf_update_filter(output_header, copy, nullptr, 0) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear variant FILTER");
                }
                if (!now_block)
                    update_mq_annotations(output_header, copy);
                remove_reblock_info_annotations(output_header, copy, options,
                                                genotype_kernel_telemetry);
                if (!options.keep_all_alts && !now_block && copy->n_allele > 2) {
                    std::vector<int> kept{0};
                    int32_t* genotypes = nullptr;
                    int genotype_count = 0;
                    if (bcf_get_genotypes(header, record, &genotypes, &genotype_count) > 0) {
                        for (int index = 0; index < genotype_count; ++index) {
                            if (bcf_gt_is_missing(genotypes[index])) continue;
                            const auto allele = bcf_gt_allele(genotypes[index]);
                            if (allele > 0 && std::find(kept.begin(), kept.end(), allele) == kept.end()) kept.push_back(allele);
                        }
                    }
                    free(genotypes);
                    if (std::find(kept.begin(), kept.end(), copy->n_allele - 1) == kept.end() &&
                        std::strcmp(copy->d.allele[copy->n_allele - 1], "<NON_REF>") == 0)
                        kept.push_back(copy->n_allele - 1);
                    if (kept.size() < static_cast<std::size_t>(copy->n_allele)) {
                        const int old_alleles = copy->n_allele;
                        for (int allele = 1; allele < old_alleles; ++allele) {
                            if (std::find(kept.begin(), kept.end(), allele) == kept.end() &&
                                copy->d.allele[allele] != nullptr && copy->d.allele[allele][0] != '<' &&
                                std::strcmp(copy->d.allele[allele], "*") != 0 &&
                                (shortest_dropped_allele < 0 ||
                                 std::strlen(copy->d.allele[allele]) <
                                 std::strlen(copy->d.allele[shortest_dropped_allele])))
                                shortest_dropped_allele = allele;
                        }
                        std::string allele_string;
                        for (std::size_t index = 0; index < kept.size(); ++index) {
                            if (index != 0) allele_string.push_back(',');
                            allele_string += copy->d.allele[kept[index]];
                        }
                        if (bcf_update_alleles_str(output_header, copy, allele_string.c_str()) != 0) {
                            bcf_destroy(copy); bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact reblocked alleles");
                        }
                        update_gt(header, output_header, copy, kept);
                        update_ad_pl(header, output_header, copy, kept, old_alleles,
                                     genotype_kernel_telemetry);
                        update_gq_from_pl(output_header, copy);
                        ++compacted_alleles;
                        compacted_this_record = true;
                    }
                }
                if (!now_block && !options.keep_all_alts && compacted_this_record &&
                    trim_common_suffix(output_header, copy, output_end)) {
                    // An explicit END describes the caller's covered span;
                    // retain it even when the normalized REF is shorter.
                    if (explicit_end >= 0) output_end = end;
                    ++genotype_kernel_telemetry.reverse_allele_trimmed_variants;
                    const int trimmed_ref_length = static_cast<int>(std::strlen(copy->d.allele[0]));
                    const int trimmed_bases = original_ref_length - trimmed_ref_length;
                    if (trimmed_bases > 0 && shortest_dropped_allele >= 0) {
                        const int original_end = explicit_end >= 0 ? end :
                            record->pos + original_ref_length;
                        int gap_start = original_end - trimmed_bases + 1;
                        gap_start = std::max(gap_start, output_end + 1);
                        if (gap_start <= original_end) {
                            auto* gap_block = make_deletion_gap_ref_block(
                                header, output_header, record, reference_index, gap_start,
                                original_end, shortest_dropped_allele, genotype_kernel_telemetry);
                            if (gap_block != nullptr) {
                                const int block_gq = format_int(output_header, gap_block, "GQ", 0);
                                const int block_min_dp = format_int(output_header, gap_block, "MIN_DP", -1);
                                RegionRecord gap_record{
                                    gap_block, gap_block->rid, static_cast<int>(gap_block->pos),
                                    original_end, std::max(0, block_gq), block_min_dp, true, {}, {}, {}};
                                gap_record.gqs = format_values(output_header, gap_block, "GQ",
                                                               sample_rows, std::max(0, block_gq));
                                gap_record.min_dps = format_values(output_header, gap_block, "MIN_DP",
                                                                  sample_rows, block_min_dp);
                                gap_record.gq = gap_record.gqs.empty() ? std::max(0, block_gq) :
                                    *std::min_element(gap_record.gqs.begin(), gap_record.gqs.end());
                                gap_record.min_dp = gap_record.min_dps.empty() ? block_min_dp :
                                    *std::min_element(gap_record.min_dps.begin(), gap_record.min_dps.end());
                                records.push_back(std::move(gap_record));
                                ++genotype_kernel_telemetry.deletion_gap_ref_blocks;
                            } else {
                                ++genotype_kernel_telemetry.deletion_gap_ref_block_fallbacks;
                            }
                        }
                    }
                }
                if (!now_block)
                    remove_non_ref_ad_kokkos(output_header, copy, genotype_kernel_telemetry);
                remove_reblock_format_annotations(output_header, copy, options,
                                                  genotype_kernel_telemetry);
                // GVCFs produced by some callers omit FORMAT/GQ and expect it
                // to be materialized from the complete triangular PL vector.
                // Do this after allele compaction so the reported GQ matches
                // the output PL ordering.
                if (format_int(output_header, copy, "GQ", -1) < 0)
                    update_gq_from_pl(output_header, copy);
                const int output_gq = format_int(output_header, copy, "GQ", gq);
                const int output_min_dp = format_int(output_header, copy, "MIN_DP",
                    min_dp >= 0 ? min_dp : depth);
                if (now_block && !min_dps.empty())
                    set_format_values(output_header, copy, "MIN_DP", min_dps);
                if (options.do_qual_approx && !now_block) {
                    long long qual_sum = 0;
                    for (int sample = 0; sample < sample_rows; ++sample) {
                        const auto value = format_pl0(output_header, copy, 0, sample);
                        if (value > 0) qual_sum += value;
                    }
                    const int qual_approx = static_cast<int>(std::clamp<long long>(
                        qual_sum, 0, std::numeric_limits<int32_t>::max()));
                    if (bcf_update_info_int32(output_header, copy, "QUALapprox", &qual_approx, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write QUALapprox");
                    const auto depth_annotations = variant_depth(output_header, copy, depth);
                    const int32_t variant_depth_value = depth_annotations.variant_depth;
                    if (bcf_update_info_int32(output_header, copy, "VarDP",
                                              &variant_depth_value, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VarDP");
                    if (depth_annotations.has_allele_depths) {
                        std::string as_depth;
                        for (std::size_t allele = 0; allele < depth_annotations.allele_depths.size(); ++allele) {
                            if (allele != 0) as_depth.push_back('|');
                            as_depth += std::to_string(depth_annotations.allele_depths[allele]);
                        }
                        if (bcf_update_info_string(output_header, copy, "AS_VarDP",
                                                   as_depth.c_str()) != 0)
                            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_VarDP");
                    }
                    ++genotype_kernel_telemetry.qual_approx_depth_annotations;
                    // GATK's AS_QUALapprox is the biallelic PL[0] after
                    // subsetting each concrete ALT; symbolic <NON_REF> and
                    // spanning-deletion alleles carry zero.  Reuse the same
                    // Kokkos Number=G remap path so the annotation remains
                    // valid for arbitrary fixed ploidy.
                    int32_t* pl_values = nullptr;
                    int pl_count = 0;
                    const auto pl_length = bcf_get_format_int32(
                        output_header, copy, "PL", &pl_values, &pl_count);
                    if (pl_length > 0 && copy->n_allele > 2 && pl_count > 0) {
                        int32_t* gt_values = nullptr;
                        int gt_count = 0;
                        const auto gt_length = bcf_get_genotypes(
                            output_header, copy, &gt_values, &gt_count);
                        const int ploidy = gt_length > 0 && gt_count > 0 &&
                            gt_count % sample_rows == 0 ? gt_count / sample_rows : 2;
                        free(gt_values);
                        const auto source_width = genotype_width(copy->n_allele, ploidy);
                        const int pl_source_width = pl_length > 0 && pl_count >= sample_rows &&
                            pl_count % sample_rows == 0 ? pl_count / sample_rows : 0;
                        if (source_width > 0 && pl_source_width >= static_cast<int>(source_width)) {
                            std::vector<int32_t> source(
                                static_cast<std::size_t>(sample_rows) * source_width,
                                bcf_int32_missing);
                            for (int sample = 0; sample < sample_rows; ++sample)
                                for (std::size_t cell = 0; cell < source_width; ++cell)
                                    source[static_cast<std::size_t>(sample) * source_width + cell] =
                                        pl_values[static_cast<std::size_t>(sample) * pl_source_width + cell];
                            // Number=A annotations reserve the leading slot
                            // for REF.  GATK serializes that non-applicable
                            // element as an empty string, so a two-ALT site
                            // is rendered as |alt1|alt2 (with symbolic
                            // <NON_REF> receiving zero below).
                            std::string as_qual("|");
                            for (int allele = 1; allele < copy->n_allele; ++allele) {
                                if (allele > 1) as_qual.push_back('|');
                                const auto* name = copy->d.allele[allele];
                                if (name != nullptr && name[0] == '<') {
                                    as_qual += "0";
                                    continue;
                                }
                                const auto projected = fastgatk::kernels::remap_genotype_pl_kokkos(
                                    source, static_cast<std::size_t>(sample_rows), copy->n_allele, 2, ploidy,
                                    {0, allele});
                                const auto target_width = genotype_width(2, ploidy);
                                long long projected_sum = 0;
                                for (int sample = 0; sample < sample_rows; ++sample) {
                                    const auto value = projected.pl.empty() ? bcf_int32_missing :
                                        projected.pl[static_cast<std::size_t>(sample) * target_width];
                                    if (value >= 0) projected_sum += value;
                                }
                                as_qual += std::to_string(std::clamp<long long>(
                                    projected_sum, 0, std::numeric_limits<int32_t>::max()));
                                ++genotype_kernel_telemetry.pl_remap_calls;
                                genotype_kernel_telemetry.pl_remap_prepare_seconds += projected.prepare_seconds;
                                genotype_kernel_telemetry.pl_remap_execute_seconds += projected.seconds;
                                genotype_kernel_telemetry.pl_remap_execution_space = projected.execution_space;
                            }
                            if (bcf_update_info_string(output_header, copy, "AS_QUALapprox",
                                                       as_qual.c_str()) != 0)
                                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_QUALapprox");
                        }
                    }
                    free(pl_values);
                }
                if (!now_block && update_raw_genotype_count(output_header, copy))
                    ++genotype_kernel_telemetry.raw_genotype_count_annotations;
                RegionRecord output_record{copy, copy->rid, static_cast<int>(copy->pos), output_end,
                                           std::max(0, output_gq), output_min_dp, now_block, {}, {}, {}};
                output_record.gqs = format_values(output_header, copy, "GQ", sample_rows,
                                                  std::max(0, output_gq));
                output_record.min_dps = format_values(output_header, copy, "MIN_DP", sample_rows,
                                                      output_min_dp);
                const auto output_dps = format_values(output_header, copy, "DP", sample_rows,
                                                      std::max(0, depth));
                output_record.dp_history.resize(static_cast<std::size_t>(sample_rows));
                for (int sample = 0; sample < sample_rows; ++sample)
                    output_record.dp_history[static_cast<std::size_t>(sample)].push_back(
                        output_dps[static_cast<std::size_t>(sample)]);
                output_record.gq = output_record.gqs.empty() ? std::max(0, output_gq) :
                    *std::min_element(output_record.gqs.begin(), output_record.gqs.end());
                output_record.min_dp = output_record.min_dps.empty() ? output_min_dp :
                    *std::min_element(output_record.min_dps.begin(), output_record.min_dps.end());
                records.push_back(std::move(output_record));
            }
            bcf_destroy(record);
            bcf_hdr_destroy(header);
            bcf_close(input);
        }
        std::sort(records.begin(), records.end(), [](const RegionRecord& left, const RegionRecord& right) {
            if (left.rid != right.rid) return left.rid < right.rid;
            if (left.pos != right.pos) return left.pos < right.pos;
            return left.end < right.end;
        });
        trim_overlapping_reference_blocks(records, output_header, reference_index,
                                          genotype_kernel_telemetry);
        std::vector<RegionRecord> merged;
        for (auto& current : records) {
            if (!merged.empty()) {
                auto& previous = merged.back();
                bool same_band = previous.gqs.size() == current.gqs.size();
                if (same_band) {
                    for (std::size_t sample = 0; sample < previous.gqs.size(); ++sample) {
                        if (band_for_gq(options, previous.gqs[sample]) !=
                            band_for_gq(options, current.gqs[sample])) {
                            same_band = false;
                            break;
                        }
                    }
                }
                const int previous_band = band_for_gq(options, previous.gq);
                if (previous.reference_block && current.reference_block && previous.rid == current.rid &&
                    same_band &&
                    // current.pos is zero-based while END is the VCF
                    // one-based inclusive coordinate.  Equality here means
                    // the records are truly adjacent; using end+1 would
                    // silently bridge a one-base uncovered hole.
                    current.pos <= previous.end) {
                    const bool merged_pl = merge_reference_block_pl(
                        output_header, previous.value, current.value,
                        static_cast<int>(previous.gqs.size()), previous.gqs);
                    if (!merged_pl && previous.gqs.size() == current.gqs.size()) {
                        for (std::size_t sample = 0; sample < previous.gqs.size(); ++sample)
                            if (current.gqs[sample] >= 0 &&
                                (previous.gqs[sample] < 0 || current.gqs[sample] < previous.gqs[sample]))
                                previous.gqs[sample] = current.gqs[sample];
                        set_format_values(output_header, previous.value, "GQ", previous.gqs);
                    }
                    append_dp_history(previous, current, output_header);
                    previous.end = std::max(previous.end, current.end);
                    if (previous.min_dps.size() == current.min_dps.size() && !previous.min_dps.empty()) {
                        for (std::size_t sample = 0; sample < previous.min_dps.size(); ++sample) {
                            if (previous.min_dps[sample] >= 0 && current.min_dps[sample] >= 0)
                                previous.min_dps[sample] = std::min(previous.min_dps[sample], current.min_dps[sample]);
                            else if (previous.min_dps[sample] < 0)
                                previous.min_dps[sample] = current.min_dps[sample];
                        }
                        set_format_values(output_header, previous.value, "MIN_DP", previous.min_dps);
                        previous.min_dp = *std::min_element(previous.min_dps.begin(), previous.min_dps.end());
                    } else if (previous.min_dp >= 0 && current.min_dp >= 0) {
                        previous.min_dp = std::min(previous.min_dp, current.min_dp);
                        set_scalar_format(output_header, previous.value, "MIN_DP", previous.min_dp);
                    } else if (previous.min_dp < 0 && current.min_dp >= 0) {
                        previous.min_dp = current.min_dp;
                        set_scalar_format(output_header, previous.value, "MIN_DP", previous.min_dp);
                    }
                    const int end = previous.end;
                    bcf_update_info_int32(output_header, previous.value, "END", &end, 1);
                    if (!previous.gqs.empty())
                        previous.gq = *std::min_element(previous.gqs.begin(), previous.gqs.end());
                    if (options.floor_blocks) {
                        if (!previous.gqs.empty()) {
                            for (auto& value : previous.gqs)
                                value = lower_bound_for_band(options, band_for_gq(options, value));
                            set_format_values(output_header, previous.value, "GQ", previous.gqs);
                            previous.gq = *std::min_element(previous.gqs.begin(), previous.gqs.end());
                        } else {
                            set_scalar_format(output_header, previous.value, "GQ",
                                lower_bound_for_band(options, previous_band));
                        }
                    }
                    bcf_destroy(current.value);
                    ++merged_blocks;
                    continue;
                }
            }
            if (current.reference_block && options.floor_blocks) {
                if (!current.gqs.empty()) {
                    for (auto& value : current.gqs)
                        value = lower_bound_for_band(options, band_for_gq(options, value));
                    set_format_values(output_header, current.value, "GQ", current.gqs);
                    current.gq = *std::min_element(current.gqs.begin(), current.gqs.end());
                } else {
                    set_scalar_format(output_header, current.value, "GQ",
                                      lower_bound_for_band(options, band_for_gq(options, current.gq)));
                }
            }
            merged.push_back(std::move(current));
        }
        records.clear();
        records = std::move(merged);
        finalize_reblock_annotation_header(output_header, options);
        if (options.floor_blocks) {
            for (auto& item : records)
                if (item.reference_block)
                    remove_floor_block_arrays(output_header, item.value);
        }
        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        htsFile* output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open reblocked output: " + options.output);
        if (bcf_hdr_write(output, output_header) != 0) { bcf_close(output); throw std::runtime_error("cannot write reblocked header"); }
        for (const auto& item : records)
            if (bcf_write(output, output_header, item.value) != 0) { bcf_close(output); throw std::runtime_error("cannot write reblocked record"); }
        bcf_close(output);
        std::string index;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build reblocked VCF index");
            } else {
                index = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index);
            }
        }
        if (!file_complete(options.output) || (!index.empty() && !file_complete(index)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: reblocked output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write ReblockGVCF manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"ReblockGVCF\",\"implementation\":\"fastgatk-reblock-gvcf\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"gvcf\","
                 << "\"compatibility\":{\"gq_bands\":true,\"reference_block_merge\":true,\"allele_compaction\":"
                 << (!options.keep_all_alts ? "true" : "false") << ",\"vcf_index\":" << (index.empty() ? "false" : "true")
                 << ",\"multi_sample\":" << (sample_count > 1 ? "true" : "false")
                 << ",\"sample_major_kokkos_remap\":true"
                 << ",\"rgq_threshold_uses_pl0\":true,\"min_dp_merge\":true,\"tree_score_threshold\":"
                 << (options.tree_score_threshold > 0.0 ? "true" : "false")
                 << ",\"qual_approx\":" << (options.do_qual_approx ? "true" : "false")
                 << ",\"qual_approx_depth\":" << (options.do_qual_approx ? "true" : "false")
                 << ",\"mq_annotations\":true"
                 << ",\"raw_genotype_count\":true"
                 << ",\"reverse_allele_trimming\":true"
                 << ",\"deletion_gap_ref_block\":true"
                 << ",\"non_ref_ad_cleanup\":true"
                 << ",\"annotation_cleanup\":true"
                 << ",\"high_confidence_homref_reblocking\":true"
                 << ",\"overlapping_ref_block_trim_split\":true"
                 << ",\"format_annotation_removal\":"
                 << (!options.format_annotations_to_remove.empty() ? "true" : "false")
                 << ",\"keep_site_filters\":" << (options.keep_filters ? "true" : "false")
                 << ",\"add_site_filters_to_genotype\":" << (options.add_filters_to_genotype ? "true" : "false")
                 << ",\"interval_subset\":" << (!regions.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"gvcf\",\"complete\":true}"
                 << (index.empty() ? "" : ",{\"path\":\"" + json_escape(index) + "\",\"kind\":\"vcf-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size() << ",\"input_records\":" << input_records
                 << ",\"output_records\":" << records.size() << ",\"merged_blocks\":" << merged_blocks
                 << ",\"converted_low_quality\":" << converted_low_quality
                 << ",\"dropped_low_quality_blocks\":" << dropped_low_quality_blocks
                 << ",\"dropped_low_quality_variants\":" << dropped_low_quality_variants
                 << ",\"converted_high_confidence_variants\":" << genotype_kernel_telemetry.converted_high_confidence_variants
                 << ",\"overlapping_ref_block_trimmed\":"
                 << genotype_kernel_telemetry.overlapping_ref_block_trimmed
                 << ",\"overlapping_ref_block_split\":"
                 << genotype_kernel_telemetry.overlapping_ref_block_split
                 << ",\"overlapping_ref_block_dropped\":"
                 << genotype_kernel_telemetry.overlapping_ref_block_dropped
                 << ",\"converted_tree_score\":" << converted_tree_score
                 << ",\"compacted_alleles\":" << compacted_alleles
                 << ",\"intervals\":" << regions.size()
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                 << ",\"pl_remap_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.pl_remap_prepare_seconds
                 << ",\"pl_remap_kernel_seconds\":"
                 << genotype_kernel_telemetry.pl_remap_execute_seconds
                 << ",\"pl_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.pl_remap_execution_space) << "\""
                 << ",\"allele_field_remap_kernel_calls\":"
                 << genotype_kernel_telemetry.allele_field_remap_calls
                 << ",\"allele_field_remap_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
                 << ",\"allele_field_remap_kernel_seconds\":"
                 << genotype_kernel_telemetry.allele_field_remap_execute_seconds
                 << ",\"allele_field_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space) << "\""
                 << ",\"raw_genotype_count_annotations\":"
                 << genotype_kernel_telemetry.raw_genotype_count_annotations
                 << ",\"qual_approx_depth_annotations\":"
                 << genotype_kernel_telemetry.qual_approx_depth_annotations
                 << ",\"reverse_allele_trimmed_variants\":"
                 << genotype_kernel_telemetry.reverse_allele_trimmed_variants
                 << ",\"deletion_gap_ref_blocks\":"
                 << genotype_kernel_telemetry.deletion_gap_ref_blocks
                 << ",\"deletion_gap_ref_block_fallbacks\":"
                 << genotype_kernel_telemetry.deletion_gap_ref_block_fallbacks
                 << ",\"non_ref_ad_zeroed\":"
                 << genotype_kernel_telemetry.non_ref_ad_zeroed
                 << ",\"non_ref_ad_kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"non_ref_ad_kernel_execution_policy\":\"RangePolicy\""
                 << ",\"non_ref_ad_kernel_calls\":"
                 << genotype_kernel_telemetry.non_ref_ad_kernel_calls
                 << ",\"non_ref_ad_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.non_ref_ad_kernel_prepare_seconds
                 << ",\"non_ref_ad_kernel_execute_seconds\":"
                 << genotype_kernel_telemetry.non_ref_ad_kernel_execute_seconds
                 << ",\"non_ref_ad_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.non_ref_ad_kernel_execution_space) << "\""
                 << ",\"info_annotations_removed\":"
                 << genotype_kernel_telemetry.info_annotations_removed
                 << ",\"format_annotations_removed\":"
                 << genotype_kernel_telemetry.format_annotations_removed
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"sample_count\":" << std::max(0, sample_count) << "}}\n";
        if (!manifest) throw std::runtime_error("cannot finalize ReblockGVCF manifest");
        std::cout << "{\"tool\":\"ReblockGVCF\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"output_records\":" << records.size()
                  << ",\"merged_blocks\":" << merged_blocks
                  << ",\"converted_tree_score\":" << converted_tree_score
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
                  << ",\"raw_genotype_count_annotations\":" << genotype_kernel_telemetry.raw_genotype_count_annotations
                  << ",\"qual_approx_depth_annotations\":" << genotype_kernel_telemetry.qual_approx_depth_annotations
                  << ",\"reverse_allele_trimmed_variants\":" << genotype_kernel_telemetry.reverse_allele_trimmed_variants
                  << ",\"deletion_gap_ref_blocks\":" << genotype_kernel_telemetry.deletion_gap_ref_blocks
                  << ",\"non_ref_ad_zeroed\":" << genotype_kernel_telemetry.non_ref_ad_zeroed
                  << ",\"converted_high_confidence_variants\":" << genotype_kernel_telemetry.converted_high_confidence_variants
                  << ",\"overlapping_ref_block_trimmed\":"
                  << genotype_kernel_telemetry.overlapping_ref_block_trimmed
                  << ",\"overlapping_ref_block_split\":"
                  << genotype_kernel_telemetry.overlapping_ref_block_split
                  << ",\"overlapping_ref_block_dropped\":"
                  << genotype_kernel_telemetry.overlapping_ref_block_dropped
                  << ",\"info_annotations_removed\":" << genotype_kernel_telemetry.info_annotations_removed
                  << ",\"format_annotations_removed\":" << genotype_kernel_telemetry.format_annotations_removed
                  << ",\"interval_skipped\":" << interval_skipped
                  << ",\"vcf_index\":" << (index.empty() ? "false" : "true")
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\"}\n";
        destroy_records(records);
        bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        destroy_records(records);
        if (output_header) bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for ReblockGVCF");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        Kokkos::initialize();
        initialized = true;
        const auto result = run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
        Kokkos::finalize();
        initialized = false;
        return result;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "fastgatk-reblock-gvcf: " << error.what() << '\n';
        return 2;
    }
}
