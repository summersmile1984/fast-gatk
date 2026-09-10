#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string reference;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::string output;
    std::string manifest;
    bool create_index = true;
    // Shared GATK writer option: retain site/INFO columns but omit all
    // sample FORMAT payload from the published VCF.  Keep processing on the
    // full header so allele/PL/AD merge semantics are resolved first.
    bool sites_only_vcf_output = false;
    // GATK CombineGVCs preserves no-call GT semantics.  Native GT/GQ
    // materialization remains an explicit opt-in for kernel experiments.
    bool native_genotype_materialization = false;
    // GATK's public CombineGVCFs switch is --call-genotypes.  Keep the
    // historical native spelling as an alias, but make the standard CLI
    // directly replaceable in Nextflow/SLURM command lines.
    bool call_genotypes = false;
    // Reference blocks can be split at deterministic genomic boundaries.  A
    // value of one is the GATK --convert-to-base-pair-resolution mode.
    bool convert_to_base_pair_resolution = false;
    int break_bands_at_multiples_of = 0;
    // Optional one-record-per-input k-way merge.  Zero keeps the historical
    // aggregate path; --stream-merge enables bounded locus state.
    bool stream_merge = false;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
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

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) &&
           std::filesystem::file_size(path, error) > 0;
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

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-combine-gvcfs (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE              input VCF/GVCF (repeatable)\n"
                         "  -R, --reference FILE            reference FASTA (accepted)\n"
                         "  -L, --intervals REGION          contig:start-end\n"
                         "      --interval-set-rule RULE  UNION (default) or INTERSECTION\n"
                         "  -O, --output FILE               combined VCF/GVCF\n"
                         "      --output-manifest FILE     OutputManifest JSON\n"
                         "      --sites-only-vcf-output [true|false]\n"
                         "                                     omit FORMAT/sample columns\n"
                         "      --call-genotypes [true|false]\n"
                         "                                     emit called GT/GQ fields\n"
                         "      --convert-to-base-pair-resolution [true|false]\n"
                         "                                     split reference bands to one base\n"
                         "      --break-bands-at-multiples-of N\n"
                         "                                     split bands before 1-based multiples of N\n"
                         "      --stream-merge                bounded k-way locus merge\n"
                         "      --fastgatk-materialize-genotypes  native experimental GT/GQ mode\n";
            std::exit(0);
        } else if (is_option(argument, "--variant") || argument == "-V")
            options.inputs.push_back(require_value(index, argc, argv, argument, "--variant", "-V"));
        else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--intervals") || is_option(argument, "--interval") ||
                 is_option(argument, "--region") || argument == "-L") {
            const char* interval_option = argument == "-L" ? "--intervals" :
                argument.rfind("--intervals", 0) == 0 ? "--intervals" :
                argument.rfind("--interval", 0) == 0 ? "--interval" : "--region";
            options.regions.push_back(require_value(index, argc, argv, argument, interval_option, "-L"));
        }
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--sites-only-vcf-output" ||
                   argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--sites-only-vcf-output");
        } else if (argument == "--fastgatk-materialize-genotypes") {
            options.native_genotype_materialization = true;
        } else if (argument == "--call-genotypes" ||
                   argument.rfind("--call-genotypes=", 0) == 0) {
            options.call_genotypes = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--call-genotypes");
            options.native_genotype_materialization =
                options.native_genotype_materialization || options.call_genotypes;
        } else if (argument == "--convert-to-base-pair-resolution" ||
                   argument.rfind("--convert-to-base-pair-resolution=", 0) == 0) {
            options.convert_to_base_pair_resolution = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--convert-to-base-pair-resolution");
        } else if (argument == "--stream-merge" || argument == "--stream-by-locus") {
            options.stream_merge = true;
        } else if (is_option(argument, "--break-bands-at-multiples-of")) {
            const auto value = require_value(index, argc, argv, argument,
                                              "--break-bands-at-multiples-of");
            try { options.break_bands_at_multiples_of = std::stoi(value); }
            catch (...) {
                throw std::invalid_argument("invalid --break-bands-at-multiples-of: " + value);
            }
            if (options.break_bands_at_multiples_of < 0)
                throw std::invalid_argument("--break-bands-at-multiples-of must be non-negative");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility flags are accepted; HTSlib validates actual contigs.
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
    if (options.convert_to_base_pair_resolution)
        options.break_bands_at_multiples_of = 1;
    // Both spellings are intentionally reflected in one execution flag so
    // telemetry and the merge path cannot diverge.
    options.native_genotype_materialization = options.native_genotype_materialization ||
                                               options.call_genotypes;
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct Region {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
    std::string contig;
};

void append_region_selector(const std::string& selector, const bcf_hdr_t* header,
                            std::vector<Region>& regions,
                            fastgatk::io::IntervalFileStats& stats) {
    std::vector<fastgatk::io::IndexedInterval> parsed;
    fastgatk::io::append_interval_selector(selector, header, parsed, stats);
    for (const auto& interval : parsed) {
        const auto* name = interval.rid >= 0 ? bcf_hdr_id2name(header, interval.rid) : nullptr;
        regions.push_back(Region{interval.rid, interval.begin, interval.end,
                                 name == nullptr ? std::string{} : std::string(name)});
    }
}

void normalize_combine_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        const auto left_name = left.contig.empty() ? std::string{} : left.contig;
        const auto right_name = right.contig.empty() ? std::string{} : right.contig;
        if (left_name != right_name) return left_name < right_name;
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.end <= region.begin) continue;
        if (!merged.empty()) {
            auto& previous = merged.back();
            const bool same_contig = !region.contig.empty() || !previous.contig.empty()
                ? region.contig == previous.contig
                : region.rid == previous.rid;
            if (same_contig && region.begin <= previous.end) {
                previous.end = std::max(previous.end, region.end);
                continue;
            }
        }
        merged.push_back(region);
    }
    regions.swap(merged);
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

// Apply GATK's interval-set rule across repeatable -L selectors.  The default
// UNION keeps historical behavior; INTERSECTION computes a half-open set
// intersection between each selector group before record traversal.
void append_combine_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_combine_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    intersection.reserve(std::min(regions.size(), incoming.size()));
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].contig < incoming[right].contig) {
            ++left;
            continue;
        }
        if (incoming[right].contig < regions[left].contig) {
            ++right;
            continue;
        }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end) {
            Region overlap = regions[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (regions[left].end < incoming[right].end) ++left;
        else ++right;
    }
    regions.swap(intersection);
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, static_cast<int>(values[0]));
    free(values);
    return end;
}

bool in_regions(const bcf_hdr_t* header, bcf1_t* record, const std::vector<Region>& regions) {
    if (regions.empty()) return true;
    const auto record_end = record_end_exclusive(header, record);
    // CombineGVCFs' Java traversal deliberately drops a reference-confidence
    // block when an interval ends inside that block (it logs "cuts in the
    // middle of one or more gVCF blocks").  This differs from ordinary
    // VariantContext overlap selection: a concrete variant is retained when
    // it overlaps the interval, while REF/<NON_REF> blocks are retained only
    // when their END is contained by the selected interval.  Keep this rule
    // at the shared input predicate so aggregate and --stream-merge paths
    // cannot diverge at an interval boundary.
    const bool reference_block = record->n_allele == 2 &&
        record->d.allele != nullptr && record->d.allele[1] != nullptr &&
        std::strcmp(record->d.allele[1], "<NON_REF>") == 0;
    const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        const bool same_contig = !region.contig.empty()
            ? (contig != nullptr && region.contig == contig)
            : record->rid == region.rid;
        const bool overlaps = same_contig && record_end > region.begin && record->pos < region.end;
        return overlaps && (!reference_block || record_end <= region.end);
    });
}

struct CombineTraversalIndex {
    hts_idx_t* index = nullptr;
    tbx_t* tabix = nullptr;

    CombineTraversalIndex() = default;

    ~CombineTraversalIndex() {
        if (tabix != nullptr) tbx_destroy(tabix);
        else if (index != nullptr) hts_idx_destroy(index);
    }
    CombineTraversalIndex(const CombineTraversalIndex&) = delete;
    CombineTraversalIndex& operator=(const CombineTraversalIndex&) = delete;

    int tid_for(const bcf_hdr_t* header, const Region& region) const {
        if (tabix != nullptr) {
            const auto* name = !region.contig.empty() ? region.contig.c_str() :
                (region.rid >= 0 ? bcf_hdr_id2name(header, region.rid) : nullptr);
            return name == nullptr ? -1 : tbx_name2id(tabix, name);
        }
        if (!region.contig.empty()) return bcf_hdr_name2id(header, region.contig.c_str());
        return region.rid;
    }
};

std::unique_ptr<CombineTraversalIndex> load_combine_traversal_index(
    const std::string& filename) {
    const auto csi = filename + ".csi";
    const auto tbi = filename + ".tbi";
    if (std::filesystem::exists(csi)) {
        if (auto* index = hts_idx_load3(filename.c_str(), csi.c_str(), HTS_FMT_CSI,
                                         HTS_IDX_SILENT_FAIL)) {
            auto result = std::make_unique<CombineTraversalIndex>();
            result->index = index;
            return result;
        }
    }
    if (std::filesystem::exists(tbi)) {
        if (auto* tabix = tbx_index_load2(filename.c_str(), tbi.c_str())) {
            auto result = std::make_unique<CombineTraversalIndex>();
            result->index = tabix->idx;
            result->tabix = tabix;
            return result;
        }
    }
    return {};
}

struct Record {
    bcf1_t* value = nullptr;
    int rid = -1;
    int pos = -1;
    int end = -1;
    int depth = -1;
    std::string alleles;
    std::vector<std::string> allele_names;
    std::string key;
    struct SampleCall {
        int ploidy = 2;
        std::vector<std::int32_t> genotype{bcf_gt_missing, bcf_gt_missing};
        std::int32_t depth = bcf_int32_missing;
        std::int32_t gq = bcf_int32_missing;
        std::vector<std::int32_t> allele_depths;
        std::vector<std::int32_t> likelihoods;

        bool has_genotype() const {
            return std::any_of(genotype.begin(), genotype.end(), [](const auto encoded) {
                return encoded != bcf_int32_vector_end && !bcf_gt_is_missing(encoded);
            });
        }

        bool operator==(const SampleCall& other) const {
            return ploidy == other.ploidy && genotype == other.genotype &&
                   depth == other.depth && gq == other.gq &&
                   allele_depths == other.allele_depths && likelihoods == other.likelihoods;
        }
    };
    std::vector<SampleCall> samples;
};

struct GenotypeKernelTelemetry {
    std::uint64_t pl_remap_calls = 0;
    double pl_remap_prepare_seconds = 0.0;
    double pl_remap_seconds = 0.0;
    std::string pl_remap_execution_space;
    std::uint64_t allele_field_remap_calls = 0;
    double allele_field_remap_prepare_seconds = 0.0;
    double allele_field_remap_seconds = 0.0;
    std::string allele_field_remap_execution_space;
    std::uint64_t genotype_calls = 0;
    double genotype_prepare_seconds = 0.0;
    double genotype_seconds = 0.0;
    std::string genotype_execution_space;
};

GenotypeKernelTelemetry genotype_kernel_telemetry;

void destroy_records(std::vector<Record>& records) {
    for (auto& record : records) bcf_destroy(record.value);
    records.clear();
}

std::string record_key(const Record& record) {
    return std::to_string(record.rid) + ":" + std::to_string(record.pos) + ":" +
           std::to_string(record.end) + ":" + record.alleles + ":" +
           std::to_string(record.depth);
}

std::vector<std::string> split_alleles(const std::string& text) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const auto end = text.find(',', begin);
        result.push_back(text.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
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

void derive_sample_gt_gq(Record::SampleCall& sample, int allele_count,
                         bool preserve_existing_gq = true) {
    if (sample.likelihoods.empty() || allele_count < 2 || sample.ploidy <= 0) return;
    const auto width = genotype_width(allele_count, sample.ploidy);
    if (width == 0 || sample.likelihoods.size() != width) return;
    const auto derived = fastgatk::kernels::derive_genotype_gt_gq_kokkos(
        sample.likelihoods, 1, allele_count, sample.ploidy);
    // CombineGVCFs --call-genotypes uses PLs to choose GT, but the Java
    // GenotypeBuilder retains an existing GQ annotation when one is present
    // on the source call.  Only synthesize GQ for PL-only records (for
    // example triploid synthetic fixtures); this preserves the standard
    // GATK command's observable FORMAT contract while keeping the native
    // kernel as the single GT/GQ implementation for missing annotations.
    const auto original_gq = sample.gq;
    ++genotype_kernel_telemetry.genotype_calls;
    genotype_kernel_telemetry.genotype_prepare_seconds += derived.prepare_seconds;
    genotype_kernel_telemetry.genotype_seconds += derived.seconds;
    genotype_kernel_telemetry.genotype_execution_space = derived.execution_space;
    sample.genotype.assign(static_cast<std::size_t>(sample.ploidy), bcf_gt_missing);
    for (int copy = 0; copy < sample.ploidy; ++copy) {
        const auto allele = derived.alleles[static_cast<std::size_t>(copy)];
        if (allele < 0) return;
        sample.genotype[static_cast<std::size_t>(copy)] = bcf_gt_unphased(allele);
    }
    if ((!preserve_existing_gq || original_gq == bcf_int32_missing) &&
        !derived.gq.empty() && derived.gq[0] >= 0)
        sample.gq = derived.gq[0];
}

bool is_reference_block(const Record& record) {
    // A variant gVCF record commonly carries `ALT=SNV,<NON_REF>`; that is a
    // variant site, not an uncovered reference block.  Only the exact
    // REF,<NON_REF> shape is a block eligible for END coalescing.
    return record.allele_names.size() == 2 && record.allele_names[1] == "<NON_REF>";
}

Record::SampleCall read_sample_call(const bcf_hdr_t* header, bcf1_t* record,
                                    int sample, int sample_count,
                                    int allele_count) {
    Record::SampleCall call;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto genotype_length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    const int genotype_width = sample_count > 0 && genotype_length > 0
        ? genotype_count / sample_count : 0;
    if (genotype_width > 0 && sample * genotype_width < genotype_count) {
        call.ploidy = genotype_width;
        call.genotype.assign(static_cast<std::size_t>(genotype_width), bcf_gt_missing);
        for (int copy = 0; copy < genotype_width; ++copy) {
            const auto encoded = genotypes[sample * genotype_width + copy];
            if (encoded != bcf_int32_vector_end && !bcf_gt_is_missing(encoded))
                call.genotype[static_cast<std::size_t>(copy)] = encoded;
        }
    }
    free(genotypes);

    int32_t* depths = nullptr;
    int depth_count = 0;
    const auto depth_length = bcf_get_format_int32(header, record, "DP", &depths, &depth_count);
    if (depth_length > 0 && depth_count >= sample_count)
        call.depth = depths[sample];
    free(depths);

    int32_t* genotype_qualities = nullptr;
    int genotype_quality_count = 0;
    const auto gq_length = bcf_get_format_int32(header, record, "GQ",
                                                &genotype_qualities, &genotype_quality_count);
    if (gq_length > 0 && genotype_quality_count >= sample_count)
        call.gq = genotype_qualities[sample];
    free(genotype_qualities);

    int32_t* allele_depths = nullptr;
    int allele_depth_count = 0;
    const auto ad_length = bcf_get_format_int32(header, record, "AD", &allele_depths,
                                                  &allele_depth_count);
    const int ad_width = sample_count > 0 && ad_length > 0
        ? allele_depth_count / sample_count : 0;
    if (ad_width > 0 && sample * ad_width + ad_width <= allele_depth_count)
        call.allele_depths.assign(allele_depths + sample * ad_width,
                                   allele_depths + (sample + 1) * ad_width);
    free(allele_depths);

    int32_t* likelihoods = nullptr;
    int likelihood_count = 0;
    const auto pl_length = bcf_get_format_int32(header, record, "PL", &likelihoods,
                                                 &likelihood_count);
    const int pl_width = sample_count > 0 && pl_length > 0
        ? likelihood_count / sample_count : 0;
    if (pl_width > 0 && sample * pl_width + pl_width <= likelihood_count)
        call.likelihoods.assign(likelihoods + sample * pl_width,
                                likelihoods + (sample + 1) * pl_width);
    free(likelihoods);
    // Some gVCFs carry a diploid GT/GQ but omit PL on reference blocks.  The
    // call remains usable for merging; AF/genotype kernels only consume a
    // complete Number=G vector and therefore leave this sample untouched.
    (void)allele_count;
    return call;
}

bool has_sample_data(const Record::SampleCall& call) {
    return call.has_genotype() || call.depth != bcf_int32_missing ||
           call.gq != bcf_int32_missing ||
           !call.allele_depths.empty() || !call.likelihoods.empty();
}

void merge_sample_call(Record::SampleCall& destination,
                       const Record::SampleCall& source) {
    if (!has_sample_data(source)) return;
    if (!has_sample_data(destination)) {
        destination = source;
        return;
    }
    if (!(destination == source))
        throw std::runtime_error(
            "OUTPUT_CONTRACT_FAILURE: conflicting FORMAT values for the same sample/site");
}

void materialize_samples(const bcf_hdr_t* header, Record& record,
                         std::size_t sample_count) {
    if (sample_count == 0) return;
    const auto allele_count = static_cast<std::size_t>(record.value->n_allele);
    int output_ploidy = 2;
    for (std::size_t sample = 0; sample < std::min(sample_count, record.samples.size()); ++sample)
        output_ploidy = std::max(output_ploidy, record.samples[sample].ploidy);
    const auto pl_count = genotype_width(static_cast<int>(allele_count), output_ploidy);
    if (pl_count == 0) throw std::runtime_error("BAD_INPUT: unsupported CombineGVCFs genotype dimensions");
    std::vector<std::int32_t> genotypes(sample_count * static_cast<std::size_t>(output_ploidy), bcf_gt_missing);
    std::vector<std::int32_t> depths(sample_count, bcf_int32_missing);
    std::vector<std::int32_t> genotype_qualities(sample_count, bcf_int32_missing);
    std::vector<std::int32_t> allele_depths(sample_count * allele_count, bcf_int32_missing);
    std::vector<std::int32_t> likelihoods(sample_count * pl_count, bcf_int32_missing);
    bool has_likelihoods = false;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        if (sample >= record.samples.size()) continue;
        const auto& call = record.samples[sample];
        for (int copy = 0; copy < call.ploidy && copy < output_ploidy; ++copy)
            if (static_cast<std::size_t>(copy) < call.genotype.size())
                genotypes[sample * static_cast<std::size_t>(output_ploidy) +
                          static_cast<std::size_t>(copy)] = call.genotype[static_cast<std::size_t>(copy)];
        depths[sample] = call.depth;
        genotype_qualities[sample] = call.gq;
        for (std::size_t i = 0; i < std::min(allele_count, call.allele_depths.size()); ++i)
            allele_depths[sample * allele_count + i] = call.allele_depths[i];
        for (std::size_t i = 0; i < std::min(pl_count, call.likelihoods.size()); ++i)
            likelihoods[sample * pl_count + i] = call.likelihoods[i];
        has_likelihoods = has_likelihoods || !call.likelihoods.empty();
    }
    // Insert fields in GATK's canonical FORMAT order.  HTSlib preserves the
    // record-local insertion order (it does not sort by header ID), so PL
    // must be updated after GT/AD/DP/GQ even though the PL-less branch is
    // decided before publication.
    if (bcf_update_genotypes(header, record.value, genotypes.data(),
                             static_cast<int>(genotypes.size())) != 0 ||
        bcf_update_format_int32(header, record.value, "AD", allele_depths.data(),
                                static_cast<int>(allele_depths.size())) != 0 ||
        bcf_update_format_int32(header, record.value, "DP", depths.data(),
                                static_cast<int>(depths.size())) != 0 ||
        bcf_update_format_int32(header, record.value, "GQ", genotype_qualities.data(),
                                static_cast<int>(genotype_qualities.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize merged FORMAT fields");
    const bool has_pl_header = bcf_hdr_id2int(header, BCF_DT_ID, "PL") >= 0;
    if (has_pl_header &&
        (has_likelihoods
            ? bcf_update_format_int32(header, record.value, "PL", likelihoods.data(),
                                      static_cast<int>(likelihoods.size()))
            : bcf_update_format_int32(header, record.value, "PL", nullptr, 0)) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize merged PL field");
}

void merge_record_samples(Record& destination, const Record& source) {
    if (destination.samples.size() < source.samples.size())
        destination.samples.resize(source.samples.size());
    for (std::size_t sample = 0; sample < source.samples.size(); ++sample)
        merge_sample_call(destination.samples[sample], source.samples[sample]);
}

void remap_record_to_allele_union(const bcf_hdr_t* output_header, Record& record,
                                  const std::vector<std::string>& union_alleles,
                                  bool native_genotype_materialization,
                                  bool gatk_call_genotypes = false) {
    const bool preserve_gq = !is_reference_block(record);
    if (record.allele_names == union_alleles) {
        if (native_genotype_materialization) {
            for (auto& sample : record.samples)
                derive_sample_gt_gq(sample, static_cast<int>(union_alleles.size()), preserve_gq);
        } else {
            for (auto& sample : record.samples)
                sample.genotype.assign(static_cast<std::size_t>(sample.ploidy), bcf_gt_missing);
        }
        return;
    }
    if (record.allele_names.empty())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: missing CombineGVCFs allele metadata");
    std::vector<int> old_to_new(record.allele_names.size(), -1);
    for (std::size_t old = 0; old < record.allele_names.size(); ++old) {
        const auto found = std::find(union_alleles.begin(), union_alleles.end(), record.allele_names[old]);
        if (found == union_alleles.end())
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: CombineGVCFs ALT missing from union");
        old_to_new[old] = static_cast<int>(found - union_alleles.begin());
    }
    const int old_alleles = static_cast<int>(record.allele_names.size());
    const int new_alleles = static_cast<int>(union_alleles.size());
    std::vector<std::int32_t> target_to_old(static_cast<std::size_t>(new_alleles), -1);
    for (int target = 0; target < new_alleles; ++target) {
        const auto found = std::find(record.allele_names.begin(), record.allele_names.end(),
                                     union_alleles[static_cast<std::size_t>(target)]);
        if (found != record.allele_names.end())
            target_to_old[static_cast<std::size_t>(target)] =
                static_cast<std::int32_t>(found - record.allele_names.begin());
    }
    for (auto& sample : record.samples) {
        for (auto& encoded : sample.genotype) {
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            const int old_index = bcf_gt_allele(encoded);
            if (old_index < 0 || old_index >= old_alleles || old_to_new[old_index] < 0) {
                encoded = bcf_gt_missing;
                continue;
            }
            const int new_index = old_to_new[old_index];
            encoded = bcf_gt_is_phased(encoded) ? bcf_gt_phased(new_index)
                                                 : bcf_gt_unphased(new_index);
        }
        if (!sample.allele_depths.empty()) {
            if (sample.allele_depths.size() == static_cast<std::size_t>(old_alleles)) {
                const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
                    sample.allele_depths, 1, old_alleles, new_alleles, target_to_old);
                ++genotype_kernel_telemetry.allele_field_remap_calls;
                genotype_kernel_telemetry.allele_field_remap_prepare_seconds += remapped.prepare_seconds;
                genotype_kernel_telemetry.allele_field_remap_seconds += remapped.seconds;
                genotype_kernel_telemetry.allele_field_remap_execution_space = remapped.execution_space;
                sample.allele_depths = remapped.values;
                if (!native_genotype_materialization || gatk_call_genotypes) {
                    const auto non_ref = std::find(record.allele_names.begin(),
                                                   record.allele_names.end(), "<NON_REF>");
                    const auto non_ref_index = non_ref == record.allele_names.end()
                        ? -1 : static_cast<int>(non_ref - record.allele_names.begin());
                    const auto non_ref_depth = non_ref_index >= 0 &&
                            non_ref_index < static_cast<int>(sample.allele_depths.size())
                        ? sample.allele_depths[static_cast<std::size_t>(old_to_new[non_ref_index])]
                        : 0;
                    for (std::size_t index = 0; index < sample.allele_depths.size(); ++index)
                        if (target_to_old[index] < 0 &&
                            sample.allele_depths[index] == bcf_int32_missing)
                            sample.allele_depths[index] = non_ref_depth;
                }
            } else {
                std::vector<std::int32_t> remapped(static_cast<std::size_t>(new_alleles), bcf_int32_missing);
                for (int old_index = 0; old_index < old_alleles && old_index < static_cast<int>(sample.allele_depths.size()); ++old_index)
                    if (old_to_new[old_index] >= 0) remapped[old_to_new[old_index]] = sample.allele_depths[old_index];
                if (!native_genotype_materialization || gatk_call_genotypes) {
                    const auto non_ref = std::find(record.allele_names.begin(),
                                                   record.allele_names.end(), "<NON_REF>");
                    const auto non_ref_index = non_ref == record.allele_names.end()
                        ? -1 : static_cast<int>(non_ref - record.allele_names.begin());
                    const auto non_ref_depth = non_ref_index >= 0 &&
                            non_ref_index < static_cast<int>(sample.allele_depths.size())
                        ? sample.allele_depths[static_cast<std::size_t>(non_ref_index)]
                        : 0;
                    for (int target = 0; target < new_alleles; ++target)
                        if (target_to_old[static_cast<std::size_t>(target)] < 0 &&
                            remapped[static_cast<std::size_t>(target)] == bcf_int32_missing)
                            remapped[static_cast<std::size_t>(target)] = non_ref_depth;
                }
                sample.allele_depths = std::move(remapped);
            }
        }
        if (!sample.likelihoods.empty()) {
            const auto old_width = genotype_width(old_alleles, sample.ploidy);
            const auto new_width = genotype_width(new_alleles, sample.ploidy);
            if (old_width > 0 && new_width > 0 && sample.likelihoods.size() == old_width) {
                const auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
                    sample.likelihoods, 1, old_alleles, new_alleles, sample.ploidy,
                    target_to_old);
                ++genotype_kernel_telemetry.pl_remap_calls;
                genotype_kernel_telemetry.pl_remap_prepare_seconds += remapped.prepare_seconds;
                genotype_kernel_telemetry.pl_remap_seconds += remapped.seconds;
                genotype_kernel_telemetry.pl_remap_execution_space = remapped.execution_space;
                sample.likelihoods = remapped.pl;
                if (!native_genotype_materialization || gatk_call_genotypes) {
                    for (auto& value : sample.likelihoods)
                        if (value == bcf_int32_missing) value = 99;
                }
            } else {
                std::vector<std::int32_t> remapped(new_width, bcf_int32_missing);
                if (!native_genotype_materialization || gatk_call_genotypes)
                    std::fill(remapped.begin(), remapped.end(), 99);
                sample.likelihoods = std::move(remapped);
            }
        }
        if (native_genotype_materialization)
            derive_sample_gt_gq(sample, new_alleles, preserve_gq);
        else
            sample.genotype.assign(static_cast<std::size_t>(sample.ploidy), bcf_gt_missing);
    }
    record.allele_names = union_alleles;
    record.alleles.clear();
    for (std::size_t index = 0; index < union_alleles.size(); ++index) {
        if (index != 0) record.alleles.push_back(',');
        record.alleles += union_alleles[index];
    }
    if (bcf_update_alleles_str(output_header, record.value, record.alleles.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update CombineGVCFs ALT union");
}

int info_int(const bcf_hdr_t* header, bcf1_t* record, const char* name, int fallback) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, name, &values, &count);
    const auto result = length > 0 && count > 0 ? values[0] : fallback;
    free(values);
    return result;
}

bcf_hdr_t* make_output_header(const bcf_hdr_t* source) {
    auto* header = bcf_hdr_init("w");
    if (!header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot create output VCF header");
    bcf_hdr_append(header, "##fileformat=VCFv4.2");
    bcf_hdr_append(header, "##source=fastgatk-combine-gvcfs");
    bcf_hdr_append(header, "##fastgatk_combine_gvcfs_status=prototype-site-level-merge");
    for (int rid = 0; rid < source->n[BCF_DT_CTG]; ++rid) {
        const char* name = bcf_hdr_id2name(source, rid);
        // Structured contig hrec values are not complete textual header lines
        // (their value pointer is null in HTSlib); reconstruct a valid line.
        bcf_hdr_printf(header, "##contig=<ID=%s>", name);
    }
    bcf_hdr_append(header, "##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>");
    bcf_hdr_append(header, "##INFO=<ID=END,Number=1,Type=Integer,Description=End position of reference block>");
    bcf_hdr_append(header, "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>");
    bcf_hdr_append(header, "##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>");
    bcf_hdr_append(header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>");
    bcf_hdr_append(header, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>");
    bcf_hdr_append(header, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>");
    bcf_hdr_append(header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>");
    bcf_hdr_append(header, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled genotype likelihoods>");
    bcf_hdr_append(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    if (bcf_hdr_sync(header) != 0) {
        bcf_hdr_destroy(header);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync output VCF header");
    }
    return header;
}

std::uint64_t combine_record_bytes(const Record& record) {
    std::uint64_t bytes = sizeof(Record) + record.alleles.size() + record.key.size();
    for (const auto& allele : record.allele_names) bytes += allele.size();
    for (const auto& sample : record.samples) {
        bytes += sizeof(Record::SampleCall) +
                 sample.genotype.size() * sizeof(std::int32_t) +
                 sample.allele_depths.size() * sizeof(std::int32_t) +
                 sample.likelihoods.size() * sizeof(std::int32_t);
    }
    // bcf1_t owns HTSlib-managed strings/INFO/FORMAT buffers which are not
    // exposed as portable fields.  Keep a conservative fixed allowance so
    // the merge queue never treats a record as zero-cost.
    return std::max<std::uint64_t>(1, bytes + 1024);
}

// Materialize only the source record envelope and sample-major FORMAT state.
// The aggregate path may turn this into a vector of reference-band segments;
// stream mode keeps this one record as the cursor's lazy split state instead.
Record materialize_combine_record(const bcf_hdr_t* input_header,
                                  bcf_hdr_t* output_header,
                                  bcf1_t* input_record,
                                  const int output_rid,
                                  const std::size_t sample_base,
                                  const std::size_t sample_count) {
    bcf_unpack(input_record, BCF_UN_ALL);
    std::string allele_string;
    for (int allele = 0; allele < input_record->n_allele; ++allele) {
        if (allele != 0) allele_string.push_back(',');
        allele_string += input_record->d.allele[allele];
    }
    const auto end = info_int(input_header, input_record, "END", input_record->pos + 1);
    const auto depth = info_int(input_header, input_record, "DP", -1);
    auto* output_record = bcf_init();
    if (!output_record || bcf_update_alleles_str(output_header, output_record,
                                                 allele_string.c_str()) != 0) {
        bcf_destroy(output_record);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize combined record");
    }
    output_record->rid = output_rid;
    output_record->pos = input_record->pos;
    bcf_float_set_missing(output_record->qual);
    if (depth >= 0) bcf_update_info_int32(output_header, output_record, "DP", &depth, 1);
    if (end > input_record->pos + 1)
        bcf_update_info_int32(output_header, output_record, "END", &end, 1);
    bcf_unpack(output_record, BCF_UN_STR);
    Record item;
    item.value = output_record;
    item.rid = output_rid;
    item.pos = input_record->pos;
    item.end = end;
    item.depth = depth;
    item.alleles = allele_string;
    item.allele_names = split_alleles(allele_string);
    item.samples.resize(sample_count);
    for (int sample = 0; sample < input_header->n[BCF_DT_SAMPLE]; ++sample)
        item.samples[sample_base + static_cast<std::size_t>(sample)] =
            read_sample_call(input_header, input_record, sample,
                             input_header->n[BCF_DT_SAMPLE], input_record->n_allele);
    if (item.depth < 0 && !is_reference_block(item)) {
        long long depth_sum = 0;
        bool has_depth = false;
        for (const auto& call : item.samples) {
            if (call.depth == bcf_int32_missing || call.depth < 0) continue;
            depth_sum += call.depth;
            has_depth = true;
        }
        if (has_depth && depth_sum <= std::numeric_limits<int>::max())
            item.depth = static_cast<int>(depth_sum);
        if (item.depth >= 0 && bcf_update_info_int32(output_header, output_record,
                                                      "DP", &item.depth, 1) != 0) {
            bcf_destroy(output_record);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize DP annotation");
        }
    }
    item.key = record_key(item);
    return item;
}

[[maybe_unused]] std::vector<Record> materialize_combine_segments(
    const bcf_hdr_t* input_header,
    bcf_hdr_t* output_header,
    bcf1_t* input_record,
    const int output_rid,
    const std::size_t sample_base,
    const std::size_t sample_count,
    faidx_t* reference_index,
    const Options& options) {
    std::vector<Record> segments;
    bcf_unpack(input_record, BCF_UN_ALL);
    std::string allele_string;
    for (int allele = 0; allele < input_record->n_allele; ++allele) {
        if (allele != 0) allele_string.push_back(',');
        allele_string += input_record->d.allele[allele];
    }
    const auto end = info_int(input_header, input_record, "END", input_record->pos + 1);
    const auto depth = info_int(input_header, input_record, "DP", -1);
    auto* output_record = bcf_init();
    if (!output_record || bcf_update_alleles_str(output_header, output_record,
                                                 allele_string.c_str()) != 0) {
        bcf_destroy(output_record);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize combined record");
    }
    output_record->rid = output_rid;
    output_record->pos = input_record->pos;
    bcf_float_set_missing(output_record->qual);
    if (depth >= 0) bcf_update_info_int32(output_header, output_record, "DP", &depth, 1);
    if (end > input_record->pos + 1)
        bcf_update_info_int32(output_header, output_record, "END", &end, 1);
    bcf_unpack(output_record, BCF_UN_STR);
    Record item;
    item.value = output_record;
    item.rid = output_rid;
    item.pos = input_record->pos;
    item.end = end;
    item.depth = depth;
    item.alleles = allele_string;
    item.allele_names = split_alleles(allele_string);
    item.samples.resize(sample_count);
    for (int sample = 0; sample < input_header->n[BCF_DT_SAMPLE]; ++sample)
        item.samples[sample_base + static_cast<std::size_t>(sample)] =
            read_sample_call(input_header, input_record, sample,
                             input_header->n[BCF_DT_SAMPLE], input_record->n_allele);
    if (item.depth < 0) {
        long long depth_sum = 0;
        bool has_depth = false;
        for (const auto& call : item.samples) {
            if (call.depth == bcf_int32_missing || call.depth < 0) continue;
            depth_sum += call.depth;
            has_depth = true;
        }
        if (has_depth && depth_sum <= std::numeric_limits<int>::max())
            item.depth = static_cast<int>(depth_sum);
        if (item.depth >= 0 && bcf_update_info_int32(output_header, output_record,
                                                      "DP", &item.depth, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize DP annotation");
    }
    const auto break_multiple = options.convert_to_base_pair_resolution
        ? 1 : options.break_bands_at_multiples_of;
    const bool split_block = is_reference_block(item) && item.end > item.pos + 1 &&
                             break_multiple > 0;
    std::vector<int> segment_starts{item.pos};
    std::vector<int> segment_ends{item.end};
    if (split_block) {
        segment_starts.clear();
        segment_ends.clear();
        int segment_start = item.pos;
        for (long long coordinate =
                 (static_cast<long long>(item.pos) + 2LL + break_multiple - 1LL) /
                 break_multiple * break_multiple;
             coordinate <= static_cast<long long>(item.end);
             coordinate += break_multiple) {
            const auto boundary = static_cast<int>(coordinate - 1LL);
            if (boundary <= segment_start || boundary >= item.end) continue;
            segment_starts.push_back(segment_start);
            segment_ends.push_back(boundary);
            segment_start = boundary;
        }
        segment_starts.push_back(segment_start);
        segment_ends.push_back(item.end);
    }
    const auto item_rid = item.rid;
    const auto item_samples = item.samples;
    const auto item_allele_names = item.allele_names;
    bcf1_t* original_value = item.value;
    item.value = nullptr;
    try {
        for (std::size_t segment = 0; segment < segment_starts.size(); ++segment) {
            Record output_segment;
            output_segment.value = segment == 0 ? original_value : bcf_dup(original_value);
            if (!output_segment.value)
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot split reference block");
            output_segment.rid = item_rid;
            output_segment.alleles = allele_string;
            output_segment.allele_names = item_allele_names;
            output_segment.samples = item_samples;
            output_segment.depth = item.depth;
            output_segment.pos = segment_starts[segment];
            output_segment.end = segment_ends[segment];
            output_segment.value->pos = output_segment.pos;
            std::string reference_allele = output_segment.allele_names.front();
            if (reference_index && split_block) {
                const char* contig_name = bcf_hdr_id2name(output_header, output_segment.rid);
                int fetched_length = 0;
                char* fetched = contig_name == nullptr ? nullptr :
                    faidx_fetch_seq(reference_index, contig_name, output_segment.pos,
                                    output_segment.pos, &fetched_length);
                if (fetched != nullptr && fetched_length == 1)
                    reference_allele.assign(1, fetched[0]);
                free(fetched);
            }
            if (split_block && !reference_allele.empty()) {
                output_segment.allele_names.front() = reference_allele;
                output_segment.alleles = reference_allele + ",<NON_REF>";
                if (bcf_update_alleles_str(output_header, output_segment.value,
                                           output_segment.alleles.c_str()) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split reference allele");
            }
            if (output_segment.end > output_segment.pos + 1) {
                const auto segment_end = output_segment.end;
                if (bcf_update_info_int32(output_header, output_segment.value,
                                          "END", &segment_end, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split END");
            } else if (bcf_update_info_int32(output_header, output_segment.value,
                                              "END", nullptr, 0) != 0) {
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear split END");
            }
            output_segment.key = record_key(output_segment);
            segments.push_back(std::move(output_segment));
        }
    } catch (...) {
        if (original_value != nullptr &&
            std::none_of(segments.begin(), segments.end(),
                         [&](const Record& segment) { return segment.value == original_value; }))
            bcf_destroy(original_value);
        destroy_records(segments);
        throw;
    }
    return segments;
}

struct CombineStreamSummary {
    std::uint64_t input_records = 0;
    std::uint64_t skipped_records = 0;
    std::uint64_t interval_skipped = 0;
    std::uint64_t output_records = 0;
    std::uint64_t peak_queued_bytes = 0;
    std::uint64_t max_inflight_records = 0;
    std::uint64_t lazy_reference_blocks = 0;
    std::uint64_t lazy_reference_segments = 0;
    std::uint64_t indexed_inputs = 0;
    std::uint64_t indexed_interval_queries = 0;
};

int run_streaming_tool(const Options& options,
                       const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto started = std::chrono::steady_clock::now();
    bcf_hdr_t* output_header = nullptr;
    bcf_hdr_t* writer_header = nullptr;
    faidx_t* reference_index = nullptr;
    fastgatk::io::IntervalFileStats interval_file_stats;
    std::vector<Region> regions;
    bool intervals_parsed = false;
    std::vector<std::string> sample_names;
    std::set<std::string> sample_name_set;
    struct Cursor {
        std::string path;
        htsFile* file = nullptr;
        bcf_hdr_t* header = nullptr;
        bcf1_t* input_record = nullptr;
        std::size_t sample_base = 0;
        std::deque<Record> pending;
        std::unique_ptr<CombineTraversalIndex> traversal_index;
        hts_itr_t* iterator = nullptr;
        kstring_t indexed_line{0, 0, nullptr};
        std::size_t region_index = 0;
        // Reference-band splitting is lazy in stream mode.  Keep the source
        // record and next coordinate only; do not allocate one Record per
        // base (or per split segment) before the merge heap consumes it.
        std::optional<Record> split_block;
        int split_cursor = -1;
        bool eof = false;
        int last_rid = -1;
        int last_pos = -1;
    };
    std::vector<Cursor> cursors;
    CombineStreamSummary summary;
        const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    std::uint64_t queued_bytes = 0;
    htsFile* output = nullptr;
    std::string index_path;
    try {
        if (!options.reference.empty()) {
            reference_index = fai_load(options.reference.c_str());
            if (!reference_index)
                throw std::runtime_error("BAD_INPUT: cannot load reference FASTA index: " + options.reference);
        }
        cursors.reserve(options.inputs.size());
        for (const auto& input_path : options.inputs) {
            Cursor cursor;
            cursor.path = input_path;
            cursor.file = bcf_open(input_path.c_str(), "r");
            if (!cursor.file) throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            cursor.header = bcf_hdr_read(cursor.file);
            if (!cursor.header)
                throw std::runtime_error("BAD_INPUT: cannot read VCF header: " + input_path);
            if (!output_header) output_header = make_output_header(cursor.header);
            cursor.sample_base = sample_names.size();
            for (int sample = 0; sample < cursor.header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(cursor.header, BCF_DT_SAMPLE, sample);
                if (name == nullptr || !sample_name_set.insert(name).second)
                    throw std::runtime_error(
                        "BAD_INPUT: duplicate or invalid sample name across CombineGVCFs inputs");
                sample_names.emplace_back(name);
                if (bcf_hdr_add_sample(output_header, name) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot add sample to output header");
            }
            if (!intervals_parsed) {
                regions.reserve(options.regions.size());
                bool first_selector = true;
                for (const auto& text : options.regions) {
                    append_combine_region_selector_with_rule(
                        text, cursor.header, regions, interval_file_stats,
                        options.interval_set_rule, first_selector);
                    first_selector = false;
                }
                intervals_parsed = true;
            }
            cursor.input_record = bcf_init();
            if (!cursor.input_record)
                throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
            cursors.push_back(std::move(cursor));
        }
        if (!regions.empty()) normalize_combine_regions(regions);
        if (!sample_names.empty() && bcf_hdr_add_sample(output_header, nullptr) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot terminate sample header");
        if (bcf_hdr_sync(output_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync merged sample header");
        writer_header = bcf_hdr_dup(output_header);
        if (!writer_header)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate CombineGVCFs writer header");
        if (options.sites_only_vcf_output && bcf_hdr_set_samples(writer_header, nullptr, 0) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create sites-only writer header");
        if (!regions.empty()) {
            for (auto& cursor : cursors) {
                cursor.traversal_index = load_combine_traversal_index(cursor.path);
                if (cursor.traversal_index != nullptr) ++summary.indexed_inputs;
            }
        }

        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open output VCF/GVCF: " + options.output);
        if (bcf_hdr_write(output, writer_header) != 0)
            throw std::runtime_error("cannot write output VCF/GVCF header");

        const auto read_next_raw = [&](Cursor& cursor) -> bool {
            if (cursor.traversal_index == nullptr)
                return bcf_read(cursor.file, cursor.header, cursor.input_record) == 0;
            while (true) {
                if (cursor.iterator != nullptr) {
                    int status = -1;
                    if (cursor.traversal_index->tabix != nullptr) {
                        status = tbx_itr_next(cursor.file, cursor.traversal_index->tabix,
                                              cursor.iterator, &cursor.indexed_line);
                        if (status >= 0) {
                            if (vcf_parse(&cursor.indexed_line, cursor.header,
                                          cursor.input_record) < 0)
                                throw std::runtime_error(
                                    "BAD_INPUT: indexed CombineGVCFs VCF record parse failed");
                            return true;
                        }
                    } else {
                        status = bcf_itr_next(cursor.file, cursor.iterator,
                                              cursor.input_record);
                        if (status >= 0) return true;
                    }
                    hts_itr_destroy(cursor.iterator);
                    cursor.iterator = nullptr;
                    if (status < -1)
                        throw std::runtime_error(
                            "BAD_INPUT: indexed CombineGVCFs traversal failed");
                }
                if (cursor.region_index >= regions.size()) return false;
                const auto& region = regions[cursor.region_index++];
                ++summary.indexed_interval_queries;
                const auto tid = cursor.traversal_index->tid_for(cursor.header, region);
                if (tid < 0) continue;
                if (cursor.traversal_index->tabix != nullptr)
                    cursor.iterator = tbx_itr_queryi(cursor.traversal_index->tabix, tid,
                                                     region.begin, region.end);
                else
                    cursor.iterator = bcf_itr_queryi(cursor.traversal_index->index, tid,
                                                     region.begin, region.end);
            }
        };

        const auto next_lazy_reference_segment = [&](Cursor& cursor)
            -> std::optional<Record> {
            if (!cursor.split_block.has_value()) return std::nullopt;
            auto& original = *cursor.split_block;
            if (cursor.split_cursor < original.pos ||
                cursor.split_cursor >= original.end) {
                bcf_destroy(original.value);
                original.value = nullptr;
                cursor.split_block.reset();
                cursor.split_cursor = -1;
                return std::nullopt;
            }
            const auto break_multiple = options.convert_to_base_pair_resolution
                ? 1 : options.break_bands_at_multiples_of;
            const int segment_begin = cursor.split_cursor;
            int segment_end = original.end;
            if (break_multiple > 0) {
                const auto coordinate =
                    (static_cast<long long>(segment_begin) + 2LL + break_multiple - 1LL) /
                    break_multiple * break_multiple;
                const auto boundary = static_cast<int>(coordinate - 1LL);
                if (boundary > segment_begin && boundary < original.end)
                    segment_end = boundary;
            }
            auto* copy = bcf_dup(original.value);
            if (!copy)
                throw std::runtime_error(
                    "RESOURCE_EXHAUSTED: cannot clone lazy CombineGVCFs reference segment");
            Record segment;
            segment.value = copy;
            segment.rid = original.rid;
            segment.pos = segment_begin;
            segment.end = segment_end;
            segment.depth = original.depth;
            segment.allele_names = original.allele_names;
            segment.alleles = original.alleles;
            segment.samples = original.samples;
            segment.value->pos = segment_begin;
            const auto reference_block = is_reference_block(original);
            if (reference_block) {
                std::string reference_allele = segment.allele_names.front();
                if (reference_index) {
                    const char* contig_name = bcf_hdr_id2name(output_header, segment.rid);
                    int fetched_length = 0;
                    char* fetched = contig_name == nullptr ? nullptr :
                        faidx_fetch_seq(reference_index, contig_name, segment.pos,
                                        segment.pos, &fetched_length);
                    if (fetched != nullptr && fetched_length == 1)
                        reference_allele.assign(1, fetched[0]);
                    free(fetched);
                }
                segment.allele_names.front() = reference_allele;
                segment.alleles = reference_allele + ",<NON_REF>";
                if (bcf_update_alleles_str(output_header, segment.value,
                                           segment.alleles.c_str()) != 0) {
                    bcf_destroy(segment.value);
                    throw std::runtime_error(
                        "OUTPUT_CONTRACT_FAILURE: cannot update lazy reference allele");
                }
            }
            if (segment.end > segment.pos + 1) {
                const auto end = segment.end;
                if (bcf_update_info_int32(output_header, segment.value, "END",
                                          &end, 1) != 0) {
                    bcf_destroy(segment.value);
                    throw std::runtime_error(
                        "OUTPUT_CONTRACT_FAILURE: cannot update lazy reference END");
                }
            } else if (bcf_update_info_int32(output_header, segment.value, "END",
                                             nullptr, 0) != 0) {
                bcf_destroy(segment.value);
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: cannot clear lazy reference END");
            }
            segment.key = record_key(segment);
            cursor.split_cursor = segment_end;
            ++summary.lazy_reference_segments;
            if (cursor.split_cursor >= original.end) {
                bcf_destroy(original.value);
                original.value = nullptr;
                cursor.split_block.reset();
                cursor.split_cursor = -1;
            }
            return segment;
        };

        const auto load_next = [&](Cursor& cursor) -> bool {
            while (cursor.pending.empty() && !cursor.eof) {
                if (cursor.split_block.has_value()) {
                    if (auto segment = next_lazy_reference_segment(cursor)) {
                        cursor.pending.push_back(std::move(*segment));
                        continue;
                    }
                }
                if (!read_next_raw(cursor)) {
                    cursor.eof = true;
                    break;
                }
                ++summary.input_records;
                bcf_unpack(cursor.input_record, BCF_UN_ALL);
                if (!in_regions(cursor.header, cursor.input_record, regions)) {
                    ++summary.interval_skipped;
                    continue;
                }
                if (cursor.input_record->rid < 0 || cursor.input_record->n_allele < 2) {
                    ++summary.skipped_records;
                    continue;
                }
                if (cursor.last_rid > cursor.input_record->rid ||
                    (cursor.last_rid == cursor.input_record->rid &&
                     cursor.last_pos > cursor.input_record->pos))
                    throw std::runtime_error(
                        "BAD_INPUT: --stream-merge requires coordinate-sorted inputs");
                cursor.last_rid = cursor.input_record->rid;
                cursor.last_pos = cursor.input_record->pos;
                const char* contig = bcf_hdr_id2name(cursor.header, cursor.input_record->rid);
                const auto output_rid = contig == nullptr
                    ? -1 : bcf_hdr_name2id(output_header, contig);
                if (output_rid < 0) {
                    ++summary.skipped_records;
                    continue;
                }
                auto materialized = materialize_combine_record(
                    cursor.header, output_header, cursor.input_record, output_rid,
                    cursor.sample_base, sample_names.size());
                const auto break_multiple = options.convert_to_base_pair_resolution
                    ? 1 : options.break_bands_at_multiples_of;
                if (is_reference_block(materialized) && materialized.end > materialized.pos + 1 &&
                    break_multiple > 0) {
                    cursor.split_cursor = materialized.pos;
                    cursor.split_block.emplace(std::move(materialized));
                    ++summary.lazy_reference_blocks;
                    continue;
                }
                cursor.pending.push_back(std::move(materialized));
            }
            return !cursor.pending.empty();
        };

        struct QueueItem {
            Record record;
            std::size_t source = 0;
        };
        const auto less = [](const QueueItem& left, const QueueItem& right) {
            if (left.record.rid != right.record.rid) return left.record.rid > right.record.rid;
            if (left.record.pos != right.record.pos) return left.record.pos > right.record.pos;
            if (left.record.end != right.record.end) return left.record.end > right.record.end;
            if (left.record.key != right.record.key) return left.record.key > right.record.key;
            return left.source > right.source;
        };
        std::priority_queue<QueueItem, std::vector<QueueItem>, decltype(less)> queue(less);
        for (std::size_t source = 0; source < cursors.size(); ++source) {
            if (load_next(cursors[source])) {
                auto record = std::move(cursors[source].pending.front());
                cursors[source].pending.pop_front();
                queued_bytes += combine_record_bytes(record);
                queue.push(QueueItem{std::move(record), source});
            }
        }
        summary.max_inflight_records = queue.size();
        summary.peak_queued_bytes = queued_bytes;
        if (safe_memory_budget > 0 && queued_bytes > safe_memory_budget)
            throw std::runtime_error("RESOURCE_EXHAUSTED: stream merge queue exceeds safe memory budget");

        std::optional<Record> pending_output;
        const auto write_record = [&](Record& record) {
            remap_record_to_allele_union(output_header, record, record.allele_names,
                                         options.native_genotype_materialization,
                                         options.call_genotypes);
            materialize_samples(output_header, record, sample_names.size());
            // GATK's sites-only CombineGVCFs writer does not promote a
            // sample-derived DP annotation onto a pure reference-confidence
            // block; END remains the only block INFO field.  Apply this
            // lexical cleanup only after merge/filter computation.  In the
            // normal writer an explicit input INFO/DP remains authoritative.
            if (options.sites_only_vcf_output && is_reference_block(record))
                bcf_update_info_int32(output_header, record.value, "DP", nullptr, 0);
            if (options.sites_only_vcf_output && bcf_subset_format(writer_header, record.value) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset sites-only FORMAT payload");
            if (bcf_write(output, writer_header, record.value) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write output VCF/GVCF record");
            ++summary.output_records;
            bcf_destroy(record.value);
            record.value = nullptr;
        };
        while (!queue.empty()) {
            auto item = std::move(const_cast<QueueItem&>(queue.top()));
            queue.pop();
            queued_bytes -= std::min<std::uint64_t>(queued_bytes, combine_record_bytes(item.record));
            if (load_next(cursors[item.source])) {
                auto next = std::move(cursors[item.source].pending.front());
                cursors[item.source].pending.pop_front();
                queued_bytes += combine_record_bytes(next);
                queue.push(QueueItem{std::move(next), item.source});
            }
            summary.max_inflight_records = std::max<std::uint64_t>(
                summary.max_inflight_records, queue.size() + 1);
            summary.peak_queued_bytes = std::max(summary.peak_queued_bytes, queued_bytes);
            if (safe_memory_budget > 0 && queued_bytes > safe_memory_budget)
                throw std::runtime_error("RESOURCE_EXHAUSTED: stream merge queue exceeds safe memory budget");
            remap_record_to_allele_union(output_header, item.record, item.record.allele_names,
                                         options.native_genotype_materialization,
                                         options.call_genotypes);
            if (!pending_output.has_value()) {
                pending_output = std::move(item.record);
                continue;
            }
            auto& previous = *pending_output;
            const auto& current = item.record;
            const bool same_site = !is_reference_block(previous) && !is_reference_block(current) &&
                previous.rid == current.rid && previous.pos == current.pos &&
                !previous.allele_names.empty() && !current.allele_names.empty() &&
                previous.allele_names.front() == current.allele_names.front();
            const bool same_block = is_reference_block(previous) && is_reference_block(current) &&
                previous.rid == current.rid && previous.alleles == current.alleles &&
                current.pos <= previous.end &&
                !((options.convert_to_base_pair_resolution || options.break_bands_at_multiples_of > 0) &&
                  ((static_cast<long long>(current.pos) + 1LL) %
                   static_cast<long long>(options.convert_to_base_pair_resolution
                                              ? 1 : options.break_bands_at_multiples_of)) == 0);
            if (same_block) {
                if (!sample_names.empty()) merge_record_samples(previous, item.record);
                previous.end = std::max(previous.end, item.record.end);
                const auto merged_end = previous.end;
                bcf_update_info_int32(output_header, previous.value, "END", &merged_end, 1);
                bcf_destroy(item.record.value);
                continue;
            }
            if (same_site) {
                std::vector<std::string> union_alleles{previous.allele_names.front()};
                bool has_non_ref = false;
                const auto append_concrete = [&](const std::vector<std::string>& alleles) {
                    for (std::size_t allele = 1; allele < alleles.size(); ++allele) {
                        if (alleles[allele] == "<NON_REF>") {
                            has_non_ref = true;
                            continue;
                        }
                        if (std::find(union_alleles.begin() + 1, union_alleles.end(),
                                      alleles[allele]) == union_alleles.end())
                            union_alleles.push_back(alleles[allele]);
                    }
                };
                if (options.native_genotype_materialization && !options.call_genotypes) {
                    append_concrete(previous.allele_names);
                    append_concrete(item.record.allele_names);
                } else {
                    append_concrete(item.record.allele_names);
                    append_concrete(previous.allele_names);
                }
                if (has_non_ref) union_alleles.emplace_back("<NON_REF>");
                remap_record_to_allele_union(output_header, previous, union_alleles,
                                             options.native_genotype_materialization,
                                             options.call_genotypes);
                remap_record_to_allele_union(output_header, item.record, union_alleles,
                                             options.native_genotype_materialization,
                                             options.call_genotypes);
                if (sample_names.empty()) {
                    bcf_destroy(item.record.value);
                    continue;
                }
                merge_record_samples(previous, item.record);
                if (previous.depth < 0) previous.depth = item.record.depth;
                else if (item.record.depth >= 0) previous.depth += item.record.depth;
                previous.end = std::max(previous.end, item.record.end);
                if (previous.depth >= 0)
                    bcf_update_info_int32(output_header, previous.value, "DP", &previous.depth, 1);
                if (previous.end > previous.pos + 1) {
                    const auto merged_end = previous.end;
                    bcf_update_info_int32(output_header, previous.value, "END", &merged_end, 1);
                }
                bcf_destroy(item.record.value);
                continue;
            }
            write_record(previous);
            pending_output = std::move(item.record);
        }
        if (pending_output.has_value()) write_record(*pending_output);
        if (bcf_close(output) != 0) throw std::runtime_error("cannot finalize output VCF/GVCF");
        output = nullptr;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0,
                                     &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const auto complete = file_complete(options.output) &&
            (index_path.empty() || file_complete(index_path));
        if (!complete) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: combined output or index is missing/empty");
        const auto manifest_path = options.manifest.empty()
            ? options.output + ".manifest.json" : options.manifest;
        const auto file_bytes = [](const std::string& path) -> std::uintmax_t {
            std::error_code error;
            const auto bytes = std::filesystem::file_size(path, error);
            return error ? 0 : bytes;
        };
        const auto wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"CombineGVCFs\","
                 << "\"implementation\":\"fastgatk-combine-gvcfs\",\"status\":\"prototype\","
                 << "\"execution_space\":\"Host\",\"determinism\":\"strict\","
                 << "\"primary_output\":\"" << json_escape(options.output)
                 << "\",\"primary_output_kind\":\"gvcf\",\"compatibility\":{"
                 << "\"site_level_merge\":true,\"k_way_merge\":true,\"stream_merge\":true,"
                 << "\"allele_union\":true,\"reference_block_merge\":true,\"format_sample_merge\":"
                 << (!sample_names.empty() ? "true" : "false") << ",\"vcf_index\":"
                 << (index_path.empty() ? "false" : "true") << ",\"interval_subset\":"
                 << (!options.regions.empty() ? "true" : "false")
                 << ",\"kokkos_pl_remap\":true,\"kokkos_allele_field_remap\":true"
                 << ",\"call_genotypes\":" << (options.call_genotypes ? "true" : "false")
                 << ",\"break_bands_at_multiples_of\":" << options.break_bands_at_multiples_of
                 << ",\"convert_to_base_pair_resolution\":"
                 << (options.convert_to_base_pair_resolution ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"lazy_reference_band_split\":true"
                 << ",\"indexed_interval_traversal\":"
                 << (!options.regions.empty() && summary.indexed_inputs > 0 ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":"
                 << ((!options.native_genotype_materialization || options.call_genotypes)
                         ? "true" : "false")
                 << "},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"gvcf\",\"complete\":"
                 << (file_complete(options.output) ? "true" : "false") << "}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) +
                     "\",\"kind\":\"vcf-index\",\"complete\":" +
                     (file_complete(index_path) ? "true}" : "false}"))
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size()
                 << ",\"input_records\":" << summary.input_records
                 << ",\"output_records\":" << summary.output_records
                 << ",\"skipped_records\":" << summary.skipped_records
                 << ",\"intervals\":" << options.regions.size()
                 << ",\"interval_skipped\":" << summary.interval_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"indexed_inputs\":" << summary.indexed_inputs
                 << ",\"indexed_interval_queries\":" << summary.indexed_interval_queries
                 << ",\"sample_count\":" << sample_names.size()
                 << ",\"stream_merge\":true,\"k_way_merge\":true"
                 << ",\"call_genotypes\":" << (options.call_genotypes ? "true" : "false")
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"lazy_reference_blocks\":" << summary.lazy_reference_blocks
                 << ",\"lazy_reference_segments\":" << summary.lazy_reference_segments
                 << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.genotype_calls
                 << ",\"genotype_kernel_prepare_seconds\":" << genotype_kernel_telemetry.genotype_prepare_seconds
                 << ",\"genotype_kernel_seconds\":" << genotype_kernel_telemetry.genotype_seconds
                 << ",\"genotype_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.genotype_execution_space) << "\""
                 << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                 << ",\"pl_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.pl_remap_prepare_seconds
                 << ",\"pl_remap_kernel_seconds\":" << genotype_kernel_telemetry.pl_remap_seconds
                 << ",\"pl_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.pl_remap_execution_space) << "\""
                 << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
                 << ",\"allele_field_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
                 << ",\"allele_field_remap_kernel_seconds\":" << genotype_kernel_telemetry.allele_field_remap_seconds
                 << ",\"allele_field_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space) << "\""
                 << ",\"max_inflight_records\":" << summary.max_inflight_records
                 << ",\"peak_queued_bytes\":" << summary.peak_queued_bytes
                 << ",\"safe_memory_budget_bytes\":" << safe_memory_budget
                 << ",\"output_bytes\":" << file_bytes(options.output)
                 << ",\"index_bytes\":" << (index_path.empty() ? 0 : file_bytes(index_path))
                 << ",\"wall_seconds\":" << wall_seconds << "}}\n";
        manifest.flush();
        if (!file_complete(manifest_path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: manifest is empty");
        std::cout << "{\"tool\":\"CombineGVCFs\",\"status\":\"prototype\","
                  << "\"stream_merge\":true,\"input_records\":" << summary.input_records
                  << ",\"output_records\":" << summary.output_records
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"sites_only_vcf_output\":"
                  << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"max_inflight_records\":" << summary.max_inflight_records
                  << ",\"peak_queued_bytes\":" << summary.peak_queued_bytes
                  << ",\"lazy_reference_blocks\":" << summary.lazy_reference_blocks
                  << ",\"lazy_reference_segments\":" << summary.lazy_reference_segments
                  << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.genotype_calls
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"allele_field_remap_kernel_calls\":"
                  << genotype_kernel_telemetry.allele_field_remap_calls << "}\n";
        for (auto& cursor : cursors) {
            if (cursor.iterator != nullptr) hts_itr_destroy(cursor.iterator);
            cursor.iterator = nullptr;
            free(cursor.indexed_line.s);
            cursor.indexed_line.s = nullptr;
            cursor.traversal_index.reset();
            for (auto& item : cursor.pending) bcf_destroy(item.value);
            cursor.pending.clear();
            if (cursor.split_block.has_value()) {
                bcf_destroy(cursor.split_block->value);
                cursor.split_block.reset();
            }
            bcf_destroy(cursor.input_record);
            cursor.input_record = nullptr;
            bcf_hdr_destroy(cursor.header);
            cursor.header = nullptr;
            bcf_close(cursor.file);
            cursor.file = nullptr;
        }
        bcf_hdr_destroy(output_header);
        output_header = nullptr;
        bcf_hdr_destroy(writer_header);
        writer_header = nullptr;
        if (reference_index) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        if (output) bcf_close(output);
        for (auto& cursor : cursors) {
            if (cursor.iterator != nullptr) hts_itr_destroy(cursor.iterator);
            cursor.iterator = nullptr;
            free(cursor.indexed_line.s);
            cursor.indexed_line.s = nullptr;
            cursor.traversal_index.reset();
            for (auto& item : cursor.pending) bcf_destroy(item.value);
            cursor.pending.clear();
            if (cursor.split_block.has_value()) {
                bcf_destroy(cursor.split_block->value);
                cursor.split_block.reset();
            }
            bcf_destroy(cursor.input_record);
            if (cursor.header) bcf_hdr_destroy(cursor.header);
            if (cursor.file) bcf_close(cursor.file);
        }
        if (output_header) bcf_hdr_destroy(output_header);
        if (writer_header) bcf_hdr_destroy(writer_header);
        if (reference_index) fai_destroy(reference_index);
        throw;
    }
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    if (options.stream_merge)
        return run_streaming_tool(options, resources);
    const auto started = std::chrono::steady_clock::now();
    bcf_hdr_t* output_header = nullptr;
    bcf_hdr_t* writer_header = nullptr;
    faidx_t* reference_index = nullptr;
    std::vector<Record> records;
    std::vector<std::string> sample_names;
    std::set<std::string> sample_name_set;
    std::uint64_t input_records = 0;
    std::uint64_t skipped = 0;
    std::uint64_t interval_skipped = 0;
    fastgatk::io::IntervalFileStats interval_file_stats;
    std::vector<Region> regions;
    bool intervals_parsed = false;
    const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    try {
        if (!options.reference.empty()) {
            reference_index = fai_load(options.reference.c_str());
            if (!reference_index)
                throw std::runtime_error("BAD_INPUT: cannot load reference FASTA index: " + options.reference);
        }
        for (const auto& input_path : options.inputs) {
            htsFile* input = bcf_open(input_path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) {
                bcf_close(input);
                throw std::runtime_error("BAD_INPUT: cannot read VCF header: " + input_path);
            }
            if (!output_header) output_header = make_output_header(header);
            const auto sample_base = sample_names.size();
            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (name == nullptr || !sample_name_set.insert(name).second) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error(
                        "BAD_INPUT: duplicate or invalid sample name across CombineGVCFs inputs");
                }
                sample_names.emplace_back(name);
                if (bcf_hdr_add_sample(output_header, name) != 0) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot add sample to output header");
                }
            }
            if (sample_names.size() != sample_base) {
                for (auto& existing : records) existing.samples.resize(sample_names.size());
            }
            if (!intervals_parsed) {
                regions.reserve(options.regions.size());
                bool first_selector = true;
                for (const auto& text : options.regions) {
                    append_combine_region_selector_with_rule(
                        text, header, regions, interval_file_stats,
                        options.interval_set_rule, first_selector);
                    first_selector = false;
                }
                intervals_parsed = true;
            }
            bcf1_t* record = bcf_init();
            if (!record) {
                bcf_hdr_destroy(header); bcf_close(input);
                throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
            }
            while (bcf_read(input, header, record) == 0) {
                ++input_records;
                bcf_unpack(record, BCF_UN_ALL);
                if (!in_regions(header, record, regions)) {
                    ++interval_skipped;
                    continue;
                }
                if (record->rid < 0 || record->n_allele < 2) { ++skipped; continue; }
                const auto* contig = bcf_hdr_id2name(header, record->rid);
                const auto output_rid = bcf_hdr_name2id(output_header, contig);
                if (output_rid < 0) { ++skipped; continue; }
                std::string allele_string;
                for (int allele = 0; allele < record->n_allele; ++allele) {
                    if (allele != 0) allele_string.push_back(',');
                    allele_string += record->d.allele[allele];
                }
                const auto end = info_int(header, record, "END", record->pos + 1);
                const auto depth = info_int(header, record, "DP", -1);
                auto* output_record = bcf_init();
                if (!output_record || bcf_update_alleles_str(output_header, output_record, allele_string.c_str()) != 0) {
                    bcf_destroy(output_record); bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize combined record");
                }
                output_record->rid = output_rid;
                output_record->pos = record->pos;
                bcf_float_set_missing(output_record->qual);
                if (depth >= 0) bcf_update_info_int32(output_header, output_record, "DP", &depth, 1);
                if (end > record->pos + 1)
                    bcf_update_info_int32(output_header, output_record, "END", &end, 1);
                bcf_unpack(output_record, BCF_UN_STR);
                const auto position = static_cast<int>(record->pos);
                Record item;
                item.value = output_record;
                item.rid = output_rid;
                item.pos = position;
                item.end = end;
                item.depth = depth;
                item.alleles = allele_string;
                item.allele_names = split_alleles(allele_string);
                item.samples.resize(sample_names.size());
                for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample)
                    item.samples[sample_base + static_cast<std::size_t>(sample)] =
                        read_sample_call(header, record, sample, header->n[BCF_DT_SAMPLE],
                                         record->n_allele);
                // GATK's CombineGVCFs annotates site DP from FORMAT/DP when
                // an input record has no INFO/DP.  Compute that fallback
                // after decoding the sample calls so --call-genotypes and
                // the default no-call path expose the same site annotation.
                if (item.depth < 0 && !is_reference_block(item)) {
                    long long depth_sum = 0;
                    bool has_depth = false;
                    for (const auto& call : item.samples) {
                        if (call.depth == bcf_int32_missing || call.depth < 0) continue;
                        depth_sum += call.depth;
                        has_depth = true;
                    }
                    if (has_depth && depth_sum <= std::numeric_limits<int>::max())
                        item.depth = static_cast<int>(depth_sum);
                    if (item.depth >= 0 &&
                        bcf_update_info_int32(output_header, output_record, "DP", &item.depth, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize DP annotation");
                }
                // GATK can split reference-confidence bands at deterministic
                // 1-based genomic boundaries.  Keep each segment as an
                // independent Record through the merge stage; the merge
                // predicate below explicitly refuses to coalesce across a
                // requested breakpoint.
                const auto break_multiple = options.convert_to_base_pair_resolution
                    ? 1 : options.break_bands_at_multiples_of;
                const bool split_block = is_reference_block(item) &&
                    item.end > item.pos + 1 && break_multiple > 0;
                std::vector<int> segment_starts{item.pos};
                std::vector<int> segment_ends{item.end};
                if (split_block) {
                    segment_starts.clear();
                    segment_ends.clear();
                    int segment_start = item.pos;
                    // A breakpoint at genomic coordinate N means that the
                    // next segment starts at zero-based N-1 (before that
                    // base), matching CombineGVCFs.getIntermediateStopSites.
                    for (long long coordinate =
                             (static_cast<long long>(item.pos) + 2LL + break_multiple - 1LL) /
                             break_multiple * break_multiple;
                         coordinate <= static_cast<long long>(item.end);
                         coordinate += break_multiple) {
                        const auto boundary = static_cast<int>(coordinate - 1LL);
                        if (boundary <= segment_start || boundary >= item.end) continue;
                        segment_starts.push_back(segment_start);
                        segment_ends.push_back(boundary);
                        segment_start = boundary;
                    }
                    segment_starts.push_back(segment_start);
                    segment_ends.push_back(item.end);
                }
                const auto item_rid = item.rid;
                const auto item_samples = item.samples;
                const auto item_allele_names = item.allele_names;
                bcf1_t* original_value = item.value;
                item.value = nullptr;
                for (std::size_t segment = 0; segment < segment_starts.size(); ++segment) {
                    Record output_segment;
                    output_segment.value = segment == 0 ? original_value : bcf_dup(original_value);
                    if (!output_segment.value)
                        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot split reference block");
                    output_segment.rid = item_rid;
                    output_segment.alleles = allele_string;
                    output_segment.allele_names = item_allele_names;
                    output_segment.samples = item_samples;
                    output_segment.depth = item.depth;
                    output_segment.pos = segment_starts[segment];
                    output_segment.end = segment_ends[segment];
                    output_segment.value->pos = output_segment.pos;
                    std::string reference_allele = output_segment.allele_names.front();
                    if (reference_index && split_block) {
                        const char* contig_name = bcf_hdr_id2name(output_header, output_segment.rid);
                        int fetched_length = 0;
                        char* fetched = contig_name == nullptr ? nullptr :
                            faidx_fetch_seq(reference_index, contig_name, output_segment.pos,
                                            output_segment.pos, &fetched_length);
                        if (fetched != nullptr && fetched_length == 1)
                            reference_allele.assign(1, fetched[0]);
                        free(fetched);
                    }
                    if (split_block && !reference_allele.empty()) {
                        output_segment.allele_names.front() = reference_allele;
                        output_segment.alleles = reference_allele + ",<NON_REF>";
                        if (bcf_update_alleles_str(output_header, output_segment.value,
                                                   output_segment.alleles.c_str()) != 0)
                            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split reference allele");
                    }
                    if (output_segment.end > output_segment.pos + 1) {
                        const auto segment_end = output_segment.end;
                        if (bcf_update_info_int32(output_header, output_segment.value,
                                                  "END", &segment_end, 1) != 0)
                            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split END");
                    } else if (bcf_update_info_int32(output_header, output_segment.value,
                                                     "END", nullptr, 0) != 0) {
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear split END");
                    }
                    output_segment.key = record_key(output_segment);
                    records.push_back(std::move(output_segment));
                    if (safe_memory_budget != 0 && records.size() > safe_memory_budget / 512) {
                        bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                        throw std::runtime_error("RESOURCE_EXHAUSTED: combined GVCF staging exceeds safe memory budget");
                    }
                }
            }
            bcf_destroy(record);
            bcf_hdr_destroy(header);
            bcf_close(input);
        }
        std::sort(records.begin(), records.end(), [](const Record& left, const Record& right) {
            if (left.rid != right.rid) return left.rid < right.rid;
            if (left.pos != right.pos) return left.pos < right.pos;
            if (left.end != right.end) return left.end < right.end;
            return left.key < right.key;
        });
        std::vector<Record> merged;
        for (auto& current : records) {
            // Materialize the source call even when a locus has no duplicate
            // record to trigger allele-union remapping.  This keeps the
            // single-input path on the same arbitrary-ploidy Kokkos GT/GQ
            // implementation as the multi-shard path.
            remap_record_to_allele_union(output_header, current, current.allele_names,
                                         options.native_genotype_materialization,
                                         options.call_genotypes);
            if (!merged.empty()) {
                auto& previous = merged.back();
                const bool same_site = !is_reference_block(previous) && !is_reference_block(current) &&
                    previous.rid == current.rid && previous.pos == current.pos &&
                    !previous.allele_names.empty() && !current.allele_names.empty() &&
                    previous.allele_names.front() == current.allele_names.front();
                const bool same_block = is_reference_block(previous) && is_reference_block(current) &&
                    previous.rid == current.rid && previous.alleles == current.alleles &&
                    // current.pos is zero-based and END is one-based
                    // inclusive; this admits only truly adjacent blocks.
                    current.pos <= previous.end &&
                    // Do not undo an explicit --break-bands-at-multiples-of
                    // split.  current.pos is zero-based, so a segment that
                    // starts at 1-based coordinate N has (pos + 1) == N.
                    !((options.convert_to_base_pair_resolution ||
                       options.break_bands_at_multiples_of > 0) &&
                      ((static_cast<long long>(current.pos) + 1LL) %
                       static_cast<long long>(options.convert_to_base_pair_resolution
                                                  ? 1 : options.break_bands_at_multiples_of)) == 0);
                if (same_block) {
                    if (!sample_names.empty()) merge_record_samples(previous, current);
                    previous.end = std::max(previous.end, current.end);
                    const auto merged_end = previous.end;
                    bcf_update_info_int32(output_header, previous.value, "END", &merged_end, 1);
                    bcf_destroy(current.value);
                    continue;
                }
                if (same_site) {
                    // Keep concrete ALTs in first-seen order and place the
                    // gVCF `<NON_REF>` sentinel last, matching GATK's
                    // allele ordering for a joint site.
                    std::vector<std::string> union_alleles{previous.allele_names.front()};
                    bool has_non_ref = false;
                    const auto append_concrete = [&](const std::vector<std::string>& alleles) {
                        for (std::size_t allele = 1; allele < alleles.size(); ++allele) {
                            if (alleles[allele] == "<NON_REF>") {
                                has_non_ref = true;
                                continue;
                            }
                            if (std::find(union_alleles.begin() + 1, union_alleles.end(),
                                          alleles[allele]) == union_alleles.end())
                                union_alleles.push_back(alleles[allele]);
                        }
                    };
                    if (options.native_genotype_materialization && !options.call_genotypes) {
                        append_concrete(previous.allele_names);
                        append_concrete(current.allele_names);
                    } else {
                        // GATK's merge tree prepends the later shard's
                        // concrete ALTs; this is observable in PL/AD order.
                        append_concrete(current.allele_names);
                        append_concrete(previous.allele_names);
                    }
                    if (has_non_ref) union_alleles.emplace_back("<NON_REF>");
                    remap_record_to_allele_union(output_header, previous, union_alleles,
                                                 options.native_genotype_materialization,
                                                 options.call_genotypes);
                    remap_record_to_allele_union(output_header, current, union_alleles,
                                                 options.native_genotype_materialization,
                                                 options.call_genotypes);
                    // With no FORMAT/sample columns this is a true duplicate
                    // site and remains deduplicated. With samples, each input
                    // contributes a disjoint sample and the calls are merged.
                    if (sample_names.empty()) {
                        bcf_destroy(current.value);
                        continue;
                    }
                    merge_record_samples(previous, current);
                    if (previous.depth < 0) previous.depth = current.depth;
                    else if (current.depth >= 0) previous.depth += current.depth;
                    previous.end = std::max(previous.end, current.end);
                    if (previous.depth >= 0)
                        bcf_update_info_int32(output_header, previous.value, "DP", &previous.depth, 1);
                    if (previous.end > previous.pos + 1) {
                        const auto merged_end = previous.end;
                        bcf_update_info_int32(output_header, previous.value, "END", &merged_end, 1);
                    }
                    bcf_destroy(current.value);
                    continue;
                }
            }
            merged.push_back(std::move(current));
        }
        records.clear();
        records = std::move(merged);
        if (!sample_names.empty() && bcf_hdr_add_sample(output_header, nullptr) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot terminate sample header");
        if (bcf_hdr_sync(output_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync merged sample header");
        writer_header = bcf_hdr_dup(output_header);
        if (!writer_header)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate CombineGVCFs writer header");
        if (options.sites_only_vcf_output && bcf_hdr_set_samples(writer_header, nullptr, 0) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create sites-only writer header");
        for (auto& record : records)
            materialize_samples(output_header, record, sample_names.size());
        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        htsFile* output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open output VCF/GVCF: " + options.output);
        if (bcf_hdr_write(output, writer_header) != 0) {
            bcf_close(output);
            throw std::runtime_error("cannot write output VCF/GVCF header");
        }
        for (const auto& item : records) {
            if (options.sites_only_vcf_output && is_reference_block(item))
                bcf_update_info_int32(output_header, item.value, "DP", nullptr, 0);
            if (options.sites_only_vcf_output && bcf_subset_format(writer_header, item.value) != 0) {
                bcf_close(output);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset sites-only FORMAT payload");
            }
            if (bcf_write(output, writer_header, item.value) != 0) {
                bcf_close(output);
                throw std::runtime_error("cannot write output VCF/GVCF record");
            }
        }
        bcf_close(output);
        bcf_hdr_destroy(writer_header);
        writer_header = nullptr;
        std::string index_path;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const bool primary_complete = file_complete(options.output);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !index_complete)
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: combined output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        const auto file_bytes = [](const std::string& path) -> std::uintmax_t {
            std::error_code error;
            const auto bytes = std::filesystem::file_size(path, error);
            return error ? 0 : bytes;
        };
        const auto wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"CombineGVCFs\",\"implementation\":\"fastgatk-combine-gvcfs\",\"status\":\"prototype\","
                 << "\"execution_space\":\"Host\",\"determinism\":\"strict\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"gvcf\","
                 << "\"compatibility\":{\"site_level_merge\":true,\"allele_union\":true,\"gq_rematerialization\":true,\"reference_block_merge\":true,\"format_sample_merge\":"
                 << (!sample_names.empty() ? "true" : "false") << ",\"vcf_index\":"
                 << (index_path.empty() ? "false" : "true") << ",\"interval_subset\":"
                 << (!options.regions.empty() ? "true" : "false")
                 << ",\"kokkos_genotype_assignment\":"
                 << (options.native_genotype_materialization ? "true" : "false")
                 << ",\"kokkos_pl_remap\":true,\"kokkos_allele_field_remap\":true"
                 << ",\"call_genotypes\":"
                 << (options.call_genotypes ? "true" : "false")
                 << ",\"break_bands_at_multiples_of\":"
                 << options.break_bands_at_multiples_of
                 << ",\"convert_to_base_pair_resolution\":"
                 << (options.convert_to_base_pair_resolution ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"gatk_format_semantics\":"
                 << ((!options.native_genotype_materialization || options.call_genotypes) ? "true" : "false")
                 << ",\"native_genotype_materialization\":"
                 << (options.native_genotype_materialization ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":"
                 << ((!options.native_genotype_materialization || options.call_genotypes) ? "true" : "false")
                 << "},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"gvcf\",\"complete\":"
                 << (primary_complete ? "true" : "false") << "}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":" +
                     (index_complete ? "true}" : "false}"))
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size() << ",\"input_records\":" << input_records
                 << ",\"output_records\":" << records.size() << ",\"skipped_records\":" << skipped
                 << ",\"intervals\":" << options.regions.size()
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"sample_count\":" << sample_names.size()
                 << ",\"call_genotypes\":" << (options.call_genotypes ? "true" : "false")
                 << ",\"break_bands_at_multiples_of\":" << options.break_bands_at_multiples_of
                 << ",\"convert_to_base_pair_resolution\":"
                 << (options.convert_to_base_pair_resolution ? "true" : "false")
                 << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.genotype_calls
                 << ",\"genotype_kernel_prepare_seconds\":" << genotype_kernel_telemetry.genotype_prepare_seconds
                 << ",\"genotype_kernel_seconds\":" << genotype_kernel_telemetry.genotype_seconds
                 << ",\"genotype_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.genotype_execution_space) << "\""
                 << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                 << ",\"pl_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.pl_remap_prepare_seconds
                 << ",\"pl_remap_kernel_seconds\":" << genotype_kernel_telemetry.pl_remap_seconds
                 << ",\"pl_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.pl_remap_execution_space) << "\""
                 << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
                 << ",\"allele_field_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
                 << ",\"allele_field_remap_kernel_seconds\":" << genotype_kernel_telemetry.allele_field_remap_seconds
                 << ",\"allele_field_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space)
                 << "\",\"output_bytes\":" << file_bytes(options.output)
                 << ",\"index_bytes\":" << (index_path.empty() ? 0 : file_bytes(index_path))
                 << ",\"wall_seconds\":" << wall_seconds
                 << ",\"native_genotype_materialization\":"
                 << (options.native_genotype_materialization ? "true" : "false") << "}}\n";
        std::cout << "{\"tool\":\"CombineGVCFs\",\"status\":\"prototype\",\"input_files\":"
                  << options.inputs.size() << ",\"input_records\":" << input_records
                  << ",\"output_records\":" << records.size()
                  << ",\"interval_skipped\":" << interval_skipped
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"sites_only_vcf_output\":"
                  << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"sample_count\":" << sample_names.size()
                  << ",\"call_genotypes\":" << (options.call_genotypes ? "true" : "false")
                  << ",\"break_bands_at_multiples_of\":" << options.break_bands_at_multiples_of
                  << ",\"convert_to_base_pair_resolution\":"
                  << (options.convert_to_base_pair_resolution ? "true" : "false")
                  << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.genotype_calls
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls << "}\n";
        destroy_records(records);
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(writer_header);
        if (reference_index) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        destroy_records(records);
        if (output_header) bcf_hdr_destroy(output_header);
        if (writer_header) bcf_hdr_destroy(writer_header);
        if (reference_index) fai_destroy(reference_index);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for CombineGVCFs");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        Kokkos::initialize();
        initialized = true;
        const auto options = parse(argc, argv);
        const auto result = run_tool(options, fastgatk::runtime::ResourceSnapshot::probe());
        Kokkos::finalize();
        initialized = false;
        return result;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "fastgatk-combine-gvcfs: " << error.what() << '\n';
        return 2;
    }
}
