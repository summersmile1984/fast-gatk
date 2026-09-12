#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/java_numeric.hpp"
#include "fastgatk/core/plan.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string variants;
    std::string output;
    std::string manifest;
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> excluded_regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    int interval_padding = 0;
    int exclusion_padding = 0;
    double min_af = 0.01;
    double max_af = 0.20;
    int minimum_mapping_quality = 50;
    bool require_read_length = false;
    std::uint32_t min_read_length = 1;
    std::uint32_t max_read_length = std::numeric_limits<std::uint32_t>::max();
    std::size_t batch_records = 1024;
    int threads = 1;
    bool include_duplicates = false;
    bool disable_tool_default_read_filters = false;
    bool disable_sequence_dictionary_validation = false;
    std::vector<std::string> read_filters;
    std::vector<std::string> inverted_read_filters;
    std::vector<std::string> disabled_read_filters;
};

std::string format_double(double value) {
    return fastgatk::io::java_double(value);
}

// GetPileupSummaries inherits the GATK LocusWalker default read-filter set.
// Keep the filter state explicit so --read-filter/--disable-read-filter never
// become silently ignored compatibility flags.
struct ReadFilterConfig {
    bool mapping_quality = true;
    bool mapping_quality_available = true;
    bool mapping_quality_not_zero = true;
    bool mapped = true;
    bool secondary = true;
    bool supplementary = true;
    bool duplicate = true;
    bool vendor_quality = true;
    bool nonzero_reference = true;
    bool mate_same_contig = true;
    bool good_cigar = true;
    bool wellformed = true;
    bool read_length = false;
};

enum class ReadFilterKind : std::uint8_t {
    MappingQuality,
    MappingQualityAvailable,
    MappingQualityNotZero,
    Mapped,
    PrimaryLine,
    NotSecondary,
    NotSupplementary,
    NotDuplicate,
    PassesVendorQuality,
    NonZeroReferenceLength,
    MateOnSameContig,
    GoodCigar,
    Wellformed,
    ReadLength,
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

bool parse_bool_text(const std::string& value, const char* name) {
    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    if (normalized == "true" || normalized == "1") return true;
    if (normalized == "false" || normalized == "0") return false;
    throw std::invalid_argument(std::string(name) + " must be true or false");
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

ReadFilterKind read_filter_kind(const std::string& name) {
    if (name == "MappingQualityReadFilter") return ReadFilterKind::MappingQuality;
    if (name == "MappingQualityAvailableReadFilter") return ReadFilterKind::MappingQualityAvailable;
    if (name == "MappingQualityNotZeroReadFilter") return ReadFilterKind::MappingQualityNotZero;
    if (name == "MappedReadFilter") return ReadFilterKind::Mapped;
    if (name == "PrimaryLineReadFilter") return ReadFilterKind::PrimaryLine;
    if (name == "NotSecondaryAlignmentReadFilter") return ReadFilterKind::NotSecondary;
    if (name == "NotSupplementaryAlignmentReadFilter") return ReadFilterKind::NotSupplementary;
    if (name == "NotDuplicateReadFilter") return ReadFilterKind::NotDuplicate;
    if (name == "PassesVendorQualityCheckReadFilter") return ReadFilterKind::PassesVendorQuality;
    if (name == "NonZeroReferenceLengthAlignmentReadFilter") return ReadFilterKind::NonZeroReferenceLength;
    if (name == "MateOnSameContigOrNoMappedMateReadFilter") return ReadFilterKind::MateOnSameContig;
    if (name == "GoodCigarReadFilter") return ReadFilterKind::GoodCigar;
    if (name == "WellformedReadFilter") return ReadFilterKind::Wellformed;
    if (name == "ReadLengthReadFilter") return ReadFilterKind::ReadLength;
    throw std::invalid_argument("UNSUPPORTED_PARAMETER: --read-filter " + name);
}

void set_read_filter(ReadFilterConfig& config, ReadFilterKind kind, bool enabled) {
    switch (kind) {
        case ReadFilterKind::MappingQuality: config.mapping_quality = enabled; break;
        case ReadFilterKind::MappingQualityAvailable: config.mapping_quality_available = enabled; break;
        case ReadFilterKind::MappingQualityNotZero: config.mapping_quality_not_zero = enabled; break;
        case ReadFilterKind::Mapped: config.mapped = enabled; break;
        case ReadFilterKind::PrimaryLine:
            config.secondary = enabled;
            config.supplementary = enabled;
            break;
        case ReadFilterKind::NotSecondary: config.secondary = enabled; break;
        case ReadFilterKind::NotSupplementary: config.supplementary = enabled; break;
        case ReadFilterKind::NotDuplicate: config.duplicate = enabled; break;
        case ReadFilterKind::PassesVendorQuality: config.vendor_quality = enabled; break;
        case ReadFilterKind::NonZeroReferenceLength: config.nonzero_reference = enabled; break;
        case ReadFilterKind::MateOnSameContig: config.mate_same_contig = enabled; break;
        case ReadFilterKind::GoodCigar: config.good_cigar = enabled; break;
        case ReadFilterKind::Wellformed: config.wellformed = enabled; break;
        case ReadFilterKind::ReadLength: config.read_length = enabled; break;
    }
}

ReadFilterConfig effective_read_filters(const Options& options) {
    ReadFilterConfig config;
    if (options.disable_tool_default_read_filters) config = ReadFilterConfig{
        false, false, false, false, false, false, false, false, false, false, false, false, false};
    if (options.include_duplicates) config.duplicate = false;
    for (const auto& name : options.read_filters)
        set_read_filter(config, read_filter_kind(name), true);
    for (const auto& name : options.disabled_read_filters)
        set_read_filter(config, read_filter_kind(name), false);
    return config;
}

Options parse(int argc, char** argv) {
    Options options;
    const auto optional_bool = [&](int& index, int argc_local, char** argv_local,
                                   const std::string& argument, const char* name) {
        const auto inline_value = option_value(argument, name);
        if (!inline_value.empty()) return parse_bool_text(inline_value, name);
        if (argument == name && index + 1 < argc_local) {
            const std::string next(argv_local[index + 1]);
            if (next == "true" || next == "false" || next == "1" || next == "0") {
                ++index;
                return parse_bool_text(next, name);
            }
        }
        return true;
    };
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-get-pileup-summaries (GATK-compatible native prototype)\n"
                         "  -I, --input FILE              BAM/CRAM/SAM input\n"
                         "  -V, --variant FILE            biallelic SNP resource with AF\n"
                         "  -L, --intervals REGION        sites/interval selector (repeatable)\n"
                         "  -XL, --exclude-intervals REGION  exclude sites/intervals (repeatable)\n"
                         "      --interval-set-rule RULE UNION (default) or INTERSECTION\n"
                         "      --interval-padding N       padding for included intervals\n"
                         "      --interval-exclusion-padding N  padding for excluded intervals\n"
                         "  -O, --output FILE             pileup summary table\n"
                         "      --min-af FLOAT            minimum population AF (default 0.01)\n"
                         "      --max-af FLOAT            maximum population AF (default 0.20)\n"
                         "      --minimum-mapping-quality N (default 50)\n"
                         "      --min-read-length N       ReadLengthReadFilter lower bound\n"
                         "      --max-read-length N       ReadLengthReadFilter upper bound\n"
                         "      --output-manifest FILE   OutputManifest JSON\n"
                         "  -RF, --read-filter NAME       add a supported GATK read filter\n"
                         "  -DF, --disable-read-filter NAME  disable one default filter\n";
            std::exit(0);
        } else if (argument == "-I" || is_option(argument, "--input"))
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        else if (argument == "-V" || is_option(argument, "--variant"))
            options.variants = require_value(index, argc, argv, argument, "--variant", "-V");
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                is_option(argument, "--region") ? "--region" :
                is_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (argument == "-XL" || is_option(argument, "--exclude-intervals"))
            options.excluded_regions.push_back(require_value(
                index, argc, argv, argument, "--exclude-intervals", "-XL"));
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-ip" || is_option(argument, "--interval-padding"))
            options.interval_padding = std::stoi(require_value(
                index, argc, argv, argument, "--interval-padding", "-ip"));
        else if (argument == "-ixp" || is_option(argument, "--interval-exclusion-padding"))
            options.exclusion_padding = std::stoi(require_value(
                index, argc, argv, argument, "--interval-exclusion-padding", "-ixp"));
        else if (is_option(argument, "--min-af") || is_option(argument, "--minimum-population-allele-frequency"))
            options.min_af = std::stod(require_value(
                index, argc, argv, argument,
                argument.rfind("--min-af", 0) == 0 ? "--min-af" : "--minimum-population-allele-frequency"));
        else if (is_option(argument, "--max-af") || is_option(argument, "--maximum-population-allele-frequency"))
            options.max_af = std::stod(require_value(
                index, argc, argv, argument,
                argument.rfind("--max-af", 0) == 0 ? "--max-af" : "--maximum-population-allele-frequency"));
        else if (is_option(argument, "--minimum-mapping-quality"))
            options.minimum_mapping_quality = std::stoi(require_value(
                index, argc, argv, argument, "--minimum-mapping-quality"));
        else if (is_option(argument, "--min-read-length")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--min-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-read-length is out of range");
            options.min_read_length = static_cast<std::uint32_t>(value);
            options.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        } else if (is_option(argument, "--max-read-length")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--max-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-read-length is out of range");
            options.max_read_length = static_cast<std::uint32_t>(value);
            options.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        }
        else if (is_option(argument, "--batch-records"))
            options.batch_records = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--batch-records")));
        else if (is_option(argument, "--threads"))
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--quiet") {
            // Progress/logging compatibility flag.
        } else if (is_option(argument, "--disable-sequence-dictionary-validation")) {
            options.disable_sequence_dictionary_validation = optional_bool(
                index, argc, argv, argument, "--disable-sequence-dictionary-validation");
        } else if (is_option(argument, "--disable-tool-default-read-filters")) {
            options.disable_tool_default_read_filters = optional_bool(
                index, argc, argv, argument, "--disable-tool-default-read-filters");
        } else if (is_option(argument, "--include-duplicates")) {
            options.include_duplicates = optional_bool(index, argc, argv, argument, "--include-duplicates");
        } else if (argument == "-RF" || is_option(argument, "--read-filter")) {
            options.read_filters.push_back(require_value(index, argc, argv, argument, "--read-filter", "-RF"));
        } else if (argument == "-DF" || is_option(argument, "--disable-read-filter")) {
            options.disabled_read_filters.push_back(require_value(
                index, argc, argv, argument, "--disable-read-filter", "-DF"));
        } else if (argument == "-XRF" || is_option(argument, "--inverted-read-filter")) {
            options.inverted_read_filters.push_back(require_value(
                index, argc, argv, argument, "--inverted-read-filter", "-XRF"));
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--seconds-between-progress-updates")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.variants.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.regions.empty()) throw std::invalid_argument("-L/--intervals is required");
    if (!std::isfinite(options.min_af) || !std::isfinite(options.max_af) ||
        options.min_af < 0.0 || options.max_af > 1.0 || options.min_af >= options.max_af)
        throw std::invalid_argument("population AF bounds must satisfy 0 <= min-af < max-af <= 1");
    if (options.minimum_mapping_quality < 0 || options.batch_records == 0 || options.threads < 1)
        throw std::invalid_argument("mapping quality, batch-records and threads must be positive");
    if (options.interval_padding < 0 || options.exclusion_padding < 0)
        throw std::invalid_argument("interval padding must be non-negative");
    if (options.require_read_length && options.min_read_length > options.max_read_length)
        throw std::invalid_argument("read length filter minimum exceeds maximum");
    // Validate names during parse, even when a filter is immediately disabled,
    // so unsupported filters fail closed instead of silently changing counts.
    (void)effective_read_filters(options);
    for (const auto& name : options.inverted_read_filters) (void)read_filter_kind(name);
    return options;
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

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

#if FASTGATK_HAS_HTSLIB

struct Site {
    std::string contig;
    std::string ref;
    std::string alt;
    int tid = -1;
    int position = -1;  // zero-based
    double allele_frequency = 0.0;
    int ref_code = -1;
    int alt_code = -1;
};

int base_code(char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

std::string bam_sample_name(const std::string& path, const std::string& reference) {
    htsFile* file = sam_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open reads: " + path);
    if (!reference.empty() && hts_set_fai_filename(file, reference.c_str()) != 0) {
        sam_close(file);
        throw std::runtime_error("BAD_INPUT: cannot configure reference for CRAM: " + reference);
    }
    sam_hdr_t* header = sam_hdr_read(file);
    if (!header) {
        sam_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read read header: " + path);
    }
    std::string sample;
    kstring_t line{0, 0, nullptr};
    const int count = sam_hdr_count_lines(header, "RG");
    for (int index = 0; index < count && sample.empty(); ++index) {
        if (sam_hdr_find_line_pos(header, "RG", index, &line) != 0) continue;
        std::istringstream fields(std::string(line.s == nullptr ? "" : line.s, line.l));
        std::string field;
        while (fields >> field) {
            if (field.rfind("SM:", 0) == 0 && field.size() > 3) {
                sample = field.substr(3);
                break;
            }
        }
    }
    free(line.s);
    sam_hdr_destroy(header);
    sam_close(file);
    if (sample.empty())
        throw std::runtime_error("BAD_INPUT: reads header has no sample (missing @RG SM tag)");
    return sample;
}

using LocusKey = std::uint64_t;

LocusKey locus_key(int tid, int position) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(tid)) << 32) |
           static_cast<std::uint32_t>(position);
}

bool is_variant_file(const std::string& path) {
    const auto lower = [&]() {
        std::string value = path;
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return value;
    }();
    return fastgatk::io::has_suffix_ci(lower, ".vcf") || fastgatk::io::has_suffix_ci(lower, ".vcf.gz") ||
           fastgatk::io::has_suffix_ci(lower, ".bcf") || fastgatk::io::has_suffix_ci(lower, ".bcf.gz");
}

void validate_variant_dictionary(const bcf_hdr_t* header,
                                 const fastgatk::io::HeaderSummary& bam_header,
                                 bool disabled) {
    if (disabled || header == nullptr || bam_header.contigs.empty()) return;
    std::unordered_map<std::string, std::int64_t> bam_lengths;
    for (std::size_t tid = 0; tid < bam_header.contigs.size(); ++tid)
        bam_lengths.emplace(bam_header.contigs[tid],
                            tid < bam_header.contig_lengths.size() ? bam_header.contig_lengths[tid] : -1);
    for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
        const auto* name = bcf_hdr_id2name(header, rid);
        if (name == nullptr) continue;
        const auto found = bam_lengths.find(name);
        if (found == bam_lengths.end()) continue;
        const auto* id_info = header->id[BCF_DT_CTG][rid].val;
        const auto variant_length = id_info == nullptr ? std::int64_t{0}
                                                        : static_cast<std::int64_t>(id_info->info[0]);
        if (variant_length > 0 && found->second > 0 && variant_length != found->second)
            throw std::runtime_error("BAD_INPUT: sequence dictionary length mismatch for contig " +
                                     std::string(name));
    }
}

std::unordered_set<std::string> load_interval_variant_keys(
    const std::vector<std::string>& regions,
    fastgatk::io::HtsIntervalSetRule rule,
    bool& active) {
    active = false;
    std::unordered_set<std::string> keys;
    bool first_variant_selector = true;
    for (const auto& selector : regions) {
        if (!is_variant_file(selector)) continue;
        active = true;
        std::unordered_set<std::string> selector_keys;
        htsFile* file = bcf_open(selector.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval VCF: " + selector);
        bcf_hdr_t* header = bcf_hdr_read(file);
        bcf1_t* record = header ? bcf_init() : nullptr;
        if (!header || !record) {
            if (record) bcf_destroy(record);
            if (header) bcf_hdr_destroy(header);
            bcf_close(file);
            throw std::runtime_error("BAD_INPUT: cannot read interval VCF: " + selector);
        }
        while (fastgatk::io::read_variant_record(file, header, record, selector) == 0) {
            const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
            if (contig) selector_keys.emplace(std::string(contig) + ":" + std::to_string(record->pos));
            bcf_clear(record);
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(file);
        if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_variant_selector) {
            keys.insert(selector_keys.begin(), selector_keys.end());
            first_variant_selector = false;
            continue;
        }
        std::unordered_set<std::string> intersection;
        for (const auto& key : keys)
            if (selector_keys.find(key) != selector_keys.end()) intersection.insert(key);
        keys.swap(intersection);
        first_variant_selector = false;
    }
    return keys;
}

struct VariantTraversalIndex {
    hts_idx_t* index = nullptr;
    tbx_t* tabix = nullptr;

    VariantTraversalIndex() = default;
    VariantTraversalIndex(const VariantTraversalIndex&) = delete;
    VariantTraversalIndex& operator=(const VariantTraversalIndex&) = delete;
    ~VariantTraversalIndex() {
        if (tabix != nullptr) tbx_destroy(tabix);
        else if (index != nullptr) hts_idx_destroy(index);
    }

    int tid_for(const bcf_hdr_t* header, int rid) const {
        if (tabix == nullptr) return rid;
        const auto* name = rid >= 0 ? bcf_hdr_id2name(header, rid) : nullptr;
        return name == nullptr ? -1 : tbx_name2id(tabix, name);
    }
};

// HTSlib's bcf_index_load() intentionally accepts CSI only.  GATK's VCF
// resources are commonly distributed with a tabix (.tbi) index.  A TBI index
// has its own compact sequence-ID namespace, so retain the tbx wrapper for
// name-to-ID translation before using the generic BCF iterator.
std::unique_ptr<VariantTraversalIndex> load_variant_index(const std::string& filename) {
    const std::filesystem::path csi_path = filename + ".csi";
    const std::filesystem::path tbi_path = filename + ".tbi";
    if (std::filesystem::exists(csi_path)) {
        if (auto* index = hts_idx_load3(filename.c_str(), csi_path.c_str(),
                                         HTS_FMT_CSI, HTS_IDX_SILENT_FAIL)) {
            auto result = std::make_unique<VariantTraversalIndex>();
            result->index = index;
            return result;
        }
    }
    if (std::filesystem::exists(tbi_path)) {
        if (auto* tabix = tbx_index_load2(filename.c_str(), tbi_path.c_str())) {
            auto result = std::make_unique<VariantTraversalIndex>();
            result->index = tabix->idx;
            result->tabix = tabix;
            return result;
        }
    }
    return {};
}

using RawAfMap = std::unordered_map<LocusKey, double>;

bool is_text_vcf(const std::string& path) {
    return fastgatk::io::has_suffix_ci(path, ".vcf") ||
           fastgatk::io::has_suffix_ci(path, ".vcf.gz");
}

// HTSlib decodes Type=Float INFO values to binary32.  GATK/htsjdk parses the
// VCF token as a binary64 value, however, so retain the textual AF for VCF
// inputs before bcf_read() loses those digits (e.g. 0.123456789).  This map is
// only used on sequential/CSI paths; tabix traversal parses the same line
// inline and therefore avoids a second scan of large indexed resources.
RawAfMap load_raw_vcf_af(const std::string& filename,
                         const std::vector<std::string>& bam_contigs) {
    RawAfMap values;
    if (!is_text_vcf(filename)) return values;
    std::unordered_map<std::string, int> tids;
    for (std::size_t tid = 0; tid < bam_contigs.size(); ++tid)
        tids.emplace(bam_contigs[tid], static_cast<int>(tid));
    htsFile* raw = hts_open(filename.c_str(), "r");
    if (!raw) return values;
    kstring_t line{0, 0, nullptr};
    while (fastgatk::io::read_text_line(raw, &line, filename) >= 0) {
        if (line.l == 0 || line.s[0] == '#') continue;
        std::string text(line.s, line.l);
        std::size_t field_start = 0;
        std::string contig;
        std::string info;
        int position = -1;
        for (int field = 0; field <= 7; ++field) {
            const auto field_end = text.find('\t', field_start);
            if (field_end == std::string::npos && field < 7) break;
            const auto value = text.substr(field_start,
                                           field_end == std::string::npos
                                               ? std::string::npos : field_end - field_start);
            if (field == 0) contig = value;
            else if (field == 1) {
                char* end = nullptr;
                const auto parsed = std::strtol(value.c_str(), &end, 10);
                if (end != value.c_str() && end != nullptr && *end == '\0' && parsed > 0 &&
                    parsed <= std::numeric_limits<int>::max())
                    position = static_cast<int>(parsed - 1);
            } else if (field == 7) info = value;
            if (field_end == std::string::npos) break;
            field_start = field_end + 1;
        }
        const auto tid = tids.find(contig);
        if (tid == tids.end() || position < 0 || info.empty()) continue;
        double af = std::numeric_limits<double>::quiet_NaN();
        std::size_t token_start = 0;
        while (token_start <= info.size()) {
            const auto token_end = info.find(';', token_start);
            const auto token = info.substr(token_start,
                                           token_end == std::string::npos
                                               ? std::string::npos : token_end - token_start);
            if (token.rfind("AF=", 0) == 0) {
                const auto first = token.substr(3, token.find(',', 3) - 3);
                char* end = nullptr;
                const auto parsed = std::strtod(first.c_str(), &end);
                if (!first.empty() && end != first.c_str() && end != nullptr && *end == '\0' &&
                    std::isfinite(parsed)) af = parsed;
                break;
            }
            if (token_end == std::string::npos) break;
            token_start = token_end + 1;
        }
        if (std::isfinite(af)) values.emplace(locus_key(tid->second, position), af);
    }
    free(line.s);
    hts_close(raw);
    return values;
}

void capture_raw_vcf_af(const char* text, std::size_t length,
                        const std::unordered_map<std::string, int>& bam_tids,
                        RawAfMap& values) {
    if (text == nullptr || length == 0) return;
    std::string line(text, length);
    std::size_t field_start = 0;
    std::string contig;
    std::string info;
    int position = -1;
    for (int field = 0; field <= 7; ++field) {
        const auto field_end = line.find('\t', field_start);
        if (field_end == std::string::npos && field < 7) return;
        const auto value = line.substr(field_start,
                                       field_end == std::string::npos
                                           ? std::string::npos : field_end - field_start);
        if (field == 0) contig = value;
        else if (field == 1) {
            char* end = nullptr;
            const auto parsed = std::strtol(value.c_str(), &end, 10);
            if (end == value.c_str() || end == nullptr || *end != '\0' || parsed <= 0 ||
                parsed > std::numeric_limits<int>::max()) return;
            position = static_cast<int>(parsed - 1);
        } else if (field == 7) info = value;
        if (field_end == std::string::npos) break;
        field_start = field_end + 1;
    }
    const auto tid = bam_tids.find(contig);
    if (tid == bam_tids.end() || position < 0) return;
    std::size_t token_start = 0;
    while (token_start <= info.size()) {
        const auto token_end = info.find(';', token_start);
        const auto token = info.substr(token_start,
                                       token_end == std::string::npos
                                           ? std::string::npos : token_end - token_start);
        if (token.rfind("AF=", 0) == 0) {
            const auto first = token.substr(3, token.find(',', 3) - 3);
            char* end = nullptr;
            const auto parsed = std::strtod(first.c_str(), &end);
            if (!first.empty() && end != first.c_str() && end != nullptr && *end == '\0' &&
                std::isfinite(parsed))
                values.emplace(locus_key(tid->second, position), parsed);
            return;
        }
        if (token_end == std::string::npos) break;
        token_start = token_end + 1;
    }
}

std::vector<Site> load_sites(const Options& options,
                             const fastgatk::io::HeaderSummary& bam_header,
                             const std::vector<fastgatk::io::IndexedInterval>& intervals,
                             const std::vector<fastgatk::io::IndexedInterval>& exclusions,
                             std::size_t& input_records, std::size_t& selected_records,
                             std::size_t& missing_af_records,
                             bool& indexed_variant_traversal,
                             std::size_t& variant_interval_queries,
                             RawAfMap& raw_af_values) {
    std::vector<Site> sites;
    std::unordered_map<std::string, int> bam_tids;
    for (std::size_t tid = 0; tid < bam_header.contigs.size(); ++tid)
        bam_tids.emplace(bam_header.contigs[tid], static_cast<int>(tid));
    bool interval_key_filter_active = false;
    const auto interval_keys = load_interval_variant_keys(
        options.regions, options.interval_set_rule, interval_key_filter_active);

    htsFile* file = bcf_open(options.variants.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open variants: " + options.variants);
    bcf_hdr_t* header = bcf_hdr_read(file);
    if (!header) {
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read variants header: " + options.variants);
    }
    validate_variant_dictionary(header, bam_header, options.disable_sequence_dictionary_validation);
    const int af_id = bcf_hdr_id2int(header, BCF_DT_ID, "AF");
    if (!bcf_hdr_idinfo_exists(header, BCF_HL_INFO, af_id)) {
        bcf_hdr_destroy(header);
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: population VCF does not have an AF INFO field");
    }
    bcf1_t* record = bcf_init();
    if (!record) {
        bcf_hdr_destroy(header);
        bcf_close(file);
        throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
    }
    std::unordered_set<LocusKey> seen;
    const auto process_record = [&]() {
        ++input_records;
        bcf_unpack(record, BCF_UN_ALL);
        const auto* contig_ptr = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
        if (!contig_ptr || record->n_allele != 2 || !record->d.allele[0] || !record->d.allele[1]) {
            bcf_clear(record);
            return;
        }
        const std::string contig(contig_ptr);
        const auto bam_tid = bam_tids.find(contig);
        if (bam_tid == bam_tids.end()) {
            bcf_clear(record);
            return;
        }
        if (interval_key_filter_active &&
            interval_keys.find(contig + ":" + std::to_string(record->pos)) == interval_keys.end()) {
            bcf_clear(record);
            return;
        }
        if (!intervals.empty() && std::none_of(intervals.begin(), intervals.end(), [&](const auto& interval) {
                return record->rid == interval.rid && record->pos >= interval.begin && record->pos < interval.end;
            })) {
            bcf_clear(record);
            return;
        }
        if (!exclusions.empty() && std::any_of(exclusions.begin(), exclusions.end(), [&](const auto& interval) {
                return record->rid == interval.rid && record->pos >= interval.begin && record->pos < interval.end;
            })) {
            bcf_clear(record);
            return;
        }
        const std::string ref(record->d.allele[0]);
        const std::string alt(record->d.allele[1]);
        if (ref.size() != 1 || alt.size() != 1) {
            bcf_clear(record);
            return;
        }
        float* af_values = nullptr;
        int af_count = 0;
        const int af_length = bcf_get_info_float(header, record, "AF", &af_values, &af_count);
        if (af_length <= 0 || af_count <= 0 || !af_values || !std::isfinite(af_values[0])) {
            ++missing_af_records;
            free(af_values);
            bcf_clear(record);
            return;
        }
        const float af_value = af_values[0];
        const auto raw_key = locus_key(bam_tid->second, static_cast<int>(record->pos));
        const auto raw_af = raw_af_values.find(raw_key);
        const double af = raw_af == raw_af_values.end()
            ? static_cast<double>(af_value) : raw_af->second;
        free(af_values);
        if (!(options.min_af < af && af < options.max_af)) {
            bcf_clear(record);
            return;
        }
        const int ref_code = base_code(ref[0]);
        const int alt_code = base_code(alt[0]);
        if (ref_code < 0 || alt_code < 0 || ref_code == alt_code) {
            bcf_clear(record);
            return;
        }
        const LocusKey key = locus_key(bam_tid->second, record->pos);
        if (!seen.insert(key).second) {
            bcf_clear(record);
            return;
        }
        if (record->pos > std::numeric_limits<int>::max()) {
            bcf_clear(record);
            return;
        }
        sites.push_back(Site{contig, ref, alt, bam_tid->second, static_cast<int>(record->pos),
                             af, ref_code, alt_code});
        ++selected_records;
        bcf_clear(record);
    };

    // A large population VCF is the dominant I/O cost of GetPileupSummaries.
    // When -L produced coordinate intervals and a tabix/CSI index is
    // available, traverse only those intervals.  VCF-file selectors still
    // use the same exact-key filter, while the indexed path avoids decoding
    // unrelated records.  If no index exists, retain the deterministic
    // sequential fallback.
    auto variant_index = intervals.empty()
        ? std::unique_ptr<VariantTraversalIndex>{}
        : load_variant_index(options.variants);
    indexed_variant_traversal = variant_index != nullptr;
    if (variant_index == nullptr || variant_index->tabix == nullptr)
        raw_af_values = load_raw_vcf_af(options.variants, bam_header.contigs);
    if (variant_index != nullptr) {
        for (const auto& interval : intervals) {
            ++variant_interval_queries;
            const int index_tid = variant_index->tid_for(header, interval.rid);
            if (index_tid < 0) continue;
            int status = 0;
            if (variant_index->tabix != nullptr) {
                // TBI sequence IDs are not required to match the full VCF
                // header's rid values.  Use tbx's name-mapped iterator and
                // parse the selected VCF lines through the public vcf_parse
                // API so the normal AF/allele/filter path remains shared.
                hts_itr_t* iterator = tbx_itr_queryi(variant_index->tabix, index_tid,
                                                     interval.begin, interval.end);
                if (iterator == nullptr) continue;
                kstring_t line{0, 0, nullptr};
                while ((status = tbx_itr_next(file, variant_index->tabix, iterator, &line)) >= 0) {
                    if (vcf_parse(&line, header, record) < 0) {
                        free(line.s);
                        tbx_itr_destroy(iterator);
                        bcf_destroy(record);
                        bcf_hdr_destroy(header);
                        bcf_close(file);
                        throw std::runtime_error("BAD_INPUT: indexed population VCF record parse failed");
                    }
                    capture_raw_vcf_af(line.s, line.l, bam_tids, raw_af_values);
                    process_record();
                }
                free(line.s);
                tbx_itr_destroy(iterator);
            } else {
                hts_itr_t* iterator = bcf_itr_queryi(variant_index->index, index_tid,
                                                     interval.begin, interval.end);
                if (iterator == nullptr) continue;
                while ((status = bcf_itr_next(file, iterator, record)) >= 0)
                    process_record();
                hts_itr_destroy(iterator);
            }
            if (status < -1) {
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(file);
                throw std::runtime_error("BAD_INPUT: indexed population VCF traversal failed");
            }
        }
    } else {
        while (fastgatk::io::read_variant_record(file, header, record,
                                               options.variants) == 0)
            process_record();
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(file);
    if (sites.empty()) {
        if (input_records != 0 && missing_af_records == input_records)
            throw std::runtime_error("BAD_INPUT: no variants in population VCF had an AF field");
        // GATK's LocusWalker succeeds when every AF-bearing record is outside
        // the requested range/interval (or is not a biallelic SNP).  Keep the
        // normal metadata + six-column header and emit zero data rows; only a
        // population resource with no usable AF at all is a bad input.
    }
    return sites;
}

void normalize_pileup_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<fastgatk::io::IndexedInterval> merged;
    merged.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (interval.end <= interval.begin) continue;
        if (!merged.empty() && merged.back().rid == interval.rid && interval.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    intervals.swap(merged);
}

void intersect_pileup_intervals(
    std::vector<fastgatk::io::IndexedInterval>& current,
    const std::vector<fastgatk::io::IndexedInterval>& incoming) {
    std::vector<fastgatk::io::IndexedInterval> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < current.size() && right < incoming.size()) {
        if (current[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < current[left].rid) { ++right; continue; }
        const int begin = std::max(current[left].begin, incoming[right].begin);
        const int end = std::min(current[left].end, incoming[right].end);
        if (begin < end) {
            auto overlap = current[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (current[left].end < incoming[right].end) ++left;
        else ++right;
    }
    current.swap(intersection);
}

std::vector<fastgatk::io::IndexedInterval> load_intervals_for_sites(
    const std::vector<std::string>& regions,
    const std::string& variants,
    fastgatk::io::HtsIntervalSetRule rule,
    int padding = 0) {
    std::vector<fastgatk::io::IndexedInterval> intervals;
    htsFile* file = bcf_open(variants.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open variants: " + variants);
    bcf_hdr_t* header = bcf_hdr_read(file);
    if (!header) {
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read variants header: " + variants);
    }
    fastgatk::io::IntervalFileStats stats;
    bool first_selector = true;
    for (const auto& selector : regions) {
        std::vector<fastgatk::io::IndexedInterval> incoming;
        if (!is_variant_file(selector)) {
            fastgatk::io::append_interval_selector(selector, header, incoming, stats);
        } else {
            // A VCF supplied as -L is an exact-locus selector.  Materialize its
            // loci as one-base indexed queries against the population resource.
            std::unordered_set<LocusKey> seen_variant_selector_loci;
            htsFile* selector_file = bcf_open(selector.c_str(), "r");
            if (!selector_file)
                throw std::runtime_error("BAD_INPUT: cannot open interval VCF: " + selector);
            bcf_hdr_t* selector_header = bcf_hdr_read(selector_file);
            bcf1_t* selector_record = selector_header ? bcf_init() : nullptr;
            if (!selector_header || !selector_record) {
                if (selector_record) bcf_destroy(selector_record);
                if (selector_header) bcf_hdr_destroy(selector_header);
                bcf_close(selector_file);
                throw std::runtime_error("BAD_INPUT: cannot read interval VCF: " + selector);
            }
            while (fastgatk::io::read_variant_record(
                       selector_file, selector_header, selector_record, variants) == 0) {
                const auto* contig = selector_record->rid >= 0
                    ? bcf_hdr_id2name(selector_header, selector_record->rid) : nullptr;
                const int rid = contig == nullptr ? -1 : bcf_hdr_name2id(header, contig);
                if (rid >= 0 && selector_record->pos >= 0 &&
                    selector_record->pos <= std::numeric_limits<int>::max() - 1 &&
                    seen_variant_selector_loci.insert(locus_key(rid, selector_record->pos)).second) {
                    const int pos = static_cast<int>(selector_record->pos);
                    const auto span_value = std::max<hts_pos_t>(1, selector_record->rlen);
                    const int span = span_value > std::numeric_limits<int>::max() - pos
                        ? std::numeric_limits<int>::max() - pos
                        : static_cast<int>(span_value);
                    incoming.push_back(fastgatk::io::IndexedInterval{rid, pos, pos + span, {}});
                }
                bcf_clear(selector_record);
            }
            bcf_destroy(selector_record);
            bcf_hdr_destroy(selector_header);
            bcf_close(selector_file);
        }
        if (padding != 0) {
            for (auto& interval : incoming) {
                interval.begin = std::max(0, interval.begin - padding);
                interval.end = interval.end > std::numeric_limits<int>::max() - padding
                    ? std::numeric_limits<int>::max() : interval.end + padding;
            }
        }
        normalize_pileup_intervals(incoming);
        if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector)
            intervals.insert(intervals.end(), incoming.begin(), incoming.end());
        else
            intersect_pileup_intervals(intervals, incoming);
        first_selector = false;
    }
    bcf_hdr_destroy(header);
    bcf_close(file);
    normalize_pileup_intervals(intervals);
    return intervals;
}

bool passes_named_read_filter(const fastgatk::io::ReadBatch& batch, std::size_t record,
                              const Options& options, ReadFilterKind filter) {
    if (record >= batch.records() || record >= batch.tids.size()) return false;
    const auto flags = batch.flags.size() == batch.records() ? batch.flags[record] : 0U;
    const auto mapq = batch.mapq.size() == batch.records() ? batch.mapq[record] : 255U;
    switch (filter) {
    case ReadFilterKind::ReadLength: {
        if (batch.offsets.size() != batch.records() + 1 ||
            batch.offsets[record + 1] < batch.offsets[record]) return false;
        const auto length = batch.offsets[record + 1] - batch.offsets[record];
        return length >= options.min_read_length && length <= options.max_read_length;
    }
    case ReadFilterKind::MappingQualityAvailable: return mapq != 255;
    case ReadFilterKind::MappingQualityNotZero: return mapq != 0;
    case ReadFilterKind::MappingQuality:
        return mapq >= static_cast<unsigned int>(options.minimum_mapping_quality);
    case ReadFilterKind::Mapped: return (flags & BAM_FUNMAP) == 0;
    case ReadFilterKind::PrimaryLine: return (flags & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
    case ReadFilterKind::NotSecondary: return (flags & BAM_FSECONDARY) == 0;
    case ReadFilterKind::NotSupplementary: return (flags & BAM_FSUPPLEMENTARY) == 0;
    case ReadFilterKind::NotDuplicate: return (flags & BAM_FDUP) == 0;
    case ReadFilterKind::PassesVendorQuality: return (flags & BAM_FQCFAIL) == 0;
    case ReadFilterKind::MateOnSameContig:
        if (!(flags & BAM_FPAIRED) || (flags & BAM_FMUNMAP)) return true;
        if (batch.mate_tids.size() != batch.records() || batch.mate_tids[record] != batch.tids[record])
            return false;
        return true;
    case ReadFilterKind::NonZeroReferenceLength:
        return batch.positions[record] >= 0 &&
               fastgatk::io::reference_end(batch, record) > batch.positions[record];
    case ReadFilterKind::GoodCigar: return batch.cigar_record_layout_valid(record);
    case ReadFilterKind::Wellformed:
        if (batch.offsets.size() != batch.records() + 1 || batch.offsets[record] > batch.offsets[record + 1] ||
            batch.offsets.back() != batch.bases.size() || batch.qualities.size() != batch.bases.size())
            return false;
        return true;
    }
    return true;
}

bool passes_read_filter(const fastgatk::io::ReadBatch& batch, std::size_t record,
                        const Options& options, const ReadFilterConfig& config) {
    if (record >= batch.records() || record >= batch.tids.size()) return false;
    const auto enabled = [&](ReadFilterKind filter, bool active) {
        return !active || passes_named_read_filter(batch, record, options, filter);
    };
    if (!enabled(ReadFilterKind::ReadLength, config.read_length) ||
        !enabled(ReadFilterKind::MappingQualityAvailable, config.mapping_quality_available) ||
        !enabled(ReadFilterKind::MappingQualityNotZero, config.mapping_quality_not_zero) ||
        !enabled(ReadFilterKind::MappingQuality, config.mapping_quality) ||
        !enabled(ReadFilterKind::Mapped, config.mapped) ||
        !enabled(ReadFilterKind::NotSecondary, config.secondary) ||
        !enabled(ReadFilterKind::NotSupplementary, config.supplementary) ||
        !enabled(ReadFilterKind::NotDuplicate, config.duplicate) ||
        !enabled(ReadFilterKind::PassesVendorQuality, config.vendor_quality) ||
        !enabled(ReadFilterKind::MateOnSameContig, config.mate_same_contig) ||
        !enabled(ReadFilterKind::NonZeroReferenceLength, config.nonzero_reference) ||
        !enabled(ReadFilterKind::GoodCigar, config.good_cigar) ||
        !enabled(ReadFilterKind::Wellformed, config.wellformed)) return false;
    for (const auto& name : options.inverted_read_filters)
        if (passes_named_read_filter(batch, record, options, read_filter_kind(name))) return false;
    return true;
}

struct CountStats {
    std::uint64_t reads_seen = 0;
    std::uint64_t reads_used = 0;
    std::uint64_t bases_used = 0;
    std::size_t count_batches = 0;
    std::size_t count_observations = 0;
    double count_prepare_seconds = 0.0;
    double count_seconds = 0.0;
    std::string count_execution_space;
    std::size_t count_buffer_capacity = 0;
    std::size_t count_buffer_allocations = 0;
    std::size_t count_buffer_reuses = 0;
};

// Host-owned result emitted by one decode/compute stage.  Keeping interval
// counts in the queued value lets the sink merge them deterministically while
// the next HTSlib batch is already being decoded.
struct PileupBatchResult {
    std::vector<std::array<std::uint64_t, 3>> counts;
    CountStats stats;
};

template<class ExecSpace>
class PileupCountBuffers {
public:
    Kokkos::View<int*> site_indices;
    Kokkos::View<int*> observed_bases;
    Kokkos::View<int*> refs;
    Kokkos::View<int*> alts;
    Kokkos::View<std::uint64_t**> device_counts;
    Kokkos::View<int*, Kokkos::HostSpace> host_refs;
    Kokkos::View<int*, Kokkos::HostSpace> host_alts;
    Kokkos::View<int*, Kokkos::HostSpace> host_site_indices;
    Kokkos::View<int*, Kokkos::HostSpace> host_observed_bases;
    std::size_t capacity = 0;
    std::size_t site_capacity = 0;
    std::size_t allocations = 0;
    std::size_t reuses = 0;

    void ensure(std::size_t observations, std::size_t sites) {
        if (sites != site_capacity) {
            refs = Kokkos::View<int*>("pileup_refs_persistent", sites);
            alts = Kokkos::View<int*>("pileup_alts_persistent", sites);
            device_counts = Kokkos::View<std::uint64_t**>("pileup_counts_persistent", sites, 3);
            host_refs = Kokkos::View<int*, Kokkos::HostSpace>("pileup_refs_host", sites);
            host_alts = Kokkos::View<int*, Kokkos::HostSpace>("pileup_alts_host", sites);
            site_capacity = sites;
            ++allocations;
        }
        if (observations <= capacity) {
            ++reuses;
            return;
        }
        const auto doubled = capacity > (std::numeric_limits<std::size_t>::max() / 2)
            ? std::numeric_limits<std::size_t>::max() : capacity * 2;
        const auto next = std::max(observations, capacity == 0 ? std::size_t{1} : doubled);
        site_indices = Kokkos::View<int*>("pileup_site_indices_persistent", next);
        observed_bases = Kokkos::View<int*>("pileup_observed_bases_persistent", next);
        host_site_indices = Kokkos::View<int*, Kokkos::HostSpace>("pileup_site_indices_host", next);
        host_observed_bases = Kokkos::View<int*, Kokkos::HostSpace>("pileup_observed_bases_host", next);
        capacity = next;
        ++allocations;
    }

    void clear() {
        site_indices = Kokkos::View<int*>{};
        observed_bases = Kokkos::View<int*>{};
        refs = Kokkos::View<int*>{};
        alts = Kokkos::View<int*>{};
        device_counts = Kokkos::View<std::uint64_t**>{};
        host_refs = Kokkos::View<int*, Kokkos::HostSpace>{};
        host_alts = Kokkos::View<int*, Kokkos::HostSpace>{};
        host_site_indices = Kokkos::View<int*, Kokkos::HostSpace>{};
        host_observed_bases = Kokkos::View<int*, Kokkos::HostSpace>{};
        capacity = 0;
        site_capacity = 0;
    }
};

PileupBatchResult count_batch(const fastgatk::io::ReadBatch& batch, const Options& options,
                              const ReadFilterConfig& filter_config,
                              const std::unordered_map<LocusKey, int>& site_lookup,
                              const std::vector<Site>& sites,
                              PileupCountBuffers<Kokkos::DefaultExecutionSpace>& buffers) {
    PileupBatchResult result;
    result.counts.assign(sites.size(), {0, 0, 0});
    auto& stats = result.stats;
    std::vector<int> observation_sites;
    std::vector<int> observation_bases;
    for (std::size_t record = 0; record < batch.records(); ++record) {
        ++stats.reads_seen;
        if (!passes_read_filter(batch, record, options, filter_config)) continue;
        if (batch.offsets.size() != batch.records() + 1 || batch.cigar_offsets.size() != batch.records() + 1 ||
            batch.positions[record] < 0 || batch.tids[record] < 0) continue;
        const auto read_begin = batch.offsets[record];
        const auto read_end = batch.offsets[record + 1];
        const auto cigar_begin = batch.cigar_offsets[record];
        const auto cigar_end = batch.cigar_offsets[record + 1];
        std::uint64_t read_cursor = 0;
        std::int64_t reference_position = batch.positions[record];
        bool valid = true;
        std::size_t read_observations = 0;
        for (std::size_t cigar_index = cigar_begin; cigar_index < cigar_end; ++cigar_index) {
            const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[cigar_index]);
            if (!operation.valid()) { valid = false; break; }
            if (operation.projects_base()) {
                if (read_cursor + operation.length > read_end - read_begin) { valid = false; break; }
                for (std::uint32_t offset = 0; offset < operation.length; ++offset) {
                    const auto key = locus_key(batch.tids[record], static_cast<int>(reference_position + offset));
                    const auto site = site_lookup.find(key);
                    if (site == site_lookup.end()) continue;
                    const auto base = base_code(static_cast<char>(batch.bases[read_begin + read_cursor + offset]));
                    if (base < 0) continue;
                    observation_sites.push_back(site->second);
                    observation_bases.push_back(base);
                    ++read_observations;
                }
                read_cursor += operation.length;
                reference_position += operation.length;
            } else {
                if (operation.consumes_read()) read_cursor += operation.length;
                if (operation.consumes_reference()) reference_position += operation.length;
                if (read_cursor > read_end - read_begin) { valid = false; break; }
            }
        }
        if (!valid || read_cursor != read_end - read_begin) continue;
        if (read_observations != 0) ++stats.reads_used;
        stats.bases_used += read_observations;
    }
    if (observation_sites.empty()) {
        stats.count_buffer_capacity = buffers.capacity;
        stats.count_buffer_allocations = buffers.allocations;
        stats.count_buffer_reuses = buffers.reuses;
        return result;
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch count_host("get-pileup-summaries-count-v1");
    count_host.records = observation_sites.size();
    count_host.bytes = observation_sites.size() * sizeof(int) +
                       observation_bases.size() * sizeof(int) +
                       sites.size() * 2 * sizeof(int) +
                       sites.size() * 3 * sizeof(std::uint64_t);
    fastgatk::core::KernelPlan<ExecSpace> count_plan("get-pileup-summaries-count");
    count_plan.begin_prepare(count_host);
    buffers.ensure(observation_sites.size(), sites.size());
    auto site_indices = Kokkos::subview(buffers.site_indices,
        std::make_pair<std::size_t, std::size_t>(0, observation_sites.size()));
    auto observed_bases = Kokkos::subview(buffers.observed_bases,
        std::make_pair<std::size_t, std::size_t>(0, observation_bases.size()));
    auto refs = buffers.refs;
    auto alts = buffers.alts;
    auto device_counts = buffers.device_counts;
    for (std::size_t index = 0; index < sites.size(); ++index) {
        buffers.host_refs(index) = sites[index].ref_code;
        buffers.host_alts(index) = sites[index].alt_code;
    }
    Kokkos::deep_copy(refs, buffers.host_refs);
    Kokkos::deep_copy(alts, buffers.host_alts);
    for (std::size_t index = 0; index < observation_sites.size(); ++index) {
        buffers.host_site_indices(index) = observation_sites[index];
        buffers.host_observed_bases(index) = observation_bases[index];
    }
    Kokkos::deep_copy(site_indices, Kokkos::subview(buffers.host_site_indices,
        std::make_pair<std::size_t, std::size_t>(0, observation_sites.size())));
    Kokkos::deep_copy(observed_bases, Kokkos::subview(buffers.host_observed_bases,
        std::make_pair<std::size_t, std::size_t>(0, observation_bases.size())));
    fastgatk::core::DeviceBatch<ExecSpace> count_device(observation_sites.size());
    count_device.bind("site_indices", buffers.site_indices, observation_sites.size());
    count_device.bind("observed_bases", buffers.observed_bases, observation_bases.size());
    count_device.bind("refs", refs);
    count_device.bind("alts", alts);
    count_device.bind("counts", device_counts);
    ExecSpace().fence();
    count_plan.end_prepare(count_device);
    count_plan.begin_execute();
    Kokkos::deep_copy(device_counts, std::uint64_t{0});
    Kokkos::parallel_for("get_pileup_summaries_count", Kokkos::RangePolicy<ExecSpace>(0, observation_sites.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            const int site = site_indices(index);
            const int base = observed_bases(index);
            const int category = base == refs(site) ? 0 : base == alts(site) ? 1 : 2;
            Kokkos::atomic_fetch_add(&device_counts(site, category), std::uint64_t{1});
        });
    Kokkos::fence();
    count_plan.end_execute();
    ++stats.count_batches;
    stats.count_observations += observation_sites.size();
    stats.count_prepare_seconds += count_plan.telemetry().prepare_seconds;
    stats.count_seconds += count_plan.telemetry().execute_seconds;
    stats.count_execution_space = ExecSpace::name();
    stats.count_buffer_capacity = buffers.capacity;
    stats.count_buffer_allocations = buffers.allocations;
    stats.count_buffer_reuses = buffers.reuses;
    auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_counts);
    for (std::size_t site = 0; site < sites.size(); ++site)
        for (int category = 0; category < 3; ++category)
            result.counts[site][category] = host_counts(site, category);
    stats.count_buffer_capacity = buffers.capacity;
    stats.count_buffer_allocations = buffers.allocations;
    stats.count_buffer_reuses = buffers.reuses;
    return result;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
    Kokkos::initialize(settings);
    bool initialized = true;
    PileupCountBuffers<Kokkos::DefaultExecutionSpace> count_buffers;
    try {
        const auto interval_selectors = [&]() {
            std::vector<std::string> selectors;
            for (const auto& selector : options.regions)
                if (!is_variant_file(selector)) selectors.push_back(selector);
            return selectors;
        }();
        const fastgatk::runtime::AdaptiveController controller;
        const auto initial_limits = controller.initial(
            resources.budget(), fastgatk::runtime::WorkEstimate{
                static_cast<std::uint32_t>(std::min<std::size_t>(
                    options.batch_records, std::numeric_limits<std::uint32_t>::max())),
                1, 2048, 0, 2048});
        const auto initial_batch_records = std::min<std::size_t>(
            options.batch_records, initial_limits.max_reads);
        fastgatk::io::HtsReader reader(options.input, options.reference, interval_selectors,
                                       initial_batch_records);
        const auto intervals = load_intervals_for_sites(
            options.regions, options.variants, options.interval_set_rule, options.interval_padding);
        const auto exclusions = load_intervals_for_sites(
            options.excluded_regions, options.variants,
            fastgatk::io::HtsIntervalSetRule::Union, options.exclusion_padding);
        std::size_t input_sites = 0, selected_sites = 0;
        std::size_t missing_af_records = 0;
        bool indexed_variant_traversal = false;
        std::size_t variant_interval_queries = 0;
        const auto filter_config = effective_read_filters(options);
        RawAfMap raw_af_values;
        auto sites = load_sites(options, reader.header(), intervals, exclusions, input_sites, selected_sites,
                                missing_af_records, indexed_variant_traversal,
                                variant_interval_queries, raw_af_values);
        std::unordered_map<LocusKey, int> site_lookup;
        site_lookup.reserve(sites.size() * 2 + 1);
        for (std::size_t index = 0; index < sites.size(); ++index)
            site_lookup.emplace(locus_key(sites[index].tid, sites[index].position), static_cast<int>(index));
        std::vector<std::array<std::uint64_t, 3>> counts(sites.size(), {0, 0, 0});
        CountStats stats;
        using ReadBatch = fastgatk::io::ReadBatch;
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<ReadBatch, PileupBatchResult, PileupBatchResult>;
        auto batch_limits = initial_limits;
        batch_limits.max_reads = static_cast<std::uint32_t>(std::min<std::size_t>(
            reader.batch_records(), std::numeric_limits<std::uint32_t>::max()));
        std::uint64_t adaptive_reductions = 0;
        const auto saturating_multiply = [](std::uint64_t left, std::uint64_t right) {
            return right != 0 && left > std::numeric_limits<std::uint64_t>::max() / right
                ? std::numeric_limits<std::uint64_t>::max() : left * right;
        };
        const auto saturating_add = [](std::uint64_t left, std::uint64_t right) {
            return right > std::numeric_limits<std::uint64_t>::max() - left
                ? std::numeric_limits<std::uint64_t>::max() : left + right;
        };
        const auto safe_budget = resources.safe_memory_budget_bytes();
        const auto decoded_capacity = safe_budget != 0
            ? safe_budget
            : std::max<std::uint64_t>(64ULL << 20,
                saturating_multiply(static_cast<std::uint64_t>(initial_batch_records), 16384));
        const auto result_capacity = std::max<std::uint64_t>(sizeof(PileupBatchResult),
            saturating_add(static_cast<std::uint64_t>(sizeof(PileupBatchResult)),
                           saturating_multiply(static_cast<std::uint64_t>(sites.size()),
                                               3U * sizeof(std::uint64_t))));
        const auto pipeline_limits = Pipeline::Limits{
            decoded_capacity, result_capacity, result_capacity};
        const auto pipeline_metrics = Pipeline(
            pipeline_limits,
            [&]() -> std::optional<ReadBatch> {
                ReadBatch batch;
                if (!reader.next(batch)) return std::nullopt;
                fastgatk::runtime::RuntimeTelemetry runtime;
                runtime.host_bytes = batch.bytes();
                runtime.inflight_bytes = runtime.host_bytes;
                runtime.compute_queue_empty = true;
                const auto pressure = safe_budget != 0 && runtime.host_bytes > safe_budget
                    ? fastgatk::runtime::Pressure::HostMemory
                    : fastgatk::runtime::Pressure::Normal;
                const auto next_limits = controller.next(runtime, batch_limits, pressure);
                if (next_limits.max_reads < reader.batch_records()) {
                    reader.set_batch_records(next_limits.max_reads);
                    ++adaptive_reductions;
                }
                batch_limits = next_limits;
                return batch;
            },
            [&](ReadBatch batch) -> std::optional<PileupBatchResult> {
                return count_batch(batch, options, filter_config, site_lookup, sites, count_buffers);
            },
            [](PileupBatchResult result) -> std::optional<PileupBatchResult> {
                return result;
            },
            [&](PileupBatchResult result) {
                stats.reads_seen += result.stats.reads_seen;
                stats.reads_used += result.stats.reads_used;
                stats.bases_used += result.stats.bases_used;
                stats.count_batches += result.stats.count_batches;
                stats.count_observations += result.stats.count_observations;
                stats.count_prepare_seconds += result.stats.count_prepare_seconds;
                stats.count_seconds += result.stats.count_seconds;
                if (!result.stats.count_execution_space.empty())
                    stats.count_execution_space = result.stats.count_execution_space;
                stats.count_buffer_capacity = result.stats.count_buffer_capacity;
                stats.count_buffer_allocations = result.stats.count_buffer_allocations;
                stats.count_buffer_reuses = result.stats.count_buffer_reuses;
                for (std::size_t site = 0; site < sites.size(); ++site)
                    for (int category = 0; category < 3; ++category)
                        counts[site][category] += result.counts[site][category];
            },
            [](const ReadBatch& batch) { return batch.bytes(); },
            [&](const PileupBatchResult& result) {
                return saturating_add(static_cast<std::uint64_t>(sizeof(PileupBatchResult)),
                                      saturating_multiply(static_cast<std::uint64_t>(result.counts.size()),
                                                          3U * sizeof(std::uint64_t)));
            },
            [&](const PileupBatchResult& result) {
                return saturating_add(static_cast<std::uint64_t>(sizeof(PileupBatchResult)),
                                      saturating_multiply(static_cast<std::uint64_t>(result.counts.size()),
                                                          3U * sizeof(std::uint64_t)));
            }).run();

        const auto sample = bam_sample_name(options.input, options.reference);
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open output: " + options.output);
        output << "#<METADATA>SAMPLE=" << sample << '\n'
               << "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n";
        for (std::size_t index = 0; index < sites.size(); ++index) {
            const auto& site = sites[index];
            output << site.contig << '\t' << site.position + 1 << '\t'
                   << counts[index][0] << '\t' << counts[index][1] << '\t' << counts[index][2]
                   << '\t' << format_double(site.allele_frequency) << '\n';
        }
        output.flush();
        if (!output || !file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: pileup table is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"GetPileupSummaries\","
                 << "\"implementation\":\"fastgatk-get-pileup-summaries\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"table\","
                 << "\"compatibility\":{\"af_site_selection\":true,\"read_filter_mask\":true,"
                 << "\"inverted_read_filter\":true,"
                 << "\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\","
                 << "\"kokkos_count_kernel\":true,\"bit_identical_to_gatk\":false,"
                 << "\"full_pileup_engine\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"table\",\"complete\":true}],"
                 << "\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_sites\":" << input_sites << ",\"selected_sites\":" << selected_sites
                 << ",\"reads_seen\":" << stats.reads_seen << ",\"reads_used\":" << stats.reads_used
                 << ",\"bases_used\":" << stats.bases_used << ",\"sample\":\"" << json_escape(sample)
                 << "\",\"count_kernel_execution_space\":\""
                 << json_escape(stats.count_execution_space)
                 << "\",\"count_kernel_batches\":" << stats.count_batches
                 << ",\"count_kernel_observations\":" << stats.count_observations
                 << ",\"count_kernel_prepare_seconds\":" << stats.count_prepare_seconds
                 << ",\"count_kernel_seconds\":" << stats.count_seconds
                 << ",\"pipeline_lifecycle\":\"decode->compute->encode->sink\""
                 << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                 << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                 << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                 << ",\"pipeline_decoded_bytes\":" << pipeline_metrics.decoded_bytes
                 << ",\"pipeline_computed_bytes\":" << pipeline_metrics.computed_bytes
                 << ",\"pipeline_encoded_bytes\":" << pipeline_metrics.encoded_bytes
                 << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes
                 << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes
                 << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes
                 << ",\"persistent_buffer_capacity_records\":" << stats.count_buffer_capacity
                 << ",\"persistent_buffer_allocations\":" << stats.count_buffer_allocations
                 << ",\"persistent_buffer_reuses\":" << stats.count_buffer_reuses
                 << ",\"initial_batch_records\":" << initial_batch_records
                 << ",\"effective_batch_records\":" << reader.batch_records()
                 << ",\"adaptive_batch_reductions\":" << adaptive_reductions
                 << ",\"min_af\":" << options.min_af << ",\"max_af\":" << options.max_af
                 << ",\"minimum_mapping_quality\":" << options.minimum_mapping_quality
                 << ",\"inverted_read_filter_count\":" << options.inverted_read_filters.size()
                 << ",\"read_filter_require_read_length\":"
                 << (options.require_read_length ? "true" : "false")
                 << ",\"read_filter_min_read_length\":" << options.min_read_length
                 << ",\"read_filter_max_read_length\":" << options.max_read_length
                 << ",\"include_duplicates\":" << (options.include_duplicates ? "true" : "false")
                 << ",\"disable_tool_default_read_filters\":"
                 << (options.disable_tool_default_read_filters ? "true" : "false")
                 << ",\"disable_sequence_dictionary_validation\":"
                 << (options.disable_sequence_dictionary_validation ? "true" : "false")
                 << ",\"interval_padding\":" << options.interval_padding
                 << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                 << ",\"excluded_intervals\":" << exclusions.size()
                 << ",\"indexed_variant_traversal\":"
                 << (indexed_variant_traversal ? "true" : "false")
                 << ",\"variant_interval_queries\":" << variant_interval_queries
                 << ",\"missing_af_records\":" << missing_af_records
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"read_filters\":[";
        for (std::size_t index = 0; index < options.read_filters.size(); ++index) {
            if (index != 0) manifest << ',';
            manifest << '"' << json_escape(options.read_filters[index]) << '"';
        }
        manifest << "],\"disabled_read_filters\":[";
        for (std::size_t index = 0; index < options.disabled_read_filters.size(); ++index) {
            if (index != 0) manifest << ',';
            manifest << '"' << json_escape(options.disabled_read_filters[index]) << '"';
        }
        manifest << "]}}\n";
        std::cout << "{\"tool\":\"GetPileupSummaries\",\"status\":\"prototype\"," 
                  << "\"sites\":" << sites.size() << ",\"reads_used\":" << stats.reads_used
                  << ",\"bases_used\":" << stats.bases_used
                  << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"excluded_intervals\":" << exclusions.size()
                  << ",\"pipeline_lifecycle\":\"decode->compute->encode->sink\""
                  << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                  << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                  << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                  << "}\n";
        count_buffers.clear();
        Kokkos::finalize();
        initialized = false;
        return 0;
    } catch (...) {
        count_buffers.clear();
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        throw;
    }
}

#else

int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for GetPileupSummaries");
}

#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-get-pileup-summaries: " << error.what() << '\n';
        return 2;
    }
}
