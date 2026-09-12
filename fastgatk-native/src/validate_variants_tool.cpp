#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string manifest;
    std::string reference;
    std::string dbsnp;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> exclusions;
    int interval_padding = 0;
    int exclusion_padding = 0;
    std::vector<std::string> excluded_types;
    bool validate_gvcf = false;
    bool skip_filtered = false;
    bool warn_on_errors = false;
    bool fail_on_overlap = false;
    bool disable_sequence_dictionary_validation = false;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option_value(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* long_name,
                          const char* short_name = nullptr) {
    auto inline_value = option_value(argument, long_name);
    if (inline_value.empty() && short_name != nullptr)
        inline_value = option_value(argument, short_name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == long_name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-validate-variants (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE              input VCF/BCF\n"
                         "  -R, --reference FILE            reference FASTA for REF validation\n"
                         "      --dbsnp, -D FILE            dbSNP VCF for ID validation\n"
                         "  -L, --intervals REGION           interval selector (repeatable)\n"
                         "      --interval-set-rule, -isr RULE  UNION (default) or INTERSECTION\n"
                         "  -XL, --exclude-intervals REGION  exclusion selector (repeatable)\n"
                         "      -ip, --interval-padding N    include interval padding\n"
                         "      -ixp, --interval-exclusion-padding N  exclusion padding\n"
                         "      --validation-type-to-exclude TYPE (repeatable; FORMAT skips INFO/FORMAT cardinality)\n"
                         "      --do-not-validate-filtered-records[=BOOL]\n"
                         "      --warn-on-errors[=BOOL]     continue and report validation errors\n"
                         "      --validate-GVCF, --gvcf[=BOOL] validate GVCF blocks and coverage\n"
                         "      --fail-gvcf-on-overlap[=BOOL] reject overlapping GVCF blocks\n"
                         "      --disable-sequence-dictionary-validation[=BOOL]\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-V" || has_option_value(argument, "--variant"))
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        else if (argument == "-R" || has_option_value(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "--dbsnp" || argument == "-D" ||
                 has_option_value(argument, "--dbsnp") || has_option_value(argument, "-D"))
            options.dbsnp = require_value(index, argc, argv, argument, "--dbsnp", "-D");
        else if (argument == "-O" || has_option_value(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "-L" || has_option_value(argument, "--intervals") ||
                 has_option_value(argument, "--interval") || has_option_value(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                has_option_value(argument, "--region") ? "--region" :
                has_option_value(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (argument == "-isr" || has_option_value(argument, "--interval-set-rule") ||
                 has_option_value(argument, "-isr")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-XL" || has_option_value(argument, "--exclude-intervals"))
            options.exclusions.push_back(require_value(
                index, argc, argv, argument, "--exclude-intervals", "-XL"));
        else if (argument == "-ip" || has_option_value(argument, "--interval-padding")) {
            try {
                options.interval_padding = std::stoi(require_value(
                    index, argc, argv, argument, "--interval-padding", "-ip"));
            } catch (...) {
                throw std::invalid_argument("invalid integer for --interval-padding");
            }
            if (options.interval_padding < 0)
                throw std::invalid_argument("--interval-padding must be non-negative");
        } else if (argument == "-ixp" || has_option_value(argument, "--interval-exclusion-padding")) {
            try {
                options.exclusion_padding = std::stoi(require_value(
                    index, argc, argv, argument, "--interval-exclusion-padding", "-ixp"));
            } catch (...) {
                throw std::invalid_argument("invalid integer for --interval-exclusion-padding");
            }
            if (options.exclusion_padding < 0)
                throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
        }
        else if (argument == "--validation-type-to-exclude" || has_option_value(argument, "--validation-type-to-exclude") ||
                 argument == "-Xtype")
            options.excluded_types.push_back(require_value(
                index, argc, argv, argument, "--validation-type-to-exclude", "-Xtype"));
        else if (argument == "--do-not-validate-filtered-records" ||
                 argument == "-do-not-validate-filtered-records" ||
                 argument == "--do-not-validate-filtered" ||
                 argument.rfind("--do-not-validate-filtered-records=", 0) == 0 ||
                 argument.rfind("-do-not-validate-filtered-records=", 0) == 0 ||
                 argument.rfind("--do-not-validate-filtered=", 0) == 0)
            options.skip_filtered = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument,
                argument == "--do-not-validate-filtered" ||
                        argument.rfind("--do-not-validate-filtered=", 0) == 0
                    ? "--do-not-validate-filtered" :
                argument == "-do-not-validate-filtered-records" ||
                        argument.rfind("-do-not-validate-filtered-records=", 0) == 0
                    ? "-do-not-validate-filtered-records" : "--do-not-validate-filtered-records");
        else if (argument == "--warn-on-errors" || argument == "-warn-on-errors" ||
                 argument.rfind("--warn-on-errors=", 0) == 0 ||
                 argument.rfind("-warn-on-errors=", 0) == 0)
            options.warn_on_errors = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument,
                argument == "-warn-on-errors" || argument.rfind("-warn-on-errors=", 0) == 0
                    ? "-warn-on-errors" : "--warn-on-errors");
        else if (argument == "--validate-GVCF" || argument == "--validate-gvcf" ||
                 argument == "--gvcf" || argument == "-gvcf" ||
                 argument.rfind("--validate-GVCF=", 0) == 0 ||
                 argument.rfind("--validate-gvcf=", 0) == 0 ||
                 argument.rfind("--gvcf=", 0) == 0)
            options.validate_gvcf = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument,
                argument == "--validate-gvcf" || argument.rfind("--validate-gvcf=", 0) == 0
                    ? "--validate-gvcf" :
                argument == "--gvcf" || argument.rfind("--gvcf=", 0) == 0 || argument == "-gvcf"
                    ? (argument == "-gvcf" ? "-gvcf" : "--gvcf") : "--validate-GVCF");
        else if (argument == "--fail-gvcf-on-overlap" || argument == "--no-overlaps" ||
                 argument.rfind("--fail-gvcf-on-overlap=", 0) == 0 ||
                 argument.rfind("--no-overlaps=", 0) == 0)
            options.fail_on_overlap = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument,
                argument == "--no-overlaps" || argument.rfind("--no-overlaps=", 0) == 0
                    ? "--no-overlaps" : "--fail-gvcf-on-overlap");
        else if (has_option_value(argument, "--output-manifest") || has_option_value(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--quiet") {
            // Progress/logging compatibility flag.
        } else if (argument == "--disable-sequence-dictionary-validation" ||
                   argument == "-disable-sequence-dictionary-validation" ||
                   argument.rfind("--disable-sequence-dictionary-validation=", 0) == 0 ||
                   argument.rfind("-disable-sequence-dictionary-validation=", 0) == 0) {
            options.disable_sequence_dictionary_validation = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument,
                argument == "-disable-sequence-dictionary-validation" ||
                        argument.rfind("-disable-sequence-dictionary-validation=", 0) == 0
                    ? "-disable-sequence-dictionary-validation"
                    : "--disable-sequence-dictionary-validation");
        } else if (has_option_value(argument, "--java-options") || has_option_value(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    // Match GATK's argument constraints: filtered-record skipping is not
    // meaningful for GVCF validation and cannot be combined with overlap
    // enforcement.  Reject before opening the input so fallback can handle
    // the same fail-closed boundary deterministically.
    if (options.skip_filtered && options.validate_gvcf)
        throw std::invalid_argument(
            "DO_NOT_VALIDATE_FILTERED cannot be used in conjunction with VALIDATE_GVCF");
    if (options.skip_filtered && options.fail_on_overlap)
        throw std::invalid_argument(
            "DO_NOT_VALIDATE_FILTERED cannot be used in conjunction with FAIL_ON_OVERLAP");
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct Issue {
    std::string type;
    std::string contig;
    int position = 0;
    std::string message;
};

struct Summary {
    std::uint64_t input_records = 0;
    std::uint64_t validated_records = 0;
    std::uint64_t skipped_filtered = 0;
    std::uint64_t skipped_interval = 0;
    std::uint64_t skipped_excluded = 0;
    std::uint64_t warnings = 0;
    std::uint64_t errors = 0;
    std::uint64_t gvcf_blocks = 0;
    std::uint64_t gvcf_nonref_records = 0;
    std::uint64_t gvcf_overlaps = 0;
    std::uint64_t gvcf_uncovered_gaps = 0;
    std::uint64_t field_cardinality_errors = 0;
    std::uint64_t unused_alternate_errors = 0;
};

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return value;
}

bool record_filtered(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && std::string(name) != "PASS" && std::string(name) != ".") return true;
    }
    return false;
}

std::string contig_name(const bcf_hdr_t* header, const bcf1_t* record) {
    const auto* name = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    return name == nullptr ? "*" : name;
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    // HTSlib positions are zero-based.  Keep the interval contract in
    // zero-based half-open coordinates and honor a gVCF END annotation when
    // it extends beyond the concrete REF span.
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, values[0]);
    free(values);
    return end;
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_validate_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals) {
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

void append_validate_interval_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<fastgatk::io::IndexedInterval>& intervals,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<fastgatk::io::IndexedInterval> incoming;
    fastgatk::io::append_interval_selector(selector, header, incoming, stats);
    normalize_validate_intervals(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        intervals.insert(intervals.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<fastgatk::io::IndexedInterval> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < intervals.size() && right < incoming.size()) {
        if (intervals[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < intervals[left].rid) { ++right; continue; }
        const int begin = std::max(intervals[left].begin, incoming[right].begin);
        const int end = std::min(intervals[left].end, incoming[right].end);
        if (begin < end) {
            auto overlap = intervals[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (intervals[left].end < incoming[right].end) ++left;
        else ++right;
    }
    intervals.swap(intersection);
}

bool in_intervals(const bcf_hdr_t* header, bcf1_t* record,
                  const std::vector<fastgatk::io::IndexedInterval>& intervals) {
    if (intervals.empty()) return true;
    const auto end = record_end_exclusive(header, record);
    return std::any_of(intervals.begin(), intervals.end(), [&](const auto& interval) {
        return interval.rid == record->rid && end > interval.begin && record->pos < interval.end;
    });
}

std::set<std::string> excluded_set(const Options& options) {
    std::set<std::string> excluded;
    for (auto type : options.excluded_types) {
        type = upper(type);
        if (type != "ALL" && type != "REF" && type != "IDS" && type != "ALLELES" &&
            type != "CHR_COUNTS" && type != "FORMAT")
            throw std::invalid_argument("BAD_INPUT: unknown validation type: " + type);
        excluded.insert(type);
    }
    return excluded;
}

void add_issue(const Options& options, Summary& summary, std::vector<Issue>& issues,
               const std::string& type, const std::string& contig, int position,
               const std::string& message);

std::size_t genotype_width(const int allele_count, const int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t width = 1;
    for (int copy = 1; copy <= ploidy; ++copy) {
        const auto numerator = static_cast<std::size_t>(allele_count + copy - 1);
        if (numerator != 0 && width > std::numeric_limits<std::size_t>::max() / numerator)
            return 0;
        width = (width * numerator) / static_cast<std::size_t>(copy);
    }
    return width;
}

struct GenotypeShape {
    int ploidy = 0;
    bool present = false;
    bool uniform = true;
};

GenotypeShape genotype_shape(const bcf_hdr_t* header, bcf1_t* record) {
    GenotypeShape shape;
    if (bcf_hdr_nsamples(header) == 0) return shape;
    int32_t* values = nullptr;
    int count = 0;
    const int length = bcf_get_genotypes(header, record, &values, &count);
    if (length <= 0 || count <= 0 || values == nullptr) {
        free(values);
        return shape;
    }
    // bcf_get_genotypes returns a rectangular sample-major buffer.  A vector
    // end terminates shorter ploidies; missing alleles do not change the
    // declared ploidy and are therefore counted as occupied GT slots.
    const int samples = bcf_hdr_nsamples(header);
    const int stride = count / samples;
    if (stride <= 0 || count % samples != 0) {
        free(values);
        return shape;
    }
    for (int sample = 0; sample < samples; ++sample) {
        int sample_ploidy = 0;
        for (int allele = 0; allele < stride; ++allele) {
            const auto value = values[sample * stride + allele];
            if (value == bcf_int32_vector_end) break;
            ++sample_ploidy;
        }
        if (sample_ploidy == 0) continue;
        if (!shape.present) {
            shape.present = true;
            shape.ploidy = sample_ploidy;
        } else if (shape.ploidy != sample_ploidy) {
            shape.uniform = false;
        }
    }
    free(values);
    return shape;
}

std::size_t expected_cardinality(const bcf_hdr_t* header, const int id, const int column_type,
                                 const int allele_count, const int ploidy) {
    if (!bcf_hdr_idinfo_exists(header, column_type, id)) return 0;
    const auto length = bcf_hdr_id2length(header, column_type, id);
    if (length == BCF_VL_VAR) return std::numeric_limits<std::size_t>::max();
    if (length == BCF_VL_A) return static_cast<std::size_t>(std::max(0, allele_count - 1));
    if (length == BCF_VL_R) return static_cast<std::size_t>(std::max(0, allele_count));
    if (length == BCF_VL_G) {
        if (ploidy <= 0) return std::numeric_limits<std::size_t>::max();
        return genotype_width(allele_count, ploidy);
    }
    return static_cast<std::size_t>(bcf_hdr_id2number(header, column_type, id));
}

void validate_field_cardinality(const Options& options, Summary& summary,
                                std::vector<Issue>& issues, const bcf_hdr_t* header,
                                bcf1_t* record, const std::string& contig, const int position,
                                const int allele_count) {
    const auto shape = genotype_shape(header, record);
    if (shape.present && !shape.uniform)
        add_issue(options, summary, issues, "FORMAT", contig, position,
                  "GT ploidy differs between samples");
    const int ploidy = shape.present && shape.uniform ? shape.ploidy : 0;

    bcf_unpack(record, BCF_UN_ALL);
    for (int index = 0; index < record->n_info; ++index) {
        const auto& info = record->d.info[index];
        const auto expected = expected_cardinality(header, info.key, BCF_HL_INFO,
                                                   allele_count, ploidy);
        if (expected == std::numeric_limits<std::size_t>::max()) continue;
        if (static_cast<std::size_t>(std::max(0, info.len)) != expected) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, info.key);
            add_issue(options, summary, issues, "FORMAT", contig, position,
                      std::string("INFO/") + (name == nullptr ? "?" : name) +
                      " has " + std::to_string(info.len) + " values; expected " +
                      std::to_string(expected));
            ++summary.field_cardinality_errors;
        }
    }

    for (int index = 0; index < record->n_fmt; ++index) {
        const auto& format = record->d.fmt[index];
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, format.id);
        // GT is a special Number=1 declaration whose actual vector length is
        // the sample ploidy, so it is checked above rather than as a scalar.
        if (name != nullptr && std::string(name) == "GT") continue;
        const auto expected = expected_cardinality(header, format.id, BCF_HL_FMT,
                                                   allele_count, ploidy);
        if (expected == std::numeric_limits<std::size_t>::max() ||
            (bcf_hdr_id2length(header, BCF_HL_FMT, format.id) == BCF_VL_G && !shape.present))
            continue;
        if (static_cast<std::size_t>(std::max(0, format.n)) != expected) {
            add_issue(options, summary, issues, "FORMAT", contig, position,
                      std::string("FORMAT/") + (name == nullptr ? "?" : name) +
                      " has " + std::to_string(format.n) + " values per sample; expected " +
                      std::to_string(expected));
            ++summary.field_cardinality_errors;
        }
    }
}

// VariantContext.validateAlternateAlleles() has a stricter meaning than the
// wire-level VCF checks above: when genotypes are present, every concrete ALT
// allele must occur in at least one called genotype.  Keep this separate from
// genotype chromosome-count validation so --validation-type-to-exclude
// CHR_COUNTS does not accidentally disable the ALLELES contract.
void validate_alternate_usage(const Options& options, Summary& summary,
                              std::vector<Issue>& issues, const bcf_hdr_t* header,
                              bcf1_t* record, const std::string& contig,
                              const int position) {
    const int samples = bcf_hdr_nsamples(header);
    if (samples <= 0 || record->n_allele <= 1) return;
    int32_t* genotypes = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
    if (length <= 0 || genotypes == nullptr || count <= 0 || count % samples != 0) {
        free(genotypes);
        return;
    }
    const int stride = count / samples;
    std::vector<bool> used(static_cast<std::size_t>(record->n_allele - 1), false);
    for (int sample = 0; sample < samples; ++sample) {
        for (int copy = 0; copy < stride; ++copy) {
            const auto encoded = genotypes[sample * stride + copy];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            const int allele = bcf_gt_allele(encoded);
            if (allele > 0 && allele < record->n_allele)
                used[static_cast<std::size_t>(allele - 1)] = true;
        }
    }
    free(genotypes);
    bcf_unpack(record, BCF_UN_STR);
    for (int alt = 1; alt < record->n_allele; ++alt) {
        // VariantContext.validateAlternateAlleles() only requires concrete
        // alternates to participate in a called genotype.  Symbolic alleles
        // such as <DEL> and <CNV> are deliberately excluded from this
        // usage check (the same rule also keeps <NON_REF> out of ordinary
        // site validation); treating them as ordinary strings would reject
        // valid symbolic records that are intentionally 0/0 or no-call.
        const auto* allele = record->d.allele != nullptr ? record->d.allele[alt] : nullptr;
        if (allele != nullptr && allele[0] == '<' &&
            std::strlen(allele) >= 2 && allele[std::strlen(allele) - 1] == '>')
            continue;
        if (used[static_cast<std::size_t>(alt - 1)]) continue;
        add_issue(options, summary, issues, "ALLELES", contig, position,
                  "ALT allele " + std::to_string(alt) + " (" +
                  (allele == nullptr ? std::string("?") : std::string(allele)) +
                  ") is not observed in any called genotype");
        ++summary.unused_alternate_errors;
    }
}

void add_issue(const Options& options, Summary& summary, std::vector<Issue>& issues,
               const std::string& type, const std::string& contig, int position,
               const std::string& message) {
    Issue issue{type, contig, position, message};
    if (options.warn_on_errors) {
        ++summary.warnings;
        issues.push_back(std::move(issue));
    } else {
        ++summary.errors;
        issues.push_back(std::move(issue));
        throw std::runtime_error("VALIDATION_FAILURE: " + type + " at " + contig + ":" +
                                 std::to_string(position) + ": " + message);
    }
}

std::map<std::string, std::set<std::string>> load_dbsnp_ids(const std::string& path) {
    std::map<std::string, std::set<std::string>> ids;
    if (path.empty()) return ids;
    htsFile* file = bcf_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open dbSNP VCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(file);
    if (!header) {
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read dbSNP VCF header: " + path);
    }
    bcf1_t* record = bcf_init();
    while (fastgatk::io::read_variant_record(file, header, record, path) >= 0) {
        bcf_unpack(record, BCF_UN_STR);
        const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
        if (contig == nullptr || record->d.id == nullptr || std::string(record->d.id) == ".") continue;
        std::ostringstream key;
        key << contig << ':' << record->pos + 1;
        for (int alt = 0; alt < record->n_allele; ++alt) key << ':' << record->d.allele[alt];
        std::string id(record->d.id);
        std::size_t start = 0;
        while (start <= id.size()) {
            const auto end = id.find(';', start);
            ids[key.str()].insert(id.substr(start, end == std::string::npos ? std::string::npos : end - start));
            if (end == std::string::npos) break;
            start = end + 1;
        }
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(file);
    return ids;
}

std::int64_t info_integer(const bcf_hdr_t* header, bcf1_t* record, const char* field, bool& present) {
    int32_t* values = nullptr;
    int count = 0;
    const int length = bcf_get_info_int32(header, record, field, &values, &count);
    if (length <= 0 || count <= 0 || values == nullptr) {
        free(values);
        present = false;
        return 0;
    }
    present = true;
    const auto result = values[0] == bcf_int32_missing ? 0 : static_cast<std::int64_t>(values[0]);
    free(values);
    return result;
}

std::vector<std::int64_t> info_integers(const bcf_hdr_t* header, bcf1_t* record, const char* field) {
    int32_t* values = nullptr;
    int count = 0;
    const int length = bcf_get_info_int32(header, record, field, &values, &count);
    if (length <= 0 || count <= 0 || values == nullptr) {
        free(values);
        return {};
    }
    std::vector<std::int64_t> result;
    result.reserve(static_cast<std::size_t>(length));
    for (int index = 0; index < length; ++index)
        result.push_back(values[index] == bcf_int32_missing ? 0 : static_cast<std::int64_t>(values[index]));
    free(values);
    return result;
}

struct AlleleCounts {
    std::uint64_t chromosomes = 0;
    std::vector<std::uint64_t> alternate;
};

AlleleCounts genotype_counts(const bcf_hdr_t* header, bcf1_t* record, int alt_count) {
    AlleleCounts counts;
    counts.alternate.assign(static_cast<std::size_t>(std::max(0, alt_count)), 0);
    if (bcf_hdr_nsamples(header) == 0) return counts;
    int32_t* genotypes = nullptr;
    int count = 0;
    const int length = bcf_get_genotypes(header, record, &genotypes, &count);
    if (length <= 0 || count <= 0 || count % bcf_hdr_nsamples(header) != 0) {
        free(genotypes);
        return counts;
    }
    const int ploidy = count / bcf_hdr_nsamples(header);
    for (int sample = 0; sample < bcf_hdr_nsamples(header); ++sample) {
        bool missing = false;
        std::vector<int> sample_alternates;
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto value = genotypes[sample * ploidy + allele];
            if (bcf_gt_is_missing(value) || value == bcf_int32_vector_end) {
                missing = true;
                continue;
            }
            const int index = bcf_gt_allele(value);
            if (index < 0 || index > alt_count) {
                missing = true;
                continue;
            }
            if (index > 0) sample_alternates.push_back(index - 1);
        }
        if (!missing) {
            counts.chromosomes += static_cast<std::uint64_t>(ploidy);
            for (const auto index : sample_alternates)
                ++counts.alternate[static_cast<std::size_t>(index)];
        }
    }
    free(genotypes);
    return counts;
}

std::string record_key(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    std::ostringstream key;
    key << contig_name(header, record) << ':' << record->pos + 1;
    for (int index = 0; index < record->n_allele; ++index) key << ':' << record->d.allele[index];
    return key.str();
}

bool validate_reference(const Options& options, Summary& summary, std::vector<Issue>& issues,
                        faidx_t* reference, const bcf_hdr_t* header, bcf1_t* record) {
    if (reference == nullptr) return true;
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele == 0 || record->d.allele[0] == nullptr) return true;
    const auto contig = contig_name(header, record);
    const int length = static_cast<int>(std::strlen(record->d.allele[0]));
    int fetched_length = 0;
    char* bases = faidx_fetch_seq(reference, contig.c_str(), record->pos,
                                   record->pos + length - 1, &fetched_length);
    if (bases == nullptr || fetched_length != length) {
        free(bases);
        add_issue(options, summary, issues, "REF", contig, record->pos + 1,
                  "reference FASTA does not cover the reported REF allele");
        return false;
    }
    const std::string observed = upper(std::string(bases, static_cast<std::size_t>(fetched_length)));
    const std::string reported = upper(record->d.allele[0]);
    free(bases);
    if (observed != reported) {
        add_issue(options, summary, issues, "REF", contig, record->pos + 1,
                  "reported REF=" + reported + " but reference has " + observed);
        return false;
    }
    return true;
}

void validate_sequence_dictionary(const Options& options, Summary& summary,
                                  std::vector<Issue>& issues, faidx_t* reference,
                                  const bcf_hdr_t* header) {
    if (options.disable_sequence_dictionary_validation || reference == nullptr || header == nullptr)
        return;
    for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
        const auto* name = bcf_hdr_id2name(header, rid);
        if (name == nullptr) continue;
        const auto reference_length = faidx_seq_len64(reference, name);
        const auto* id_info = header->id[BCF_DT_CTG][rid].val;
        const auto declared_length = id_info == nullptr ? std::int64_t{0}
                                                         : static_cast<std::int64_t>(id_info->info[0]);
        if (reference_length < 0) {
            add_issue(options, summary, issues, "DICTIONARY", name, 0,
                      "VCF contig is absent from reference sequence dictionary");
        } else if (declared_length > 0 && declared_length != reference_length) {
            add_issue(options, summary, issues, "DICTIONARY", name, 0,
                      "VCF contig length does not match reference sequence dictionary");
        }
    }
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto excluded = excluded_set(options);
    const bool validate_all = excluded.empty();
    const bool exclude_all = excluded.count("ALL") != 0;
    const bool validate_ref = !exclude_all && (validate_all || excluded.count("REF") == 0);
    const bool validate_ids = !exclude_all && (validate_all || excluded.count("IDS") == 0);
    // GATK deliberately disables strict alternate-allele validation for GVCF
    // traversal because <NON_REF> and reference-confidence records are not
    // ordinary variant records (ValidateVariants.calculateValidationTypesToApply).
    const bool validate_alleles = !options.validate_gvcf &&
                                  !exclude_all && (validate_all || excluded.count("ALLELES") == 0);
    const bool validate_chr_counts = !exclude_all && (validate_all || excluded.count("CHR_COUNTS") == 0);
    const bool validate_fields = !exclude_all && (validate_all || excluded.count("FORMAT") == 0);
    if (validate_ref && options.reference.empty() && !validate_all)
        throw std::invalid_argument("BAD_INPUT: REF validation requires -R/--reference");

    htsFile* input = bcf_open(options.input.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/BCF: " + options.input);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (!header) {
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + options.input);
    }
    std::vector<fastgatk::io::IndexedInterval> intervals;
    fastgatk::io::IntervalFileStats interval_stats;
    bool first_interval_selector = true;
    for (const auto& selector : options.regions) {
        append_validate_interval_selector_with_rule(
            selector, header, intervals, interval_stats, options.interval_set_rule,
            first_interval_selector);
        first_interval_selector = false;
    }
    normalize_validate_intervals(intervals);
    std::vector<fastgatk::io::IndexedInterval> exclusions;
    fastgatk::io::IntervalFileStats exclusion_stats;
    for (const auto& selector : options.exclusions)
        fastgatk::io::append_interval_selector(selector, header, exclusions, exclusion_stats);
    const auto apply_padding = [](std::vector<fastgatk::io::IndexedInterval>& values, int padding) {
        if (padding <= 0) return;
        for (auto& interval : values) {
            interval.begin = std::max(0, interval.begin - padding);
            if (interval.end > std::numeric_limits<int>::max() - padding)
                interval.end = std::numeric_limits<int>::max();
            else interval.end += padding;
        }
    };
    apply_padding(intervals, options.interval_padding);
    apply_padding(exclusions, options.exclusion_padding);

    faidx_t* reference = nullptr;
    if (validate_ref && !options.reference.empty()) {
        reference = fai_load(options.reference.c_str());
        if (!reference) {
            bcf_hdr_destroy(header);
            bcf_close(input);
            throw std::runtime_error("BAD_INPUT: cannot load reference FASTA/index: " + options.reference);
        }
    }
    const auto dbsnp_ids = validate_ids ? load_dbsnp_ids(options.dbsnp) : std::map<std::string, std::set<std::string>>{};
    Summary summary;
    std::vector<Issue> issues;
    validate_sequence_dictionary(options, summary, issues, reference, header);
    bcf1_t* record = bcf_init();
    std::map<int, int> previous_end;
    std::map<int, int> previous_start;
    std::map<int, int> first_start;
    std::map<int, bool> previous_reference_block;
    std::map<int, std::vector<std::pair<int, int>>> gvcf_spans;
    while (fastgatk::io::read_variant_record(input, header, record,
                                           options.input) >= 0) {
        ++summary.input_records;
        bcf_unpack(record, BCF_UN_ALL);
        if (!in_intervals(header, record, intervals)) {
            ++summary.skipped_interval;
            bcf_clear(record);
            continue;
        }
        if (!exclusions.empty() && in_intervals(header, record, exclusions)) {
            ++summary.skipped_excluded;
            bcf_clear(record);
            continue;
        }
        if (options.skip_filtered && record_filtered(header, record)) {
            ++summary.skipped_filtered;
            bcf_clear(record);
            continue;
        }
        ++summary.validated_records;
        const auto contig = contig_name(header, record);
        const int position = record->pos + 1;
        const int alt_count = std::max(0, record->n_allele - 1);
        if (validate_alleles) {
            std::set<std::string> seen;
            bcf_unpack(record, BCF_UN_STR);
            if (record->n_allele < 2 || record->d.allele[0] == nullptr || std::strlen(record->d.allele[0]) == 0)
                add_issue(options, summary, issues, "ALLELES", contig, position, "record has no valid REF/ALT allele");
            for (int alt = 1; alt < record->n_allele; ++alt) {
                const std::string value = record->d.allele[alt] == nullptr ? "" : record->d.allele[alt];
                if (value.empty() || value == record->d.allele[0] || !seen.insert(value).second)
                    add_issue(options, summary, issues, "ALLELES", contig, position, "duplicate, empty, or REF-equal ALT allele");
            }
            validate_alternate_usage(options, summary, issues, header, record, contig, position);
        }
        if (validate_ref) validate_reference(options, summary, issues, reference, header, record);
        if (validate_ids && !options.dbsnp.empty() && record->d.id != nullptr && std::string(record->d.id) != ".") {
            const auto key = record_key(header, record);
            const auto found = dbsnp_ids.find(key);
            if (found == dbsnp_ids.end())
                add_issue(options, summary, issues, "IDS", contig, position, "record ID is absent from dbSNP track");
        }
        if (validate_chr_counts) {
            bool an_present = false;
            bool ac_present = false;
            const auto an = info_integer(header, record, "AN", an_present);
            const auto ac = info_integers(header, record, "AC");
            ac_present = !ac.empty();
            const auto counts = genotype_counts(header, record, alt_count);
            if (an_present && an != static_cast<std::int64_t>(counts.chromosomes))
                add_issue(options, summary, issues, "CHR_COUNTS", contig, position,
                          "AN does not match called genotype chromosome count");
            if (ac_present) {
                if (ac.size() != counts.alternate.size())
                    add_issue(options, summary, issues, "CHR_COUNTS", contig, position,
                              "AC allele count length does not match ALT count");
                const auto common = std::min(ac.size(), counts.alternate.size());
                for (std::size_t index = 0; index < common; ++index) {
                    if (ac[index] < 0 || ac[index] != static_cast<std::int64_t>(counts.alternate[index]))
                        add_issue(options, summary, issues, "CHR_COUNTS", contig, position,
                                  "AC does not match called genotype allele count");
                }
            }
        }
        if (validate_fields)
            validate_field_cardinality(options, summary, issues, header, record,
                                       contig, position, record->n_allele);
        if (options.validate_gvcf) {
            bcf_unpack(record, BCF_UN_STR | BCF_UN_INFO);
            bool has_nonref = false;
            for (int alt = 1; alt < record->n_allele; ++alt)
                if (record->d.allele[alt] != nullptr && std::string(record->d.allele[alt]) == "<NON_REF>") has_nonref = true;
            if (!has_nonref)
                add_issue(options, summary, issues, "GVCF", contig, position, "record lacks <NON_REF> allele");
            const bool reference_block = record->n_allele == 2 && has_nonref;
            bool end_present = false;
            const auto end_value = info_integer(header, record, "END", end_present);
            const int end = end_present && end_value >= position ? static_cast<int>(end_value) : position;
            if (reference != nullptr) {
                const auto reference_length = faidx_seq_len64(reference, contig.c_str());
                if (reference_length >= 0 && static_cast<std::int64_t>(end) > reference_length)
                    add_issue(options, summary, issues, "GVCF_END", contig, position,
                              "GVCF END exceeds the reference sequence length");
            }
            if (reference_block) ++summary.gvcf_blocks;
            else ++summary.gvcf_nonref_records;
            if (first_start.find(record->rid) == first_start.end())
                first_start[record->rid] = position;
            const auto prior_contig = previous_start.find(record->rid);
            if (prior_contig != previous_start.end() && position < prior_contig->second)
                add_issue(options, summary, issues, "GVCF_ORDER", contig, position, "records are not coordinate sorted");
            const auto prior_end = previous_end.find(record->rid);
            if (prior_end != previous_end.end() && position <= prior_end->second) {
                ++summary.gvcf_overlaps;
                if (options.fail_on_overlap && (reference_block || previous_reference_block[record->rid]))
                    add_issue(options, summary, issues, "GVCF_OVERLAP", contig, position, "overlapping GVCF blocks");
            }
            // Every GVCF record contributes coverage, not only <NON_REF>
            // reference blocks.  Keep all record spans for the final union
            // check; restricting gap detection to adjacent reference blocks
            // would miss uncovered sequence between two concrete variants.
            gvcf_spans[record->rid].emplace_back(record->pos, end);
            previous_start[record->rid] = position;
            previous_end[record->rid] = std::max(previous_end[record->rid], end);
            previous_reference_block[record->rid] = reference_block;
        }
        bcf_clear(record);
    }
    // Validate complete GVCF coverage against the requested interval union or
    // the entire sequence dictionary.  Java GATK treats concrete variants and
    // reference blocks alike for coverage; merging half-open spans here also
    // handles adjacent records and overlapping haplotypes without double
    // counting.  The coverage metric counts uncovered segments (rather than
    // bases) to keep telemetry bounded on malformed large files.
    if (options.validate_gvcf) {
        std::map<int, std::vector<std::pair<int, int>>> targets;
        if (!intervals.empty()) {
            for (const auto& interval : intervals)
                if (interval.end > interval.begin)
                    targets[interval.rid].emplace_back(interval.begin, interval.end);
        } else {
            for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
                const auto* name = bcf_hdr_id2name(header, rid);
                if (name == nullptr) continue;
                std::int64_t length = -1;
                if (reference != nullptr)
                    length = faidx_seq_len64(reference, name);
                if (length < 0 && header->id[BCF_DT_CTG][rid].val != nullptr)
                    length = static_cast<std::int64_t>(header->id[BCF_DT_CTG][rid].val->info[0]);
                if (length > 0 && length <= std::numeric_limits<int>::max())
                    targets[rid].emplace_back(0, static_cast<int>(length));
            }
        }

        // GATK applies -XL after -L (and applies interval padding before the
        // set operation).  Subtract exclusions from the coverage targets so
        // a deliberately omitted span does not create a false GVCF gap.
        if (!exclusions.empty()) {
            std::map<int, std::vector<std::pair<int, int>>> excluded_spans;
            for (const auto& exclusion : exclusions)
                if (exclusion.rid >= 0 && exclusion.end > exclusion.begin)
                    excluded_spans[exclusion.rid].emplace_back(exclusion.begin, exclusion.end);
            for (auto& [rid, target_spans] : targets) {
                auto exclusion_it = excluded_spans.find(rid);
                if (exclusion_it == excluded_spans.end()) continue;
                auto& spans = exclusion_it->second;
                std::sort(spans.begin(), spans.end());
                std::vector<std::pair<int, int>> merged_exclusions;
                for (const auto span : spans) {
                    if (merged_exclusions.empty() || span.first > merged_exclusions.back().second)
                        merged_exclusions.push_back(span);
                    else
                        merged_exclusions.back().second = std::max(merged_exclusions.back().second, span.second);
                }
                std::vector<std::pair<int, int>> remaining;
                for (const auto target : target_spans) {
                    int cursor = target.first;
                    for (const auto exclusion : merged_exclusions) {
                        if (exclusion.second <= cursor) continue;
                        if (exclusion.first >= target.second) break;
                        if (exclusion.first > cursor)
                            remaining.emplace_back(cursor, std::min(exclusion.first, target.second));
                        cursor = std::max(cursor, exclusion.second);
                        if (cursor >= target.second) break;
                    }
                    if (cursor < target.second) remaining.emplace_back(cursor, target.second);
                }
                target_spans = std::move(remaining);
            }
        }

        std::string first_gap_contig = "*";
        int first_gap_position = 0;
        for (auto& [rid, target_spans] : targets) {
            auto& observed = gvcf_spans[rid];
            std::sort(target_spans.begin(), target_spans.end());
            std::vector<std::pair<int, int>> merged_targets;
            for (const auto span : target_spans) {
                if (merged_targets.empty() || span.first > merged_targets.back().second)
                    merged_targets.push_back(span);
                else
                    merged_targets.back().second = std::max(merged_targets.back().second, span.second);
            }
            std::sort(observed.begin(), observed.end());
            std::vector<std::pair<int, int>> merged_observed;
            for (const auto span : observed) {
                if (span.second <= span.first) continue;
                if (merged_observed.empty() || span.first > merged_observed.back().second)
                    merged_observed.push_back(span);
                else
                    merged_observed.back().second = std::max(merged_observed.back().second, span.second);
            }
            const auto* contig_name = bcf_hdr_id2name(header, rid);
            for (const auto target : merged_targets) {
                int cursor = target.first;
                for (const auto span : merged_observed) {
                    if (span.second <= cursor) continue;
                    if (span.first > cursor) {
                        ++summary.gvcf_uncovered_gaps;
                        if (first_gap_position == 0) {
                            first_gap_contig = contig_name == nullptr ? "*" : contig_name;
                            first_gap_position = cursor + 1;
                        }
                    }
                    cursor = std::max(cursor, span.second);
                    if (cursor >= target.second) break;
                }
                if (cursor < target.second) {
                    ++summary.gvcf_uncovered_gaps;
                    if (first_gap_position == 0) {
                        first_gap_contig = contig_name == nullptr ? "*" : contig_name;
                        first_gap_position = cursor + 1;
                    }
                }
            }
        }
        if (summary.gvcf_uncovered_gaps > 0)
            add_issue(options, summary, issues, "GVCF_COVERAGE", first_gap_contig,
                      first_gap_position,
                      "GVCF coverage does not reach the end of a requested sequence or interval; "
                      "selected GVCF records contain uncovered coordinate gaps");
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(input);
    if (reference) fai_destroy(reference);

    if (!options.manifest.empty()) {
        std::ofstream manifest(options.manifest);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create manifest: " + options.manifest);
        manifest << "{\"schema_version\":1,\"tool\":\"ValidateVariants\","
                 << "\"implementation\":\"fastgatk-validate-variants\",\"status\":\""
                 << (summary.errors == 0 ? "pass" : "failed") << "\","
                 << "\"compatibility\":{\"vcf_parse\":true,\"alleles\":" << (validate_alleles ? "true" : "false")
                 << ",\"ref\":" << (validate_ref ? "true" : "false")
                 << ",\"ids\":" << (validate_ids && !options.dbsnp.empty() ? "true" : "false")
                 << ",\"chr_counts\":" << (validate_chr_counts ? "true" : "false")
                 << ",\"field_cardinality\":" << (validate_fields ? "true" : "false")
                 << ",\"gvcf\":" << (options.validate_gvcf ? "true" : "false")
                 << ",\"interval_exclusion\":" << (!exclusions.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_padding\":" << options.interval_padding
                 << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                 << ",\"sequence_dictionary_validation\":"
                 << (!options.disable_sequence_dictionary_validation ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":false},\"telemetry\":{"
                 << "\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << summary.input_records
                 << ",\"validated_records\":" << summary.validated_records
                 << ",\"skipped_filtered\":" << summary.skipped_filtered
                 << ",\"skipped_interval\":" << summary.skipped_interval
                 << ",\"skipped_excluded\":" << summary.skipped_excluded
                 << ",\"warnings\":" << summary.warnings
                 << ",\"errors\":" << summary.errors
                 << ",\"gvcf_blocks\":" << summary.gvcf_blocks
                 << ",\"gvcf_overlaps\":" << summary.gvcf_overlaps
                 << ",\"gvcf_uncovered_gaps\":" << summary.gvcf_uncovered_gaps
                 << ",\"field_cardinality_errors\":" << summary.field_cardinality_errors
                 << ",\"unused_alternate_errors\":" << summary.unused_alternate_errors
                 << ",\"interval_list_inputs\":" << interval_stats.files
                 << ",\"interval_list_records\":" << interval_stats.records
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << "},\"issues\":[";
        for (std::size_t index = 0; index < issues.size(); ++index) {
            if (index != 0) manifest << ',';
            manifest << "{\"type\":\"" << issues[index].type << "\",\"contig\":\""
                     << issues[index].contig << "\",\"position\":" << issues[index].position
                     << ",\"message\":\"";
            for (const char character : issues[index].message) {
                if (character == '"' || character == '\\') manifest << '\\';
                manifest << character;
            }
            manifest << "\"}";
        }
        manifest << "]}\n";
    }
    if (summary.errors != 0 && !options.warn_on_errors)
        throw std::runtime_error("VALIDATION_FAILURE: validation failed");
    if (!options.output.empty()) {
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create validation report");
        output << "tool\tstatus\tinput_records\tvalidated_records\twarnings\terrors\n"
               << "ValidateVariants\t" << (summary.errors == 0 ? "PASS" : "WARN") << '\t'
               << summary.input_records << '\t' << summary.validated_records << '\t'
               << summary.warnings << '\t' << summary.errors << '\n';
    }
    std::cout << "{\"tool\":\"ValidateVariants\",\"status\":\""
              << (summary.errors == 0 ? "pass" : "warn") << "\",\"input_records\":"
              << summary.input_records << ",\"validated_records\":" << summary.validated_records
              << ",\"skipped_interval\":" << summary.skipped_interval
              << ",\"skipped_excluded\":" << summary.skipped_excluded
              << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
              << ",\"warnings\":" << summary.warnings << ",\"errors\":" << summary.errors << "}\n";
    return 0;
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for ValidateVariants");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse(argc, argv);
        return run_tool(options, fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-validate-variants: " << error.what() << '\n';
        return 2;
    }
}
