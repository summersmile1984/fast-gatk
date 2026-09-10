#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/java_numeric.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/output.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string manifest;
    std::string reference;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    // GATK's interval merging rule controls whether abutting half-open
    // selectors are coalesced.  ALL is the default; OVERLAPPING_ONLY keeps
    // adjacent selectors as distinct intervals for output matching.
    bool merge_abutting_intervals = true;
    enum class VariantOutputFiltering { Overlaps, StartsIn, EndsIn, Contained, Anywhere };
    VariantOutputFiltering variant_output_filtering = VariantOutputFiltering::Overlaps;
    bool variant_output_filtering_explicit = false;
    std::vector<std::string> exclusions;
    int interval_padding = 0;
    int exclusion_padding = 0;
    std::vector<std::string> fields;
    std::vector<std::string> genotype_fields;
    std::vector<std::string> allele_fields;
    std::vector<std::string> allele_genotype_fields;
    bool show_filtered = false;
    bool split_multi_allelic = false;
    bool moltenize = false;
    bool error_if_missing = false;
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

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-variants-to-table (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE              input VCF/BCF\n"
                         "  -O, --output FILE               output tab-delimited table\n"
                         "  -F, --fields FIELD              site/INFO field (repeatable)\n"
                         "  -GF, --genotype-fields FIELD    FORMAT field (repeatable)\n"
                         "  -ASF FIELD                      allele-specific INFO field\n"
                         "  -ASGF FIELD                     allele-specific FORMAT field\n"
                         "  -SMA, --split-multi-allelic     one row per ALT allele\n"
                         "  --show-filtered, --raw          include filtered records\n"
                         "  --moltenize                      RecordID/Sample/Variable/Value output\n"
                         "  --error-if-missing-data, -EMD    fail when a requested value is absent\n"
                         "  -L, --intervals REGION           interval selector (repeatable)\n"
                         "      --interval-set-rule RULE     UNION (default) or INTERSECTION\n"
                         "      --interval-merging-rule RULE ALL (default) or OVERLAPPING_ONLY\n"
                         "      --variant-output-filtering MODE STARTS_IN/ENDS_IN/OVERLAPS/CONTAINED/ANYWHERE\n"
                         "  -XL, --exclude-intervals REGION  exclusion selector (repeatable)\n"
                         "  -ip, --interval-padding N        include interval padding\n"
                         "  -ixp, --interval-exclusion-padding N  exclusion padding\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-V" || is_option(argument, "--variant"))
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
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
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (is_option(argument, "--interval-merging-rule") || argument == "-imr") {
            auto value = require_value(index, argc, argv, argument, "--interval-merging-rule", "-imr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "ALL") options.merge_abutting_intervals = true;
            else if (value == "OVERLAPPING_ONLY") options.merge_abutting_intervals = false;
            else throw std::invalid_argument("invalid --interval-merging-rule: " + value);
        }
        else if (is_option(argument, "--variant-output-filtering")) {
            auto value = require_value(index, argc, argv, argument, "--variant-output-filtering");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "OVERLAPS") options.variant_output_filtering = Options::VariantOutputFiltering::Overlaps;
            else if (value == "STARTS_IN") options.variant_output_filtering = Options::VariantOutputFiltering::StartsIn;
            else if (value == "ENDS_IN") options.variant_output_filtering = Options::VariantOutputFiltering::EndsIn;
            else if (value == "CONTAINED") options.variant_output_filtering = Options::VariantOutputFiltering::Contained;
            else if (value == "ANYWHERE") options.variant_output_filtering = Options::VariantOutputFiltering::Anywhere;
            else throw std::invalid_argument("invalid --variant-output-filtering: " + value);
            options.variant_output_filtering_explicit = true;
        }
        else if (argument == "-XL" || is_option(argument, "--exclude-intervals"))
            options.exclusions.push_back(require_value(
                index, argc, argv, argument, "--exclude-intervals", "-XL"));
        else if (argument == "-ip" || is_option(argument, "--interval-padding"))
            options.interval_padding = std::stoi(require_value(
                index, argc, argv, argument, "--interval-padding", "-ip"));
        else if (argument == "-ixp" || is_option(argument, "--interval-exclusion-padding"))
            options.exclusion_padding = std::stoi(require_value(
                index, argc, argv, argument, "--interval-exclusion-padding", "-ixp"));
        else if (argument == "-F" || is_option(argument, "--fields"))
            options.fields.push_back(require_value(index, argc, argv, argument, "--fields", "-F"));
        else if (argument == "-GF" || is_option(argument, "--genotype-fields"))
            options.genotype_fields.push_back(require_value(index, argc, argv, argument, "--genotype-fields", "-GF"));
        else if (argument == "-ASF" || is_option(argument, "--allele-specific-fields") ||
                 is_option(argument, "--asFieldsToTake"))
            options.allele_fields.push_back(require_value(
                index, argc, argv, argument,
                argument.rfind("--asFieldsToTake", 0) == 0
                    ? "--asFieldsToTake" : "--allele-specific-fields", "-ASF"));
        else if (argument == "-ASGF" || is_option(argument, "--allele-specific-genotype-fields") ||
                 is_option(argument, "--asGenotypeFieldsToTake"))
            options.allele_genotype_fields.push_back(require_value(
                index, argc, argv, argument,
                argument.rfind("--asGenotypeFieldsToTake", 0) == 0
                    ? "--asGenotypeFieldsToTake" : "--allele-specific-genotype-fields", "-ASGF"));
        else if (argument == "-SMA" || argument == "--split-multi-allelic") options.split_multi_allelic = true;
        else if (argument == "--show-filtered" || argument == "--raw") options.show_filtered = true;
        else if (argument == "--moltenize") options.moltenize = true;
        else if (argument == "-EMD" || argument == "--error-if-missing-data") options.error_if_missing = true;
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Accepted compatibility flags. Table extraction itself does not
            // require a sequence dictionary, but never drops records silently.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.interval_padding < 0 || options.exclusion_padding < 0)
        throw std::invalid_argument("interval padding must be non-negative");
    return options;
}

#if FASTGATK_HAS_HTSLIB

constexpr const char* kMissing = "NA";

// HTSlib decodes numeric INFO values to binary32 when reading VCF/BCF.  For
// text VCF, htsjdk/GATK keeps the original INFO token for VariantsToTable
// rendering (including long decimals, exponent spelling, and a literal `.`).
// Keep that lexical value at the Host/file boundary; BCF has no equivalent
// source token and intentionally continues through the typed decoder.
struct RawInfoRecord {
    std::unordered_map<std::string, std::vector<std::string>> values;
};

bool is_text_vcf(const std::string& path) {
    return fastgatk::io::has_suffix_ci(path, ".vcf") ||
           fastgatk::io::has_suffix_ci(path, ".vcf.gz");
}

std::vector<std::string> split_raw_csv(const std::string& text) {
    std::vector<std::string> values;
    std::size_t start = 0;
    while (true) {
        const auto end = text.find(',', start);
        values.push_back(text.substr(start, end == std::string::npos
            ? std::string::npos : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return values;
}

class RawInfoReader {
public:
    explicit RawInfoReader(const std::string& path) {
        if (is_text_vcf(path) && path != "-") {
            file_ = hts_open(path.c_str(), "r");
            if (file_ != nullptr && hts_get_format(file_)->format != vcf) {
                hts_close(file_);
                file_ = nullptr;
            }
        }
    }

    RawInfoReader(const RawInfoReader&) = delete;
    RawInfoReader& operator=(const RawInfoReader&) = delete;

    ~RawInfoReader() {
        free(line_.s);
        if (file_ != nullptr) hts_close(file_);
    }

    bool enabled() const noexcept { return file_ != nullptr; }

    bool next(RawInfoRecord& record) {
        if (file_ == nullptr) return false;
        while (hts_getline(file_, '\n', &line_) >= 0) {
            if (line_.l == 0 || line_.s[0] == '#') continue;
            record = RawInfoRecord{};
            std::size_t field_start = 0;
            std::string info;
            const std::string_view line(line_.s, line_.l);
            for (int field = 0; field <= 7; ++field) {
                const auto field_end = line.find('\t', field_start);
                if (field_end == std::string_view::npos && field < 7) {
                    // Preserve record alignment with bcf_read even for a
                    // malformed data line; HTSlib reports the parse error.
                    break;
                }
                if (field == 7) {
                    info.assign(line.data() + field_start,
                                (field_end == std::string_view::npos ? line.size() : field_end) - field_start);
                }
                if (field_end == std::string_view::npos) break;
                field_start = field_end + 1;
            }
            if (!info.empty() && info != ".") {
                std::size_t token_start = 0;
                while (true) {
                    const auto token_end = info.find(';', token_start);
                    const auto token = info.substr(token_start,
                        token_end == std::string::npos ? std::string::npos : token_end - token_start);
                    if (!token.empty() && token != ".") {
                        const auto equals = token.find('=');
                        if (equals == std::string::npos) {
                            record.values.emplace(token, std::vector<std::string>{"true"});
                        } else if (equals > 0) {
                            record.values[token.substr(0, equals)] =
                                split_raw_csv(token.substr(equals + 1));
                        }
                    }
                    if (token_end == std::string::npos) break;
                    token_start = token_end + 1;
                }
            }
            return true;
        }
        return false;
    }

private:
    htsFile* file_ = nullptr;
    kstring_t line_{0, 0, nullptr};
};

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    // HTSlib POS is zero-based and VCF END is one-based inclusive; the
    // resulting integer is a zero-based half-open end suitable for the
    // shared interval selector.  This keeps gVCF reference blocks visible
    // when a requested locus lies inside the block rather than at its POS.
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, static_cast<int>(values[0]));
    free(values);
    return end;
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_table_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals,
                               bool merge_abutting = true) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<fastgatk::io::IndexedInterval> merged;
    merged.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (interval.end <= interval.begin) continue;
        const bool joins = !merged.empty() && (merge_abutting
            ? interval.begin <= merged.back().end
            : interval.begin < merged.back().end);
        if (!merged.empty() && merged.back().rid == interval.rid && joins)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    intervals.swap(merged);
}

void append_table_interval_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<fastgatk::io::IndexedInterval>& intervals,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector,
    bool merge_abutting) {
    std::vector<fastgatk::io::IndexedInterval> incoming;
    fastgatk::io::append_interval_selector(selector, header, incoming, stats);
    normalize_table_intervals(incoming, merge_abutting);
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

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return value;
}

std::vector<std::string> split_csv(const std::string& text) {
    std::vector<std::string> values;
    std::string current;
    std::istringstream stream(text);
    while (std::getline(stream, current, ',')) values.push_back(current.empty() ? kMissing : current);
    if (values.empty()) values.push_back(kMissing);
    return values;
}

std::string join(const std::vector<std::string>& values, const std::string& delimiter = ",") {
    std::ostringstream output;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) output << delimiter;
        output << values[index];
    }
    return output.str();
}

const char* variant_output_filtering_name(Options::VariantOutputFiltering mode) {
    switch (mode) {
        case Options::VariantOutputFiltering::StartsIn: return "STARTS_IN";
        case Options::VariantOutputFiltering::EndsIn: return "ENDS_IN";
        case Options::VariantOutputFiltering::Contained: return "CONTAINED";
        case Options::VariantOutputFiltering::Anywhere: return "ANYWHERE";
        case Options::VariantOutputFiltering::Overlaps: return "OVERLAPS";
    }
    return "OVERLAPS";
}

std::string format_float(float value) {
    if (bcf_float_is_missing(value) || bcf_float_is_vector_end(value) || !std::isfinite(value))
        return kMissing;
    return fastgatk::io::java_float(value);
}

std::string format_int(int32_t value) {
    if (value == bcf_int32_missing || value == bcf_int32_vector_end) return kMissing;
    return std::to_string(value);
}

int number_code(const bcf_hdr_t* header, const char* kind, const std::string& field) {
    const auto id = bcf_hdr_id2int(header, BCF_DT_ID, field.c_str());
    const int type = kind == std::string("INFO") ? BCF_HL_INFO : BCF_HL_FMT;
    if (id < 0 || !bcf_hdr_idinfo_exists(header, type, id)) return BCF_VL_FIXED;
    return bcf_hdr_id2length(header, type, id);
}

std::size_t genotype_count(int allele_count, int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t result = 1;
    for (int position = 1; position <= ploidy; ++position) {
        const auto numerator = static_cast<std::size_t>(allele_count + position - 1);
        if (numerator != 0 && result > std::numeric_limits<std::size_t>::max() / numerator) return 0;
        result *= numerator;
        result /= static_cast<std::size_t>(position);
    }
    return result;
}

std::size_t rank_genotype(const std::vector<int>& alleles) {
    if (alleles.empty()) return 0;
    const int max_allele = *std::max_element(alleles.begin(), alleles.end());
    std::size_t rank = 0;
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        rank += genotype_count(max_allele + 1, static_cast<int>(index + 1));
    }
    // The compact diploid form is by far the common case.  Use the VCF
    // triangular rank directly; the generic fallback above is only used to
    // size malformed/non-diploid vectors and is never used for projection.
    if (alleles.size() == 2) {
        const int low = std::min(alleles[0], alleles[1]);
        const int high = std::max(alleles[0], alleles[1]);
        return static_cast<std::size_t>(high * (high + 1) / 2 + low);
    }
    return rank;
}

std::vector<std::string> project_values(std::vector<std::string> values, int length_code,
                                        int alt_index, int alt_count, int ploidy = 2) {
    if (alt_index < 0 || values.empty() || values.size() == 1) return values;
    if (length_code == BCF_VL_A) {
        return alt_index < static_cast<int>(values.size())
            ? std::vector<std::string>{values[static_cast<std::size_t>(alt_index)]}
            : std::vector<std::string>{kMissing};
    }
    if (length_code == BCF_VL_R) {
        std::vector<std::string> result;
        result.push_back(values.front());
        if (alt_index + 1 < static_cast<int>(values.size())) result.push_back(values[static_cast<std::size_t>(alt_index + 1)]);
        else result.push_back(kMissing);
        return result;
    }
    if (length_code == BCF_VL_G && ploidy == 2 && alt_count > 0) {
        std::vector<std::string> result;
        const std::vector<std::vector<int>> genotypes{
            {0, 0}, {0, alt_index + 1}, {alt_index + 1, alt_index + 1}};
        for (const auto& genotype : genotypes) {
            const auto index = rank_genotype(genotype);
            result.push_back(index < values.size() ? values[index] : kMissing);
        }
        return result;
    }
    return values;
}

std::vector<std::string> info_values(const bcf_hdr_t* header, bcf1_t* record,
                                     const std::string& field,
                                     const RawInfoRecord* raw_info = nullptr) {
    if (raw_info != nullptr) {
        const auto raw = raw_info->values.find(field);
        if (raw != raw_info->values.end()) return raw->second;
    }
    char* text = nullptr;
    int text_size = 0;
    const auto string_length = bcf_get_info_string(header, record, field.c_str(), &text, &text_size);
    if (string_length >= 0 && text != nullptr) {
        const auto result = split_csv(std::string(text, static_cast<std::size_t>(string_length)));
        free(text);
        return result;
    }
    free(text);
    int32_t* integers = nullptr;
    int integer_count = 0;
    const auto integer_length = bcf_get_info_int32(header, record, field.c_str(), &integers, &integer_count);
    if (integer_length > 0 && integer_count > 0) {
        std::vector<std::string> result;
        result.reserve(static_cast<std::size_t>(integer_length));
        for (int index = 0; index < integer_length; ++index) result.push_back(format_int(integers[index]));
        free(integers);
        return result;
    }
    free(integers);
    float* reals = nullptr;
    int real_count = 0;
    const auto real_length = bcf_get_info_float(header, record, field.c_str(), &reals, &real_count);
    if (real_length > 0 && real_count > 0) {
        std::vector<std::string> result;
        result.reserve(static_cast<std::size_t>(real_length));
        for (int index = 0; index < real_length; ++index) result.push_back(format_float(reals[index]));
        free(reals);
        return result;
    }
    free(reals);
    int flag = 0;
    int flag_count = 0;
    if (bcf_get_info_flag(header, record, field.c_str(), &flag, &flag_count) > 0) return {"true"};
    return {};
}

std::vector<std::string> genotype_values(const bcf_hdr_t* header, bcf1_t* record,
                                         const std::string& field, int sample_index) {
    const int sample_count = bcf_hdr_nsamples(header);
    if (sample_index < 0 || sample_index >= sample_count) return {};
    if (field == "GT") {
        int32_t* genotypes = nullptr;
        int count = 0;
        const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
        if (length <= 0 || count < sample_count || count % sample_count != 0) {
            free(genotypes);
            return {};
        }
        const int ploidy = count / sample_count;
        std::ostringstream output;
        bool wrote = false;
        bool phased = false;
        for (int copy = 0; copy < ploidy; ++copy) {
            const auto encoded = genotypes[sample_index * ploidy + copy];
            if (encoded == bcf_int32_vector_end) break;
            if (wrote) output << (phased ? '|' : '/');
            if (encoded == bcf_int32_missing || bcf_gt_is_missing(encoded)) {
                output << '.';
            } else {
                const int allele = bcf_gt_allele(encoded);
                if (allele >= 0 && allele < record->n_allele && record->d.allele[allele] != nullptr)
                    output << record->d.allele[allele];
                else output << '.';
            }
            phased = bcf_gt_is_phased(encoded);
            wrote = true;
        }
        free(genotypes);
        return {wrote ? output.str() : kMissing};
    }
    const auto id = bcf_hdr_id2int(header, BCF_DT_ID, field.c_str());
    const int type = id >= 0 ? bcf_hdr_id2type(header, BCF_HL_FMT, id) : BCF_HT_INT;
    if (type == BCF_HT_STR) {
        char** values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_string(header, record, field.c_str(), &values, &count);
        std::vector<std::string> result;
        if (length > sample_index && values != nullptr && values[sample_index] != nullptr)
            result = split_csv(values[sample_index]);
        if (values != nullptr) {
            free(values[0]);
            free(values);
        }
        return result;
    }
    int32_t* integers = nullptr;
    int count = 0;
    const auto integer_length = bcf_get_format_int32(header, record, field.c_str(), &integers, &count);
    if (integer_length > 0 && count >= sample_count && count % sample_count == 0) {
        const int width = count / sample_count;
        std::vector<std::string> result;
        result.reserve(static_cast<std::size_t>(width));
        for (int index = 0; index < width; ++index)
            result.push_back(format_int(integers[sample_index * width + index]));
        free(integers);
        return result;
    }
    free(integers);
    float* reals = nullptr;
    count = 0;
    const auto real_length = bcf_get_format_float(header, record, field.c_str(), &reals, &count);
    if (real_length > 0 && count >= sample_count && count % sample_count == 0) {
        const int width = count / sample_count;
        std::vector<std::string> result;
        result.reserve(static_cast<std::size_t>(width));
        for (int index = 0; index < width; ++index)
            result.push_back(format_float(reals[sample_index * width + index]));
        free(reals);
        return result;
    }
    free(reals);
    return {};
}

std::string variant_type(bcf1_t* record) {
    // VariantContext.Type (used by GATK VariantsToTable) treats symbolic
    // alleles as SYMBOLIC even when htslib reports the record as VCF_REF.
    // Classify from the actual allele strings so <NON_REF>, spanning '*', and
    // breakend notation are not mistaken for reference-only records.
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele <= 1 || record->d.allele == nullptr || record->d.allele[0] == nullptr)
        return "NO_VARIATION";
    const std::string ref(record->d.allele[0]);
    bool snp = false;
    bool mnp = false;
    bool indel = false;
    bool symbolic = false;
    for (int index = 1; index < record->n_allele; ++index) {
        const char* alt_ptr = record->d.allele[index];
        // htsjdk's VariantContext treats the spanning-deletion allele `*`
        // as an ordinary one-base allele for type classification.  It is
        // therefore SNP-like when REF is one base (A>* and A>G,*), while it
        // still participates in an INDEL when another ALT changes length
        // (AT>A,*).  Symbolic `<...>` and breakend alleles are the only
        // categories that force a symbolic/mixed result here.
        if (alt_ptr == nullptr || *alt_ptr == '\0' || *alt_ptr == '<' ||
            std::strchr(alt_ptr, '[') != nullptr || std::strchr(alt_ptr, ']') != nullptr) {
            symbolic = true;
            continue;
        }
        const std::string alt(alt_ptr);
        if (ref.size() == 1 && alt.size() == 1) snp = true;
        else if (ref.size() == alt.size() && ref.size() > 1) mnp = true;
        else indel = true;
    }
    const int concrete = static_cast<int>(snp) + static_cast<int>(mnp) + static_cast<int>(indel);
    if (symbolic && concrete > 0) return "MIXED";
    if (symbolic) return "SYMBOLIC";
    if (concrete > 1) return "MIXED";
    if (snp) return "SNP";
    if (mnp) return "MNP";
    if (indel) return "INDEL";
    return "NO_VARIATION";
}

std::string format_qual(float value) {
    // HTSJDK represents a missing QUAL as VariantContext.NO_LOG10_PERROR,
    // whose table rendering is -10.0 rather than the generic missing token.
    // htslib stores a VCF '.' QUAL as bcf_float_missing, so translate that
    // sentinel before applying the normal float formatting path.
    if (bcf_float_is_missing(value) || bcf_float_is_vector_end(value)) return "-10.0";
    return format_float(value);
}

bool is_transition(bcf1_t* record) {
    if (record->n_allele != 2 || record->d.allele == nullptr || record->d.allele[0] == nullptr ||
        record->d.allele[1] == nullptr || std::strlen(record->d.allele[0]) != 1 ||
        std::strlen(record->d.allele[1]) != 1) return false;
    const char ref = static_cast<char>(std::toupper(record->d.allele[0][0]));
    const char alt = static_cast<char>(std::toupper(record->d.allele[1][0]));
    return (ref == 'A' && alt == 'G') || (ref == 'G' && alt == 'A') ||
           (ref == 'C' && alt == 'T') || (ref == 'T' && alt == 'C');
}

struct GenotypeSummary {
    std::uint64_t hom_ref = 0;
    std::uint64_t het = 0;
    std::uint64_t hom_var = 0;
    std::uint64_t no_call = 0;
};

GenotypeSummary summarize_genotypes(const bcf_hdr_t* header, bcf1_t* record) {
    GenotypeSummary summary;
    const int samples = bcf_hdr_nsamples(header);
    if (samples <= 0) return summary;
    int32_t* genotypes = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
    if (length <= 0 || count < samples || count % samples != 0) {
        free(genotypes);
        summary.no_call = static_cast<std::uint64_t>(samples);
        return summary;
    }
    const int ploidy = count / samples;
    for (int sample = 0; sample < samples; ++sample) {
        bool missing = false;
        bool any_alt = false;
        bool all_same = true;
        int first = -1;
        for (int copy = 0; copy < ploidy; ++copy) {
            const auto encoded = genotypes[sample * ploidy + copy];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                missing = true;
                continue;
            }
            const int allele = bcf_gt_allele(encoded);
            if (first < 0) first = allele;
            else if (first != allele) all_same = false;
            if (allele > 0) any_alt = true;
        }
        if (missing || first < 0) ++summary.no_call;
        else if (!any_alt) ++summary.hom_ref;
        else if (all_same) ++summary.hom_var;
        else ++summary.het;
    }
    free(genotypes);
    return summary;
}

std::string sample_name(const bcf_hdr_t* header, int index) {
    if (index < 0 || index >= bcf_hdr_nsamples(header)) return kMissing;
    const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
    return name == nullptr ? kMissing : std::string(name);
}

bool filtered(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && std::strcmp(name, "PASS") != 0 && std::strcmp(name, ".") != 0) return true;
    }
    return false;
}

std::vector<std::string> default_site_fields(const bcf_hdr_t* header) {
    std::vector<std::string> fields{"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"};
    const int count = header->n[BCF_DT_ID];
    for (int id = 0; id < count; ++id) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, id);
        if (name != nullptr && bcf_hdr_idinfo_exists(header, BCF_HL_INFO, id))
            fields.emplace_back(name);
    }
    return fields;
}

std::vector<std::string> default_genotype_fields(const bcf_hdr_t* header) {
    std::vector<std::string> fields;
    const int count = header->n[BCF_DT_ID];
    for (int id = 0; id < count; ++id) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, id);
        if (name != nullptr && bcf_hdr_idinfo_exists(header, BCF_HL_FMT, id)) {
            if (std::strcmp(name, "GT") == 0) fields.insert(fields.begin(), name);
            else fields.emplace_back(name);
        }
    }
    return fields;
}

std::string site_value(const bcf_hdr_t* header, bcf1_t* record, const std::string& raw_field,
                       int alt_index, int alt_count, bool split,
                       bool allele_specific = false,
                       const RawInfoRecord* raw_info = nullptr) {
    const auto field = trim(raw_field);
    const auto upper_field = upper(field);
    bcf_unpack(record, BCF_UN_STR | BCF_UN_FLT);
    if (!field.empty() && field.back() == '*') {
        // VariantsToTable wildcard fields are a single column.  GATK walks
        // present INFO attributes, converts multi-valued attributes through
        // Java List.toString() (`[a, b]`), sorts those rendered values, and
        // joins them with commas.  Do the same from HTSlib's typed INFO
        // records instead of expanding the wildcard into header columns.
        const auto prefix = field.substr(0, field.size() - 1);
        std::set<std::string> rendered;
        bcf_unpack(record, BCF_UN_INFO);
        for (int index = 0; index < record->n_info; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.info[index].key);
            if (name == nullptr || std::string(name).rfind(prefix, 0) != 0) continue;
            const auto values = info_values(header, record, name, raw_info);
            if (values.empty()) continue;
            if (values.size() == 1) rendered.insert(values.front());
            else rendered.insert("[" + join(values, ", ") + "]");
        }
        if (rendered.empty()) return kMissing;
        return join(std::vector<std::string>(rendered.begin(), rendered.end()));
    }
    if (upper_field == "CHROM") return bcf_hdr_id2name(header, record->rid) == nullptr ? kMissing : bcf_hdr_id2name(header, record->rid);
    if (upper_field == "POS") return std::to_string(record->pos + 1);
    if (upper_field == "ID") return record->d.id == nullptr ? "." : record->d.id;
    if (upper_field == "REF") return record->n_allele > 0 && record->d.allele[0] ? record->d.allele[0] : kMissing;
    if (upper_field == "ALT") {
        if (split && alt_index >= 0 && alt_index + 1 < record->n_allele) return record->d.allele[alt_index + 1];
        std::vector<std::string> alleles;
        for (int index = 1; index < record->n_allele; ++index) alleles.emplace_back(record->d.allele[index]);
        return alleles.empty() ? kMissing : join(alleles);
    }
    if (upper_field == "QUAL") return format_qual(record->qual);
    if (upper_field == "FILTER") {
        // VariantContext renders multiple filters with commas (the VCF
        // FILTER column itself is semicolon-delimited in the wire format).
        bcf_unpack(record, BCF_UN_FLT);
        if (record->d.n_flt == 0) return "PASS";
        std::ostringstream filters;
        for (int index = 0; index < record->d.n_flt; ++index) {
            if (index != 0) filters << ',';
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
            filters << (name == nullptr ? kMissing : name);
        }
        return filters.str();
    }
    if (upper_field == "TYPE") return variant_type(record);
    if (upper_field == "TRANSITION") {
        // GATK returns -1 for non-biallelic/non-SNP records, not false.
        if (record->n_allele != 2 || variant_type(record) != "SNP") return "-1";
        return is_transition(record) ? "1" : "0";
    }
    // VariantsToTable exposes TRANSITION but does not define TRANSVERSION,
    // INDEL, or N_ALLELES as standard table fields.  They consequently flow
    // through the regular INFO lookup and render as NA, matching GATK's
    // Field extraction contract instead of inventing native-only columns.
    if (upper_field == "MULTI-ALLELIC" || upper_field == "MULTIALLELIC") return record->n_allele > 2 ? "true" : "false";
    if (upper_field == "HET" || upper_field == "HOM-REF" || upper_field == "HOM-VAR" ||
        upper_field == "NO-CALL" || upper_field == "VAR" || upper_field == "NSAMPLES" ||
        upper_field == "NCALLED" || upper_field == "SAMPLE_NAME") {
        const auto summary = summarize_genotypes(header, record);
        if (upper_field == "HET") return std::to_string(summary.het);
        if (upper_field == "HOM-REF") return std::to_string(summary.hom_ref);
        if (upper_field == "HOM-VAR") return std::to_string(summary.hom_var);
        if (upper_field == "NO-CALL") return std::to_string(summary.no_call);
        if (upper_field == "VAR") return std::to_string(summary.het + summary.hom_var);
        if (upper_field == "NSAMPLES") return std::to_string(bcf_hdr_nsamples(header));
        if (upper_field == "NCALLED") {
            return std::to_string(static_cast<std::uint64_t>(bcf_hdr_nsamples(header)) -
                                  std::min<std::uint64_t>(summary.no_call, bcf_hdr_nsamples(header)));
        }
        return sample_name(header, 0);
    }
    if (upper_field == "EVENTLENGTH") {
        if (record->n_allele <= 1 || record->d.allele[0] == nullptr) return "0";
        const auto ref_length = static_cast<long long>(std::strlen(record->d.allele[0]));
        long long maximum = 0;
        bool saw_symbolic = false;
        for (int index = 1; index < record->n_allele; ++index) {
            if (record->d.allele[index] == nullptr) continue;
            const char* alt = record->d.allele[index];
            if (*alt == '\0' || *alt == '<' ||
                std::strchr(alt, '[') != nullptr || std::strchr(alt, ']') != nullptr) {
                saw_symbolic = true;
                continue;
            }
            const auto delta = static_cast<long long>(std::strlen(record->d.allele[index])) - ref_length;
            if (std::llabs(delta) > std::llabs(maximum)) maximum = delta;
        }
        if (saw_symbolic && maximum == 0) return "-1";
        return std::to_string(maximum);
    }
    // GATK accepts bare INFO keys (for example `DP`), not the htslib-style
    // `INFO/DP` namespace prefix.  Keep the requested spelling intact so an
    // explicit prefix follows the normal missing-field path and renders NA.
    const auto& info_field = field;
    auto values = info_values(header, record, info_field, raw_info);
    if (values.empty()) return kMissing;
    if (!split || alt_index < 0) return join(values);
    // Java's ordinary -F path uses addFieldValue(): a List is split only
    // when its actual length equals the number of ALT rows.  In particular a
    // Number=R INFO value stays as ref,alt1,alt2 on every row.  -ASF is
    // different only for Number=R, where the ref element is first removed;
    // the remaining ALT list is then distributed one value per row.
    if (allele_specific) {
        // GATK's split -ASF implementation obtains the value through
        // VariantContext.getAttributeAsString(), whose List rendering is
        // `[first, second, ...]`, then only removes brackets before splitting
        // on commas.  The space before every non-first element is therefore
        // observable output and must survive the native typed decoder.
        for (std::size_t index = 1; index < values.size(); ++index)
            values[index] = " " + values[index];
    }
    if (allele_specific && number_code(header, "INFO", info_field) == BCF_VL_R) {
        if (values.size() > 1) values.erase(values.begin());
        else values = {kMissing};
    }
    if (values.size() == static_cast<std::size_t>(alt_count))
        return values[static_cast<std::size_t>(alt_index)];
    return join(values);
}

std::string genotype_value(const bcf_hdr_t* header, bcf1_t* record, const std::string& raw_field,
                           int sample_index, int alt_index, int alt_count, bool split,
                           bool allele_specific) {
    std::string field = trim(raw_field);
    auto values = genotype_values(header, record, field, sample_index);
    if (values.empty()) return kMissing;
    if (split && allele_specific) {
        const auto length_code = number_code(header, "FORMAT", field);
        if (length_code == BCF_VL_R && field != "AD") {
            // GATK's ASGF path drops the reference entry for generic
            // Number=R fields; AD is the one special case rendered as
            // ref,selected-alt.
            if (alt_index >= 0 && alt_index + 1 < static_cast<int>(values.size()))
                values = {values[static_cast<std::size_t>(alt_index + 1)]};
            else values = {kMissing};
        } else if (length_code == BCF_VL_A || (length_code == BCF_VL_R && field == "AD")) {
            values = project_values(std::move(values), length_code, alt_index, alt_count);
        }
        // Number=G values are deliberately not projected.  This matches
        // VariantsToTable's addAlleleSpecificFieldValue: only A/R fields are
        // split; genotype likelihood vectors remain intact on every ALT row.
    }
    return join(values);
}

// `NA` is the rendering for an absent value, but it is also a legal literal
// String INFO/FORMAT value.  Keep presence separate from rendering so -EMD
// follows VariantsToTable's hasAttribute()/hasAnyAttribute() checks rather
// than rejecting a genuine value whose spelling happens to be "NA".
bool is_standard_site_field(const std::string& upper_field) {
    return upper_field == "CHROM" || upper_field == "POS" || upper_field == "ID" ||
           upper_field == "REF" || upper_field == "ALT" || upper_field == "QUAL" ||
           upper_field == "FILTER" || upper_field == "TYPE" ||
           upper_field == "TRANSITION" || upper_field == "MULTI-ALLELIC" ||
           upper_field == "MULTIALLELIC" || upper_field == "HET" ||
           upper_field == "HOM-REF" || upper_field == "HOM-VAR" ||
           upper_field == "NO-CALL" || upper_field == "VAR" ||
           upper_field == "NSAMPLES" || upper_field == "NCALLED" ||
           upper_field == "SAMPLE_NAME" || upper_field == "EVENTLENGTH";
}

bool site_value_missing(const bcf_hdr_t* header, bcf1_t* record, const std::string& raw_field,
                        const RawInfoRecord* raw_info = nullptr) {
    const auto field = trim(raw_field);
    const auto upper_field = upper(field);
    if (is_standard_site_field(upper_field)) return false;
    bcf_unpack(record, BCF_UN_INFO);
    if (!field.empty() && field.back() == '*') {
        const auto prefix = field.substr(0, field.size() - 1);
        for (int index = 0; index < record->n_info; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.info[index].key);
            if (name != nullptr && std::string(name).rfind(prefix, 0) == 0) return false;
        }
        return true;
    }
    return info_values(header, record, field, raw_info).empty();
}

bool genotype_value_missing(const bcf_hdr_t* header, bcf1_t* record,
                            const std::string& field, int sample_index) {
    return genotype_values(header, record, trim(field), sample_index).empty();
}

struct Sample {
    std::string name;
    int index = -1;
};

struct GenotypeField {
    std::string name;
    bool allele_specific = false;
};

struct TableCell {
    std::string sample;
    std::string field;
    std::string value;
    bool missing = false;
    // GATK's molten writer intentionally iterates fieldsToTake and
    // genotypeFieldsToTake only.  -ASF/-ASGF values are extracted for normal
    // table rows, but are not emitted in the four-column molten projection.
    bool molten_visible = true;
};

std::string table_sample_name(const char* name) {
    if (name == nullptr) return kMissing;
    std::string result(name);
    std::replace(result.begin(), result.end(), ' ', '_');
    return result;
}

std::vector<Sample> sorted_samples(const bcf_hdr_t* header) {
    std::vector<Sample> samples;
    for (int index = 0; index < bcf_hdr_nsamples(header); ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
        if (name != nullptr) samples.push_back(Sample{table_sample_name(name), index});
    }
    std::sort(samples.begin(), samples.end(), [](const Sample& left, const Sample& right) {
        return left.name < right.name;
    });
    return samples;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    htsFile* input = nullptr;
    bcf_hdr_t* header = nullptr;
    bcf1_t* record = nullptr;
    std::ofstream file_output;
    std::ostream* output = nullptr;
    try {
        input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/BCF: " + options.input);
        header = bcf_hdr_read(input);
        if (!header) throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + options.input);
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate VCF record");

        auto site_fields = options.fields.empty() && options.genotype_fields.empty() &&
                           options.allele_fields.empty() && options.allele_genotype_fields.empty()
            ? default_site_fields(header) : options.fields;
        auto genotype_fields = options.fields.empty() && options.genotype_fields.empty() &&
                               options.allele_fields.empty() && options.allele_genotype_fields.empty()
            ? default_genotype_fields(header) : options.genotype_fields;
        // Keep wildcard fields as a single requested column.  GATK's
        // VariantsToTable renders `PREFIX*` from present INFO attributes;
        // expanding against the header would change both the column contract
        // and missing-record behavior.
        const auto allele_fields = options.allele_fields;
        const auto allele_genotype_fields = options.allele_genotype_fields;
        site_fields.insert(site_fields.end(), allele_fields.begin(), allele_fields.end());
        std::vector<GenotypeField> genotype_specs;
        genotype_specs.reserve(genotype_fields.size() + allele_genotype_fields.size());
        for (const auto& field : genotype_fields) genotype_specs.push_back(GenotypeField{field, false});
        for (const auto& field : allele_genotype_fields) genotype_specs.push_back(GenotypeField{field, true});
        const auto samples = sorted_samples(header);
        if (!genotype_specs.empty() && samples.empty())
            throw std::invalid_argument("BAD_INPUT: genotype fields requested but VCF has no samples");
        RawInfoReader raw_info_reader(options.input);

        if (options.output == "-") output = &std::cout;
        else {
            file_output.open(options.output);
            if (!file_output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create table: " + options.output);
            output = &file_output;
        }
        if (options.moltenize) {
            *output << "RecordID\tSample\tVariable\tValue\n";
        } else {
            bool first = true;
            for (const auto& field : site_fields) {
                if (!first) *output << '\t';
                *output << field;
                first = false;
            }
            for (const auto& sample : samples) {
                for (const auto& field : genotype_specs) {
                    if (!first) *output << '\t';
                    *output << sample.name << '.' << field.name;
                    first = false;
                }
            }
            *output << '\n';
        }

        std::vector<fastgatk::io::IndexedInterval> intervals;
        fastgatk::io::IntervalFileStats interval_stats;
        bool first_selector = true;
        for (const auto& selector : options.regions) {
            append_table_interval_selector_with_rule(
                selector, header, intervals, interval_stats,
                options.interval_set_rule, first_selector,
                options.merge_abutting_intervals);
            first_selector = false;
        }
        normalize_table_intervals(intervals, options.merge_abutting_intervals);
        // GATK's IntervalFilteringVcfWriter requires a non-empty interval
        // set for every mode other than ANYWHERE.  Keep the failure explicit
        // instead of silently treating a misspelled/under-specified selector
        // as a whole-file extraction.
        if (options.variant_output_filtering_explicit && intervals.empty() &&
            options.variant_output_filtering != Options::VariantOutputFiltering::Anywhere)
            throw std::invalid_argument(
                "--variant-output-filtering requires at least one -L/--intervals selector");
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
        std::uint64_t input_records = 0;
        std::uint64_t output_rows = 0;
        std::uint64_t table_record_id = 0;
        std::uint64_t skipped_filtered = 0;
        std::uint64_t skipped_interval = 0;
        std::uint64_t skipped_excluded = 0;
        while (bcf_read(input, header, record) == 0) {
            ++input_records;
            RawInfoRecord raw_info_record;
            const RawInfoRecord* raw_info = raw_info_reader.next(raw_info_record)
                ? &raw_info_record : nullptr;
            bcf_unpack(record, BCF_UN_ALL);
            if (!intervals.empty()) {
                const auto record_end = record_end_exclusive(header, record);
                // VariantsToTable writes a text table directly and therefore
                // never constructs GATK's IntervalFilteringVcfWriter (that
                // decorator is used only by VCF-producing tools).  In GATK,
                // --variant-output-filtering is consequently an inherited,
                // validation-only option here: traversal remains the normal
                // OVERLAPS query for every mode, including ANYWHERE.
                const bool selected = std::any_of(intervals.begin(), intervals.end(), [&](const auto& interval) {
                    return record->rid == interval.rid && record_end > interval.begin &&
                           record->pos < interval.end;
                });
                if (!selected) {
                    ++skipped_interval;
                    bcf_clear(record);
                    continue;
                }
            }
            if (!exclusions.empty()) {
                const auto record_end = record_end_exclusive(header, record);
                const bool excluded = std::any_of(exclusions.begin(), exclusions.end(), [&](const auto& interval) {
                    return record->rid == interval.rid && record_end > interval.begin &&
                           record->pos < interval.end;
                });
                if (excluded) {
                    ++skipped_excluded;
                    bcf_clear(record);
                    continue;
                }
            }
            if (!options.show_filtered && filtered(header, record)) {
                ++skipped_filtered;
                bcf_clear(record);
                continue;
            }
            const int alt_count = std::max(0, record->n_allele - 1);
            const int row_count = options.split_multi_allelic && alt_count > 0 ? alt_count : 1;
            const std::uint64_t record_id = ++table_record_id;
            for (int row = 0; row < row_count; ++row) {
                const int alt_index = options.split_multi_allelic && alt_count > 0 ? row : -1;
                std::vector<TableCell> cells;
                cells.reserve(site_fields.size() + samples.size() * genotype_specs.size());
                const std::size_t ordinary_site_fields = site_fields.size() - allele_fields.size();
                for (std::size_t index = 0; index < site_fields.size(); ++index) {
                    const auto& field = site_fields[index];
                    cells.push_back(TableCell{
                        "site", field,
                        site_value(header, record, field, alt_index, alt_count, alt_index >= 0,
                                   index >= ordinary_site_fields, raw_info),
                        site_value_missing(header, record, field, raw_info),
                        index < ordinary_site_fields,
                    });
                }
                for (const auto& sample : samples) {
                    for (const auto& field : genotype_specs) {
                        cells.push_back(TableCell{
                            sample.name, field.name,
                            genotype_value(header, record, field.name, sample.index, alt_index, alt_count,
                                           alt_index >= 0, field.allele_specific),
                            genotype_value_missing(header, record, field.name, sample.index),
                            !field.allele_specific,
                        });
                    }
                }
                for (const auto& cell : cells) {
                    if (cell.missing && options.error_if_missing)
                        throw std::runtime_error("BAD_INPUT: missing VariantsToTable value for " + cell.field +
                                                 " in sample " + cell.sample);
                }
                if (options.moltenize) {
                    // Preserve GATK 4.6.2.0's public molten writer exactly:
                    // extraction first appends -F, -ASF, -GF, then -ASGF
                    // values, while emitMoltenizedOutput iterates only -F
                    // and -GF and consumes the extracted list sequentially.
                    // Thus an -ASF/-ASGF request shifts subsequent emitted
                    // values instead of producing its own molten rows.  It
                    // is quirky, but pipelines relying on direct replacement
                    // observe this concrete behavior.
                    std::size_t emission_index = 0;
                    for (std::size_t index = 0; index < ordinary_site_fields; ++index) {
                        const auto& cell = cells.at(emission_index++);
                        *output << record_id << "\tsite\t" << site_fields[index] << '\t'
                                << cell.value << '\n';
                    }
                    for (const auto& sample : samples) {
                        for (const auto& field : genotype_specs) {
                            if (field.allele_specific) continue;
                            const auto& cell = cells.at(emission_index++);
                            *output << record_id << '\t' << sample.name << '\t' << field.name << '\t'
                                    << cell.value << '\n';
                        }
                    }
                } else {
                    for (std::size_t index = 0; index < cells.size(); ++index) {
                        if (index != 0) *output << '\t';
                        *output << cells[index].value;
                    }
                    *output << '\n';
                }
                ++output_rows;
            }
            // For split multi-allelic records all rows share the original
            // RecordID, matching GATK's table contract.  The count itself is
            // tracked separately from output_rows below.
            bcf_clear(record);
        }
        if (options.output != "-") file_output.flush();
        bcf_destroy(record); record = nullptr;
        bcf_hdr_destroy(header); header = nullptr;
        bcf_close(input); input = nullptr;
        const auto primary_complete = options.output == "-" || file_complete(options.output);
        if (!primary_complete) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: table is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        if (options.output != "-") {
            std::ofstream manifest(manifest_path);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create manifest: " + manifest_path);
            manifest << "{\"schema_version\":1,\"tool\":\"VariantsToTable\",\"implementation\":\"fastgatk-variants-to-table\",\"status\":\"prototype\","
                     << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"table\","
                     << "\"compatibility\":{\"field_extraction\":true,\"sample_name_sort\":true,\"missing_values_na\":true,"
                     << "\"show_filtered\":" << (options.show_filtered ? "true" : "false")
                     << ",\"split_multi_allelic\":" << (options.split_multi_allelic ? "true" : "false")
                     << ",\"moltenize\":" << (options.moltenize ? "true" : "false")
                     << ",\"interval_subset\":" << (!intervals.empty() ? "true" : "false")
                     << ",\"interval_set_rule\":\""
                     << interval_set_rule_name(options.interval_set_rule) << "\""
                     << ",\"interval_merging_rule\":\""
                     << (options.merge_abutting_intervals ? "ALL" : "OVERLAPPING_ONLY") << "\""
                     << ",\"variant_output_filtering\":\""
                     << variant_output_filtering_name(options.variant_output_filtering) << "\""
                     << ",\"interval_list_inputs\":" << interval_stats.files
                     << ",\"interval_list_records\":" << interval_stats.records
                     << ",\"interval_padding\":" << options.interval_padding
                     << ",\"interval_exclusion\":" << (!exclusions.empty() ? "true" : "false")
                     << ",\"interval_exclusion_list_inputs\":" << exclusion_stats.files
                     << ",\"interval_exclusion_list_records\":" << exclusion_stats.records
                     << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                     << ",\"cloud_uri_staging\":false,\"raw_info_tokens\":"
                     << (raw_info_reader.enabled() ? "true" : "false")
                     << ",\"bit_identical_to_gatk\":false},"
                     << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"table\",\"complete\":" << (primary_complete ? "true" : "false")
                     << "}],\"telemetry\":{\"resources\":" << resources.to_json()
                     << ",\"input_records\":" << input_records << ",\"output_rows\":" << output_rows
                     << ",\"skipped_filtered\":" << skipped_filtered
                     << ",\"skipped_interval\":" << skipped_interval
                     << ",\"interval_set_rule\":\""
                     << interval_set_rule_name(options.interval_set_rule) << "\""
                     << ",\"skipped_excluded\":" << skipped_excluded << "}}\n";
            manifest.flush();
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize manifest: " + manifest_path);
            manifest.close();
            fastgatk::runtime::require_complete_output({
                std::filesystem::path(options.output), {}, std::filesystem::path(manifest_path), false, true});
        }
        std::cout << "{\"tool\":\"VariantsToTable\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"output_rows\":" << output_rows
                  << ",\"skipped_excluded\":" << skipped_excluded
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\",\"interval_merging_rule\":\""
                  << (options.merge_abutting_intervals ? "ALL" : "OVERLAPPING_ONLY")
                  << "\",\"variant_output_filtering\":\""
                  << variant_output_filtering_name(options.variant_output_filtering) << "\"}\n";
        return 0;
    } catch (...) {
        if (record) bcf_destroy(record);
        if (header) bcf_hdr_destroy(header);
        if (input) bcf_close(input);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for VariantsToTable");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-variants-to-table: " << error.what() << '\n';
        return 2;
    }
}
