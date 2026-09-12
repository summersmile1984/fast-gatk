#include <algorithm>
#include <chrono>
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
#include <sstream>
#include <stdexcept>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "optional_boolean.hpp"

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string recal_file;
    std::string tranches_file;
    std::string output;
    std::string manifest;
    std::string mode = "SNP";
    double truth_sensitivity = std::numeric_limits<double>::quiet_NaN();
    double lod_cutoff = std::numeric_limits<double>::quiet_NaN();
    bool ignore_all_filters = false;
    bool exclude_filtered = false;
    bool allele_specific = false;
    // GATK's output index switch is an optional boolean and defaults to true.
    // Keep the choice in the output contract so a false request never leaves
    // an apparently valid but stale Tabix sidecar behind.
    bool create_output_variant_index = true;
    // GATK's writer-only switch: retain all INFO/FILTER recalibration work,
    // but omit FORMAT/sample columns from the published VCF.  Keep the full
    // input header active until the final writer boundary.
    bool sites_only_vcf_output = false;
    std::vector<std::string> ignore_filters;
    std::vector<std::string> regions;
    std::vector<std::string> excluded_regions;
    int exclusion_padding = 0;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
};

struct Tranche {
    double target = 0.0;
    double min_vqslod = 0.0;
    std::string filter_name;
    std::string model;
};

std::string tranche_filter_for_score(const std::vector<Tranche>& tranches, double score) {
    // GATK reverses the tranches selected by --truth-sensitivity-filter-level
    // (best/highest sensitivity first), then walks from the least sensitive
    // tranche toward the best score.  Preserve that exact boundary behavior:
    // a score below the weakest selected tranche receives the `+` filter,
    // scores between tranche minima receive the preceding tranche name, and a
    // score at/above the best selected tranche passes.
    if (tranches.empty()) return "LOW_VQSLOD";
    for (std::size_t index = tranches.size(); index-- > 0;) {
        if (score >= tranches[index].min_vqslod)
            return index == tranches.size() - 1 ? "PASS" : tranches[index].filter_name;
    }
    return tranches.front().filter_name + "+";
}

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = [&]() {
        const auto primary = option_value(argument, name);
        if (!primary.empty() || short_name == nullptr) return primary;
        return option_value(argument, short_name);
    }();
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find(',', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
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

bool complete_file(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() && path.compare(path.size() - value.size(), value.size(), value) == 0;
}

std::uintmax_t file_size_or_zero(const std::string& path) {
    if (path.empty()) return 0;
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : size;
}

std::int64_t file_mtime_ticks_or_zero(const std::string& path) {
    if (path.empty()) return 0;
    std::error_code error;
    const auto time = std::filesystem::last_write_time(path, error);
    if (error) return 0;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        time.time_since_epoch()).count();
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-apply-vqsr (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE                 input VCF/BCF\n"
                         "      --recal-file FILE              recal VCF containing INFO/VQSLOD\n"
                         "      --tranches-file FILE           GATK tranches CSV\n"
                         "      --truth-sensitivity-filter-level F  tranche target in percent\n"
                         "      --lod-score-cutoff F            explicit VQSLOD cutoff\n"
                         "      --mode SNP|INDEL|BOTH           recalibration mode\n"
                         "  -L, --intervals REGION             interval selector (repeatable)\n"
                         "      --interval-set-rule RULE       UNION (default) or INTERSECTION\n"
                         "      --use-allele-specific-annotations[=true|false]  AS mode\n"
                         "      --create-output-variant-index[=true|false]     write .tbi (default true)\n"
                         "      --sites-only-vcf-output[=true|false]           omit FORMAT/sample columns\n"
                         "  -XL, --exclude-intervals REGION    exclude loci after -L selection (repeatable)\n"
                         "      --interval-exclusion-padding N pad excluded intervals (default 0)\n"
                         "      --ignore-filter NAME            apply to records carrying NAME\n"
                         "  -O, --output FILE                  filtered VCF/BCF\n"
                         "      --exclude-filtered             omit records below cutoff\n"
                         "      --output-manifest FILE         OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-V" || has_option(argument, "--variant")) {
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        } else if (has_option(argument, "--recal-file")) {
            options.recal_file = require_value(index, argc, argv, argument, "--recal-file");
        } else if (has_option(argument, "--tranches-file")) {
            options.tranches_file = require_value(index, argc, argv, argument, "--tranches-file");
        } else if (has_option(argument, "--truth-sensitivity-filter-level") ||
                   has_option(argument, "--ts-filter-level")) {
            const auto name = argument.rfind("--ts-filter-level", 0) == 0
                ? "--ts-filter-level" : "--truth-sensitivity-filter-level";
            options.truth_sensitivity = std::stod(require_value(index, argc, argv, argument, name));
        } else if (has_option(argument, "--lod-score-cutoff")) {
            options.lod_cutoff = std::stod(require_value(index, argc, argv, argument, "--lod-score-cutoff"));
        } else if (has_option(argument, "--mode")) {
            options.mode = require_value(index, argc, argv, argument, "--mode");
            if (options.mode != "SNP" && options.mode != "INDEL" && options.mode != "BOTH")
                throw std::invalid_argument("--mode must be SNP, INDEL, or BOTH");
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                (has_option(argument, "--region") ? "--region" :
                 (has_option(argument, "--intervals") ? "--intervals" : "--interval"));
            options.regions.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals") ||
                   has_option(argument, "--exclude-interval") || has_option(argument, "--exclude-region")) {
            const char* name = argument == "-XL" ? "--exclude-intervals" :
                (has_option(argument, "--exclude-region") ? "--exclude-region" :
                 (has_option(argument, "--exclude-intervals") ? "--exclude-intervals" : "--exclude-interval"));
            options.excluded_regions.push_back(require_value(index, argc, argv, argument, name, "-XL"));
        } else if (argument == "-ixp" || has_option(argument, "--interval-exclusion-padding")) {
            const auto value = require_value(index, argc, argv, argument,
                                             "--interval-exclusion-padding", "-ixp");
            options.exclusion_padding = std::stoi(value);
            if (options.exclusion_padding < 0)
                throw std::invalid_argument("interval-exclusion-padding must be non-negative");
        } else if (has_option(argument, "--interval-set-rule") || has_option(argument, "-isr")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
                return static_cast<char>(std::toupper(character));
            });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--use-allele-specific-annotations" ||
                   argument.rfind("--use-allele-specific-annotations=", 0) == 0 ||
                   argument == "--AS" || argument.rfind("--AS=", 0) == 0) {
            const char* name = argument.rfind("--AS", 0) == 0
                ? "--AS" : "--use-allele-specific-annotations";
            options.allele_specific = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (argument == "--create-output-variant-index" ||
                   argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_output_variant_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--sites-only-vcf-output" ||
                   argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--sites-only-vcf-output");
        } else if (argument == "--ignore-all-filters" ||
                   argument.rfind("--ignore-all-filters=", 0) == 0) {
            // Barclay exposes this as an optional Boolean, not a presence
            // flag: the separated false form must preserve the input FILTER
            // just like its default.  Keep it on the shared parser used by
            // the other writer switches so malformed values fail closed.
            options.ignore_all_filters = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--ignore-all-filters");
        } else if (argument == "--exclude-filtered" ||
                   argument.rfind("--exclude-filtered=", 0) == 0) {
            options.exclude_filtered = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--exclude-filtered");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility switches.
        } else if (has_option(argument, "--ignore-filter")) {
            options.ignore_filters.push_back(require_value(index, argc, argv, argument, "--ignore-filter"));
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.recal_file.empty()) throw std::invalid_argument("--recal-file is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (std::isfinite(options.truth_sensitivity) &&
        (!std::isfinite(options.truth_sensitivity) || options.truth_sensitivity < 0.0 || options.truth_sensitivity > 100.0))
        throw std::invalid_argument("truth-sensitivity-filter-level must be in [0,100]");
    if (!std::isfinite(options.lod_cutoff) && !options.tranches_file.empty()) return options;
    return options;
}

std::map<std::string, std::size_t> csv_indices(const std::vector<std::string>& header) {
    std::map<std::string, std::size_t> result;
    for (std::size_t index = 0; index < header.size(); ++index) result[header[index]] = index;
    return result;
}

template <typename Fn>
void for_each_tranche_line(const std::string& path, Fn&& callback) {
#if FASTGATK_HAS_HTSLIB
    if (suffix(path, ".gz")) {
        htsFile* input = hts_open(path.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open compressed tranches file: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (fastgatk::io::read_text_line(input, &line, path) >= 0)
                callback(line.s == nullptr ? std::string{} : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(input);
            throw;
        }
        const auto close_status = hts_close(input);
        free(line.s);
        if (close_status != 0)
            throw std::runtime_error("BAD_INPUT: cannot read compressed tranches file: " + path);
        return;
    }
#endif
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open tranches file: " + path);
    for (std::string line; std::getline(input, line);) callback(std::move(line));
    if (input.bad()) throw std::runtime_error("BAD_INPUT: cannot read tranches file: " + path);
}

std::vector<Tranche> read_tranches(const std::string& path, const std::string& mode) {
    std::vector<std::string> header;
    std::vector<Tranche> rows;
    for_each_tranche_line(path, [&](std::string line) {
        if (line.empty() || line[0] == '#') return;
        const auto fields = split_csv(line);
        if (header.empty()) { header = fields; return; }
        const auto index = csv_indices(header);
        for (const auto* name : {"targetTruthSensitivity", "minVQSLod", "filterName", "model"})
            if (!index.count(name)) throw std::runtime_error("BAD_INPUT: tranches header missing " + std::string(name));
        if (fields.size() != header.size()) throw std::runtime_error("BAD_INPUT: malformed tranches row");
        Tranche row;
        row.target = std::stod(fields.at(index.at("targetTruthSensitivity")));
        row.min_vqslod = std::stod(fields.at(index.at("minVQSLod")));
        row.filter_name = fields.at(index.at("filterName"));
        row.model = fields.at(index.at("model"));
        if (mode == "BOTH" || row.model == mode) rows.push_back(std::move(row));
    });
    if (rows.empty()) throw std::runtime_error("BAD_INPUT: no tranches match mode " + mode);
    return rows;
}

#if FASTGATK_HAS_HTSLIB

std::string variant_key(const bcf_hdr_t* header, bcf1_t* record, int alt_index) {
    bcf_unpack(record, BCF_UN_STR);
    const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    if (!contig || !record->d.allele || !record->d.allele[0] ||
        alt_index < 1 || alt_index >= record->n_allele || !record->d.allele[alt_index])
        return {};
    return std::string(contig) + ":" + std::to_string(record->pos + 1) + ":" +
           record->d.allele[0] + ":" + record->d.allele[alt_index];
}

std::string variant_type(bcf1_t* record) {
    const auto types = bcf_get_variant_types(record);
    const bool snp = (types & VCF_SNP) != 0;
    const bool mnp = (types & VCF_MNP) != 0;
    const bool indel = (types & VCF_INDEL) != 0;
    if ((snp || mnp) && indel) return "MIXED";
    if (indel) return "INDEL";
    if (snp) return "SNP";
    if (mnp) return "MNP";
    bcf_unpack(record, BCF_UN_STR);
    if (record->d.allele != nullptr)
        for (int allele = 1; allele < record->n_allele; ++allele)
            if (record->d.allele[allele] != nullptr && record->d.allele[allele][0] == '<')
                return "SYMBOLIC";
    return "OTHER";
}

bool record_matches_mode(const std::string& type, const std::string& mode) {
    if (mode == "BOTH") return true;
    if (mode == "SNP") return type == "SNP" || type == "MNP";
    if (mode == "INDEL") return type == "INDEL" || type == "MIXED" || type == "SYMBOLIC";
    return false;
}

bool allele_matches_mode(const bcf1_t* record, int alt_index, const std::string& mode) {
    if (mode == "BOTH") return true;
    if (!record || !record->d.allele || !record->d.allele[0] ||
        alt_index < 1 || alt_index >= record->n_allele || !record->d.allele[alt_index])
        return false;
    const std::string ref(record->d.allele[0]);
    const std::string alt(record->d.allele[alt_index]);
    if (alt.empty()) return false;
    // VariantDataManager.checkVariationClass is allele based: a spanning
    // deletion '*' has the same one-base length as a one-base REF and is
    // therefore considered SNP, while symbolic structural alleles belong to
    // INDEL mode.  Do not discard both categories before applying that rule.
    const bool symbolic_structural = alt.front() == '<';
    const bool snp = !symbolic_structural && ref.size() == alt.size() && ref.size() == 1;
    const bool indel = symbolic_structural || ref.size() != alt.size();
    return mode == "SNP" ? snp : (mode == "INDEL" ? indel : false);
}

struct RecalScores {
    std::unordered_map<std::string, float> scores;
    std::unordered_map<std::string, std::string> culprits;
    std::unordered_map<std::string, bool> positive_labels;
    std::unordered_map<std::string, bool> negative_labels;
    std::unordered_map<std::string, float> site_scores;
    std::unordered_map<std::string, std::string> site_culprits;
    std::unordered_map<std::string, bool> site_positive_labels;
    std::unordered_map<std::string, bool> site_negative_labels;
    bool allele_specific = false;
};

std::string variant_location_key(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    if (!contig || !record->d.allele || !record->d.allele[0]) return {};
    int32_t* end = nullptr;
    int end_count = 0;
    const int end_length = bcf_get_info_int32(header, record, "END", &end, &end_count);
    int64_t end_position = static_cast<int64_t>(record->pos) + 1 +
        static_cast<int64_t>(std::max<std::size_t>(1, std::strlen(record->d.allele[0]))) - 1;
    if (end_length > 0 && end != nullptr && end_count > 0 && end[0] > 0)
        end_position = end[0];
    free(end);
    return std::string(contig) + ":" + std::to_string(record->pos + 1) + ":" +
           std::to_string(end_position);
}

RecalScores read_recal_vcf(const std::string& path, bool require_allele_specific) {
    htsFile* input = bcf_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open recal VCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    bcf1_t* record = bcf_init();
    if (!header || !record) throw std::runtime_error("BAD_INPUT: cannot read recal VCF header");
    RecalScores result;
    while (fastgatk::io::read_variant_record(input, header, record, path) == 0) {
        bcf_unpack(record, BCF_UN_ALL);
        float* values = nullptr;
        int count = 0;
        const char* info_name = require_allele_specific ? "AS_VQSLOD" : "VQSLOD";
        const int length = bcf_get_info_float(header, record, info_name, &values, &count);
        if (length <= 0 || values == nullptr) {
            free(values);
            bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
            throw std::runtime_error(std::string("BACKEND_UNAVAILABLE: recal VCF must contain INFO/") + info_name);
        }
        char* as_culprit_value = nullptr;
        int as_culprit_count = 0;
        const int as_culprit_length = bcf_get_info_string(header, record, "AS_culprit",
                                                          &as_culprit_value, &as_culprit_count);
        char* site_culprit = nullptr;
        int site_culprit_count = 0;
        char* site_culprit_upper = nullptr;
        int site_culprit_upper_count = 0;
        if (as_culprit_length <= 0) {
            // GATK's recalibration VCF uses lower-case INFO/culprit.  Accept
            // the historical upper-case CULPRIT spelling emitted by older
            // native reports as an input compatibility fallback.
            const int lower_length = bcf_get_info_string(
                header, record, "culprit", &site_culprit, &site_culprit_count);
            if (lower_length <= 0)
                (void)bcf_get_info_string(header, record, "CULPRIT",
                                          &site_culprit_upper, &site_culprit_upper_count);
        }
        int flag_count = 0;
        const int positive_flag = bcf_get_info_flag(header, record, "POSITIVE_TRAIN_SITE", nullptr, &flag_count);
        flag_count = 0;
        const int negative_flag = bcf_get_info_flag(header, record, "NEGATIVE_TRAIN_SITE", nullptr, &flag_count);
        const bool positive = positive_flag > 0;
        const bool negative = negative_flag > 0;
        const auto location = variant_location_key(header, record);
        if (!location.empty() && length > 0) {
            result.site_scores[location] = values[0];
            if (as_culprit_length <= 0) {
                if (site_culprit != nullptr)
                    result.site_culprits[location] = site_culprit;
                else if (site_culprit_upper != nullptr)
                    result.site_culprits[location] = site_culprit_upper;
            }
            result.site_positive_labels[location] = positive;
            result.site_negative_labels[location] = negative;
        }
        for (int alt = 1; alt < record->n_allele; ++alt) {
            const auto key = variant_key(header, record, alt);
            if (key.empty()) continue;
            result.scores[key] = values[std::min(alt - 1, length - 1)];
            if (as_culprit_value != nullptr && as_culprit_length > 0) {
                // Older native AS fixtures used a Number=A AS_culprit vector.
                // Keep accepting that spelling, although GATK's actual AS
                // recalibration VCF stores one scalar culprit on each
                // per-ALT recalibration record.
                std::string encoded(as_culprit_value);
                std::size_t begin = 0;
                for (int selected = 0; selected < alt; ++selected) {
                    const auto end = encoded.find(',', begin);
                    if (selected == alt - 1) {
                        result.culprits[key] = encoded.substr(
                            begin, end == std::string::npos ? std::string::npos : end - begin);
                        break;
                    }
                    if (end == std::string::npos) break;
                    begin = end + 1;
                }
            } else if (as_culprit_length <= 0) {
                if (site_culprit != nullptr)
                    result.culprits[key] = site_culprit;
                else if (site_culprit_upper != nullptr)
                    result.culprits[key] = site_culprit_upper;
            }
            result.positive_labels[key] = positive;
            result.negative_labels[key] = negative;
        }
        result.allele_specific = result.allele_specific || require_allele_specific;
        free(as_culprit_value);
        free(site_culprit);
        free(site_culprit_upper);
        free(values);
        bcf_clear(record);
    }
    bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
    if (result.scores.empty() && result.site_scores.empty())
        throw std::runtime_error("BAD_INPUT: recal VCF contains no VQSLOD records");
    return result;
}

bool has_filters(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name && std::string(name) != "PASS" && std::string(name) != ".") return true;
    }
    return false;
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection ? "INTERSECTION" : "UNION";
}

void normalize_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.rid != rhs.rid) return lhs.rid < rhs.rid;
        if (lhs.begin != rhs.begin) return lhs.begin < rhs.begin;
        return lhs.end < rhs.end;
    });
    std::vector<fastgatk::io::IndexedInterval> merged;
    for (const auto& interval : intervals) {
        if (interval.end <= interval.begin) continue;
        if (!merged.empty() && merged.back().rid == interval.rid &&
            interval.begin <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }
    intervals.swap(merged);
}

void pad_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals, int padding) {
    if (padding <= 0) return;
    for (auto& interval : intervals) {
        interval.begin = std::max(0, interval.begin - padding);
        if (interval.end != std::numeric_limits<int>::max()) {
            const auto end = static_cast<std::int64_t>(interval.end) + padding;
            interval.end = static_cast<int>(std::min<std::int64_t>(
                std::numeric_limits<int>::max(), end));
        }
    }
    normalize_intervals(intervals);
}

void append_interval_selector_with_rule(
    const std::string& selector, const bcf_hdr_t* header,
    std::vector<fastgatk::io::IndexedInterval>& intervals,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule, bool first_selector) {
    std::vector<fastgatk::io::IndexedInterval> incoming;
    fastgatk::io::append_interval_selector(selector, header, incoming, stats);
    normalize_intervals(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        intervals.insert(intervals.end(), incoming.begin(), incoming.end());
        normalize_intervals(intervals);
        return;
    }
    std::vector<fastgatk::io::IndexedInterval> intersection;
    for (const auto& lhs : intervals) {
        for (const auto& rhs : incoming) {
            if (lhs.rid != rhs.rid) continue;
            const int begin = std::max(lhs.begin, rhs.begin);
            const int end = std::min(lhs.end, rhs.end);
            if (end > begin)
                intersection.push_back({lhs.rid, begin, end, lhs.contig});
        }
    }
    intervals.swap(intersection);
    normalize_intervals(intervals);
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    const auto reference_span = std::max<hts_pos_t>(1, record->rlen);
    int end = static_cast<int>(record->pos + reference_span);
    int32_t* values = nullptr;
    int count = 0;
    if (bcf_get_info_int32(header, record, "END", &values, &count) > 0 &&
        values != nullptr && count > 0 && values[0] > 0)
        end = values[0];
    free(values);
    return end;
}

bool record_in_intervals(const bcf_hdr_t* header, bcf1_t* record,
                         const std::vector<fastgatk::io::IndexedInterval>& intervals) {
    if (intervals.empty()) return false;
    const int begin = record->pos;
    const int end = record_end_exclusive(header, record);
    for (const auto& interval : intervals) {
        if (interval.rid != record->rid) continue;
        if (begin < interval.end && end > interval.begin) return true;
    }
    return false;
}

#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse(argc, argv);
#if !FASTGATK_HAS_HTSLIB
        (void)options;
        throw std::runtime_error("BACKEND_UNAVAILABLE: ApplyVQSR requires HTSlib");
#else
        const auto recal_scores = read_recal_vcf(options.recal_file, options.allele_specific);
        const auto& scores = recal_scores.scores;
        const auto& culprits = recal_scores.culprits;
        double cutoff = options.lod_cutoff;
        std::string filter_name = "LOW_VQSLOD";
        std::vector<Tranche> tranches;
        bool tranche_mode = false;
        if (!std::isfinite(cutoff) && !options.tranches_file.empty() &&
            std::isfinite(options.truth_sensitivity)) {
            const auto all_tranches = read_tranches(options.tranches_file, options.mode);
            if (all_tranches.empty()) throw std::runtime_error("BAD_INPUT: no tranches match mode " + options.mode);
            for (const auto& tranche : all_tranches)
                if (tranche.target + 1e-9 >= options.truth_sensitivity) tranches.push_back(tranche);
            if (tranches.empty())
                throw std::runtime_error("BAD_INPUT: no tranches at or above requested truth sensitivity");
            // ApplyVQSR's Java implementation reverses the selected rows:
            // [highest truth sensitivity ... lowest] for score walking.
            std::reverse(tranches.begin(), tranches.end());
            tranche_mode = true;
        } else if (!std::isfinite(cutoff)) {
            // ApplyVQSR does not enter tranche-walking mode merely because a
            // tranches file was supplied.  Without an explicit truth
            // sensitivity target GATK keeps VQSLOD >= 0 and emits LOW_VQSLOD;
            // this is also the default when no tranches file is supplied.
            cutoff = 0.0;
        }

        htsFile* input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open input VCF: " + options.input);
        bcf_hdr_t* input_header = bcf_hdr_read(input);
        if (input_header == nullptr) throw std::runtime_error("BAD_INPUT: cannot read input VCF header");
        std::vector<fastgatk::io::IndexedInterval> intervals;
        std::vector<fastgatk::io::IndexedInterval> excluded_intervals;
        fastgatk::io::IntervalFileStats interval_stats;
        bool first_interval_selector = true;
        for (const auto& selector : options.regions) {
            append_interval_selector_with_rule(
                selector, input_header, intervals, interval_stats,
                options.interval_set_rule, first_interval_selector);
            first_interval_selector = false;
        }
        normalize_intervals(intervals);
        fastgatk::io::IntervalFileStats exclusion_interval_stats;
        for (const auto& selector : options.excluded_regions) {
            // Exclusions are always unioned by GATK, independent of the
            // include interval set rule.  Applying them after include
            // normalization also keeps the traversal boundary explicit.
            fastgatk::io::append_interval_selector(selector, input_header,
                                                   excluded_intervals,
                                                   exclusion_interval_stats);
        }
        normalize_intervals(excluded_intervals);
        pad_intervals(excluded_intervals, options.exclusion_padding);
        // Keep the source header aware of every field/filter that can be
        // added below.  bcf_translate() remaps IDs using the source header;
        // registering these definitions there prevents newly-added training
        // flags and tranche filters from being discarded during translation.
        const auto ensure_input_header_line = [&](const std::string& line, const char* id,
                                                  int dictionary) {
            if (bcf_hdr_id2int(input_header, dictionary, id) < 0 && bcf_hdr_append(input_header, line.c_str()) != 0)
                throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot extend input header with ") + id);
        };
        ensure_input_header_line("##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=fastgatk VQSR score>", "VQSLOD", BCF_DT_ID);
        ensure_input_header_line("##INFO=<ID=culprit,Number=1,Type=String,Description=The annotation which was the worst performing in the Gaussian mixture model>", "culprit", BCF_DT_ID);
        ensure_input_header_line("##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the positive training set of good variants>", "POSITIVE_TRAIN_SITE", BCF_DT_ID);
        ensure_input_header_line("##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the negative training set of bad variants>", "NEGATIVE_TRAIN_SITE", BCF_DT_ID);
        if (options.allele_specific) {
            ensure_input_header_line("##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=fastgatk allele-specific VQSR score>", "AS_VQSLOD", BCF_DT_ID);
            ensure_input_header_line("##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=fastgatk allele-specific VQSR filter status>", "AS_FilterStatus", BCF_DT_ID);
            ensure_input_header_line("##INFO=<ID=AS_culprit,Number=A,Type=String,Description=fastgatk allele-specific VQSR culprit>", "AS_culprit", BCF_DT_ID);
        }
        if (tranche_mode) {
            for (const auto& tranche : tranches) {
                if (!tranche.filter_name.empty())
                    ensure_input_header_line("##FILTER=<ID=" + tranche.filter_name + ",Description=Truth sensitivity tranche>", tranche.filter_name.c_str(), BCF_DT_ID);
            }
            if (!tranches.empty()) {
                const auto plus_name = tranches.front().filter_name + "+";
                ensure_input_header_line("##FILTER=<ID=" + plus_name + ",Description=VQSLOD below selected tranche>", plus_name.c_str(), BCF_DT_ID);
            }
        } else {
            ensure_input_header_line("##FILTER=<ID=LOW_VQSLOD,Description=VQSR tranche below VQSLOD cutoff>", "LOW_VQSLOD", BCF_DT_ID);
        }
        if (bcf_hdr_sync(input_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync input header");
        bcf_hdr_t* output_header = input_header ? bcf_hdr_dup(input_header) : nullptr;
        bcf1_t* record = bcf_init();
        if (!input_header || !output_header || !record)
            throw std::runtime_error("BAD_INPUT: cannot read input VCF header");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "VQSLOD") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=fastgatk VQSR score>");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "culprit") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=culprit,Number=1,Type=String,Description=The annotation which was the worst performing in the Gaussian mixture model>");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "POSITIVE_TRAIN_SITE") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the positive training set of good variants>");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "NEGATIVE_TRAIN_SITE") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the negative training set of bad variants>");
        if (options.allele_specific) {
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AS_VQSLOD") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=fastgatk allele-specific VQSR score>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AS_FilterStatus") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=fastgatk allele-specific VQSR filter status>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AS_culprit") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AS_culprit,Number=A,Type=String,Description=fastgatk allele-specific VQSR culprit>");
        }
        if (tranche_mode) {
            for (const auto& tranche : tranches) {
                if (tranche.filter_name.empty()) continue;
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, tranche.filter_name.c_str()) < 0)
                    bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=Truth sensitivity tranche at VQSLOD %.4f>",
                                   tranche.filter_name.c_str(), tranche.min_vqslod);
            }
            if (!tranches.empty()) {
                const auto plus_name = tranches.front().filter_name + "+";
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, plus_name.c_str()) < 0)
                    bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=VQSLOD below the lowest selected tranche>", plus_name.c_str());
            }
        } else if (std::isfinite(cutoff) && bcf_hdr_id2int(output_header, BCF_DT_ID, filter_name.c_str()) < 0) {
            bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=VQSR tranche below VQSLOD cutoff>", filter_name.c_str());
        }
        bcf_hdr_append(output_header, "##source=fastgatk-apply-vqsr");
        bcf_hdr_append(output_header, "##fastgatk_apply_vqsr_status=prototype-recal-vcf");
        if (bcf_hdr_sync(output_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync VCF header");
        bcf_hdr_t* writer_header = bcf_hdr_dup(output_header);
        if (!writer_header)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate ApplyVQSR writer header");
        if (options.sites_only_vcf_output &&
            bcf_hdr_set_samples(writer_header, nullptr, 0) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot configure sites-only VCF header");
        htsFile* output = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
        if (!output || bcf_hdr_write(output, writer_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VCF header");
        const auto write_record = [&]() {
            if (bcf_translate(output_header, input_header, record) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VQSR record");
            // Filtering/scoring above consumes the complete input record.
            // Strip FORMAT/sample payload only after translation and all INFO/
            // FILTER updates, matching GATK's writer-only option.
            if (options.sites_only_vcf_output && bcf_subset_format(writer_header, record) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sites-only VCF record");
            if (bcf_write(output, writer_header, record) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VQSR record");
        };
        const auto filter_is_ignored = [&](bcf1_t* candidate) {
            if (options.ignore_all_filters) return true;
            if (options.ignore_filters.empty()) return false;
            bcf_unpack(candidate, BCF_UN_FLT);
            for (int index = 0; index < candidate->d.n_flt; ++index) {
                const auto* name = bcf_hdr_int2id(input_header, BCF_DT_ID, candidate->d.flt[index]);
                if (name != nullptr && std::find(options.ignore_filters.begin(), options.ignore_filters.end(), name) != options.ignore_filters.end())
                    return true;
            }
            return false;
        };
        std::uint64_t input_records = 0, scored_records = 0, filtered_records = 0, excluded_records = 0;
        std::uint64_t pre_filtered_records = 0, ignored_filtered_records = 0;
        std::uint64_t applicable_alleles = 0, mode_skipped_alleles = 0;
        std::uint64_t interval_skipped_records = 0;
        std::uint64_t interval_excluded_records = 0;
        while (fastgatk::io::read_variant_record(input, input_header, record,
                                               options.input) == 0) {
            if (!options.regions.empty() && !record_in_intervals(input_header, record, intervals)) {
                ++interval_skipped_records;
                bcf_clear(record);
                continue;
            }
            if (!excluded_intervals.empty() && record_in_intervals(input_header, record, excluded_intervals)) {
                ++interval_excluded_records;
                bcf_clear(record);
                continue;
            }
            ++input_records;
            bcf_unpack(record, BCF_UN_ALL);
            const auto type = variant_type(record);
            const bool apply_mode = record_matches_mode(type, options.mode);
            const bool record_has_filters = has_filters(input_header, record);
            const bool ignored_record_filters = record_has_filters &&
                (options.ignore_all_filters || filter_is_ignored(record));
            if (record_has_filters && !ignored_record_filters) {
                // GATK does not annotate a pre-filtered record unless the
                // caller explicitly asks to ignore that filter.
                ++pre_filtered_records;
                write_record();
                bcf_clear(record);
                continue;
            }
            if (ignored_record_filters) ++ignored_filtered_records;
            float site_score = std::numeric_limits<float>::quiet_NaN();
            bool score_found = false;
            double best_score = -std::numeric_limits<double>::infinity();
            std::string site_filter = "PASS";
            std::vector<float> allele_scores;
            std::vector<std::string> allele_statuses;
            std::vector<std::string> allele_culprits;
            bool applied_allele = false;
            bool positive_label = false;
            bool negative_label = false;
            if (apply_mode || options.allele_specific) {
                for (int alt = 1; alt < record->n_allele; ++alt) {
                    const bool applies_to_allele = !options.allele_specific ||
                        options.mode == "BOTH" || allele_matches_mode(record, alt, options.mode);
                    if (!applies_to_allele) {
                        if (options.allele_specific) {
                            allele_scores.push_back(std::numeric_limits<float>::quiet_NaN());
                            allele_statuses.emplace_back("NA");
                            allele_culprits.emplace_back("NA");
                            ++mode_skipped_alleles;
                        }
                        continue;
                    }
                    const auto key = variant_key(input_header, record, alt);
                    const auto location = variant_location_key(input_header, record);
                    auto iterator = scores.end();
                    auto site_iterator = recal_scores.site_scores.end();
                    if (options.allele_specific)
                        iterator = scores.find(key);
                    else
                        site_iterator = recal_scores.site_scores.find(location);
                    const bool found = options.allele_specific ? iterator != scores.end() : site_iterator != recal_scores.site_scores.end();
                    // Scalar recalibration tables use dummy N/<VQSR>
                    // alleles.  Accept a matching allele key as a useful
                    // compatibility fallback for native recal VCFs.
                    if (!found && !options.allele_specific) iterator = scores.find(key);
                    const bool found_after_fallback = options.allele_specific ? iterator != scores.end() :
                        (site_iterator != recal_scores.site_scores.end() || iterator != scores.end());
                    if (!found_after_fallback)
                        throw std::runtime_error("BACKEND_UNAVAILABLE: recal VCF lacks input variant key");
                    const float score = options.allele_specific ? iterator->second :
                        (site_iterator != recal_scores.site_scores.end() ? site_iterator->second : iterator->second);
                    if (options.allele_specific) allele_scores.push_back(score);
                    if (!score_found || score > site_score) site_score = score;
                    best_score = std::max(best_score, static_cast<double>(score));
                    score_found = true;
                    applied_allele = true;
                    if (options.allele_specific) {
                        ++applicable_alleles;
                        const auto culprit = culprits.find(key);
                        allele_culprits.push_back(culprit == culprits.end() ? "." : culprit->second);
                        const auto allele_filter = tranche_mode ? tranche_filter_for_score(tranches, score) :
                            (score < cutoff ? filter_name : "PASS");
                        allele_statuses.push_back(allele_filter);
                        const auto positive = recal_scores.positive_labels.find(key);
                        const auto negative = recal_scores.negative_labels.find(key);
                        positive_label = positive_label || (positive != recal_scores.positive_labels.end() && positive->second);
                        negative_label = negative_label || (negative != recal_scores.negative_labels.end() && negative->second);
                    }
                }
            }
            if (score_found) {
                ++scored_records;
                site_filter = tranche_mode ? tranche_filter_for_score(tranches, best_score) :
                    (best_score < cutoff ? filter_name : "PASS");
                if (options.allele_specific) {
                    if (bcf_update_info_float(output_header, record, "AS_VQSLOD",
                                              allele_scores.data(), static_cast<int>(allele_scores.size())) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_VQSLOD");
                    std::ostringstream status;
                    for (std::size_t index = 0; index < allele_statuses.size(); ++index) {
                        if (index) status << ',';
                        status << allele_statuses[index];
                    }
                    if (bcf_update_info_string(output_header, record, "AS_FilterStatus", status.str().c_str()) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_FilterStatus");
                    std::ostringstream culprit;
                    for (std::size_t index = 0; index < allele_culprits.size(); ++index) {
                        if (index) culprit << ',';
                        culprit << allele_culprits[index];
                    }
                    if (bcf_update_info_string(output_header, record, "AS_culprit", culprit.str().c_str()) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_culprit");
                } else {
                    if (bcf_update_info_float(output_header, record, "VQSLOD", &site_score, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VQSLOD");
                    const auto location = variant_location_key(input_header, record);
                    auto culprit = recal_scores.site_culprits.find(location);
                    const auto allele_culprit = culprits.find(variant_key(input_header, record, 1));
                    const std::string culprit_value = culprit != recal_scores.site_culprits.end() ? culprit->second :
                        (allele_culprit != culprits.end() ? allele_culprit->second : ".");
                    if (bcf_update_info_string(output_header, record, "culprit", culprit_value.c_str()) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write culprit");
                    const auto positive = recal_scores.site_positive_labels.find(location);
                    const auto negative = recal_scores.site_negative_labels.find(location);
                    positive_label = positive != recal_scores.site_positive_labels.end() && positive->second;
                    negative_label = negative != recal_scores.site_negative_labels.end() && negative->second;
                }
                if (positive_label) {
                    if (bcf_update_info_flag(output_header, record, "POSITIVE_TRAIN_SITE", nullptr, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write POSITIVE_TRAIN_SITE");
                }
                if (negative_label) {
                    if (bcf_update_info_flag(output_header, record, "NEGATIVE_TRAIN_SITE", nullptr, 1) != 0)
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write NEGATIVE_TRAIN_SITE");
                }
            }
            const bool filtered = score_found && site_filter != "PASS" && site_filter != "UNFILTERED";
            if (filtered) {
                ++filtered_records;
                if (options.exclude_filtered) { ++excluded_records; bcf_clear(record); continue; }
                if (bcf_update_filter(output_header, record, nullptr, 0) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear input FILTER");
                const int filter_id = bcf_hdr_id2int(output_header, BCF_DT_ID, site_filter.c_str());
                if (filter_id < 0 || bcf_add_filter(output_header, record, filter_id) < 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot add VQSR FILTER");
            } else if (score_found && site_filter == "PASS") {
                int pass_id = 0;  // HTSlib's reserved FILTER id for PASS.
                if (bcf_update_filter(output_header, record, &pass_id, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set PASS FILTER");
            }
            if (options.allele_specific && !applied_allele && !score_found) {
                // No allele in this mode was recalibrated; retain the input
                // record without manufacturing an empty AS annotation.
            }
            write_record();
            bcf_clear(record);
        }
        bcf_destroy(record); bcf_close(output); bcf_close(input);
        bcf_hdr_destroy(writer_header); bcf_hdr_destroy(output_header); bcf_hdr_destroy(input_header);
        std::string index_path;
        if (options.create_output_variant_index) {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VQSR Tabix index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        if (!complete_file(options.output) || (!index_path.empty() && !complete_file(index_path)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: VQSR output or index is incomplete");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest");
        manifest << "{\"schema_version\":1,\"tool\":\"ApplyVQSR\",\"implementation\":\"fastgatk-apply-vqsr\","
                 << "\"status\":\"prototype\",\"input_records\":" << input_records
                 << ",\"scored_records\":" << scored_records << ",\"filtered_records\":" << filtered_records
                 << ",\"excluded_records\":" << excluded_records << ",\"mode\":\"" << json_escape(options.mode)
                 << "\",\"allele_specific\":" << (options.allele_specific ? "true" : "false")
                 << ",\"recal_allele_specific\":" << (recal_scores.allele_specific ? "true" : "false")
                 << ",\"applicable_alleles\":" << applicable_alleles
                 << ",\"mode_skipped_alleles\":" << mode_skipped_alleles
                 << ",\"pre_filtered_records\":" << pre_filtered_records
                 << ",\"create_output_variant_index\":"
                 << (options.create_output_variant_index ? "true" : "false")
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"ignored_filtered_records\":" << ignored_filtered_records
                 << ",\"interval_skipped_records\":" << interval_skipped_records
                 << ",\"interval_excluded_records\":" << interval_excluded_records
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_stats.files
                 << ",\"interval_list_records\":" << interval_stats.records
                 << ",\"excluded_interval_list_inputs\":" << exclusion_interval_stats.files
                 << ",\"excluded_interval_list_records\":" << exclusion_interval_stats.records
                 << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                 << ",\"cutoff\":" << (std::isfinite(cutoff) ? std::to_string(cutoff) : "null")
                 << ",\"cutoff_source\":\""
                 << (tranche_mode ? "tranches" : (std::isfinite(options.lod_cutoff) ? "lod-score-cutoff" : "default-zero"))
                 << "\",\"tranche_count\":" << tranches.size()
                 << ",\"provenance\":{\"recal_file_size\":" << file_size_or_zero(options.recal_file)
                 << ",\"recal_file_mtime_ns\":" << file_mtime_ticks_or_zero(options.recal_file)
                 << ",\"tranches_file_size\":" << file_size_or_zero(options.tranches_file)
                 << ",\"tranches_file_mtime_ns\":" << file_mtime_ticks_or_zero(options.tranches_file)
                 << "}"
                 << ",\"recal_file\":\"" << json_escape(options.recal_file) << "\",\"output\":\""
                 << json_escape(options.output) << "\",\"bit_identical_to_gatk\":false}\n";
        manifest.close();
        if (!manifest || !complete_file(manifest_path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest");
        std::cout << "{\"tool\":\"ApplyVQSR\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"scored_records\":" << scored_records << ",\"filtered_records\":"
                  << filtered_records << ",\"excluded_records\":" << excluded_records
                  << ",\"allele_specific\":" << (options.allele_specific ? "true" : "false")
                  << ",\"create_output_variant_index\":"
                  << (options.create_output_variant_index ? "true" : "false")
                  << ",\"sites_only_vcf_output\":"
                  << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"interval_skipped_records\":" << interval_skipped_records
                  << ",\"interval_excluded_records\":" << interval_excluded_records
                  << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"applicable_alleles\":" << applicable_alleles
                  << ",\"mode_skipped_alleles\":" << mode_skipped_alleles << "}\n";
        return 0;
#endif
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
