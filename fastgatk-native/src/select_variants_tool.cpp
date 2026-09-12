#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/output.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
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
    std::vector<std::string> include_types;
    std::vector<std::string> exclude_types;
    std::vector<std::string> include_samples;
    std::vector<std::string> exclude_samples;
    std::vector<std::string> select_expressions;
    std::vector<std::string> select_genotype_expressions;
    std::string concordance;
    std::string discordance;
    std::string restrict_alleles;
    int min_indel_size = 0;
    int max_indel_size = std::numeric_limits<int>::max();
    bool exclude_filtered = false;
    bool exclude_non_variants = false;
    bool set_filtered_gt_to_nocall = false;
    bool keep_original_ac = false;
    bool keep_original_dp = false;
    bool remove_unused_alternates = false;
    bool compare_genotypes = false;
    bool create_index = true;
    bool sites_only_vcf_output = false;
};

struct GenotypeKernelTelemetry {
    std::uint64_t pl_remap_calls = 0;
    double pl_remap_prepare_seconds = 0.0;
    double pl_remap_execute_seconds = 0.0;
    std::uint64_t gt_gq_calls = 0;
    double gt_gq_prepare_seconds = 0.0;
    double gt_gq_execute_seconds = 0.0;
    std::uint64_t allele_field_remap_calls = 0;
    double allele_field_remap_prepare_seconds = 0.0;
    double allele_field_remap_execute_seconds = 0.0;
    std::string allele_field_remap_execution_space;
    std::string execution_space;
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
    return path.size() >= value.size() && path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-select-variants (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE              input VCF/BCF\n"
                         "  -O, --output FILE               output VCF/BCF\n"
                         "  --select-type-to-include TYPE   SNP/INDEL/MNP/MIXED (repeatable)\n"
                         "  --select-type-to-exclude TYPE   exclude variant type (repeatable)\n"
                         "  --restrict-alleles-to MODE      ALL/BI/MULTI\n"
                         "  --min-indel-size N              minimum absolute indel length\n"
                         "  --max-indel-size N              maximum absolute indel length\n"
                         "  --sample-name NAME              keep sample (repeatable)\n"
                         "  --exclude-sample-name NAME      remove sample (repeatable)\n"
                         "  --select-expression EXPR       JEXL subset (repeatable)\n"
                         "  --select-genotype EXPR        genotype JEXL subset over any sample (repeatable)\n"
                         "  --concordance FILE             retain records with matching REF/ALT sites\n"
                         "  --discordance FILE             retain records absent from REF/ALT sites\n"
                         "  --concordance-genotypes        compare GTs for shared samples\n"
                         "  --interval-set-rule RULE       UNION (default) or INTERSECTION\n"
                         "  --exclude-filtered              drop non-PASS records\n"
                         "  --exclude-non-variants[=BOOL]    drop sites without a called ALT genotype (default false)\n"
                         "  --set-filtered-gt-to-nocall[=BOOL]  set filtered called genotypes to no-call\n"
                         "  --sites-only-vcf-output[=BOOL]  omit FORMAT/sample columns at final writer\n"
                         "  --keep-original-ac              preserve AC/AF/AN before subsetting\n"
                         "  --keep-original-dp              preserve DP before subsetting\n"
                         "  --remove-unused-alternates      subset ALT/GT/AD/PL\n"
                         "      --output-manifest FILE     OutputManifest JSON\n";
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
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (is_option(argument, "--select-type-to-include"))
            options.include_types.push_back(require_value(index, argc, argv, argument, "--select-type-to-include"));
        else if (is_option(argument, "--select-type-to-exclude"))
            options.exclude_types.push_back(require_value(index, argc, argv, argument, "--select-type-to-exclude"));
        else if (is_option(argument, "--restrict-alleles-to"))
            options.restrict_alleles = require_value(index, argc, argv, argument, "--restrict-alleles-to");
        else if (is_option(argument, "--min-indel-size")) {
            const auto value = require_value(index, argc, argv, argument, "--min-indel-size");
            try { options.min_indel_size = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --min-indel-size: " + value); }
            if (options.min_indel_size < 0)
                throw std::invalid_argument("--min-indel-size must be non-negative");
        } else if (is_option(argument, "--max-indel-size")) {
            const auto value = require_value(index, argc, argv, argument, "--max-indel-size");
            try { options.max_indel_size = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --max-indel-size: " + value); }
            if (options.max_indel_size < 0)
                throw std::invalid_argument("--max-indel-size must be non-negative");
        }
        else if (is_option(argument, "--sample-name") || is_option(argument, "--sample-name-to-include"))
            options.include_samples.push_back(require_value(index, argc, argv, argument,
                argument.rfind("--sample-name-to-include", 0) == 0 ? "--sample-name-to-include" : "--sample-name"));
        else if (is_option(argument, "--exclude-sample-name") || is_option(argument, "--sample-name-to-exclude"))
            options.exclude_samples.push_back(require_value(index, argc, argv, argument,
                argument.rfind("--sample-name-to-exclude", 0) == 0 ? "--sample-name-to-exclude" : "--exclude-sample-name"));
        else if (is_option(argument, "--select-expression") || is_option(argument, "--select"))
            options.select_expressions.push_back(require_value(index, argc, argv, argument,
                (argument == "--select" || !option_value(argument, "--select").empty())
                    ? "--select" : "--select-expression"));
        else if (is_option(argument, "--select-genotype") || is_option(argument, "--select-genotype-expressions"))
            options.select_genotype_expressions.push_back(require_value(
                index, argc, argv, argument,
                argument.rfind("--select-genotype-expressions", 0) == 0
                    ? "--select-genotype-expressions" : "--select-genotype"));
        else if (is_option(argument, "--concordance"))
            options.concordance = require_value(index, argc, argv, argument, "--concordance");
        else if (is_option(argument, "--discordance"))
            options.discordance = require_value(index, argc, argv, argument, "--discordance");
        else if (argument == "--concordance-genotypes" || argument == "--discordance-genotypes" ||
                 argument == "--sample-level-concordance")
            options.compare_genotypes = true;
        else if (argument == "--exclude-filtered" ||
                 argument.rfind("--exclude-filtered=", 0) == 0)
            options.exclude_filtered = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--exclude-filtered");
        else if (argument == "--exclude-filtered-variants" ||
                 argument.rfind("--exclude-filtered-variants=", 0) == 0)
            options.exclude_filtered = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--exclude-filtered-variants");
        else if (argument == "--exclude-non-variants" ||
                 argument.rfind("--exclude-non-variants=", 0) == 0)
            options.exclude_non_variants = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--exclude-non-variants");
        else if (argument == "--exclude-non-variant-sites" ||
                 argument.rfind("--exclude-non-variant-sites=", 0) == 0)
            options.exclude_non_variants = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--exclude-non-variant-sites");
        else if (argument == "--set-filtered-gt-to-nocall" ||
                 argument.rfind("--set-filtered-gt-to-nocall=", 0) == 0) {
            options.set_filtered_gt_to_nocall = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--set-filtered-gt-to-nocall");
        }
        else if (argument == "--keep-original-ac" ||
                 argument.rfind("--keep-original-ac=", 0) == 0) {
            options.keep_original_ac = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--keep-original-ac");
        } else if (argument == "--keep-original-dp" ||
                   argument.rfind("--keep-original-dp=", 0) == 0) {
            options.keep_original_dp = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--keep-original-dp");
        }
        else if (argument == "--remove-unused-alternates") options.remove_unused_alternates = true;
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--sites-only-vcf-output") {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--sites-only-vcf-output");
        } else if (argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            // Barclay rejects an embedded '=' spelling for this optional
            // Boolean (the value is accepted only as a following token).
            throw std::invalid_argument("unknown option: " + argument);
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // accepted compatibility flags
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
    if (!options.restrict_alleles.empty() && options.restrict_alleles != "ALL" &&
        options.restrict_alleles != "BI" && options.restrict_alleles != "MULTI")
        throw std::invalid_argument("--restrict-alleles-to must be ALL, BI or MULTI");
    if (!options.concordance.empty() && !options.discordance.empty())
        throw std::invalid_argument("--concordance and --discordance are mutually exclusive");
    if (options.min_indel_size > options.max_indel_size)
        throw std::invalid_argument("--min-indel-size cannot exceed --max-indel-size");
    if (options.compare_genotypes && options.concordance.empty() && options.discordance.empty())
        throw std::invalid_argument("--concordance-genotypes requires --concordance or --discordance");
    return options;
}

#if FASTGATK_HAS_HTSLIB

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
    return value;
}

bool select_allele_is_breakpoint(const std::string& allele) {
    return allele.size() > 1 && allele.find_first_of("[]") != std::string::npos;
}

bool select_allele_is_single_breakend(const std::string& allele) {
    return allele.size() > 1 && (allele.front() == '.' || allele.back() == '.');
}

bool select_allele_is_symbolic(const std::string& allele) {
    return allele == "*" || (allele.size() > 1 &&
        (allele.front() == '<' || allele.back() == '>' ||
         select_allele_is_breakpoint(allele) || select_allele_is_single_breakend(allele)));
}

std::string variant_type(bcf1_t* record) {
    const auto types = bcf_get_variant_types(record);
    const bool snp = (types & VCF_SNP) != 0;
    const bool indel = (types & VCF_INDEL) != 0;
    const bool mnp = (types & VCF_MNP) != 0;
    if (types == VCF_REF) return "NO_VARIATION";
    if (snp && indel) return "MIXED";
    if (snp) return "SNP";
    if (indel) return "INDEL";
    if (mnp) return "MNP";
    for (int allele = 1; allele < record->n_allele; ++allele)
        if (record->d.allele[allele] != nullptr &&
            select_allele_is_symbolic(record->d.allele[allele])) return "SYMBOLIC";
    return "OTHER";
}

bool concrete_or_symbolic_alt_present(bcf1_t* record) {
    if (record->n_allele <= 1 || record->d.allele == nullptr) return false;
    bcf_unpack(record, BCF_UN_STR);
    for (int index = 1; index < record->n_allele; ++index) {
        const auto* allele = record->d.allele[index];
        if (allele != nullptr && allele[0] != '\0' && std::string(allele) != ".") return true;
    }
    return false;
}

// SelectVariants --exclude-non-variants uses VariantContext.isVariant().
// With samples, a site with only hom-ref/no-call genotypes is non-variant even
// when an ALT is listed in the VCF.  A sample-free record is judged from its
// ALT column itself, which keeps annotation-only VCFs selectable.
bool is_non_variant(const bcf_hdr_t* header, bcf1_t* record) {
    if (!concrete_or_symbolic_alt_present(record)) return true;
    if (record->n_sample <= 0) return false;
    int32_t* genotypes = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
    if (length <= 0 || genotypes == nullptr || count <= 0) {
        free(genotypes);
        return false;
    }
    bool called_alt = false;
    for (int index = 0; index < count; ++index) {
        const auto encoded = genotypes[index];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
        const int allele = bcf_gt_allele(encoded);
        if (allele > 0 && allele < record->n_allele) {
            called_alt = true;
            break;
        }
    }
    free(genotypes);
    return !called_alt;
}

bool indel_size_range(bcf1_t* record, int& minimum, int& maximum) {
    minimum = std::numeric_limits<int>::max();
    maximum = 0;
    if (record->n_allele <= 1 || record->d.allele == nullptr) return false;
    bcf_unpack(record, BCF_UN_STR);
    const auto* reference = record->d.allele[0];
    if (reference == nullptr || reference[0] == '\0') return false;
    const std::size_t reference_length = std::strlen(reference);
    bool found = false;
    for (int index = 1; index < record->n_allele; ++index) {
        const auto* alternate = record->d.allele[index];
        if (alternate == nullptr || alternate[0] == '\0' || alternate[0] == '<' ||
            std::string(alternate) == "*" || std::string(alternate).find_first_of("[]") != std::string::npos)
            continue;
        const auto alternate_length = std::strlen(alternate);
        const auto difference = reference_length > alternate_length
            ? reference_length - alternate_length : alternate_length - reference_length;
        if (difference == 0) continue;
        const auto size = difference > static_cast<std::size_t>(std::numeric_limits<int>::max())
            ? std::numeric_limits<int>::max() : static_cast<int>(difference);
        minimum = std::min(minimum, size);
        maximum = std::max(maximum, size);
        found = true;
    }
    return found;
}

std::string variant_key(const bcf_hdr_t* header, bcf1_t* record) {
    const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    if (contig == nullptr || record->n_allele < 2 || record->d.allele == nullptr ||
        record->d.allele[0] == nullptr)
        return {};
    std::ostringstream key;
    key << contig << ':' << record->pos << ':' << record->d.allele[0] << ':';
    for (int allele = 1; allele < record->n_allele; ++allele) {
        if (allele != 1) key << ',';
        key << (record->d.allele[allele] == nullptr ? "" : record->d.allele[allele]);
    }
    return key.str();
}

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

void normalize_select_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.end <= region.begin) continue;
        if (!merged.empty() && merged.back().rid == region.rid && region.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, region.end);
        else merged.push_back(region);
    }
    regions.swap(merged);
}

void append_select_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_select_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < regions[left].rid) { ++right; continue; }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end) intersection.push_back(Region{regions[left].rid, begin, end});
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
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        return record->rid == region.rid && record_end > region.begin && record->pos < region.end;
    });
}

std::string genotype_string(const bcf_hdr_t* header, bcf1_t* record, int sample) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &values, &count);
    if (length <= 0 || record->n_sample <= 0 || count < record->n_sample ||
        count % record->n_sample != 0 || sample < 0 || sample >= record->n_sample) {
        free(values);
        return {};
    }
    const int width = count / record->n_sample;
    std::vector<int> alleles;
    alleles.reserve(static_cast<std::size_t>(width));
    for (int index = 0; index < width; ++index) {
        const auto encoded = values[sample * width + index];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
            free(values);
            return {};
        }
        alleles.push_back(bcf_gt_allele(encoded));
    }
    free(values);
    std::sort(alleles.begin(), alleles.end());
    std::ostringstream result;
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        if (index != 0) result << '/';
        result << alleles[index];
    }
    return result.str();
}

struct ComparisonData {
    std::set<std::string> site_keys;
    std::set<std::string> genotype_keys;
    std::set<std::string> sample_names;
};

ComparisonData load_comparison_data(const std::string& path) {
    ComparisonData data;
    htsFile* input = bcf_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open comparison VCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (!header) {
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read comparison VCF header: " + path);
    }
    bcf1_t* record = bcf_init();
    if (!record) {
        bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate comparison VCF record");
    }
    for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
        if (name != nullptr) data.sample_names.emplace(name);
    }
    while (fastgatk::io::read_variant_record(input, header, record, path) == 0) {
        bcf_unpack(record, BCF_UN_STR);
        const auto key = variant_key(header, record);
        if (!key.empty()) {
            data.site_keys.insert(key);
            for (int sample = 0; sample < record->n_sample; ++sample) {
                const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                const auto genotype = genotype_string(header, record, sample);
                if (name != nullptr && !genotype.empty())
                    data.genotype_keys.insert(key + "\t" + name + "\t" + genotype);
            }
        }
        bcf_clear(record);
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(input);
    return data;
}

bool contains(const std::vector<std::string>& values, const std::string& needle) {
    const auto wanted = upper(needle);
    return std::any_of(values.begin(), values.end(), [&](const std::string& value) {
        return upper(value) == wanted;
    });
}

bool filtered(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->d.n_flt == 0) return false;
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && std::strcmp(name, "PASS") != 0 && std::strcmp(name, ".") != 0) return true;
    }
    return false;
}

void canonicalize_site_filters(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->d.n_flt <= 1) return;
    std::vector<int> filter_ids;
    filter_ids.reserve(static_cast<std::size_t>(record->d.n_flt));
    for (int index = 0; index < record->d.n_flt; ++index)
        filter_ids.push_back(record->d.flt[index]);
    std::stable_sort(filter_ids.begin(), filter_ids.end(), [&](int left, int right) {
        const auto* left_name = bcf_hdr_int2id(header, BCF_DT_ID, left);
        const auto* right_name = bcf_hdr_int2id(header, BCF_DT_ID, right);
        if (left_name == nullptr || right_name == nullptr)
            return left < right;
        return std::string(left_name) < std::string(right_name);
    });
    if (bcf_update_filter(header, record, filter_ids.data(),
                          static_cast<int>(filter_ids.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot canonicalize SelectVariants FILTER labels");
}

std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool wrapped_by_parentheses(const std::string& value) {
    if (value.size() < 2 || value.front() != '(' || value.back() != ')') return false;
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const auto character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') ++depth;
        else if (character == ')' && --depth == 0 && index + 1 != value.size()) return false;
    }
    return depth == 0;
}

std::size_t top_level_token(const std::string& value, const std::string& token) {
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    for (std::size_t index = 0; index + token.size() <= value.size(); ++index) {
        const auto character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') { ++depth; continue; }
        if (character == ')') { --depth; continue; }
        if (depth == 0 && value.compare(index, token.size(), token) == 0) return index;
    }
    return std::string::npos;
}

struct SelectNumericExpr {
    enum class Kind { Literal, Field, Add, Subtract, Multiply, Divide, Modulo, Negate };
    Kind kind = Kind::Field;
    double literal = 0.0;
    std::string field;
    std::shared_ptr<SelectNumericExpr> left;
    std::shared_ptr<SelectNumericExpr> right;
};

std::pair<std::size_t, char> select_top_level_arithmetic(const std::string& value) {
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    std::size_t additive = std::string::npos;
    std::size_t multiplicative = std::string::npos;
    char additive_op = 0;
    char multiplicative_op = 0;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const char character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') { ++depth; continue; }
        if (character == ')') { --depth; continue; }
        if (depth != 0 || (character != '+' && character != '-' &&
                           character != '*' && character != '/' && character != '%'))
            continue;
        if ((character == '+' || character == '-') &&
            (index == 0 || value[index - 1] == '(' || value[index - 1] == '+' ||
             value[index - 1] == '-' || value[index - 1] == '*' ||
             value[index - 1] == '/' || value[index - 1] == '%' ||
             ((value[index - 1] == 'e' || value[index - 1] == 'E') && index > 1 &&
              std::isdigit(static_cast<unsigned char>(value[index - 2])))))
            continue;
        if (character == '+' || character == '-') {
            additive = index;
            additive_op = character;
        } else {
            multiplicative = index;
            multiplicative_op = character;
        }
    }
    if (additive != std::string::npos) return {additive, additive_op};
    if (multiplicative != std::string::npos) return {multiplicative, multiplicative_op};
    return {std::string::npos, 0};
}

std::shared_ptr<SelectNumericExpr> parse_select_numeric_expression(const std::string& source) {
    auto value = trim_copy(source);
    if (value.empty()) throw std::invalid_argument("empty numeric SelectVariants expression");
    while (wrapped_by_parentheses(value)) value = trim_copy(value.substr(1, value.size() - 2));
    const auto split = select_top_level_arithmetic(value);
    if (split.first != std::string::npos) {
        auto node = std::make_shared<SelectNumericExpr>();
        if (split.second == '+') node->kind = SelectNumericExpr::Kind::Add;
        else if (split.second == '-') node->kind = SelectNumericExpr::Kind::Subtract;
        else if (split.second == '*') node->kind = SelectNumericExpr::Kind::Multiply;
        else if (split.second == '/') node->kind = SelectNumericExpr::Kind::Divide;
        else node->kind = SelectNumericExpr::Kind::Modulo;
        node->left = parse_select_numeric_expression(value.substr(0, split.first));
        node->right = parse_select_numeric_expression(value.substr(split.first + 1));
        return node;
    }
    if (!value.empty() && value.front() == '-') {
        auto node = std::make_shared<SelectNumericExpr>();
        node->kind = SelectNumericExpr::Kind::Negate;
        node->left = parse_select_numeric_expression(value.substr(1));
        return node;
    }
    std::size_t consumed = 0;
    try {
        const auto literal = std::stod(value, &consumed);
        if (consumed == value.size()) {
            auto node = std::make_shared<SelectNumericExpr>();
            node->kind = SelectNumericExpr::Kind::Literal;
            node->literal = literal;
            return node;
        }
    } catch (...) {
    }
    auto node = std::make_shared<SelectNumericExpr>();
    node->kind = SelectNumericExpr::Kind::Field;
    node->field = value;
    return node;
}

bool select_numeric_expression_has_operator(const SelectNumericExpr& expression) {
    return expression.kind != SelectNumericExpr::Kind::Literal &&
           expression.kind != SelectNumericExpr::Kind::Field;
}

std::string unquote(std::string value) {
    value = trim_copy(std::move(value));
    if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                              (value.front() == '\'' && value.back() == '\'')))
        return value.substr(1, value.size() - 2);
    return value;
}

bool parse_call(const std::string& raw, std::string& name, std::string& argument) {
    const auto value = trim_copy(raw);
    const auto open = value.find('(');
    if (open == std::string::npos || value.empty() || value.back() != ')' || open == 0) return false;
    name = trim_copy(value.substr(0, open));
    argument = unquote(value.substr(open + 1, value.size() - open - 2));
    return true;
}

bool select_compare(double value, const std::string& op, double threshold) {
    if (op == "<") return value < threshold;
    if (op == "<=") return value <= threshold;
    if (op == ">") return value > threshold;
    if (op == ">=") return value >= threshold;
    if (op == "==") return value == threshold;
    if (op == "!=") return value != threshold;
    return false;
}

// Accept both the compact INFO-vector spelling (AF[0]) and the explicit
// VariantContext form (vc.getAttribute("AF")[0]).  Returning the base field
// lets the existing HTSlib readers preserve the declared numeric type while
// selecting the requested zero-based element.
bool parse_select_indexed_info_field(const std::string& raw, std::string& base, int& index) {
    const auto value = trim_copy(raw);
    if (value.size() < 4) return false;
    std::string index_text;
    std::size_t base_end = std::string::npos;
    if (value.back() == ']') {
        const auto open = value.rfind('[');
        if (open == std::string::npos || open == 0 || open + 1 >= value.size() - 1) return false;
        index_text = value.substr(open + 1, value.size() - open - 2);
        base_end = open;
    } else if (value.back() == ')') {
        // HTSJDK's JEXL collection accessor is
        // vc.getAttribute("TAG").get(i).
        const auto token = value.rfind(".get(");
        if (token == std::string::npos || token == 0 || token + 5 >= value.size() - 1) return false;
        index_text = value.substr(token + 5, value.size() - token - 6);
        base_end = token;
    } else {
        // HTSJDK/JEXL also supports a dot-indexed collection element such as
        // `vc.getAttribute("AF").1`.  Restrict the suffix to decimal digits
        // so a normal dotted method/field is not misclassified.
        const auto dot = value.rfind('.');
        if (dot == std::string::npos || dot == 0 || dot + 1 >= value.size()) return false;
        index_text = value.substr(dot + 1);
        if (!std::all_of(index_text.begin(), index_text.end(),
                         [](unsigned char character) { return std::isdigit(character) != 0; }))
            return false;
        base_end = dot;
    }
    try {
        std::size_t consumed = 0;
        const auto parsed = std::stoi(index_text, &consumed);
        if (consumed != index_text.size() || parsed < 0) return false;
        base = trim_copy(value.substr(0, base_end));
        if (base.empty()) return false;
        index = parsed;
        return true;
    } catch (...) {
        return false;
    }
}

bool select_has_attribute(const bcf_hdr_t* header, bcf1_t* record, const std::string& tag) {
    if (tag.empty() || bcf_hdr_id2int(header, BCF_DT_ID, tag.c_str()) < 0) return false;
    bcf_unpack(record, BCF_UN_INFO);
    const auto id = bcf_hdr_id2int(header, BCF_DT_ID, tag.c_str());
    for (int index = 0; index < record->n_info; ++index)
        if (record->d.info[index].key == id && record->d.info[index].vptr != nullptr) return true;
    return false;
}

bool select_filter_set_contains(const bcf_hdr_t* header, bcf1_t* record, const std::string& wanted) {
    if (wanted.empty()) return false;
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && wanted == name) return true;
    }
    return false;
}

bool select_transition_state(bcf1_t* record, bool want_transition) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2 || record->d.allele == nullptr || record->d.allele[0] == nullptr ||
        std::strlen(record->d.allele[0]) != 1)
        return false;
    const char reference = static_cast<char>(std::toupper(record->d.allele[0][0]));
    bool seen = false;
    for (int index = 1; index < record->n_allele; ++index) {
        const char* alternate = record->d.allele[index];
        if (alternate == nullptr || std::strlen(alternate) != 1) return false;
        const char allele = static_cast<char>(std::toupper(alternate[0]));
        const bool transition = (reference == 'A' && allele == 'G') ||
                                (reference == 'G' && allele == 'A') ||
                                (reference == 'C' && allele == 'T') ||
                                (reference == 'T' && allele == 'C');
        if (transition != want_transition) return false;
        seen = true;
    }
    return seen;
}

// VariantContext's sample-aware polymorphism predicates are defined over
// called genotype alleles, rather than over the presence of ALT strings in
// the record.  Keep this lowering explicit so 0/0-only records remain
// selectable as monomorphic even when the VCF carries unused alternates.
bool select_polymorphic_in_samples(const bcf_hdr_t* header, bcf1_t* record) {
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || record->n_sample <= 0 || genotype_count <= 0 ||
        genotype_count % record->n_sample != 0) {
        free(genotypes);
        return false;
    }
    const int ploidy = genotype_count / record->n_sample;
    bool polymorphic = false;
    for (int sample = 0; sample < record->n_sample && !polymorphic; ++sample) {
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto encoded = genotypes[sample * ploidy + allele];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            if (bcf_gt_allele(encoded) > 0) {
                polymorphic = true;
                break;
            }
        }
    }
    free(genotypes);
    return polymorphic;
}

struct SelectGenotypeCounts {
    int called_chromosomes = 0;
    int no_call_count = 0;
    int hom_ref_count = 0;
    int het_count = 0;
    int hom_var_count = 0;
    bool has_genotypes = false;
};

bool select_genotype_counts(const bcf_hdr_t* header, bcf1_t* record,
                            SelectGenotypeCounts& counts) {
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || record->n_sample <= 0 || genotype_count <= 0 ||
        genotype_count % record->n_sample != 0) {
        free(genotypes);
        return false;
    }
    counts.has_genotypes = true;
    const int ploidy = genotype_count / record->n_sample;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        bool called = true;
        bool any_alt = false;
        bool all_same = true;
        int first_allele = -1;
        for (int copy = 0; copy < ploidy; ++copy) {
            const auto encoded = genotypes[sample * ploidy + copy];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                called = false;
                continue;
            }
            ++counts.called_chromosomes;
            const int allele = bcf_gt_allele(encoded);
            if (first_allele < 0) first_allele = allele;
            else if (allele != first_allele) all_same = false;
            if (allele > 0) any_alt = true;
        }
        if (!called) {
            ++counts.no_call_count;
        } else if (all_same && first_allele == 0) {
            ++counts.hom_ref_count;
        } else if (all_same && any_alt) {
            ++counts.hom_var_count;
        } else {
            ++counts.het_count;
        }
    }
    free(genotypes);
    return true;
}

bool select_site_method(const std::string& raw, bcf1_t* record, const bcf_hdr_t* header, bool& value) {
    if (raw == "vc.getFilters().isEmpty()" || raw == "getFilters().isEmpty()") {
        bcf_unpack(record, BCF_UN_FLT);
        value = true;
        for (int index = 0; index < record->d.n_flt; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
            if (name != nullptr && std::strcmp(name, "PASS") != 0 && std::strcmp(name, ".") != 0) {
                value = false;
                break;
            }
        }
        return true;
    }
    if (raw == "vc.getGenotypes().isEmpty()" || raw == "getGenotypes().isEmpty()") {
        value = record->n_sample == 0;
        return true;
    }
    std::string name;
    std::string argument;
    if (!parse_call(raw, name, argument)) return false;
    if (name.rfind("vc.", 0) == 0) name.erase(0, 3);
    const auto type = variant_type(record);
    if (name == "isSNP" && argument.empty()) value = type == "SNP";
    else if (name == "isIndel" && argument.empty()) value = type == "INDEL";
    else if (name == "isMNP" && argument.empty()) value = type == "MNP";
    else if (name == "isTransition" && argument.empty()) value = select_transition_state(record, true);
    else if (name == "isTransversion" && argument.empty()) value = select_transition_state(record, false);
    else if (name == "isMixed" && argument.empty()) {
        std::set<std::string> categories;
        if (record->n_allele > 1 && record->d.allele != nullptr && record->d.allele[0] != nullptr) {
            const auto reference_length = std::strlen(record->d.allele[0]);
            for (int index = 1; index < record->n_allele; ++index) {
                const char* alternate = record->d.allele[index];
                if (alternate == nullptr) continue;
                if (alternate[0] == '<' || alternate[0] == '*') categories.insert("SYMBOLIC");
                else if (std::strlen(alternate) == reference_length && reference_length > 1) categories.insert("MNP");
                else if (std::strlen(alternate) == reference_length) categories.insert("SNP");
                else categories.insert("INDEL");
            }
        }
        value = categories.size() > 1;
    } else if (name == "isSymbolic" && argument.empty()) value = type == "SYMBOLIC";
    else if (name == "isVariant" && argument.empty()) value = record->n_allele > 1 && type != "NO_VARIATION";
    else if (name == "isBiallelic" && argument.empty()) value = record->n_allele == 2;
    else if (name == "isMultiallelic" && argument.empty()) value = record->n_allele > 2;
    else if (name == "isFiltered" && argument.empty()) value = filtered(header, record);
    else if ((name == "isPass" || name == "isNotFiltered") && argument.empty()) value = !filtered(header, record);
    else if (name == "isNoVariation" && argument.empty()) value = type == "NO_VARIATION";
    else if (name == "isPolymorphicInSamples" && argument.empty())
        value = select_polymorphic_in_samples(header, record);
    else if (name == "isMonomorphicInSamples" && argument.empty())
        value = !select_polymorphic_in_samples(header, record);
    else if (name == "hasGenotypes" && argument.empty()) {
        SelectGenotypeCounts counts;
        value = select_genotype_counts(header, record, counts) && counts.has_genotypes;
    }
    else if (name == "hasAlternateAllele" && !argument.empty()) {
        try {
            const auto index = std::stoi(argument);
            value = index >= 0 && index + 1 < record->n_allele;
        } catch (...) {
            value = false;
        }
    }
    else if (name == "hasAttribute" && !argument.empty()) value = select_has_attribute(header, record, argument);
    else return false;
    return true;
}

bool select_read_string_field(const bcf_hdr_t* header, bcf1_t* record,
                              const std::string& raw_field, std::string& value) {
    auto field = trim_copy(raw_field);
    if (field == "vc.getType()" || field == "getType()") {
        value = variant_type(record);
        return true;
    }
    if (field == "vc.getID()" || field == "getID()") {
        bcf_unpack(record, BCF_UN_STR);
        if (record->d.id == nullptr || std::strcmp(record->d.id, ".") == 0) return false;
        value = record->d.id;
        return true;
    }
    if (field == "vc.getContig()" || field == "getContig()") {
        const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
        if (contig == nullptr) return false;
        value = contig;
        return true;
    }
    if (field == "vc.getFilters()" || field == "getFilters()") {
        bcf_unpack(record, BCF_UN_FLT);
        value.clear();
        for (int index = 0; index < record->d.n_flt; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
            if (name == nullptr || std::strcmp(name, "PASS") == 0 || std::strcmp(name, ".") == 0) continue;
            if (!value.empty()) value.push_back(';');
            value += name;
        }
        return true;
    }
    const std::string reference_prefix = "vc.getReference().getBaseString()";
    if (field == reference_prefix || field == "getReference().getBaseString()") {
        bcf_unpack(record, BCF_UN_STR);
        if (record->d.allele == nullptr || record->n_allele == 0 || record->d.allele[0] == nullptr)
            return false;
        value = record->d.allele[0];
        return true;
    }
    const std::string alternate_prefix = "vc.getAlternateAllele(";
    const std::string alternate_suffix = ").getBaseString()";
    if (field.rfind(alternate_prefix, 0) == 0 && field.size() > alternate_prefix.size() + alternate_suffix.size() &&
        field.compare(field.size() - alternate_suffix.size(), alternate_suffix.size(), alternate_suffix) == 0) {
        try {
            const auto index_text = field.substr(alternate_prefix.size(),
                field.size() - alternate_prefix.size() - alternate_suffix.size());
            const auto index = std::stoi(index_text);
            bcf_unpack(record, BCF_UN_STR);
            if (index < 0 || index + 1 >= record->n_allele || record->d.allele[index + 1] == nullptr)
                return false;
            value = record->d.allele[index + 1];
            return true;
        } catch (...) {
            return false;
        }
    }
    const std::string prefix = "vc.getAttribute(";
    if (field.rfind(prefix, 0) == 0 && !field.empty() && field.back() == ')')
        field = unquote(field.substr(prefix.size(), field.size() - prefix.size() - 1));
    if (field.empty()) return false;
    char* text = nullptr;
    int text_size = 0;
    const auto length = bcf_get_info_string(header, record, field.c_str(), &text, &text_size);
    if (length >= 0 && text != nullptr) {
        value.assign(text, static_cast<std::size_t>(length));
        free(text);
        return true;
    }
    free(text);
    // JEXL frequently calls String methods on numeric INFO attributes.  Keep
    // the conversion deterministic and locale-independent instead of rejecting
    // expressions such as vc.getAttribute("DP").contains("2").
    int32_t* integers = nullptr;
    int integer_count = 0;
    const auto integer_length = bcf_get_info_int32(header, record, field.c_str(), &integers, &integer_count);
    if (integer_length > 0 && integer_count > 0 && integers[0] != bcf_int32_missing &&
        integers[0] != bcf_int32_vector_end) {
        value = std::to_string(integers[0]);
        free(integers);
        return true;
    }
    free(integers);
    return false;
}

bool parse_select_allele_method(const std::string& raw, bool& reference_allele,
                                int& allele_index, std::string& method) {
    const auto value = trim_copy(raw);
    const std::string reference_prefixes[] = {"vc.getReference().", "getReference()."};
    for (const auto& prefix : reference_prefixes) {
        if (value.rfind(prefix, 0) != 0) continue;
        const auto suffix = value.substr(prefix.size());
        if (suffix == "isSymbolic()" || suffix == "isReference()" ||
            suffix == "isNoCall()" || suffix == "isCalled()" ||
            suffix == "isBreakpoint()" || suffix == "isSingleBreakend()" ||
            suffix == "isNonReference()" || suffix == "isNonRefAllele()" ||
            suffix == "length()") {
            reference_allele = true;
            allele_index = -1;
            method = suffix;
            return true;
        }
    }
    const std::string prefixes[] = {"vc.getAlternateAllele(", "getAlternateAllele("};
    const std::string* selected = nullptr;
    for (const auto& prefix : prefixes) {
        if (value.rfind(prefix, 0) == 0) {
            selected = &prefix;
            break;
        }
    }
    if (selected == nullptr) return false;
    const auto close = value.find(").", selected->size());
    if (close == std::string::npos || close + 2 >= value.size()) return false;
    try {
        std::size_t consumed = 0;
        const auto index_text = value.substr(selected->size(), close - selected->size());
        const auto index = std::stoi(index_text, &consumed);
        if (consumed != index_text.size() || index < 0) return false;
        const auto suffix = value.substr(close + 2);
        if (suffix != "isSymbolic()" && suffix != "isReference()" &&
            suffix != "isNoCall()" && suffix != "isCalled()" &&
            suffix != "isBreakpoint()" && suffix != "isSingleBreakend()" &&
            suffix != "isNonReference()" && suffix != "isNonRefAllele()" &&
            suffix != "length()")
            return false;
        reference_allele = false;
        allele_index = index;
        method = suffix;
        return true;
    } catch (...) {
        return false;
    }
}

bool select_read_allele_method(const std::string& raw, bcf1_t* record,
                               double& numeric_value, bool& boolean_value,
                               bool& numeric) {
    bool reference_allele = false;
    int allele_index = -1;
    std::string method;
    if (!parse_select_allele_method(raw, reference_allele, allele_index, method)) return false;
    bcf_unpack(record, BCF_UN_STR);
    const int index = reference_allele ? 0 : allele_index + 1;
    if (record->d.allele == nullptr || index < 0 || index >= record->n_allele ||
        record->d.allele[index] == nullptr) return false;
    const std::string allele(record->d.allele[index]);
    if (method == "length()") {
        const bool symbolic = select_allele_is_symbolic(allele);
        numeric_value = static_cast<double>(allele == "." || symbolic || allele == "*" ? 0 : allele.size());
        numeric = true;
        return allele != ".";
    }
    numeric = false;
    if (method == "isSymbolic()")
        boolean_value = select_allele_is_symbolic(allele);
    else if (method == "isReference()") boolean_value = reference_allele;
    else if (method == "isNoCall()") boolean_value = allele == ".";
    else if (method == "isCalled()") boolean_value = allele != ".";
    else if (method == "isBreakpoint()") boolean_value = select_allele_is_breakpoint(allele);
    else if (method == "isSingleBreakend()") boolean_value = select_allele_is_single_breakend(allele);
    else if (method == "isNonReference()") boolean_value = !reference_allele && allele != ".";
    else if (method == "isNonRefAllele()") boolean_value = allele == "<NON_REF>";
    else return false;
    return true;
}

bool select_string_method(const std::string& raw, const bcf_hdr_t* header,
                          bcf1_t* record, bool& value) {
    const auto source = trim_copy(raw);
    struct Method { const char* token; int kind; };
    constexpr Method methods[] = {
        {".contains(", 0}, {".startsWith(", 1}, {".endsWith(", 2}, {".matches(", 3},
    };
    for (const auto& method : methods) {
        const auto position = source.rfind(method.token);
        if (position == std::string::npos || source.empty() || source.back() != ')') continue;
        const auto argument_begin = position + std::strlen(method.token);
        if (argument_begin >= source.size() - 1) return false;
        const auto argument = unquote(source.substr(argument_begin, source.size() - argument_begin - 1));
        std::string observed;
        if ((source.substr(0, position) == "vc.getFilters()" || source.substr(0, position) == "getFilters()") &&
            method.kind == 0)
            return value = select_filter_set_contains(header, record, argument), true;
        if (!select_read_string_field(header, record, source.substr(0, position), observed)) return false;
        if (method.kind == 0) value = observed.find(argument) != std::string::npos;
        else if (method.kind == 1) value = observed.rfind(argument, 0) == 0;
        else if (method.kind == 2) value = observed.size() >= argument.size() &&
            observed.compare(observed.size() - argument.size(), argument.size(), argument) == 0;
        else {
            try { value = std::regex_search(observed, std::regex(argument)); }
            catch (const std::regex_error&) { throw std::invalid_argument("invalid regex in SelectVariants expression: " + raw); }
        }
        return true;
    }
    return false;
}

bool select_read_field(const bcf_hdr_t* header, bcf1_t* record,
                       const std::string& raw_field, double& value) {
    auto field = trim_copy(raw_field);
    int requested_index = -1;
    std::string indexed_base;
    if (parse_select_indexed_info_field(field, indexed_base, requested_index))
        field = std::move(indexed_base);
    if (field == "vc.getNAlleles()" || field == "getNAlleles()") {
        value = record->n_allele;
        return true;
    }
    if (field == "vc.getAlleles().size()" || field == "getAlleles().size()") {
        value = record->n_allele;
        return true;
    }
    if (field == "vc.getFilters().size()" || field == "getFilters().size()") {
        bcf_unpack(record, BCF_UN_FLT);
        value = record->d.n_flt;
        return true;
    }
    if (field == "vc.getGenotypes().size()" || field == "getGenotypes().size()") {
        value = record->n_sample;
        return true;
    }
    if (field == "vc.getNSamples()" || field == "getNSamples()") {
        value = record->n_sample;
        return true;
    }
    if (field == "vc.getCalledChrCount()" || field == "getCalledChrCount()" ||
        field == "vc.getNoCallCount()" || field == "getNoCallCount()" ||
        field == "vc.getHomRefCount()" || field == "getHomRefCount()" ||
        field == "vc.getHetCount()" || field == "getHetCount()" ||
        field == "vc.getHomVarCount()" || field == "getHomVarCount()") {
        SelectGenotypeCounts counts;
        if (!select_genotype_counts(header, record, counts)) return false;
        if (field.find("CalledChrCount") != std::string::npos) value = counts.called_chromosomes;
        else if (field.find("NoCallCount") != std::string::npos) value = counts.no_call_count;
        else if (field.find("HomRefCount") != std::string::npos) value = counts.hom_ref_count;
        else if (field.find("HetCount") != std::string::npos) value = counts.het_count;
        else value = counts.hom_var_count;
        return true;
    }
    if (field == "vc.getStart()" || field == "getStart()") {
        value = static_cast<double>(record->pos + 1);
        return record->pos >= 0;
    }
    if (field == "vc.getEnd()" || field == "getEnd()") {
        int32_t* end = nullptr;
        int count = 0;
        const auto length = bcf_get_info_int32(header, record, "END", &end, &count);
        if (length > 0 && count > 0 && end[0] != bcf_int32_missing && end[0] != bcf_int32_vector_end) {
            value = end[0];
            free(end);
            return true;
        }
        free(end);
        bcf_unpack(record, BCF_UN_STR);
        if (record->d.allele == nullptr || record->d.allele[0] == nullptr) return false;
        value = static_cast<double>(record->pos + std::strlen(record->d.allele[0]));
        return record->pos >= 0;
    }
    const std::string prefix = "vc.getAttribute(";
    if (field.rfind(prefix, 0) == 0 && !field.empty() && field.back() == ')')
        field = unquote(field.substr(prefix.size(), field.size() - prefix.size() - 1));
    if (field == "QUAL") {
        if (bcf_float_is_missing(record->qual) || bcf_float_is_vector_end(record->qual)) return false;
        value = record->qual;
        return true;
    }
    float* floats = nullptr;
    int float_count = 0;
    const auto float_length = bcf_get_info_float(header, record, field.c_str(), &floats, &float_count);
    const int float_index = requested_index >= 0 ? requested_index : 0;
    if (float_length > float_index && float_count > float_index &&
        !bcf_float_is_missing(floats[float_index]) &&
        !bcf_float_is_vector_end(floats[float_index])) {
        value = floats[float_index];
        free(floats);
        return std::isfinite(value);
    }
    free(floats);
    int32_t* integers = nullptr;
    int integer_count = 0;
    const auto integer_length = bcf_get_info_int32(header, record, field.c_str(), &integers, &integer_count);
    const int integer_index = requested_index >= 0 ? requested_index : 0;
    if (integer_length <= integer_index || integer_count <= integer_index ||
        integers[integer_index] == bcf_int32_missing ||
        integers[integer_index] == bcf_int32_vector_end) {
        free(integers);
        return false;
    }
    value = integers[integer_index];
    free(integers);
    return true;
}

bool evaluate_select_numeric_expression(const SelectNumericExpr& expression,
                                        const bcf_hdr_t* header, bcf1_t* record,
                                        double& value) {
    if (expression.kind == SelectNumericExpr::Kind::Literal) {
        value = expression.literal;
        return true;
    }
    if (expression.kind == SelectNumericExpr::Kind::Field)
        return select_read_field(header, record, expression.field, value);
    if (expression.kind == SelectNumericExpr::Kind::Negate) {
        if (!expression.left || !evaluate_select_numeric_expression(*expression.left, header, record, value))
            return false;
        value = -value;
        return true;
    }
    if (!expression.left || !expression.right) return false;
    double left = 0.0;
    double right = 0.0;
    if (!evaluate_select_numeric_expression(*expression.left, header, record, left) ||
        !evaluate_select_numeric_expression(*expression.right, header, record, right))
        return false;
    switch (expression.kind) {
        case SelectNumericExpr::Kind::Add: value = left + right; break;
        case SelectNumericExpr::Kind::Subtract: value = left - right; break;
        case SelectNumericExpr::Kind::Multiply: value = left * right; break;
        case SelectNumericExpr::Kind::Divide: value = left / right; break;
        case SelectNumericExpr::Kind::Modulo: value = std::fmod(left, right); break;
        default: return false;
    }
    return true;
}

// Forward declaration: ordinary --select expressions may reuse the
// sample-aware genotype evaluator defined below.
bool select_genotype_expression(const std::string& source, const bcf_hdr_t* header,
                                bcf1_t* record, int sample);

bool select_expression(const std::string& source, const bcf_hdr_t* header, bcf1_t* record) {
    auto value = trim_copy(source);
    while (wrapped_by_parentheses(value)) value = trim_copy(value.substr(1, value.size() - 2));
    if (value.empty()) throw std::invalid_argument("empty SelectVariants expression");
    const auto or_position = top_level_token(value, "||");
    if (or_position != std::string::npos)
        return select_expression(value.substr(0, or_position), header, record) ||
               select_expression(value.substr(or_position + 2), header, record);
    const auto and_position = top_level_token(value, "&&");
    if (and_position != std::string::npos)
        return select_expression(value.substr(0, and_position), header, record) &&
               select_expression(value.substr(and_position + 2), header, record);
    if (value.front() == '!') return !select_expression(value.substr(1), header, record);

    // `--select` accepts the same explicit sample receiver as GATK's JEXL
    // surface, independently of the any-sample `--select-genotype` option.
    // Strip the receiver once, resolve the named sample through the VCF
    // header, and reuse the compact genotype evaluator so comparisons,
    // indexed AD/PL and boolean predicates cannot diverge between options.
    const std::string genotype_prefix = "vc.getGenotype(";
    if (value.rfind(genotype_prefix, 0) == 0) {
        const auto close = value.find(')', genotype_prefix.size());
        if (close == std::string::npos)
            throw std::invalid_argument("malformed SelectVariants genotype receiver: " + source);
        const auto sample_name = unquote(value.substr(genotype_prefix.size(), close - genotype_prefix.size()));
        const auto sample_index = bcf_hdr_id2int(header, BCF_DT_SAMPLE, sample_name.c_str());
        if (sample_index < 0) return false;
        auto suffix = trim_copy(value.substr(close + 1));
        if (!suffix.empty() && suffix.front() == '.') suffix.erase(0, 1);
        if (suffix.empty())
            throw std::invalid_argument("missing SelectVariants genotype method: " + source);
        return select_genotype_expression(suffix, header, record, sample_index);
    }

    double allele_numeric = 0.0;
    bool allele_boolean = false;
    bool allele_is_numeric = false;
    if (select_read_allele_method(value, record, allele_numeric, allele_boolean, allele_is_numeric)) {
        if (allele_is_numeric)
            throw std::invalid_argument("allele length requires a numeric comparison: " + source);
        return allele_boolean;
    }

    bool string_method_value = false;
    if (select_string_method(value, header, record, string_method_value))
        return string_method_value;

    if (value.back() == ')' && value.find('(') != std::string::npos) {
        bool method_value = false;
        if (select_site_method(value, record, header, method_value)) return method_value;
    }

    const std::vector<std::string> operators{"<=", ">=", "==", "!=", "=~", "!~", "<", ">"};
    std::size_t position = std::string::npos;
    std::string op;
    for (const auto& candidate : operators) {
        const auto found = top_level_token(value, candidate);
        if (found != std::string::npos && (position == std::string::npos || found < position)) {
            position = found;
            op = candidate;
        }
    }
    if (position == std::string::npos)
        throw std::invalid_argument("unsupported SelectVariants expression: " + source);
    const auto lhs = trim_copy(value.substr(0, position));
    const auto rhs = trim_copy(value.substr(position + op.size()));
    // Lower compound numeric operands before the string/boolean special
    // cases.  Simple fields continue through the established readers, while
    // expressions such as `QUAL / 2 > 10` and `DP + 5 == 25` use one
    // deterministic evaluator for both sides of the comparison.
    if (rhs != "null" && rhs != "true" && rhs != "false" &&
        !(rhs.size() >= 2 && ((rhs.front() == '"' && rhs.back() == '"') ||
                              (rhs.front() == '\'' && rhs.back() == '\'')))) {
        auto lhs_numeric = parse_select_numeric_expression(lhs);
        auto rhs_numeric = parse_select_numeric_expression(rhs);
        if (select_numeric_expression_has_operator(*lhs_numeric) ||
            select_numeric_expression_has_operator(*rhs_numeric)) {
            // SelectVariants' HTSJDK JEXL evaluator rejects collection
            // elements as arithmetic operands (the indexed value must be the
            // direct comparison operand).  Fail closed instead of accepting
            // a native expression whose result would diverge from GATK.
            const auto contains_indexed_collection = [](const std::string& operand) {
                return operand.find(".get(") != std::string::npos ||
                       operand.find('[') != std::string::npos ||
                       operand.find(").") != std::string::npos;
            };
            if (contains_indexed_collection(lhs) || contains_indexed_collection(rhs))
                throw std::invalid_argument("SelectVariants does not support arithmetic on indexed INFO/FORMAT vectors: " + source);
            double observed = 0.0;
            double threshold = 0.0;
            if (!evaluate_select_numeric_expression(*lhs_numeric, header, record, observed) ||
                !evaluate_select_numeric_expression(*rhs_numeric, header, record, threshold))
                return false;
            return select_compare(observed, op, threshold);
        }
    }
    std::string null_lhs = lhs;
    int null_index = -1;
    std::string null_base;
    const bool lhs_is_indexed = parse_select_indexed_info_field(lhs, null_base, null_index);
    if (lhs_is_indexed) null_lhs = null_base;
    if (rhs == "null" && (null_lhs.rfind("vc.getAttribute(", 0) == 0 || null_lhs == "QUAL")) {
        if (lhs == "QUAL") {
            const bool present = !bcf_float_is_missing(record->qual) &&
                !bcf_float_is_vector_end(record->qual);
            return op == "==" ? !present : present;
        }
        if (lhs_is_indexed) {
            double ignored = 0.0;
            const bool present = select_read_field(header, record, lhs, ignored);
            return op == "==" ? !present : present;
        }
        const auto attribute_begin = null_lhs.find('(');
        const auto attribute_end = null_lhs.rfind(')');
        if (null_lhs.rfind("vc.getAttribute(", 0) == 0 && attribute_begin != std::string::npos &&
            attribute_end > attribute_begin) {
            const auto attribute = unquote(null_lhs.substr(attribute_begin + 1,
                attribute_end - attribute_begin - 1));
            const bool present = select_has_attribute(header, record, attribute);
            return op == "==" ? !present : present;
        }
    }
    // Boolean VariantContext methods can be compared explicitly, e.g.
    // vc.isSNP() == true.  The direct method path above only handles a
    // predicate with no comparison, so lower the left-hand method here.
    if (rhs == "true" || rhs == "false") {
        if (op != "==" && op != "!=")
            throw std::invalid_argument("unsupported boolean SelectVariants expression: " + source);
        bool method_value = false;
        if (select_site_method(lhs, record, header, method_value))
            return op == "==" ? method_value == (rhs == "true")
                               : method_value != (rhs == "true");
    }
    if ((rhs.size() >= 2 && ((rhs.front() == '"' && rhs.back() == '"') ||
                             (rhs.front() == '\'' && rhs.back() == '\'')))) {
        if (op != "==" && op != "!=" && op != "=~" && op != "!~")
            throw std::invalid_argument("unsupported string SelectVariants expression: " + source);
        std::string observed;
        if (!select_read_string_field(header, record, lhs, observed))
            // Match GATK's JEXL null coercion for regex operators: a missing
            // String operand cannot match (=~ false), but it does satisfy the
            // negated form (!~ true).
            return op == "!~";
        if (op == "==") return observed == unquote(rhs);
        if (op == "!=") return observed != unquote(rhs);
        auto pattern = unquote(rhs);
        if (pattern.size() >= 2 && pattern.front() == '/' && pattern.back() == '/')
            pattern = pattern.substr(1, pattern.size() - 2);
        try {
            const bool matched = std::regex_search(observed, std::regex(pattern));
            return op == "=~" ? matched : !matched;
        }
        catch (const std::regex_error&) {
            throw std::invalid_argument("invalid regex in SelectVariants expression: " + source);
        }
    }
    // GATK expressions commonly use the enum spelling without quotes, e.g.
    // vc.getType() == VariantContext.Type.SNP.
    if (op == "==" || op == "!=") {
        const std::string enum_prefix = "VariantContext.Type.";
        if (rhs.rfind(enum_prefix, 0) == 0) {
            std::string observed;
            if (!select_read_string_field(header, record, lhs, observed)) return false;
            const auto expected = rhs.substr(enum_prefix.size());
            return op == "==" ? observed == expected : observed != expected;
        }
    }
    bool ignored_boolean = false;
    bool numeric = false;
    double observed = 0.0;
    if (select_read_allele_method(lhs, record, observed, ignored_boolean, numeric)) {
        if (!numeric) {
            if ((rhs == "true" || rhs == "false") && (op == "==" || op == "!="))
                return op == "==" ? ignored_boolean == (rhs == "true")
                                   : ignored_boolean != (rhs == "true");
            throw std::invalid_argument("only allele length supports numeric comparison: " + source);
        }
        std::size_t allele_consumed = 0;
        double threshold = 0.0;
        try { threshold = std::stod(rhs, &allele_consumed); }
        catch (...) { throw std::invalid_argument("unsupported allele comparison: " + source); }
        if (allele_consumed != rhs.size())
            throw std::invalid_argument("unsupported allele comparison: " + source);
        return select_compare(observed, op, threshold);
    }
    std::size_t consumed = 0;
    double threshold = 0.0;
    try { threshold = std::stod(rhs, &consumed); }
    catch (...) { throw std::invalid_argument("unsupported non-numeric SelectVariants expression: " + source); }
    if (consumed != rhs.size()) throw std::invalid_argument("unsupported SelectVariants expression: " + source);
    double ordinary_observed = 0.0;
    return select_read_field(header, record, lhs, ordinary_observed) &&
           select_compare(ordinary_observed, op, threshold);
}

bool genotype_predicate(const bcf_hdr_t* header, bcf1_t* record, int sample,
                        std::string value) {
    value = trim_copy(std::move(value));
    const auto open = value.find('(');
    if (open == std::string::npos || value.back() != ')') return false;
    const auto name = value.substr(0, open);
    if (value.substr(open + 1, value.size() - open - 2) != "") return false;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || record->n_sample <= 0 || genotype_count % record->n_sample != 0 ||
        sample < 0 || sample >= record->n_sample) {
        free(genotypes);
        return false;
    }
    if (name == "isAvailable" || name == "g.isAvailable") {
        free(genotypes);
        return true;
    }
    if (name == "hasDP" || name == "g.hasDP" || name == "hasGQ" || name == "g.hasGQ" ||
        name == "hasAD" || name == "g.hasAD" || name == "hasPL" || name == "g.hasPL") {
        free(genotypes);
        std::string tag;
        if (name.find("DP") != std::string::npos) tag = "DP";
        else if (name.find("GQ") != std::string::npos) tag = "GQ";
        else if (name.find("AD") != std::string::npos) tag = "AD";
        else tag = "PL";
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(header, record, tag.c_str(), &values, &count);
        bool available = false;
        if (length > 0 && count >= record->n_sample && count % record->n_sample == 0) {
            const int width = count / record->n_sample;
            const int offset = sample * width;
            available = width > 0 && offset < count &&
                values[offset] != bcf_int32_missing && values[offset] != bcf_int32_vector_end;
        }
        free(values);
        return available;
    }
    const int ploidy = genotype_count / record->n_sample;
    bool called = true;
    bool any_alt = false;
    bool all_same = true;
    int first_allele = -1;
    for (int copy = 0; copy < ploidy; ++copy) {
        const auto encoded = genotypes[sample * ploidy + copy];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
            called = false;
            continue;
        }
        const int allele = bcf_gt_allele(encoded);
        if (first_allele < 0) first_allele = allele;
        else if (allele != first_allele) all_same = false;
        if (allele > 0) any_alt = true;
    }
    free(genotypes);
    if (name == "isNoCall" || name == "g.isNoCall") return !called;
    if (name == "isCalled" || name == "g.isCalled") return called;
    if (name == "isHet" || name == "g.isHet") return called && !all_same;
    if (name == "isHom" || name == "g.isHom") return called && all_same;
    if (name == "isHomRef" || name == "g.isHomRef") return called && all_same && first_allele == 0;
    if (name == "isHomVar" || name == "g.isHomVar") return called && all_same && any_alt;
    return false;
}

bool genotype_read_field(const bcf_hdr_t* header, bcf1_t* record, int sample,
                         std::string raw_field, double& value) {
    auto field = trim_copy(std::move(raw_field));
    int requested_index = -1;
    const auto open = field.find('[');
    if (open != std::string::npos && field.back() == ']') {
        try {
            requested_index = std::stoi(field.substr(open + 1, field.size() - open - 2));
        } catch (...) {
            return false;
        }
        field = field.substr(0, open);
    } else {
        // JEXL's collection shorthand uses `.N` after getAD()/getPL().
        // Parse it before stripping any receiver prefix below.
        const auto dot_index = field.rfind('.');
        if (dot_index != std::string::npos && dot_index + 1 < field.size() &&
            std::all_of(field.begin() + static_cast<std::ptrdiff_t>(dot_index + 1), field.end(),
                        [](unsigned char character) { return std::isdigit(character) != 0; })) {
            try {
                std::size_t consumed = 0;
                requested_index = std::stoi(field.substr(dot_index + 1), &consumed);
                if (consumed != field.size() - dot_index - 1 || requested_index < 0) return false;
                field = field.substr(0, dot_index);
            } catch (...) {
                return false;
            }
        }
    }
    const auto dot = field.rfind(".");
    if (dot != std::string::npos && dot + 1 < field.size()) field = field.substr(dot + 1);
    if (field == "getDP()") field = "DP";
    else if (field == "getGQ()") field = "GQ";
    else if (field == "getAD()") field = "AD";
    else if (field == "getPL()") field = "PL";
    else if (field == "getMIN_DP()") field = "MIN_DP";
    else if (field == "getPhredScaledQual()") field = "GQ";
    else if (field == "getPloidy()") {
        int32_t* genotypes = nullptr;
        int genotype_count = 0;
        const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
        if (length <= 0 || record->n_sample <= 0 || genotype_count % record->n_sample != 0) {
            free(genotypes);
            return false;
        }
        value = genotype_count / record->n_sample;
        free(genotypes);
        return true;
    }
    if (record->n_sample <= 0 || sample < 0 || sample >= record->n_sample) return false;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, field.c_str(), &values, &count);
    if (length <= 0 || count < record->n_sample || count % record->n_sample != 0) {
        free(values);
        return false;
    }
    const int width = count / record->n_sample;
    const int index = requested_index < 0 ? 0 : requested_index;
    if (index >= width) {
        free(values);
        return false;
    }
    const auto observed = values[sample * width + index];
    const bool missing = observed == bcf_int32_missing || observed == bcf_int32_vector_end;
    if (!missing) value = observed;
    free(values);
    return !missing;
}

bool select_genotype_expression(const std::string& source, const bcf_hdr_t* header,
                                bcf1_t* record, int sample) {
    auto value = trim_copy(source);
    while (wrapped_by_parentheses(value)) value = trim_copy(value.substr(1, value.size() - 2));
    if (value.empty()) throw std::invalid_argument("empty SelectVariants genotype expression");
    // Accept both SelectVariants' compact genotype syntax (GQ/AD/isHet())
    // and the explicit VariantContext form used by GATK JEXL expressions.
    const std::string genotype_prefix = "vc.getGenotype(";
    if (value.rfind(genotype_prefix, 0) == 0) {
        const auto close = value.find(')', genotype_prefix.size());
        if (close == std::string::npos)
            throw std::invalid_argument("malformed SelectVariants genotype expression: " + source);
        const auto requested_sample = unquote(value.substr(genotype_prefix.size(), close - genotype_prefix.size()));
        const auto* current_sample = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
        if (current_sample == nullptr || requested_sample != current_sample) return false;
        value = trim_copy(value.substr(close + 1));
        if (!value.empty() && value.front() == '.') value.erase(0, 1);
        if (value.empty()) throw std::invalid_argument("missing genotype method: " + source);
    }
    const auto or_position = top_level_token(value, "||");
    if (or_position != std::string::npos)
        return select_genotype_expression(value.substr(0, or_position), header, record, sample) ||
               select_genotype_expression(value.substr(or_position + 2), header, record, sample);
    const auto and_position = top_level_token(value, "&&");
    if (and_position != std::string::npos)
        return select_genotype_expression(value.substr(0, and_position), header, record, sample) &&
               select_genotype_expression(value.substr(and_position + 2), header, record, sample);
    if (value.front() == '!') return !select_genotype_expression(value.substr(1), header, record, sample);
    if (value.back() == ')' && value.find('(') != std::string::npos) {
        const auto open = value.find('(');
        const auto name = value.substr(0, open);
        static const std::set<std::string> predicates{
            "isAvailable", "g.isAvailable", "isNoCall", "g.isNoCall",
            "isCalled", "g.isCalled", "isHet", "g.isHet", "isHom", "g.isHom",
            "isHomRef", "g.isHomRef", "isHomVar", "g.isHomVar",
            "hasDP", "g.hasDP", "hasGQ", "g.hasGQ", "hasAD", "g.hasAD",
            "hasPL", "g.hasPL"};
        if (predicates.count(name) != 0)
            return genotype_predicate(header, record, sample, value);
    }
    const std::vector<std::string> operators{"<=", ">=", "==", "!=", "<", ">"};
    std::size_t position = std::string::npos;
    std::string op;
    for (const auto& candidate : operators) {
        const auto found = top_level_token(value, candidate);
        if (found != std::string::npos && (position == std::string::npos || found < position)) {
            position = found;
            op = candidate;
        }
    }
    if (position == std::string::npos)
        throw std::invalid_argument("unsupported SelectVariants genotype expression: " + source);
    const auto lhs = trim_copy(value.substr(0, position));
    const auto rhs = trim_copy(value.substr(position + op.size()));
    const bool boolean_rhs = rhs == "true" || rhs == "false";
    if (boolean_rhs && op != "==" && op != "!=")
        throw std::invalid_argument("unsupported boolean SelectVariants genotype expression: " + source);
    if (!boolean_rhs && rhs != "null") {
        auto lhs_numeric = parse_select_numeric_expression(lhs);
        auto rhs_numeric = parse_select_numeric_expression(rhs);
        if (select_numeric_expression_has_operator(*lhs_numeric) ||
            select_numeric_expression_has_operator(*rhs_numeric)) {
            // Genotype arithmetic is evaluated against the current sample;
            // use a tiny wrapper so each leaf retains AD/PL indexing and the
            // same missing-value behavior as the scalar genotype path.
            auto evaluate_genotype_numeric = [&](const SelectNumericExpr& expression,
                                                  auto&& self, double& observed) -> bool {
                if (expression.kind == SelectNumericExpr::Kind::Literal) {
                    observed = expression.literal;
                    return true;
                }
                if (expression.kind == SelectNumericExpr::Kind::Field)
                    return genotype_read_field(header, record, sample, expression.field, observed);
                if (expression.kind == SelectNumericExpr::Kind::Negate) {
                    if (!expression.left || !self(*expression.left, self, observed)) return false;
                    observed = -observed;
                    return true;
                }
                if (!expression.left || !expression.right) return false;
                double left = 0.0;
                double right = 0.0;
                if (!self(*expression.left, self, left) || !self(*expression.right, self, right)) return false;
                switch (expression.kind) {
                    case SelectNumericExpr::Kind::Add: observed = left + right; break;
                    case SelectNumericExpr::Kind::Subtract: observed = left - right; break;
                    case SelectNumericExpr::Kind::Multiply: observed = left * right; break;
                    case SelectNumericExpr::Kind::Divide: observed = left / right; break;
                    case SelectNumericExpr::Kind::Modulo: observed = std::fmod(left, right); break;
                    default: return false;
                }
                return true;
            };
            double observed = 0.0;
            double threshold_value = 0.0;
            if (!evaluate_genotype_numeric(*lhs_numeric, evaluate_genotype_numeric, observed) ||
                !evaluate_genotype_numeric(*rhs_numeric, evaluate_genotype_numeric, threshold_value))
                return false;
            return select_compare(observed, op, threshold_value);
        }
    }
    std::size_t consumed = 0;
    double threshold = 0.0;
    try {
        threshold = boolean_rhs ? (rhs == "true" ? 1.0 : 0.0) : std::stod(rhs, &consumed);
        if (boolean_rhs) consumed = rhs.size();
    }
    catch (...) { throw std::invalid_argument("unsupported non-numeric SelectVariants genotype expression: " + source); }
    if (consumed != rhs.size()) throw std::invalid_argument("unsupported SelectVariants genotype expression: " + source);
    // VariantFiltration/SelectVariants share GATK's compact numeric genotype
    // predicates (`isHet == 1`, `isHomRef == 1`, ...).  Lower the boolean
    // result to 0/1 before applying the requested numeric comparison.
    static const std::set<std::string> boolean_predicates{
        "isAvailable", "g.isAvailable", "isNoCall", "g.isNoCall",
        "isCalled", "g.isCalled", "isHet", "g.isHet", "isHom", "g.isHom",
        "isHomRef", "g.isHomRef", "isHomVar", "g.isHomVar",
        "hasDP", "g.hasDP", "hasGQ", "g.hasGQ", "hasAD", "g.hasAD",
        "hasPL", "g.hasPL"
    };
    auto predicate_name = lhs;
    if (predicate_name.size() > 2 && predicate_name.ends_with("()"))
        predicate_name.erase(predicate_name.size() - 2);
    if (boolean_predicates.contains(predicate_name)) {
        const bool predicate = genotype_predicate(header, record, sample, predicate_name + "()");
        return select_compare(predicate ? 1.0 : 0.0, op, threshold);
    }
    double observed = 0.0;
    return genotype_read_field(header, record, sample, lhs, observed) &&
           select_compare(observed, op, threshold);
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

void compact_fields(const bcf_hdr_t* header, bcf1_t* record, const std::vector<int>& kept,
                    int old_alleles, GenotypeKernelTelemetry& telemetry) {
    const int sample_count = record->n_sample;
    if (sample_count <= 0) return;
    int ploidy = 2;
    bool genotype_ploidy_observed = false;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    if (bcf_get_genotypes(header, record, &genotypes, &genotype_count) > 0 &&
        genotype_count >= sample_count && genotype_count % sample_count == 0) {
        const int width = genotype_count / sample_count;
        if (width > 0) {
            ploidy = width;
            genotype_ploidy_observed = true;
        }
        for (int index = 0; index < genotype_count; ++index) {
            if (bcf_gt_is_missing(genotypes[index]) || genotypes[index] == bcf_int32_vector_end) {
                genotypes[index] = bcf_gt_missing;
                continue;
            }
            const int old_allele = bcf_gt_allele(genotypes[index]);
            const auto found = std::find(kept.begin(), kept.end(), old_allele);
            const bool phased = bcf_gt_is_phased(genotypes[index]);
            genotypes[index] = found == kept.end() ? bcf_gt_missing :
                (phased ? bcf_gt_phased(static_cast<int>(found - kept.begin())) :
                          bcf_gt_unphased(static_cast<int>(found - kept.begin())));
        }
        if (bcf_update_genotypes(header, record, genotypes, genotype_count) != 0) {
            free(genotypes);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset GT FORMAT");
        }
        (void)width;
    }
    free(genotypes);
    // PL is Number=G and remains useful even when a VCF omits GT.  Infer the
    // unique bounded ploidy from the triangular width instead of silently
    // assuming diploid, so allele subsetting preserves arbitrary fixed
    // ploidy in PL-only records.
    if (!genotype_ploidy_observed && bcf_hdr_id2int(header, BCF_DT_ID, "PL") >= 0) {
        int32_t* pls = nullptr;
        int pl_count = 0;
        const auto pl_length = bcf_get_format_int32(header, record, "PL", &pls, &pl_count);
        if (pl_length > 0 && pl_count >= sample_count && pl_count % sample_count == 0) {
            const auto width = static_cast<std::size_t>(pl_count / sample_count);
            for (int candidate = 1; candidate <= 32; ++candidate) {
                if (genotype_width(old_alleles, candidate) == width) {
                    ploidy = candidate;
                    break;
                }
            }
        }
        free(pls);
    }
    const auto old_pl_width = genotype_width(old_alleles, ploidy);
    const int new_alleles = static_cast<int>(kept.size());
    const auto new_pl_width = genotype_width(new_alleles, ploidy);
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AD") >= 0) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
        if (length > 0 && count >= sample_count * old_alleles && count % sample_count == 0) {
            std::vector<int32_t> kept_indices;
            kept_indices.reserve(kept.size());
            for (const auto allele : kept) kept_indices.push_back(allele);
            const std::vector<int32_t> source(
                values, values + static_cast<std::size_t>(sample_count) * old_alleles);
            const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
                source, static_cast<std::size_t>(sample_count), old_alleles,
                new_alleles, kept_indices);
            if (bcf_update_format_int32(header, record, "AD", remapped.values.data(),
                                        static_cast<int>(remapped.values.size())) != 0) {
                free(values);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset AD FORMAT");
            }
            ++telemetry.allele_field_remap_calls;
            telemetry.allele_field_remap_prepare_seconds += remapped.prepare_seconds;
            telemetry.allele_field_remap_execute_seconds += remapped.seconds;
            telemetry.allele_field_remap_execution_space = remapped.execution_space;
        }
        free(values);
    }
    if (bcf_hdr_id2int(header, BCF_DT_ID, "PL") >= 0 && old_pl_width > 0 && new_pl_width > 0) {
        int32_t* pls = nullptr;
        int pl_count = 0;
        const auto length = bcf_get_format_int32(header, record, "PL", &pls, &pl_count);
        if (length > 0 && pl_count >= static_cast<int>(sample_count * old_pl_width)) {
            if (new_alleles == 1) {
                // A reference-only allele space has exactly one Number=G cell:
                // the all-reference likelihood, which is old PL index 0 because
                // AlleleSubsettingUtils.subsettedPLIndices maps the new
                // (ref,ref,..) genotype onto oldAlleleCounts (0,0,..).  The
                // Kokkos remap kernel is defined for two or more target alleles,
                // so collapse the vector directly.  GQ is left to the caller:
                // AlleleSubsettingUtils.java:98-99 keeps the original GQ for
                // this subset, and compact_fields must not re-derive it here.
                std::vector<std::int32_t> collapsed(static_cast<std::size_t>(sample_count));
                for (int sample = 0; sample < sample_count; ++sample)
                    collapsed[static_cast<std::size_t>(sample)] =
                        pls[static_cast<std::size_t>(sample) * static_cast<std::size_t>(old_pl_width)];
                if (bcf_update_format_int32(header, record, "PL", collapsed.data(), sample_count) != 0) {
                    free(pls);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot collapse SelectVariants PL FORMAT");
                }
            } else {
                std::vector<std::int32_t> source(
                    pls, pls + sample_count * old_pl_width);
                std::vector<std::int32_t> kept_indices;
                kept_indices.reserve(kept.size());
                for (const auto allele : kept)
                    kept_indices.push_back(static_cast<std::int32_t>(allele));
                const auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
                    source, static_cast<std::size_t>(sample_count), old_alleles,
                    new_alleles, ploidy, kept_indices);
                if (bcf_update_format_int32(header, record, "PL", remapped.pl.data(),
                                            static_cast<int>(remapped.pl.size())) != 0) {
                    free(pls);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset PL FORMAT");
                }
                ++telemetry.pl_remap_calls;
                telemetry.pl_remap_prepare_seconds += remapped.prepare_seconds;
                telemetry.pl_remap_execute_seconds += remapped.seconds;
                telemetry.execution_space = remapped.execution_space;
                if (bcf_hdr_id2int(header, BCF_DT_ID, "GQ") >= 0) {
                    const auto derived = fastgatk::kernels::derive_genotype_gt_gq_kokkos(
                        remapped.pl, static_cast<std::size_t>(sample_count),
                        new_alleles, ploidy);
                    if (bcf_update_format_int32(header, record, "GQ", derived.gq.data(), sample_count) != 0) {
                        free(pls);
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot recompute SelectVariants GQ");
                    }
                    ++telemetry.gt_gq_calls;
                    telemetry.gt_gq_prepare_seconds += derived.prepare_seconds;
                    telemetry.gt_gq_execute_seconds += derived.seconds;
                    telemetry.execution_space = derived.execution_space;
                }
            }
        }
        free(pls);
    }
}

// Pinned GATK 4.6.2.0 does *not* drop an all-hom-ref record when
// --remove-unused-alternates is enabled: SelectVariants.subsetGenotypesBySampleNames
// trims the unused ALTs and returns a single-allele VariantContext
// (SelectVariants.java:1214-1233 -> GATKVariantContextUtils.trimAlleles at
// GATKVariantContextUtils.java:1455), which the writer emits with ALT='.'.
// The non-variant record is only removed when --exclude-non-variants is set
// (SelectVariants.java:708-715).
//
// Two further observable details of that ref-only subset, both measured against
// the pinned jar:
//   * AlleleSubsettingUtils.java:94-100 keeps the *original* GQ when the new
//     allele subset has a single PL cell ("if we subset to just ref allele, keep
//     the GQ"), instead of re-deriving it from the collapsed PL vector, and it
//     only emits GQ when the input genotype had one.
//   * htsjdk's chromosome-count refresh then writes AN (including AN=0 for an
//     all-no-call genotype set) and no AC/AF, because the record has no ALT.
void trim_record_to_reference_only(const bcf_hdr_t* header, bcf1_t* record, int old_alleles,
                                   GenotypeKernelTelemetry& telemetry) {
    std::vector<std::int32_t> original_gq;
    int32_t* gq = nullptr;
    int gq_count = 0;
    const bool had_gq = bcf_get_format_int32(header, record, "GQ", &gq, &gq_count) > 0 &&
        gq != nullptr && gq_count >= record->n_sample;
    if (had_gq) original_gq.assign(gq, gq + record->n_sample);
    free(gq);

    bcf_unpack(record, BCF_UN_STR);
    const std::string reference_allele = record->d.allele[0] == nullptr ? "" : record->d.allele[0];
    if (reference_allele.empty() || reference_allele == "." ||
        bcf_update_alleles_str(header, record, reference_allele.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot trim SelectVariants record to the reference allele");

    const std::vector<int> kept{0};
    compact_fields(header, record, kept, old_alleles, telemetry);

    if (had_gq) {
        if (bcf_update_format_int32(header, record, "GQ", original_gq.data(),
                                    static_cast<int>(original_gq.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot preserve SelectVariants GQ for a reference-only record");
    } else if (bcf_hdr_id2int(header, BCF_DT_ID, "GQ") >= 0 &&
               bcf_update_format_int32(header, record, "GQ", nullptr, 0) != 0) {
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot drop unrequested SelectVariants GQ");
    }

    int32_t allele_number = 0;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    if (bcf_get_genotypes(header, record, &genotypes, &genotype_count) > 0) {
        for (int index = 0; index < genotype_count; ++index) {
            const auto encoded = genotypes[index];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            ++allele_number;
        }
    }
    free(genotypes);
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AC") >= 0 &&
        bcf_update_info_int32(header, record, "AC", nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear SelectVariants AC for a reference-only record");
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AF") >= 0 &&
        bcf_update_info_float(header, record, "AF", nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear SelectVariants AF for a reference-only record");
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AN") >= 0 &&
        bcf_update_info_int32(header, record, "AN", &allele_number, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot refresh SelectVariants AN for a reference-only record");
}

// SelectVariants applies --set-filtered-gt-to-nocall after sample subsetting
// and before --remove-unused-alternates.  Only called genotypes are eligible:
// an already no-call genotype must retain its existing payload and a PASS/.
// FT must not be treated as filtered.  This ordering is observable for a
// filtered genotype that is the sole carrier of an ALT, because the subsequent
// allele compaction then removes that ALT and remaps Number=G/R FORMAT fields.
std::uint64_t set_filtered_genotypes_to_nocall(const bcf_hdr_t* header,
                                                bcf1_t* record) {
    if (record->n_sample <= 0) return 0;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto genotype_length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (genotype_length <= 0 || genotypes == nullptr || genotype_count <= 0 ||
        genotype_count < record->n_sample || genotype_count % record->n_sample != 0) {
        free(genotypes);
        return 0;
    }
    const int ploidy = genotype_count / record->n_sample;
    if (ploidy <= 0) {
        free(genotypes);
        return 0;
    }

    char** filter_values = nullptr;
    int filter_count = 0;
    const auto filter_length = bcf_get_format_string(
        header, record, "FT", &filter_values, &filter_count);
    if (filter_length <= 0 || filter_values == nullptr || filter_count < record->n_sample) {
        if (filter_values != nullptr) {
            if (filter_count > 0) free(filter_values[0]);
            free(filter_values);
        }
        free(genotypes);
        return 0;
    }

    std::uint64_t changed = 0;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        const char* observed = filter_values[sample];
        const bool genotype_filtered = observed != nullptr && observed[0] != '\0' &&
            std::strcmp(observed, ".") != 0 && std::strcmp(observed, "PASS") != 0;
        if (!genotype_filtered) continue;
        bool called = true;
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto encoded = genotypes[sample * ploidy + allele];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                called = false;
                break;
            }
        }
        if (!called) continue;
        for (int allele = 0; allele < ploidy; ++allele) {
            // GenotypeBuilder(g).alleles(noCallAlleles) keeps the original
            // phased bit.  HTSlib represents that as a phased missing allele
            // (.|.) rather than flattening it to an unphased ./. 
            const auto encoded = genotypes[sample * ploidy + allele];
            genotypes[sample * ploidy + allele] = bcf_gt_is_phased(encoded)
                ? bcf_gt_phased(-1) : bcf_gt_missing;
        }
        ++changed;
    }

    if (changed > 0 && bcf_update_genotypes(header, record, genotypes, genotype_count) != 0) {
        if (filter_count > 0) free(filter_values[0]);
        free(filter_values);
        free(genotypes);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set filtered GT to no-call");
    }
    if (filter_count > 0) free(filter_values[0]);
    free(filter_values);
    free(genotypes);
    return changed;
}

// A VariantContext writer omits FORMAT/FT when every selected genotype has
// an unfiltered status (PASS or missing).  HTSlib otherwise carries a literal
// PASS vector through bcf_subset(), which is observably different in a
// sample-subset SelectVariants output.  Keep FT only when at least one
// non-PASS status is present in this record.
void clear_all_pass_genotype_filters(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->n_sample <= 0) return;
    char** filter_values = nullptr;
    int filter_count = 0;
    const auto filter_length = bcf_get_format_string(
        header, record, "FT", &filter_values, &filter_count);
    if (filter_length <= 0 || filter_values == nullptr || filter_count < record->n_sample) {
        if (filter_values != nullptr) {
            if (filter_count > 0) free(filter_values[0]);
            free(filter_values);
        }
        return;
    }
    bool has_filtered_genotype = false;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        const char* observed = filter_values[sample];
        if (observed != nullptr && observed[0] != '\0' &&
            std::strcmp(observed, ".") != 0 && std::strcmp(observed, "PASS") != 0) {
            has_filtered_genotype = true;
            break;
        }
    }
    if (!has_filtered_genotype && bcf_update_format_string(
            header, record, "FT", nullptr, 0) != 0) {
        if (filter_count > 0) free(filter_values[0]);
        free(filter_values);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear unfiltered genotype FT");
    }
    if (filter_count > 0) free(filter_values[0]);
    free(filter_values);
}

void recompute_site_annotations(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->n_allele < 2 || record->n_sample <= 0) return;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || genotype_count <= 0) {
        free(genotypes);
        return;
    }
    std::vector<int32_t> ac(static_cast<std::size_t>(record->n_allele - 1), 0);
    int32_t an = 0;
    for (int index = 0; index < genotype_count; ++index) {
        const auto encoded = genotypes[index];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
        const int allele = bcf_gt_allele(encoded);
        if (allele < 0 || allele >= record->n_allele) continue;
        ++an;
        if (allele > 0) ++ac[static_cast<std::size_t>(allele - 1)];
    }
    free(genotypes);
    // HTSJDK's chromosome-count refresh removes AC/AN/AF when the resulting
    // genotype set contains no called chromosomes (for example, a sole
    // phased ALT genotype converted to no-call).  Do not emit synthetic zero
    // annotations in that case.
    if (an == 0) {
        if (bcf_update_info_int32(header, record, "AC", nullptr, 0) != 0 ||
            bcf_update_info_int32(header, record, "AN", nullptr, 0) != 0 ||
            bcf_update_info_float(header, record, "AF", nullptr, 0) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear SelectVariants AC/AN/AF");
        return;
    }
    if (bcf_update_info_int32(header, record, "AC", ac.data(), static_cast<int>(ac.size())) != 0 ||
        bcf_update_info_int32(header, record, "AN", &an, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update SelectVariants AC/AN");
    std::vector<float> af(ac.size(), 0.0F);
    for (std::size_t index = 0; index < ac.size(); ++index)
        af[index] = an == 0 ? 0.0F : static_cast<float>(ac[index]) / static_cast<float>(an);
    if (bcf_update_info_float(header, record, "AF", af.data(), static_cast<int>(af.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update SelectVariants AF");
}

std::vector<int> sample_indices(const bcf_hdr_t* header, const Options& options,
                                std::string& sample_spec) {
    const int count = header->n[BCF_DT_SAMPLE];
    std::set<std::string> excluded(options.exclude_samples.begin(), options.exclude_samples.end());
    std::vector<int> indices;
    for (int index = 0; index < count; ++index) {
        const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
        if (name == nullptr) continue;
        const bool explicitly_included = options.include_samples.empty() ||
            std::find(options.include_samples.begin(), options.include_samples.end(), name) != options.include_samples.end();
        if (explicitly_included && !excluded.contains(name)) indices.push_back(index);
    }
    for (const auto& name : options.include_samples)
        if (bcf_hdr_id2int(header, BCF_DT_SAMPLE, name.c_str()) < 0)
            throw std::runtime_error("BAD_INPUT: sample not present in VCF header: " + name);
    for (const auto& name : options.exclude_samples)
        if (bcf_hdr_id2int(header, BCF_DT_SAMPLE, name.c_str()) < 0)
            throw std::runtime_error("BAD_INPUT: sample not present in VCF header: " + name);
    for (std::size_t position = 0; position < indices.size(); ++position) {
        if (position != 0) sample_spec.push_back(',');
        sample_spec += bcf_hdr_int2id(header, BCF_DT_SAMPLE, indices[position]);
    }
    return indices;
}

struct OriginalAnnotations {
    std::vector<std::int32_t> ac;
    std::vector<float> af;
    std::int32_t an = 0;
    std::int32_t dp = 0;
    bool has_ac = false;
    bool has_af = false;
    bool has_an = false;
    bool has_dp = false;
};

OriginalAnnotations capture_original_annotations(const bcf_hdr_t* header, bcf1_t* record,
                                                 bool keep_ac, bool keep_dp) {
    OriginalAnnotations result;
    if (keep_ac) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_info_int32(header, record, "AC", &values, &count);
        if (length > 0 && values != nullptr) {
            result.ac.assign(values, values + length);
            result.has_ac = true;
        }
        free(values);
        float* frequencies = nullptr;
        count = 0;
        const auto frequency_length = bcf_get_info_float(header, record, "AF", &frequencies, &count);
        if (frequency_length > 0 && frequencies != nullptr) {
            result.af.assign(frequencies, frequencies + frequency_length);
            result.has_af = true;
        }
        free(frequencies);
        values = nullptr;
        count = 0;
        const auto allele_number_length = bcf_get_info_int32(header, record, "AN", &values, &count);
        if (allele_number_length > 0 && values != nullptr) {
            result.an = values[0];
            result.has_an = true;
        }
        free(values);
    }
    if (keep_dp) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_info_int32(header, record, "DP", &values, &count);
        if (length > 0 && values != nullptr) {
            result.dp = values[0];
            result.has_dp = true;
        }
        free(values);
    }
    return result;
}

void restore_original_annotations(const bcf_hdr_t* header, bcf1_t* record,
                                  const Options& options, const OriginalAnnotations& original,
                                  const std::vector<int>& retained_alt_indices) {
    if (options.keep_original_ac) {
        if (original.has_ac) {
            std::vector<std::int32_t> values;
            values.reserve(retained_alt_indices.size());
            for (const auto allele : retained_alt_indices)
                if (allele > 0 && static_cast<std::size_t>(allele - 1) < original.ac.size())
                    values.push_back(original.ac[static_cast<std::size_t>(allele - 1)]);
            if (!values.empty() && bcf_update_info_int32(header, record, "AC_Orig",
                                                         values.data(), static_cast<int>(values.size())) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AC_Orig");
        }
        if (original.has_af) {
            std::vector<float> values;
            values.reserve(retained_alt_indices.size());
            for (const auto allele : retained_alt_indices)
                if (allele > 0 && static_cast<std::size_t>(allele - 1) < original.af.size())
                    values.push_back(original.af[static_cast<std::size_t>(allele - 1)]);
            if (!values.empty() && bcf_update_info_float(header, record, "AF_Orig",
                                                         values.data(), static_cast<int>(values.size())) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AF_Orig");
        }
        if (original.has_an && bcf_update_info_int32(header, record, "AN_Orig", &original.an, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AN_Orig");
    }
    if (options.keep_original_dp && original.has_dp &&
        bcf_update_info_int32(header, record, "DP_Orig", &original.dp, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write DP_Orig");
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    htsFile* input = nullptr;
    htsFile* output = nullptr;
    bcf_hdr_t* input_header = nullptr;
    bcf_hdr_t* output_header = nullptr;
    bcf_hdr_t* writer_header = nullptr;
    bcf1_t* record = nullptr;
    std::uint64_t input_records = 0, output_records = 0, filtered_records = 0, type_skipped = 0;
    std::uint64_t non_variant_skipped = 0, indel_size_skipped = 0;
    std::uint64_t filtered_gt_nocall = 0;
    std::uint64_t expression_skipped = 0;
    std::uint64_t comparison_skipped = 0;
    std::uint64_t interval_skipped = 0;
    GenotypeKernelTelemetry genotype_kernel_telemetry;
    fastgatk::io::IntervalFileStats interval_file_stats;
    std::string sample_spec;
    try {
        input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open variant input: " + options.input);
        input_header = bcf_hdr_read(input);
        if (!input_header) throw std::runtime_error("BAD_INPUT: cannot read variant header");
        std::vector<Region> regions;
        regions.reserve(options.regions.size());
        bool first_selector = true;
        for (const auto& text : options.regions) {
            append_select_region_selector_with_rule(
                text, input_header, regions, interval_file_stats,
                options.interval_set_rule, first_selector);
            first_selector = false;
        }
        normalize_select_regions(regions);
        const bool comparison_mode = !options.concordance.empty() || !options.discordance.empty();
        const auto comparison_data = comparison_mode
            ? load_comparison_data(options.concordance.empty() ? options.discordance : options.concordance)
            : ComparisonData{};
        if (comparison_mode && options.compare_genotypes) {
            bool shared_sample = false;
            for (int sample = 0; sample < input_header->n[BCF_DT_SAMPLE]; ++sample) {
                const auto* name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, sample);
                if (name != nullptr && comparison_data.sample_names.contains(name)) {
                    shared_sample = true;
                    break;
                }
            }
            if (!shared_sample)
                throw std::invalid_argument("BAD_INPUT: sample-level concordance requires a shared sample");
        }
        output_header = bcf_hdr_dup(input_header);
        if (!output_header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate variant header");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=AC,Number=A,Type=Integer,Description=Alternate allele count>");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=AN,Number=1,Type=Integer,Description=Total called allele count>");
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=AF,Number=A,Type=Float,Description=Alternate allele frequency>");
        if (options.keep_original_ac) {
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AC_Orig,Number=A,Type=Integer,Description=Original AC>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AF_Orig,Number=A,Type=Float,Description=Original AF>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AN_Orig,Number=1,Type=Integer,Description=Original AN>");
        }
        if (options.keep_original_dp && bcf_hdr_id2int(output_header, BCF_DT_ID, "DP_Orig") < 0)
            bcf_hdr_append(output_header, "##INFO=<ID=DP_Orig,Number=1,Type=Integer,Description=Original DP>");
        const auto indices = sample_indices(input_header, options, sample_spec);
        const bool subset_samples = !options.include_samples.empty() || !options.exclude_samples.empty();
        if (subset_samples) {
            const char* selection = sample_spec.empty() ? nullptr : sample_spec.c_str();
            if (bcf_hdr_set_samples(output_header, selection, 0) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset sample header");
            if (bcf_hdr_sync(output_header) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync subset sample header");
        }
        bcf_hdr_append(output_header, "##fastgatk_select_variants_status=prototype-filter-subset-jexl");
        if (bcf_hdr_sync(output_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync output header");
        writer_header = bcf_hdr_dup(output_header);
        if (!writer_header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate SelectVariants writer header");
        if (options.sites_only_vcf_output) {
            if (bcf_hdr_set_samples(writer_header, nullptr, 0) != 0 ||
                bcf_hdr_sync(writer_header) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot configure sites-only VCF header");
        }
        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open SelectVariants output: " + options.output);
        if (bcf_hdr_write(output, writer_header) != 0) throw std::runtime_error("cannot write SelectVariants header");
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
        while (fastgatk::io::read_variant_record(input, input_header, record,
                                               options.input) == 0) {
            ++input_records;
            bcf_unpack(record, BCF_UN_ALL);
            canonicalize_site_filters(input_header, record);
            if (!in_regions(input_header, record, regions)) {
                ++interval_skipped;
                bcf_clear(record);
                continue;
            }
            if (comparison_mode) {
                const auto key = variant_key(input_header, record);
                bool present = comparison_data.site_keys.contains(key);
                if (options.compare_genotypes && present) {
                    bool genotype_match = false;
                    for (int sample = 0; sample < record->n_sample; ++sample) {
                        const auto* name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, sample);
                        if (name == nullptr || !comparison_data.sample_names.contains(name)) continue;
                        const auto genotype = genotype_string(input_header, record, sample);
                        if (!genotype.empty() && comparison_data.genotype_keys.contains(
                                key + "\t" + name + "\t" + genotype)) {
                            genotype_match = true;
                            break;
                        }
                    }
                    present = genotype_match;
                }
                const bool keep = !options.concordance.empty() ? present : !present;
                if (!keep) {
                    ++comparison_skipped;
                    bcf_clear(record);
                    continue;
                }
            }
            const auto type = variant_type(record);
            if ((!options.include_types.empty() && !contains(options.include_types, type)) ||
                contains(options.exclude_types, type)) {
                ++type_skipped;
                bcf_clear(record);
                continue;
            }
            if (options.exclude_non_variants && is_non_variant(input_header, record)) {
                ++non_variant_skipped;
                bcf_clear(record);
                continue;
            }
            if (options.min_indel_size > 0 || options.max_indel_size < std::numeric_limits<int>::max()) {
                int minimum_indel = 0;
                int maximum_indel = 0;
                if (indel_size_range(record, minimum_indel, maximum_indel) &&
                    (minimum_indel < options.min_indel_size || maximum_indel > options.max_indel_size)) {
                    ++indel_size_skipped;
                    bcf_clear(record);
                    continue;
                }
            }
            if (options.exclude_filtered && filtered(input_header, record)) {
                ++filtered_records;
                bcf_clear(record);
                continue;
            }
            // GATK combines repeated -select and -select-genotype arguments
            // with logical OR; &&/|| inside one argument are evaluated by the
            // recursive expression parser.  Genotype expressions match when
            // at least one sample satisfies the predicate.
            bool selected_by_expression = options.select_expressions.empty() &&
                options.select_genotype_expressions.empty();
            for (const auto& expression : options.select_expressions)
                selected_by_expression = selected_by_expression ||
                    select_expression(expression, input_header, record);
            for (const auto& expression : options.select_genotype_expressions) {
                for (int sample = 0; sample < record->n_sample; ++sample) {
                    if (select_genotype_expression(expression, input_header, record, sample)) {
                        selected_by_expression = true;
                        break;
                    }
                }
                if (selected_by_expression) break;
            }
            if (!selected_by_expression) {
                ++expression_skipped;
                bcf_clear(record);
                continue;
            }
            const int concrete_alts = std::count_if(record->d.allele + 1, record->d.allele + record->n_allele,
                [](const char* allele) { return allele != nullptr && allele[0] != '<'; });
            if (options.restrict_alleles == "BI" && concrete_alts > 1) { ++type_skipped; bcf_clear(record); continue; }
            if (options.restrict_alleles == "MULTI" && concrete_alts <= 1) { ++type_skipped; bcf_clear(record); continue; }
            const auto original_annotations = capture_original_annotations(
                input_header, record, options.keep_original_ac, options.keep_original_dp);
            std::vector<int> retained_alt_indices;
            if (record->n_allele > 1)
                for (int allele = 1; allele < record->n_allele; ++allele)
                    retained_alt_indices.push_back(allele);
            if (subset_samples) {
                std::vector<int> map = indices;
                if (bcf_subset(input_header, record, static_cast<int>(map.size()), map.data()) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset variant FORMAT fields");
                if (bcf_translate(output_header, input_header, record) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot translate variant header");
            }
            std::uint64_t record_filtered_gt_nocall = 0;
            if (options.set_filtered_gt_to_nocall)
                record_filtered_gt_nocall = set_filtered_genotypes_to_nocall(
                    output_header, record);
            clear_all_pass_genotype_filters(output_header, record);
            filtered_gt_nocall += record_filtered_gt_nocall;
            if (options.remove_unused_alternates && record->n_allele > 1 && record->n_sample > 0) {
                int32_t* genotypes = nullptr;
                int genotype_count = 0;
                std::vector<int> kept{0};
                if (bcf_get_genotypes(output_header, record, &genotypes, &genotype_count) > 0) {
                    for (int index = 0; index < genotype_count; ++index) {
                        if (bcf_gt_is_missing(genotypes[index])) continue;
                        const auto allele = bcf_gt_allele(genotypes[index]);
                        if (allele > 0 && std::find(kept.begin(), kept.end(), allele) == kept.end()) kept.push_back(allele);
                    }
                }
                free(genotypes);
                if (kept.size() < static_cast<std::size_t>(record->n_allele)) {
                    const int old_alleles = record->n_allele;
                    // A record whose samples use no ALT keeps the reference
                    // allele alone; GATK emits it with ALT='.'.  The
                    // --exclude-non-variants test of SelectVariants.java:708-715
                    // runs *after* the trimming (and after
                    // --set-filtered-gt-to-nocall), so a record trimmed down to
                    // the reference allele is non-variant by definition and is
                    // dropped here even though the pre-trim check above, which
                    // sees the original record, could not tell.
                    if (kept.size() == 1) {
                        if (options.exclude_non_variants) {
                            ++non_variant_skipped;
                            bcf_clear(record);
                            continue;
                        }
                        retained_alt_indices.clear();
                        trim_record_to_reference_only(output_header, record, old_alleles,
                                                      genotype_kernel_telemetry);
                    } else {
                        std::string alleles;
                        retained_alt_indices.clear();
                        for (std::size_t index = 0; index < kept.size(); ++index) {
                            if (index != 0) alleles.push_back(',');
                            alleles += record->d.allele[kept[index]];
                            if (index != 0) retained_alt_indices.push_back(kept[index]);
                        }
                        if (bcf_update_alleles_str(output_header, record, alleles.c_str()) != 0)
                            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove unused alternate alleles");
                        compact_fields(output_header, record, kept, old_alleles,
                                       genotype_kernel_telemetry);
                    }
                }
            }
            // GATK refreshes AC/AN/AF only when the selected sample/allele
            // space changes.  A pure site-selection pass preserves existing
            // INFO annotations byte-for-byte; this matters for annotation
            // pipelines that intentionally carry precomputed values.
            if (subset_samples || options.remove_unused_alternates || record_filtered_gt_nocall > 0)
                recompute_site_annotations(output_header, record);
            if (options.keep_original_ac || options.keep_original_dp)
                restore_original_annotations(output_header, record, options, original_annotations,
                                             retained_alt_indices);
            if (options.sites_only_vcf_output && bcf_subset_format(writer_header, record) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sites-only SelectVariants record");
            if (bcf_write(output, writer_header, record) != 0)
                throw std::runtime_error("cannot write SelectVariants record");
            ++output_records;
            bcf_clear(record);
        }
        bcf_destroy(record); record = nullptr;
        bcf_close(output); output = nullptr;
        bcf_close(input); input = nullptr;
        std::string index;
        if (options.create_index) {
            if (suffix(options.output, ".gz")) {
                index = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build SelectVariants index");
            } else {
                index = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index);
            }
        }
        if (!file_complete(options.output) || (!index.empty() && !file_complete(index)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: SelectVariants output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write SelectVariants manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"SelectVariants\",\"implementation\":\"fastgatk-select-variants\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"type_filter\":true,\"sample_subset\":"
                 << ((!options.include_samples.empty() || !options.exclude_samples.empty()) ? "true" : "false")
                 << ",\"exclude_filtered\":" << (options.exclude_filtered ? "true" : "false")
                 << ",\"remove_unused_alternates\":" << (options.remove_unused_alternates ? "true" : "false")
                 << ",\"exclude_non_variants\":" << (options.exclude_non_variants ? "true" : "false")
                 << ",\"set_filtered_gt_to_nocall\":" << (options.set_filtered_gt_to_nocall ? "true" : "false")
                 << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"indel_size_filter\":" << ((options.min_indel_size > 0 || options.max_indel_size < std::numeric_limits<int>::max()) ? "true" : "false")
                 << ",\"keep_original_ac\":" << (options.keep_original_ac ? "true" : "false")
                 << ",\"keep_original_dp\":" << (options.keep_original_dp ? "true" : "false")
                 << ",\"site_annotation_recompute\":true"
                 << ",\"jexl_expression_subset\":" << (!options.select_expressions.empty() ? "true" : "false")
                 << ",\"jexl_info_vector_indexing\":true,\"jexl_boolean_comparisons\":true,\"jexl_regex_operators\":true,\"jexl_arithmetic\":true"
                 << ",\"genotype_expression_subset\":" << (!options.select_genotype_expressions.empty() ? "true" : "false")
                 << ",\"concordance_discordance_subset\":" << (comparison_mode ? "true" : "false")
                 << ",\"sample_level_concordance\":" << (options.compare_genotypes ? "true" : "false")
                 << ",\"interval_subset\":" << (!regions.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"vcf_index\":" << (index.empty() ? "false" : "true")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":true}"
                 << (index.empty() ? "" : ",{\"path\":\"" + json_escape(index) + "\",\"kind\":\"vcf-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << input_records << ",\"output_records\":" << output_records
                 << ",\"filtered_records\":" << filtered_records << ",\"type_skipped\":" << type_skipped
                 << ",\"non_variant_skipped\":" << non_variant_skipped
                 << ",\"indel_size_skipped\":" << indel_size_skipped
                 << ",\"expression_skipped\":" << expression_skipped
                 << ",\"comparison_skipped\":" << comparison_skipped
                 << ",\"filtered_gt_nocall\":" << filtered_gt_nocall
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                 << ",\"pl_remap_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.pl_remap_prepare_seconds
                 << ",\"pl_remap_kernel_seconds\":"
                 << genotype_kernel_telemetry.pl_remap_execute_seconds
                 << ",\"allele_field_remap_kernel_calls\":"
                 << genotype_kernel_telemetry.allele_field_remap_calls
                 << ",\"allele_field_remap_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
                 << ",\"allele_field_remap_kernel_seconds\":"
                 << genotype_kernel_telemetry.allele_field_remap_execute_seconds
                 << ",\"allele_field_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space) << "\""
                 << ",\"gt_gq_kernel_calls\":" << genotype_kernel_telemetry.gt_gq_calls
                 << ",\"gt_gq_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.gt_gq_prepare_seconds
                 << ",\"gt_gq_kernel_seconds\":"
                 << genotype_kernel_telemetry.gt_gq_execute_seconds
                 << ",\"genotype_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.execution_space) << "\""
                 << ",\"sample_count\":" << writer_header->n[BCF_DT_SAMPLE] << "}}\n";
        manifest.flush();
        if (!manifest) throw std::runtime_error("cannot finalize SelectVariants manifest");
        manifest.close();
        fastgatk::runtime::require_complete_output({
            std::filesystem::path(options.output), std::filesystem::path(index),
            std::filesystem::path(manifest_path), !index.empty(), true});
        std::cout << "{\"tool\":\"SelectVariants\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"output_records\":" << output_records
                  << ",\"expression_skipped\":" << expression_skipped
                  << ",\"comparison_skipped\":" << comparison_skipped
                  << ",\"filtered_gt_nocall\":" << filtered_gt_nocall
                  << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"interval_skipped\":" << interval_skipped << "}\n";
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(writer_header);
        return 0;
    } catch (...) {
        bcf_destroy(record);
        if (output) bcf_close(output);
        if (input) bcf_close(input);
        if (output_header) bcf_hdr_destroy(output_header);
        if (writer_header) bcf_hdr_destroy(writer_header);
        if (input_header) bcf_hdr_destroy(input_header);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for SelectVariants");
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
        std::cerr << "fastgatk-select-variants: " << error.what() << '\n';
        return 2;
    }
}
