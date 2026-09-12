#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
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
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> excluded_regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::string interval_merging_rule = "ALL";
    int interval_padding = 0;
    int exclusion_padding = 0;
    std::string output;
    std::string manifest;
    bool dont_trim = false;
    bool split_multi_allelics = false;
    bool keep_original_ac = false;
    bool sites_only_vcf_output = false;
    bool create_index = true;
    int max_indel_length = 200;
    int max_leading_bases = 1000;
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
    int max_ploidy = 0;
    std::uint64_t symbolic_split_records = 0;
    std::uint64_t symbolic_info_fields_cleared = 0;
    std::uint64_t symbolic_format_fields_cleared = 0;
    std::uint64_t oversized_indel_records = 0;
    std::uint64_t left_shift_bases = 0;
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
    const bool regular = std::filesystem::is_regular_file(path, error);
    if (error || !regular) return false;
    const auto size = std::filesystem::file_size(path, error);
    return !error && size > 0;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-left-align-trim (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE             input VCF/BCF\n"
                         "  -R, --reference FILE           reference FASTA\n"
                         "  -L, --intervals REGION         interval subset (repeatable)\n"
                         "  -XL, --exclude-intervals REGION  exclude interval (repeatable)\n"
                         "      --interval-set-rule RULE  UNION (default) or INTERSECTION\n"
                         "  -imr, --interval-merging-rule RULE ALL (default) or OVERLAPPING_ONLY\n"
                         "  -ip, --interval-padding INT    include padding (default 0)\n"
                         "  -ixp, --interval-exclusion-padding INT exclude padding (default 0)\n"
                         "  -O, --output FILE              normalized VCF/BCF\n"
                         "      --dont-trim-alleles        left-align without common prefix/suffix trimming\n"
                         "      --sites-only-vcf-output BOOL  omit FORMAT/sample columns (default false)\n"
                         "      --max-indel-length INT     maximum indel size to realign (default 200)\n"
                         "      --max-leading-bases INT    maximum reference bases to look back (default 1000)\n"
                         "      --output-manifest FILE    OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-V" || is_option(argument, "--variant"))
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region")) {
            const char* interval_option = argument.rfind("--intervals", 0) == 0 ? "--intervals" :
                (argument.rfind("--interval", 0) == 0 ? "--interval" :
                 (argument.rfind("--region", 0) == 0 ? "--region" : "--intervals"));
            options.regions.push_back(require_value(index, argc, argv, argument, interval_option, "-L"));
        }
        else if (argument == "-XL" || is_option(argument, "--exclude-intervals") ||
                 is_option(argument, "--exclude-interval") || is_option(argument, "--exclude-region")) {
            const char* interval_option = argument == "-XL" ? "--exclude-intervals" :
                (argument.rfind("--exclude-region", 0) == 0 ? "--exclude-region" :
                 (argument.rfind("--exclude-interval", 0) == 0 ? "--exclude-interval" : "--exclude-intervals"));
            options.excluded_regions.push_back(require_value(index, argc, argv, argument, interval_option, "-XL"));
        }
        else if (argument == "-ip" || is_option(argument, "--interval-padding")) {
            options.interval_padding = std::stoi(require_value(index, argc, argv, argument,
                                                                "--interval-padding", "-ip"));
            if (options.interval_padding < 0)
                throw std::invalid_argument("--interval-padding must be non-negative");
        }
        else if (argument == "-ixp" || is_option(argument, "--interval-exclusion-padding")) {
            options.exclusion_padding = std::stoi(require_value(index, argc, argv, argument,
                                                                 "--interval-exclusion-padding", "-ixp"));
            if (options.exclusion_padding < 0)
                throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
        }
        else if (argument == "-isr" || is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-imr" || is_option(argument, "--interval-merging-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-merging-rule", "-imr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value != "ALL" && value != "OVERLAPPING_ONLY")
                throw std::invalid_argument("invalid --interval-merging-rule: " + value);
            options.interval_merging_rule = std::move(value);
        }
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--dont-trim-alleles" || argument == "--no-trim" || argument == "-no-trim") {
            const char* name = argument == "-no-trim" ? "-no-trim" :
                (argument == "--no-trim" ? "--no-trim" : "--dont-trim-alleles");
            options.dont_trim = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        }
        else if (is_option(argument, "--max-indel-length")) {
            options.max_indel_length = std::stoi(require_value(index, argc, argv, argument,
                                                                "--max-indel-length"));
            if (options.max_indel_length < 0)
                throw std::invalid_argument("--max-indel-length must be non-negative");
        } else if (is_option(argument, "--max-leading-bases")) {
            options.max_leading_bases = std::stoi(require_value(index, argc, argv, argument,
                                                                "--max-leading-bases"));
            if (options.max_leading_bases < 0)
                throw std::invalid_argument("--max-leading-bases must be non-negative");
        }
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0 ||
                 argument == "-OVI" || argument.rfind("-OVI=", 0) == 0) {
            const char* name = argument.rfind("-OVI", 0) == 0
                ? "-OVI" : "--create-output-variant-index";
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (argument == "--split-multi-allelics" ||
                   argument.rfind("--split-multi-allelics=", 0) == 0) {
            options.split_multi_allelics = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--split-multi-allelics");
        } else if (argument == "--keep-original-ac" ||
                   argument.rfind("--keep-original-ac=", 0) == 0) {
            options.keep_original_ac = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--keep-original-ac");
        } else if (argument == "--sites-only-vcf-output" ||
                   argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--sites-only-vcf-output");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Accepted compatibility switch.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
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

void normalize_left_align_regions(std::vector<Region>& regions) {
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

void apply_region_padding(std::vector<Region>& regions, const faidx_t* reference,
                          const bcf_hdr_t* header, int padding) {
    if (padding == 0) return;
    for (auto& region : regions) {
        region.begin = std::max(0, region.begin - padding);
        const char* contig = region.rid >= 0 ? bcf_hdr_id2name(header, region.rid) : nullptr;
        const auto length = contig == nullptr ? 0 : faidx_seq_len64(reference, contig);
        const auto limit = length > 0 ? static_cast<std::int64_t>(length)
                                      : static_cast<std::int64_t>(std::numeric_limits<int>::max());
        const auto padded_end = static_cast<std::int64_t>(region.end) + padding;
        region.end = static_cast<int>(std::min(limit, padded_end));
    }
}

void append_left_align_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_left_align_regions(incoming);
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
        end = std::max(end, values[0]);
    free(values);
    return end;
}

bool in_regions(const bcf_hdr_t* header, bcf1_t* record,
                const std::vector<Region>& regions) {
    if (regions.empty()) return true;
    const auto end = record_end_exclusive(header, record);
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        return record->rid == region.rid && end > region.begin && record->pos < region.end;
    });
}

bool overlaps_regions(const bcf_hdr_t* header, bcf1_t* record,
                      const std::vector<Region>& regions) {
    if (regions.empty()) return false;
    const auto end = record_end_exclusive(header, record);
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        return record->rid == region.rid && end > region.begin && record->pos < region.end;
    });
}

std::string reference_base(const faidx_t* reference, const std::string& contig, int one_based_position) {
    if (one_based_position < 1) return {};
    const std::string region = contig + ":" + std::to_string(one_based_position) + "-" +
                               std::to_string(one_based_position);
    hts_pos_t length = 0;
    char* bases = fai_fetch64(reference, region.c_str(), &length);
    if (!bases || length != 1) {
        free(bases);
        return {};
    }
    const std::string result(1, bases[0]);
    free(bases);
    return result;
}

bool symbolic(const std::string& allele) {
    return allele.empty() || allele == "." || allele == "*" || allele.front() == '<' ||
           allele.find('[') != std::string::npos || allele.find(']') != std::string::npos ||
           (allele.size() > 1 && (allele.front() == '.' || allele.back() == '.' || allele.back() == '>'));
}

bool trim_common(std::vector<std::string>& alleles, int& one_based_position) {
    if (alleles.size() < 2 || symbolic(alleles.front())) return false;
    bool changed = false;
    while (true) {
        if (std::any_of(alleles.begin(), alleles.end(), [](const std::string& allele) { return allele.size() <= 1; })) break;
        const char last = alleles.front().back();
        if (!std::all_of(alleles.begin(), alleles.end(), [last](const std::string& allele) { return allele.back() == last; })) break;
        for (auto& allele : alleles) allele.pop_back();
        changed = true;
    }
    while (true) {
        if (std::any_of(alleles.begin(), alleles.end(), [](const std::string& allele) { return allele.size() <= 1; })) break;
        const char first = alleles.front().front();
        if (!std::all_of(alleles.begin(), alleles.end(), [first](const std::string& allele) { return allele.front() == first; })) break;
        for (auto& allele : alleles) allele.erase(allele.begin());
        ++one_based_position;
        changed = true;
    }
    return changed;
}

// splitVariantContextToBiallelics always right-trims each newly created pair,
// even when LeftAlignAndTrimVariants was called with --dont-trim-alleles.  The
// full left-align pass below performs the optional forward trim separately.
// Keep this operation concrete-only: symbolic alleles are not compared by
// HTSJDK's trimAlleles helper, but remain in the pair unchanged.
bool trim_common_suffix(std::vector<std::string>& alleles) {
    if (alleles.size() < 2 || symbolic(alleles.front())) return false;
    std::vector<std::size_t> comparable;
    comparable.reserve(alleles.size());
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        if (symbolic(alleles[index])) continue;
        // A one-base concrete allele is the VCF anchor and disables the
        // helper's common trimming for this pair.
        if (alleles[index].size() <= 1) return false;
        comparable.push_back(index);
    }
    if (comparable.size() < 2) return false;
    bool changed = false;
    while (true) {
        const auto first = comparable.front();
        if (alleles[first].size() <= 1) break;
        const char last = alleles[first].back();
        if (!std::all_of(comparable.begin(), comparable.end(), [&](const std::size_t index) {
                return alleles[index].size() > 1 && alleles[index].back() == last;
            })) break;
        for (const auto index : comparable) alleles[index].pop_back();
        changed = true;
    }
    return changed;
}

bool right_trim_split_record(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2 || record->d.allele == nullptr) return false;
    std::vector<std::string> alleles;
    alleles.reserve(record->n_allele);
    for (int index = 0; index < record->n_allele; ++index) {
        if (!record->d.allele[index]) return false;
        alleles.emplace_back(record->d.allele[index]);
    }
    if (!trim_common_suffix(alleles)) return false;
    std::string allele_string;
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        if (index != 0) allele_string.push_back(',');
        allele_string += alleles[index];
    }
    if (bcf_update_alleles_str(header, record, allele_string.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot right-trim split alleles");
    return true;
}

bool normalize_record(const faidx_t* reference, const bcf_hdr_t* header, bcf1_t* record,
                      bool dont_trim, int max_indel_length, int max_leading_bases,
                      GenotypeKernelTelemetry& telemetry) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2 || record->rid < 0) return false;
    const char* contig_name = bcf_hdr_id2name(header, record->rid);
    if (!contig_name) return false;
    std::vector<std::string> alleles;
    alleles.reserve(record->n_allele);
    for (int index = 0; index < record->n_allele; ++index) {
        if (!record->d.allele[index]) return false;
        alleles.emplace_back(record->d.allele[index]);
    }
    // GATK's workhorse is only called for biallelic, concrete indels.  In
    // particular a multiallelic record without --split-multi-allelics is
    // copied byte-for-byte, and a non-indel is not common-trimmed.
    if (record->n_allele != 2 || symbolic(alleles.front()) ||
        symbolic(alleles[1]) || alleles.front().size() == alleles[1].size()) return false;
    if (max_leading_bases <= 0) return false;
    const auto oversized = std::any_of(alleles.begin() + 1, alleles.end(), [&](const std::string& allele) {
        const auto difference = alleles.front().size() > allele.size()
            ? alleles.front().size() - allele.size() : allele.size() - alleles.front().size();
        return difference > static_cast<std::size_t>(max_indel_length);
    });
    if (oversized) {
        ++telemetry.oversized_indel_records;
        return false;
    }
    int one_based_position = record->pos + 1;
    bool changed = false;
    if (!dont_trim) changed = trim_common(alleles, one_based_position) || changed;

    const bool indel = std::any_of(alleles.begin() + 1, alleles.end(), [&](const std::string& allele) {
        return allele.size() != alleles.front().size();
    });
    // Rotate all alleles together through a repeated reference base. This is
    // the deterministic left-shift rule for normalized primitive indels; a
    // multiallelic site is shifted only when every ALT can make the same move.
    if (indel) {
        int shifted_bases = 0;
        while (one_based_position > 1 && shifted_bases < max_leading_bases) {
            const auto previous = reference_base(reference, contig_name, one_based_position - 1);
            if (previous.empty() || std::any_of(alleles.begin(), alleles.end(), [&](const std::string& allele) {
                    return allele.empty() || allele.back() != previous.front();
                })) break;
            for (auto& allele : alleles) allele = previous + allele.substr(0, allele.size() - 1);
            --one_based_position;
            ++shifted_bases;
            ++telemetry.left_shift_bases;
            changed = true;
            if (!dont_trim) changed = trim_common(alleles, one_based_position) || changed;
        }
    }
    if (!changed) return false;
    std::string allele_string;
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        if (index != 0) allele_string.push_back(',');
        allele_string += alleles[index];
    }
    if (bcf_update_alleles_str(header, record, allele_string.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update normalized alleles");
    record->pos = one_based_position - 1;
    return true;
}

std::size_t genotype_count(int allele_count, int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t result = 1;
    for (int position = 1; position <= ploidy; ++position) {
        const auto numerator = static_cast<std::size_t>(allele_count + position - 1);
        if (numerator != 0 && result > std::numeric_limits<std::size_t>::max() / numerator)
            return 0;
        result *= numerator;
        result /= static_cast<std::size_t>(position);
    }
    return result;
}

int infer_ploidy_from_width(int allele_count, std::size_t width) {
    for (int ploidy = 1; ploidy <= 32; ++ploidy)
        if (genotype_count(allele_count, ploidy) == width) return ploidy;
    return 0;
}

int gq_from_normalized_pl(const std::int32_t* values, int width) {
    int best_index = -1;
    std::int32_t best = std::numeric_limits<std::int32_t>::max();
    std::int32_t runner_up = std::numeric_limits<std::int32_t>::max();
    for (int index = 0; index < width; ++index) {
        const auto value = values[index];
        if (value < 0) continue;
        if (value < best) {
            runner_up = best;
            best = value;
            best_index = index;
        } else if (value < runner_up) {
            runner_up = value;
        }
    }
    if (best_index < 0 || runner_up == std::numeric_limits<std::int32_t>::max()) return -1;
    // GenotypeLikelihoods.getGQLog10FromLikelihoods() returns the simple
    // best-vs-runner-up log-space difference whenever the winner is valid.
    // Since PL is -10*log10(L), this is the runner-up PL after rescaling;
    // Math.round is the conversion used by GenotypeBuilder.log10PError.
    return std::min(99, std::max(0, static_cast<int>(runner_up - best)));
}

void compact_split_fields(const bcf_hdr_t* header, bcf1_t* record,
                          int selected_alt, int original_alleles,
                          GenotypeKernelTelemetry& telemetry) {
    const int samples = record->n_sample;
    if (samples <= 0) return;
    int32_t* source_gq_values = nullptr;
    int source_gq_count = 0;
    const auto source_gq_length = bcf_get_format_int32(header, record, "GQ",
                                                       &source_gq_values, &source_gq_count);
    const bool source_has_gq = source_gq_length > 0 && source_gq_count >= samples &&
                               source_gq_count % samples == 0;
    int32_t* source_pl_values = nullptr;
    int source_pl_count = 0;
    const auto source_pl_length = bcf_get_format_int32(header, record, "PL",
                                                        &source_pl_values, &source_pl_count);
    const auto source_pl_width = (source_pl_length > 0 && source_pl_count >= samples &&
                                  source_pl_count % samples == 0)
        ? static_cast<std::size_t>(source_pl_count / samples) : 0;
    int ploidy = 0;
    int32_t* genotypes = nullptr;
    int genotype_value_count = 0;
    const auto genotype_length = bcf_get_genotypes(header, record, &genotypes, &genotype_value_count);
    if (genotype_length > 0 && genotype_value_count >= samples && genotype_value_count % samples == 0) {
        ploidy = genotype_value_count / samples;
        for (int index = 0; index < genotype_value_count; ++index) {
            if (bcf_gt_is_missing(genotypes[index]) || genotypes[index] == bcf_int32_vector_end) {
                genotypes[index] = bcf_gt_missing;
                continue;
            }
            const int sample = index / ploidy;
            const bool no_call_zero_gq = source_has_gq && source_gq_values[sample] == 0 &&
                (source_pl_width == 0 || source_pl_values[sample * source_pl_width] == 0);
            if (no_call_zero_gq) {
                genotypes[index] = bcf_gt_missing;
                continue;
            }
            const auto allele = bcf_gt_allele(genotypes[index]);
            const bool phased = bcf_gt_is_phased(genotypes[index]);
            if (allele == 0) genotypes[index] = phased ? bcf_gt_phased(0) : bcf_gt_unphased(0);
            else if (allele == selected_alt)
                genotypes[index] = phased ? bcf_gt_phased(1) : bcf_gt_unphased(1);
            // GATK's allele-subset path projects removed alternate alleles
            // onto REF (rather than introducing a no-call), preserving
            // ploidy and AN for split records.
            else genotypes[index] = phased ? bcf_gt_phased(0) : bcf_gt_unphased(0);
        }
        if (bcf_update_genotypes(header, record, genotypes, genotype_value_count) != 0) {
            free(genotypes);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact split GT");
        }
    }
    free(genotypes);

    // Number=G is defined by the sample ploidy, not always by diploid
    // triangular width.  If GT is missing, infer the unique bounded ploidy
    // from the source PL width below.
    if (ploidy <= 0 && source_pl_width > 0)
        ploidy = infer_ploidy_from_width(original_alleles, source_pl_width);
    if (ploidy > 0) {
        telemetry.max_ploidy = std::max(telemetry.max_ploidy, ploidy);
    }

    if (bcf_hdr_id2int(header, BCF_DT_ID, "AD") >= 0) {
        int32_t* source_ad_values = nullptr;
        int source_ad_count = 0;
        const auto source_ad_length = bcf_get_format_int32(header, record, "AD",
                                                            &source_ad_values, &source_ad_count);
        if (source_ad_length > 0 && source_ad_count == samples * original_alleles) {
            const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
                std::vector<std::int32_t>(source_ad_values, source_ad_values + source_ad_count),
                static_cast<std::size_t>(samples), original_alleles, 2,
                {0, selected_alt});
            telemetry.allele_field_remap_calls++;
            telemetry.allele_field_remap_prepare_seconds += remapped.prepare_seconds;
            telemetry.allele_field_remap_execute_seconds += remapped.seconds;
            telemetry.allele_field_remap_execution_space = remapped.execution_space;
            if (bcf_update_format_int32(header, record, "AD", remapped.values.data(),
                                        static_cast<int>(remapped.values.size())) != 0) {
                free(source_ad_values);
                free(source_pl_values);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact split AD");
            }
        }
        free(source_ad_values);
    }
    if (ploidy > 0 && source_pl_width > 0 &&
        source_pl_width == genotype_count(original_alleles, ploidy)) {
        auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
            std::vector<std::int32_t>(source_pl_values, source_pl_values + source_pl_count),
            static_cast<std::size_t>(samples), original_alleles, 2, ploidy,
            {0, selected_alt});
        telemetry.pl_remap_calls++;
        telemetry.pl_remap_prepare_seconds += remapped.prepare_seconds;
        telemetry.pl_remap_execute_seconds += remapped.seconds;
        telemetry.pl_remap_execution_space = remapped.execution_space;
        // HTSJDK normalizes each sample's PL vector by subtracting its
        // minimum finite likelihood after allele projection.
        const int target_width = static_cast<int>(genotype_count(2, ploidy));
        for (int sample = 0; sample < samples; ++sample) {
            auto begin = remapped.pl.begin() + static_cast<std::ptrdiff_t>(sample * target_width);
            auto end = begin + target_width;
            auto minimum = std::numeric_limits<std::int32_t>::max();
            for (auto value = begin; value != end; ++value)
                if (*value >= 0) minimum = std::min(minimum, *value);
            if (minimum != std::numeric_limits<std::int32_t>::max()) {
                for (auto value = begin; value != end; ++value)
                    if (*value >= 0) *value -= minimum;
            }
        }
        if (bcf_update_format_int32(header, record, "PL", remapped.pl.data(),
                                    static_cast<int>(remapped.pl.size())) != 0) {
            free(source_pl_values);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact split PL");
        }
        if (source_has_gq && bcf_hdr_id2int(header, BCF_DT_ID, "GQ") >= 0 && target_width >= 2) {
            std::vector<int32_t> gq(static_cast<std::size_t>(samples), bcf_int32_missing);
            for (int sample = 0; sample < samples; ++sample) {
                const auto source_gq = source_gq_values[sample];
                if (source_gq == bcf_int32_missing || source_gq == bcf_int32_vector_end) continue;
                gq[static_cast<std::size_t>(sample)] = gq_from_normalized_pl(
                    remapped.pl.data() + static_cast<std::ptrdiff_t>(sample * target_width),
                    target_width);
            }
            if (bcf_update_format_int32(header, record, "GQ", gq.data(), samples) != 0) {
                free(source_pl_values);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot recompute split GQ");
            }
        }
    }
    free(source_gq_values);
    free(source_pl_values);
}

void recompute_split_annotations(const bcf_hdr_t* header, bcf1_t* record) {
    // The splitter only retains AC/AF/AN as candidate site annotations.  A
    // record with no called alleles loses them during GATK's chromosome-count
    // recomputation rather than carrying stale source values forward.
    for (const auto* name : {"AC", "AF", "AN"}) {
        if (bcf_hdr_id2int(header, BCF_DT_ID, name) >= 0 &&
            bcf_update_info(header, record, name, nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot clear split ") + name);
    }
    if (record->n_sample <= 0) return;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || genotype_count < record->n_sample || genotype_count % record->n_sample != 0) {
        free(genotypes);
        return;
    }
    const auto ploidy = genotype_count / record->n_sample;
    int32_t an = 0;
    int32_t ac = 0;
    for (int index = 0; index < genotype_count; ++index) {
        if (bcf_gt_is_missing(genotypes[index]) || genotypes[index] == bcf_int32_vector_end) continue;
        ++an;
        if (bcf_gt_allele(genotypes[index]) == 1) ++ac;
    }
    free(genotypes);
    if (an == 0) return;
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AN") >= 0 &&
        bcf_update_info_int32(header, record, "AN", &an, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split AN");
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AC") >= 0 &&
        bcf_update_info_int32(header, record, "AC", &ac, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split AC");
    if (bcf_hdr_id2int(header, BCF_DT_ID, "AF") >= 0) {
        const float af = an == 0 ? 0.0F : static_cast<float>(ac) / static_cast<float>(an);
        if (bcf_update_info_float(header, record, "AF", &af, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update split AF");
    }
    (void)ploidy;
}

void preserve_original_counts(const bcf_hdr_t* header, const bcf1_t* source,
                              bcf1_t* destination, int selected_alt) {
    int32_t* ac = nullptr;
    int ac_count = 0;
    const auto ac_length = bcf_get_info_int32(header, const_cast<bcf1_t*>(source), "AC", &ac, &ac_count);
    if (ac_length > 0 && selected_alt > 0 && selected_alt <= ac_count) {
        const int32_t value = ac[selected_alt - 1];
        if (bcf_update_info_int32(header, destination, "AC_Orig", &value, 1) != 0) {
            free(ac);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot preserve AC_Orig");
        }
    }
    free(ac);
    float* af = nullptr;
    int af_count = 0;
    const auto af_length = bcf_get_info_float(header, const_cast<bcf1_t*>(source), "AF", &af, &af_count);
    if (af_length > 0 && selected_alt > 0 && selected_alt <= af_count) {
        const float value = af[selected_alt - 1];
        if (bcf_update_info_float(header, destination, "AF_Orig", &value, 1) != 0) {
            free(af);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot preserve AF_Orig");
        }
    }
    free(af);
    int32_t* an = nullptr;
    int an_count = 0;
    const auto an_length = bcf_get_info_int32(header, const_cast<bcf1_t*>(source), "AN", &an, &an_count);
    if (an_length > 0 && an_count > 0) {
        if (bcf_update_info_int32(header, destination, "AN_Orig", an, 1) != 0) {
            free(an);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot preserve AN_Orig");
        }
    }
    free(an);
}

bool record_has_symbolic_alt(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele <= 2 || record->d.allele == nullptr) return false;
    for (int index = 1; index < record->n_allele; ++index) {
        if (record->d.allele[index] != nullptr && symbolic(record->d.allele[index])) return true;
    }
    (void)header;
    return false;
}

// GATK switches the whole split operation to SET_TO_NO_CALL_NO_ANNOTATIONS
// when any sample is a diploid heterozygous non-reference call (1/2).  The
// switch is independent of whether an ALT is symbolic.  Typed DP survives
// this assignment; AD/PL/GQ and extended FORMAT attributes do not.
bool has_het_non_ref(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->n_sample <= 0) return false;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || genotype_count < record->n_sample ||
        genotype_count % record->n_sample != 0) {
        free(genotypes);
        return false;
    }
    const int ploidy = genotype_count / record->n_sample;
    bool result = false;
    // Genotype.isHetNonRef() is the diploid 1/2 case.  Polyploid calls are
    // not classified as HET_NON_REF by HTSJDK's GenotypeType.
    if (ploidy == 2) {
        for (int sample = 0; sample < record->n_sample; ++sample) {
            const auto first = genotypes[2 * sample];
            const auto second = genotypes[2 * sample + 1];
            if (first == bcf_int32_vector_end || second == bcf_int32_vector_end ||
                bcf_gt_is_missing(first) || bcf_gt_is_missing(second)) continue;
            const auto first_allele = bcf_gt_allele(first);
            const auto second_allele = bcf_gt_allele(second);
            if (first_allele > 0 && second_allele > 0 && first_allele != second_allele) {
                result = true;
                break;
            }
        }
    }
    free(genotypes);
    return result;
}

void reset_no_call_split_fields(const bcf_hdr_t* header, bcf1_t* record,
                                GenotypeKernelTelemetry& telemetry) {
    bcf_unpack(record, BCF_UN_FMT | BCF_UN_INFO);

    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto genotype_length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (genotype_length > 0 && genotype_count > 0) {
        for (int index = 0; index < genotype_count; ++index)
            genotypes[index] = bcf_gt_missing;
        if (bcf_update_genotypes(header, record, genotypes, genotype_count) != 0) {
            free(genotypes);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot reset symbolic split GT");
        }
    }
    free(genotypes);

    std::vector<std::string> format_names;
    format_names.reserve(record->n_fmt);
    for (int index = 0; index < record->n_fmt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.fmt[index].id);
        if (name != nullptr && std::strcmp(name, "GT") != 0) format_names.emplace_back(name);
    }
    for (const auto& name : format_names) {
        // GenotypeBuilder.noAttributes().noGQ().noAD().noPL() leaves DP in
        // place.  This matters for depth-only records produced by GATK.
        if (name == "DP") continue;
        if (bcf_update_format(header, record, name.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear no-call FORMAT/" + name);
        ++telemetry.symbolic_format_fields_cleared;
    }

    std::vector<std::string> info_names;
    info_names.reserve(record->n_info);
    for (int index = 0; index < record->n_info; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.info[index].key);
        if (name != nullptr) info_names.emplace_back(name);
    }
    for (const auto& name : info_names) {
        if (bcf_update_info(header, record, name.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear no-call INFO/" + name);
        ++telemetry.symbolic_info_fields_cleared;
    }
}

void clear_split_info_annotations(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_INFO);
    std::vector<std::string> info_names;
    info_names.reserve(record->n_info);
    for (int index = 0; index < record->n_info; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.info[index].key);
        if (name != nullptr && std::strcmp(name, "AC") != 0 &&
            std::strcmp(name, "AF") != 0 && std::strcmp(name, "AN") != 0)
            info_names.emplace_back(name);
    }
    for (const auto& name : info_names) {
        if (bcf_update_info(header, record, name.c_str(), nullptr, 0, BCF_HT_STR) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear split INFO/" + name);
    }
}

std::vector<bcf1_t*> split_record(const bcf_hdr_t* header, const bcf1_t* source,
                                  bool keep_original_ac,
                                  GenotypeKernelTelemetry& telemetry) {
    std::vector<bcf1_t*> result;
    auto* mutable_source = const_cast<bcf1_t*>(source);
    bcf_unpack(mutable_source, BCF_UN_STR | BCF_UN_FMT | BCF_UN_INFO);
    if (source->n_allele <= 1) {
        // GATK's splitVariantContextToBiallelics treats a non-variant as an
        // empty list, so --split-multi-allelics drops it from the output.
        return result;
    }
    if (source->n_allele == 2) {
        auto* copy = bcf_dup(const_cast<bcf1_t*>(source));
        if (!copy) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate variant record");
        result.push_back(copy);
        return result;
    }
    const int original_alleles = source->n_allele;
    const bool symbolic_alt = record_has_symbolic_alt(header, mutable_source);
    const bool no_annotations = has_het_non_ref(header, mutable_source);
    for (int selected = 1; selected < original_alleles; ++selected) {
        auto* copy = bcf_dup(const_cast<bcf1_t*>(source));
        if (!copy) {
            for (auto* record : result) bcf_destroy(record);
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate split variant record");
        }
        std::string alleles = std::string(source->d.allele[0]) + "," + source->d.allele[selected];
        if (bcf_update_alleles_str(header, copy, alleles.c_str()) != 0) {
            bcf_destroy(copy);
            for (auto* record : result) bcf_destroy(record);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot split multiallelic alleles");
        }
        // The split helper always right-trims, independently of the later
        // --dont-trim-alleles alignment option.  Left trimming is left to the
        // traversal because its window is bounded by the previous output.
        (void)right_trim_split_record(header, copy);
        if (no_annotations) {
            reset_no_call_split_fields(header, copy, telemetry);
            if (symbolic_alt) ++telemetry.symbolic_split_records;
        } else {
            clear_split_info_annotations(header, copy);
            compact_split_fields(header, copy, selected, original_alleles, telemetry);
            recompute_split_annotations(header, copy);
            if (keep_original_ac) preserve_original_counts(header, source, copy, selected);
        }
        result.push_back(copy);
    }
    return result;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    htsFile* input = nullptr;
    htsFile* output = nullptr;
    bcf_hdr_t* input_header = nullptr;
    bcf_hdr_t* output_header = nullptr;
    bcf_hdr_t* writer_header = nullptr;
    bcf1_t* record = nullptr;
    faidx_t* reference = nullptr;
    std::uint64_t input_records = 0;
    std::uint64_t normalized_records = 0;
    std::uint64_t interval_skipped = 0;
    std::uint64_t exclusion_skipped = 0;
    GenotypeKernelTelemetry genotype_kernel_telemetry;
    fastgatk::io::IntervalFileStats interval_file_stats;
    fastgatk::io::IntervalFileStats exclusion_file_stats;
    try {
        reference = fai_load(options.reference.c_str());
        if (!reference) throw std::runtime_error("BAD_INPUT: cannot load FASTA index: " + options.reference);
        input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open variant input: " + options.input);
        input_header = bcf_hdr_read(input);
        if (!input_header) throw std::runtime_error("BAD_INPUT: cannot read variant header");
        output_header = bcf_hdr_dup(input_header);
        if (!output_header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate variant header");
        if (options.split_multi_allelics) {
            // GATK adds the standard chromosome-count declarations whenever
            // splitting is enabled, even when the input header omitted them.
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC") < 0)
                bcf_hdr_append(output_header,
                    "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count in genotypes>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF") < 0)
                bcf_hdr_append(output_header,
                    "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele Frequency>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN") < 0)
                bcf_hdr_append(output_header,
                    "##INFO=<ID=AN,Number=1,Type=Integer,Description=Total number of alleles>");
        }
        if (options.keep_original_ac) {
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AC_Orig,Number=A,Type=Integer,Description=Original AC>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AF_Orig,Number=A,Type=Float,Description=Original AF>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN_Orig") < 0)
                bcf_hdr_append(output_header, "##INFO=<ID=AN_Orig,Number=1,Type=Integer,Description=Original AN>");
        }
        bcf_hdr_append(output_header, "##source=fastgatk-left-align-trim");
        bcf_hdr_append(output_header, "##fastgatk_left_align_status=prototype-trim-left-shift");
        if (bcf_hdr_sync(output_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync output header");
        // Keep the input-side FORMAT/sample payload available while records
        // are normalized and split.  The duplicate writer header is reduced
        // only at the final boundary, matching GATK's sites-only writer.
        writer_header = bcf_hdr_dup(output_header);
        if (!writer_header)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate sites-only writer header");
        if (options.sites_only_vcf_output &&
            bcf_hdr_set_samples(writer_header, nullptr, 0) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot configure sites-only VCF header");
        std::vector<Region> regions;
        regions.reserve(options.regions.size());
        bool first_selector = true;
        for (const auto& text : options.regions) {
            append_left_align_region_selector_with_rule(
                text, input_header, regions, interval_file_stats,
                options.interval_set_rule, first_selector);
            first_selector = false;
        }
        std::vector<Region> excluded_regions;
        excluded_regions.reserve(options.excluded_regions.size());
        bool first_exclusion_selector = true;
        for (const auto& text : options.excluded_regions) {
            append_left_align_region_selector_with_rule(
                text, input_header, excluded_regions, exclusion_file_stats,
                fastgatk::io::HtsIntervalSetRule::Union, first_exclusion_selector);
            first_exclusion_selector = false;
        }
        apply_region_padding(regions, reference, input_header, options.interval_padding);
        apply_region_padding(excluded_regions, reference, input_header, options.exclusion_padding);
        normalize_left_align_regions(regions);
        normalize_left_align_regions(excluded_regions);
        output = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
        if (!output) throw std::runtime_error("cannot open normalized output: " + options.output);
        if (bcf_hdr_write(output, writer_header) != 0) throw std::runtime_error("cannot write normalized VCF header");
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
        int previous_rid = -1;
        int previous_pos = -1;
        int previous_end = -1;  // one-based inclusive end of the last written record
        while (fastgatk::io::read_variant_record(input, input_header, record,
                                               options.input) == 0) {
            ++input_records;
            if (!in_regions(input_header, record, regions)) {
                ++interval_skipped;
                bcf_clear(record);
                continue;
            }
            if (overlaps_regions(input_header, record, excluded_regions)) {
                ++exclusion_skipped;
                bcf_clear(record);
                continue;
            }
            auto split_records = options.split_multi_allelics
                ? split_record(output_header, record, options.keep_original_ac,
                               genotype_kernel_telemetry)
                : std::vector<bcf1_t*>{bcf_dup(record)};
            // In split mode a non-variant record intentionally yields no
            // output.  A failed duplication remains an actual resource
            // error, while an empty split is the normal GATK contract.
            if (options.split_multi_allelics && split_records.empty()) {
                bcf_clear(record);
                continue;
            }
            if (split_records.empty() || split_records.front() == nullptr) {
                for (auto* value : split_records) bcf_destroy(value);
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate normalized VCF record");
            }
            for (auto* split : split_records) {
                if (!split) {
                    for (auto* value : split_records) if (value) bcf_destroy(value);
                    throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate normalized VCF record");
                }
                bcf_unpack(split, BCF_UN_STR | BCF_UN_FMT | BCF_UN_INFO);
                const int split_rid = split->rid;
                const int split_start = split->pos + 1;
                std::int64_t window = options.max_leading_bases;
                if (split_rid == previous_rid && previous_end >= 1) {
                    // This is exactly VariantWalker.apply's
                    // min(maxLeadingBases, distanceToLastVariant - 1), with
                    // one-based inclusive coordinates.  The bound uses the
                    // record as written, including a record skipped for
                    // oversized indel length.
                    const std::int64_t distance = static_cast<std::int64_t>(split_start) - previous_end;
                    window = std::min<std::int64_t>(window, distance - 1);
                }
                if (window > std::numeric_limits<int>::max()) window = std::numeric_limits<int>::max();
                if (window < std::numeric_limits<int>::min()) window = std::numeric_limits<int>::min();
                const bool normalize = !options.split_multi_allelics || split->n_allele == 2;
                const bool changed = normalize && normalize_record(
                    reference, output_header, split, options.dont_trim,
                    options.max_indel_length, static_cast<int>(window),
                    genotype_kernel_telemetry);
                if (changed) ++normalized_records;
                const char* contig = split->rid >= 0 ? bcf_hdr_id2name(input_header, split->rid) : nullptr;
                const int output_rid = contig ? bcf_hdr_name2id(output_header, contig) : -1;
                if (output_rid < 0 || (output_rid < previous_rid) ||
                    (output_rid == previous_rid && split->pos < previous_pos)) {
                    for (auto* value : split_records) bcf_destroy(value);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: normalized records are not coordinate-sorted");
                }
                previous_rid = output_rid;
                previous_pos = split->pos;
                previous_end = record_end_exclusive(output_header, split);
                if (options.sites_only_vcf_output && bcf_subset_format(writer_header, split) != 0) {
                    for (auto* value : split_records) if (value && value != split) bcf_destroy(value);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot subset sites-only FORMAT payload");
                }
                if (bcf_translate(writer_header, input_header, split) != 0 ||
                    bcf_write(output, writer_header, split) != 0) {
                    for (auto* value : split_records) if (value && value != split) bcf_destroy(value);
                    throw std::runtime_error("cannot write normalized VCF record");
                }
                bcf_destroy(split);
            }
            bcf_clear(record);
        }
        bcf_destroy(record); record = nullptr;
        bcf_close(output); output = nullptr;
        bcf_close(input); input = nullptr;
        bcf_hdr_destroy(writer_header); writer_header = nullptr;
        std::string index_path;
        if (options.create_index) {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build normalized VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const bool primary_complete = file_complete(options.output);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !index_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: normalized output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"LeftAlignAndTrimVariants\",\"implementation\":\"fastgatk-left-align-trim\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"common_trim\":" << (options.dont_trim ? "false" : "true")
                 << ",\"max_indel_length\":" << options.max_indel_length
                 << ",\"max_leading_bases\":" << options.max_leading_bases
                 << ",\"split_multiallelics\":" << (options.split_multi_allelics ? "true" : "false")
                 << ",\"keep_original_ac\":" << (options.keep_original_ac ? "true" : "false")
                 << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"interval_subset\":" << (!regions.empty() ? "true" : "false")
                 << ",\"interval_exclusion\":" << (!excluded_regions.empty() ? "true" : "false")
                 << ",\"interval_padding\":" << options.interval_padding
                 << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                 << ",\"interval_merging_rule\":\"" << options.interval_merging_rule << "\""
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"repeat_left_shift\":true,\"vcf_index\":" << (index_path.empty() ? "false" : "true")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\"" << json_escape(options.output)
                 << "\",\"kind\":\"vcf\",\"complete\":true}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << input_records << ",\"normalized_records\":" << normalized_records
                 << ",\"intervals\":" << regions.size() << ",\"interval_skipped\":" << interval_skipped
                 << ",\"excluded_intervals\":" << excluded_regions.size()
                 << ",\"exclusion_skipped\":" << exclusion_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"exclusion_interval_list_inputs\":" << exclusion_file_stats.files
                 << ",\"exclusion_interval_list_records\":" << exclusion_file_stats.records
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
                 << ",\"max_ploidy\":" << genotype_kernel_telemetry.max_ploidy
                 << ",\"symbolic_split_records\":" << genotype_kernel_telemetry.symbolic_split_records
                 << ",\"symbolic_info_fields_cleared\":"
                 << genotype_kernel_telemetry.symbolic_info_fields_cleared
                 << ",\"symbolic_format_fields_cleared\":"
                 << genotype_kernel_telemetry.symbolic_format_fields_cleared
                 << ",\"oversized_indel_records\":" << genotype_kernel_telemetry.oversized_indel_records
                 << ",\"left_shift_bases\":" << genotype_kernel_telemetry.left_shift_bases << "}}\n";
        std::cout << "{\"tool\":\"LeftAlignAndTrimVariants\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"normalized_records\":" << normalized_records
                  << ",\"interval_skipped\":" << interval_skipped
                  << ",\"exclusion_skipped\":" << exclusion_skipped
                  << ",\"interval_padding\":" << options.interval_padding
                  << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                  << ",\"interval_merging_rule\":\"" << options.interval_merging_rule << "\""
                  << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"allele_field_remap_kernel_calls\":"
                  << genotype_kernel_telemetry.allele_field_remap_calls
                  << ",\"max_ploidy\":" << genotype_kernel_telemetry.max_ploidy
                  << ",\"symbolic_split_records\":" << genotype_kernel_telemetry.symbolic_split_records
                  << ",\"symbolic_info_fields_cleared\":"
                  << genotype_kernel_telemetry.symbolic_info_fields_cleared
                  << ",\"symbolic_format_fields_cleared\":"
                  << genotype_kernel_telemetry.symbolic_format_fields_cleared
                  << ",\"oversized_indel_records\":" << genotype_kernel_telemetry.oversized_indel_records
                  << ",\"left_shift_bases\":" << genotype_kernel_telemetry.left_shift_bases << "}\n";
        fai_destroy(reference);
        bcf_hdr_destroy(output_header);
        return 0;
    } catch (...) {
        bcf_destroy(record);
        if (output) bcf_close(output);
        if (input) bcf_close(input);
        if (output_header) bcf_hdr_destroy(output_header);
        if (writer_header) bcf_hdr_destroy(writer_header);
        if (input_header) bcf_hdr_destroy(input_header);
        if (reference) fai_destroy(reference);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for LeftAlignAndTrimVariants");
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
        std::cerr << "fastgatk-left-align-trim: " << error.what() << '\n';
        return 2;
    }
}
