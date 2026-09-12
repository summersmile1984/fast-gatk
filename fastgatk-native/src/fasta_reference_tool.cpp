#include "fastgatk/core/plan.hpp"
#include "fastgatk/reference_io.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include "fastgatk/io/hts_read_guard.hpp"
#ifdef FASTGATK_FASTA_ALTERNATE
#include <htslib/vcf.h>
#endif

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Region {
    int tid = -1;
    std::string contig;
    hts_pos_t start = 1;
    hts_pos_t end = 0;
};

struct Options {
    std::string reference;
    std::string output;
    std::vector<std::string> selectors;
    std::string manifest;
    int line_width = 60;
    int threads = 1;
#ifdef FASTGATK_FASTA_ALTERNATE
    std::string variants;
    std::string snp_mask;
    bool snp_mask_priority = false;
    std::string iupac_sample;
#endif
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string value_for(int& index, int argc, char** argv, const std::string& argument,
                     std::string_view name, const char* short_name = nullptr) {
    if (const auto value = inline_value(argument, name); !value.empty()) return value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument("missing value for " + std::string(name));
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
#ifdef FASTGATK_FASTA_ALTERNATE
            std::cout << "fastgatk-fasta-alternate-reference-maker (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE       reference FASTA (+ .fai)\n"
                         "  -V, --variant FILE         input VCF/BCF\n"
                         "  -L, --intervals REGION     interval or interval-list (repeatable)\n"
                         "      --snp-mask FILE       optional SNP-mask VCF\n"
                         "      --snp-mask-priority  mask takes priority over variants\n"
                         "      --use-iupac-sample S  encode diploid heterozygous SNPs\n"
                         "  -O, --output FILE          output FASTA\n"
                         "      --line-width N         bases per sequence line (default 60)\n"
                         "      --output-manifest FILE OutputManifest JSON\n";
#else
            std::cout << "fastgatk-fasta-reference-maker (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE       reference FASTA (+ .fai)\n"
                         "  -L, --intervals REGION     interval or interval-list (repeatable)\n"
                         "  -O, --output FILE          output FASTA\n"
                         "      --line-width N         bases per sequence line (default 60)\n"
                         "      --output-manifest FILE OutputManifest JSON\n";
#endif
            std::exit(0);
        } else if (argument == "-R" || argument == "--reference" ||
                   !inline_value(argument, "--reference").empty()) {
            options.reference = value_for(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || argument == "--intervals" || argument == "--interval" ||
                   argument == "--region" || !inline_value(argument, "--intervals").empty() ||
                   !inline_value(argument, "--interval").empty() ||
                   !inline_value(argument, "--region").empty()) {
            const std::string_view name = argument == "-L" ? "--intervals" :
                (argument == "--region" || !inline_value(argument, "--region").empty() ? "--region" :
                 (argument == "--interval" || !inline_value(argument, "--interval").empty() ? "--interval" : "--intervals"));
            options.selectors.push_back(value_for(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-O" || argument == "--output" || !inline_value(argument, "--output").empty()) {
            options.output = value_for(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--line-width" || !inline_value(argument, "--line-width").empty()) {
            options.line_width = std::stoi(value_for(index, argc, argv, argument, "--line-width"));
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(value_for(index, argc, argv, argument, "--threads"));
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() ||
                   !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = value_for(index, argc, argv, argument, name);
#ifdef FASTGATK_FASTA_ALTERNATE
        } else if (argument == "-V" || argument == "--variant" || argument == "--variants" ||
                   !inline_value(argument, "--variant").empty() ||
                   !inline_value(argument, "--variants").empty()) {
            const auto name = argument == "-V" || argument == "--variant" ||
                    !inline_value(argument, "--variant").empty()
                ? "--variant" : "--variants";
            options.variants = value_for(index, argc, argv, argument, name, "-V");
        } else if (argument == "--snp-mask" || !inline_value(argument, "--snp-mask").empty()) {
            options.snp_mask = value_for(index, argc, argv, argument, "--snp-mask");
        } else if (argument == "--snp-mask-priority") {
            options.snp_mask_priority = true;
        } else if (argument == "--use-iupac-sample" || !inline_value(argument, "--use-iupac-sample").empty()) {
            options.iupac_sample = value_for(index, argc, argv, argument, "--use-iupac-sample");
#endif
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--seconds-between-progress-updates", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                     "--seconds-between-progress-updates");
                (void)value_for(index, argc, argv, argument, name);
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.line_width < 1) throw std::invalid_argument("--line-width must be positive");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
#ifdef FASTGATK_FASTA_ALTERNATE
    if (options.variants.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.snp_mask_priority && options.snp_mask.empty())
        throw std::invalid_argument("--snp-mask-priority requires --snp-mask");
#endif
    return options;
}

std::string without_commas(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), ','), value.end());
    return value;
}

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while (stream >> field) fields.push_back(field);
    return fields;
}

std::vector<std::string> read_text_lines(const std::string& path) {
    std::vector<std::string> lines;
    if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        kstring_t line{0, 0, nullptr};
        while (fastgatk::io::read_text_line(file, &line, path) >= 0)
            lines.emplace_back(line.s == nullptr ? "" : std::string(line.s, line.l));
        free(line.s);
        if (hts_close(file) != 0) throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
        return lines;
    }
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
    for (std::string line; std::getline(input, line);) lines.push_back(std::move(line));
    return lines;
}

int contig_id(const std::string& contig, const faidx_t* fai) {
    for (int tid = 0; tid < faidx_nseq(fai); ++tid)
        if (contig == faidx_iseq(fai, tid)) return tid;
    return -1;
}

Region parse_interval(const std::string& token, const faidx_t* fai, bool bed = false) {
    const std::string text = without_commas(token);
    const auto colon = text.find(':');
    const std::string contig = colon == std::string::npos ? text : text.substr(0, colon);
    const int tid = contig_id(contig, fai);
    if (tid < 0) throw std::invalid_argument("BAD_INPUT: interval contig is absent from reference: " + contig);
    const auto length = static_cast<hts_pos_t>(faidx_seq_len64(fai, contig.c_str()));
    if (colon == std::string::npos) return {tid, contig, 1, length};
    const auto dash = text.find('-', colon + 1);
    const auto begin_text = text.substr(colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1);
    const auto end_text = dash == std::string::npos ? begin_text : text.substr(dash + 1);
    if (begin_text.empty() || end_text.empty()) throw std::invalid_argument("BAD_INPUT: malformed interval: " + token);
    const auto begin = static_cast<hts_pos_t>(std::stoll(begin_text));
    const auto end = static_cast<hts_pos_t>(std::stoll(end_text));
    if (bed) {
        if (begin < 0 || end <= begin || end > length)
            throw std::invalid_argument("BAD_INPUT: interval outside reference: " + token);
        return {tid, contig, begin + 1, end};
    }
    if (begin < 1 || end < begin || end > length)
        throw std::invalid_argument("BAD_INPUT: interval outside reference: " + token);
    return {tid, contig, begin, end};
}

void append_interval_selector(const std::string& selector, const faidx_t* fai,
                              std::vector<Region>& regions, std::size_t& interval_files,
                              std::size_t& interval_records) {
    std::error_code error;
    const bool regular = std::filesystem::is_regular_file(selector, error) && !error;
    if (!regular) {
        regions.push_back(parse_interval(selector, fai));
        return;
    }
    ++interval_files;
    const bool bed = selector.size() >= 4 &&
        (selector.substr(selector.size() - 4) == ".bed" ||
         (selector.size() >= 7 && selector.substr(selector.size() - 7) == ".bed.gz"));
    std::size_t records_this_file = 0;
    for (const auto& line : read_text_lines(selector)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '@') continue;
        const auto fields = split_fields(line.substr(first));
        if (fields.empty() || fields.front() == "track" || fields.front() == "browser") continue;
        if (bed) {
            if (fields.size() < 3) throw std::invalid_argument("BAD_INPUT: malformed BED line in " + selector);
            regions.push_back(parse_interval(fields[0] + ":" + fields[1] + "-" + fields[2], fai, true));
        } else if (fields.size() >= 3 && fields[0].find(':') == std::string::npos) {
            regions.push_back(parse_interval(fields[0] + ":" + fields[1] + "-" + fields[2], fai));
        } else {
            regions.push_back(parse_interval(fields.front(), fai));
        }
        ++interval_records;
        ++records_this_file;
    }
    if (records_this_file == 0) throw std::invalid_argument("BAD_INPUT: interval file is empty: " + selector);
}

std::vector<Region> selected_regions(const Options& options, const faidx_t* fai,
                                     std::size_t& interval_files, std::size_t& interval_records) {
    std::vector<Region> regions;
    if (options.selectors.empty()) {
        regions.reserve(static_cast<std::size_t>(faidx_nseq(fai)));
        for (int tid = 0; tid < faidx_nseq(fai); ++tid) {
            const char* contig = faidx_iseq(fai, tid);
            regions.push_back({tid, contig, 1, faidx_seq_len64(fai, contig)});
        }
    } else {
        for (const auto& selector : options.selectors)
            append_interval_selector(selector, fai, regions, interval_files, interval_records);
    }
    std::stable_sort(regions.begin(), regions.end(), [](const Region& lhs, const Region& rhs) {
        if (lhs.tid != rhs.tid) return lhs.tid < rhs.tid;
        if (lhs.start != rhs.start) return lhs.start < rhs.start;
        return lhs.end < rhs.end;
    });
    std::vector<Region> merged;
    for (const auto& region : regions) {
        if (!merged.empty() && merged.back().tid == region.tid &&
            region.start <= merged.back().end + 1) {
            merged.back().end = std::max(merged.back().end, region.end);
        } else {
            merged.push_back(region);
        }
    }
    return merged;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char character : value) {
        if (character == '"' || character == '\\') out << '\\';
        if (character == '\n') out << "\\n";
        else if (character == '\r') out << "\\r";
        else if (character == '\t') out << "\\t";
        else out << character;
    }
    return out.str();
}

struct KernelStats {
    std::size_t records = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

std::string copy_reference_bytes_kokkos(const std::string& sequence, int threads, KernelStats& stats) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<unsigned char*> device_input("fasta_reference_input", sequence.size());
    Kokkos::View<unsigned char*> device_output("fasta_reference_output", sequence.size());
    auto host = Kokkos::create_mirror_view(device_input);
    for (std::size_t index = 0; index < sequence.size(); ++index)
        host(index) = static_cast<unsigned char>(std::toupper(static_cast<unsigned char>(sequence[index])));
    fastgatk::core::HostBatch batch("fasta-reference-maker-v1");
    batch.records = sequence.size();
    batch.bytes = sequence.size();
    fastgatk::core::KernelPlan<ExecSpace> plan("fasta-reference-copy");
    plan.begin_prepare(batch);
    Kokkos::deep_copy(device_input, host);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(sequence.size());
    device_batch.bind("input", device_input);
    device_batch.bind("output", device_output);
    plan.end_prepare(device_batch);
    plan.begin_execute();
    Kokkos::parallel_for("fasta_reference_copy", Kokkos::RangePolicy<ExecSpace>(0, sequence.size()),
        KOKKOS_LAMBDA(const std::size_t index) { device_output(index) = device_input(index); });
    ExecSpace().fence();
    plan.end_execute();
    auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_output);
    std::string copied(sequence.size(), '\0');
    for (std::size_t index = 0; index < sequence.size(); ++index)
        copied[index] = static_cast<char>(result(index));
    stats.records += sequence.size();
    stats.prepare_seconds += plan.telemetry().prepare_seconds;
    stats.execute_seconds += plan.telemetry().execute_seconds;
    (void)threads;
    return copied;
}

std::string md5_for_bytes(const std::string& sequence) {
    hts_md5_context* context = hts_md5_init();
    if (!context) throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib MD5 context allocation failed");
    hts_md5_update(context, sequence.data(), static_cast<unsigned long>(sequence.size()));
    unsigned char digest[16]{};
    hts_md5_final(digest, context);
    hts_md5_destroy(context);
    char hex[33]{};
    hts_md5_hex(hex, digest);
    return std::string(hex);
}

void write_dictionary(const std::string& fasta,
                      const std::vector<std::pair<std::string, std::pair<std::int64_t, std::string>>>& sequences) {
    std::ofstream dictionary(fastgatk::reference::dictionary_path(fasta));
    if (!dictionary) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sequence dictionary: " + fasta);
    dictionary << "@HD\tVN:1.6\n";
    for (const auto& [name, info] : sequences)
        dictionary << "@SQ\tSN:" << name << "\tLN:" << info.first << "\tM5:" << info.second << "\n";
}

#ifdef FASTGATK_FASTA_ALTERNATE
struct Variant {
    std::string ref;
    std::string alt;
    std::string iupac;
    bool filtered = false;
    // HTSJDK classifies a VariantContext as SNP/ simple indel at the record
    // level.  Keep that information while flattening ALT alleles so a mixed
    // or multi-allelic record is not accidentally applied as an independent
    // simple allele.
    bool record_is_snp = false;
    bool record_is_simple_indel = false;
    bool simple_snp = false;
    bool simple_insertion = false;
    bool simple_deletion = false;
};

using VariantMap = std::unordered_map<std::string, std::map<hts_pos_t, std::vector<Variant>>>;

char iupac_code(char first, char second) {
    first = static_cast<char>(std::toupper(static_cast<unsigned char>(first)));
    second = static_cast<char>(std::toupper(static_cast<unsigned char>(second)));
    if (first > second) std::swap(first, second);
    if (first == second) return first;
    const std::string pair{first, second};
    if (pair == "AC") return 'M';
    if (pair == "AG") return 'R';
    if (pair == "AT") return 'W';
    if (pair == "CG") return 'S';
    if (pair == "CT") return 'Y';
    if (pair == "GT") return 'K';
    return 'N';
}

std::string sample_iupac(const bcf_hdr_t* header, bcf1_t* record, int allele_count,
                         const std::string& sample) {
    // IUPAC is restricted to diploid genotypes, not to biallelic records.
    // A diploid 1/2 genotype at a multi-ALT SNP is still encoded from the two
    // concrete one-base alleles by GATK's getIUPACBase().
    if (sample.empty() || allele_count < 2) return {};
    const int sample_index = bcf_hdr_id2int(header, BCF_DT_SAMPLE, sample.c_str());
    if (sample_index < 0) throw std::invalid_argument("BAD_INPUT: IUPAC sample is absent from VCF: " + sample);
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const int values = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    std::string result;
    if (values >= (sample_index + 1) * 2 && genotype_count >= (sample_index + 1) * 2) {
        const int first = bcf_gt_allele(genotypes[sample_index * 2]);
        const int second = bcf_gt_allele(genotypes[sample_index * 2 + 1]);
        if (first >= 0 && second >= 0 && first < allele_count && second < allele_count) {
            const char* first_allele = record->d.allele[first];
            const char* second_allele = record->d.allele[second];
            if (first_allele && second_allele && std::strlen(first_allele) == 1 && std::strlen(second_allele) == 1) {
                // GATK's getIUPACBase() returns the selected allele for a
                // homozygous diploid genotype (including hom-ref/hom-alt),
                // rather than falling back to the first ALT.  This matters
                // for multi-ALT records where the selected homozygous ALT is
                // not ALT[0].  Heterozygous calls continue through the
                // canonical IUPAC lookup below.
                result.push_back(first == second ? first_allele[0] :
                                 iupac_code(first_allele[0], second_allele[0]));
            }
        }
    }
    free(genotypes);
    return result;
}

VariantMap read_variants(const std::string& path, const std::string& iupac_sample,
                         std::unordered_map<std::string, bool>* samples_seen = nullptr) {
    htsFile* file = bcf_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BACKEND_UNAVAILABLE: cannot open VCF/BCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(file);
    if (!header) {
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + path);
    }
    if (samples_seen && !iupac_sample.empty()) (*samples_seen)[iupac_sample] = bcf_hdr_id2int(header, BCF_DT_SAMPLE, iupac_sample.c_str()) >= 0;
    VariantMap result;
    bcf1_t* record = bcf_init();
    while (record && fastgatk::io::read_variant_record(file, header, record, path) == 0) {
        bcf_unpack(record, BCF_UN_STR | BCF_UN_FLT);
        if (record->rid < 0 || record->n_allele < 2) {
            bcf_clear(record);
            continue;
        }
        const char* contig = bcf_hdr_id2name(header, record->rid);
        const char* ref = record->d.allele[0];
        if (!contig || !ref) {
            bcf_clear(record);
            continue;
        }
        const std::string reference_allele(ref);
        // Mirror htsjdk's VariantContext type classification for the subset
        // consumed by FastaAlternateReferenceMaker.  A record with a SNP and
        // an indel (or a symbolic ALT) is MIXED and must not be treated as a
        // simple SNP merely because one flattened ALT is one base long.
        bool record_is_snp = reference_allele.size() == 1;
        bool record_is_simple_indel = record->n_allele == 2 && !reference_allele.empty();
        if (record->n_allele < 2) record_is_snp = false;
        for (int allele_index = 1; allele_index < record->n_allele; ++allele_index) {
            const char* alt = record->d.allele[allele_index];
            const std::string alternate = alt == nullptr ? std::string{} : std::string(alt);
            const bool symbolic = alternate.empty() || alternate[0] == '<' || alternate == "*";
            if (symbolic || alternate.size() != reference_allele.size() || alternate.size() != 1)
                record_is_snp = false;
            if (record->n_allele != 2 || symbolic || alternate.empty() ||
                alternate[0] != reference_allele[0] ||
                alternate.size() == reference_allele.size() ||
                (reference_allele.size() != 1 && alternate.size() != 1))
                record_is_simple_indel = false;
        }
        for (int allele_index = 1; allele_index < record->n_allele; ++allele_index) {
            const char* alt = record->d.allele[allele_index];
            if (!alt || alt[0] == '<' || std::string(alt) == "*") continue;
            Variant variant;
            variant.ref = ref;
            variant.alt = alt;
            variant.record_is_snp = record_is_snp;
            variant.record_is_simple_indel = record_is_simple_indel;
            // HTSlib represents the VCF PASS token as one filter ID with
            // numeric value zero; actual FILTER labels have positive IDs.
            variant.filtered = record->d.n_flt > 0 &&
                !(record->d.n_flt == 1 && record->d.flt[0] == 0);
            variant.simple_snp = record_is_snp && variant.ref.size() == 1 && variant.alt.size() == 1;
            variant.simple_insertion = record_is_simple_indel &&
                variant.alt.size() > variant.ref.size() &&
                variant.alt.compare(0, variant.ref.size(), variant.ref) == 0;
            variant.simple_deletion = record_is_simple_indel &&
                variant.ref.size() > variant.alt.size() &&
                variant.ref.compare(0, variant.alt.size(), variant.alt) == 0;
            if (!(variant.simple_snp || variant.simple_insertion || variant.simple_deletion)) continue;
            if (variant.simple_snp && !iupac_sample.empty())
                variant.iupac = sample_iupac(header, record, record->n_allele, iupac_sample);
            result[contig][static_cast<hts_pos_t>(record->pos + 1)].push_back(std::move(variant));
        }
        bcf_clear(record);
    }
    if (record) bcf_destroy(record);
    bcf_hdr_destroy(header);
    const int close_status = bcf_close(file);
    if (close_status != 0) throw std::runtime_error("BAD_INPUT: failed reading VCF/BCF: " + path);
    return result;
}

bool has_snp(const VariantMap& masks, const std::string& contig, hts_pos_t position) {
    const auto contig_it = masks.find(contig);
    if (contig_it == masks.end()) return false;
    const auto position_it = contig_it->second.find(position);
    if (position_it == contig_it->second.end()) return false;
    return std::any_of(position_it->second.begin(), position_it->second.end(),
                       [](const Variant& variant) {
                           // GATK's FastaAlternateReferenceMaker deliberately
                           // does not inspect FILTER for the SNP mask: a
                           // filtered SNP is still a mask feature.  It also
                           // asks VariantContext::isSNP(), not whether one
                           // flattened ALT happens to be one base long.
                           return variant.record_is_snp;
                       });
}

std::string apply_alternate(const Region& region, const std::string& reference, const VariantMap& variants,
                            const VariantMap& mask, const Options& options) {
    std::string output;
    output.reserve(reference.size());
    std::size_t deletion_bases_remaining = 0;
    const auto variant_contig = variants.find(region.contig);
    const auto mask_contig = mask.find(region.contig);
    for (std::size_t offset = 0; offset < reference.size(); ++offset) {
        const hts_pos_t position = region.start + static_cast<hts_pos_t>(offset);
        if (deletion_bases_remaining > 0) {
            --deletion_bases_remaining;
            continue;
        }
        const bool masked = mask_contig != mask.end() && has_snp(mask, region.contig, position);
        if (options.snp_mask_priority && masked) {
            output.push_back('N');
            continue;
        }
        const std::vector<Variant>* calls = nullptr;
        if (variant_contig != variants.end()) {
            const auto call_it = variant_contig->second.find(position);
            if (call_it != variant_contig->second.end()) calls = &call_it->second;
        }
        bool replaced = false;
        if (calls) {
            for (const auto& call : *calls) {
                if (call.filtered) continue;
                if (call.simple_deletion) {
                    deletion_bases_remaining = call.ref.size() - 1;
                    output.push_back(reference[offset]);
                    replaced = true;
                    break;
                }
                if (call.simple_insertion) {
                    output.append(call.alt);
                    replaced = true;
                    break;
                }
                if (call.simple_snp) {
                    output.append(call.iupac.empty() ? call.alt : call.iupac);
                    replaced = true;
                    break;
                }
            }
        }
        if (!replaced) {
            if (!options.snp_mask_priority && masked) output.push_back('N');
            else output.push_back(reference[offset]);
        }
    }
    return output;
}
#endif

void write_fasta(const Options& options, const std::vector<Region>& regions,
                 faidx_t* fai, KernelStats& kernel_stats,
                 std::vector<std::pair<std::string, std::pair<std::int64_t, std::string>>>& dictionary_entries
#ifdef FASTGATK_FASTA_ALTERNATE
                 , const VariantMap& variants, const VariantMap& mask
#endif
                 ) {
    std::ofstream output(options.output, std::ios::binary);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write FASTA: " + options.output);
    std::size_t sequence_id = 0;
    for (const auto& region : regions) {
        hts_pos_t fetched = 0;
        char* raw = faidx_fetch_seq64(fai, region.contig.c_str(), region.start - 1, region.end - 1, &fetched);
        if (!raw || fetched != region.end - region.start + 1) {
            free(raw);
            throw std::runtime_error("BAD_INPUT: failed to fetch reference interval: " + region.contig);
        }
        std::string sequence(raw, static_cast<std::size_t>(fetched));
        free(raw);
#ifdef FASTGATK_FASTA_ALTERNATE
        sequence = apply_alternate(region, sequence, variants, mask, options);
        sequence = copy_reference_bytes_kokkos(sequence, options.threads, kernel_stats);
#else
        sequence = copy_reference_bytes_kokkos(sequence, options.threads, kernel_stats);
#endif
        if (sequence.empty()) continue;
        ++sequence_id;
        const std::string sequence_name = std::to_string(sequence_id);
        output << ">" << sequence_name << " " << region.contig << ":" << region.start << "-" << region.end << "\n";
        for (std::size_t offset = 0; offset < sequence.size(); offset += static_cast<std::size_t>(options.line_width))
            output.write(sequence.data() + static_cast<std::streamoff>(offset),
                         static_cast<std::streamsize>(std::min<std::size_t>(options.line_width, sequence.size() - offset))) << '\n';
        dictionary_entries.push_back({sequence_name, {static_cast<std::int64_t>(sequence.size()), md5_for_bytes(sequence)}});
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed finalizing FASTA: " + options.output);
    if (fai_build(options.output.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed creating FASTA index: " + options.output + ".fai");
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    faidx_t* fai = fai_load(options.reference.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
    try {
        std::size_t interval_files = 0;
        std::size_t interval_records = 0;
        const auto regions = selected_regions(options, fai, interval_files, interval_records);
#ifdef FASTGATK_FASTA_ALTERNATE
        std::unordered_map<std::string, bool> samples_seen;
        const auto variants = read_variants(options.variants, options.iupac_sample, &samples_seen);
        if (!options.iupac_sample.empty() && !samples_seen[options.iupac_sample])
            throw std::invalid_argument("BAD_INPUT: IUPAC sample is absent from VCF: " + options.iupac_sample);
        const auto mask = options.snp_mask.empty() ? VariantMap{} : read_variants(options.snp_mask, "");
#endif
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(options.threads);
        Kokkos::initialize(settings);
        bool initialized = true;
        KernelStats kernel_stats;
        std::vector<std::pair<std::string, std::pair<std::int64_t, std::string>>> dictionary_entries;
        try {
            write_fasta(options, regions, fai, kernel_stats, dictionary_entries
#ifdef FASTGATK_FASTA_ALTERNATE
                        , variants, mask
#endif
                        );
            write_dictionary(options.output, dictionary_entries);
            if (!options.manifest.empty()) {
                std::ofstream manifest(options.manifest);
                if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
#ifdef FASTGATK_FASTA_ALTERNATE
                manifest << "{\"schema_version\":1,\"tool\":\"FastaAlternateReferenceMaker\",\"implementation\":\"fastgatk-fasta-alternate-reference-maker\",\"status\":\"prototype\",\"execution_space\":\""
                         << Kokkos::DefaultExecutionSpace::name() << "/Host/HTSlib\",\"records\":" << kernel_stats.records
                         << ",\"prepare_seconds\":" << kernel_stats.prepare_seconds << ",\"execute_seconds\":" << kernel_stats.execute_seconds << ",\"reference\":\""
                         << json_escape(options.reference) << "\",\"variants\":\"" << json_escape(options.variants)
                         << "\",\"regions\":" << regions.size() << ",\"interval_list_inputs\":" << interval_files
                         << ",\"interval_list_records\":" << interval_records << ",\"output\":\"" << json_escape(options.output)
                         << "\",\"dictionary\":\"" << json_escape(fastgatk::reference::dictionary_path(options.output)) << "\"}\n";
#else
                manifest << "{\"schema_version\":1,\"tool\":\"FastaReferenceMaker\",\"implementation\":\"fastgatk-fasta-reference-maker\",\"status\":\"prototype\",\"execution_space\":\""
                         << Kokkos::DefaultExecutionSpace::name() << "\",\"reference\":\"" << json_escape(options.reference)
                         << "\",\"regions\":" << regions.size() << ",\"interval_list_inputs\":" << interval_files
                         << ",\"interval_list_records\":" << interval_records << ",\"records\":" << kernel_stats.records
                         << ",\"prepare_seconds\":" << kernel_stats.prepare_seconds << ",\"execute_seconds\":" << kernel_stats.execute_seconds
                         << ",\"output\":\"" << json_escape(options.output) << "\",\"dictionary\":\""
                         << json_escape(fastgatk::reference::dictionary_path(options.output)) << "\"}\n";
#endif
            }
            Kokkos::finalize();
            initialized = false;
        } catch (...) {
            if (initialized) Kokkos::finalize();
            throw;
        }
        fai_destroy(fai);
        return 0;
    } catch (...) {
        fai_destroy(fai);
        throw;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return main_impl(argc, argv);
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
