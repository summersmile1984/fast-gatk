#include "fastgatk/reference_io.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::string variant;
    std::vector<std::string> references;
    std::string output;
    std::string manifest;
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv, const std::string& argument,
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
            std::cout << "fastgatk-check-reference-compatibility (GATK-compatible native contract-compatible)\n"
                         "  -I, --input FILE                 BAM/CRAM/SAM sequence dictionary\n"
                         "  -V, --variant FILE               VCF/BCF sequence dictionary\n"
                         "  -refcomp, --references-to-compare FILE  reference FASTA (repeatable)\n"
                         "  -O, --output FILE                compatibility table (optional)\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || argument == "--input" || !inline_value(argument, "--input").empty()) {
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-V" || argument == "--variant" || !inline_value(argument, "--variant").empty()) {
            options.variant = require_value(index, argc, argv, argument, "--variant", "-V");
        } else if (argument == "-refcomp" || argument == "--refcomp" || argument == "--references-to-compare" ||
                   !inline_value(argument, "--refcomp").empty() || !inline_value(argument, "--references-to-compare").empty()) {
            const auto name = argument.rfind("--refcomp", 0) == 0 ? "--refcomp" : "--references-to-compare";
            options.references.push_back(require_value(index, argc, argv, argument, name, "-refcomp"));
        } else if (argument == "-O" || argument == "--output" || !inline_value(argument, "--output").empty()) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() || !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = require_value(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0 || argument.rfind("--seconds-between-progress-updates", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" : "--seconds-between-progress-updates");
                (void)require_value(index, argc, argv, argument, name);
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty() == options.variant.empty())
        throw std::invalid_argument("provide exactly one of -I/--input or -V/--variant");
    if (options.references.empty()) throw std::invalid_argument("-refcomp/--references-to-compare is required");
    return options;
}

struct Dictionary {
    std::string name;
    std::vector<fastgatk::reference::SequenceInfo> sequences;
    bool md5s_present = false;
};

std::vector<fastgatk::reference::SequenceInfo> parse_vcf_text_contigs(const std::string& path) {
    std::vector<fastgatk::reference::SequenceInfo> result;
    htsFile* stream = hts_open(path.c_str(), "r");
    if (!stream) return result;
    kstring_t line{0, 0, nullptr};
    while (fastgatk::io::read_text_line(stream, &line, path) >= 0) {
        if (!line.s || std::string_view(line.s, line.l).rfind("##contig=<", 0) != 0) continue;
        const std::string text(line.s, line.l);
        const auto begin = text.find('<');
        const auto end = text.rfind('>');
        if (begin == std::string::npos || end <= begin) continue;
        fastgatk::reference::SequenceInfo sequence;
        std::istringstream fields(text.substr(begin + 1, end - begin - 1));
        for (std::string field; std::getline(fields, field, ',');) {
            const auto equal = field.find('=');
            if (equal == std::string::npos) continue;
            const auto key = field.substr(0, equal);
            const auto value = field.substr(equal + 1);
            if (key == "ID") sequence.name = value;
            else if (key == "length") sequence.length = std::stoll(value);
            else if (key == "md5" || key == "M5") sequence.md5 = value;
        }
        if (!sequence.name.empty() && sequence.length >= 0) result.push_back(std::move(sequence));
    }
    free(line.s);
    hts_close(stream);
    return result;
}

Dictionary bam_dictionary(const std::string& path) {
    samFile* file = sam_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open reads input: " + path);
    bam_hdr_t* header = sam_hdr_read(file);
    sam_close(file);
    if (!header) throw std::runtime_error("BAD_INPUT: reads input has no SAM header: " + path);
    Dictionary result;
    result.name = fastgatk::reference::basename(path);
    result.md5s_present = true;
    for (int index = 0; index < header->n_targets; ++index) {
        fastgatk::reference::SequenceInfo sequence{
            header->target_name[index], static_cast<std::int64_t>(header->target_len[index]), {}};
        kstring_t value{0, 0, nullptr};
        if (sam_hdr_find_tag_pos(header, "SQ", index, "M5", &value) == 0 && value.s)
            sequence.md5 = value.s;
        free(value.s);
        if (sequence.md5.empty()) result.md5s_present = false;
        result.sequences.push_back(std::move(sequence));
    }
    sam_hdr_destroy(header);
    return result;
}

Dictionary vcf_dictionary(const std::string& path) {
    htsFile* file = bcf_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open variant input: " + path);
    bcf_hdr_t* header = bcf_hdr_read(file);
    bcf_close(file);
    if (!header) throw std::runtime_error("BAD_INPUT: variant input has no VCF header: " + path);
    Dictionary result;
    result.name = fastgatk::reference::basename(path);
    const auto raw_contigs = parse_vcf_text_contigs(path);
    if (!raw_contigs.empty()) {
        result.sequences = raw_contigs;
        result.md5s_present = std::all_of(result.sequences.begin(), result.sequences.end(),
            [](const auto& sequence) { return !sequence.md5.empty(); });
        bcf_hdr_destroy(header);
        return result;
    }
    result.md5s_present = true;
    for (int index = 0; index < header->n[BCF_DT_CTG]; ++index) {
        const char* name = bcf_hdr_int2id(header, BCF_DT_CTG, index);
        fastgatk::reference::SequenceInfo sequence{name ? name : "", -1, {}};
        bcf_hrec_t* hrec = bcf_hdr_get_hrec(header, BCF_HL_CTG, "ID", sequence.name.c_str(), nullptr);
        if (hrec) {
            for (int key = 0; key < hrec->nkeys; ++key) {
                if (std::string_view(hrec->keys[key]) == "length") sequence.length = std::stoll(hrec->vals[key]);
                else if (std::string_view(hrec->keys[key]) == "md5") sequence.md5 = hrec->vals[key];
            }
        }
        if (sequence.length < 0) {
            bcf_hdr_destroy(header);
            throw std::invalid_argument("BAD_INPUT: VCF contig has no length: " + sequence.name);
        }
        if (sequence.md5.empty()) result.md5s_present = false;
        result.sequences.push_back(std::move(sequence));
    }
    // HTSlib intentionally drops malformed structured lines.  GATK's
    // VCF reader still exposes a dictionary for the common legacy spelling
    // ``M5:...`` (and treats that field as absent), so recover contig ID/LN
    // from the raw header when the parsed dictionary is empty.
    if (result.sequences.empty()) {
        kstring_t formatted{0, 0, nullptr};
        if (bcf_hdr_format(header, 0, &formatted) == 0 && formatted.s) {
            const char* raw = formatted.s;
            std::istringstream lines(raw);
            for (std::string line; std::getline(lines, line);) {
                if (line.rfind("##contig=<", 0) != 0) continue;
                const auto begin = line.find('<');
                const auto end = line.rfind('>');
                if (begin == std::string::npos || end <= begin) continue;
                fastgatk::reference::SequenceInfo sequence;
                std::string fields = line.substr(begin + 1, end - begin - 1);
                std::istringstream field_stream(fields);
                for (std::string field; std::getline(field_stream, field, ',');) {
                    const auto equal = field.find('=');
                    if (equal == std::string::npos) continue;
                    const auto key = field.substr(0, equal);
                    const auto value = field.substr(equal + 1);
                    if (key == "ID") sequence.name = value;
                    else if (key == "length") sequence.length = std::stoll(value);
                    else if (key == "md5" || key == "M5") sequence.md5 = value;
                }
                if (!sequence.name.empty() && sequence.length >= 0) {
                    if (sequence.md5.empty()) result.md5s_present = false;
                    result.sequences.push_back(std::move(sequence));
                }
            }
        }
        free(formatted.s);
    }
    if (result.sequences.empty()) result.md5s_present = false;
    bcf_hdr_destroy(header);
    return result;
}

Dictionary load_query(const Options& options) {
    const auto& path = options.input.empty() ? options.variant : options.input;
    // GATK accepts VCF/BCF paths whose compression suffix is not the
    // conventional `.vcf.gz`/`.bcf` spelling (for example `.vcf.bgz` or an
    // object-store staged name).  Dispatch by HTSlib's content format rather
    // than the filename so the native adapter preserves that input contract.
    htsFile* probe = hts_open(path.c_str(), "r");
    if (!probe) throw std::runtime_error("BAD_INPUT: cannot open query input: " + path);
    const auto* format = hts_get_format(probe);
    const bool is_variant = format != nullptr &&
        (format->format == vcf || format->format == bcf);
    if (hts_close(probe) != 0)
        throw std::runtime_error("BAD_INPUT: cannot close query input: " + path);
    if (is_variant) return vcf_dictionary(path);
    return bam_dictionary(path);
}

std::unordered_map<std::string, fastgatk::reference::SequenceInfo> by_name(const Dictionary& dictionary) {
    std::unordered_map<std::string, fastgatk::reference::SequenceInfo> result;
    for (const auto& sequence : dictionary.sequences) result.emplace(sequence.name, sequence);
    return result;
}

struct Compatibility {
    std::string reference;
    std::string status;
    std::string summary;
};

Compatibility compare(const Dictionary& query, const fastgatk::reference::FastaInfo& reference,
                      const std::string& query_name, bool use_md5) {
    Dictionary target{fastgatk::reference::basename(reference.path), reference.sequences, true};
    const auto query_map = by_name(query);
    const auto target_map = by_name(target);
    bool identical = query.sequences.size() == target.sequences.size();
    bool query_subset = true;
    bool names_lengths_match = true;
    bool sequence_mismatch = false;
    bool target_subset = true;
    std::vector<std::string> missing;
    for (const auto& sequence : query.sequences) {
        const auto it = target_map.find(sequence.name);
        if (it == target_map.end()) {
            query_subset = false;
            identical = false;
            names_lengths_match = false;
            continue;
        }
        if (it->second.length != sequence.length) {
            names_lengths_match = false;
            identical = false;
        }
        if (use_md5 && !sequence.md5.empty() && !it->second.md5.empty() && sequence.md5 != it->second.md5) {
            sequence_mismatch = true;
            identical = false;
        }
    }
    for (const auto& sequence : target.sequences) {
        if (query_map.find(sequence.name) == query_map.end()) {
            missing.push_back(sequence.name);
            identical = false;
            target_subset = false;
        }
    }
    if (use_md5 && query.md5s_present && sequence_mismatch) {
        return {target.name, "NOT_COMPATIBLE", "Status: [DIFFER_IN_SEQUENCE]. Run CompareReferences tool for more information on reference differences."};
    }
    if (identical && names_lengths_match) {
        if (use_md5 && query.md5s_present)
            return {target.name, "COMPATIBLE", "The sequence dictionaries exactly match"};
        return {target.name, "COMPATIBLE", "All sequence names and lengths match in the sequence dictionaries. Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references."};
    }
    if (query_subset && !missing.empty() && names_lengths_match) {
        std::ostringstream summary;
        summary << "The sequence dictionary in " << query_name << " is a subset of the " << target.name
                << " reference sequence dictionary. Missing sequence(s): [";
        for (std::size_t index = 0; index < missing.size(); ++index) {
            if (index) summary << ", ";
            summary << missing[index];
        }
        summary << ']';
        if (!use_md5 || !query.md5s_present)
            summary << ". Since the MD5s are lacking, we can't confirm there aren't mismatching bases in the references.";
        return {target.name, "COMPATIBLE_SUBSET", summary.str()};
    }
    const bool query_is_superset = target_subset && target.sequences.size() < query.sequences.size();
    const auto status = query_is_superset ? "SUPERSET" : "DIFFER_IN_SEQUENCES_PRESENT";
    return {target.name, "NOT_COMPATIBLE", "Status: [" + std::string(status) + "]. Run CompareReferences tool for more information on reference differences."};
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

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    const auto query = load_query(options);
    const auto query_path = options.input.empty() ? options.variant : options.input;
    const auto query_name = fastgatk::reference::basename(query_path);
    std::vector<Compatibility> results;
    for (const auto& path : options.references) {
        const auto reference = fastgatk::reference::load_fasta_info(path, "RECALCULATE_IF_MISSING");
        results.push_back(compare(query, reference, query_name, query.md5s_present));
    }
    std::ostringstream table;
    table << "#Current Reference: " << query_name << '\n'
          << "Reference\tCompatibility\tSummary\n";
    for (const auto& result : results)
        table << result.reference << '\t' << result.status << '\t' << result.summary << '\n';
    if (options.output.empty()) std::cout << table.str();
    else {
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compatibility table: " + options.output);
        output << table.str();
    }
    if (!options.manifest.empty()) {
        std::ofstream manifest(options.manifest);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
        manifest << "{\"schema_version\":1,\"tool\":\"CheckReferenceCompatibility\",\"implementation\":\"fastgatk-check-reference-compatibility\",\"status\":\"contract-compatible\",\"query\":\""
                 << json_escape(query_path) << "\",\"query_md5s_present\":" << (query.md5s_present ? "true" : "false")
                 << ",\"references\":" << results.size() << ",\"output\":\"" << json_escape(options.output)
                 << "\",\"telemetry\":{\"host_reference_io\":true,\"dictionary_records\":" << query.sequences.size() << "}}\n";
    }
    std::cerr << "{\"tool\":\"CheckReferenceCompatibility\",\"status\":\"contract-compatible\",\"references\":"
              << results.size() << ",\"query_md5s_present\":" << (query.md5s_present ? "true" : "false") << "}\n";
    return 0;
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
