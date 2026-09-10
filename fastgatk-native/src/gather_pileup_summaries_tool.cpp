#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/io/java_numeric.hpp"

#include <array>
#include <algorithm>
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
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#ifndef FASTGATK_HAS_ZLIB
#define FASTGATK_HAS_ZLIB 0
#endif

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string sequence_dictionary;
    std::string output;
    std::string manifest;
    bool allow_overlaps = false;
    bool reject_overlaps = false;
    bool disable_sequence_dictionary_validation = false;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
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

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-gather-pileup-summaries (GATK-compatible native prototype)\n"
                         "  -I/--I, --input FILE              pileup table (repeatable)\n"
                         "  -SD, --sequence-dictionary FILE  dictionary or FASTA/FAI\n"
                         "  -O/--O, --output FILE             gathered pileup table\n"
                         "      --allow-overlaps             retain duplicate loci (GATK default)\n"
                         "      --reject-overlaps            fail on duplicate loci (strict extension)\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-I" || is_option(argument, "--I") || is_option(argument, "--input")) {
            const char* name = argument.rfind("--input", 0) == 0 ? "--input" : "--I";
            options.inputs.push_back(require_value(index, argc, argv, argument, name, "-I"));
        }
        else if (argument == "-SD" || is_option(argument, "--sequence-dictionary"))
            options.sequence_dictionary = require_value(index, argc, argv, argument, "--sequence-dictionary", "-SD");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.sequence_dictionary = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-O" || is_option(argument, "--O") || is_option(argument, "--output")) {
            const char* name = argument.rfind("--output", 0) == 0 ? "--output" : "--O";
            options.output = require_value(index, argc, argv, argument, name, "-O");
        }
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--allow-overlaps")
            options.allow_overlaps = true;
        else if (argument == "--reject-overlaps")
            options.reject_overlaps = true;
        else if (argument == "--quiet") {
            // Progress/logging compatibility flag.
        } else if (argument == "--disable-sequence-dictionary-validation") {
            options.disable_sequence_dictionary_validation = true;
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("at least one -I/--input is required");
    if (options.sequence_dictionary.empty()) throw std::invalid_argument("-SD/--sequence-dictionary is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
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

std::string format_double(double value) {
    return fastgatk::io::java_double(value);
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

struct Site {
    std::string contig;
    int position = 0;
    int ref_count = 0;
    int alt_count = 0;
    int other_count = 0;
    double allele_frequency = 0.0;
    std::size_t input_file = 0;

    int total() const { return ref_count + alt_count + other_count; }
};

struct Table {
    std::string sample = "UNKNOWN";
    std::vector<Site> sites;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const auto end = line.find('\t', start);
        fields.push_back(line.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return fields;
}

bool compressed_path(const std::string& path) {
    if (path.size() < 3) return false;
    const auto lower = [](unsigned char value) {
        return static_cast<char>(std::tolower(value));
    };
    const auto offset = path.size() - 3;
    return lower(static_cast<unsigned char>(path[offset])) == '.' &&
           lower(static_cast<unsigned char>(path[offset + 1])) == 'g' &&
           lower(static_cast<unsigned char>(path[offset + 2])) == 'z';
}

template<class Callback>
void for_each_text_line(const std::string& path, Callback&& callback) {
    if (!compressed_path(path)) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open pileup table: " + path);
        for (std::string line; std::getline(input, line);)
            callback(std::move(line));
        if (input.bad()) throw std::runtime_error("BAD_INPUT: failed reading pileup table: " + path);
        return;
    }
#if FASTGATK_HAS_ZLIB
    gzFile input = gzopen(path.c_str(), "rb");
    if (input == nullptr) throw std::runtime_error("BAD_INPUT: cannot open compressed pileup table: " + path);
    try {
        char buffer[64 * 1024];
        while (true) {
            std::string line;
            while (true) {
                char* result = gzgets(input, buffer, static_cast<int>(sizeof(buffer)));
                if (result == nullptr) break;
                const auto length = std::strlen(buffer);
                line.append(buffer, length);
                if (length == 0 || buffer[length - 1] == '\n') break;
            }
            if (line.empty()) {
                int error = Z_OK;
                (void)gzerror(input, &error);
                if (error != Z_OK && error != Z_STREAM_END)
                    throw std::runtime_error("BAD_INPUT: failed reading compressed pileup table: " + path);
                break;
            }
            if (!line.empty() && line.back() == '\n') line.pop_back();
            if (!line.empty() && line.back() == '\r') line.pop_back();
            callback(std::move(line));
        }
    } catch (...) {
        gzclose(input);
        throw;
    }
    int close_status = gzclose(input);
    if (close_status != Z_OK)
        throw std::runtime_error("BAD_INPUT: failed closing compressed pileup table: " + path);
#else
    throw std::runtime_error("BACKEND_UNAVAILABLE: zlib is required for compressed pileup table: " + path);
#endif
}

class TextWriter {
public:
    explicit TextWriter(const std::string& path) : path_(path), compressed_(compressed_path(path)) {
        if (!compressed_) {
            plain_.open(path);
            if (!plain_) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open output: " + path);
        } else {
#if FASTGATK_HAS_ZLIB
            compressed_file_ = gzopen(path.c_str(), "wb");
            if (compressed_file_ == nullptr)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open compressed output: " + path);
#else
            throw std::runtime_error("BACKEND_UNAVAILABLE: zlib is required for compressed output: " + path);
#endif
        }
    }

    TextWriter(const TextWriter&) = delete;
    TextWriter& operator=(const TextWriter&) = delete;

    ~TextWriter() {
#if FASTGATK_HAS_ZLIB
        if (compressed_file_ != nullptr) gzclose(compressed_file_);
#endif
    }

    void write(const std::string& text) {
        if (!compressed_) {
            plain_ << text;
            if (!plain_) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write output: " + path_);
            return;
        }
#if FASTGATK_HAS_ZLIB
        const auto expected = static_cast<unsigned int>(text.size());
        if (expected != 0 && gzwrite(compressed_file_, text.data(), expected) != static_cast<int>(expected))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed output: " + path_);
#endif
    }

    void close() {
        if (!compressed_) {
            plain_.flush();
            if (!plain_) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize output: " + path_);
            plain_.close();
            return;
        }
#if FASTGATK_HAS_ZLIB
        if (compressed_file_ != nullptr) {
            const auto status = gzclose(compressed_file_);
            compressed_file_ = nullptr;
            if (status != Z_OK)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize compressed output: " + path_);
        }
#endif
    }

private:
    std::string path_;
    bool compressed_ = false;
    std::ofstream plain_;
#if FASTGATK_HAS_ZLIB
    gzFile compressed_file_ = nullptr;
#endif
};

Table read_table(const std::string& path, std::size_t input_file) {
    Table table;
    bool header_seen = false;
    std::size_t line_number = 0;
    for_each_text_line(path, [&](std::string line) {
        ++line_number;
        if (line.empty()) return;
        if (line.rfind("#<METADATA>SAMPLE=", 0) == 0) {
            table.sample = line.substr(std::string("#<METADATA>SAMPLE=").size());
            return;
        }
        if (line[0] == '#') return;
        const auto fields = split_tab(line);
        if (!header_seen) {
            if (fields.size() != 6 || fields[0] != "contig" || fields[1] != "position")
                throw std::runtime_error("BAD_INPUT: invalid pileup table header at " + path + ":" + std::to_string(line_number));
            header_seen = true;
            return;
        }
        if (fields.size() != 6)
            throw std::runtime_error("BAD_INPUT: invalid pileup table row at " + path + ":" + std::to_string(line_number));
        Site site;
        site.contig = fields[0];
        site.position = std::stoi(fields[1]);
        site.ref_count = std::stoi(fields[2]);
        site.alt_count = std::stoi(fields[3]);
        site.other_count = std::stoi(fields[4]);
        site.allele_frequency = std::stod(fields[5]);
        site.input_file = input_file;
        if (site.contig.empty() || site.position < 1 || site.ref_count < 0 || site.alt_count < 0 ||
            site.other_count < 0 || !std::isfinite(site.allele_frequency) ||
            site.allele_frequency < 0.0 || site.allele_frequency > 1.0)
            throw std::runtime_error("BAD_INPUT: invalid pileup table values at " + path + ":" + std::to_string(line_number));
        table.sites.push_back(std::move(site));
    });
    if (!header_seen) throw std::runtime_error("BAD_INPUT: pileup table has no header: " + path);
    return table;
}

struct DictionaryEntry {
    std::string name;
    std::int64_t length = -1;
};

std::vector<DictionaryEntry> dictionary_entries(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        const std::string fai = path + ".fai";
        input.open(fai);
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open sequence dictionary/FAI: " + path);
    }
    std::vector<std::string> lines;
    for (std::string line; std::getline(input, line);) lines.push_back(std::move(line));
    const bool fasta = std::any_of(lines.begin(), lines.end(), [](const std::string& line) {
        return !line.empty() && line[0] == '>';
    });
    std::vector<DictionaryEntry> contigs;
    if (fasta) {
        DictionaryEntry current;
        for (const auto& line : lines) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!current.name.empty()) contigs.push_back(std::move(current));
                current = DictionaryEntry{};
                const auto end = line.find_first_of(" \t", 1);
                current.name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
                current.length = 0;
            } else if (!current.name.empty()) {
                current.length += static_cast<std::int64_t>(
                    std::count_if(line.begin(), line.end(), [](const char c) {
                        return !std::isspace(static_cast<unsigned char>(c));
                    }));
            }
        }
        if (!current.name.empty()) contigs.push_back(std::move(current));
        if (contigs.empty()) throw std::runtime_error("BAD_INPUT: sequence dictionary has no contigs: " + path);
        return contigs;
    }
    for (const auto& line : lines) {
        if (line.empty()) continue;
        if (line.rfind("@SQ", 0) == 0) {
            DictionaryEntry entry;
            for (const auto& field : split_tab(line)) {
                if (field.rfind("SN:", 0) == 0) entry.name = field.substr(3);
                else if (field.rfind("LN:", 0) == 0) entry.length = std::stoll(field.substr(3));
            }
            if (!entry.name.empty()) contigs.push_back(std::move(entry));
        } else if (line[0] != '@') {
            const auto fields = split_tab(line);
            if (!fields.empty() && !fields[0].empty()) {
                const auto length = fields.size() > 1 ? std::stoll(fields[1]) : -1;
                contigs.push_back(DictionaryEntry{fields[0], length});
            }
        }
    }
    if (contigs.empty()) throw std::runtime_error("BAD_INPUT: sequence dictionary has no contigs: " + path);
    return contigs;
}

std::string site_key(const Site& site) {
    return site.contig + ":" + std::to_string(site.position);
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto dictionary = dictionary_entries(options.sequence_dictionary);
    std::unordered_map<std::string, int> contig_order;
    std::unordered_map<std::string, std::int64_t> contig_lengths;
    contig_order.reserve(dictionary.size() * 2 + 1);
    contig_lengths.reserve(dictionary.size() * 2 + 1);
    for (std::size_t index = 0; index < dictionary.size(); ++index) {
        contig_order.emplace(dictionary[index].name, static_cast<int>(index));
        contig_lengths.emplace(dictionary[index].name, dictionary[index].length);
    }

    std::vector<Table> tables;
    std::size_t all_site_count = 0;
    std::string sample = "UNKNOWN";
    bool nonempty_sample_seen = false;
    std::size_t empty_inputs = 0;
    std::unordered_set<std::string> unknown_contigs;
    for (std::size_t input_index = 0; input_index < options.inputs.size(); ++input_index) {
        auto table = read_table(options.inputs[input_index], input_index);
        if (table.sites.empty()) ++empty_inputs;
        // GATK drops empty PileupSummaryTable inputs before it validates the
        // sample metadata.  Consequently empty shards may carry different
        // SAMPLE values without affecting a non-empty gather.  Validate only
        // tables that contribute rows, and select the sample from the first
        // contributing table (rather than from an empty shard encountered
        // earlier in the input list).
        if (!table.sites.empty()) {
            if (!nonempty_sample_seen) {
                sample = table.sample;
                nonempty_sample_seen = true;
            } else if (table.sample != "UNKNOWN" && table.sample != sample)
                throw std::runtime_error("BAD_INPUT: pileup tables contain different samples: " + sample + " vs " + table.sample);
        }
        for (const auto& site : table.sites) {
            const auto order = contig_order.find(site.contig);
            if (order == contig_order.end()) {
                if (!options.disable_sequence_dictionary_validation)
                    throw std::runtime_error("BAD_INPUT: contig is absent from sequence dictionary: " + site.contig);
                unknown_contigs.insert(site.contig);
            } else {
                const auto length = contig_lengths.at(site.contig);
                if (!options.disable_sequence_dictionary_validation && length > 0 &&
                    static_cast<std::int64_t>(site.position) > length)
                    throw std::runtime_error("BAD_INPUT: pileup position exceeds sequence dictionary length for contig " +
                                             site.contig);
            }
            ++all_site_count;
        }
        if (!table.sites.empty()) tables.push_back(std::move(table));
    }
    if (!unknown_contigs.empty()) {
        std::vector<std::string> sorted_unknown(unknown_contigs.begin(), unknown_contigs.end());
        std::sort(sorted_unknown.begin(), sorted_unknown.end());
        for (const auto& contig : sorted_unknown) {
            const auto order = static_cast<int>(contig_order.size());
            contig_order.emplace(contig, order);
            contig_lengths.emplace(contig, -1);
        }
    }
    if (sample.empty()) sample = "UNKNOWN";
    // GATK's GatherPileupSummaries sorts non-empty input files by the first
    // pileup record only, then writes each file's rows in their original
    // order.  It does not globally sort loci or deduplicate overlaps.  Keep
    // that file-level ordering so this native replacement has identical
    // behavior for scattered tables that are already internally sorted.
    std::stable_sort(tables.begin(), tables.end(), [&](const Table& left, const Table& right) {
        const auto& left_site = left.sites.front();
        const auto& right_site = right.sites.front();
        const int left_contig = contig_order.at(left_site.contig);
        const int right_contig = contig_order.at(right_site.contig);
        if (left_contig != right_contig) return left_contig < right_contig;
        return left_site.position < right_site.position;
    });

    std::unordered_set<std::string> seen;
    seen.reserve(all_site_count * 2 + 1);
    std::size_t duplicate_sites = 0;
    std::vector<Site> gathered;
    gathered.reserve(all_site_count);
    for (const auto& table : tables) {
        for (const auto& site : table.sites) {
            if (!seen.insert(site_key(site)).second) {
                ++duplicate_sites;
                if (options.reject_overlaps)
                    throw std::runtime_error("BAD_INPUT: overlapping pileup summaries at " + site_key(site));
            }
            gathered.push_back(site);
        }
    }
    TextWriter output(options.output);
    // An all-empty gather is a valid GATK result, but its writer emits only
    // the column header and no SAMPLE metadata line.
    if (!gathered.empty()) output.write("#<METADATA>SAMPLE=" + sample + "\n");
    output.write("contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n");
    for (const auto& site : gathered) {
        std::ostringstream row;
        row << site.contig << '\t' << site.position << '\t'
            << site.ref_count << '\t' << site.alt_count << '\t' << site.other_count
            << '\t' << format_double(site.allele_frequency) << '\n';
        output.write(row.str());
    }
    output.close();
    if (!file_complete(options.output))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: gathered table is incomplete");
    const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create manifest: " + manifest_path);
    manifest << "{\"schema_version\":1,\"tool\":\"GatherPileupSummaries\","
             << "\"implementation\":\"fastgatk-gather-pileup-summaries\",\"status\":\"prototype\","
             << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"table\","
             << "\"compatibility\":{\"dictionary_order\":true,\"sample_validation\":true,"
             << "\"file_first_record_order\":true,\"overlap_detection\":"
             << (options.reject_overlaps ? "true" : "false")
             << ",\"reject_overlaps\":" << (options.reject_overlaps ? "true" : "false")
             << ",\"allow_overlaps\":true"
             << ",\"compressed_input_count\":" << std::count_if(
                    options.inputs.begin(), options.inputs.end(), compressed_path)
             << ",\"compressed_output\":" << (compressed_path(options.output) ? "true" : "false")
             << ",\"sequence_dictionary_validation\":"
             << (!options.disable_sequence_dictionary_validation ? "true" : "false")
             << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
             << json_escape(options.output) << "\",\"kind\":\"table\",\"complete\":true}],"
             << "\"telemetry\":{\"resources\":" << resources.to_json()
             << ",\"input_files\":" << options.inputs.size() << ",\"empty_inputs\":" << empty_inputs
             << ",\"input_sites\":" << all_site_count << ",\"output_sites\":" << gathered.size()
             << ",\"duplicate_sites\":" << duplicate_sites << ",\"sample\":\"" << json_escape(sample) << "\"}}\n";
    std::cout << "{\"tool\":\"GatherPileupSummaries\",\"status\":\"prototype\","
              << "\"input_files\":" << options.inputs.size() << ",\"output_sites\":" << gathered.size()
              << ",\"duplicate_sites\":" << duplicate_sites << "}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-gather-pileup-summaries: " << error.what() << '\n';
        return 2;
    }
}
