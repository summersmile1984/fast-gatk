#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/core/plan.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string output;
    std::string manifest;
    std::string sample;
    double convergence_threshold = 1e-4;
    int max_iterations = 20;
    int max_depth = 200;
    int threads = 1;
    std::string tmp_dir;
    std::string verbosity = "INFO";
    bool quiet = false;
};

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

struct OrientationRecord {
    std::uint8_t alt = 0;
    std::uint64_t f1r2 = 0;
    std::uint64_t r1f2 = 0;
};

// The standard GATK CollectF1R2Counts hand-off contains three metrics/table
// members.  Keep this representation flat and value-owned at the Host
// boundary; the EM below only consumes integer counts and fixed-size arrays.
struct AltSiteRecord {
    std::string context;
    int ref_count = 0;
    int alt_count = 0;
    int ref_f1r2 = 0;
    int alt_f1r2 = 0;
    char alt = 'N';
};

struct StandardOrientationInput {
    std::unordered_map<std::string, std::vector<std::uint64_t>> ref_hist;
    std::unordered_map<std::string, std::vector<std::uint64_t>> alt_hist;
    std::vector<AltSiteRecord> alt_sites;
    std::string sample;
    std::size_t ref_records = 0;
    std::size_t alt_records = 0;
};

struct ParsedOrientationInput {
    bool standard = false;
    std::vector<OrientationRecord> legacy_records;
    StandardOrientationInput standard_input;
    // A single CollectF1R2Counts archive may contain one ref/alt histogram
    // triplet for every read-group sample.  GATK preserves those sample
    // boundaries and LearnReadOrientationModel writes one prior member per
    // sample; do not key the standard tables only by context.
    std::map<std::string, StandardOrientationInput> standard_inputs;
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
            std::cout << "fastgatk-learn-read-orientation-model (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 Mutect2 F1R2 TSV sidecar (repeatable)\n"
                         "  -O, --output FILE                artifact-prior .tar.gz\n"
                         "      --sample NAME                sample metadata (default UNKNOWN)\n"
                         "      --convergence-threshold F    EM compatibility option\n"
                         "      --num-em-iterations N        EM compatibility option\n"
                         "      --max-depth N                histogram compatibility option\n"
                         "      --threads N                  Kokkos execution threads\n"
                         "      --QUIET[=true|false]         suppress job-summary logging\n"
                         "      --tmp-dir DIR                temporary directory\n"
                         "      --verbosity LEVEL            ERROR/WARNING/INFO/DEBUG\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || is_option(argument, "--input"))
            options.inputs.push_back(require_value(index, argc, argv, argument, "--input", "-I"));
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--sample"))
            options.sample = require_value(index, argc, argv, argument, "--sample");
        else if (is_option(argument, "--convergence-threshold"))
            options.convergence_threshold = std::stod(require_value(
                index, argc, argv, argument, "--convergence-threshold"));
        else if (is_option(argument, "--num-em-iterations") || is_option(argument, "--max-em-iterations"))
            options.max_iterations = std::stoi(require_value(
                index, argc, argv, argument,
                argument.rfind("--max-em-iterations", 0) == 0
                    ? "--max-em-iterations" : "--num-em-iterations"));
        else if (is_option(argument, "--max-depth"))
            options.max_depth = std::stoi(require_value(index, argc, argv, argument, "--max-depth"));
        else if (is_option(argument, "--threads"))
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--QUIET" || argument == "--quiet" ||
                 argument.rfind("--QUIET=", 0) == 0 || argument.rfind("--quiet=", 0) == 0)
            options.quiet = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, argument.rfind("--QUIET", 0) == 0 ? "--QUIET" : "--quiet");
        else if (argument == "--tmp-dir" || is_option(argument, "--tmp-dir"))
            options.tmp_dir = require_value(index, argc, argv, argument, "--tmp-dir");
        else if (argument == "--use-jdk-deflater" || argument == "-jdk-deflater" ||
                 argument == "--use-jdk-inflater" || argument == "-jdk-inflater" ||
                 argument.rfind("--use-jdk-deflater=", 0) == 0 || argument.rfind("--use-jdk-inflater=", 0) == 0) {
            const bool deflater = argument.rfind("--use-jdk-deflater", 0) == 0 || argument == "-jdk-deflater";
            (void)fastgatk::native::parse_optional_boolean(index, argc, argv, argument,
                deflater ? (argument == "-jdk-deflater" ? "-jdk-deflater" : "--use-jdk-deflater")
                         : (argument == "-jdk-inflater" ? "-jdk-inflater" : "--use-jdk-inflater"));
        } else if (argument == "--disable-sequence-dictionary-validation") {
            // Accepted launcher compatibility flag.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--VERBOSITY")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--VERBOSITY", 0) == 0 ? "--VERBOSITY" : "--verbosity");
            if (argument == "--verbosity" || argument.rfind("--verbosity=", 0) == 0 ||
                argument == "--VERBOSITY" || argument.rfind("--VERBOSITY=", 0) == 0)
                options.verbosity = argument.find('=') == std::string::npos ? argv[index] :
                    option_value(argument, argument.rfind("--VERBOSITY", 0) == 0 ? "--VERBOSITY" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("at least one -I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.output.size() < 7 || options.output.substr(options.output.size() - 7) != ".tar.gz")
        throw std::invalid_argument("-O/--output must end in .tar.gz");
    if (options.sample.empty()) options.sample = "UNKNOWN";
    if (!std::isfinite(options.convergence_threshold) || options.convergence_threshold <= 0.0 ||
        options.max_iterations < 1 || options.max_depth < 1 || options.threads < 1)
        throw std::invalid_argument("orientation model parameters must be positive");
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

std::string trim_left(std::string value) {
    const auto first = value.find_first_not_of(" \t");
    return first == std::string::npos ? std::string{} : value.substr(first);
}

int base_index(const std::string& value) {
    if (value.size() != 1) return -1;
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(value[0])))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

std::vector<OrientationRecord> read_input_stream(std::istream& input, const std::string& path,
                                                 std::string& sample) {
    std::vector<OrientationRecord> records;
    bool header_seen = false;
    std::size_t line_number = 0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty()) continue;
        if (line.rfind("#<METADATA>SAMPLE=", 0) == 0) {
            const auto value = line.substr(std::string("#<METADATA>SAMPLE=").size());
            if (sample == "UNKNOWN") sample = value;
            else if (sample != value && value != "UNKNOWN")
                throw std::runtime_error("BAD_INPUT: F1R2 inputs have different samples");
            continue;
        }
        std::string data = line;
        if (!data.empty() && data[0] == '#') data = trim_left(data.substr(1));
        const auto fields = split_tab(data);
        if (!header_seen) {
            if (fields.size() >= 5 && fields[0] == "contig" && fields[1] == "position" &&
                fields[2] == "alt" && fields[3] == "f1r2" && fields[4] == "r1f2") {
                header_seen = true;
                continue;
            }
            if (line[0] == '#') continue;
            throw std::runtime_error("BAD_INPUT: missing F1R2 TSV header at " + path + ":" +
                                     std::to_string(line_number));
        }
        if (fields.size() < 5)
            throw std::runtime_error("BAD_INPUT: malformed F1R2 row at " + path + ":" +
                                     std::to_string(line_number));
        const int alt = base_index(fields[2]);
        if (alt < 0) throw std::runtime_error("BAD_INPUT: ALT must be A/C/G/T at " + path + ":" +
                                             std::to_string(line_number));
        const auto f1r2 = std::stoull(fields[3]);
        const auto r1f2 = std::stoull(fields[4]);
        if (fields.size() >= 6 && std::stoull(fields[5]) != f1r2 + r1f2)
            throw std::runtime_error("BAD_INPUT: total does not match F1R2/R1F2 at " + path + ":" +
                                     std::to_string(line_number));
        records.push_back(OrientationRecord{static_cast<std::uint8_t>(alt), f1r2, r1f2});
    }
    if (!header_seen) throw std::runtime_error("BAD_INPUT: F1R2 TSV has no header: " + path);
    return records;
}

std::vector<OrientationRecord> read_input(const std::string& path, std::string& sample) {
#if FASTGATK_HAS_ZLIB
    if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
        gzFile compressed = gzopen(path.c_str(), "rb");
        if (!compressed) throw std::runtime_error("BAD_INPUT: cannot open compressed F1R2 TSV: " + path);
        std::string text;
        std::array<char, 1 << 16> buffer{};
        int read = 0;
        while ((read = gzread(compressed, buffer.data(), static_cast<unsigned>(buffer.size()))) > 0)
            text.append(buffer.data(), static_cast<std::size_t>(read));
        const auto close_status = gzclose(compressed);
        if (read < 0 || close_status != Z_OK)
            throw std::runtime_error("BAD_INPUT: cannot read compressed F1R2 TSV: " + path);
        std::istringstream input(text);
        return read_input_stream(input, path, sample);
    }
#else
    if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz")
        throw std::runtime_error("BACKEND_UNAVAILABLE: compressed F1R2 TSV requires zlib: " + path);
#endif
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open F1R2 TSV: " + path);
    return read_input_stream(input, path, sample);
}

#if FASTGATK_HAS_ZLIB
std::uint64_t tar_octal(const char* field, std::size_t length) {
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < length && field[index] != '\0'; ++index) {
        if (field[index] >= '0' && field[index] <= '7') value = value * 8 + (field[index] - '0');
    }
    return value;
}

std::string read_gzip_file(const std::string& path) {
    gzFile stream = gzopen(path.c_str(), "rb");
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open F1R2 tar.gz: " + path);
    std::string bytes;
    std::array<char, 1 << 16> buffer{};
    int read = 0;
    while ((read = gzread(stream, buffer.data(), static_cast<unsigned>(buffer.size()))) > 0)
        bytes.append(buffer.data(), static_cast<std::size_t>(read));
    const int close_status = gzclose(stream);
    if (read < 0 || close_status != Z_OK)
        throw std::runtime_error("BAD_INPUT: cannot decompress F1R2 tar.gz: " + path);
    return bytes;
}

std::string tar_member_name(const unsigned char* header) {
    std::size_t length = 0;
    while (length < 100 && header[length] != 0) ++length;
    return std::string(reinterpret_cast<const char*>(header), length);
}

std::uint64_t parse_metric_count(const std::string& value, const std::string& path) {
    try {
        std::size_t consumed = 0;
        const auto result = std::stoull(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("trailing characters");
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("BAD_INPUT: invalid histogram count in " + path + ": " + value);
    }
}

void parse_standard_metrics(const std::string& content, const std::string& path,
                            bool reference_histogram, std::size_t max_depth,
                            StandardOrientationInput& result) {
    std::vector<std::string> columns;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty()) continue;
        if (line.rfind("# ", 0) == 0 && result.sample.empty()) {
            const auto candidate = trim_left(line.substr(2));
            if (!candidate.empty()) result.sample = candidate;
        }
        const auto fields = split_tab(line);
        if (!fields.empty() && fields[0] == "depth") {
            columns = fields;
            break;
        }
    }
    if (columns.size() < 2)
        throw std::runtime_error("BAD_INPUT: GATK histogram header is missing in " + path);
    auto& destination = reference_histogram ? result.ref_hist : result.alt_hist;
    for (std::size_t column = 1; column < columns.size(); ++column) {
        if (columns[column].empty()) continue;
        // `--max-depth` is an EM observation bound in LearnReadOrientationModel,
        // not a request to re-bin an already materialized CollectF1R2Counts
        // histogram.  GATK retains source bins above the bound for
        // numExamples/numAltExamples, while the engine only visits depths
        // 1..maxDepth during each E/M step.  Seed the requested range here,
        // then grow for any higher bins encountered in the input.
        destination.try_emplace(columns[column], max_depth + 1, 0);
    }
    bool data_started = false;
    stream.clear();
    stream.seekg(0);
    while (std::getline(stream, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.empty() || fields[0] == "depth") { data_started = true; continue; }
        if (!data_started || columns.empty() || fields.size() < columns.size()) continue;
        bool numeric = !fields[0].empty() && std::all_of(fields[0].begin(), fields[0].end(),
                                                         [](const char c) { return std::isdigit(static_cast<unsigned char>(c)); });
        if (!numeric) continue;
        const auto raw_depth = parse_metric_count(fields[0], path);
        if (raw_depth > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() - 1))
            throw std::runtime_error("BAD_INPUT: histogram depth is too large in " + path);
        const auto depth = static_cast<std::size_t>(raw_depth);
        for (std::size_t column = 1; column < columns.size(); ++column) {
            const auto count = parse_metric_count(fields[column], path);
            if (count == 0 || columns[column].empty()) continue;
            auto& bins = destination[columns[column]];
            if (bins.size() <= depth) bins.resize(depth + 1, 0);
            bins[depth] += count;
            if (reference_histogram) result.ref_records += count;
            else result.alt_records += count;
        }
    }
}

void parse_standard_alt_table(const std::string& content, const std::string& path,
                              StandardOrientationInput& result) {
    std::vector<std::string> columns;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty()) continue;
        if (line.rfind("#<METADATA>SAMPLE=", 0) == 0) {
            const auto sample = line.substr(std::string("#<METADATA>SAMPLE=").size());
            if (result.sample.empty() || result.sample == "UNKNOWN") result.sample = sample;
            else if (result.sample != sample && sample != "UNKNOWN")
                throw std::runtime_error("BAD_INPUT: F1R2 table samples disagree");
            continue;
        }
        if (line[0] == '#') continue;
        columns = split_tab(line);
        break;
    }
    const std::array<std::string, 7> expected{
        "context", "ref_count", "alt_count", "ref_f1r2", "alt_f1r2", "depth", "alt"};
    if (columns.size() < expected.size())
        throw std::runtime_error("BAD_INPUT: GATK alt table header is missing in " + path);
    for (std::size_t i = 0; i < expected.size(); ++i)
        if (columns[i] != expected[i])
            throw std::runtime_error("BAD_INPUT: unexpected GATK alt table column in " + path);

    stream.clear();
    stream.seekg(0);
    while (std::getline(stream, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() < expected.size() || fields[0] == "context") continue;
        const int alt = base_index(fields[6]);
        if (fields[0].size() != 3 || alt < 0)
            throw std::runtime_error("BAD_INPUT: malformed GATK alt table row in " + path);
        try {
            AltSiteRecord record;
            record.context = fields[0];
            record.ref_count = std::stoi(fields[1]);
            record.alt_count = std::stoi(fields[2]);
            record.ref_f1r2 = std::stoi(fields[3]);
            record.alt_f1r2 = std::stoi(fields[4]);
            const int depth = std::stoi(fields[5]);
            record.alt = static_cast<char>(std::toupper(static_cast<unsigned char>(fields[6][0])));
            if (record.ref_count < 0 || record.alt_count <= 0 || record.ref_f1r2 < 0 ||
                record.alt_f1r2 < 0 || depth != record.ref_count + record.alt_count ||
                record.ref_f1r2 > record.ref_count || record.alt_f1r2 > record.alt_count)
                throw std::invalid_argument("invalid count relationship");
            result.alt_sites.push_back(std::move(record));
        } catch (const std::exception&) {
            throw std::runtime_error("BAD_INPUT: malformed GATK alt table row in " + path);
        }
    }
}

void merge_standard_orientation_input(StandardOrientationInput& destination,
                                      const StandardOrientationInput& source) {
    if (destination.sample.empty()) destination.sample = source.sample;
    for (const auto& [label, counts] : source.ref_hist) {
        auto& output = destination.ref_hist[label];
        if (output.size() < counts.size()) output.resize(counts.size(), 0);
        for (std::size_t index = 0; index < counts.size(); ++index) output[index] += counts[index];
    }
    for (const auto& [label, counts] : source.alt_hist) {
        auto& output = destination.alt_hist[label];
        if (output.size() < counts.size()) output.resize(counts.size(), 0);
        for (std::size_t index = 0; index < counts.size(); ++index) output[index] += counts[index];
    }
    destination.alt_sites.insert(destination.alt_sites.end(), source.alt_sites.begin(), source.alt_sites.end());
    destination.ref_records += source.ref_records;
    destination.alt_records += source.alt_records;
}

ParsedOrientationInput read_standard_tar_input(const std::string& path, std::string& sample,
                                               std::size_t max_depth) {
    const auto archive = read_gzip_file(path);
    ParsedOrientationInput parsed;
    std::size_t offset = 0;
    while (offset + 512 <= archive.size()) {
        const auto* header = reinterpret_cast<const unsigned char*>(archive.data() + offset);
        bool empty = true;
        for (std::size_t index = 0; index < 512; ++index)
            if (header[index] != 0) { empty = false; break; }
        if (empty) break;
        const std::string member = tar_member_name(header);
        const auto size = tar_octal(reinterpret_cast<const char*>(header + 124), 12);
        const char type = static_cast<char>(header[156]);
        offset += 512;
        if (offset + size > archive.size())
            throw std::runtime_error("BAD_INPUT: truncated F1R2 tar member: " + member);
        if (type == '0' || type == '\0') {
            const std::string content = archive.substr(offset, size);
            if (member.size() >= std::string(".ref_histogram").size() &&
                member.ends_with(".ref_histogram")) {
                parsed.standard = true;
                StandardOrientationInput member_input;
                parse_standard_metrics(content, member, true, max_depth, member_input);
                const auto member_sample = member_input.sample.empty() ? "UNKNOWN" : member_input.sample;
                merge_standard_orientation_input(parsed.standard_inputs[member_sample], member_input);
            } else if (member.size() >= std::string(".alt_histogram").size() &&
                       member.ends_with(".alt_histogram")) {
                parsed.standard = true;
                StandardOrientationInput member_input;
                parse_standard_metrics(content, member, false, max_depth, member_input);
                const auto member_sample = member_input.sample.empty() ? "UNKNOWN" : member_input.sample;
                merge_standard_orientation_input(parsed.standard_inputs[member_sample], member_input);
            } else if (member.size() >= std::string(".alt_table").size() &&
                member.ends_with(".alt_table")) {
                parsed.standard = true;
                StandardOrientationInput member_input;
                parse_standard_alt_table(content, member, member_input);
                const auto member_sample = member_input.sample.empty() ? "UNKNOWN" : member_input.sample;
                merge_standard_orientation_input(parsed.standard_inputs[member_sample], member_input);
            } else {
                std::istringstream input(content);
                auto shard = read_input_stream(input, member.empty() ? path : member, sample);
                parsed.legacy_records.insert(parsed.legacy_records.end(), shard.begin(), shard.end());
            }
        }
        offset += ((size + 511) / 512) * 512;
    }
    if (parsed.standard) {
        for (auto& [member_sample, input] : parsed.standard_inputs) {
            if (input.sample.empty()) input.sample = member_sample;
            if (member_sample == "UNKNOWN" && sample != "UNKNOWN") input.sample = sample;
            else if (sample != "UNKNOWN" && member_sample != "UNKNOWN" && member_sample != sample)
                throw std::runtime_error("BAD_INPUT: F1R2 tar samples disagree");
        }
        if (sample == "UNKNOWN" && parsed.standard_inputs.size() == 1) {
            const auto& input = parsed.standard_inputs.begin()->second;
            if (!input.sample.empty() && input.sample != "UNKNOWN") sample = input.sample;
        }
        return parsed;
    }
    if (parsed.legacy_records.empty())
        throw std::runtime_error("BAD_INPUT: F1R2 tar.gz has no recognized members: " + path);
    return parsed;
}

#else
ParsedOrientationInput read_standard_tar_input(const std::string& path, std::string&, std::size_t) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: zlib is required for .tar.gz orientation input: " + path);
}
#endif

std::string reverse_complement(const std::string& context) {
    std::string result(context.size(), 'N');
    for (std::size_t i = 0; i < context.size(); ++i) {
        const char base = static_cast<char>(std::toupper(static_cast<unsigned char>(context[context.size() - i - 1])));
        result[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
    }
    return result;
}

char complement_base(const char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::vector<std::string> all_contexts() {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    std::vector<std::string> contexts;
    contexts.reserve(64);
    for (const auto left : bases)
        for (const auto middle : bases)
            for (const auto right : bases)
                contexts.emplace_back(std::string{left, middle, right});
    return contexts;
}

// GATK's ArtifactPriorCollection stores the 64 contexts in a Java HashMap and
// writes map.values() directly.  That is not lexical order: with the pinned
// Java 17 runtime the table grows from 64 to 128 buckets while the contexts
// are inserted, and the final iteration order is the bucket order.  Recreate
// the small, observable part of HashMap here so a native prior table has the
// same row order without making the Host/HTSlib path depend on Java.
std::uint32_t java_string_hash(const std::string& value) {
    std::uint32_t hash = 0;
    for (const auto character : value)
        hash = hash * 31U + static_cast<unsigned char>(character);
    return hash;
}

std::vector<std::string> gatk_context_order() {
    constexpr std::size_t bucket_count = 128;
    std::array<std::vector<std::string>, bucket_count> buckets;
    for (const auto& context : all_contexts()) {
        const auto hash = java_string_hash(context);
        const auto spread = hash ^ (hash >> 16U);
        buckets[spread & (bucket_count - 1U)].push_back(context);
    }
    std::vector<std::string> result;
    result.reserve(64);
    for (const auto& bucket : buckets)
        result.insert(result.end(), bucket.begin(), bucket.end());
    return result;
}

// TableWriter ultimately calls Java Double.toString().  C++'s shortest
// round-trip conversion has the same significant digits, but chooses a
// different fixed/scientific cut-over and emits `0`/`1` without `.0`.  Keep a
// tiny formatter at this file boundary so the numerical model can stay in
// binary64 while the serialized prior table follows the GATK text contract.
std::string java_double(double value) {
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return value < 0.0 ? "-Infinity" : "Infinity";
    if (value == 0.0) return std::signbit(value) ? "-0.0" : "0.0";

    char buffer[128]{};
    const auto converted = std::to_chars(buffer, buffer + sizeof(buffer), value,
                                         std::chars_format::general);
    if (converted.ec != std::errc{})
        throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: cannot format orientation prior");
    std::string text(buffer, converted.ptr);

    // Split the shortest representation into decimal digits and a decimal
    // point position.  This lets us apply Java's exponent cut-over while
    // preserving exactly the digits selected by to_chars.
    const bool negative = !text.empty() && text.front() == '-';
    if (negative) text.erase(text.begin());
    int exponent = 0;
    const auto exponent_position = text.find_first_of("eE");
    if (exponent_position != std::string::npos) {
        exponent = std::stoi(text.substr(exponent_position + 1));
        text.erase(exponent_position);
    }
    const auto point = text.find('.');
    int decimal_position = point == std::string::npos
        ? static_cast<int>(text.size()) : static_cast<int>(point);
    std::string digits = text;
    if (point != std::string::npos) digits.erase(point, 1);
    while (!digits.empty() && digits.front() == '0') {
        digits.erase(digits.begin());
        --decimal_position;
    }
    if (digits.empty()) return "0.0";
    decimal_position += exponent;

    // Java uses fixed notation for -3 <= exponent < 7, scientific notation
    // otherwise.  (The sign is kept outside the digit stream.)
    const int scientific_exponent = decimal_position - 1;
    const int abs_exponent = scientific_exponent;
    const bool scientific = abs_exponent < -3 || abs_exponent >= 7;
    std::string result;
    if (!scientific) {
        if (negative) result.push_back('-');
        if (decimal_position <= 0) {
            result += "0.";
            result.append(static_cast<std::size_t>(-decimal_position), '0');
            result += digits;
        } else if (decimal_position >= static_cast<int>(digits.size())) {
            result += digits;
            result.append(static_cast<std::size_t>(decimal_position - digits.size()), '0');
            result += ".0";
        } else {
            result.append(digits.data(), static_cast<std::size_t>(decimal_position));
            result.push_back('.');
            result.append(digits.data() + decimal_position,
                          digits.size() - static_cast<std::size_t>(decimal_position));
        }
        if (result.find('.') == std::string::npos) result += ".0";
        return result;
    }

    if (negative) result.push_back('-');
    result.push_back(digits.front());
    if (digits.size() > 1) {
        result.push_back('.');
        result.append(digits.data() + 1, digits.size() - 1);
    } else {
        result += ".0";
    }
    result.push_back('E');
    result += scientific_exponent < 0 ? "-" : "";
    result += std::to_string(std::abs(scientific_exponent));
    return result;
}

double log_beta_binomial_probability(const int successes, const int trials,
                                     const double alpha, const double beta) {
    if (successes < 0 || trials < 0 || successes > trials) return -std::numeric_limits<double>::infinity();
    const double log_choose = std::lgamma(static_cast<double>(trials) + 1.0) -
                              std::lgamma(static_cast<double>(successes) + 1.0) -
                              std::lgamma(static_cast<double>(trials - successes) + 1.0);
    const double log_beta = std::lgamma(static_cast<double>(successes) + alpha) +
                            std::lgamma(static_cast<double>(trials - successes) + beta) -
                            std::lgamma(static_cast<double>(trials) + alpha + beta);
    const double prior_beta = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
    return log_choose + log_beta - prior_beta;
}

std::pair<double, double> artifact_fraction_pseudocounts(const int state) {
    if (state < 8) return {1.0, 9.0};
    if (state == 8) return {3.0, 10000.0};
    if (state == 9) return {5.0, 5.0};
    if (state == 10) return {2.0, 5.0};
    return {10000.0, 3.0};
}

std::pair<double, double> orientation_fraction_pseudocounts(const int state) {
    if (state < 4) return {100.0, 1.0};
    if (state < 8) return {1.0, 100.0};
    return {10.0, 10.0};
}

std::array<double, 12> flat_artifact_prior(const char ref) {
    std::array<double, 12> result{};
    const int ref_index = base_index(std::string(1, ref));
    if (ref_index < 0) throw std::runtime_error("BAD_INPUT: invalid reference base in orientation context");
    constexpr double denominator = 1.0 / 10.0;
    result.fill(denominator);
    result[ref_index] = 0.0;
    result[4 + ref_index] = 0.0;
    return result;
}

[[maybe_unused]] std::array<double, 12> compute_orientation_responsibilities(
    const char ref, const char alt, const int alt_depth, const int alt_f1r2,
    const int depth, const std::array<double, 12>& prior) {
    std::array<double, 12> log_values{};
    log_values.fill(-std::numeric_limits<double>::infinity());
    const int ref_index = base_index(std::string(1, ref));
    const int alt_index = base_index(std::string(1, alt));
    if (ref_index < 0 || alt_index < 0 || depth < 0 || alt_depth < 0 || alt_f1r2 < 0)
        throw std::runtime_error("BAD_INPUT: invalid orientation EM observation");
    for (int state = 0; state < 12; ++state) {
        if (state == ref_index || state == 4 + ref_index) continue; // ref-to-ref artifacts
        if (state < 8 && state % 4 != alt_index) continue; // incompatible artifact allele
        if (!(prior[state] > 0.0) || !std::isfinite(prior[state])) continue;
        const auto af = artifact_fraction_pseudocounts(state);
        const auto orientation = orientation_fraction_pseudocounts(state);
        log_values[state] = std::log(prior[state]) +
            log_beta_binomial_probability(alt_depth, depth, af.first, af.second) +
            log_beta_binomial_probability(alt_f1r2, alt_depth, orientation.first, orientation.second);
    }
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto value : log_values) maximum = std::max(maximum, value);
    if (!std::isfinite(maximum)) throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: orientation posterior has no finite state");
    double normalizer = 0.0;
    for (auto& value : log_values) {
        value = std::isfinite(value) ? std::exp(value - maximum) : 0.0;
        normalizer += value;
    }
    if (!(normalizer > 0.0) || !std::isfinite(normalizer))
        throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: orientation posterior normalization failed");
    for (auto& value : log_values) value /= normalizer;
    return log_values;
}

struct OrientationPriorRow {
    std::string context;
    std::array<double, 12> probabilities{};
    std::uint64_t examples = 0;
    std::uint64_t alt_examples = 0;
};

struct FullOrientationResult {
    std::vector<OrientationPriorRow> rows;
    std::size_t records = 0;
    std::uint64_t f1r2_alt_support = 0;
    std::uint64_t r1f2_alt_support = 0;
    std::size_t contexts_with_data = 0;
    std::size_t em_iterations = 0;
    std::size_t em_kokkos_observations = 0;
    std::size_t em_kokkos_batches = 0;
};

struct OrientationEmObservation {
    int alt_index = 0;
    int alt_depth = 0;
    int alt_f1r2 = 0;
    int depth = 0;
    std::uint64_t weight = 0;
};

KOKKOS_INLINE_FUNCTION double orientation_em_log_beta_binomial(
    const int successes, const int trials, const double alpha, const double beta) {
    if (successes < 0 || trials < 0 || successes > trials || !(alpha > 0.0) || !(beta > 0.0))
        return -1.0e300;
    const auto log_choose = Kokkos::lgamma(static_cast<double>(trials) + 1.0) -
                            Kokkos::lgamma(static_cast<double>(successes) + 1.0) -
                            Kokkos::lgamma(static_cast<double>(trials - successes) + 1.0);
    return log_choose + Kokkos::lgamma(static_cast<double>(successes) + alpha) +
           Kokkos::lgamma(static_cast<double>(trials - successes) + beta) -
           Kokkos::lgamma(static_cast<double>(trials) + alpha + beta) +
           Kokkos::lgamma(alpha + beta) - Kokkos::lgamma(alpha) - Kokkos::lgamma(beta);
}

// The standard CollectF1R2Counts path can contain thousands of histogram bins
// per context.  Evaluate all 12 latent states for an observation in one
// Kokkos work item, then perform the final weighted sum on Host in input order.
// Keeping the reduction order on Host is intentional: it preserves the Java
// EM's deterministic update order while moving the expensive beta-binomial
// state evaluation to the selected Kokkos execution space.
std::vector<std::array<double, 12>> orientation_em_responsibilities_kokkos(
    int ref_index, const std::array<double, 12>& prior,
    const std::vector<OrientationEmObservation>& observations,
    KernelStats* kernel_stats = nullptr) {
    if (observations.empty()) return {};
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<int*> alt("orientation_em_alt", observations.size());
    Kokkos::View<int*> alt_depth("orientation_em_alt_depth", observations.size());
    Kokkos::View<int*> alt_f1r2("orientation_em_alt_f1r2", observations.size());
    Kokkos::View<int*> depth("orientation_em_depth", observations.size());
    Kokkos::View<double*> prior_device("orientation_em_prior", 12);
    Kokkos::View<double**> log_values("orientation_em_log_values", observations.size(), 12);
    fastgatk::core::HostBatch host_batch("learn-orientation-em-v1");
    host_batch.records = observations.size();
    host_batch.bytes = observations.size() * sizeof(OrientationEmObservation) +
                       prior.size() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("learn-orientation-em");
    kernel_plan.begin_prepare(host_batch);
    auto alt_host = Kokkos::create_mirror_view(alt);
    auto alt_depth_host = Kokkos::create_mirror_view(alt_depth);
    auto alt_f1r2_host = Kokkos::create_mirror_view(alt_f1r2);
    auto depth_host = Kokkos::create_mirror_view(depth);
    auto prior_host = Kokkos::create_mirror_view(prior_device);
    for (std::size_t i = 0; i < observations.size(); ++i) {
        alt_host(i) = observations[i].alt_index;
        alt_depth_host(i) = observations[i].alt_depth;
        alt_f1r2_host(i) = observations[i].alt_f1r2;
        depth_host(i) = observations[i].depth;
    }
    for (int state = 0; state < 12; ++state) prior_host(state) = prior[static_cast<std::size_t>(state)];
    Kokkos::deep_copy(alt, alt_host);
    Kokkos::deep_copy(alt_depth, alt_depth_host);
    Kokkos::deep_copy(alt_f1r2, alt_f1r2_host);
    Kokkos::deep_copy(depth, depth_host);
    Kokkos::deep_copy(prior_device, prior_host);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(observations.size());
    device_batch.bind("alt", alt);
    device_batch.bind("alt_depth", alt_depth);
    device_batch.bind("alt_f1r2", alt_f1r2);
    device_batch.bind("depth", depth);
    device_batch.bind("prior", prior_device);
    device_batch.bind("log_values", log_values);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "fastgatk_orientation_em_responsibilities",
        Kokkos::RangePolicy<ExecSpace>(0, observations.size()),
        KOKKOS_LAMBDA(const std::size_t observation) {
            const int observed_alt = alt(observation);
            const int observed_depth = depth(observation);
            const int observed_alt_depth = alt_depth(observation);
            const int observed_alt_f1r2 = alt_f1r2(observation);
            double maximum = -1.0e300;
            for (int state = 0; state < 12; ++state) {
                double value = -1.0e300;
                if (state != ref_index && state != 4 + ref_index &&
                    !(state < 8 && state % 4 != observed_alt)) {
                    const double alpha_af = state < 8 ? 1.0 : state == 8 ? 3.0 : state == 9 ? 5.0 : state == 10 ? 2.0 : 10000.0;
                    const double beta_af = state < 8 ? 9.0 : state == 8 ? 10000.0 : state == 9 ? 5.0 : state == 10 ? 5.0 : 3.0;
                    const double alpha_orientation = state < 4 ? 100.0 : state < 8 ? 1.0 : 10.0;
                    const double beta_orientation = state < 4 ? 1.0 : state < 8 ? 100.0 : 10.0;
                    // The merged design matrix is already strand-canonicalized
                    // (reverse-complement rows swap F1R2/F2R1).  Therefore the
                    // EM observation is the canonical F1R2 count for every
                    // latent state, matching GATK's DesignMatrix convention.
                    const int orientation_successes = observed_alt_f1r2;
                    const double prior_value = Kokkos::fmax(prior_device(state), 1.0e-300);
                    value = Kokkos::log(prior_value) +
                        orientation_em_log_beta_binomial(observed_alt_depth, observed_depth,
                                                         alpha_af, beta_af) +
                        orientation_em_log_beta_binomial(orientation_successes, observed_alt_depth,
                                                         alpha_orientation, beta_orientation);
                }
                log_values(observation, state) = value;
                maximum = Kokkos::fmax(maximum, value);
            }
            double normalizer = 0.0;
            for (int state = 0; state < 12; ++state) {
                const double value = log_values(observation, state) < -1.0e299
                    ? 0.0 : Kokkos::exp(log_values(observation, state) - maximum);
                log_values(observation, state) = value;
                normalizer += value;
            }
            for (int state = 0; state < 12; ++state)
                log_values(observation, state) = normalizer > 0.0
                    ? log_values(observation, state) / normalizer : 0.0;
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += observations.size();
        kernel_stats->prepare_seconds += kernel_plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += kernel_plan.telemetry().execute_seconds;
    }
    auto log_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, log_values);
    std::vector<std::array<double, 12>> result(observations.size());
    for (std::size_t observation = 0; observation < observations.size(); ++observation)
        for (int state = 0; state < 12; ++state)
            result[observation][static_cast<std::size_t>(state)] = log_host(observation, state);
    return result;
}

std::vector<std::uint64_t> histogram_for(
    const std::unordered_map<std::string, std::vector<std::uint64_t>>& histograms,
    const std::string& label, const std::size_t max_depth) {
    const auto found = histograms.find(label);
    if (found == histograms.end()) return std::vector<std::uint64_t>(max_depth + 1, 0);
    // Preserve all source bins.  The Java engine applies maxDepth when
    // materializing EM observations, not while parsing or summing histograms.
    return found->second;
}

std::vector<std::uint64_t> add_histograms(const std::vector<std::uint64_t>& left,
                                          const std::vector<std::uint64_t>& right) {
    const std::size_t size = std::max(left.size(), right.size());
    std::vector<std::uint64_t> result(size, 0);
    for (std::size_t i = 0; i < size; ++i) {
        if (i < left.size()) result[i] += left[i];
        if (i < right.size()) result[i] += right[i];
    }
    return result;
}

AltSiteRecord reverse_complement_record(const AltSiteRecord& input) {
    AltSiteRecord result = input;
    result.context = reverse_complement(input.context);
    result.alt = complement_base(input.alt);
    result.ref_f1r2 = input.ref_count - input.ref_f1r2;
    result.alt_f1r2 = input.alt_count - input.alt_f1r2;
    return result;
}

std::string make_full_prior_table(const std::string& sample,
                                  const std::vector<OrientationPriorRow>& rows) {
    std::ostringstream table;
    table << "#<METADATA>SAMPLE=" << sample << '\n';
    table << "context\trev_comp\tf1r2_a\tf1r2_c\tf1r2_g\tf1r2_t\tf2r1_a\tf2r1_c\tf2r1_g\tf2r1_t\t"
             "hom_ref\tgermline_het\tsomatic_het\thom_var\tnum_examples\tnum_alt_examples\n";
    std::map<std::string, const OrientationPriorRow*> by_context;
    for (const auto& row : rows) by_context[row.context] = &row;
    for (const auto& context : gatk_context_order()) {
        const auto found = by_context.find(context);
        if (found == by_context.end()) continue;
        const auto& row = *found->second;
        table << row.context << '\t' << reverse_complement(row.context);
        for (const auto value : row.probabilities) table << '\t' << java_double(value);
        table << '\t' << row.examples << '\t' << row.alt_examples << '\n';
    }
    return table.str();
}

FullOrientationResult learn_full_orientation_model(const StandardOrientationInput& input,
                                                   const Options& options,
                                                   KernelStats* kernel_stats = nullptr) {
    const auto contexts = all_contexts();
    const std::size_t max_depth = static_cast<std::size_t>(options.max_depth);
    FullOrientationResult result;
    result.rows.reserve(contexts.size());
    for (const auto& context : contexts) {
        const auto reverse = reverse_complement(context);
        if (context > reverse) continue; // ArtifactPriorCollection writes RC from the canonical row.
        const char ref = context[1];

        const auto ref_hist = add_histograms(
            histogram_for(input.ref_hist, context, max_depth),
            histogram_for(input.ref_hist, reverse, max_depth));

        std::vector<AltSiteRecord> design;
        if (context == reverse) {
            // Java's mergeDesignMatrices intentionally appends the reverse
            // complement of the same list for palindromic contexts.
            for (const auto& record : input.alt_sites) {
                if (record.context == context) design.push_back(record);
            }
            const auto original = design;
            for (const auto& record : original) design.push_back(reverse_complement_record(record));
        } else {
            for (const auto& record : input.alt_sites) {
                if (record.context == context) design.push_back(record);
                if (record.context == reverse) design.push_back(reverse_complement_record(record));
            }
        }

        struct AltHistogram {
            char alt;
            bool f1r2;
            std::vector<std::uint64_t> counts;
        };
        std::vector<AltHistogram> alt_histograms;
        for (const char alt : std::array<char, 4>{'A', 'C', 'G', 'T'}) {
            if (alt == ref) continue;
            for (const bool f1r2 : {true, false}) {
                const char orientation = f1r2 ? 'F' : 'R';
                const char reverse_orientation = f1r2 ? 'R' : 'F';
                const std::string label = context + "_" + alt + "_" +
                    (f1r2 ? "F1R2" : "F2R1");
                const std::string reverse_label = reverse + "_" + complement_base(alt) + "_" +
                    (reverse_orientation == 'F' ? "F1R2" : "F2R1");
                alt_histograms.push_back(AltHistogram{
                    alt, f1r2,
                    add_histograms(histogram_for(input.alt_hist, label, max_depth),
                                   histogram_for(input.alt_hist, reverse_label, max_depth))});
                (void)orientation;
            }
        }

        std::uint64_t ref_examples = 0;
        for (std::size_t depth = 1; depth < ref_hist.size(); ++depth) ref_examples += ref_hist[depth];
        std::uint64_t alt_hist_examples = 0;
        for (const auto& histogram : alt_histograms)
            for (std::size_t depth = 1; depth < histogram.counts.size(); ++depth)
                alt_hist_examples += histogram.counts[depth];
        const std::uint64_t alt_examples = static_cast<std::uint64_t>(design.size()) + alt_hist_examples;

        OrientationPriorRow row;
        row.context = context;
        row.probabilities = flat_artifact_prior(ref);
        if (ref_examples == 0 || design.empty()) {
            result.rows.push_back(std::move(row));
            auto reverse_row = result.rows.back();
            reverse_row.context = reverse;
            const auto canonical = result.rows.back().probabilities;
            for (int state = 0; state < 8; ++state) reverse_row.probabilities[state] = canonical[7 - state];
            for (int state = 8; state < 12; ++state) reverse_row.probabilities[state] = canonical[state];
            result.rows.push_back(std::move(reverse_row));
            continue;
        }
        ++result.contexts_with_data;
        std::array<double, 12> prior = row.probabilities;
        std::vector<OrientationEmObservation> observations;
        observations.reserve(ref_hist.size() + design.size() +
                             alt_histograms.size() * max_depth);
        const int ref_index = base_index(std::string(1, ref));
        const auto em_depth_limit = std::min(ref_hist.size(), max_depth + 1);
        for (std::size_t depth = 1; depth < em_depth_limit; ++depth) {
            if (ref_hist[depth] == 0) continue;
            observations.push_back(OrientationEmObservation{
                ref_index, 0, 0, static_cast<int>(depth), ref_hist[depth]});
        }
        for (const auto& record : design) {
            observations.push_back(OrientationEmObservation{
                base_index(std::string(1, record.alt)), record.alt_count, record.alt_f1r2,
                record.ref_count + record.alt_count, 1});
        }
        for (const auto& histogram : alt_histograms) {
            const int alt_index = base_index(std::string(1, histogram.alt));
            const auto em_depth_limit = std::min(histogram.counts.size(), max_depth + 1);
            for (std::size_t depth = 1; depth < em_depth_limit; ++depth) {
                if (histogram.counts[depth] == 0) continue;
                observations.push_back(OrientationEmObservation{
                    alt_index, 1, histogram.f1r2 ? 1 : 0,
                    static_cast<int>(depth), histogram.counts[depth]});
            }
        }
        result.em_kokkos_observations += observations.size();
        std::size_t iterations = 0;
        for (; iterations < static_cast<std::size_t>(options.max_iterations); ++iterations) {
            std::array<double, 12> effective{};
            const auto responsibilities = orientation_em_responsibilities_kokkos(
                ref_index, prior, observations, kernel_stats);
            ++result.em_kokkos_batches;
            for (std::size_t observation = 0; observation < observations.size(); ++observation)
                for (int state = 0; state < 12; ++state)
                    effective[state] += static_cast<double>(observations[observation].weight) *
                        responsibilities[observation][static_cast<std::size_t>(state)];
            const auto old_prior = prior;
            double normalizer = 0.0;
            for (int state = 0; state < 12; ++state) {
                prior[state] = effective[state] + row.probabilities[state];
                normalizer += prior[state];
            }
            for (auto& value : prior) value /= normalizer;
            double l2 = 0.0;
            for (int state = 0; state < 12; ++state) {
                const double delta = old_prior[state] - prior[state];
                l2 += delta * delta;
            }
            if (std::sqrt(l2) <= options.convergence_threshold) {
                ++iterations;
                break;
            }
        }
        result.em_iterations = std::max(result.em_iterations, iterations);
        row.probabilities = prior;
        row.examples = ref_examples + alt_examples;
        row.alt_examples = alt_examples;
        result.rows.push_back(std::move(row));
        auto reverse_row = result.rows.back();
        reverse_row.context = reverse;
        const auto canonical = result.rows.back().probabilities;
        for (int state = 0; state < 8; ++state) reverse_row.probabilities[state] = canonical[7 - state];
        for (int state = 8; state < 12; ++state) reverse_row.probabilities[state] = canonical[state];
        result.rows.push_back(std::move(reverse_row));
    }
    for (const auto& record : input.alt_sites) {
        result.f1r2_alt_support += static_cast<std::uint64_t>(record.alt_f1r2);
        result.r1f2_alt_support += static_cast<std::uint64_t>(record.alt_count - record.alt_f1r2);
    }
    for (const auto& [label, counts] : input.alt_hist) {
        const auto separator = label.rfind('_');
        if (separator == std::string::npos) continue;
        const bool f1r2 = label.substr(separator + 1) == "F1R2";
        std::uint64_t total = 0;
        for (std::size_t depth = 1; depth < counts.size(); ++depth) total += counts[depth];
        if (f1r2) result.f1r2_alt_support += total;
        else result.r1f2_alt_support += total;
    }
    result.records = input.alt_records + input.alt_sites.size();
    return result;
}

std::string url_encode(const std::string& value) {
    std::ostringstream encoded;
    encoded << std::uppercase << std::hex;
    for (const unsigned char character : value) {
        if (std::isalnum(character) || character == '.' || character == '-' || character == '_')
            encoded << static_cast<char>(character);
        else
            encoded << '%' << std::setw(2) << std::setfill('0') << static_cast<int>(character)
                    << std::setfill(' ');
    }
    return encoded.str();
}

std::string make_prior_table(const std::string& sample,
                             const std::array<std::uint64_t, 8>& counts,
                             std::uint64_t records) {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    std::uint64_t total_alt = 0;
    for (const auto count : counts) total_alt += count;
    const double scale = static_cast<double>(std::max<std::uint64_t>(total_alt, 1));
    std::ostringstream table;
    table << "#<METADATA>SAMPLE=" << sample << '\n';
    table << "context\trev_comp\tf1r2_a\tf1r2_c\tf1r2_g\tf1r2_t\tf2r1_a\tf2r1_c\tf2r1_g\tf2r1_t\t"
             "hom_ref\tgermline_het\tsomatic_het\thom_var\tnum_examples\tnum_alt_examples\n";
    for (const auto& context : gatk_context_order()) {
        const char middle = context[1];
        std::array<double, 12> mass{};
        for (int orientation = 0; orientation < 2; ++orientation) {
            for (int alt = 0; alt < 4; ++alt) {
                const auto state = orientation * 4 + alt;
                if (bases[alt] == middle) mass[state] = 0.0;
                else mass[state] = 1.0 + 8.0 * static_cast<double>(counts[state]) / scale;
            }
        }
        mass[8] = mass[9] = mass[10] = mass[11] = 1.0;
        double normalizer = 0.0;
        for (const auto value : mass) normalizer += value;
        table << context << '\t' << reverse_complement(context);
        for (const auto value : mass) table << '\t' << java_double(value / normalizer);
        table << '\t' << records << '\t' << records << '\n';
    }
    return table.str();
}

#if FASTGATK_HAS_ZLIB
void append_tar_field(std::array<unsigned char, 512>& header, std::size_t offset,
                      std::size_t length, const std::string& value) {
    const auto count = std::min(length, value.size());
    std::copy_n(value.data(), count, header.data() + offset);
}

void append_tar_octal(std::array<unsigned char, 512>& header, std::size_t offset,
                      std::size_t length, std::uint64_t value) {
    std::ostringstream octal;
    octal << std::oct << value;
    const auto text = octal.str();
    const auto begin = text.size() < length - 1 ? 0 : text.size() - (length - 1);
    std::fill(header.begin() + offset, header.begin() + offset + length, '0');
    std::copy(text.begin() + begin, text.end(), header.begin() + offset + (length - 1) - (text.size() - begin));
    header[offset + length - 1] = '\0';
}

void write_tar_gz(const std::string& output,
                  const std::vector<std::pair<std::string, std::string>>& members) {
    if (members.empty())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: orientation archive has no members");
    gzFile stream = gzopen(output.c_str(), "wb");
    if (!stream) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create archive: " + output);
    auto write = [&](const void* data, std::size_t size) {
        const auto* bytes = static_cast<const unsigned char*>(data);
        while (size != 0) {
            const auto chunk = static_cast<unsigned int>(std::min<std::size_t>(size, std::numeric_limits<unsigned int>::max()));
            if (gzwrite(stream, bytes, chunk) != static_cast<int>(chunk)) {
                gzclose(stream);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write archive: " + output);
            }
            bytes += chunk;
            size -= chunk;
        }
    };
    const std::array<unsigned char, 512> padding{};
    for (const auto& [member, content] : members) {
    std::array<unsigned char, 512> header{};
    append_tar_field(header, 0, 100, member);
    append_tar_field(header, 100, 8, "0000644\0");
    append_tar_field(header, 108, 8, "0000000\0");
    append_tar_field(header, 116, 8, "0000000\0");
    append_tar_octal(header, 124, 12, content.size());
    append_tar_field(header, 136, 12, "00000000000\0");
    std::fill(header.begin() + 148, header.begin() + 156, ' ');
    header[156] = '0';
    append_tar_field(header, 257, 6, "ustar\0");
    append_tar_field(header, 263, 2, "00");
    std::uint32_t checksum = 0;
    for (const auto byte : header) checksum += byte;
    append_tar_octal(header, 148, 8, checksum);

    write(header.data(), header.size());
    write(content.data(), content.size());
    const auto remainder = content.size() % 512;
    if (remainder != 0) write(padding.data(), 512 - remainder);
    }
    write(padding.data(), padding.size());
    write(padding.data(), padding.size());
    if (gzclose(stream) != Z_OK) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize archive: " + output);
}

void write_tar_gz(const std::string& output, const std::string& member, const std::string& content) {
    write_tar_gz(output, std::vector<std::pair<std::string, std::string>>{{member, content}});
}
#else
void write_tar_gz(const std::string&, const std::vector<std::pair<std::string, std::string>>&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: zlib is required for .tar.gz orientation priors");
}
void write_tar_gz(const std::string&, const std::string&, const std::string&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: zlib is required for .tar.gz orientation priors");
}
#endif

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        KernelStats kernel_stats;

        std::vector<OrientationRecord> records;
        bool standard_input = false;
        StandardOrientationInput standard;
        std::map<std::string, StandardOrientationInput> standards;
        std::string sample = options.sample;
        for (const auto& input : options.inputs) {
            const bool tar_input = input.size() >= 7 && input.substr(input.size() - 7) == ".tar.gz";
            if (tar_input) {
                auto parsed = read_standard_tar_input(input, sample, static_cast<std::size_t>(options.max_depth));
                if (parsed.standard) {
                    if (!records.empty())
                        throw std::runtime_error("BAD_INPUT: cannot mix legacy F1R2 TSV and CollectF1R2Counts tar inputs");
                    standard_input = true;
                    for (const auto& [parsed_sample, parsed_input] : parsed.standard_inputs) {
                        auto& destination = standards[parsed_sample];
                        merge_standard_orientation_input(destination, parsed_input);
                    }
                } else {
                    if (standard_input)
                        throw std::runtime_error("BAD_INPUT: cannot mix CollectF1R2Counts tar and legacy F1R2 TSV inputs");
                    records.insert(records.end(), parsed.legacy_records.begin(), parsed.legacy_records.end());
                }
            } else {
                if (standard_input)
                    throw std::runtime_error("BAD_INPUT: cannot mix CollectF1R2Counts tar and legacy F1R2 TSV inputs");
                auto shard = read_input(input, sample);
                records.insert(records.end(), shard.begin(), shard.end());
            }
        }
        if (standard_input) {
            if (standards.empty())
                throw std::runtime_error("BAD_INPUT: CollectF1R2Counts tar has no reference histogram");
            if (sample != "UNKNOWN" && standards.size() == 1 && standards.begin()->first == "UNKNOWN") {
                auto input = std::move(standards.begin()->second);
                standards.clear();
                input.sample = sample;
                standards.emplace(sample, std::move(input));
            }
            std::vector<std::pair<std::string, std::string>> members;
            members.reserve(standards.size());
            std::uint64_t total_records = 0;
            std::uint64_t total_contexts_with_data = 0;
            std::uint64_t total_f1r2 = 0;
            std::uint64_t total_r1f2 = 0;
            for (auto& [output_sample, input] : standards) {
                if (input.sample.empty()) input.sample = output_sample == "UNKNOWN" ? sample : output_sample;
                if (input.ref_hist.empty())
                    throw std::runtime_error("BAD_INPUT: CollectF1R2Counts tar has no reference histogram for sample " + output_sample);
                const auto effective_sample = input.sample.empty() ? output_sample : input.sample;
                const auto full = learn_full_orientation_model(input, options, &kernel_stats);
                total_records += full.records;
                total_contexts_with_data += full.contexts_with_data;
                total_f1r2 += full.f1r2_alt_support;
                total_r1f2 += full.r1f2_alt_support;
                members.emplace_back(url_encode(effective_sample) + ".orientation_priors",
                                     make_full_prior_table(effective_sample, full.rows));
            }
            write_tar_gz(options.output, members);
            if (!file_complete(options.output))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: orientation archive is missing or empty");

            const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
            std::ofstream manifest(manifest_path);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + manifest_path);
            manifest << "{\"schema_version\":1,\"tool\":\"LearnReadOrientationModel\","
                     << "\"implementation\":\"fastgatk-orientation-beta-binomial-em-v1\",\"status\":\"prototype\","
                     << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\","
                     << "\"primary_output\":\"" << json_escape(options.output) << "\","
                     << "\"primary_output_kind\":\"artifact-prior-tar-gz\",\"sample\":\""
                     << json_escape(standards.size() == 1 ? standards.begin()->first : "MULTIPLE") << "\",\"samples\":"
                     << standards.size() << ",\"compatibility\":{" 
                     << "\"tar_gz\":true,\"artifact_prior_table\":true,\"full_gatk_em\":true,"
                     << "\"collect_f1r2_tar_input\":true,\"em_kokkos\":true,"
                     << "\"quiet\":" << (options.quiet ? "true" : "false") << ",\"tmp_dir\":\""
                     << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                     << json_escape(options.verbosity) << "\",\"bit_identical_to_gatk\":false},"
                     << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"artifact-prior-tar-gz\",\"complete\":true}],"
                     << "\"telemetry\":{" << "\"resources\":" << resources.to_json()
                     << ",\"input_files\":" << options.inputs.size()
                     << ",\"records\":" << total_records
                     << ",\"contexts\":64"
                     << ",\"contexts_with_data\":" << total_contexts_with_data
                     << ",\"em_kokkos_observations\":" << kernel_stats.observations
                     << ",\"em_kokkos_batches\":" << kernel_stats.batches
                     << ",\"em_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                     << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"kernel_execution_policy\":\"RangePolicy\""
                     << ",\"kernel_batches\":" << kernel_stats.batches
                     << ",\"kernel_observations\":" << kernel_stats.observations
                     << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds
                     << ",\"f1r2_alt_support\":" << total_f1r2
                     << ",\"r1f2_alt_support\":" << total_r1f2
                     << ",\"convergence_threshold\":" << options.convergence_threshold
                     << ",\"max_iterations\":" << options.max_iterations
                     << ",\"max_depth\":" << options.max_depth
                     << ",\"quiet\":" << (options.quiet ? "true" : "false")
                     << ",\"tmp_dir\":\"" << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                     << json_escape(options.verbosity) << "\"}}\n";
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize manifest");
            const auto summary_sample = standards.size() == 1 ? standards.begin()->first : "MULTIPLE";
            std::cout << "{\"tool\":\"LearnReadOrientationModel\",\"status\":\"prototype\","
                      << "\"sample\":\"" << json_escape(summary_sample) << "\",\"samples\":" << standards.size()
                      << ",\"records\":" << total_records << ",\"contexts\":64"
                      << ",\"contexts_with_data\":" << total_contexts_with_data
                      << ",\"f1r2_alt_support\":" << total_f1r2
                      << ",\"r1f2_alt_support\":" << total_r1f2 << "}\n";
            Kokkos::finalize();
            return 0;
        }
        if (records.empty()) throw std::runtime_error("BAD_INPUT: no F1R2 records found");

        std::array<std::uint64_t, 8> aggregate{};
        {
            Kokkos::View<int*> alt("orientation_alt", records.size());
            Kokkos::View<std::uint64_t*> f1r2("orientation_f1r2", records.size());
            Kokkos::View<std::uint64_t*> r1f2("orientation_r1f2", records.size());
            auto alt_host = Kokkos::create_mirror_view(alt);
            auto f1_host = Kokkos::create_mirror_view(f1r2);
            auto r1_host = Kokkos::create_mirror_view(r1f2);
            for (std::size_t i = 0; i < records.size(); ++i) {
                alt_host(i) = records[i].alt;
                f1_host(i) = records[i].f1r2;
                r1_host(i) = records[i].r1f2;
            }
            Kokkos::deep_copy(alt, alt_host);
            Kokkos::deep_copy(f1r2, f1_host);
            Kokkos::deep_copy(r1f2, r1_host);
            Kokkos::View<std::uint64_t*> counts("orientation_counts", 8);
            fastgatk::core::HostBatch host_batch("learn-orientation-aggregate-v1");
            host_batch.records = records.size();
            host_batch.bytes = records.size() * sizeof(OrientationRecord);
            fastgatk::core::KernelPlan<Kokkos::DefaultExecutionSpace> kernel_plan(
                "learn-orientation-aggregate");
            kernel_plan.begin_prepare(host_batch);
            Kokkos::deep_copy(counts, std::uint64_t{0});
            fastgatk::core::DeviceBatch<Kokkos::DefaultExecutionSpace> device_batch(records.size());
            device_batch.bind("alt", alt);
            device_batch.bind("f1r2", f1r2);
            device_batch.bind("r1f2", r1f2);
            device_batch.bind("counts", counts);
            Kokkos::DefaultExecutionSpace().fence();
            kernel_plan.end_prepare(device_batch);
            kernel_plan.begin_execute();
            Kokkos::parallel_for("learn_orientation_model_aggregate",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, records.size()),
                KOKKOS_LAMBDA(const std::size_t i) {
                const auto base = alt(i);
                Kokkos::atomic_add(&counts[base], f1r2(i));
                Kokkos::atomic_add(&counts[4 + base], r1f2(i));
            });
            Kokkos::DefaultExecutionSpace().fence();
            kernel_plan.end_execute();
            ++kernel_stats.batches;
            kernel_stats.observations += records.size();
            kernel_stats.prepare_seconds += kernel_plan.telemetry().prepare_seconds;
            kernel_stats.execute_seconds += kernel_plan.telemetry().execute_seconds;
            auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, counts);
            for (int i = 0; i < 8; ++i) aggregate[i] = host_counts(i);
        }

        const auto table = make_prior_table(sample, aggregate, records.size());
        const auto member = url_encode(sample) + ".orientation_priors";
        write_tar_gz(options.output, member, table);
        if (!file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: orientation archive is missing or empty");

        const auto total_f1r2 = aggregate[0] + aggregate[1] + aggregate[2] + aggregate[3];
        const auto total_r1f2 = aggregate[4] + aggregate[5] + aggregate[6] + aggregate[7];
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"LearnReadOrientationModel\","
                 << "\"implementation\":\"fastgatk-orientation-strand-mixture-v1\",\"status\":\"prototype\","
                 << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\","
                 << "\"primary_output_kind\":\"artifact-prior-tar-gz\",\"sample\":\""
                 << json_escape(sample) << "\",\"compatibility\":{"
                 << "\"tar_gz\":true,"
                 << "\"artifact_prior_table\":true,\"full_gatk_em\":false,"
                 << "\"collect_f1r2_tar_input\":true,"
                 << "\"quiet\":" << (options.quiet ? "true" : "false") << ",\"tmp_dir\":\""
                 << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                 << json_escape(options.verbosity) << "\",\"bit_identical_to_gatk\":false},"
                 << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                 << "\",\"kind\":\"artifact-prior-tar-gz\",\"complete\":true}],"
                 << "\"telemetry\":{" << "\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size() << ",\"records\":" << records.size()
                 << ",\"f1r2_alt_support\":" << total_f1r2 << ",\"r1f2_alt_support\":" << total_r1f2
                 << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                 << ",\"kernel_execution_policy\":\"RangePolicy\""
                 << ",\"kernel_batches\":" << kernel_stats.batches
                 << ",\"kernel_observations\":" << kernel_stats.observations
                 << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                 << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds
                 << ",\"convergence_threshold\":" << options.convergence_threshold
                 << ",\"max_iterations\":" << options.max_iterations
                 << ",\"max_depth\":" << options.max_depth
                 << ",\"quiet\":" << (options.quiet ? "true" : "false")
                 << ",\"tmp_dir\":\"" << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                 << json_escape(options.verbosity) << "\"}}\n";
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize manifest");
        std::cout << "{\"tool\":\"LearnReadOrientationModel\",\"status\":\"prototype\","
                  << "\"sample\":\"" << json_escape(sample) << "\",\"records\":" << records.size()
                  << ",\"contexts\":64,\"f1r2_alt_support\":" << total_f1r2
                  << ",\"r1f2_alt_support\":" << total_r1f2 << "}\n";
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
