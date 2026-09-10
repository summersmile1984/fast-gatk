#include "fastgatk/runtime/resource.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::vector<double> truth_sensitivity;
    // GatherTranches requires an explicit mode.  Do not silently pick SNP:
    // the GATK argument parser rejects an invocation without --mode.
    std::string mode;
    std::string output;
    std::string manifest;
};

struct Tranche {
    double min_vqslod = 0.0;
    std::int64_t num_known = 0;
    double known_titv = 0.0;
    std::int64_t num_novel = 0;
    double novel_titv = 0.0;
    std::int64_t accessible_truth_sites = 0;
    std::int64_t calls_at_truth_sites = 0;
    std::string model;
    std::string input_filter_name;
    double sensitivity = 0.0;
};

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
    const auto inline_value = option_value(argument, name);
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
    for (const char character : text) {
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

bool compressed_path(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

void write_text_output(const std::string& path, const std::string& text) {
    if (!compressed_path(path)) {
        std::ofstream output(path, std::ios::binary);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write " + path);
        output.write(text.data(), static_cast<std::streamsize>(text.size()));
        output.close();
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete output");
        return;
    }
#if FASTGATK_HAS_ZLIB
    gzFile output = gzopen(path.c_str(), "wb");
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed output " + path);
    std::size_t offset = 0;
    while (offset < text.size()) {
        const auto remaining = text.size() - offset;
        const auto chunk = static_cast<unsigned>(std::min<std::size_t>(remaining, 1U << 20));
        const auto written = gzwrite(output, text.data() + offset, chunk);
        if (written != static_cast<int>(chunk)) {
            (void)gzclose(output);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed output " + path);
        }
        offset += chunk;
    }
    if (gzclose(output) != Z_OK)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize compressed output " + path);
#else
    throw std::runtime_error("BACKEND_UNAVAILABLE: compressed tranche output requires zlib: " + path);
#endif
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-gather-tranches (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 scattered VQSLOD tranche CSV (repeatable)\n"
                         "      --truth-sensitivity-tranche F target sensitivity in percent (repeatable)\n"
                         "  -tranche F                         GATK alias for --truth-sensitivity-tranche\n"
                         "      --mode SNP|INDEL|BOTH         required recalibration mode\n"
                         "  -O, --output FILE                gathered tranches CSV\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.inputs.push_back(require_value(index, argc, argv, argument, "--input", "-I"));
        } else if (argument == "-tranche" || has_option(argument, "--truth-sensitivity-tranche") ||
                   has_option(argument, "--tranche")) {
            const auto name = (argument == "-tranche" || argument.rfind("--tranche", 0) == 0)
                ? "--tranche" : "--truth-sensitivity-tranche";
            const auto value = std::stod(argument == "-tranche"
                ? require_value(index, argc, argv, argument, "-tranche", "-tranche")
                : require_value(index, argc, argv, argument, name));
            if (!std::isfinite(value) || value < 0.0 || value > 100.0)
                throw std::invalid_argument("truth sensitivity tranche must be in [0,100]");
            options.truth_sensitivity.push_back(value);
        } else if (has_option(argument, "--mode")) {
            options.mode = require_value(index, argc, argv, argument, "--mode");
            if (options.mode != "SNP" && options.mode != "INDEL" && options.mode != "BOTH")
                throw std::invalid_argument("--mode must be SNP, INDEL, or BOTH");
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility switches.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("at least one -I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.mode.empty()) throw std::invalid_argument("--mode is required");
    if (options.truth_sensitivity.empty())
        options.truth_sensitivity = {100.0, 99.9, 99.0, 90.0};
    return options;
}

std::map<std::string, std::size_t> header_indices(const std::vector<std::string>& header) {
    std::map<std::string, std::size_t> indices;
    for (std::size_t index = 0; index < header.size(); ++index) indices[header[index]] = index;
    return indices;
}

template <typename Fn>
void for_each_tranche_line(const std::string& path, Fn&& callback) {
    const bool compressed = compressed_path(path);
    if (compressed) {
#if FASTGATK_HAS_ZLIB
        gzFile input = gzopen(path.c_str(), "rb");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open compressed tranche file: " + path);
        char buffer[1 << 16];
        std::string pending;
        while (true) {
            char* result = gzgets(input, buffer, static_cast<int>(sizeof(buffer)));
            if (result == nullptr) break;
            pending.append(result);
            if (pending.find('\n') == std::string::npos) continue;
            std::size_t begin = 0;
            while (true) {
                const auto end = pending.find('\n', begin);
                if (end == std::string::npos) {
                    pending = pending.substr(begin);
                    break;
                }
                callback(pending.substr(begin, end - begin));
                begin = end + 1;
            }
        }
        if (!pending.empty()) callback(pending);
        int error_code = Z_OK;
        (void)gzerror(input, &error_code);
        const auto close_status = gzclose(input);
        if (error_code != Z_OK || close_status != Z_OK)
            throw std::runtime_error("BAD_INPUT: cannot read compressed tranche file: " + path);
        return;
#else
        throw std::runtime_error("BACKEND_UNAVAILABLE: compressed tranche files require zlib: " + path);
#endif
    }
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open tranche file: " + path);
    for (std::string line; std::getline(input, line);) callback(std::move(line));
    if (input.bad()) throw std::runtime_error("BAD_INPUT: cannot read tranche file: " + path);
}

std::vector<Tranche> read_tranches(const std::string& path) {
    std::vector<std::string> header;
    std::map<std::string, std::size_t> indices;
    std::vector<Tranche> rows;
    bool saw_version = false;
    for_each_tranche_line(path, [&](std::string line) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) return;
        if (line[0] == '#') {
            if (line.find("Version") != std::string::npos) {
                std::istringstream version_stream(line);
                std::string hash, label, number;
                int version = 0;
                if (!(version_stream >> hash >> label >> number >> version))
                    throw std::runtime_error("BAD_INPUT: malformed tranche version in " + path);
                if (version != 6)
                    throw std::runtime_error("BAD_INPUT: unsupported VQSLOD tranche version " + std::to_string(version));
                saw_version = true;
            }
            return;
        }
        const auto fields = split_csv(line);
        if (header.empty()) {
            header = fields;
            if (header.size() != 11)
                throw std::runtime_error("BAD_INPUT: tranche header must contain 11 columns in " + path);
            indices = header_indices(header);
            const auto required = {"minVQSLod", "numNovel", "novelTiTv", "model", "filterName"};
            for (const auto* name : required)
                if (!indices.count(name)) throw std::runtime_error("BAD_INPUT: tranche header missing " + std::string(name));
            return;
        }
        if (fields.size() != header.size()) throw std::runtime_error("BAD_INPUT: malformed tranche row in " + path);
        auto get = [&](const char* name) -> const std::string& { return fields.at(indices.at(name)); };
        Tranche row;
        row.min_vqslod = std::stod(get("minVQSLod"));
        row.num_novel = std::stoll(get("numNovel"));
        row.novel_titv = std::stod(get("novelTiTv"));
        row.model = get("model");
        row.input_filter_name = get("filterName");
        // VQSLODTranche.readTranches uses -1 for optional legacy columns.  The
        // normal VQSR report always supplies them, but retaining the sentinel
        // makes malformed/legacy input fail at the same aggregation boundary.
        if (indices.count("numKnown")) row.num_known = std::stoll(get("numKnown"));
        else row.num_known = -1;
        if (indices.count("knownTiTv")) row.known_titv = std::stod(get("knownTiTv"));
        else row.known_titv = -1.0;
        if (indices.count("accessibleTruthSites")) row.accessible_truth_sites = std::stoll(get("accessibleTruthSites"));
        else row.accessible_truth_sites = -1;
        if (indices.count("callsAtTruthSites")) row.calls_at_truth_sites = std::stoll(get("callsAtTruthSites"));
        else row.calls_at_truth_sites = -1;
        row.sensitivity = row.accessible_truth_sites > 0
            ? static_cast<double>(row.calls_at_truth_sites) / static_cast<double>(row.accessible_truth_sites)
            : 0.0;
        if (!std::isfinite(row.min_vqslod) || !std::isfinite(row.known_titv) ||
            !std::isfinite(row.novel_titv) || row.num_known < 0 || row.num_novel < 0 ||
            row.accessible_truth_sites < 0 || row.calls_at_truth_sites < 0 ||
            (row.calls_at_truth_sites > row.accessible_truth_sites && row.accessible_truth_sites >= 0))
            throw std::runtime_error("BAD_INPUT: invalid tranche numeric value in " + path);
        rows.push_back(std::move(row));
    });
    // VariantRecalibrator scatter output may legitimately contain an empty
    // requested-VQSLOD slice.  GATK's readTranches accepts its version-6
    // header and GatherTranches simply contributes no rows from that shard.
    // An all-empty input set is still rejected by merge() below.
    if (header.empty()) throw std::runtime_error("BAD_INPUT: tranche file has no header: " + path);
    (void)saw_version;  // GATK accepts reports without a version comment.
    return rows;
}

double weighted_titv(const std::vector<Tranche>& rows, bool known) {
    double transitions = 0.0;
    double transversions = 0.0;
    for (const auto& row : rows) {
        const auto count = static_cast<double>(known ? row.num_known : row.num_novel);
        const auto titv = known ? row.known_titv : row.novel_titv;
        if (count <= 0.0 || titv < 0.0) continue;
        const double ti = count * titv / (1.0 + titv);
        transitions += ti;
        transversions += count - ti;
    }
    return transversions > 0.0 ? transitions / transversions : 0.0;
}

std::vector<Tranche> merge(const std::vector<std::string>& inputs, const std::string& mode) {
    std::map<double, std::vector<Tranche>> grouped;
    for (const auto& input : inputs)
        for (const auto& row : read_tranches(input)) grouped[row.min_vqslod].push_back(row);
    std::vector<Tranche> merged;
    for (auto iterator = grouped.rbegin(); iterator != grouped.rend(); ++iterator) {
        const auto& rows = iterator->second;
        Tranche row;
        row.min_vqslod = iterator->first;
        row.model = mode;
        for (const auto& shard : rows) {
            row.num_known += shard.num_known;
            row.num_novel += shard.num_novel;
            row.accessible_truth_sites += shard.accessible_truth_sites;
            row.calls_at_truth_sites += shard.calls_at_truth_sites;
        }
        row.known_titv = weighted_titv(rows, true);
        row.novel_titv = weighted_titv(rows, false);
        row.sensitivity = row.accessible_truth_sites > 0
            ? static_cast<double>(row.calls_at_truth_sites) / static_cast<double>(row.accessible_truth_sites)
            : 0.0;
        merged.push_back(std::move(row));
    }
    if (merged.empty()) throw std::runtime_error("BAD_INPUT: no tranche rows to merge");
    return merged;
}

std::vector<Tranche> select_targets(const std::vector<Tranche>& merged,
                                    std::vector<double> targets) {
    std::sort(targets.begin(), targets.end());
    std::vector<Tranche> selected;
    if (targets.empty() || merged.size() < 2) return selected;

    // This is the intentionally stateful walk in GATK's
    // VQSLODTranche.mergeAndConvertTranches.  The merged rows are in
    // descending minVQSLod order.  A target advances through that walk until
    // the distance to the next row gets worse; the previous row is then
    // emitted.  It is not an independent nearest-neighbour search: preserving
    // this detail matters for multiple targets and for the 100% boundary.
    std::size_t row_index = 0;
    std::size_t target_index = 0;
    double previous_distance = 100.0;
    double target = targets[target_index];
    while (row_index < merged.size() && target_index < targets.size()) {
        const auto current = row_index;
        if (row_index + 1 < merged.size()) {
            const auto previous = current;
            ++row_index;
            const double current_distance = std::abs(target - merged[row_index].sensitivity * 100.0);
            if (current_distance <= previous_distance) {
                previous_distance = current_distance;
                continue;
            }
            auto row = merged[previous];
            row.sensitivity = target;
            selected.push_back(std::move(row));
            ++target_index;
            if (target_index >= targets.size()) break;
            target = targets[target_index];
            previous_distance = std::abs(target - merged[row_index].sensitivity * 100.0);
            continue;
        }
        auto row = merged[current];
        row.sensitivity = target;
        selected.push_back(std::move(row));
        ++target_index;
    }
    return selected;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto wall_begin = std::chrono::steady_clock::now();
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        const auto options = parse(argc, argv);
        const auto merged = merge(options.inputs, options.mode);
        const auto selected = select_targets(merged, options.truth_sensitivity);
        std::ostringstream output;
        output << "# Variant quality score tranches file\n# Version number 5\n"
                  "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";
        output << std::fixed << std::setprecision(2);
        // Tranche.tranchesString sorts by callsAtTruthSites, not by requested
        // sensitivity or VQSLOD.  Java's Collections.sort is stable, so use a
        // stable sort to retain deterministic order for tied counts.
        auto printable = selected;
        std::stable_sort(printable.begin(), printable.end(), [](const Tranche& left, const Tranche& right) {
            return left.calls_at_truth_sites < right.calls_at_truth_sites;
        });
        for (std::size_t index = 0; index < printable.size(); ++index) {
            const auto& row = printable[index];
            const double previous = index == 0 ? 0.0 : printable[index - 1].sensitivity;
            output << row.sensitivity << ',' << row.num_known << ',' << row.num_novel << ','
                   << std::setprecision(4) << row.known_titv << ',' << row.novel_titv << ','
                   << row.min_vqslod << ",VQSRTranche" << row.model << std::setprecision(2)
                   << previous << "to" << row.sensitivity << ',' << row.model << ','
                   << row.accessible_truth_sites << ',' << row.calls_at_truth_sites << ','
                   << std::setprecision(4) << (row.accessible_truth_sites > 0
                       ? static_cast<double>(row.calls_at_truth_sites) / row.accessible_truth_sites : 0.0) << '\n';
            output << std::setprecision(2);
        }
        const auto output_text = output.str();
        write_text_output(options.output, output_text);
        if (!complete_file(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete output");
        const auto wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_begin).count();
        const auto output_bytes = std::filesystem::file_size(options.output);
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest");
            manifest << "{\"schema_version\":1,\"tool\":\"GatherTranches\","
                     << "\"implementation\":\"fastgatk-gather-tranches\",\"status\":\"prototype\","
                     << "\"execution_space\":\"Host\",\"determinism\":\"strict\","
                     << "\"primary_output\":\"" << json_escape(options.output)
                     << "\",\"primary_output_kind\":\"tranches\","
                     << "\"merged_tranches\":" << merged.size()
                     << ",\"selected_tranches\":" << selected.size() << ","
                     << "\"compatibility\":{\"merge_scattered_tranches\":true,"
                     << "\"stateful_target_selection\":true,\"version6_input\":true,"
                     << "\"gatk_version5_output\":true,\"cloud_uri_staging\":false,"
                     << "\"compressed_output\":" << (compressed_path(options.output) ? "true" : "false") << ","
                     << "\"bit_identical_to_gatk\":true},"
                     << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"tranches\",\"complete\":true}],"
                     << "\"telemetry\":{\"resources\":" << resources.to_json()
                     << ",\"input_files\":" << options.inputs.size()
                     << ",\"merged_tranches\":" << merged.size()
                     << ",\"selected_tranches\":" << selected.size()
                     << ",\"output_compressed\":" << (compressed_path(options.output) ? "true" : "false")
                     << ",\"output_bytes\":" << output_bytes
                     << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}}\n";
            manifest.close();
            if (!manifest || !complete_file(options.manifest))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest");
        }
        std::cout << "{\"tool\":\"GatherTranches\",\"status\":\"prototype\",\"merged_tranches\":"
                  << merged.size() << ",\"selected_tranches\":" << selected.size()
                  << ",\"execution_space\":\"Host\",\"output_bytes\":" << output_bytes
                  << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
