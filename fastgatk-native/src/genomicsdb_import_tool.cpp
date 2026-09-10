#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/vcf.h>
#endif

#if defined(__unix__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {

struct Options {
    std::vector<std::string> forwarded;
    std::vector<std::string> input_paths;
    std::vector<std::string> sample_name_maps;
    std::string workspace;
    bool update_workspace = false;
    bool workspace_path_seen = false;
    bool update_workspace_path_seen = false;
    std::string manifest;
    std::string tmp_dir;
    std::uint32_t requested_batch = 0;
    std::uint32_t effective_batch = 0;
    std::size_t requested_reader_threads = 0;
    std::size_t effective_reader_threads = 0;
    bool has_batch = false;
    bool has_reader_threads = false;
    bool has_tmp_dir = false;
    bool has_input = false;
    bool overwrite = false;
    bool native_workspace = false;
    bool resume_native_workspace = false;
    bool checkpoint_resumed = false;
    std::size_t checkpoint_inputs_reused = 0;
    // GATK's sample-name-map is observable in callset.json and in the
    // downstream GenotypeGVCFs sample columns.  Keep the mapping with the
    // normalized input path so the native sparse workspace can materialize a
    // private header with the requested sample name.
    std::map<std::string, std::string> sample_name_by_input;
    std::size_t sample_name_map_entries = 0;
    // GATK GenomicsDBImport restricts the imported array to records whose
    // reference span overlaps the include interval set.  Keep the raw
    // selectors so interval-list/BED files can be resolved against each
    // input header (shards may have disjoint dictionaries).
    std::vector<std::string> interval_selectors;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::uint64_t interval_records_kept = 0;
    std::uint64_t interval_records_skipped = 0;
};

std::string env(const char* name) {
    const char* value = std::getenv(name);
    return value == nullptr ? std::string{} : std::string(value);
}

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
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

std::uint64_t parse_uint(const std::string& value, const char* option) {
    if (value.empty()) throw std::invalid_argument(std::string("missing value for ") + option);
    std::size_t consumed = 0;
    const auto parsed = std::stoull(value, &consumed);
    if (consumed != value.size() || parsed == 0 || parsed > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(std::string("invalid positive integer for ") + option + ": " + value);
    return parsed;
}

std::string require_next(int& index, int argc, char** argv,
                         const std::string& argument, const char* long_name) {
    const auto inline_value = option_value(argument, long_name);
    if (!inline_value.empty()) return inline_value;
    if (argument == long_name && index + 1 < argc) return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

void replace_or_append(std::vector<std::string>& arguments, const char* long_name,
                       const std::string& value, bool had_option) {
    if (!had_option) {
        arguments.emplace_back(long_name);
        arguments.push_back(value);
        return;
    }
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        if (!is_option(arguments[index], long_name)) continue;
        if (!option_value(arguments[index], long_name).empty()) {
            arguments[index] = std::string(long_name) + "=" + value;
        } else if (index + 1 < arguments.size()) {
            arguments[index + 1] = value;
        }
        return;
    }
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-genomicsdb-import (resource-aware external GenomicsDB adapter)\n"
                         "  -V, --variant FILE                         input GVCF (repeatable)\n"
                         "  --sample-name-map FILE                     sample map input\n"
                         "  -L, --intervals REGION                    imported span selector (repeatable)\n"
                         "      --interval-set-rule RULE              UNION (default) or INTERSECTION\n"
                         "  --genomicsdb-workspace-path DIRECTORY      output workspace\n"
                         "  --genomicsdb-update-workspace-path DIR    update an existing workspace\n"
                         "  --batch-size INTEGER                       reader batch cap\n"
                         "  --reader-threads INTEGER                   reader thread cap\n"
                         "  --fastgatk-native-workspace                create a native sparse workspace index\n"
                         "  --resume-native-workspace                  resume a native workspace input checkpoint\n"
                         "  --output-manifest FILE                     adapter manifest\n";
            std::exit(0);
        }
        if (is_option(argument, "--genomicsdb-workspace-path")) {
            if (options.update_workspace_path_seen)
                throw std::invalid_argument(
                    "--genomicsdb-workspace-path cannot be combined with --genomicsdb-update-workspace-path");
            options.workspace = require_next(index, argc, argv, argument, "--genomicsdb-workspace-path");
            options.update_workspace = false;
            options.workspace_path_seen = true;
            options.forwarded.push_back(argument);
            if (option_value(argument, "--genomicsdb-workspace-path").empty())
                options.forwarded.push_back(options.workspace);
            continue;
        }
        if (is_option(argument, "--genomicsdb-update-workspace-path")) {
            if (options.workspace_path_seen)
                throw std::invalid_argument(
                    "--genomicsdb-update-workspace-path cannot be combined with --genomicsdb-workspace-path");
            options.workspace = require_next(index, argc, argv, argument,
                                             "--genomicsdb-update-workspace-path");
            options.update_workspace = true;
            options.update_workspace_path_seen = true;
            options.forwarded.push_back(argument);
            if (option_value(argument, "--genomicsdb-update-workspace-path").empty())
                options.forwarded.push_back(options.workspace);
            continue;
        }
        if (argument == "-V" || argument == "--variant" || is_option(argument, "--variant")) {
            const auto value = argument == "-V"
                ? (index + 1 < argc ? std::string(argv[++index]) :
                   throw std::invalid_argument("missing value for --variant"))
                : require_next(index, argc, argv, argument, "--variant");
            options.has_input = true;
            options.input_paths.push_back(value);
            options.forwarded.push_back(argument);
            if (option_value(argument, "--variant").empty()) options.forwarded.push_back(value);
            continue;
        }
        if (argument == "--sample-name-map" || is_option(argument, "--sample-name-map")) {
            const auto value = require_next(index, argc, argv, argument, "--sample-name-map");
            options.has_input = true;
            options.sample_name_maps.push_back(value);
            options.forwarded.push_back(argument);
            if (option_value(argument, "--sample-name-map").empty()) options.forwarded.push_back(value);
            continue;
        }
        if (argument == "-L" || is_option(argument, "--intervals") ||
            is_option(argument, "--interval") || is_option(argument, "--region")) {
            const char* option_name = argument == "-L" ? "--intervals" :
                argument.rfind("--intervals", 0) == 0 ? "--intervals" :
                argument.rfind("--interval", 0) == 0 ? "--interval" : "--region";
            const auto value = argument == "-L"
                ? (index + 1 < argc ? std::string(argv[++index])
                                    : throw std::invalid_argument("missing value for --intervals"))
                : require_next(index, argc, argv, argument, option_name);
            options.interval_selectors.push_back(value);
            options.forwarded.push_back(argument);
            if (option_value(argument, option_name).empty()) options.forwarded.push_back(value);
            continue;
        }
        if (is_option(argument, "--interval-set-rule")) {
            const auto value = require_next(index, argc, argv, argument, "--interval-set-rule");
            std::string normalized = value;
            std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (normalized == "UNION") {
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            } else if (normalized == "INTERSECTION") {
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            } else {
                throw std::invalid_argument("BAD_INPUT: --interval-set-rule must be UNION or INTERSECTION");
            }
            options.forwarded.push_back(argument);
            if (option_value(argument, "--interval-set-rule").empty()) options.forwarded.push_back(value);
            continue;
        }
        if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest")) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = require_next(index, argc, argv, argument, name);
            continue;
        }
        if (is_option(argument, "--batch-size")) {
            const auto value = require_next(index, argc, argv, argument, "--batch-size");
            options.requested_batch = static_cast<std::uint32_t>(parse_uint(value, "--batch-size"));
            options.has_batch = true;
            options.forwarded.push_back(argument);
            if (option_value(argument, "--batch-size").empty()) options.forwarded.push_back(value);
            continue;
        }
        if (is_option(argument, "--reader-threads")) {
            const auto value = require_next(index, argc, argv, argument, "--reader-threads");
            options.requested_reader_threads = static_cast<std::size_t>(parse_uint(value, "--reader-threads"));
            options.has_reader_threads = true;
            options.forwarded.push_back(argument);
            if (option_value(argument, "--reader-threads").empty()) options.forwarded.push_back(value);
            continue;
        }
        if (is_option(argument, "--tmp-dir")) {
            options.tmp_dir = require_next(index, argc, argv, argument, "--tmp-dir");
            options.has_tmp_dir = true;
            options.forwarded.push_back(argument);
            if (option_value(argument, "--tmp-dir").empty()) options.forwarded.push_back(options.tmp_dir);
            continue;
        }
        if (argument == "--fastgatk-native-workspace") {
            options.native_workspace = true;
            continue;
        }
        if (argument == "--resume-native-workspace") {
            options.native_workspace = true;
            options.resume_native_workspace = true;
            continue;
        }
        if (argument == "--overwrite-existing-genomicsdb-workspace") options.overwrite = true;
        if (argument == "-V" || argument == "--variant" || is_option(argument, "--variant") ||
            argument == "--sample-name-map" || is_option(argument, "--sample-name-map"))
            options.has_input = true;
        // Keep unsupported/release-specific arguments verbatim for the
        // external backend; direct replacement depends on argv preservation.
        options.forwarded.push_back(argument);
    }
    if (options.workspace.empty())
        throw std::invalid_argument(
            "BAD_INPUT: --genomicsdb-workspace-path or --genomicsdb-update-workspace-path is required");
    if (!options.has_input)
        throw std::invalid_argument("BAD_INPUT: at least one -V/--variant or --sample-name-map is required");
    if (options.resume_native_workspace && options.overwrite)
        throw std::invalid_argument("--resume-native-workspace cannot be combined with --overwrite-existing-genomicsdb-workspace");
    if (options.update_workspace && options.overwrite)
        throw std::invalid_argument(
            "--genomicsdb-update-workspace-path cannot be combined with --overwrite-existing-genomicsdb-workspace");
    return options;
}

std::uint32_t choose_batch(const Options& options,
                           const fastgatk::runtime::ResourceSnapshot& resources) {
    // GenomicsDB keeps one reader and bookkeeping state per sample in a batch.
    // Reserve a fixed control-plane allowance, then cap the backend batch at a
    // conservative 16 MiB per reader. Unknown limits leave an explicit user
    // value untouched and otherwise use GATK's own default.
    const auto budget = resources.safe_memory_budget_bytes();
    if (budget == 0) return options.has_batch ? options.requested_batch : 0;
    constexpr std::uint64_t control_plane = 8ULL * 1024ULL * 1024ULL;
    constexpr std::uint64_t bytes_per_reader = 16ULL * 1024ULL * 1024ULL;
    const auto usable = budget > control_plane ? budget - control_plane : bytes_per_reader;
    const auto computed = static_cast<std::uint32_t>(std::clamp<std::uint64_t>(
        usable / bytes_per_reader, 1, 512));
    return options.has_batch ? std::min(options.requested_batch, computed) : computed;
}

std::vector<std::string> build_backend_argv(Options& options,
                                             const fastgatk::runtime::ResourceSnapshot& resources) {
    options.effective_batch = choose_batch(options, resources);
    if (options.effective_batch != 0)
        replace_or_append(options.forwarded, "--batch-size",
                          std::to_string(options.effective_batch), options.has_batch);
    options.effective_reader_threads = resources.effective_threads(
        options.has_reader_threads ? options.requested_reader_threads : 0);
    if (options.has_reader_threads)
        replace_or_append(options.forwarded, "--reader-threads",
                          std::to_string(options.effective_reader_threads), true);
    else if (options.effective_reader_threads > 1)
        replace_or_append(options.forwarded, "--reader-threads",
                          std::to_string(options.effective_reader_threads), false);
    if (!options.has_tmp_dir && !resources.scratch_directory.empty()) {
        options.tmp_dir = resources.scratch_directory;
        options.has_tmp_dir = true;
        replace_or_append(options.forwarded, "--tmp-dir", options.tmp_dir, false);
    }
    const auto gatk = env("FASTGATK_GATK_BINARY").empty() ? std::string("gatk") : env("FASTGATK_GATK_BINARY");
    std::vector<std::string> command{gatk, "GenomicsDBImport"};
    command.insert(command.end(), options.forwarded.begin(), options.forwarded.end());
    return command;
}

#if defined(__unix__) || defined(__APPLE__)
int run_backend(const std::vector<std::string>& command) {
    std::vector<char*> argv;
    argv.reserve(command.size() + 1);
    for (const auto& value : command) argv.push_back(const_cast<char*>(value.c_str()));
    argv.push_back(nullptr);
    const auto pid = fork();
    if (pid < 0) throw std::runtime_error("BACKEND_UNAVAILABLE: cannot fork GenomicsDB backend");
    if (pid == 0) {
        execvp(argv[0], argv.data());
        _exit(127);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0)
        throw std::runtime_error("BACKEND_UNAVAILABLE: cannot wait for GenomicsDB backend");
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    if (WIFSIGNALED(status)) return 128 + WTERMSIG(status);
    return 1;
}
#else
int run_backend(const std::vector<std::string>&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: external GenomicsDB adapter requires POSIX process support");
}
#endif

bool workspace_complete(const std::string& path) {
    std::error_code error;
    const std::filesystem::path workspace(path);
    if (!std::filesystem::is_directory(workspace, error) || error) return false;
    for (const auto& entry : std::filesystem::directory_iterator(workspace, error)) {
        if (error) return false;
        if (entry.is_regular_file(error) && !error && entry.file_size(error) > 0 && !error) return true;
    }
    return false;
}

std::filesystem::path native_workspace_checkpoint_path(const std::string& workspace) {
    return std::filesystem::path(workspace) / "fastgatk-native-workspace.checkpoint";
}

std::uintmax_t input_size_or_zero(const std::string& path) {
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : size;
}

std::intmax_t input_mtime_or_zero(const std::string& path) {
    std::error_code error;
    const auto stamp = std::filesystem::last_write_time(path, error);
    return error ? 0 : static_cast<std::intmax_t>(stamp.time_since_epoch().count());
}

struct NativeWorkspaceCheckpoint {
    std::vector<std::string> inputs;
    std::vector<std::uintmax_t> sizes;
    std::vector<std::intmax_t> mtimes;
    std::vector<bool> completed;
    std::map<std::string, std::string> sample_name_by_input;
    std::vector<std::string> interval_selectors;
    std::string interval_set_rule = "UNION";
};

NativeWorkspaceCheckpoint make_native_checkpoint(
    const std::vector<std::string>& inputs,
    const std::map<std::string, std::string>& sample_name_by_input,
    const std::vector<std::string>& interval_selectors,
    fastgatk::io::HtsIntervalSetRule interval_set_rule) {
    NativeWorkspaceCheckpoint checkpoint;
    checkpoint.inputs = inputs;
    checkpoint.sample_name_by_input = sample_name_by_input;
    checkpoint.interval_selectors = interval_selectors;
    checkpoint.interval_set_rule = interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
    checkpoint.sizes.reserve(inputs.size());
    checkpoint.mtimes.reserve(inputs.size());
    checkpoint.completed.assign(inputs.size(), false);
    for (const auto& input : inputs) {
        const auto size = input_size_or_zero(input);
        if (size == 0)
            throw std::runtime_error("BAD_INPUT: cannot fingerprint native GenomicsDB input: " + input);
        checkpoint.sizes.push_back(size);
        checkpoint.mtimes.push_back(input_mtime_or_zero(input));
    }
    return checkpoint;
}

void write_native_workspace_checkpoint(const std::filesystem::path& path,
                                       const NativeWorkspaceCheckpoint& checkpoint) {
    if (checkpoint.inputs.size() != checkpoint.sizes.size() ||
        checkpoint.inputs.size() != checkpoint.mtimes.size() ||
        checkpoint.inputs.size() != checkpoint.completed.size())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: malformed native GenomicsDB checkpoint state");
    const auto temporary = path.string() + ".tmp";
    {
        std::ofstream stream(temporary, std::ios::trunc);
        if (!stream)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot write native GenomicsDB checkpoint: " + temporary);
        stream << "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V3\n"
               << "input_count " << checkpoint.inputs.size() << "\n";
        for (std::size_t index = 0; index < checkpoint.inputs.size(); ++index) {
            stream << "input " << std::quoted(checkpoint.inputs[index]) << ' '
                   << checkpoint.sizes[index] << ' ' << checkpoint.mtimes[index] << ' '
                   << (checkpoint.completed[index] ? 1 : 0) << "\n";
        }
        stream << "sample_map_count " << checkpoint.sample_name_by_input.size() << "\n";
        for (const auto& mapping : checkpoint.sample_name_by_input)
            stream << "sample_map " << std::quoted(mapping.first) << ' '
                   << std::quoted(mapping.second) << "\n";
        stream << "interval_set_rule " << checkpoint.interval_set_rule << "\n"
               << "interval_count " << checkpoint.interval_selectors.size() << "\n";
        for (const auto& selector : checkpoint.interval_selectors)
            stream << "interval " << std::quoted(selector) << "\n";
        stream.flush();
        if (!stream)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot finalize native GenomicsDB checkpoint: " + temporary);
    }
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) {
        std::filesystem::remove(path, error);
        error.clear();
        std::filesystem::rename(temporary, path, error);
    }
    if (error)
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot install native GenomicsDB checkpoint: " + path.string());
}

NativeWorkspaceCheckpoint read_native_workspace_checkpoint(
    const std::filesystem::path& path, const std::vector<std::string>& expected_inputs,
    const std::vector<std::string>& expected_intervals = {},
    const std::string& expected_interval_set_rule = "UNION") {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("BAD_INPUT: cannot open native GenomicsDB checkpoint: " + path.string());
    std::string version;
    if (!(stream >> version) ||
        (version != "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V1" &&
         version != "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V2" &&
         version != "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V3"))
        throw std::runtime_error("BAD_INPUT: unsupported native GenomicsDB checkpoint version");
    const bool has_sample_map = version != "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V1";
    const bool has_intervals = version == "FASTGATK_NATIVE_GENOMICSDB_CHECKPOINT_V3";
    std::string key;
    std::size_t input_count = 0;
    std::size_t sample_map_count = 0;
    std::size_t interval_count = 0;
    NativeWorkspaceCheckpoint checkpoint;
    while (stream >> key) {
        if (key == "input_count") {
            stream >> input_count;
            if (!stream || input_count > (1ULL << 20))
                throw std::runtime_error("BAD_INPUT: invalid native GenomicsDB checkpoint input count");
            checkpoint.inputs.reserve(input_count);
            checkpoint.sizes.reserve(input_count);
            checkpoint.mtimes.reserve(input_count);
            checkpoint.completed.reserve(input_count);
        } else if (key == "input") {
            std::string value;
            std::uintmax_t size = 0;
            std::intmax_t mtime = 0;
            int complete = 0;
            stream >> std::quoted(value) >> size >> mtime >> complete;
            if (!stream || (complete != 0 && complete != 1))
                throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB checkpoint input");
            checkpoint.inputs.push_back(std::move(value));
            checkpoint.sizes.push_back(size);
            checkpoint.mtimes.push_back(mtime);
            checkpoint.completed.push_back(complete != 0);
        } else if (key == "sample_map") {
            if (!has_sample_map)
                throw std::runtime_error("BAD_INPUT: sample-map field is not supported by V1 checkpoint");
            std::string input;
            std::string sample;
            stream >> std::quoted(input) >> std::quoted(sample);
            if (!stream || input.empty() || sample.empty() ||
                !checkpoint.sample_name_by_input.emplace(std::move(input), std::move(sample)).second)
                throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB checkpoint sample map");
        } else if (key == "sample_map_count") {
            if (!has_sample_map)
                throw std::runtime_error("BAD_INPUT: sample-map count is not supported by V1 checkpoint");
            stream >> sample_map_count;
            if (!stream || sample_map_count > (1ULL << 20))
                throw std::runtime_error("BAD_INPUT: invalid native GenomicsDB checkpoint sample-map count");
        } else if (key == "interval_set_rule") {
            if (!has_intervals)
                throw std::runtime_error("BAD_INPUT: interval fields are not supported by legacy checkpoint");
            stream >> checkpoint.interval_set_rule;
            if (!stream || (checkpoint.interval_set_rule != "UNION" &&
                            checkpoint.interval_set_rule != "INTERSECTION"))
                throw std::runtime_error("BAD_INPUT: invalid native GenomicsDB checkpoint interval set rule");
        } else if (key == "interval_count") {
            if (!has_intervals)
                throw std::runtime_error("BAD_INPUT: interval fields are not supported by legacy checkpoint");
            stream >> interval_count;
            if (!stream || interval_count > (1ULL << 20))
                throw std::runtime_error("BAD_INPUT: invalid native GenomicsDB checkpoint interval count");
            checkpoint.interval_selectors.reserve(interval_count);
        } else if (key == "interval") {
            if (!has_intervals)
                throw std::runtime_error("BAD_INPUT: interval fields are not supported by legacy checkpoint");
            std::string selector;
            stream >> std::quoted(selector);
            if (!stream || selector.empty())
                throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB checkpoint interval");
            checkpoint.interval_selectors.push_back(std::move(selector));
        } else {
            throw std::runtime_error("BAD_INPUT: unknown native GenomicsDB checkpoint field: " + key);
        }
    }
    if (input_count == 0 || checkpoint.inputs.size() != input_count ||
        checkpoint.sizes.size() != input_count || checkpoint.mtimes.size() != input_count ||
        checkpoint.completed.size() != input_count)
        throw std::runtime_error("BAD_INPUT: incomplete native GenomicsDB checkpoint");
    if (has_sample_map && checkpoint.sample_name_by_input.size() != sample_map_count)
        throw std::runtime_error("BAD_INPUT: incomplete native GenomicsDB checkpoint sample map");
    if (has_intervals && checkpoint.interval_selectors.size() != interval_count)
        throw std::runtime_error("BAD_INPUT: incomplete native GenomicsDB checkpoint intervals");
    if (!has_intervals && (!expected_intervals.empty() || expected_interval_set_rule != "UNION"))
        throw std::runtime_error("BAD_INPUT: legacy native GenomicsDB checkpoint lacks interval selection metadata");
    if (has_intervals && (checkpoint.interval_selectors != expected_intervals ||
                          checkpoint.interval_set_rule != expected_interval_set_rule))
        throw std::runtime_error("BAD_INPUT: native GenomicsDB checkpoint interval selection mismatch");
    if (checkpoint.inputs != expected_inputs)
        throw std::runtime_error("BAD_INPUT: native GenomicsDB checkpoint input set mismatch");
    for (std::size_t index = 0; index < checkpoint.inputs.size(); ++index) {
        if (input_size_or_zero(checkpoint.inputs[index]) != checkpoint.sizes[index] ||
            input_mtime_or_zero(checkpoint.inputs[index]) != checkpoint.mtimes[index])
            throw std::runtime_error("BAD_INPUT: native GenomicsDB input changed since checkpoint: " +
                                     checkpoint.inputs[index]);
    }
    for (const auto& mapping : checkpoint.sample_name_by_input) {
        if (std::find(checkpoint.inputs.begin(), checkpoint.inputs.end(), mapping.first) == checkpoint.inputs.end())
            throw std::runtime_error("BAD_INPUT: checkpoint sample-map input is not in input set: " + mapping.first);
    }
    return checkpoint;
}

std::size_t native_checkpoint_completed_count(const NativeWorkspaceCheckpoint& checkpoint) {
    return static_cast<std::size_t>(std::count(
        checkpoint.completed.begin(), checkpoint.completed.end(), true));
}

bool is_remote_uri(const std::string& path) {
    return path.rfind("http://", 0) == 0 || path.rfind("https://", 0) == 0 ||
           path.rfind("s3://", 0) == 0 || path.rfind("gs://", 0) == 0 ||
           path.rfind("az://", 0) == 0 || path.rfind("file://", 0) == 0;
}

// A GenomicsDB workspace is a directory tree.  The non-empty preflight alone
// does not prevent two retried SLURM/Nextflow tasks from writing it at the same
// time, so take an atomic sibling-directory lock before inspecting or
// publishing any workspace files.  `create_directory` is atomic on the local
// filesystems used by scratch and is portable across POSIX and Windows.
class WorkspaceLock {
public:
    explicit WorkspaceLock(const std::string& workspace)
        : path_(workspace + ".fastgatk.lock") {
        std::error_code error;
        if (!std::filesystem::create_directory(path_, error)) {
            if (error) {
                throw std::runtime_error(
                    "RESOURCE_EXHAUSTED: cannot acquire GenomicsDB workspace lock " +
                    path_.string() + ": " + error.message());
            }
            throw std::runtime_error(
                "RESOURCE_EXHAUSTED: GenomicsDB workspace is locked by another task: " +
                path_.string() +
                "; wait for the active writer or remove the lock only after validating that it is stale");
        }
        owns_ = true;
        std::ofstream owner(path_ / "owner", std::ios::trunc);
        if (!owner) {
            release();
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: cannot write GenomicsDB workspace lock owner: " +
                path_.string());
        }
#if defined(__unix__) || defined(__APPLE__)
        owner << "pid=" << static_cast<long long>(::getpid()) << '\n';
#else
        owner << "pid=unknown\n";
#endif
        owner << "thread=" << std::this_thread::get_id() << '\n';
        owner.flush();
        if (!owner) {
            release();
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: cannot finalize GenomicsDB workspace lock owner: " +
                path_.string());
        }
    }

    WorkspaceLock(const WorkspaceLock&) = delete;
    WorkspaceLock& operator=(const WorkspaceLock&) = delete;

    ~WorkspaceLock() { release(); }

    const std::filesystem::path& path() const noexcept { return path_; }

private:
    void release() noexcept {
        if (!owns_) return;
        std::error_code error;
        std::filesystem::remove_all(path_, error);
        owns_ = false;
    }

    std::filesystem::path path_;
    bool owns_ = false;
};

void validate_input_paths(const Options& options) {
    for (const auto& path : options.input_paths) {
        if (is_remote_uri(path)) continue;
        std::error_code error;
        if (!std::filesystem::is_regular_file(path, error) || error)
            throw std::runtime_error("BAD_INPUT: GenomicsDB variant input is missing: " + path);
    }
    for (const auto& path : options.sample_name_maps) {
        std::error_code error;
        if (!std::filesystem::is_regular_file(path, error) || error)
            throw std::runtime_error("BAD_INPUT: GenomicsDB sample-name-map is missing: " + path);
    }
}

void validate_workspace_preflight(const Options& options) {
    std::error_code error;
    const auto workspace = std::filesystem::path(options.workspace);
    const auto exists = std::filesystem::exists(workspace, error);
    if (error)
        throw std::runtime_error("BAD_INPUT: cannot inspect GenomicsDB workspace: " + options.workspace);
    if (!exists) {
        if (options.resume_native_workspace)
            throw std::runtime_error(
                "BAD_INPUT: --resume-native-workspace requested but native GenomicsDB workspace is missing: " +
                options.workspace);
        if (options.update_workspace)
            throw std::runtime_error(
                "BAD_INPUT: --genomicsdb-update-workspace-path requires an existing GenomicsDB workspace: " +
                options.workspace);
        return;
    }
    if (!std::filesystem::is_directory(workspace, error) || error)
        throw std::runtime_error("BAD_INPUT: GenomicsDB workspace path is not a directory: " + options.workspace);
    bool nonempty = false;
    for (const auto& entry : std::filesystem::directory_iterator(workspace, error)) {
        (void)entry;
        if (error) throw std::runtime_error("BAD_INPUT: cannot inspect GenomicsDB workspace: " + options.workspace);
        nonempty = true;
        break;
    }
    const auto checkpoint = native_workspace_checkpoint_path(options.workspace);
    std::error_code checkpoint_error;
    const bool resumable = options.resume_native_workspace &&
        std::filesystem::is_regular_file(checkpoint, checkpoint_error) && !checkpoint_error;
    if (options.resume_native_workspace && !resumable)
        throw std::runtime_error(
            "BAD_INPUT: --resume-native-workspace requested but no native GenomicsDB checkpoint exists: " +
            checkpoint.string());
    // GATK's incremental mode intentionally accepts the populated workspace
    // and lets GenomicsDB validate its schema/callset before appending.  It
    // also ignores this invocation's intervals in favor of the initial
    // import's interval set; preserve those semantics in the external mode.
    if (options.update_workspace)
        return;
    if (nonempty && !options.overwrite && !resumable)
        throw std::runtime_error(
            "BAD_INPUT: GenomicsDB workspace is non-empty; pass --overwrite-existing-genomicsdb-workspace: " +
            options.workspace);
}

std::string normalize_workspace_input_path(std::string input) {
    if (is_remote_uri(input)) return input;
    std::error_code error;
    const auto absolute = std::filesystem::absolute(input, error);
    return error ? input : absolute.string();
}

std::vector<std::string> collect_workspace_inputs(Options& options) {
    std::vector<std::string> inputs = options.input_paths;
    std::set<std::string> sample_names;
    for (const auto& map_path : options.sample_name_maps) {
        std::ifstream stream(map_path);
        if (!stream) throw std::runtime_error("BAD_INPUT: cannot read --sample-name-map: " + map_path);
        std::string line;
        while (std::getline(stream, line)) {
            const auto first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#') continue;
            std::istringstream fields(line);
            std::string sample;
            std::string path;
            std::string extra;
            fields >> sample >> path >> extra;
            if (sample.empty() || path.empty() || !extra.empty())
                throw std::runtime_error(
                    "BAD_INPUT: --sample-name-map requires exactly two columns (sample and path): " +
                    map_path);
            if (!sample_names.insert(sample).second)
                throw std::runtime_error("BAD_INPUT: duplicate sample name in --sample-name-map: " + sample);
            const auto normalized = normalize_workspace_input_path(path);
            const auto existing = options.sample_name_by_input.find(normalized);
            if (existing != options.sample_name_by_input.end() && existing->second != sample)
                throw std::runtime_error(
                    "BAD_INPUT: one GenomicsDB input is assigned multiple sample names: " + normalized);
            options.sample_name_by_input[normalized] = sample;
            ++options.sample_name_map_entries;
            inputs.push_back(path);
        }
    }
    // GenomicsDB callset/sample order is observable downstream (the merged
    // VCF header and genotype columns follow import order).  Normalize local
    // paths without sorting: GATK preserves -V/sample-map order, so sorting
    // by filename would silently change a valid workspace's sample columns.
    std::vector<std::string> ordered;
    ordered.reserve(inputs.size());
    std::set<std::string> seen;
    for (auto& input : inputs) input = normalize_workspace_input_path(input);
    for (const auto& input : inputs) {
        if (seen.insert(input).second) ordered.push_back(input);
    }
    return ordered;
}

void write_input_index(const std::filesystem::path& target,
                       const std::vector<std::string>& inputs) {
    const auto temporary = target.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::trunc);
    if (!stream) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write workspace input index");
    stream << "# fastgatk GenomicsDB adapter input index v1\n";
    for (const auto& input : inputs) stream << input << '\n';
    stream.flush();
    if (!stream) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize workspace input index");
    std::error_code error;
    std::filesystem::rename(temporary, target, error);
    if (error) {
        std::filesystem::remove(target, error);
        error.clear();
        std::filesystem::rename(temporary, target, error);
    }
    if (error) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish workspace input index");
}

#if FASTGATK_HAS_HTSLIB
struct RecordSpan {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::uint64_t records = 0;
};

using NativeImportIntervals = std::vector<fastgatk::io::IndexedInterval>;

NativeImportIntervals parse_native_import_intervals(
    const bcf_hdr_t* header, const std::vector<std::string>& selectors) {
    NativeImportIntervals intervals;
    fastgatk::io::IntervalFileStats file_stats;
    for (const auto& selector : selectors) {
        // Different input shards are allowed to have disjoint dictionaries.
        // Keep an unknown contig by name here; it simply has no overlap with
        // this shard, while rejecting malformed coordinates still happens in
        // the shared Host parser.
        fastgatk::io::append_interval_selector(
            selector, header, intervals, file_stats, true);
    }
    return intervals;
}

std::int64_t native_record_end(const bcf1_t* record) {
    const auto start = static_cast<std::int64_t>(record->pos);
    auto end = start + std::max<std::int64_t>(1, static_cast<std::int64_t>(record->rlen));
    return end;
}

bool native_record_overlaps_import_intervals(
    const bcf_hdr_t* header, const bcf1_t* record,
    const NativeImportIntervals& intervals,
    fastgatk::io::HtsIntervalSetRule interval_set_rule) {
    if (intervals.empty()) return true;
    const auto* contig = bcf_hdr_id2name(header, record->rid);
    if (contig == nullptr) return false;
    const auto start = static_cast<std::int64_t>(record->pos);
    auto end = native_record_end(record);
    int32_t* end_value = nullptr;
    int end_count = 0;
    const auto end_length = bcf_get_info_int32(
        header, const_cast<bcf1_t*>(record), "END", &end_value, &end_count);
    if (end_length > 0 && end_count > 0 && end_value[0] > start + 1)
        end = static_cast<std::int64_t>(end_value[0]);
    free(end_value);
    const auto overlaps = [&](const fastgatk::io::IndexedInterval& interval) {
        return interval.contig == contig &&
               start < static_cast<std::int64_t>(interval.end) &&
               end > static_cast<std::int64_t>(interval.begin);
    };
    if (interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection)
        return std::all_of(intervals.begin(), intervals.end(), overlaps);
    return std::any_of(intervals.begin(), intervals.end(), overlaps);
}

std::vector<RecordSpan> scan_record_spans(const std::string& path) {
    htsFile* input = bcf_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot scan native GenomicsDB input: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    bcf1_t* record = bcf_init();
    if (!header || !record) {
        if (record) bcf_destroy(record);
        if (header) bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB input header: " + path);
    }
    std::map<std::string, RecordSpan> spans;
    int status = 0;
    while ((status = bcf_read(input, header, record)) == 0) {
        bcf_unpack(record, BCF_UN_STR | BCF_UN_INFO);
        const char* contig = bcf_hdr_id2name(header, record->rid);
        if (contig == nullptr) continue;
        auto& span = spans[contig];
        span.contig = contig;
        const auto start = static_cast<std::int64_t>(record->pos);
        auto end = start + 1;
        int32_t* end_value = nullptr;
        int end_count = 0;
        const auto end_length = bcf_get_info_int32(header, record, "END", &end_value, &end_count);
        if (end_length > 0 && end_count > 0 && end_value[0] > start + 1)
            end = static_cast<std::int64_t>(end_value[0]);
        free(end_value);
        if (span.records == 0) {
            span.start = start;
            span.end = end;
        } else {
            span.start = std::min(span.start, start);
            span.end = std::max(span.end, end);
        }
        ++span.records;
    }
    if (status < -1) {
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: failed scanning native GenomicsDB input: " + path);
    }
    const int close_status = bcf_close(input);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    if (close_status != 0)
        throw std::runtime_error("BAD_INPUT: failed closing native GenomicsDB input: " + path);
    std::vector<RecordSpan> result;
    result.reserve(spans.size());
    for (auto& entry : spans) result.push_back(std::move(entry.second));
    return result;
}
#endif

void write_record_index(const std::filesystem::path& target,
                        const std::vector<std::string>& inputs) {
#if FASTGATK_HAS_HTSLIB
    const auto temporary = target.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::trunc);
    if (!stream) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write native record index");
    stream << "# fastgatk native sparse record index v1\n"
           << "# input\\tcontig\\tstart0\\tend0\\trecords\n";
    for (const auto& input : inputs) {
        for (const auto& span : scan_record_spans(input))
            stream << input << '\t' << span.contig << '\t' << span.start << '\t'
                   << span.end << '\t' << span.records << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize native record index");
    std::error_code error;
    std::filesystem::rename(temporary, target, error);
    if (error) {
        std::filesystem::remove(target, error);
        error.clear();
        std::filesystem::rename(temporary, target, error);
    }
    if (error) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish native record index");
#else
    (void)target;
    (void)inputs;
    throw std::runtime_error(
        "BACKEND_UNAVAILABLE: native sparse record index requires HTSlib support");
#endif
}

void write_workspace_input_index(const std::string& workspace,
                                 const std::vector<std::string>& inputs) {
    write_input_index(std::filesystem::path(workspace) / "fastgatk-inputs.tsv", inputs);
}

std::vector<std::string> merge_workspace_input_index(
    const std::string& workspace, const std::vector<std::string>& inputs) {
    std::vector<std::string> merged;
    std::ifstream stream(std::filesystem::path(workspace) / "fastgatk-inputs.tsv");
    if (stream) {
        std::string line;
        while (std::getline(stream, line)) {
            const auto first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#') continue;
            const auto value = line.substr(first);
            if (!value.empty() && std::find(merged.begin(), merged.end(), value) == merged.end())
                merged.push_back(value);
        }
    }
    for (const auto& input : inputs) {
        if (std::find(merged.begin(), merged.end(), input) == merged.end())
            merged.push_back(input);
    }
    return merged;
}

// Native sparse workspaces do not have a TileDB callset map, but the metadata
// sidecar still records sample-name-map rewrites made during the initial
// import.  Preserve those private header rewrites when an incremental native
// update rebuilds the materialized index.  This is intentionally a tiny
// parser for the deterministic JSON emitted above; it is not exposed as a
// general JSON API and malformed entries fail closed at the call site.
std::string json_unescape_string(const std::string& value) {
    std::string result;
    result.reserve(value.size());
    bool escaped = false;
    for (const char character : value) {
        if (escaped) {
            switch (character) {
                case 'n': result.push_back('\n'); break;
                case 'r': result.push_back('\r'); break;
                case 't': result.push_back('\t'); break;
                case '\\': result.push_back('\\'); break;
                case '"': result.push_back('"'); break;
                default: result.push_back(character); break;
            }
            escaped = false;
        } else if (character == '\\') {
            escaped = true;
        } else {
            result.push_back(character);
        }
    }
    if (escaped) result.push_back('\\');
    return result;
}

std::optional<std::pair<std::string, std::size_t>> json_string_at(
    const std::string& text, std::size_t quote) {
    if (quote >= text.size() || text[quote] != '"') return std::nullopt;
    std::string encoded;
    bool escaped = false;
    for (std::size_t index = quote + 1; index < text.size(); ++index) {
        const auto character = text[index];
        if (!escaped && character == '"')
            return std::make_pair(json_unescape_string(encoded), index + 1);
        encoded.push_back(character);
        if (escaped) escaped = false;
        else if (character == '\\') escaped = true;
    }
    return std::nullopt;
}

std::map<std::string, std::string> read_native_workspace_sample_name_map(
    const std::string& workspace) {
    const auto metadata_path = std::filesystem::path(workspace) / "fastgatk-workspace.json";
    std::ifstream stream(metadata_path);
    if (!stream) return {};
    std::ostringstream text;
    text << stream.rdbuf();
    if (!stream.good() && !stream.eof())
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB workspace metadata: " +
                                 metadata_path.string());
    const auto source = text.str();
    const auto array_marker = source.find("\"sample_name_map\":[");
    if (array_marker == std::string::npos) return {};
    const auto array_end = source.find(']', array_marker);
    if (array_end == std::string::npos)
        throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map");
    std::map<std::string, std::string> result;
    std::size_t cursor = array_marker + std::string("\"sample_name_map\":[").size();
    while (cursor < array_end) {
        const auto input_marker = source.find("\"input\":", cursor);
        if (input_marker == std::string::npos || input_marker >= array_end) break;
        const auto input_quote = source.find('"', input_marker + std::string("\"input\":").size());
        if (input_quote == std::string::npos || input_quote >= array_end)
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map input");
        const auto input = json_string_at(source, input_quote);
        if (!input.has_value() || input->second > array_end)
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map input");
        const auto sample_marker = source.find("\"sample_name\":", input->second);
        if (sample_marker == std::string::npos || sample_marker >= array_end)
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map name");
        const auto sample_quote = source.find('"', sample_marker + std::string("\"sample_name\":").size());
        if (sample_quote == std::string::npos || sample_quote >= array_end)
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map name");
        const auto sample = json_string_at(source, sample_quote);
        if (!sample.has_value() || sample->second > array_end || input->first.empty() || sample->first.empty() ||
            !result.emplace(input->first, sample->first).second)
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace sample map entry");
        cursor = sample->second;
    }
    return result;
}

std::vector<std::string> read_native_workspace_interval_selectors(
    const std::string& workspace) {
    const auto metadata_path = std::filesystem::path(workspace) / "fastgatk-workspace.json";
    std::ifstream stream(metadata_path);
    if (!stream) return {};
    std::ostringstream text;
    text << stream.rdbuf();
    if (!stream.good() && !stream.eof())
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB workspace metadata: " +
                                 metadata_path.string());
    const auto source = text.str();
    const auto marker = source.find("\"import_intervals\":[");
    if (marker == std::string::npos) return {};
    const auto begin = marker + std::string("\"import_intervals\":[").size();
    const auto end = source.find(']', begin);
    if (end == std::string::npos)
        throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace intervals");
    std::vector<std::string> result;
    std::size_t cursor = begin;
    while (cursor < end) {
        while (cursor < end && (source[cursor] == ',' ||
                                std::isspace(static_cast<unsigned char>(source[cursor])))) ++cursor;
        if (cursor >= end) break;
        if (source[cursor] != '"')
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace interval selector");
        const auto parsed = json_string_at(source, cursor);
        if (!parsed.has_value() || parsed->second > end || parsed->first.empty())
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB workspace interval selector");
        result.push_back(parsed->first);
        cursor = parsed->second;
    }
    return result;
}

fastgatk::io::HtsIntervalSetRule read_native_workspace_interval_set_rule(
    const std::string& workspace) {
    const auto metadata_path = std::filesystem::path(workspace) / "fastgatk-workspace.json";
    std::ifstream stream(metadata_path);
    if (!stream) return fastgatk::io::HtsIntervalSetRule::Union;
    std::ostringstream text;
    text << stream.rdbuf();
    if (!stream.good() && !stream.eof())
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB workspace metadata: " +
                                 metadata_path.string());
    const auto source = text.str();
    const auto marker = source.find("\"interval_set_rule\":");
    if (marker == std::string::npos) return fastgatk::io::HtsIntervalSetRule::Union;
    const auto quote = source.find('"', marker + std::string("\"interval_set_rule\":").size());
    if (quote == std::string::npos) throw std::runtime_error(
        "BAD_INPUT: malformed native GenomicsDB workspace interval set rule");
    const auto parsed = json_string_at(source, quote);
    if (!parsed.has_value()) throw std::runtime_error(
        "BAD_INPUT: malformed native GenomicsDB workspace interval set rule");
    if (parsed->first == "INTERSECTION") return fastgatk::io::HtsIntervalSetRule::Intersection;
    if (parsed->first == "UNION") return fastgatk::io::HtsIntervalSetRule::Union;
    throw std::runtime_error("BAD_INPUT: unknown native GenomicsDB workspace interval set rule: " +
                             parsed->first);
}

bool native_workspace_has_interval_metadata(const std::string& workspace) {
    const auto metadata_path = std::filesystem::path(workspace) / "fastgatk-workspace.json";
    std::ifstream stream(metadata_path);
    if (!stream) return false;
    std::ostringstream text;
    text << stream.rdbuf();
    if (!stream.good() && !stream.eof())
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB workspace metadata: " +
                                 metadata_path.string());
    return text.str().find("\"import_intervals\":[") != std::string::npos;
}

std::uint64_t parse_failure_after_native_inputs() {
    const auto value = env("FASTGATK_GENOMICSDB_FAIL_AFTER_INPUTS");
    if (value.empty()) return 0;
    std::size_t consumed = 0;
    try {
        const auto parsed = std::stoull(value, &consumed);
        if (consumed != value.size() || parsed == 0)
            throw std::invalid_argument("non-positive");
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            "invalid positive FASTGATK_GENOMICSDB_FAIL_AFTER_INPUTS: " + value);
    }
}

std::filesystem::path native_input_target(const std::filesystem::path& directory,
                                          std::size_t index,
                                          const std::filesystem::path& source) {
    std::ostringstream basename;
    basename << "input-" << std::setw(8) << std::setfill('0') << index
             << source.extension().string();
    return directory / basename.str();
}

#if FASTGATK_HAS_HTSLIB
void materialize_native_input_rewritten(
    const std::string& source, const std::filesystem::path& target,
    const std::optional<std::string>& sample_name,
    const std::vector<std::string>& interval_selectors,
    fastgatk::io::HtsIntervalSetRule interval_set_rule,
    std::uint64_t& interval_records_kept,
    std::uint64_t& interval_records_skipped) {
    htsFile* input = hts_open(source.c_str(), "r");
    if (!input)
        throw std::runtime_error("BAD_INPUT: cannot open native GenomicsDB input for sample-map rewrite: " + source);
    bcf_hdr_t* input_header = bcf_hdr_read(input);
    if (!input_header) {
        hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read native GenomicsDB input header: " + source);
    }
    const int input_samples = bcf_hdr_nsamples(input_header);
    if (sample_name.has_value() && input_samples != 1) {
        bcf_hdr_destroy(input_header);
        hts_close(input);
        throw std::runtime_error(
            "BAD_INPUT: --sample-name-map input must contain exactly one sample: " + source);
    }
    bcf_hdr_t* output_header = bcf_hdr_dup(input_header);
    if (!output_header) {
        bcf_hdr_destroy(input_header);
        hts_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate native GenomicsDB input header: " + source);
    }
    if (sample_name.has_value()) {
        // Removing all samples resets the sample dictionary.  Clear
        // keep_samples afterwards: this header is a writer header, so
        // bcf_write must retain the one input sample payload rather than
        // apply the read-side subset.
        if (bcf_hdr_set_samples(output_header, nullptr, 0) != 0 ||
            bcf_hdr_add_sample(output_header, sample_name->c_str()) != 0) {
            bcf_hdr_destroy(output_header);
            bcf_hdr_destroy(input_header);
            hts_close(input);
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: cannot install --sample-name-map sample name: " +
                *sample_name);
        }
        free(output_header->keep_samples);
        output_header->keep_samples = nullptr;
        output_header->nsamples_ori = bcf_hdr_nsamples(output_header);
    }
    const auto extension = target.extension().string();
    const char* mode = extension == ".gz" || extension == ".bgz" ? "wz" :
                       extension == ".bcf" ? "wb" : "w";
    htsFile* output = hts_open(target.string().c_str(), mode);
    if (!output || bcf_hdr_write(output, output_header) != 0) {
        if (output) hts_close(output);
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(input_header);
        hts_close(input);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write rewritten native input: " +
                                     target.string());
    }
    bcf1_t* record = bcf_init();
    if (!record) {
        hts_close(output);
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(input_header);
        hts_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate native GenomicsDB input record");
    }
    int status = 0;
    const auto intervals = parse_native_import_intervals(input_header, interval_selectors);
    while ((status = bcf_read(input, input_header, record)) == 0) {
        if (!native_record_overlaps_import_intervals(
                input_header, record, intervals, interval_set_rule)) {
            ++interval_records_skipped;
            continue;
        }
        if (!interval_selectors.empty()) ++interval_records_kept;
        if (bcf_write(output, output_header, record) != 0) {
            bcf_destroy(record);
            hts_close(output);
            bcf_hdr_destroy(output_header);
            bcf_hdr_destroy(input_header);
            hts_close(input);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write rewritten native input record: " +
                                     target.string());
        }
    }
    const int output_status = hts_close(output);
    const int input_status = hts_close(input);
    if (status < -1 || output_status != 0 || input_status != 0) {
        bcf_destroy(record);
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(input_header);
        throw std::runtime_error("BAD_INPUT: failed reading rewritten native GenomicsDB input: " + source);
    }
    bcf_destroy(record);
    bcf_hdr_destroy(output_header);
    bcf_hdr_destroy(input_header);
}
#endif

std::vector<std::string> materialize_native_inputs(
    const std::filesystem::path& workspace,
    const std::vector<std::string>& inputs,
    NativeWorkspaceCheckpoint& checkpoint,
    const std::filesystem::path& checkpoint_path,
    bool allow_failure_injection,
    const std::map<std::string, std::string>& sample_name_by_input,
    const std::vector<std::string>& interval_selectors,
    fastgatk::io::HtsIntervalSetRule interval_set_rule,
    std::uint64_t& interval_records_kept,
    std::uint64_t& interval_records_skipped) {
    if (checkpoint.inputs != inputs || checkpoint.completed.size() != inputs.size())
        throw std::runtime_error("BAD_INPUT: native GenomicsDB checkpoint input set mismatch");
    const auto directory = workspace / "native-inputs";
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create native GenomicsDB input directory: " +
                                 directory.string());
    std::vector<std::string> materialized;
    materialized.reserve(inputs.size());
    const auto fail_after = allow_failure_injection ? parse_failure_after_native_inputs() : 0;
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        const auto source = std::filesystem::path(inputs[index]);
        const auto target = native_input_target(directory, index, source);
        const auto mapped_sample = sample_name_by_input.find(inputs[index]);
        const bool rewritten_input = !interval_selectors.empty() ||
            mapped_sample != sample_name_by_input.end();
        if (rewritten_input) {
            // An overwrite/resume may leave an index from an earlier
            // unfiltered materialization.  It is unsafe after rewriting the
            // record stream, so remove it before publishing the replacement.
            for (const char* suffix : {".tbi", ".csi"}) {
                std::error_code stale_error;
                std::filesystem::remove(target.string() + suffix, stale_error);
                if (stale_error)
                    throw std::runtime_error(
                        "OUTPUT_CONTRACT_FAILURE: cannot remove stale native GenomicsDB index: " +
                        (target.string() + suffix));
            }
        }
        if (checkpoint.completed[index]) {
            error.clear();
            if (!std::filesystem::is_regular_file(target, error) || error ||
                std::filesystem::file_size(target, error) == 0 || error)
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: checkpointed native GenomicsDB input is missing: " +
                    target.string());
        } else {
#if FASTGATK_HAS_HTSLIB
            if (rewritten_input) {
                const std::optional<std::string> sample = mapped_sample == sample_name_by_input.end()
                    ? std::nullopt : std::optional<std::string>(mapped_sample->second);
                materialize_native_input_rewritten(
                    inputs[index], target, sample, interval_selectors, interval_set_rule,
                    interval_records_kept, interval_records_skipped);
            } else {
                std::filesystem::copy_file(source, target,
                                            std::filesystem::copy_options::overwrite_existing, error);
            }
#else
            std::filesystem::copy_file(source, target,
                                        std::filesystem::copy_options::overwrite_existing, error);
#endif
            if (error)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize native GenomicsDB input: " +
                                         source.string());
            error.clear();
            if (!std::filesystem::is_regular_file(target, error) || error ||
                std::filesystem::file_size(target, error) == 0 || error)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: materialized native GenomicsDB input is empty: " +
                                         target.string());
        }
        // Indexed iteration is optional for the native reader, but preserving
        // a sidecar makes the materialized file usable by downstream HTSlib
        // tools without reopening the original path.
        // A copied index describes the original byte stream and is invalid
        // after interval filtering or a sample-name rewrite.  Keep source
        // indexes only for byte-for-byte materialization; rewritten inputs
        // intentionally fall back to sequential HTSlib traversal unless a
        // future index-builder path is enabled.
        if (!rewritten_input) for (const char* suffix : {".tbi", ".csi"}) {
            const auto source_index = source.string() + suffix;
            error.clear();
            if (!std::filesystem::is_regular_file(source_index, error) || error) {
                error.clear();
                continue;
            }
            const auto target_index = target.string() + suffix;
            if (!checkpoint.completed[index]) {
                std::filesystem::copy_file(source_index, target_index,
                                            std::filesystem::copy_options::overwrite_existing, error);
                if (error)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize native GenomicsDB index: " +
                                             source_index);
            } else if (!std::filesystem::is_regular_file(target_index, error) || error) {
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: checkpointed native GenomicsDB index is missing: " +
                    target_index);
            }
        }
        materialized.push_back(std::filesystem::absolute(target, error).string());
        if (error) {
            error.clear();
            materialized.back() = target.string();
        }
        if (!checkpoint.completed[index]) {
            checkpoint.completed[index] = true;
            write_native_workspace_checkpoint(checkpoint_path, checkpoint);
            if (fail_after != 0 && native_checkpoint_completed_count(checkpoint) >= fail_after)
                throw std::runtime_error(
                    "RESOURCE_EXHAUSTED: injected native GenomicsDB checkpoint interruption");
        }
    }
    write_input_index(workspace / "fastgatk-native-inputs.tsv", materialized);
    return materialized;
}

void write_native_workspace_metadata(Options& options,
                                     const std::vector<std::string>& inputs) {
    std::error_code error;
    const auto workspace = std::filesystem::path(options.workspace);
    std::filesystem::create_directories(workspace, error);
    if (error)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create native GenomicsDB workspace: " +
                                 options.workspace);
    for (const auto& input : inputs) {
        if (is_remote_uri(input))
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: native GenomicsDB workspace requires local VCF/GVCF inputs: " + input);
    }
    const auto checkpoint_path = native_workspace_checkpoint_path(options.workspace);
    NativeWorkspaceCheckpoint checkpoint;
    if (options.resume_native_workspace) {
        const auto expected_interval_rule =
            options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
                ? "INTERSECTION" : "UNION";
        checkpoint = read_native_workspace_checkpoint(
            checkpoint_path, inputs, options.interval_selectors, expected_interval_rule);
        if (checkpoint.sample_name_by_input != options.sample_name_by_input)
            throw std::runtime_error(
                "BAD_INPUT: native GenomicsDB checkpoint sample-name-map does not match resume inputs");
        options.checkpoint_resumed = true;
        options.checkpoint_inputs_reused = native_checkpoint_completed_count(checkpoint);
    } else {
        checkpoint = make_native_checkpoint(
            inputs, options.sample_name_by_input, options.interval_selectors,
            options.interval_set_rule);
        write_native_workspace_checkpoint(checkpoint_path, checkpoint);
        options.checkpoint_resumed = false;
        options.checkpoint_inputs_reused = 0;
    }
    write_workspace_input_index(options.workspace, inputs);
    const auto materialized = materialize_native_inputs(
        workspace, inputs, checkpoint, checkpoint_path, !options.resume_native_workspace,
        options.sample_name_by_input, options.interval_selectors,
        options.interval_set_rule, options.interval_records_kept,
        options.interval_records_skipped);
    write_record_index(workspace / "fastgatk-record-index.tsv", materialized);
    const auto metadata_path = workspace / "fastgatk-workspace.json";
    const auto temporary = metadata_path.string() + ".tmp";
    std::ofstream metadata(temporary, std::ios::trunc);
    if (!metadata)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write native GenomicsDB metadata");
    metadata << "{\"schema_version\":1,\"backend\":\"fastgatk-sparse-index\"," 
             << "\"storage_format\":\"fastgatk-portable-sparse-v1\"," 
             << "\"tiledb_schema\":false,\"genomicsdb_query_compatible\":false,"
             << "\"query_engine\":\"fastgatk-native-genotype-gvcf\"," 
             << "\"input_count\":" << inputs.size() << ",\"batch_size\":"
             << options.effective_batch << ",\"reader_threads\":"
             << options.effective_reader_threads << ",\"materialized_index\":\"fastgatk-native-inputs.tsv\",\"record_index\":\"fastgatk-record-index.tsv\",\"record_index_schema\":1," 
             << "\"incremental_update\":" << (options.update_workspace ? "true" : "false") << ","
             << "\"checkpoint_schema\":1,\"checkpoint_resumed\":"
             << (options.checkpoint_resumed ? "true" : "false")
             << ",\"checkpoint_inputs_reused\":" << options.checkpoint_inputs_reused
             << ",\"sample_name_map_entries\":" << options.sample_name_map_entries
             << ",\"sample_name_map\":[";
    bool first_sample_mapping = true;
    for (const auto& input : inputs) {
        const auto mapping = options.sample_name_by_input.find(input);
        if (mapping == options.sample_name_by_input.end()) continue;
        if (!first_sample_mapping) metadata << ',';
        first_sample_mapping = false;
        metadata << "{\"input\":\"" << json_escape(input)
                 << "\",\"sample_name\":\"" << json_escape(mapping->second)
                 << "\"}";
    }
    metadata << "],\"interval_set_rule\":\""
             << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
                     ? "INTERSECTION" : "UNION")
             << "\",\"import_intervals\":[";
    for (std::size_t index = 0; index < options.interval_selectors.size(); ++index) {
        if (index != 0) metadata << ',';
        metadata << '"' << json_escape(options.interval_selectors[index]) << '"';
    }
    metadata << "],\"inputs\":[";
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        if (index != 0) metadata << ',';
        metadata << '"' << json_escape(inputs[index]) << '"';
    }
    metadata << "],\"materialized_inputs\":[";
    for (std::size_t index = 0; index < materialized.size(); ++index) {
        if (index != 0) metadata << ',';
        metadata << '"' << json_escape(materialized[index]) << '"';
    }
    metadata << "]}\n";
    metadata.flush();
    if (!metadata)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize native GenomicsDB metadata");
    std::filesystem::rename(temporary, metadata_path, error);
    if (error) {
        std::filesystem::remove(metadata_path, error);
        error.clear();
        std::filesystem::rename(temporary, metadata_path, error);
    }
    if (error)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish native GenomicsDB metadata");
    std::filesystem::remove(checkpoint_path, error);
    if (error)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove native GenomicsDB checkpoint: " +
                                 checkpoint_path.string());
}

int run_tool(Options options, const fastgatk::runtime::ResourceSnapshot& resources) {
    validate_input_paths(options);
    WorkspaceLock workspace_lock(options.workspace);
    validate_workspace_preflight(options);
    const auto command = build_backend_argv(options, resources);
    std::cerr << "fastgatk-genomicsdb-import: "
              << ((options.native_workspace || env("FASTGATK_GENOMICSDB_BACKEND") == "native")
                      ? "native sparse workspace" : "external GenomicsDB backend")
              << " batch_size=" << options.effective_batch
              << " reader_threads=" << options.effective_reader_threads;
    if (!options.tmp_dir.empty()) std::cerr << " tmp_dir=" << options.tmp_dir;
    std::cerr << '\n';
    auto workspace_inputs = collect_workspace_inputs(options);
    if (workspace_inputs.empty())
        throw std::runtime_error("BAD_INPUT: no VCF/GVCF inputs found for workspace index");
    const auto configured_backend = env("FASTGATK_GENOMICSDB_BACKEND");
    const bool use_native_workspace = options.native_workspace || configured_backend == "native";
    if (use_native_workspace && options.update_workspace) {
        // A native sparse workspace is an input-index publication rather than
        // a TileDB array. It can nevertheless provide a safe incremental
        // boundary by atomically rebuilding the deterministic materialized
        // index from the existing inputs plus the new shards. Require the
        // prior index and preserve sample-map rewrites from its sidecar; do
        // not guess at opaque/partial workspaces.
        const auto prior_index = std::filesystem::path(options.workspace) / "fastgatk-inputs.tsv";
        std::error_code prior_error;
        if (!std::filesystem::is_regular_file(prior_index, prior_error) || prior_error)
            throw std::runtime_error(
                "BAD_INPUT: native incremental update requires fastgatk-inputs.tsv: " +
                prior_index.string());
        // GATK ignores -L supplied on an update and reuses the initial
        // workspace interval set.  Preserve that boundary for the native
        // sparse mode as well; an initial whole-genome import is represented
        // by an explicit empty selector array in the current metadata.
        if (native_workspace_has_interval_metadata(options.workspace)) {
            options.interval_selectors =
                read_native_workspace_interval_selectors(options.workspace);
            options.interval_set_rule =
                read_native_workspace_interval_set_rule(options.workspace);
        }
        const auto prior_sample_map = read_native_workspace_sample_name_map(options.workspace);
        for (const auto& mapping : prior_sample_map) {
            const auto existing = options.sample_name_by_input.find(mapping.first);
            if (existing != options.sample_name_by_input.end() && existing->second != mapping.second)
                throw std::runtime_error(
                    "BAD_INPUT: native incremental update sample-name-map conflicts with existing workspace: " +
                    mapping.first);
            options.sample_name_by_input.emplace(mapping.first, mapping.second);
        }
        options.sample_name_map_entries = options.sample_name_by_input.size();
        workspace_inputs = merge_workspace_input_index(options.workspace, workspace_inputs);
        if (workspace_inputs.empty())
            throw std::runtime_error("BAD_INPUT: native incremental update has no workspace inputs");
    }
    if (use_native_workspace) {
        write_native_workspace_metadata(options, workspace_inputs);
    } else {
        const auto status = run_backend(command);
        if (status != 0)
            throw std::runtime_error("BACKEND_FAILURE: GenomicsDB backend exited with status " + std::to_string(status));
        if (!workspace_complete(options.workspace))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GenomicsDB workspace is missing or empty");
        const auto indexed_inputs = options.update_workspace
            ? merge_workspace_input_index(options.workspace, workspace_inputs)
            : workspace_inputs;
        write_workspace_input_index(options.workspace, indexed_inputs);
    }
    const auto manifest_path = options.manifest.empty()
        ? options.workspace + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write adapter manifest: " + manifest_path);
    const auto complete = workspace_complete(options.workspace);
    manifest << "{\"schema_version\":1,\"tool\":\"GenomicsDBImport\","
             << "\"implementation\":\"fastgatk-genomicsdb-import\","
             << "\"status\":\"adapter\",\"execution_mode\":\""
             << (use_native_workspace ? "native-sparse-index" : "external-genomicsdb") << "\","
             << "\"primary_output\":\"" << json_escape(options.workspace) << "\","
             << "\"primary_output_kind\":\"genomicsdb-workspace\","
             << "\"compatibility\":{\"external_backend\":"
             << (use_native_workspace ? "false" : "true")
             << ",\"native_workspace_index\":" << (use_native_workspace ? "true" : "false")
             << ",\"native_workspace_materialized\":" << (use_native_workspace ? "true" : "false")
             << ",\"native_workspace_checkpoint\":" << (use_native_workspace ? "true" : "false")
             << ",\"external_update_workspace\":" << (options.update_workspace ? "true" : "false")
             << ",\"native_incremental_update\":"
             << (use_native_workspace && options.update_workspace ? "true" : "false")
             << ",\"checkpoint_resumed\":" << (options.checkpoint_resumed ? "true" : "false")
             << ",\"resource_adaptation\":true,"
             << "\"workspace_nonempty\":" << (complete ? "true" : "false")
             << ",\"workspace_input_index\":true,\"workspace_lock\":true,\"bit_identical_to_gatk\":"
             << (use_native_workspace ? "false" : "true")
             << ",\"tiledb_storage\":" << (use_native_workspace ? "false" : "true")
             << ",\"genomicsdb_query_compatible\":"
             << (use_native_workspace ? "false" : "true")
             << ",\"storage_format\":\""
             << (use_native_workspace ? "fastgatk-portable-sparse-v1" : "external-genomicsdb")
             << "\""
             << ",\"sample_name_map\":" << (options.sample_name_map_entries != 0 ? "true" : "false")
             << "},\"outputs\":[{\"path\":\""
             << json_escape(options.workspace) << "\",\"kind\":\"genomicsdb-workspace\",\"complete\":"
             << (complete ? "true" : "false") << "}],\"telemetry\":{\"resources\":"
             << resources.to_json() << ",\"requested_batch_size\":"
             << (options.has_batch ? std::to_string(options.requested_batch) : "0")
             << ",\"effective_batch_size\":" << options.effective_batch
             << ",\"requested_reader_threads\":" << options.requested_reader_threads
             << ",\"effective_reader_threads\":" << options.effective_reader_threads
             << ",\"tmp_dir\":\"" << json_escape(options.tmp_dir)
             << "\",\"workspace_input_count\":"
             << ((use_native_workspace || !options.update_workspace)
                     ? workspace_inputs.size()
                     : merge_workspace_input_index(options.workspace, workspace_inputs).size())
             << ",\"sample_name_map_entries\":" << options.sample_name_map_entries
             << ",\"checkpoint_inputs_reused\":" << options.checkpoint_inputs_reused
             << ",\"checkpoint_inputs_completed\":" << (use_native_workspace ? workspace_inputs.size() : 0)
             << ",\"interval_selectors\":" << options.interval_selectors.size()
             << ",\"interval_set_rule\":\""
             << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
                     ? "INTERSECTION" : "UNION")
             << "\",\"interval_records_kept\":" << options.interval_records_kept
             << ",\"interval_records_skipped\":" << options.interval_records_skipped
             << ",\"checkpoint_path\":\""
             << json_escape(native_workspace_checkpoint_path(options.workspace).string()) << "\""
             << ",\"workspace_lock_path\":\"" << json_escape(workspace_lock.path().string())
             << "\"}}\n";
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize adapter manifest");
    std::cout << "{\"tool\":\"GenomicsDBImport\",\"status\":\"adapter\","
              << "\"workspace\":\"" << json_escape(options.workspace)
              << "\",\"effective_batch_size\":" << options.effective_batch
              << ",\"effective_reader_threads\":" << options.effective_reader_threads
              << ",\"sample_name_map_entries\":" << options.sample_name_map_entries
              << ",\"interval_selectors\":" << options.interval_selectors.size()
              << ",\"interval_records_kept\":" << options.interval_records_kept
              << ",\"interval_records_skipped\":" << options.interval_records_skipped
              << ",\"checkpoint_resumed\":" << (options.checkpoint_resumed ? "true" : "false")
              << ",\"external_update_workspace\":" << (options.update_workspace ? "true" : "false")
              << ",\"native_incremental_update\":"
              << (use_native_workspace && options.update_workspace ? "true" : "false")
              << ",\"execution_mode\":\""
              << (use_native_workspace ? "native-sparse-index" : "external-genomicsdb") << "\"}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-genomicsdb-import: " << error.what() << '\n';
        return 2;
    }
}
