#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

#include <algorithm>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::string reference;
    std::string output;
    std::string manifest;
    std::vector<std::string> intervals;
    std::vector<std::string> exclude_intervals;
    std::string mode = "INTERVAL_SUBDIVISION";
    std::string extension = "-scattered.interval_list";
    std::string prefix;
    int scatter_count = 1;
    int digits = 4;
    std::int64_t min_contig_size = 0;
    bool dont_mix_contigs = false;
};

struct Contig {
    std::string name;
    std::int64_t length = 0;
};

struct Interval {
    int tid = -1;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

struct KernelStats {
    std::uint64_t bases = 0;
    std::size_t records = 0;
    std::size_t prepare_bytes = 0;
    std::size_t device_bytes = 0;
    std::size_t batches = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space;
};

KernelStats summarize_intervals_kokkos(const std::vector<Interval>& intervals) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host("split-intervals-v2");
    host.records = intervals.size();
    host.bytes = intervals.size() * sizeof(Interval);
    fastgatk::core::KernelPlan<ExecSpace> plan("split-intervals-bases");
    plan.begin_prepare(host);

    Kokkos::View<std::int64_t*> lengths("interval_lengths", intervals.size());
    auto host_lengths = Kokkos::create_mirror_view(lengths);
    for (std::size_t index = 0; index < intervals.size(); ++index)
        host_lengths(index) = intervals[index].end - intervals[index].start;
    Kokkos::deep_copy(lengths, host_lengths);
    fastgatk::core::DeviceBatch<ExecSpace> device(intervals.size());
    device.bind("interval_lengths", lengths);
    ExecSpace().fence();
    plan.end_prepare(device);

    plan.begin_execute();
    std::uint64_t bases = 0;
    Kokkos::parallel_reduce("split_intervals_base_count",
                            Kokkos::RangePolicy<ExecSpace>(0, intervals.size()),
                            KOKKOS_LAMBDA(const std::size_t index, std::uint64_t& value) {
                                value += static_cast<std::uint64_t>(lengths(index));
                            },
                            bases);
    ExecSpace().fence();
    plan.end_execute();

    const auto& telemetry = plan.telemetry();
    KernelStats stats;
    stats.bases = bases;
    stats.records = intervals.size();
    stats.prepare_bytes = telemetry.prepare_bytes;
    stats.device_bytes = telemetry.device_bytes;
    stats.batches = 1;
    stats.prepare_seconds = telemetry.prepare_seconds;
    stats.execute_seconds = telemetry.execute_seconds;
    stats.execution_space = ExecSpace::name();
    return stats;
}

std::string option_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, std::string_view name,
                          std::string_view short_name = {}) {
    if (const auto inline_value = option_value(argument, name); !inline_value.empty())
        return inline_value;
    if ((argument == name || (!short_name.empty() && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument("missing value for " + std::string(name));
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool suffix_ci(const std::string& value, std::string_view ending) {
    if (value.size() < ending.size()) return false;
    const auto offset = value.size() - ending.size();
    for (std::size_t index = 0; index < ending.size(); ++index) {
        const auto lhs = static_cast<unsigned char>(value[offset + index]);
        const auto rhs = static_cast<unsigned char>(ending[index]);
        if (std::tolower(lhs) != std::tolower(rhs)) return false;
    }
    return true;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-split-intervals (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE            reference FASTA with .fai\n"
                         "  -L, --intervals SELECTOR        repeatable interval or file\n"
                         "  -XL, --exclude-intervals SELECTOR\n"
                         "  -O, --output DIRECTORY           output directory\n"
                         "      --scatter-count N             number of shards (default 1)\n"
                         "      --subdivision-mode MODE      INTERVAL_SUBDIVISION, INTERVAL_COUNT,\n"
                         "                                     BALANCING_WITHOUT_INTERVAL_SUBDIVISION\n"
                         "      --min-contig-size N          filter reference-derived contigs shorter than N\n"
                         "      --dont-mix-contigs           split each shard by contig\n"
                         "      --extension SUFFIX           output suffix\n"
                         "      --interval-file-prefix PFX   output filename prefix\n"
                         "      --interval-file-num-digits N zero-padding width\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-R" || argument == "--reference" ||
            !option_value(argument, "--reference").empty())
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || argument == "--intervals" ||
                 !option_value(argument, "--intervals").empty())
            options.intervals.push_back(require_value(index, argc, argv, argument, "--intervals", "-L"));
        else if (argument == "-XL" || argument == "--exclude-intervals" ||
                 !option_value(argument, "--exclude-intervals").empty())
            options.exclude_intervals.push_back(require_value(index, argc, argv, argument,
                                                               "--exclude-intervals", "-XL"));
        else if (argument == "-O" || argument == "--output" ||
                 !option_value(argument, "--output").empty())
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "--scatter-count" || argument == "--scatter" ||
                 !option_value(argument, "--scatter-count").empty() ||
                 !option_value(argument, "--scatter").empty()) {
            const auto name = argument.rfind("--scatter=", 0) == 0 ? "--scatter" : "--scatter-count";
            options.scatter_count = std::stoi(require_value(index, argc, argv, argument, name));
        } else if (argument == "--subdivision-mode" || argument == "--mode" ||
                   !option_value(argument, "--subdivision-mode").empty() ||
                   !option_value(argument, "--mode").empty()) {
            const auto name = argument.rfind("--mode", 0) == 0 ? "--mode" : "--subdivision-mode";
            options.mode = require_value(index, argc, argv, argument, name);
            std::transform(options.mode.begin(), options.mode.end(), options.mode.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
        } else if (argument == "--dont-mix-contigs") {
            options.dont_mix_contigs = true;
        } else if (argument == "--extension" || !option_value(argument, "--extension").empty()) {
            options.extension = require_value(index, argc, argv, argument, "--extension");
        } else if (argument == "--interval-file-prefix" ||
                   !option_value(argument, "--interval-file-prefix").empty()) {
            options.prefix = require_value(index, argc, argv, argument, "--interval-file-prefix");
        } else if (argument == "--interval-file-num-digits" ||
                   !option_value(argument, "--interval-file-num-digits").empty()) {
            options.digits = std::stoi(require_value(index, argc, argv, argument,
                                                      "--interval-file-num-digits"));
        } else if (argument == "--min-contig-size" ||
                   !option_value(argument, "--min-contig-size").empty()) {
            options.min_contig_size = std::stoll(require_value(index, argc, argv, argument,
                                                                "--min-contig-size"));
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !option_value(argument, "--output-manifest").empty() ||
                   !option_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = require_value(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--lenient") {
            // Accepted launcher-compatible flags; sequence dictionaries are
            // generated from the requested reference and remain deterministic.
        } else if (argument == "--java-options" || argument == "--verbosity") {
            (void)require_value(index, argc, argv, argument, argument);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.scatter_count < 1) throw std::invalid_argument("--scatter-count must be positive");
    if (options.min_contig_size < 0)
        throw std::invalid_argument("--min-contig-size must be non-negative");
    if (options.digits < 1 || options.digits > 12)
        throw std::invalid_argument("--interval-file-num-digits must be in [1,12]");
    if (options.mode != "INTERVAL_SUBDIVISION" && options.mode != "INTERVAL_COUNT" &&
        options.mode != "BALANCING_WITHOUT_INTERVAL_SUBDIVISION" &&
        options.mode != "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW" &&
        options.mode != "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER")
        throw std::invalid_argument("unsupported --subdivision-mode: " + options.mode);
    return options;
}

std::vector<Contig> read_fai(const std::string& reference) {
    std::ifstream stream(reference + ".fai");
    if (!stream) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + reference);
    std::vector<Contig> contigs;
    std::string line;
    while (std::getline(stream, line)) {
        line = trim(std::move(line));
        if (line.empty()) continue;
        std::istringstream fields(line);
        Contig contig;
        std::string ignored;
        fields >> contig.name >> contig.length >> ignored;
        if (contig.name.empty() || contig.length <= 0)
            throw std::runtime_error("BAD_INPUT: malformed FASTA index line");
        contigs.push_back(std::move(contig));
    }
    if (contigs.empty()) throw std::runtime_error("BAD_INPUT: FASTA index has no contigs");
    return contigs;
}

int find_contig(const std::vector<Contig>& contigs, const std::string& name) {
    for (std::size_t index = 0; index < contigs.size(); ++index)
        if (contigs[index].name == name) return static_cast<int>(index);
    return -1;
}

std::int64_t coordinate(const std::string& value, const std::string& source) {
    if (value.empty()) throw std::runtime_error("BAD_INPUT: empty interval coordinate: " + source);
    std::size_t consumed = 0;
    const auto parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed < 1)
        throw std::runtime_error("BAD_INPUT: invalid interval coordinate: " + source);
    return parsed;
}

Interval parse_literal(const std::vector<Contig>& contigs, const std::string& raw) {
    const auto source = trim(raw);
    const auto colon = source.find(':');
    const auto name = source.substr(0, colon);
    const auto tid = find_contig(contigs, name);
    if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in reference: " + name);
    if (colon == std::string::npos) return Interval{tid, 0, contigs[tid].length};
    const auto coordinates = source.substr(colon + 1);
    const auto dash = coordinates.find('-');
    const auto start = coordinate(coordinates.substr(0, dash), source);
    const auto end = dash == std::string::npos ? start : coordinate(coordinates.substr(dash + 1), source);
    if (end < start || end > contigs[tid].length)
        throw std::runtime_error("BAD_INPUT: interval outside reference: " + source);
    return Interval{tid, start - 1, end};
}

void append_interval_file(const std::vector<Contig>& contigs, const std::string& path,
                          std::vector<Interval>& result) {
    const bool bed = suffix_ci(path, ".bed") || suffix_ci(path, ".bed.gz");
    std::size_t parsed = 0;
    const auto process_line = [&](std::string line) {
        line = trim(std::move(line));
        if (line.empty() || line.front() == '#' || line.front() == '@' ||
            line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) return;
        std::istringstream fields(line);
        std::string first;
        fields >> first;
        if (first.find(':') != std::string::npos) {
            result.push_back(parse_literal(contigs, first));
            ++parsed;
            return;
        }
        std::string second, third;
        if (!(fields >> second >> third))
            throw std::runtime_error("BAD_INPUT: malformed interval file line: " + line);
        const auto tid = find_contig(contigs, first);
        if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in reference: " + first);
        std::size_t consumed_start = 0, consumed_end = 0;
        const auto start = std::stoll(second, &consumed_start);
        const auto end = std::stoll(third, &consumed_end);
        if (consumed_start != second.size() || consumed_end != third.size())
            throw std::runtime_error("BAD_INPUT: malformed interval coordinates: " + line);
        const auto begin = bed ? start : start - 1;
        const auto finish = end;
        if (begin < 0 || finish <= begin || finish > contigs[tid].length)
            throw std::runtime_error("BAD_INPUT: interval outside reference: " + line);
        result.push_back(Interval{tid, begin, finish});
        ++parsed;
    };
    if (suffix_ci(path, ".gz")) {
#if FASTGATK_HAS_HTSLIB
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file)
            throw std::runtime_error("BAD_INPUT: cannot open compressed interval file: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (fastgatk::io::read_text_line(file, &line, path) >= 0)
                process_line(line.s == nullptr ? std::string{} : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(file);
            throw;
        }
        const auto close_status = hts_close(file);
        free(line.s);
        if (close_status != 0)
            throw std::runtime_error("BAD_INPUT: failed reading compressed interval file: " + path);
#else
        throw std::runtime_error(
            "BACKEND_UNAVAILABLE: compressed interval files require HTSlib: " + path);
#endif
    } else {
        std::ifstream stream(path);
        if (!stream) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        for (std::string line; std::getline(stream, line);)
            process_line(std::move(line));
        if (stream.bad()) throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
    }
    if (parsed == 0) throw std::runtime_error("BAD_INPUT: interval file contains no intervals: " + path);
}

std::vector<Interval> materialize(const std::vector<Contig>& contigs,
                                  const std::vector<std::string>& selectors) {
    std::vector<Interval> intervals;
    for (const auto& selector : selectors) {
        std::error_code error;
        if (std::filesystem::is_regular_file(selector, error) && !error)
            append_interval_file(contigs, selector, intervals);
        else
            intervals.push_back(parse_literal(contigs, selector));
    }
    if (intervals.empty()) {
        for (std::size_t tid = 0; tid < contigs.size(); ++tid)
            intervals.push_back(Interval{static_cast<int>(tid), 0, contigs[tid].length});
    }
    std::sort(intervals.begin(), intervals.end(), [](const Interval& lhs, const Interval& rhs) {
        if (lhs.tid != rhs.tid) return lhs.tid < rhs.tid;
        if (lhs.start != rhs.start) return lhs.start < rhs.start;
        return lhs.end < rhs.end;
    });
    std::vector<Interval> merged;
    for (const auto& interval : intervals) {
        if (!merged.empty() && merged.back().tid == interval.tid && interval.start <= merged.back().end)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    return merged;
}

std::vector<Interval> subtract(const std::vector<Interval>& input,
                               const std::vector<Interval>& excludes) {
    std::vector<Interval> result;
    for (const auto& source : input) {
        std::vector<Interval> pieces{source};
        for (const auto& exclusion : excludes) {
            if (exclusion.tid != source.tid) continue;
            std::vector<Interval> next;
            for (const auto& piece : pieces) {
                if (exclusion.end <= piece.start || exclusion.start >= piece.end) {
                    next.push_back(piece);
                    continue;
                }
                if (piece.start < exclusion.start)
                    next.push_back(Interval{piece.tid, piece.start, exclusion.start});
                if (exclusion.end < piece.end)
                    next.push_back(Interval{piece.tid, exclusion.end, piece.end});
            }
            pieces.swap(next);
            if (pieces.empty()) break;
        }
        result.insert(result.end(), pieces.begin(), pieces.end());
    }
    return result;
}

std::uint64_t scatter_weight(const Interval& interval, const std::string& mode) {
    return mode == "INTERVAL_COUNT" || mode == "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER"
        ? 1ULL : static_cast<std::uint64_t>(interval.end - interval.start);
}

std::uint64_t scatter_list_weight(const std::vector<Interval>& intervals, const std::string& mode) {
    std::uint64_t weight = 0;
    for (const auto& interval : intervals) weight += scatter_weight(interval, mode);
    return weight;
}

Interval interval_prefix(const Interval& interval, std::int64_t length) {
    return Interval{interval.tid, interval.start, interval.start + length};
}

Interval interval_suffix(const Interval& interval, std::int64_t length) {
    return Interval{interval.tid, interval.start + length, interval.end};
}

// Mirrors Picard IntervalListScatter's state machine: each output computes a
// dynamic target from the remaining weight and remaining lists, while the
// final requested list drains all intervals left in the queue.  Keeping this
// in one Host routine makes the boundary semantics identical across modes;
// the only device-side work is the independent Kokkos base-count telemetry.
std::vector<std::vector<Interval>> scatter_picard(const std::vector<Interval>& intervals,
                                                  int scatter_count,
                                                  const std::string& mode) {
    std::deque<Interval> queue(intervals.begin(), intervals.end());
    std::vector<std::vector<Interval>> outputs;
    if (queue.empty()) return outputs;

    std::uint64_t weight_remaining = scatter_list_weight(intervals, mode);
    std::uint64_t ideal = std::max<std::uint64_t>(1, weight_remaining / static_cast<std::uint64_t>(scatter_count));
    if (mode == "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW") {
        for (const auto& interval : intervals)
            ideal = std::max<std::uint64_t>(ideal, scatter_weight(interval, mode));
    }

    for (int returned = 1; !queue.empty(); ++returned) {
        std::vector<Interval> output;
        if (returned >= scatter_count) {
            while (!queue.empty()) {
                output.push_back(queue.front());
                queue.pop_front();
            }
            outputs.push_back(std::move(output));
            break;
        }

        while (!queue.empty()) {
            const auto interval = queue.front();
            queue.pop_front();
            const auto current_weight = scatter_list_weight(output, mode);
            const auto interval_weight = scatter_weight(interval, mode);
            const auto remaining_lists = static_cast<std::uint64_t>(scatter_count - returned);
            const auto remaining_after_current = weight_remaining >= current_weight
                ? weight_remaining - current_weight : 0ULL;
            const double target_remaining = remaining_lists == 0
                ? 0.0 : static_cast<double>(remaining_after_current) /
                        static_cast<double>(remaining_lists);

            bool include = false;
            bool split = false;
            std::int64_t split_length = 0;
            if (mode == "INTERVAL_SUBDIVISION") {
                const auto remaining_target = ideal > current_weight ? ideal - current_weight : 0ULL;
                if (remaining_target >= interval_weight) {
                    include = true;
                } else if (remaining_target > 0 && remaining_target < interval_weight) {
                    include = true;
                    split = true;
                    split_length = static_cast<std::int64_t>(remaining_target);
                }
            } else if (mode == "INTERVAL_COUNT") {
                include = current_weight + interval_weight <= ideal;
            } else if (mode == "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER") {
                // Picard's distributed-remainder scatterer includes the next
                // interval only when the dynamic target is strictly larger
                // than the current list weight.
                include = target_remaining > static_cast<double>(current_weight);
            } else {
                const auto new_weight = current_weight + interval_weight;
                include = new_weight <= ideal ||
                          (mode == "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW" &&
                           target_remaining > static_cast<double>(ideal));
            }

            // Picard can produce an empty list for a pathological oversized
            // interval in the non-overflow mode.  Preserve data and avoid an
            // empty output artifact by assigning that interval to the list.
            if (!include && output.empty()) include = true;
            if (include) {
                if (split) {
                    output.push_back(interval_prefix(interval, split_length));
                    queue.push_front(interval_suffix(interval, split_length));
                } else {
                    output.push_back(interval);
                }
                continue;
            }
            queue.push_front(interval);
            break;
        }

        weight_remaining = weight_remaining >= scatter_list_weight(output, mode)
            ? weight_remaining - scatter_list_weight(output, mode) : 0ULL;
        if (!output.empty()) outputs.push_back(std::move(output));
    }
    return outputs;
}

std::vector<std::vector<Interval>> separate_contigs(const std::vector<std::vector<Interval>>& shards) {
    std::vector<std::vector<Interval>> output;
    for (const auto& shard : shards) {
        std::vector<std::vector<Interval>> by_contig;
        for (const auto& interval : shard) {
            auto it = std::find_if(by_contig.begin(), by_contig.end(), [&](const auto& group) {
                return !group.empty() && group.front().tid == interval.tid;
            });
            if (it == by_contig.end()) by_contig.push_back({interval});
            else it->push_back(interval);
        }
        output.insert(output.end(), by_contig.begin(), by_contig.end());
    }
    return output;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const auto character : value) {
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
    const auto contigs = read_fai(options.reference);
    auto intervals = materialize(contigs, options.intervals);
    // GATK applies min-contig-size only when the interval argument is omitted;
    // explicit -L selectors must remain authoritative even for short contigs.
    if (options.intervals.empty() && options.min_contig_size > 0) {
        intervals.erase(std::remove_if(intervals.begin(), intervals.end(), [&](const Interval& interval) {
            return contigs[static_cast<std::size_t>(interval.tid)].length < options.min_contig_size;
        }), intervals.end());
    }
    const auto excludes = materialize(contigs, options.exclude_intervals);
    if (!options.exclude_intervals.empty()) intervals = subtract(intervals, excludes);

    auto shards = scatter_picard(intervals, options.scatter_count, options.mode);
    if (options.dont_mix_contigs) shards = separate_contigs(shards);
    const auto kernel = summarize_intervals_kokkos(intervals);

    std::error_code error;
    std::filesystem::create_directories(options.output, error);
    if (error) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create output directory: " + options.output);
    std::vector<std::string> outputs;
    auto header_path = options.reference + ".dict";
    if (!std::filesystem::is_regular_file(header_path)) {
        auto dictionary = std::filesystem::path(options.reference);
        dictionary.replace_extension(".dict");
        header_path = dictionary.string();
    }
    std::vector<std::string> header;
    std::ifstream dict(header_path);
    header.emplace_back("@HD\tVN:1.6");
    std::string line;
    while (dict && std::getline(dict, line)) {
        line = trim(std::move(line));
        if (line.rfind("@SQ", 0) == 0) header.push_back(line);
    }
    if (header.size() == 1) {
        for (const auto& contig : contigs)
            header.push_back("@SQ\tSN:" + contig.name + "\tLN:" + std::to_string(contig.length));
    }
    for (int index = 0; index < static_cast<int>(shards.size()); ++index) {
        std::ostringstream name;
        name << options.prefix << std::setfill('0') << std::setw(options.digits) << index << options.extension;
        const auto path = (std::filesystem::path(options.output) / name.str()).string();
        std::ofstream out(path);
        if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write interval file: " + path);
        for (const auto& hline : header) out << hline << '\n';
        for (const auto& interval : shards[static_cast<std::size_t>(index)])
            out << contigs[static_cast<std::size_t>(interval.tid)].name << '\t'
                << interval.start + 1 << '\t' << interval.end << "\t+\t.\n";
        if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize interval file: " + path);
        outputs.push_back(path);
    }
    if (!options.manifest.empty()) {
        std::ofstream out(options.manifest);
        if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
        out << "{\"schema_version\":1,\"tool\":\"SplitIntervals\",\"implementation\":\"fastgatk-split-intervals\","
               "\"status\":\"prototype\",\"reference\":\"" << json_escape(options.reference)
            << "\",\"mode\":\"" << json_escape(options.mode) << "\",\"scatter_count\":"
            << options.scatter_count << ",\"output_shards\":" << outputs.size() << ",\"intervals\":"
            << intervals.size() << ",\"bases\":" << kernel.bases << ",\"min_contig_size\":"
            << options.min_contig_size << ",\"dont_mix_contigs\":"
            << (options.dont_mix_contigs ? "true" : "false")
            << ",\"kernel\":{\"execution_space\":\"" << json_escape(kernel.execution_space)
            << "\",\"policy\":\"RangePolicy\",\"batches\":" << kernel.batches
            << ",\"records\":" << kernel.records << ",\"prepare_bytes\":" << kernel.prepare_bytes
            << ",\"device_bytes\":" << kernel.device_bytes << ",\"prepare_seconds\":"
            << kernel.prepare_seconds << ",\"execute_seconds\":" << kernel.execute_seconds
            << "},\"outputs\":[";
        for (std::size_t index = 0; index < outputs.size(); ++index) {
            if (index != 0) out << ',';
            out << "{\"path\":\"" << json_escape(outputs[index]) << "\",\"complete\":true}";
        }
        out << "]}\n";
    }
    std::cout << "{\"tool\":\"SplitIntervals\",\"status\":\"prototype\",\"mode\":\""
              << json_escape(options.mode) << "\",\"output_shards\":" << outputs.size()
              << ",\"intervals\":" << intervals.size() << ",\"bases\":" << kernel.bases
              << ",\"min_contig_size\":" << options.min_contig_size
              << ",\"execution_space\":\"" << json_escape(kernel.execution_space)
              << "\",\"policy\":\"RangePolicy\",\"prepare_seconds\":"
              << kernel.prepare_seconds << ",\"execute_seconds\":" << kernel.execute_seconds << "}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Kokkos::ScopeGuard guard(argc, argv);
        return main_impl(argc, argv);
    }
    catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
