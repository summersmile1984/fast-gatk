#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr double DEFAULT_LOWER = 0.9;
constexpr double DEFAULT_UPPER = 1.1;
constexpr double DEFAULT_OUTLIER_Z = 2.0;
constexpr double DEFAULT_CALLING_Z = 2.0;

struct Options {
    std::string input;
    std::string output;
    std::string legacy_output;
    std::string manifest;
    double lower = DEFAULT_LOWER;
    double upper = DEFAULT_UPPER;
    double outlier_z = DEFAULT_OUTLIER_Z;
    double calling_z = DEFAULT_CALLING_Z;
    int threads = 1;
};

struct Segment {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::int64_t num_points = 0;
    double mean_log2 = 0.0;
    double copy_ratio = 0.0;
};

struct Table {
    std::vector<std::string> header_lines;
    std::vector<Segment> segments;
    std::string sample = "UNKNOWN";
};

struct Statistics {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double standard_deviation = std::numeric_limits<double>::quiet_NaN();
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string default_legacy_output(const std::string& output) {
    const std::filesystem::path path(output);
    const auto parent = path.parent_path();
    const auto stem = path.stem().string();
    const auto legacy_name = stem + ".igv.seg";
    return (parent / legacy_name).string();
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

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-call-copy-ratio-segments (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                         segment TSV\n"
                         "  -O, --output FILE                        called segment TSV\n"
                         "      --legacy-output FILE                override derived IGV .seg output\n"
                         "      --neutral-segment-copy-ratio-lower-bound F (default 0.9)\n"
                         "      --neutral-segment-copy-ratio-upper-bound F (default 1.1)\n"
                         "      --outlier-neutral-segment-copy-ratio-z-score-threshold F (default 2.0)\n"
                         "      --calling-copy-ratio-z-score-threshold F             (default 2.0)\n"
                         "      --output-manifest FILE               OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--legacy-output") || has_option(argument, "--output-legacy")) {
            options.legacy_output = require_value(
                index, argc, argv, argument,
                argument.rfind("--output-legacy", 0) == 0 ? "--output-legacy" : "--legacy-output");
        } else if (has_option(argument, "--neutral-segment-copy-ratio-lower-bound")) {
            options.lower = std::stod(require_value(index, argc, argv, argument,
                "--neutral-segment-copy-ratio-lower-bound"));
        } else if (has_option(argument, "--neutral-segment-copy-ratio-upper-bound")) {
            options.upper = std::stod(require_value(index, argc, argv, argument,
                "--neutral-segment-copy-ratio-upper-bound"));
        } else if (has_option(argument, "--outlier-z-score") ||
                   has_option(argument, "--outlier-neutral-segment-copy-ratio-z-score-threshold")) {
            const char* name = argument.rfind("--outlier-neutral-segment-copy-ratio-z-score-threshold", 0) == 0
                ? "--outlier-neutral-segment-copy-ratio-z-score-threshold" : "--outlier-z-score";
            options.outlier_z = std::stod(require_value(index, argc, argv, argument, name));
        } else if (has_option(argument, "--calling-z-score") ||
                   has_option(argument, "--calling-copy-ratio-z-score-threshold")) {
            const char* name = argument.rfind("--calling-copy-ratio-z-score-threshold", 0) == 0
                ? "--calling-copy-ratio-z-score-threshold" : "--calling-z-score";
            options.calling_z = std::stod(require_value(index, argc, argv, argument, name));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--disable-tool-default-read-filters") {
            // Compatibility switches with no effect on an already-materialized segment table.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (!std::isfinite(options.lower) || options.lower < 0.0)
        throw std::invalid_argument("java.lang.IllegalArgumentException: Copy-neutral lower bound must be non-negative.");
    if (!std::isfinite(options.upper) || options.lower >= options.upper)
        throw std::invalid_argument("java.lang.IllegalArgumentException: Copy-neutral lower bound must be less than upper bound.");
    if (!std::isfinite(options.outlier_z) || options.outlier_z <= 0.0)
        throw std::invalid_argument("java.lang.IllegalArgumentException: Outlier z-score threshold must be positive.");
    if (!std::isfinite(options.calling_z) || options.calling_z <= 0.0)
        throw std::invalid_argument("java.lang.IllegalArgumentException: Calling z-score threshold must be positive.");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    return options;
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

Table read_table(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open segment table: " + path);
    Table table;
    bool columns_seen = false;
    int contig_column = -1;
    int start_column = -1;
    int end_column = -1;
    int points_column = -1;
    int mean_column = -1;
    std::string previous_contig;
    std::int64_t previous_end = 0;
    std::unordered_map<std::string, std::size_t> sequence_rank;
    std::size_t line_number = 0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty()) continue;
        if (line.rfind("@RG", 0) == 0) {
            for (const auto& field : split_tab(line))
                if (field.rfind("SM:", 0) == 0) table.sample = field.substr(3);
            table.header_lines.push_back(line);
            continue;
        }
        if (line.rfind("@SQ", 0) == 0) {
            for (const auto& field : split_tab(line)) {
                if (field.rfind("SN:", 0) == 0)
                    sequence_rank.emplace(field.substr(3), sequence_rank.size());
            }
            table.header_lines.push_back(line);
            continue;
        }
        if (line[0] == '@') {
            table.header_lines.push_back(line);
            continue;
        }
        const auto fields = split_tab(line);
        if (!columns_seen) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "CONTIG") contig_column = static_cast<int>(i);
                else if (fields[i] == "START") start_column = static_cast<int>(i);
                else if (fields[i] == "END") end_column = static_cast<int>(i);
                else if (fields[i] == "NUM_POINTS_COPY_RATIO") points_column = static_cast<int>(i);
                else if (fields[i] == "MEAN_LOG2_COPY_RATIO") mean_column = static_cast<int>(i);
            }
            if (contig_column < 0 || start_column < 0 || end_column < 0 ||
                points_column < 0 || mean_column < 0)
                throw std::runtime_error("BAD_INPUT: segment table requires CONTIG/START/END/NUM_POINTS_COPY_RATIO/MEAN_LOG2_COPY_RATIO at " +
                                         path + ":" + std::to_string(line_number));
            columns_seen = true;
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max({contig_column, start_column, end_column,
                                                                  points_column, mean_column})))
            throw std::runtime_error("BAD_INPUT: invalid segment row at " + path + ":" + std::to_string(line_number));
        Segment segment;
        segment.contig = fields[contig_column];
        segment.start = std::stoll(fields[start_column]);
        segment.end = std::stoll(fields[end_column]);
        segment.num_points = std::stoll(fields[points_column]);
        segment.mean_log2 = std::stod(fields[mean_column]);
        // CopyRatioSegment permits non-finite MEAN_LOG2_COPY_RATIO values
        // from the TSV decoder.  SimpleCopyRatioCaller deliberately lets
        // Math.pow(2, NaN/+-Infinity) flow through its IEEE-754 statistics
        // and emits a neutral call.  Rejecting these rows here would make
        // the native tool stricter than GATK and would prevent a faithful
        // pass-through of an otherwise parseable segment table.
        if (segment.contig.empty() || segment.start < 1 || segment.end < segment.start ||
            segment.num_points < 0)
            throw std::runtime_error("BAD_INPUT: invalid segment values at " + path + ":" + std::to_string(line_number));
        // CopyRatioSegmentCollection validates the interval partition before
        // SimpleCopyRatioCaller runs: rows must be in sequence-dictionary
        // order and may not overlap.  Accepting an overlapping/unsorted table
        // changes the length-weighted statistics while Java rejects it.
        if (!sequence_rank.empty() && !sequence_rank.contains(segment.contig))
            throw std::runtime_error("BAD_INPUT: contig is absent from sequence dictionary at " + path + ":" + std::to_string(line_number));
        if (!previous_contig.empty()) {
            const auto previous_rank = sequence_rank.empty() ? 0 : sequence_rank.at(previous_contig);
            const auto current_rank = sequence_rank.empty() ? 0 : sequence_rank.at(segment.contig);
            if (current_rank < previous_rank ||
                (segment.contig == previous_contig && segment.start <= previous_end))
                throw std::runtime_error("BAD_INPUT: copy-ratio segments must be sorted and non-overlapping at " + path + ":" + std::to_string(line_number));
        }
        previous_contig = segment.contig;
        previous_end = segment.end;
        table.segments.push_back(std::move(segment));
    }
    if (!columns_seen || table.segments.empty())
        throw std::runtime_error("BAD_INPUT: segment table has no records: " + path);
    return table;
}

double kahan_sum(const std::vector<double>& values) {
    double sum = 0.0;
    double compensation = 0.0;
    double simple_sum = 0.0;
    for (const auto value : values) {
        const double corrected = value - compensation;
        const double next = sum + corrected;
        compensation = (next - sum) - corrected;
        sum = next;
        simple_sum += value;
    }
    // java.util.stream.DoublePipeline.sum delegates to
    // Collectors.computeFinalSum: the compensation word is the negated
    // low-order error, so Java returns high - compensation rather than the
    // running high word alone.  The distinction is observable at an exact
    // calling-z boundary.  Keep Java's same-sign infinity fallback as well;
    // otherwise compensated arithmetic can spuriously turn a valid infinite
    // sum into NaN.
    const double compensated_sum = sum - compensation;
    return std::isnan(compensated_sum) && std::isinf(simple_sum)
        ? simple_sum
        : compensated_sum;
}

Statistics statistics(const std::vector<const Segment*>& segments) {
    // Preserve Java's IEEE-754 result for an empty neutral set.  GATK's
    // DoubleStream sums divide by a zero total length, yielding NaN mean/sd;
    // calls then fall through to NEUTRAL because both comparisons with NaN
    // are false.  Treating this valid case as an input error changes output.
    if (segments.empty()) return Statistics{};
    std::vector<double> lengths;
    std::vector<double> weighted;
    lengths.reserve(segments.size());
    weighted.reserve(segments.size());
    for (const auto* segment : segments) {
        const double length = static_cast<double>(segment->end - segment->start + 1);
        lengths.push_back(length);
        weighted.push_back(length * segment->copy_ratio);
    }
    const double total_length = kahan_sum(lengths);
    const double mean = kahan_sum(weighted) / total_length;
    std::vector<double> squares;
    squares.reserve(segments.size());
    for (const auto* segment : segments) {
        const double length = static_cast<double>(segment->end - segment->start + 1);
        squares.push_back(length * (segment->copy_ratio - mean) * (segment->copy_ratio - mean));
    }
    const double count = static_cast<double>(segments.size());
    if (!(total_length > 0.0) || !std::isfinite(total_length))
        throw std::runtime_error("BAD_INPUT: copy-ratio segment lengths are not positive");
    // The Java caller uses a weighted sample variance.  For one neutral
    // segment its denominator is exactly zero, so 0/0 is NaN (not zero).
    const double denominator = ((count - 1.0) / count) * total_length;
    const double variance = kahan_sum(squares) / denominator;
    // Clamp only finite negative round-off; std::max(0, NaN) would
    // accidentally turn Java's NaN into zero on common libstdc++ builds.
    const double nonnegative_variance =
        std::isfinite(variance) && variance < 0.0 ? 0.0 : variance;
    return Statistics{mean, std::sqrt(nonnegative_variance)};
}

std::vector<std::string> calls(const std::vector<Segment>& segments, double lower, double upper,
                               double outlier_z, double calling_z, Statistics& final_statistics) {
    std::vector<const Segment*> neutral;
    for (const auto& segment : segments)
        if (lower <= segment.copy_ratio && segment.copy_ratio <= upper) neutral.push_back(&segment);
    const auto unfiltered = statistics(neutral);
    std::vector<const Segment*> filtered;
    for (const auto* segment : neutral) {
        const bool inlier = unfiltered.standard_deviation == 0.0
            ? std::abs(segment->copy_ratio - unfiltered.mean) <= 1.0e-12
            : std::abs(segment->copy_ratio - unfiltered.mean) <=
                  unfiltered.standard_deviation * outlier_z;
        if (inlier) filtered.push_back(segment);
    }
    // Do not replace an empty filtered set with the unfiltered set.  Java's
    // second statistics pass legitimately returns NaN for singleton/empty
    // neutral sets; that NaN is part of the no-neutral calling contract.
    final_statistics = statistics(filtered);
    std::vector<std::string> result;
    result.reserve(segments.size());
    for (const auto& segment : segments) {
        if (lower <= segment.copy_ratio && segment.copy_ratio <= upper) {
            result.emplace_back("0");
            continue;
        }
        const double deviation = segment.copy_ratio - final_statistics.mean;
        if (deviation < -final_statistics.standard_deviation * calling_z) result.emplace_back("-");
        else if (deviation > final_statistics.standard_deviation * calling_z) result.emplace_back("+");
        else result.emplace_back("0");
    }
    return result;
}

struct CopyRatioKernelResult {
    std::vector<double> values;
    std::size_t simd_width = 1;
    std::size_t simd_groups = 0;
    std::size_t batches = 0;
    std::size_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

CopyRatioKernelResult compute_copy_ratios(const Table& table) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using Simd = Kokkos::Experimental::simd<double>;
    constexpr std::size_t width = Simd::size();
    Kokkos::View<double*> means("copy_ratio_means", table.segments.size());
    auto host_means = Kokkos::create_mirror_view(means);
    for (std::size_t i = 0; i < table.segments.size(); ++i) host_means(i) = table.segments[i].mean_log2;
    fastgatk::core::HostBatch host_batch("call-copy-ratio-v1");
    host_batch.records = table.segments.size();
    host_batch.bytes = table.segments.size() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("call-copy-ratio");
    kernel_plan.begin_prepare(host_batch);
    Kokkos::deep_copy(means, host_means);
    Kokkos::View<double*> ratios("copy_ratios", table.segments.size());
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(table.segments.size());
    device_batch.bind("means", means);
    device_batch.bind("ratios", ratios);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    const std::size_t groups = (table.segments.size() + width - 1) / width;
    kernel_plan.begin_execute();
    Kokkos::parallel_for("call_copy_ratio_segments_copy_ratio",
                         Kokkos::RangePolicy<ExecSpace>(0, groups),
                         KOKKOS_LAMBDA(const std::size_t group) {
                             const std::size_t base = group * width;
                             const Simd input([=](auto lane) -> double {
                                 const std::size_t index = base + static_cast<std::size_t>(lane);
                                 return index < means.extent(0) ? means(index) : 0.0;
                             });
                             // Keep the scalar base and operation identical to
                             // the historical path.  Kokkos SIMD's portable
                             // math wrapper lowers this to AVX2/AVX-512/NEON
                             // where available and remains scalar elsewhere.
                             const Simd output = Kokkos::pow(Simd(2.0), input);
                             for (std::size_t lane = 0; lane < width; ++lane) {
                                 const std::size_t index = base + lane;
                                 if (index < ratios.extent(0)) ratios(index) = output[lane];
                             }
                         });
    ExecSpace().fence();
    kernel_plan.end_execute();
    auto host_ratios = Kokkos::create_mirror_view(ratios);
    Kokkos::deep_copy(host_ratios, ratios);
    CopyRatioKernelResult result;
    result.values.resize(table.segments.size());
    result.simd_width = width;
    result.simd_groups = groups;
    result.batches = kernel_plan.telemetry().execute_calls;
    result.observations = table.segments.size();
    result.prepare_seconds = kernel_plan.telemetry().prepare_seconds;
    result.execute_seconds = kernel_plan.telemetry().execute_seconds;
    for (std::size_t i = 0; i < result.values.size(); ++i) result.values[i] = host_ratios(i);
    return result;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto c : text) {
        if (c == '"' || c == '\\') escaped << '\\';
        if (c == '\n') escaped << "\\n";
        else if (c == '\r') escaped << "\\r";
        else if (c == '\t') escaped << "\\t";
        else escaped << c;
    }
    return escaped.str();
}

std::string json_double(double value) {
    // IEEE-754 NaN/Inf are valid internal results for a degenerate GATK
    // neutral set, but JSON has no such numeric literals.  Keep the manifest
    // machine-readable and represent non-finite statistics as null.
    if (!std::isfinite(value)) return "null";
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return output.str();
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

void write_called(const std::string& path, const Table& table, const std::vector<std::string>& called) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write called segments: " + path);
    for (const auto& line : table.header_lines) output << line << '\n';
    output << "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL\n";
    output << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < table.segments.size(); ++i) {
        const auto& segment = table.segments[i];
        const auto format_mean = [](double value) {
            if (std::isnan(value)) return std::string("NaN");
            if (std::isinf(value)) return value > 0.0 ? std::string("Infinity") : std::string("-Infinity");
            std::ostringstream formatted;
            formatted << std::fixed << std::setprecision(6) << value;
            return formatted.str();
        };
        output << segment.contig << '\t' << segment.start << '\t' << segment.end << '\t'
               << segment.num_points << '\t' << format_mean(segment.mean_log2) << '\t' << called[i] << '\n';
    }
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize called segments: " + path);
}

void write_legacy(const std::string& path, const Table& table, const std::vector<std::string>& called) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write legacy segments: " + path);
    output << "Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n";
    output << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < table.segments.size(); ++i) {
        const auto& segment = table.segments[i];
        const auto format_mean = [](double value) {
            if (std::isnan(value)) return std::string("NaN");
            if (std::isinf(value)) return value > 0.0 ? std::string("Infinity") : std::string("-Infinity");
            std::ostringstream formatted;
            formatted << std::fixed << std::setprecision(6) << value;
            return formatted.str();
        };
        output << table.sample << '\t' << segment.contig << '\t' << segment.start << '\t'
               << segment.end << '\t' << segment.num_points << '\t' << called[i] << '\t'
               << format_mean(segment.mean_log2) << '\n';
    }
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize legacy segments: " + path);
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto started = std::chrono::steady_clock::now();
        auto options = parse(argc, argv);
        // GATK always writes the IGV-compatible sidecar next to -O.  Keep
        // --legacy-output as a compatibility override, but derive the
        // default using the same remove-extension + ".igv.seg" convention.
        if (options.legacy_output.empty())
            options.legacy_output = default_legacy_output(options.output);
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(static_cast<unsigned>(std::max(1, options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        auto table = read_table(options.input);

        const auto copy_ratio_kernel = compute_copy_ratios(table);
        const auto& copy_ratios = copy_ratio_kernel.values;
        for (std::size_t i = 0; i < table.segments.size(); ++i)
            table.segments[i].copy_ratio = copy_ratios[i];

        Statistics final_statistics;
        const auto called_values = calls(table.segments, options.lower, options.upper,
                                         options.outlier_z, options.calling_z, final_statistics);
        write_called(options.output, table, called_values);
        if (!options.legacy_output.empty()) write_legacy(options.legacy_output, table, called_values);
        if (!file_complete(options.output) || (!options.legacy_output.empty() && !file_complete(options.legacy_output)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: called segment output is incomplete");
        if (!options.manifest.empty()) {
            const auto file_bytes = [](const std::string& path) -> std::uintmax_t {
                if (path.empty()) return 0;
                std::error_code error;
                const auto bytes = std::filesystem::file_size(path, error);
                return error ? 0 : bytes;
            };
            const auto wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"CallCopyRatioSegments\","
                     << "\"implementation\":\"fastgatk-call-copy-ratio-segments\","
                     << "\"status\":\"prototype\",\"input\":\"" << json_escape(options.input)
                     << "\",\"output\":\"" << json_escape(options.output)
                     << "\",\"legacy_output\":\"" << json_escape(options.legacy_output)
                     << "\",\"segments\":" << table.segments.size()
                     << ",\"determinism\":\"strict\""
                     << ",\"neutral_lower\":" << options.lower << ",\"neutral_upper\":" << options.upper
                     << ",\"outlier_z\":" << options.outlier_z << ",\"calling_z\":" << options.calling_z
                     << ",\"statistics_mean\":" << json_double(final_statistics.mean)
                     << ",\"statistics_sd\":" << json_double(final_statistics.standard_deviation)
                     << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"called-segments\",\"complete\":"
                     << (file_complete(options.output) ? "true" : "false") << "}"
                     << (options.legacy_output.empty() ? "" : ",{\"path\":\"" + json_escape(options.legacy_output)
                         + "\",\"kind\":\"legacy-segments\",\"complete\":"
                         + (file_complete(options.legacy_output) ? "true}" : "false}"))
                     << "],\"telemetry\":{\"output_bytes\":" << file_bytes(options.output)
                     << ",\"legacy_output_bytes\":" << file_bytes(options.legacy_output)
                     << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                     << ",\"kernel_execution_policy\":\"RangePolicy\""
                     << ",\"kernel_batches\":" << copy_ratio_kernel.batches
                     << ",\"kernel_observations\":" << copy_ratio_kernel.observations
                     << ",\"kernel_prepare_seconds\":" << copy_ratio_kernel.prepare_seconds
                     << ",\"copy_ratio_kernel_execution_space\":\""
                     << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"copy_ratio_kernel_simd_width\":" << copy_ratio_kernel.simd_width
                     << ",\"copy_ratio_kernel_simd_groups\":" << copy_ratio_kernel.simd_groups
                     << ",\"copy_ratio_kernel_execute_seconds\":" << copy_ratio_kernel.execute_seconds
                     << ",\"wall_seconds\":" << wall_seconds << "}}\n";
            manifest.close();
            if (!manifest || !file_complete(options.manifest))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest: " + options.manifest);
        }
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"CallCopyRatioSegments\",\"status\":\"prototype\","
                  << "\"segments\":" << table.segments.size()
                  << ",\"output\":\"" << json_escape(options.output) << "\"}\n";
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << error.what() << '\n';
        return 2;
    }
}
