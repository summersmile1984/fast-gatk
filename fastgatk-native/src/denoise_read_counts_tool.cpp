#include <Kokkos_Core.hpp>
#include "fastgatk/hdf5_count_collection.hpp"
#include "fastgatk/core/plan.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double EPSILON = 1.0e-9;
constexpr double INV_LOG_2 = 1.44269504088896340735992468100189214;

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

struct Options {
    std::string input;
    std::string output;
    std::string standardized_output;
    std::string manifest;
    std::string panel_of_normals;
    std::string annotated_intervals;
    std::string format = "TSV";
    std::string sample;
    int svd_rank = -1;
    int threads = 1;
};

struct Interval {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    double count = 0.0;
};

struct CountTable {
    std::vector<std::string> header_lines;
    std::vector<Interval> intervals;
    std::string sample = "UNKNOWN";
};

struct PanelTable {
    std::vector<std::string> header_lines;
    std::vector<Interval> intervals;
    std::vector<std::string> samples;
    std::vector<std::vector<double>> counts;
    bool wide = false;
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

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-denoise-read-counts (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                    SimpleCountCollection TSV\n"
                         "  -O, --output FILE                   denoised copy-ratio TSV\n"
                         "      --standardized-copy-ratios FILE standardized copy-ratio TSV\n"
                         "      --panel-of-normals FILE          TSV/legacy PoN (interval-matched normalization)\n"
                         "      --count-panel-of-normals FILE    GATK HDF5 SVD PoN v7\n"
                         "      --number-of-eigensamples N       remove N deterministic PoN SVD components\n"
                         "      --svd-rank N                     alias for --number-of-eigensamples\n"
                         "      --format TSV|HDF5                input count collection format\n"
                         "      --sample NAME                    sample metadata override\n"
                         "      --output-manifest FILE           OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-O" || has_option(argument, "--output") ||
                   has_option(argument, "--denoised-copy-ratios")) {
            const char* name = argument.rfind("--denoised-copy-ratios", 0) == 0
                ? "--denoised-copy-ratios" : "--output";
            options.output = require_value(index, argc, argv, argument, name,
                                           argument == "-O" ? "-O" : nullptr);
        } else if (has_option(argument, "--standardized-copy-ratios")) {
            options.standardized_output = require_value(
                index, argc, argv, argument, "--standardized-copy-ratios");
        } else if (has_option(argument, "--panel-of-normals") || has_option(argument, "--count-panel-of-normals")) {
            const char* name = argument.rfind("--count-panel-of-normals", 0) == 0
                ? "--count-panel-of-normals" : "--panel-of-normals";
            options.panel_of_normals = require_value(index, argc, argv, argument, name);
        } else if (has_option(argument, "--format")) {
            options.format = require_value(index, argc, argv, argument, "--format");
        } else if (has_option(argument, "--sample")) {
            options.sample = require_value(index, argc, argv, argument, "--sample");
        } else if (has_option(argument, "--number-of-eigensamples") || has_option(argument, "--svd-rank")) {
            const char* name = argument.rfind("--svd-rank", 0) == 0
                ? "--svd-rank" : "--number-of-eigensamples";
            options.svd_rank = std::stoi(require_value(index, argc, argv, argument, name));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(
                index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--disable-tool-default-read-filters" ||
                   argument == "--dont-trim-intervals") {
            // Accepted compatibility switches with no effect on a count table.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates") ||
                   has_option(argument, "--normalization-target-coverage") ||
                   has_option(argument, "--annotated-intervals") ||
                   has_option(argument, "--sequence-dictionary")) {
            // Keep the parser deterministic, but make it clear that rich GATK metadata
            // is not silently consumed by this TSV-only path.
            const auto metadata_value = require_value(index, argc, argv, argument,
                argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                argument.rfind("--seconds-between-progress-updates", 0) == 0
                    ? "--seconds-between-progress-updates" :
                argument.rfind("--annotated-intervals", 0) == 0 ? "--annotated-intervals" :
                argument.rfind("--sequence-dictionary", 0) == 0
                    ? "--sequence-dictionary" : "--normalization-target-coverage");
            if (argument.rfind("--annotated-intervals", 0) == 0)
                options.annotated_intervals = metadata_value;
            if (argument.rfind("--panel-of-normals", 0) != 0 &&
                argument.rfind("--normalization-target-coverage", 0) == 0)
                throw std::invalid_argument("BACKEND_UNAVAILABLE: normalization target coverage requires GATK fallback");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty() || options.standardized_output.empty())
        throw std::invalid_argument(
            "-O/--output (denoised copy ratios) and --standardized-copy-ratios are required");
    for (auto& c : options.format) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (options.format == "HDF5" && !fastgatk::hdf5::available())
        throw std::invalid_argument("BACKEND_UNAVAILABLE: HDF5 support was not built; use --format TSV or --fallback");
    if (!options.panel_of_normals.empty() &&
        (options.panel_of_normals.ends_with(".h5") || options.panel_of_normals.ends_with(".hdf5")) &&
        !fastgatk::hdf5::available())
        throw std::invalid_argument("BACKEND_UNAVAILABLE: HDF5 panel-of-normals support was not built; use --fallback");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    if (options.svd_rank < -1) throw std::invalid_argument("--number-of-eigensamples must be non-negative");
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

// SimpleCountCollection's TSV codec decodes COUNT with DataLine.getInt(),
// rather than getDouble().  Keeping the value as a double after accepting a
// fractional token would silently change the fractional-coverage denominator
// relative to GATK (and would accept files that Java rejects).  Parse the
// lexical integer here so values such as "10.0", "1e1", and out-of-range
// counts fail at the same input boundary.
double parse_integer_count(const std::string& text, const std::string& path,
                           std::size_t line_number, const char* column = "COUNT") {
    if (text.empty())
        throw std::runtime_error("BAD_INPUT: " + std::string(column) +
                                 " must be a non-negative integer at " + path + ":" +
                                 std::to_string(line_number));
    std::size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(text, &consumed, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("BAD_INPUT: " + std::string(column) +
                                 " must be a non-negative integer at " + path + ":" +
                                 std::to_string(line_number));
    }
    if (consumed != text.size() || value < 0 || value > INT_MAX)
        throw std::runtime_error("BAD_INPUT: " + std::string(column) +
                                 " must be a non-negative integer at " + path + ":" +
                                 std::to_string(line_number));
    return static_cast<double>(value);
}

CountTable read_counts(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open count table: " + path);
    CountTable table;
    bool columns_seen = false;
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
        if (line[0] == '@') {
            table.header_lines.push_back(line);
            continue;
        }
        const auto fields = split_tab(line);
        if (!columns_seen) {
            if (fields.size() != 4 || fields[0] != "CONTIG" || fields[1] != "START" ||
                fields[2] != "END" || fields[3] != "COUNT")
                throw std::runtime_error("BAD_INPUT: expected CONTIG\\tSTART\\tEND\\tCOUNT header at " +
                                         path + ":" + std::to_string(line_number));
            columns_seen = true;
            continue;
        }
        if (fields.size() != 4)
            throw std::runtime_error("BAD_INPUT: invalid count row at " + path + ":" + std::to_string(line_number));
        Interval interval;
        interval.contig = fields[0];
        interval.start = std::stoll(fields[1]);
        interval.end = std::stoll(fields[2]);
        interval.count = parse_integer_count(fields[3], path, line_number);
        if (interval.contig.empty() || interval.start < 1 || interval.end < interval.start)
            throw std::runtime_error("BAD_INPUT: invalid count values at " + path + ":" + std::to_string(line_number));
        table.intervals.push_back(std::move(interval));
    }
    if (!columns_seen || table.intervals.empty())
        throw std::runtime_error("BAD_INPUT: count table has no intervals: " + path);
    return table;
}

CountTable read_counts_any(const std::string& path, bool hdf5) {
    if (!hdf5) return read_counts(path);
    const auto collection = fastgatk::hdf5::read_simple_count_collection(path);
    if (collection.intervals.empty() || collection.intervals.size() != collection.counts.size())
        throw std::runtime_error("BAD_INPUT: HDF5 count collection is empty or misaligned");
    CountTable table;
    table.sample = collection.sample;
    // SimpleCountCollection HDF5 stores the complete SAM sequence dictionary
    // separately from the interval/count matrices.  DenoiseReadCounts writes
    // CopyRatioCollection metadata from the input collection, so carry this
    // header through the HDF5 path instead of silently dropping all @SQ
    // records from the two TSV outputs.  This is observable downstream: GATK
    // readers use the dictionary for interval ordering and metadata checks.
    std::size_t begin = 0;
    while (begin < collection.sequence_dictionary.size()) {
        const auto newline = collection.sequence_dictionary.find('\n', begin);
        const auto limit = newline == std::string::npos
            ? collection.sequence_dictionary.size() : newline;
        auto line = collection.sequence_dictionary.substr(begin, limit - begin);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) table.header_lines.push_back(std::move(line));
        if (newline == std::string::npos) break;
        begin = newline + 1;
    }
    for (std::size_t index = 0; index < collection.intervals.size(); ++index) {
        const auto& source = collection.intervals[index];
        table.intervals.push_back({source.contig, source.start, source.end, collection.counts[index]});
    }
    return table;
}

struct AnnotatedIntervals {
    std::vector<Interval> intervals;
    std::vector<double> gc_content;
};

AnnotatedIntervals read_annotated_intervals(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open annotated intervals: " + path);
    AnnotatedIntervals result;
    bool header_seen = false;
    std::size_t line_number = 0;
    std::size_t gc_column = 3;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty() || line[0] == '@') continue;
        const auto fields = split_tab(line);
        if (!header_seen) {
            if (fields.size() < 4 || fields[0] != "CONTIG" || fields[1] != "START" || fields[2] != "END")
                throw std::runtime_error("BAD_INPUT: annotated intervals must contain CONTIG/START/END/GC_CONTENT: " + path);
            gc_column = 0;
            while (gc_column < fields.size() && fields[gc_column] != "GC_CONTENT") ++gc_column;
            if (gc_column == fields.size()) throw std::runtime_error("BAD_INPUT: annotated intervals lack GC_CONTENT: " + path);
            header_seen = true;
            continue;
        }
        if (fields.size() <= gc_column) throw std::runtime_error("BAD_INPUT: malformed annotated interval at " + path + ":" + std::to_string(line_number));
        const auto start = std::stoll(fields[1]);
        const auto end = std::stoll(fields[2]);
        const auto gc = std::stod(fields[gc_column]);
        if (fields[0].empty() || start < 1 || end < start ||
            !(std::isnan(gc) || (std::isfinite(gc) && gc >= 0.0 && gc <= 1.0)))
            throw std::runtime_error("BAD_INPUT: invalid annotated interval at " + path + ":" + std::to_string(line_number));
        result.intervals.push_back({fields[0], start, end, 0.0});
        result.gc_content.push_back(gc);
    }
    if (!header_seen || result.intervals.empty()) throw std::runtime_error("BAD_INPUT: annotated intervals contain no rows: " + path);
    return result;
}

PanelTable read_panel(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open panel-of-normals: " + path);
    PanelTable table;
    bool columns_seen = false;
    std::size_t line_number = 0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty()) continue;
        if (line[0] == '@') {
            table.header_lines.push_back(line);
            continue;
        }
        const auto fields = split_tab(line);
        if (!columns_seen) {
            if (fields.size() < 4 || fields[0] != "CONTIG" || fields[1] != "START" ||
                fields[2] != "END")
                throw std::runtime_error("BAD_INPUT: expected CONTIG\\tSTART\\tEND\\tCOUNT or sample columns at " +
                                         path + ":" + std::to_string(line_number));
            if (fields.size() == 4 && fields[3] == "COUNT") {
                table.samples = {"PANEL"};
                table.wide = false;
            } else {
                table.wide = true;
                for (std::size_t index = 3; index < fields.size(); ++index) {
                    if (fields[index].empty())
                        throw std::runtime_error("BAD_INPUT: empty panel sample name at " + path);
                    table.samples.push_back(fields[index]);
                }
            }
            columns_seen = true;
            continue;
        }
        if (fields.size() != table.samples.size() + 3)
            throw std::runtime_error("BAD_INPUT: invalid panel row at " + path + ":" +
                                     std::to_string(line_number));
        Interval interval;
        interval.contig = fields[0];
        interval.start = std::stoll(fields[1]);
        interval.end = std::stoll(fields[2]);
        if (interval.contig.empty() || interval.start < 1 || interval.end < interval.start)
            throw std::runtime_error("BAD_INPUT: invalid panel interval at " + path + ":" +
                                     std::to_string(line_number));
        interval.count = std::stod(fields[3]);
        if (!std::isfinite(interval.count) || interval.count < 0.0)
            throw std::runtime_error("BAD_INPUT: invalid panel count at " + path + ":" +
                                     std::to_string(line_number));
        table.intervals.push_back(interval);
        if (table.counts.empty()) table.counts.resize(table.samples.size());
        for (std::size_t sample = 0; sample < table.samples.size(); ++sample) {
            const double value = std::stod(fields[sample + 3]);
            if (!std::isfinite(value) || value < 0.0)
                throw std::runtime_error("BAD_INPUT: invalid panel count at " + path + ":" +
                                         std::to_string(line_number));
            table.counts[sample].push_back(value);
        }
    }
    if (!columns_seen || table.intervals.empty() || table.samples.empty())
        throw std::runtime_error("BAD_INPUT: panel-of-normals has no intervals: " + path);
    return table;
}

double median(std::vector<double> values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    if (values.size() % 2 == 0) return 0.5 * (values[middle - 1] + values[middle]);
    return values[middle];
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

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

int gc_bin(double gc_content) {
    if (!std::isfinite(gc_content)) return 0;
    return std::clamp(static_cast<int>(std::llround(gc_content * 100.0)), 0, 100);
}

void correct_gc_bias(std::vector<double>& values, const std::vector<double>& gc_content) {
    if (gc_content.empty()) return;
    if (values.size() != gc_content.size())
        throw std::invalid_argument("BAD_INPUT: GC annotations do not match read-count intervals");
    constexpr int bins = 101;
    constexpr double decay_rate = 1.0 / (0.02 * bins);
    const double total_before = std::accumulate(values.begin(), values.end(), 0.0);
    std::vector<std::vector<double>> by_gc(bins);
    for (std::size_t index = 0; index < values.size(); ++index)
        by_gc[gc_bin(gc_content[index])].push_back(values[index]);
    std::vector<double> medians(bins, 1.0), factors(bins, 1.0);
    for (int bin = 0; bin < bins; ++bin)
        if (!by_gc[bin].empty()) medians[bin] = median(by_gc[bin]);
    for (int bin = 0; bin < bins; ++bin) {
        double weighted_sum = 0.0, weight_total = 0.0;
        for (int source = 0; source < bins; ++source) {
            const double weight = static_cast<double>(by_gc[source].size()) *
                std::exp(-std::abs(bin - source) * decay_rate);
            weighted_sum += weight * medians[source];
            weight_total += weight;
        }
        if (weight_total > 0.0 && std::isfinite(weighted_sum / weight_total) && weighted_sum != 0.0)
            factors[bin] = 1.0 / (weighted_sum / weight_total);
    }
    double total_after = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] *= factors[gc_bin(gc_content[index])];
        total_after += values[index];
    }
    if (total_before > 0.0 && total_after > 0.0) {
        const double normalization = total_before / total_after;
        for (auto& value : values) value *= normalization;
    }
}

void write_table(const std::string& path, const CountTable& input, const std::string& sample,
                 const std::vector<double>& values) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write copy ratios: " + path);
    output << "@HD\tVN:1.6\n";
    for (const auto& line : input.header_lines)
        if (line.rfind("@SQ", 0) == 0) output << line << '\n';
    output << "@RG\tID:GATKCopyNumber\tSM:" << sample << '\n';
    output << "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n";
    output << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < input.intervals.size(); ++i) {
        const auto& interval = input.intervals[i];
        output << interval.contig << '\t' << interval.start << '\t' << interval.end << '\t'
               << values[i] << '\n';
    }
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize copy ratios: " + path);
}

std::vector<double> standardize(const std::vector<double>& counts,
                                const std::vector<double>* panel_counts,
                                int requested_threads,
                                const std::vector<double>* gc_content = nullptr,
                                KernelStats* kernel_stats = nullptr) {
    double total = 0.0;
    for (const auto count : counts) total += count;
    if (!(total > 0.0))
        throw std::invalid_argument("java.lang.IllegalArgumentException: Sample does not have a positive sample median.");
    std::vector<double> fractional(counts.size());
    if (panel_counts == nullptr) {
        for (std::size_t i = 0; i < counts.size(); ++i) fractional[i] = counts[i] / total;
    } else {
        if (panel_counts->size() != counts.size())
            throw std::invalid_argument("BAD_INPUT: panel-of-normals interval count does not match sample");
        double panel_total = 0.0;
        for (const auto count : *panel_counts) panel_total += count;
        if (!(panel_total > 0.0))
            throw std::invalid_argument("BAD_INPUT: panel-of-normals has no positive counts");
        std::vector<double> panel_fractional(panel_counts->size());
        for (std::size_t i = 0; i < counts.size(); ++i) {
            fractional[i] = counts[i] / total;
            panel_fractional[i] = (*panel_counts)[i] / panel_total;
        }
        if (gc_content != nullptr) {
            correct_gc_bias(fractional, *gc_content);
            correct_gc_bias(panel_fractional, *gc_content);
        }
        for (std::size_t i = 0; i < counts.size(); ++i)
            fractional[i] /= std::max(EPSILON, panel_fractional[i]);
    }
    if (gc_content != nullptr && panel_counts == nullptr) correct_gc_bias(fractional, *gc_content);
    const double sample_median = median(fractional);
    if (!(sample_median > 0.0))
        throw std::invalid_argument("java.lang.IllegalArgumentException: Sample does not have a positive sample median.");

    Kokkos::View<double*> device_fractional("denoise_fractional", counts.size());
    Kokkos::View<double*> device_logs("denoise_logs", counts.size());
    Kokkos::View<double*> device_values("denoise_values", counts.size());
    auto host_fractional = Kokkos::create_mirror_view(device_fractional);
    for (std::size_t i = 0; i < fractional.size(); ++i) host_fractional(i) = fractional[i];
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch standardize_host("denoise-standardize-v1");
    standardize_host.records = counts.size();
    standardize_host.bytes = counts.size() * 3 * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> standardize_plan("denoise-standardize");
    standardize_plan.begin_prepare(standardize_host);
    Kokkos::deep_copy(device_fractional, host_fractional);
    fastgatk::core::DeviceBatch<ExecSpace> standardize_device(counts.size());
    standardize_device.bind("fractional", device_fractional);
    standardize_device.bind("logs", device_logs);
    standardize_device.bind("values", device_values);
    ExecSpace().fence();
    standardize_plan.end_prepare(standardize_device);
    const double log_floor = std::log(EPSILON) * INV_LOG_2;
    const auto policy = Kokkos::RangePolicy<ExecSpace>(0, counts.size());
    standardize_plan.begin_execute();
    Kokkos::parallel_for("denoise_read_counts_log2", policy, KOKKOS_LAMBDA(const std::size_t i) {
        const double ratio = device_fractional(i) / sample_median;
        device_logs(i) = ratio < EPSILON ? log_floor : Kokkos::log(ratio) * INV_LOG_2;
    });
    ExecSpace().fence();
    standardize_plan.end_execute();
    auto host_logs = Kokkos::create_mirror_view(device_logs);
    Kokkos::deep_copy(host_logs, device_logs);
    std::vector<double> logs(counts.size());
    for (std::size_t i = 0; i < logs.size(); ++i) logs[i] = host_logs(i);
    const double log_median = median(logs);
    standardize_plan.begin_execute();
    Kokkos::parallel_for("denoise_read_counts_center", policy, KOKKOS_LAMBDA(const std::size_t i) {
        device_values(i) = device_logs(i) - log_median;
    });
    ExecSpace().fence();
    standardize_plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += counts.size();
        kernel_stats->prepare_seconds += standardize_plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += standardize_plan.telemetry().execute_seconds;
    }
    auto host_values = Kokkos::create_mirror_view(device_values);
    Kokkos::deep_copy(host_values, device_values);
    std::vector<double> result(counts.size());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = host_values(i);
    (void)requested_threads;
    return result;
}

std::vector<double> denoise_svd(const std::vector<double>& target,
                                const std::vector<std::vector<double>>& panel_values,
                                int requested_rank,
                                int requested_threads,
                                KernelStats* kernel_stats = nullptr) {
    if (panel_values.empty()) return target;
    const std::size_t intervals = target.size();
    const std::size_t samples = panel_values.size();
    for (const auto& column : panel_values)
        if (column.size() != intervals)
            throw std::invalid_argument("BAD_INPUT: panel-of-normals matrix dimensions do not match sample");
    const int rank_limit = static_cast<int>(std::min(intervals, samples));
    const int rank = requested_rank < 0 ? std::min(20, rank_limit)
                                       : std::min(requested_rank, rank_limit);
    if (rank == 0) return target;

    // X is an interval x panel-sample matrix.  This is the dense GEMM boundary:
    // the covariance X^T X and all projections are evaluated through Kokkos
    // Views/parallel_for so the same code maps to Serial/OpenMP/CUDA/HIP.
    Kokkos::View<double**> device_x("denoise_pon_matrix", intervals, samples);
    Kokkos::View<double**> device_gram("denoise_pon_gram", samples, samples);
    Kokkos::View<double*> device_target("denoise_target", intervals);
    auto host_x = Kokkos::create_mirror_view(device_x);
    auto host_target = Kokkos::create_mirror_view(device_target);
    for (std::size_t i = 0; i < intervals; ++i) {
        host_target(i) = target[i];
        for (std::size_t j = 0; j < samples; ++j) host_x(i, j) = panel_values[j][i];
    }
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch gram_host("denoise-pon-gram-v1");
    gram_host.records = samples * samples;
    gram_host.bytes = device_x.span() * sizeof(double) + device_gram.span() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> gram_plan("denoise-pon-gram");
    gram_plan.begin_prepare(gram_host);
    Kokkos::deep_copy(device_x, host_x);
    Kokkos::deep_copy(device_target, host_target);
    fastgatk::core::DeviceBatch<ExecSpace> gram_device(samples * samples);
    gram_device.bind("matrix", device_x);
    gram_device.bind("gram", device_gram);
    gram_device.bind("target", device_target);
    ExecSpace().fence();
    gram_plan.end_prepare(gram_device);
    gram_plan.begin_execute();
    Kokkos::parallel_for("denoise_pon_gemm", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>(
        {0, 0}, {samples, samples}), KOKKOS_LAMBDA(const std::size_t row, const std::size_t col) {
        double sum = 0.0;
        for (std::size_t i = 0; i < intervals; ++i) sum += device_x(i, row) * device_x(i, col);
        device_gram(row, col) = sum;
    });
    ExecSpace().fence();
    gram_plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += samples * samples;
        kernel_stats->prepare_seconds += gram_plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += gram_plan.telemetry().execute_seconds;
    }
    auto host_gram = Kokkos::create_mirror_view(device_gram);
    Kokkos::deep_copy(host_gram, device_gram);
    std::vector<double> gram(samples * samples, 0.0);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t col = 0; col < samples; ++col) gram[row * samples + col] = host_gram(row, col);

    std::vector<double> denoised = target;
    std::vector<double> vector(samples, 0.0), next(samples, 0.0);
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    for (int component = 0; component < rank; ++component) {
        for (std::size_t j = 0; j < samples; ++j)
            vector[j] = 1.0 + 0.1732050807568877 * static_cast<double>((j + 1) * (component + 1));
        for (const auto& prior : eigenvectors) {
            double projection = 0.0;
            for (std::size_t j = 0; j < samples; ++j) projection += vector[j] * prior[j];
            for (std::size_t j = 0; j < samples; ++j) vector[j] -= projection * prior[j];
        }
        double norm = 0.0;
        for (const auto value : vector) norm += value * value;
        norm = std::sqrt(norm);
        if (!(norm > EPSILON)) break;
        for (auto& value : vector) value /= norm;
        double eigenvalue = 0.0;
        for (int iteration = 0; iteration < 64; ++iteration) {
            std::fill(next.begin(), next.end(), 0.0);
            for (std::size_t row = 0; row < samples; ++row)
                for (std::size_t col = 0; col < samples; ++col)
                    next[row] += gram[row * samples + col] * vector[col];
            for (const auto& prior : eigenvectors) {
                double projection = 0.0;
                for (std::size_t j = 0; j < samples; ++j) projection += next[j] * prior[j];
                for (std::size_t j = 0; j < samples; ++j) next[j] -= projection * prior[j];
            }
            norm = 0.0;
            for (const auto value : next) norm += value * value;
            norm = std::sqrt(norm);
            if (!(norm > EPSILON)) break;
            for (std::size_t j = 0; j < samples; ++j) vector[j] = next[j] / norm;
        }
        for (std::size_t row = 0; row < samples; ++row)
            for (std::size_t col = 0; col < samples; ++col)
                eigenvalue += vector[row] * gram[row * samples + col] * vector[col];
        if (!(eigenvalue > EPSILON) || !std::isfinite(eigenvalue)) break;
        eigenvalues.push_back(eigenvalue);
        eigenvectors.push_back(vector);
        for (std::size_t row = 0; row < samples; ++row)
            for (std::size_t col = 0; col < samples; ++col)
                gram[row * samples + col] -= eigenvalue * vector[row] * vector[col];
    }

    Kokkos::View<double*> device_denoised("denoise_svd_result", intervals);
    Kokkos::deep_copy(device_denoised, device_target);
    auto host_denoised = Kokkos::create_mirror_view(device_denoised);
    for (std::size_t component = 0; component < eigenvectors.size(); ++component) {
        const auto& right = eigenvectors[component];
        const double scale = std::sqrt(eigenvalues[component]);
        Kokkos::View<double*> device_right("denoise_svd_right", samples);
        Kokkos::View<double*> device_left("denoise_svd_left", intervals);
        Kokkos::View<double*> device_score("denoise_svd_score", 1);
        auto host_right = Kokkos::create_mirror_view(device_right);
        for (std::size_t j = 0; j < samples; ++j) host_right(j) = right[j];
        fastgatk::core::HostBatch component_host("denoise-pon-component-v1");
        component_host.records = intervals;
        component_host.bytes = (samples + intervals + 1) * sizeof(double);
        fastgatk::core::KernelPlan<ExecSpace> component_plan("denoise-pon-component");
        component_plan.begin_prepare(component_host);
        Kokkos::deep_copy(device_right, host_right);
        fastgatk::core::DeviceBatch<ExecSpace> component_device(intervals);
        component_device.bind("matrix", device_x);
        component_device.bind("right", device_right);
        component_device.bind("left", device_left);
        component_device.bind("score", device_score);
        component_device.bind("denoised", device_denoised);
        ExecSpace().fence();
        component_plan.end_prepare(component_device);
        component_plan.begin_execute();
        Kokkos::parallel_for("denoise_pon_left_projection", Kokkos::RangePolicy<ExecSpace>(0, intervals), KOKKOS_LAMBDA(const std::size_t i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < samples; ++j) sum += device_x(i, j) * device_right(j);
            device_left(i) = sum / scale;
        });
        component_plan.end_execute();
        component_plan.begin_execute();
        Kokkos::parallel_for("denoise_pon_score", Kokkos::RangePolicy<ExecSpace>(0, 1), KOKKOS_LAMBDA(const std::size_t) {
            double sum = 0.0;
            for (std::size_t i = 0; i < intervals; ++i) sum += device_left(i) * device_denoised(i);
            device_score(0) = sum;
        });
        ExecSpace().fence();
        component_plan.end_execute();
        auto host_score = Kokkos::create_mirror_view(device_score);
        Kokkos::deep_copy(host_score, device_score);
        const double score = host_score(0);
        component_plan.begin_execute();
        Kokkos::parallel_for("denoise_pon_reconstruct", Kokkos::RangePolicy<ExecSpace>(0, intervals), KOKKOS_LAMBDA(const std::size_t i) {
            device_denoised(i) -= device_left(i) * score;
        });
        ExecSpace().fence();
        component_plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += intervals;
            kernel_stats->prepare_seconds += component_plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += component_plan.telemetry().execute_seconds;
        }
    }
    Kokkos::deep_copy(host_denoised, device_denoised);
    for (std::size_t i = 0; i < intervals; ++i) denoised[i] = host_denoised(i);
    (void)requested_threads;
    return denoised;
}

std::vector<double> standardize_with_hdf5_panel(const CountTable& input,
                                                const fastgatk::hdf5::PanelOfNormals& panel,
                                                CountTable& output_input,
                                                const std::vector<double>* annotated_gc_content = nullptr,
                                                KernelStats* kernel_stats = nullptr) {
    if (panel.original_intervals.empty() || panel.panel_intervals.empty() ||
        panel.panel_interval_fractional_medians.empty())
        throw std::invalid_argument("BAD_INPUT: HDF5 panel of normals has no usable intervals");
    if (panel.original_intervals.size() != panel.original_read_counts.front().size())
        throw std::invalid_argument("BAD_INPUT: HDF5 panel original interval metadata is misaligned");

    // GATK's SVDReadCountPanelOfNormals validates the case collection against
    // the *original* interval list before it filters the panel intervals.  It
    // requires SimpleInterval list equality (same contig/start/end and order),
    // rather than accepting a case that merely contains every surviving panel
    // interval.  This is a critical workflow boundary: accepting a subset or
    // an extra interval would change the sample fractional-coverage denominator
    // and produce a successful-looking, but non-GATK, denoised copy-ratio file.
    if (input.intervals.size() != panel.original_intervals.size())
        throw std::invalid_argument(
            "BAD_INPUT: sample read-count intervals must be identical to HDF5 panel original intervals");
    for (std::size_t index = 0; index < input.intervals.size(); ++index) {
        const auto& actual = input.intervals[index];
        const auto& expected = panel.original_intervals[index];
        if (actual.contig != expected.contig || actual.start != expected.start ||
            actual.end != expected.end)
            throw std::invalid_argument(
                "BAD_INPUT: sample read-count intervals must be identical to HDF5 panel original intervals");
    }
    std::vector<double> counts(panel.panel_intervals.size(), 0.0);
    std::vector<std::size_t> input_indices;
    input_indices.reserve(panel.panel_intervals.size());
    for (const auto& wanted : panel.panel_intervals) {
        const auto found = std::find_if(input.intervals.begin(), input.intervals.end(),
            [&](const Interval& actual) {
                return actual.contig == wanted.contig && actual.start == wanted.start && actual.end == wanted.end;
            });
        if (found == input.intervals.end())
            throw std::invalid_argument("BAD_INPUT: case read-count intervals do not cover HDF5 panel intervals");
        input_indices.push_back(static_cast<std::size_t>(found - input.intervals.begin()));
    }
    double total = 0.0;
    for (const auto& interval : input.intervals) total += interval.count;
    if (!(total > 0.0)) throw std::invalid_argument("java.lang.IllegalArgumentException: Sample does not have a positive sample median.");
    std::vector<double> fractional(input.intervals.size());
    for (std::size_t index = 0; index < input.intervals.size(); ++index)
        fractional[index] = input.intervals[index].count / total;
    // GATK deliberately ignores --annotated-intervals whenever a panel of
    // normals is supplied: the panel carries the GC correction model that was
    // fitted during CreateReadCountPanelOfNormals.  Keep the pointer in the
    // signature for the no-panel/legacy call sites, but never let an ad-hoc
    // annotation override the panel metadata here.
    (void)annotated_gc_content;
    const auto& panel_gc = panel.original_interval_gc_content;
    if (!panel_gc.empty()) {
        if (panel.original_intervals.size() != input.intervals.size() ||
            panel_gc.size() != input.intervals.size())
            throw std::invalid_argument("BAD_INPUT: HDF5 panel GC annotations do not match original case intervals");
        for (std::size_t index = 0; index < input.intervals.size(); ++index) {
            if (input.intervals[index].contig != panel.original_intervals[index].contig ||
                input.intervals[index].start != panel.original_intervals[index].start ||
                input.intervals[index].end != panel.original_intervals[index].end)
                throw std::invalid_argument("BAD_INPUT: case intervals do not match HDF5 panel original intervals");
        }
        correct_gc_bias(fractional, panel_gc);
    }
    for (std::size_t index = 0; index < input_indices.size(); ++index) {
        const auto source = input_indices[index];
        const double interval_median = panel.panel_interval_fractional_medians[index];
        if (!(interval_median > 0.0) || !std::isfinite(interval_median))
            throw std::invalid_argument("BAD_INPUT: HDF5 panel interval median is not positive");
        counts[index] = fractional[source] / interval_median;
    }
    const double sample_median = median(counts);
    if (!(sample_median > 0.0))
        throw std::invalid_argument("java.lang.IllegalArgumentException: Sample does not have a positive sample median.");
    std::vector<double> values(counts.size());
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<double*> device_counts("denoise-hdf5-counts", counts.size());
    Kokkos::View<double*> device_values("denoise-hdf5-values", counts.size());
    auto host_counts = Kokkos::create_mirror_view(device_counts);
    for (std::size_t index = 0; index < counts.size(); ++index) host_counts(index) = counts[index];
    fastgatk::core::HostBatch standardize_host("denoise-hdf5-standardize-v1");
    standardize_host.records = counts.size();
    standardize_host.bytes = 2 * counts.size() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> standardize_plan("denoise-hdf5-standardize");
    standardize_plan.begin_prepare(standardize_host);
    Kokkos::deep_copy(device_counts, host_counts);
    fastgatk::core::DeviceBatch<ExecSpace> standardize_device(counts.size());
    standardize_device.bind("counts", device_counts);
    standardize_device.bind("values", device_values);
    ExecSpace().fence();
    standardize_plan.end_prepare(standardize_device);
    const double log_floor = std::log(EPSILON) * INV_LOG_2;
    standardize_plan.begin_execute();
    Kokkos::parallel_for("denoise_hdf5_standardize_log2", Kokkos::RangePolicy<ExecSpace>(0, counts.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            const double ratio = device_counts(index) / sample_median;
            device_values(index) = ratio < EPSILON ? log_floor : Kokkos::log(ratio) * INV_LOG_2;
        });
    ExecSpace().fence();
    standardize_plan.end_execute();
    auto host_values = Kokkos::create_mirror_view(device_values);
    Kokkos::deep_copy(host_values, device_values);
    for (std::size_t index = 0; index < counts.size(); ++index) values[index] = host_values(index);
    const double center = median(values);
    standardize_plan.begin_execute();
    Kokkos::parallel_for("denoise_hdf5_standardize_center", Kokkos::RangePolicy<ExecSpace>(0, counts.size()),
        KOKKOS_LAMBDA(const std::size_t index) { device_values(index) -= center; });
    ExecSpace().fence();
    standardize_plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += counts.size();
        kernel_stats->prepare_seconds += standardize_plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += standardize_plan.telemetry().execute_seconds;
    }
    Kokkos::deep_copy(host_values, device_values);
    for (std::size_t index = 0; index < counts.size(); ++index) values[index] = host_values(index);
    output_input.header_lines = input.header_lines;
    output_input.sample = input.sample;
    output_input.intervals.clear();
    for (const auto& interval : panel.panel_intervals)
        output_input.intervals.push_back({interval.contig, interval.start, interval.end, 0.0});
    return values;
}

std::vector<double> denoise_hdf5_panel(const std::vector<double>& target,
                                       const fastgatk::hdf5::PanelOfNormals& panel,
                                       int requested_rank,
                                       KernelStats* kernel_stats = nullptr) {
    if (panel.eigensample_vectors.empty() || target.empty()) return target;
    const std::size_t intervals = target.size();
    if (panel.eigensample_vectors.size() != intervals)
        throw std::invalid_argument("BAD_INPUT: HDF5 panel eigensample interval dimension does not match case");
    const std::size_t available = panel.eigensample_vectors.front().size();
    const std::size_t rank = requested_rank < 0
        ? available
        : std::min<std::size_t>(available, static_cast<std::size_t>(requested_rank));
    if (rank == 0) return target;
    Kokkos::View<double*> device_values("pon_case_values", intervals);
    Kokkos::View<double*> device_eigen("pon_eigen", intervals);
    Kokkos::View<double*> device_score("pon_score", 1);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t index = 0; index < intervals; ++index) host_values(index) = target[index];
    Kokkos::deep_copy(device_values, host_values);
    for (std::size_t component = 0; component < rank; ++component) {
        auto host_eigen = Kokkos::create_mirror_view(device_eigen);
        for (std::size_t index = 0; index < intervals; ++index) host_eigen(index) = panel.eigensample_vectors[index][component];
        fastgatk::core::HostBatch component_host("denoise-hdf5-pon-component-v1");
        component_host.records = intervals;
        component_host.bytes = 2 * intervals * sizeof(double) + sizeof(double);
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::KernelPlan<ExecSpace> component_plan("denoise-hdf5-pon-component");
        component_plan.begin_prepare(component_host);
        Kokkos::deep_copy(device_eigen, host_eigen);
        fastgatk::core::DeviceBatch<ExecSpace> component_device(intervals);
        component_device.bind("values", device_values);
        component_device.bind("eigen", device_eigen);
        component_device.bind("score", device_score);
        ExecSpace().fence();
        component_plan.end_prepare(component_device);
        component_plan.begin_execute();
        Kokkos::parallel_for("denoise_hdf5_pon_projection_score", Kokkos::RangePolicy<ExecSpace>(0, 1), KOKKOS_LAMBDA(const std::size_t) {
            double score = 0.0;
            for (std::size_t index = 0; index < intervals; ++index) score += device_values(index) * device_eigen(index);
            device_score(0) = score;
        });
        ExecSpace().fence();
        component_plan.end_execute();
        auto host_score = Kokkos::create_mirror_view(device_score);
        Kokkos::deep_copy(host_score, device_score);
        const double score = host_score(0);
        component_plan.begin_execute();
        Kokkos::parallel_for("denoise_hdf5_pon_projection_subtract", Kokkos::RangePolicy<ExecSpace>(0, intervals),
            KOKKOS_LAMBDA(const std::size_t index) { device_values(index) -= score * device_eigen(index); });
        ExecSpace().fence();
        component_plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += intervals;
            kernel_stats->prepare_seconds += component_plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += component_plan.telemetry().execute_seconds;
        }
    }
    Kokkos::deep_copy(host_values, device_values);
    std::vector<double> result(intervals);
    for (std::size_t index = 0; index < intervals; ++index) result[index] = host_values(index);
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto started = std::chrono::steady_clock::now();
        const auto options = parse(argc, argv);
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(static_cast<unsigned>(std::max(1, options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        const bool input_hdf5 = options.format == "HDF5" || options.input.ends_with(".h5") || options.input.ends_with(".hdf5");
        const auto input = read_counts_any(options.input, input_hdf5);
        std::vector<double> annotated_gc_content;
        if (!options.annotated_intervals.empty()) {
            const auto annotated = read_annotated_intervals(options.annotated_intervals);
            if (annotated.intervals.size() != input.intervals.size())
                throw std::invalid_argument("BAD_INPUT: annotated intervals do not match case interval count");
            for (std::size_t index = 0; index < input.intervals.size(); ++index) {
                if (annotated.intervals[index].contig != input.intervals[index].contig ||
                    annotated.intervals[index].start != input.intervals[index].start ||
                    annotated.intervals[index].end != input.intervals[index].end)
                    throw std::invalid_argument("BAD_INPUT: annotated intervals are not aligned with case intervals");
            }
            annotated_gc_content = annotated.gc_content;
        }
        PanelTable panel;
        fastgatk::hdf5::PanelOfNormals hdf5_panel;
        bool hdf5_panel_of_normals = false;
        std::vector<double> panel_counts;
        const bool has_panel = !options.panel_of_normals.empty();
        bool svd_used = false;
        int svd_rank = 0;
        bool gc_content_used = false;
        if (has_panel) {
            const bool panel_hdf5 = options.panel_of_normals.ends_with(".h5") || options.panel_of_normals.ends_with(".hdf5");
            if (panel_hdf5) {
                try {
                    hdf5_panel = fastgatk::hdf5::read_panel_of_normals(options.panel_of_normals);
                    hdf5_panel_of_normals = true;
                    panel.wide = false;
                    panel.samples = hdf5_panel.panel_sample_filenames;
                    for (const auto& interval : hdf5_panel.panel_intervals)
                        panel.intervals.push_back({interval.contig, interval.start, interval.end, 0.0});
                } catch (const std::exception&) {
                    const auto panel_collection = fastgatk::hdf5::read_simple_count_collection(options.panel_of_normals);
                    panel.wide = false;
                    panel.samples = {panel_collection.sample};
                    for (std::size_t index = 0; index < panel_collection.intervals.size(); ++index) {
                        panel.intervals.push_back({panel_collection.intervals[index].contig,
                                                   panel_collection.intervals[index].start,
                                                   panel_collection.intervals[index].end,
                                                   panel_collection.counts[index]});
                    }
                }
            } else {
                panel = read_panel(options.panel_of_normals);
            }
            if (!hdf5_panel_of_normals && panel.intervals.size() != input.intervals.size())
                throw std::invalid_argument("BAD_INPUT: panel-of-normals interval count does not match sample");
            if (!hdf5_panel_of_normals) {
                for (std::size_t index = 0; index < input.intervals.size(); ++index) {
                    const auto& lhs = input.intervals[index];
                    const auto& rhs = panel.intervals[index];
                    if (lhs.contig != rhs.contig || lhs.start != rhs.start || lhs.end != rhs.end)
                        throw std::invalid_argument("BAD_INPUT: panel-of-normals intervals are not aligned");
                    if (!panel.wide) panel_counts.push_back(rhs.count);
                }
            }
        }
        std::string sample = options.sample.empty() ? input.sample : options.sample;
        if (sample.empty()) sample = "UNKNOWN";
        std::vector<double> counts;
        counts.reserve(input.intervals.size());
        for (const auto& interval : input.intervals) counts.push_back(interval.count);
        CountTable output_input = input;
        KernelStats kernel_stats;
        std::vector<double> standardized_values;
        std::vector<double> denoised_values;
        if (hdf5_panel_of_normals) {
            standardized_values = standardize_with_hdf5_panel(input, hdf5_panel, output_input,
                                                              nullptr,
                                                              &kernel_stats);
            gc_content_used = !hdf5_panel.original_interval_gc_content.empty();
            denoised_values = denoise_hdf5_panel(standardized_values, hdf5_panel, options.svd_rank,
                                                 &kernel_stats);
            const auto available = hdf5_panel.eigensample_vectors.empty()
                ? 0 : static_cast<int>(hdf5_panel.eigensample_vectors.front().size());
            svd_rank = options.svd_rank < 0 ? available : std::min(options.svd_rank, available);
            svd_used = svd_rank > 0;
        } else if (has_panel && panel.wide && panel.samples.size() > 1) {
            std::vector<std::vector<double>> panel_values;
            panel_values.reserve(panel.samples.size());
            for (const auto& panel_sample : panel.counts)
                panel_values.push_back(standardize(panel_sample, nullptr, options.threads,
                                                   nullptr,
                                                   &kernel_stats));
            standardized_values = standardize(counts, nullptr, options.threads,
                                              nullptr,
                                              &kernel_stats);
            denoised_values = denoise_svd(standardized_values, panel_values,
                                          options.svd_rank, options.threads, &kernel_stats);
            svd_rank = options.svd_rank < 0
                ? std::min(20, static_cast<int>(std::min(panel.samples.size(), input.intervals.size())))
                : std::min(options.svd_rank,
                           static_cast<int>(std::min(panel.samples.size(), input.intervals.size())));
            svd_used = svd_rank > 0;
        } else {
            standardized_values = standardize(counts, has_panel ? &panel_counts : nullptr, options.threads,
                                              has_panel ? nullptr
                                                        : (annotated_gc_content.empty() ? nullptr : &annotated_gc_content),
                                              &kernel_stats);
            gc_content_used = !has_panel && !annotated_gc_content.empty();
            denoised_values = standardized_values;
        }
        const std::string standardized = options.standardized_output;
        if (!standardized.empty()) write_table(standardized, output_input, sample, standardized_values);
        if (!options.output.empty()) write_table(options.output, output_input, sample, denoised_values);
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
            manifest << "{\"schema_version\":1,\"tool\":\"DenoiseReadCounts\","
                     << "\"implementation\":\"fastgatk-denoise-read-counts\","
                     << "\"status\":\"prototype\",\"input\":\"" << json_escape(options.input)
                     << "\",\"output\":\"" << json_escape(options.output)
                     << "\",\"standardized_output\":\"" << json_escape(standardized)
                     << "\",\"format\":\"" << options.format << "\",\"sample\":\"" << json_escape(sample)
                     << "\",\"intervals\":" << output_input.intervals.size()
                     << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"determinism\":\"strict\""
                     << ",\"panel_of_normals\":" << (has_panel ? "true" : "false")
                     << ",\"panel_samples\":" << (has_panel ? panel.samples.size() : 0)
                     << ",\"gc_content\":" << (gc_content_used ? "true" : "false")
                     << ",\"hdf5\":" << ((input_hdf5 || hdf5_panel_of_normals) ? "true" : "false") << ",\"svd\":" << (svd_used ? "true" : "false")
                     << ",\"svd_rank\":" << svd_rank
                     << ",\"semantics\":\"fractional-coverage/median/log2/median\""
                     << ",\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"denoised-copy-ratios\",\"complete\":"
                     << (file_complete(options.output) ? "true" : "false") << "}"
                     << (standardized.empty() ? "" : ",{\"path\":\"" + json_escape(standardized)
                         + "\",\"kind\":\"standardized-copy-ratios\",\"complete\":"
                         + (file_complete(standardized) ? "true}" : "false}"))
                     << "]"
                     << ",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                     << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"kernel_execution_policy\":\"MDRangePolicy+RangePolicy\""
                     << ",\"kernel_batches\":" << kernel_stats.batches
                     << ",\"kernel_observations\":" << kernel_stats.observations
                     << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds
                     << ",\"output_bytes\":" << file_bytes(options.output)
                     << ",\"standardized_output_bytes\":" << file_bytes(standardized)
                     << ",\"wall_seconds\":" << wall_seconds << "}}\n";
            manifest.close();
            if (!manifest || !file_complete(options.manifest))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest: " + options.manifest);
        }
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"DenoiseReadCounts\",\"status\":\"prototype\","
                  << "\"format\":\"" << options.format << "\",\"intervals\":" << output_input.intervals.size()
                  << ",\"panel_of_normals\":" << (has_panel ? "true" : "false")
                  << ",\"panel_samples\":" << (has_panel ? panel.samples.size() : 0)
                  << ",\"gc_content\":" << (gc_content_used ? "true" : "false")
                  << ",\"svd\":" << (svd_used ? "true" : "false")
                  << ",\"svd_rank\":" << svd_rank
                  << ",\"lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                  << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                  << ",\"policy\":\"MDRangePolicy+RangePolicy\""
                  << ",\"kernel_batches\":" << kernel_stats.batches
                  << ",\"kernel_observations\":" << kernel_stats.observations
                  << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                  << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds << "}\n";
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << error.what() << '\n';
        return 2;
    }
}
