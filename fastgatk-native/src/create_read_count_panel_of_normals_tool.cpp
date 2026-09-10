#include <Kokkos_Core.hpp>
#include "fastgatk/core/plan.hpp"
#include "fastgatk/hdf5_count_collection.hpp"

#include <algorithm>
#include <cmath>
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

struct Options {
    std::vector<std::string> inputs;
    std::string output;
    std::string manifest;
    std::string annotated_intervals;
    int eigensamples = 20;
    double minimum_interval_median_percentile = 10.0;
    double maximum_zeros_in_sample_percentage = 5.0;
    double maximum_zeros_in_interval_percentage = 5.0;
    double extreme_sample_median_percentile = 2.5;
    bool do_impute_zeros = true;
    double extreme_outlier_truncation_percentile = 0.1;
    int maximum_chunk_size = 134217727;
    int threads = 1;
};

struct CountTable {
    std::string sample = "UNKNOWN";
    std::string sequence_dictionary;
    std::vector<fastgatk::hdf5::CountInterval> intervals;
    std::vector<double> counts;
};

bool same_interval(const fastgatk::hdf5::CountInterval& lhs,
                   const fastgatk::hdf5::CountInterval& rhs) {
    return lhs.contig == rhs.contig && lhs.start == rhs.start && lhs.end == rhs.end;
}

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv, const std::string& argument,
                          const char* name, const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

bool parse_bool(const std::string& value, const char* name) {
    if (value == "true" || value == "TRUE" || value == "1") return true;
    if (value == "false" || value == "FALSE" || value == "0") return false;
    throw std::invalid_argument(std::string("BAD_INPUT: ") + name + " expects true or false");
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-create-read-count-panel-of-normals (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                         input SimpleCountCollection TSV/HDF5 (repeatable)\n"
                         "  -O, --output FILE                        GATK HDF5 panel-of-normals output\n"
                         "      --number-of-eigensamples N           truncated SVD rank (default 20)\n"
                         "      --minimum-interval-median-percentile P\n"
                         "      --maximum-zeros-in-sample-percentage P\n"
                         "      --maximum-zeros-in-interval-percentage P\n"
                         "      --extreme-sample-median-percentile P\n"
                         "      --do-impute-zeros true|false\n"
                         "      --extreme-outlier-truncation-percentile P\n"
                         "      --maximum-chunk-size N               maximum matrix cells per HDF5 row chunk\n"
                         "      --annotated-intervals FILE           GC_CONTENT annotated intervals TSV\n"
                         "      --threads N                          Kokkos execution threads\n"
                         "      --output-manifest FILE\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.inputs.push_back(require_value(index, argc, argv, argument, "--input", "-I"));
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--number-of-eigensamples") || has_option(argument, "--svd-rank")) {
            const char* name = argument.rfind("--svd-rank", 0) == 0 ? "--svd-rank" : "--number-of-eigensamples";
            options.eigensamples = std::stoi(require_value(index, argc, argv, argument, name));
        } else if (has_option(argument, "--minimum-interval-median-percentile")) {
            options.minimum_interval_median_percentile = std::stod(require_value(index, argc, argv, argument, "--minimum-interval-median-percentile"));
        } else if (has_option(argument, "--maximum-zeros-in-sample-percentage")) {
            options.maximum_zeros_in_sample_percentage = std::stod(require_value(index, argc, argv, argument, "--maximum-zeros-in-sample-percentage"));
        } else if (has_option(argument, "--maximum-zeros-in-interval-percentage")) {
            options.maximum_zeros_in_interval_percentage = std::stod(require_value(index, argc, argv, argument, "--maximum-zeros-in-interval-percentage"));
        } else if (has_option(argument, "--extreme-sample-median-percentile")) {
            options.extreme_sample_median_percentile = std::stod(require_value(index, argc, argv, argument, "--extreme-sample-median-percentile"));
        } else if (has_option(argument, "--do-impute-zeros")) {
            options.do_impute_zeros = parse_bool(require_value(index, argc, argv, argument, "--do-impute-zeros"), "--do-impute-zeros");
        } else if (has_option(argument, "--extreme-outlier-truncation-percentile")) {
            options.extreme_outlier_truncation_percentile = std::stod(require_value(index, argc, argv, argument, "--extreme-outlier-truncation-percentile"));
        } else if (has_option(argument, "--maximum-chunk-size")) {
            options.maximum_chunk_size = std::stoi(require_value(index, argc, argv, argument, "--maximum-chunk-size"));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (has_option(argument, "--annotated-intervals") || has_option(argument, "--sequence-dictionary")) {
            const char* name = argument.rfind("--annotated-intervals", 0) == 0
                ? "--annotated-intervals" : "--sequence-dictionary";
            const auto value = require_value(index, argc, argv, argument, name);
            if (argument.rfind("--annotated-intervals", 0) == 0) options.annotated_intervals = value;
            else throw std::invalid_argument("BACKEND_UNAVAILABLE: explicit sequence dictionary requires the GATK fallback");
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates") || argument == "--quiet") {
            if (argument != "--quiet")
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" : "--seconds-between-progress-updates");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (!fastgatk::hdf5::available())
        throw std::invalid_argument("BACKEND_UNAVAILABLE: HDF5 support was not built; use --fallback");
    if (options.inputs.empty()) throw std::invalid_argument("-I/--input is required at least once");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.eigensamples < 0 || options.threads < 1 || options.maximum_chunk_size < 1)
        throw std::invalid_argument("BAD_INPUT: rank, threads and maximum chunk size must be positive or zero as applicable");
    const auto in_range = [](double value, double low, double high) {
        return std::isfinite(value) && value >= low && value <= high;
    };
    if (!in_range(options.minimum_interval_median_percentile, 0.0, 100.0) ||
        !in_range(options.maximum_zeros_in_sample_percentage, 0.0, 100.0) ||
        !in_range(options.maximum_zeros_in_interval_percentage, 0.0, 100.0) ||
        !in_range(options.extreme_sample_median_percentile, 0.0, 50.0) ||
        !in_range(options.extreme_outlier_truncation_percentile, 0.0, 50.0))
        throw std::invalid_argument("BAD_INPUT: filtering percentiles must be in their GATK ranges");
    return options;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

CountTable read_tsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open read-counts file: " + path);
    CountTable table;
    bool columns_seen = false;
    std::size_t line_number = 0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty()) continue;
        if (line[0] == '@') {
            if (line.rfind("@RG", 0) == 0) {
                for (const auto& field : split_tab(line))
                    if (field.rfind("SM:", 0) == 0) table.sample = field.substr(3);
            }
            if (line.rfind("@HD", 0) == 0 || line.rfind("@SQ", 0) == 0)
                table.sequence_dictionary += line + "\n";
            continue;
        }
        const auto fields = split_tab(line);
        if (!columns_seen) {
            if (fields.size() != 4 || fields[0] != "CONTIG" || fields[1] != "START" ||
                fields[2] != "END" || fields[3] != "COUNT")
                throw std::runtime_error("BAD_INPUT: expected SimpleCountCollection TSV header at " + path + ":" + std::to_string(line_number));
            columns_seen = true;
            continue;
        }
        if (fields.size() != 4) throw std::runtime_error("BAD_INPUT: invalid count row at " + path);
        const auto start = std::stoll(fields[1]);
        const auto end = std::stoll(fields[2]);
        const auto count = std::stod(fields[3]);
        if (fields[0].empty() || start < 1 || end < start || !std::isfinite(count) || count < 0.0 ||
            std::abs(count - std::round(count)) > 1.0e-8)
            throw std::runtime_error("BAD_INPUT: invalid count row at " + path + ":" + std::to_string(line_number));
        table.intervals.push_back({fields[0], start, end});
        table.counts.push_back(count);
    }
    if (!columns_seen || table.intervals.empty()) throw std::runtime_error("BAD_INPUT: read-counts file has no intervals: " + path);
    if (table.sample.empty()) table.sample = "UNKNOWN";
    return table;
}

CountTable read_any(const std::string& path) {
    const bool hdf5 = std::filesystem::path(path).extension() == ".h5" ||
                      std::filesystem::path(path).extension() == ".hdf5";
    if (!hdf5) return read_tsv(path);
    const auto collection = fastgatk::hdf5::read_simple_count_collection(path);
    CountTable table;
    table.sample = collection.sample;
    table.sequence_dictionary = collection.sequence_dictionary;
    table.intervals = collection.intervals;
    table.counts = collection.counts;
    if (table.intervals.empty() || table.intervals.size() != table.counts.size())
        throw std::runtime_error("BAD_INPUT: HDF5 read-counts file is empty or misaligned: " + path);
    return table;
}

struct AnnotatedIntervals {
    std::vector<fastgatk::hdf5::CountInterval> intervals;
    std::vector<double> gc_content;
};

AnnotatedIntervals read_annotated_intervals(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open annotated intervals: " + path);
    AnnotatedIntervals result;
    bool header_seen = false;
    std::size_t gc_column = 0;
    std::size_t line_number = 0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        if (line.empty() || line[0] == '@') continue;
        const auto fields = split_tab(line);
        if (!header_seen) {
            if (fields.size() < 4 || fields[0] != "CONTIG" || fields[1] != "START" ||
                fields[2] != "END")
                throw std::runtime_error("BAD_INPUT: annotated intervals must contain CONTIG/START/END/GC_CONTENT: " + path);
            while (gc_column < fields.size() && fields[gc_column] != "GC_CONTENT") ++gc_column;
            if (gc_column == fields.size()) throw std::runtime_error("BAD_INPUT: annotated intervals lack GC_CONTENT: " + path);
            result.gc_content.clear();
            result.intervals.clear();
            header_seen = true;
            continue;
        }
        if (fields.size() <= gc_column || fields.size() < 4)
            throw std::runtime_error("BAD_INPUT: malformed annotated interval at " + path + ":" + std::to_string(line_number));
        const auto start = std::stoll(fields[1]);
        const auto end = std::stoll(fields[2]);
        const auto gc = std::stod(fields[gc_column]);
        if (fields[0].empty() || start < 1 || end < start ||
            !(std::isnan(gc) || (std::isfinite(gc) && gc >= 0.0 && gc <= 1.0)))
            throw std::runtime_error("BAD_INPUT: invalid annotated interval at " + path + ":" + std::to_string(line_number));
        result.intervals.push_back({fields[0], start, end});
        result.gc_content.push_back(gc);
    }
    if (!header_seen || result.intervals.empty() || result.intervals.size() != result.gc_content.size())
        throw std::runtime_error("BAD_INPUT: annotated intervals contain no rows: " + path);
    return result;
}

double median(std::vector<double> values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 == 0 ? 0.5 * (values[middle - 1] + values[middle]) : values[middle];
}

double percentile(std::vector<double> values, double percent) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    if (values.size() == 1 || percent <= 0.0) return values.front();
    if (percent >= 100.0) return values.back();
    const double position = (percent / 100.0) * static_cast<double>(values.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(position));
    const auto upper = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return values[lower] + fraction * (values[upper] - values[lower]);
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const char value : text) {
        if (value == '"' || value == '\\') escaped << '\\';
        if (value == '\n') escaped << "\\n";
        else if (value == '\r') escaped << "\\r";
        else if (value == '\t') escaped << "\\t";
        else escaped << value;
    }
    return escaped.str();
}

struct PreparedPanel {
    std::vector<std::vector<double>> original;
    std::vector<std::string> sample_files;
    std::vector<fastgatk::hdf5::CountInterval> original_intervals;
    std::vector<double> original_interval_gc_content;
    std::string sequence_dictionary;
    std::vector<std::string> panel_sample_files;
    std::vector<fastgatk::hdf5::CountInterval> panel_intervals;
    std::vector<double> panel_interval_fractional_medians;
    std::vector<std::vector<double>> standardized;
    std::vector<double> singular_values;
    std::vector<std::vector<double>> eigensample_vectors;
    std::size_t filtered_samples = 0;
    std::size_t filtered_intervals = 0;
    std::size_t imputed_values = 0;
};

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

int gc_bin(double gc_content) {
    if (!std::isfinite(gc_content)) return 0;
    return std::clamp(static_cast<int>(std::llround(gc_content * 100.0)), 0, 100);
}

void correct_gc_bias(std::vector<std::vector<double>>& values,
                     const std::vector<double>& gc_content) {
    if (gc_content.empty()) return;
    if (values.empty() || values.front().size() != gc_content.size())
        throw std::invalid_argument("BAD_INPUT: GC annotations do not match read-count intervals");
    constexpr int bins = 101;
    constexpr double decay_rate = 1.0 / (0.02 * bins);
    for (auto& row : values) {
        std::vector<std::vector<double>> by_gc(bins);
        const double total_before = std::accumulate(row.begin(), row.end(), 0.0);
        for (std::size_t interval = 0; interval < row.size(); ++interval)
            by_gc[gc_bin(gc_content[interval])].push_back(row[interval]);
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
        for (std::size_t interval = 0; interval < row.size(); ++interval) {
            row[interval] *= factors[gc_bin(gc_content[interval])];
            total_after += row[interval];
        }
        if (total_after > 0.0 && total_before > 0.0) {
            const double normalization = total_before / total_after;
            for (auto& value : row) value *= normalization;
        }
    }
}

PreparedPanel preprocess(const std::vector<CountTable>& tables, const Options& options,
                         const std::vector<double>& interval_gc_content,
                         const std::vector<std::string>& input_sample_files) {
    PreparedPanel result;
    const std::size_t sample_count = tables.size();
    if (input_sample_files.size() != sample_count)
        throw std::invalid_argument("BAD_INPUT: sample-file metadata does not match count tables");
    const std::size_t interval_count = tables.front().intervals.size();
    result.original_intervals = tables.front().intervals;
    result.original_interval_gc_content = interval_gc_content;
    result.sequence_dictionary = tables.front().sequence_dictionary;
    result.sample_files.reserve(sample_count);
    result.original.resize(sample_count, std::vector<double>(interval_count));
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        // GATK's HDF5SVDReadCountPanelOfNormals stores the input count-file
        // paths (File#getAbsolutePath), not the @RG SM/sample identifier.
        // Keeping this metadata from the CLI inputs also ensures that the
        // filtered / panel sample list is aligned with the original rows.
        result.sample_files.push_back(input_sample_files[sample]);
        result.original[sample] = tables[sample].counts;
    }
    std::vector<std::vector<double>> values = result.original;
    std::vector<bool> filter_samples(sample_count, false), filter_intervals(interval_count, false);
    std::vector<double> interval_medians(interval_count);
    for (std::size_t interval = 0; interval < interval_count; ++interval) {
        std::vector<double> column(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample) column[sample] = values[sample][interval];
        double total = 0.0;
        for (const auto count : column) total += count;
        if (!(total > 0.0)) throw std::runtime_error("BAD_INPUT: panel sample has no positive counts");
        interval_medians[interval] = median(column); // replaced below after fractional conversion
    }
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const double total = std::accumulate(values[sample].begin(), values[sample].end(), 0.0);
        if (!(total > 0.0)) throw std::runtime_error("BAD_INPUT: sample has no positive read counts");
        for (auto& value : values[sample]) value /= total;
    }
    correct_gc_bias(values, interval_gc_content);
    for (std::size_t interval = 0; interval < interval_count; ++interval) {
        std::vector<double> column(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample) column[sample] = values[sample][interval];
        interval_medians[interval] = median(column);
    }
    if (options.minimum_interval_median_percentile > 0.0) {
        const double threshold = percentile(interval_medians, options.minimum_interval_median_percentile);
        for (std::size_t interval = 0; interval < interval_count; ++interval)
            if (interval_medians[interval] <= threshold) filter_intervals[interval] = true;
    }
    for (std::size_t interval = 0; interval < interval_count; ++interval) {
        if (filter_intervals[interval]) continue;
        if (!(interval_medians[interval] > 0.0)) throw std::runtime_error("BAD_INPUT: retained interval has zero fractional median");
        for (std::size_t sample = 0; sample < sample_count; ++sample) values[sample][interval] /= interval_medians[interval];
    }
    auto passing_intervals = [&]() {
        return static_cast<std::size_t>(std::count(filter_intervals.begin(), filter_intervals.end(), false));
    };
    if (passing_intervals() == 0) throw std::runtime_error("BAD_INPUT: interval filtering removed all intervals");
    if (options.maximum_zeros_in_sample_percentage < 100.0) {
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            std::size_t zeros = 0;
            for (std::size_t interval = 0; interval < interval_count; ++interval)
                if (!filter_intervals[interval] && values[sample][interval] == 0.0) ++zeros;
            if (static_cast<double>(zeros) / static_cast<double>(passing_intervals()) >= options.maximum_zeros_in_sample_percentage / 100.0)
                filter_samples[sample] = true;
        }
    }
    auto passing_samples = [&]() {
        return static_cast<std::size_t>(std::count(filter_samples.begin(), filter_samples.end(), false));
    };
    if (passing_samples() == 0) throw std::runtime_error("BAD_INPUT: sample filtering removed all samples");
    if (options.maximum_zeros_in_interval_percentage < 100.0) {
        for (std::size_t interval = 0; interval < interval_count; ++interval) {
            if (filter_intervals[interval]) continue;
            std::size_t zeros = 0;
            for (std::size_t sample = 0; sample < sample_count; ++sample)
                if (!filter_samples[sample] && values[sample][interval] == 0.0) ++zeros;
            if (static_cast<double>(zeros) / static_cast<double>(passing_samples()) >= options.maximum_zeros_in_interval_percentage / 100.0)
                filter_intervals[interval] = true;
        }
    }
    if (passing_intervals() == 0) throw std::runtime_error("BAD_INPUT: interval filtering removed all intervals");
    if (options.extreme_sample_median_percentile > 0.0) {
        std::vector<double> sample_medians(sample_count, 0.0);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            std::vector<double> row;
            for (std::size_t interval = 0; interval < interval_count; ++interval)
                if (!filter_intervals[interval]) row.push_back(values[sample][interval]);
            sample_medians[sample] = median(row);
        }
        const double low = percentile(sample_medians, options.extreme_sample_median_percentile);
        const double high = percentile(sample_medians, 100.0 - options.extreme_sample_median_percentile);
        for (std::size_t sample = 0; sample < sample_count; ++sample)
            if (sample_medians[sample] < low || sample_medians[sample] > high) filter_samples[sample] = true;
    }
    if (passing_samples() == 0 || passing_intervals() == 0) throw std::runtime_error("BAD_INPUT: filtering removed all panel rows or columns");
    std::vector<std::size_t> kept_samples, kept_intervals;
    for (std::size_t sample = 0; sample < sample_count; ++sample) if (!filter_samples[sample]) kept_samples.push_back(sample);
    for (std::size_t interval = 0; interval < interval_count; ++interval) if (!filter_intervals[interval]) kept_intervals.push_back(interval);
    result.filtered_samples = sample_count - kept_samples.size();
    result.filtered_intervals = interval_count - kept_intervals.size();
    result.panel_sample_files.reserve(kept_samples.size());
    for (const auto sample : kept_samples) result.panel_sample_files.push_back(result.sample_files[sample]);
    result.panel_intervals.reserve(kept_intervals.size());
    for (const auto interval : kept_intervals) result.panel_intervals.push_back(result.original_intervals[interval]);
    for (const auto interval : kept_intervals) result.panel_interval_fractional_medians.push_back(interval_medians[interval]);
    std::vector<std::vector<double>> panel(kept_samples.size(), std::vector<double>(kept_intervals.size()));
    for (std::size_t row = 0; row < kept_samples.size(); ++row)
        for (std::size_t column = 0; column < kept_intervals.size(); ++column)
            panel[row][column] = values[kept_samples[row]][kept_intervals[column]];
    if (options.do_impute_zeros) {
        for (std::size_t column = 0; column < kept_intervals.size(); ++column) {
            std::vector<double> nonzero;
            for (const auto& row : panel) if (row[column] > 0.0) nonzero.push_back(row[column]);
            if (nonzero.empty()) throw std::runtime_error("BAD_INPUT: cannot impute an all-zero retained interval");
            const double replacement = median(nonzero);
            for (auto& row : panel) if (row[column] == 0.0) { row[column] = replacement; ++result.imputed_values; }
        }
    }
    if (options.extreme_outlier_truncation_percentile > 0.0) {
        std::vector<double> all;
        for (const auto& row : panel) all.insert(all.end(), row.begin(), row.end());
        const double low = percentile(all, options.extreme_outlier_truncation_percentile);
        const double high = percentile(all, 100.0 - options.extreme_outlier_truncation_percentile);
        for (auto& row : panel) for (auto& value : row) value = std::clamp(value, low, high);
    }
    std::vector<double> row_medians(kept_samples.size());
    for (std::size_t row = 0; row < panel.size(); ++row) {
        row_medians[row] = median(panel[row]);
        if (!(row_medians[row] > 0.0)) throw std::runtime_error("BAD_INPUT: panel sample has non-positive median");
        for (auto& value : panel[row]) {
            const double ratio = value / row_medians[row];
            value = ratio < EPSILON ? std::log(EPSILON) * INV_LOG_2 : std::log(ratio) * INV_LOG_2;
        }
        row_medians[row] = median(panel[row]);
    }
    const double median_of_medians = median(row_medians);
    // Persist the actual post-log2 matrix.  Copying `panel` before the
    // per-sample log2 transform would leave the SVD operating on fractional
    // coverage ratios instead of GATK's standardized values, and could make
    // an identical-normal panel appear to contain a non-zero component.
    result.standardized = panel;
    for (auto& row : result.standardized) for (auto& value : row) value -= median_of_medians;
    return result;
}

void compute_svd(PreparedPanel& panel, int requested_rank, KernelStats* kernel_stats = nullptr) {
    const std::size_t samples = panel.standardized.size();
    const std::size_t intervals = panel.standardized.front().size();
    const std::size_t rank = std::min<std::size_t>({static_cast<std::size_t>(std::max(0, requested_rank)), samples, intervals});
    if (samples <= 1 || rank == 0) return;
    Kokkos::View<double**> matrix("pon_standardized", samples, intervals);
    Kokkos::View<double**> gram("pon_gram", samples, samples);
    auto host_matrix = Kokkos::create_mirror_view(matrix);
    for (std::size_t sample = 0; sample < samples; ++sample)
        for (std::size_t interval = 0; interval < intervals; ++interval) host_matrix(sample, interval) = panel.standardized[sample][interval];
    Kokkos::deep_copy(matrix, host_matrix);
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch gram_host("create-read-count-pon-gram-v1");
    gram_host.records = samples * samples;
    gram_host.bytes = matrix.span() * sizeof(double) + gram.span() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> gram_plan("create-pon-gram");
    gram_plan.begin_prepare(gram_host);
    fastgatk::core::DeviceBatch<ExecSpace> gram_device(samples * samples);
    gram_device.bind("matrix", matrix);
    gram_device.bind("gram", gram);
    ExecSpace().fence();
    gram_plan.end_prepare(gram_device);
    gram_plan.begin_execute();
    Kokkos::parallel_for("create_pon_gram", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0, 0}, {samples, samples}),
                         KOKKOS_LAMBDA(const std::size_t row, const std::size_t column) {
        double sum = 0.0;
        for (std::size_t interval = 0; interval < intervals; ++interval) sum += matrix(row, interval) * matrix(column, interval);
        gram(row, column) = sum;
    });
    ExecSpace().fence();
    gram_plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += samples * samples;
        kernel_stats->prepare_seconds += gram_plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += gram_plan.telemetry().execute_seconds;
    }
    auto host_gram = Kokkos::create_mirror_view(gram);
    Kokkos::deep_copy(host_gram, gram);
    std::vector<double> covariance(samples * samples);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t column = 0; column < samples; ++column) covariance[row * samples + column] = host_gram(row, column);
    std::vector<std::vector<double>> right_vectors;
    for (std::size_t component = 0; component < rank; ++component) {
        std::vector<double> vector(samples), next(samples);
        for (std::size_t index = 0; index < samples; ++index)
            vector[index] = 1.0 + 0.1732050807568877 * static_cast<double>((index + 1) * (component + 1));
        for (const auto& prior : right_vectors) {
            double projection = std::inner_product(vector.begin(), vector.end(), prior.begin(), 0.0);
            for (std::size_t index = 0; index < samples; ++index) vector[index] -= projection * prior[index];
        }
        double norm = std::sqrt(std::inner_product(vector.begin(), vector.end(), vector.begin(), 0.0));
        if (!(norm > EPSILON)) break;
        for (auto& value : vector) value /= norm;
        for (int iteration = 0; iteration < 96; ++iteration) {
            std::fill(next.begin(), next.end(), 0.0);
            for (std::size_t row = 0; row < samples; ++row)
                for (std::size_t column = 0; column < samples; ++column) next[row] += covariance[row * samples + column] * vector[column];
            for (const auto& prior : right_vectors) {
                const double projection = std::inner_product(next.begin(), next.end(), prior.begin(), 0.0);
                for (std::size_t index = 0; index < samples; ++index) next[index] -= projection * prior[index];
            }
            norm = std::sqrt(std::inner_product(next.begin(), next.end(), next.begin(), 0.0));
            if (!(norm > EPSILON)) break;
            for (std::size_t index = 0; index < samples; ++index) vector[index] = next[index] / norm;
        }
        double eigenvalue = 0.0;
        for (std::size_t row = 0; row < samples; ++row)
            for (std::size_t column = 0; column < samples; ++column) eigenvalue += vector[row] * covariance[row * samples + column] * vector[column];
        if (!(eigenvalue > EPSILON) || !std::isfinite(eigenvalue)) break;
        right_vectors.push_back(vector);
        panel.singular_values.push_back(std::sqrt(eigenvalue));
        for (std::size_t row = 0; row < samples; ++row)
            for (std::size_t column = 0; column < samples; ++column) covariance[row * samples + column] -= eigenvalue * vector[row] * vector[column];
    }
    panel.eigensample_vectors.assign(intervals, std::vector<double>(right_vectors.size(), 0.0));
    for (std::size_t component = 0; component < right_vectors.size(); ++component) {
        const double singular = panel.singular_values[component];
        // Use a View-backed left projection so the dense projection remains a Kokkos kernel.
        Kokkos::View<double*> device_right("pon_right", samples);
        Kokkos::View<double*> device_left("pon_left", intervals);
        auto host_right = Kokkos::create_mirror_view(device_right);
        for (std::size_t sample = 0; sample < samples; ++sample) host_right(sample) = right_vectors[component][sample];
        Kokkos::deep_copy(device_right, host_right);
        fastgatk::core::HostBatch left_host("create-read-count-pon-left-v1");
        left_host.records = intervals;
        left_host.bytes = (samples + intervals) * sizeof(double);
        fastgatk::core::KernelPlan<ExecSpace> left_plan("create-pon-left-singular-vector");
        left_plan.begin_prepare(left_host);
        fastgatk::core::DeviceBatch<ExecSpace> left_device(intervals);
        left_device.bind("matrix", matrix);
        left_device.bind("right", device_right);
        left_device.bind("left", device_left);
        ExecSpace().fence();
        left_plan.end_prepare(left_device);
        left_plan.begin_execute();
        Kokkos::parallel_for("create_pon_left_singular_vector", Kokkos::RangePolicy<ExecSpace>(0, intervals), KOKKOS_LAMBDA(const std::size_t interval) {
            double sum = 0.0;
            for (std::size_t sample = 0; sample < samples; ++sample) sum += matrix(sample, interval) * device_right(sample);
            device_left(interval) = sum / singular;
        });
        ExecSpace().fence();
        left_plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += intervals;
            kernel_stats->prepare_seconds += left_plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += left_plan.telemetry().execute_seconds;
        }
        auto host_left = Kokkos::create_mirror_view(device_left);
        Kokkos::deep_copy(host_left, device_left);
        for (std::size_t interval = 0; interval < intervals; ++interval) panel.eigensample_vectors[interval][component] = host_left(interval);
    }
}

std::string command_line(int argc, char** argv) {
    std::ostringstream command;
    for (int index = 0; index < argc; ++index) {
        if (index) command << ' ';
        command << argv[index];
    }
    return command.str();
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto options = parse(argc, argv);
        std::vector<CountTable> tables;
        tables.reserve(options.inputs.size());
        for (const auto& path : options.inputs) tables.push_back(read_any(path));
        const auto& first = tables.front();
        for (std::size_t sample = 1; sample < tables.size(); ++sample) {
            if (tables[sample].intervals.size() != first.intervals.size() ||
                !std::equal(tables[sample].intervals.begin(), tables[sample].intervals.end(),
                            first.intervals.begin(), same_interval))
                throw std::invalid_argument("BAD_INPUT: all read-count files must have identical intervals in identical order");
        }
        std::vector<double> interval_gc_content;
        if (!options.annotated_intervals.empty()) {
            const auto annotated = read_annotated_intervals(options.annotated_intervals);
            if (annotated.intervals.size() != first.intervals.size() ||
                !std::equal(annotated.intervals.begin(), annotated.intervals.end(), first.intervals.begin(), same_interval))
                throw std::invalid_argument("BAD_INPUT: annotated intervals must match read-count intervals in identical order");
            interval_gc_content = annotated.gc_content;
        }
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(static_cast<unsigned>(std::max(1, options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        std::vector<std::string> input_sample_files;
        input_sample_files.reserve(options.inputs.size());
        for (const auto& input : options.inputs)
            input_sample_files.push_back(std::filesystem::absolute(input).string());
        auto panel = preprocess(tables, options, interval_gc_content, input_sample_files);
        KernelStats kernel_stats;
        compute_svd(panel, options.eigensamples, &kernel_stats);
        // GATK's HDF5SVDReadCountPanelOfNormals rejects a multi-sample panel
        // when SVD was requested but every singular value is numerically zero
        // (the Java path checks for at least one value > EPSILON before
        // writing the HDF5 container).  Do this check before creating the
        // output so a degenerate panel cannot look like a valid PoN and later
        // silently disable denoising.
        if (panel.standardized.size() > 1 && options.eigensamples > 0 &&
            panel.singular_values.empty()) {
            throw std::runtime_error(
                "java.lang.IllegalArgumentException: No non-zero singular values were found. "
                "It may be necessary to use stricter parameters for filtering.");
        }
        fastgatk::hdf5::PanelOfNormals output;
        output.maximum_chunk_size = static_cast<std::size_t>(options.maximum_chunk_size);
        output.version = 7.0;
        output.command_line = command_line(argc, argv);
        output.sequence_dictionary = panel.sequence_dictionary;
        output.original_read_counts = std::move(panel.original);
        output.original_sample_filenames = std::move(panel.sample_files);
        output.original_intervals = std::move(panel.original_intervals);
        output.original_interval_gc_content = std::move(panel.original_interval_gc_content);
        output.panel_sample_filenames = std::move(panel.panel_sample_files);
        output.panel_intervals = std::move(panel.panel_intervals);
        output.panel_interval_fractional_medians = std::move(panel.panel_interval_fractional_medians);
        output.singular_values = std::move(panel.singular_values);
        output.eigensample_vectors = std::move(panel.eigensample_vectors);
        fastgatk::hdf5::write_panel_of_normals(options.output, output);
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"CreateReadCountPanelOfNormals\",\"implementation\":\"fastgatk-create-read-count-panel-of-normals\",\"status\":\"prototype\",\"output\":\""
                     << json_escape(options.output) << "\",\"input_samples\":" << tables.size()
                     << ",\"original_intervals\":" << first.intervals.size()
                     << ",\"panel_samples\":" << output.panel_sample_filenames.size()
                     << ",\"panel_intervals\":" << output.panel_intervals.size()
                     << ",\"filtered_samples\":" << panel.filtered_samples
                     << ",\"filtered_intervals\":" << panel.filtered_intervals
                     << ",\"imputed_values\":" << panel.imputed_values
                     << ",\"gc_content\":" << (!output.original_interval_gc_content.empty() ? "true" : "false")
                     << ",\"eigensamples\":" << (output.eigensample_vectors.empty() ? 0 : output.eigensample_vectors.front().size())
                     << ",\"format\":\"HDF5-SVD-ReadCountPanelOfNormals-v7\",\"execution_space\":\""
                     << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                     << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"kernel_execution_policy\":\"MDRangePolicy+RangePolicy\""
                     << ",\"kernel_batches\":" << kernel_stats.batches
                     << ",\"kernel_observations\":" << kernel_stats.observations
                     << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds << "}}\n";
        }
        const auto eigensamples = output.eigensample_vectors.empty() ? 0 : output.eigensample_vectors.front().size();
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"CreateReadCountPanelOfNormals\",\"status\":\"prototype\",\"format\":\"HDF5-SVD-ReadCountPanelOfNormals-v7\",\"input_samples\":"
                  << tables.size() << ",\"original_intervals\":" << first.intervals.size()
                  << ",\"panel_samples\":" << output.panel_sample_filenames.size()
                  << ",\"panel_intervals\":" << output.panel_intervals.size()
                  << ",\"gc_content\":" << (!output.original_interval_gc_content.empty() ? "true" : "false")
                  << ",\"eigensamples\":" << eigensamples
                  << ",\"kernel_batches\":" << kernel_stats.batches
                  << ",\"kernel_observations\":" << kernel_stats.observations
                  << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                  << "\"}\n";
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << error.what() << '\n';
        return 2;
    }
}
