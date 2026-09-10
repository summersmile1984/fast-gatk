#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string output;
    std::string manifest;
    std::string tmp_dir;
    std::string verbosity = "INFO";
    bool quiet = false;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if (argument == name && index + 1 < argc) return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto ch : text) {
        if (ch == '"' || ch == '\\') escaped << '\\';
        if (ch == '\n') escaped << "\\n";
        else if (ch == '\r') escaped << "\\r";
        else if (ch == '\t') escaped << "\\t";
        else escaped << ch;
    }
    return escaped.str();
}

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-gather-bqsr-reports\n"
                         "  -I, --input REPORT             input FASTGATK-BQSR report (repeatable)\n"
                         "  -O, --output REPORT            merged report\n"
                         "      --QUIET[=true|false]       Picard/GATK quiet switch\n"
                         "      --tmp-dir DIR              temporary directory\n"
                         "      --verbosity LEVEL          ERROR/WARNING/INFO/DEBUG\n"
                         "      --output-manifest FILE     manifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || argument == "--input" || !option_value(argument, "--input").empty())
            options.inputs.push_back(argument == "-I" ?
                (i + 1 < argc ? argv[++i] : throw std::invalid_argument("missing value for -I")) :
                require_value(i, argc, argv, argument, "--input"));
        else if (argument == "-O" || argument == "--output" || !option_value(argument, "--output").empty())
            options.output = argument == "-O" ?
                (i + 1 < argc ? argv[++i] : throw std::invalid_argument("missing value for -O")) :
                require_value(i, argc, argv, argument, "--output");
        else if (argument == "--output-manifest" || !option_value(argument, "--output-manifest").empty())
            options.manifest = require_value(i, argc, argv, argument, "--output-manifest");
        else if (argument == "--QUIET" || argument == "--quiet" ||
                 argument.rfind("--QUIET=", 0) == 0 || argument.rfind("--quiet=", 0) == 0)
            options.quiet = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, argument.rfind("--QUIET", 0) == 0 ? "--QUIET" : "--quiet");
        else if (argument == "--tmp-dir" || !option_value(argument, "--tmp-dir").empty())
            options.tmp_dir = require_value(i, argc, argv, argument, "--tmp-dir");
        else if (argument == "--verbosity" || !option_value(argument, "--verbosity").empty()) {
            options.verbosity = require_value(i, argc, argv, argument, "--verbosity");
            std::transform(options.verbosity.begin(), options.verbosity.end(), options.verbosity.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
            if (options.verbosity != "ERROR" && options.verbosity != "WARNING" &&
                options.verbosity != "INFO" && options.verbosity != "DEBUG")
                throw std::invalid_argument("--verbosity must be ERROR, WARNING, INFO, or DEBUG");
        } else if (argument == "--use-jdk-deflater" || argument == "--jdk-deflater" ||
                   argument == "--use-jdk-inflater" || argument == "--jdk-inflater" ||
                   argument.rfind("--use-jdk-deflater=", 0) == 0 || argument.rfind("--jdk-deflater=", 0) == 0 ||
                   argument.rfind("--use-jdk-inflater=", 0) == 0 || argument.rfind("--jdk-inflater=", 0) == 0) {
            const bool deflater = argument.rfind("--use-jdk-deflater", 0) == 0 ||
                                  argument.rfind("--jdk-deflater", 0) == 0;
            (void)fastgatk::native::parse_optional_boolean(i, argc, argv, argument,
                deflater ? (argument.rfind("--jdk-", 0) == 0 ? "--jdk-deflater" : "--use-jdk-deflater")
                         : (argument.rfind("--jdk-", 0) == 0 ? "--jdk-inflater" : "--use-jdk-inflater"));
        }
        else
            throw std::invalid_argument("unknown option: " + argument);
    }
    if (options.inputs.empty()) throw std::invalid_argument("at least one --input is required");
    if (options.output.empty()) throw std::invalid_argument("--output is required");
    return options;
}

struct Bin {
    std::uint64_t count = 0;
    std::uint64_t mismatches = 0;
};

struct CovariateKey {
    std::string read_group;
    int cycle = 0;
    std::string context;
    int quality = 0;
    char event = 'M';
    // A zero value denotes the native sidecar representation, where one row
    // carries both the context and cycle covariates and the report writer
    // expands it into the two GATK RecalTable2 rows.  Java GATK reports carry
    // one row per covariate and use C (Context) or Y (Cycle) here.
    char covariate = '\0';
    bool operator<(const CovariateKey& other) const {
        if (read_group != other.read_group) return read_group < other.read_group;
        if (quality != other.quality) return quality < other.quality;
        if (covariate != other.covariate) return covariate < other.covariate;
        if (event != other.event) return event < other.event;
        if (cycle != other.cycle) return cycle < other.cycle;
        if (context != other.context) return context < other.context;
        return false;
    }
};

using CovariateTable = std::map<CovariateKey, Bin>;

struct ReportedDatum {
    Bin bin;
    double reported_quality = 0.0;
    bool initialized = false;
};

using ReportedDatumTable = std::map<std::pair<std::string, char>, ReportedDatum>;

double empirical_quality(const Bin& bin, int default_quality) {
    if (bin.count == 0) return static_cast<double>(default_quality);
    return -10.0 * std::log10(std::max(
        static_cast<double>(bin.mismatches) / static_cast<double>(bin.count), 1.0e-6));
}

int bayesian_empirical_quality(const Bin& bin, double prior_quality) {
    const auto observations = static_cast<long double>(bin.count);
    const auto errors = static_cast<long double>(bin.mismatches);
    if (observations == 0.0L)
        return std::clamp(static_cast<int>(std::lround(prior_quality)), 0, 40);
    const auto rounded_errors = static_cast<long double>(std::llround(errors)) + 1.0L;
    const auto smoothed_observations = observations + 2.0L;
    const auto prior = std::clamp(prior_quality, 0.0, 40.0);
    double best_score = -std::numeric_limits<double>::infinity();
    int best_quality = 0;
    for (int quality = 0; quality <= 40; ++quality) {
        const auto error_probability = std::pow(10.0L, -static_cast<long double>(quality) / 10.0L);
        if (!(error_probability > 0.0L && error_probability < 1.0L)) continue;
        // RecalDatum.getLogPrior casts before taking abs; retain its first
        // maximum tie-break for deterministic Java compatibility.
        const auto difference = std::min(std::abs(static_cast<int>(quality - prior)), 40);
        const auto log_prior = -0.5L * static_cast<long double>(difference * difference) / 0.25L;
        const auto log_binomial = std::lgammal(smoothed_observations + 1.0L) -
            std::lgammal(rounded_errors + 1.0L) -
            std::lgammal(smoothed_observations - rounded_errors + 1.0L) +
            rounded_errors * std::log(error_probability) +
            (smoothed_observations - rounded_errors) * std::log1pl(-error_probability);
        const auto score = static_cast<double>(log_prior + log_binomial);
        if (score > best_score) { best_score = score; best_quality = quality; }
    }
    return best_quality;
}

struct QuantizedResult {
    std::array<std::uint64_t, 94> counts{};
    std::array<int, 94> values{};
};

struct QuantizerLeaf {
    int quality = 0;
    std::uint64_t observations = 0;
};

struct QuantizerInterval {
    int begin = 0;
    int end = 0;
    std::uint64_t observations = 0;
    std::uint64_t errors = 0;
    int fixed_quality = -1;
    std::vector<QuantizerLeaf> leaves;
};

double quantizer_error_rate(const QuantizerInterval& interval) {
    if (interval.fixed_quality >= 0)
        return std::pow(10.0, -static_cast<double>(interval.fixed_quality) / 10.0);
    if (interval.observations == 0) return 0.0;
    return (static_cast<double>(interval.errors) + 1.0) /
           (static_cast<double>(interval.observations) + 1.0);
}

double quantizer_penalty(const QuantizerInterval& interval, int min_interesting_quality) {
    const auto global_error = quantizer_error_rate(interval);
    if (!(global_error > 0.0)) return 0.0;
    double penalty = 0.0;
    for (const auto& leaf : interval.leaves) {
        if (leaf.quality <= min_interesting_quality || leaf.observations == 0) continue;
        const auto leaf_error = std::pow(10.0, -static_cast<double>(leaf.quality) / 10.0);
        penalty += std::abs(std::log10(leaf_error) - std::log10(global_error)) *
                   static_cast<double>(leaf.observations);
    }
    return penalty;
}

QuantizedResult quantize_gatk(const std::array<Bin, 94>& reported_bins,
                              int quantization_levels = 16) {
    QuantizedResult result;
    std::vector<QuantizerInterval> intervals;
    intervals.reserve(94);
    // QuantizationInfo builds its histogram from rounded empirical qualities,
    // not raw reported qualities.  The RecalDatum Bayesian estimate above is
    // the same integer-valued statistic used by GATK's table writer.
    for (std::size_t quality = 0; quality < reported_bins.size(); ++quality) {
        const auto& bin = reported_bins[quality];
        const auto empirical = bayesian_empirical_quality(bin, static_cast<double>(quality));
        if (bin.count != 0) result.counts[static_cast<std::size_t>(empirical)] += bin.count;
    }
    for (int quality = 0; quality < 94; ++quality) {
        QuantizerInterval interval;
        interval.begin = interval.end = quality;
        interval.observations = result.counts[static_cast<std::size_t>(quality)];
        interval.errors = static_cast<std::uint64_t>(std::floor(
            static_cast<double>(interval.observations) *
            std::pow(10.0, -static_cast<double>(quality) / 10.0)));
        interval.fixed_quality = quality;
        interval.leaves.push_back({quality, interval.observations});
        intervals.push_back(std::move(interval));
    }
    quantization_levels = std::max(1, std::min(quantization_levels, 94));
    while (static_cast<int>(intervals.size()) > quantization_levels) {
        std::size_t best_index = 0;
        double best_penalty = std::numeric_limits<double>::infinity();
        for (std::size_t index = 0; index + 1 < intervals.size(); ++index) {
            QuantizerInterval merged;
            merged.begin = intervals[index].begin;
            merged.end = intervals[index + 1].end;
            merged.observations = intervals[index].observations + intervals[index + 1].observations;
            merged.errors = intervals[index].errors + intervals[index + 1].errors;
            merged.leaves = intervals[index].leaves;
            merged.leaves.insert(merged.leaves.end(), intervals[index + 1].leaves.begin(),
                                 intervals[index + 1].leaves.end());
            const auto penalty = quantizer_penalty(merged, 6);
            if (penalty < best_penalty) {
                best_penalty = penalty;
                best_index = index;
            }
        }
        QuantizerInterval merged;
        merged.begin = intervals[best_index].begin;
        merged.end = intervals[best_index + 1].end;
        merged.observations = intervals[best_index].observations + intervals[best_index + 1].observations;
        merged.errors = intervals[best_index].errors + intervals[best_index + 1].errors;
        merged.leaves = intervals[best_index].leaves;
        merged.leaves.insert(merged.leaves.end(), intervals[best_index + 1].leaves.begin(),
                             intervals[best_index + 1].leaves.end());
        intervals[best_index] = std::move(merged);
        intervals.erase(intervals.begin() + static_cast<std::ptrdiff_t>(best_index + 1));
    }
    result.values.fill(0);
    for (const auto& interval : intervals) {
        int value = interval.fixed_quality;
        if (value < 0) {
            const auto error = quantizer_error_rate(interval);
            value = error == 0.0 ? 93 : std::clamp(static_cast<int>(std::lround(
                -10.0 * std::log10(error))), 0, 93);
        }
        for (int quality = interval.begin; quality <= interval.end; ++quality)
            result.values[static_cast<std::size_t>(quality)] = value;
    }
    return result;
}

std::string display_read_group(const std::string& read_group) {
    return read_group.empty() ? "UNKNOWN" : read_group;
}

struct ReadGroupSummary { std::array<Bin, 94> bins{}; };

using ReadGroupTable = std::map<std::string, ReadGroupSummary>;

struct GatherInput {
    std::array<Bin, 94> bins{};
    CovariateTable covariates;
    ReadGroupTable groups;
    ReportedDatumTable reported;
    // GATK RecalibrationReport::isEmpty() is based on the materialized
    // RecalTable0/1/2 data, not merely on the presence of an Arguments table.
    // Keep this explicit so an all-empty BaseRecalibrator report fails at the
    // same gather boundary instead of producing a misleading zero-row report.
    bool has_recalibration_rows = false;
    // Keep the first report's argument rows in their source order.  GATK's
    // GatherBQSRReports preserves this table (and uses its quantizing level)
    // instead of silently reverting to the BaseRecalibrator defaults.
    std::vector<std::pair<std::string, std::string>> arguments;
};

const std::vector<std::pair<std::string, std::string>>& default_arguments() {
    static const std::vector<std::pair<std::string, std::string>> defaults = {
        {"binary_tag_name", "null"},
        {"covariate", "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate"},
        {"default_platform", "null"},
        {"deletions_default_quality", "45"},
        {"force_platform", "null"},
        {"indels_context_size", "3"},
        {"insertions_default_quality", "45"},
        {"low_quality_tail", "2"},
        {"maximum_cycle_value", "500"},
        {"mismatches_context_size", "2"},
        {"mismatches_default_quality", "-1"},
        {"no_standard_covs", "false"},
        {"quantizing_levels", "16"},
        {"recalibration_report", "null"},
        {"run_without_dbsnp", "false"},
        {"solid_nocall_strategy", "THROW_EXCEPTION"},
        {"solid_recal_mode", "SET_Q_ZERO"},
    };
    return defaults;
}

std::string argument_value(const std::vector<std::pair<std::string, std::string>>& arguments,
                           const std::string& name) {
    for (const auto& [key, value] : arguments)
        if (key == name) return value;
    return {};
}

int quantizing_levels(const std::vector<std::pair<std::string, std::string>>& arguments) {
    const auto text = argument_value(arguments, "quantizing_levels");
    if (text.empty()) return 16;
    try {
        const auto value = std::stoi(text);
        if (value < 1 || value > 94)
            throw std::runtime_error("BAD_INPUT: quantizing_levels must be in [1,94]");
        return value;
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("BAD_INPUT: malformed quantizing_levels argument");
    } catch (const std::out_of_range&) {
        throw std::runtime_error("BAD_INPUT: malformed quantizing_levels argument");
    }
}

void validate_structural_arguments(
    const std::vector<std::pair<std::string, std::string>>& first,
    const std::vector<std::pair<std::string, std::string>>& candidate,
    const std::string& path) {
    // These values define the dimensions of RecalibrationTables.  Java's
    // RecalUtils.combineTables rejects reports whose dimensions differ;
    // accepting them would merge incompatible covariate keys and produce a
    // report that ApplyBQSR interprets differently.
    for (const char* key : {"covariate", "no_standard_covs", "mismatches_context_size",
                            "indels_context_size", "maximum_cycle_value"}) {
        const auto left = argument_value(first, key);
        const auto right = argument_value(candidate, key);
        if (!left.empty() && !right.empty() && left != right)
            throw std::runtime_error("BAD_INPUT: incompatible BQSR report argument " +
                                     std::string(key) + " in " + path);
    }
}

std::map<std::string, ReadGroupSummary> summarize_read_groups(
    const CovariateTable& covariates, const std::array<Bin, 94>& bins,
    const ReadGroupTable* supplied_groups = nullptr) {
    if (supplied_groups != nullptr && !supplied_groups->empty()) return *supplied_groups;
    std::map<std::string, ReadGroupSummary> groups;
    for (const auto& [key, value] : covariates) {
        if (key.event != 'M') continue;
        // A Java report has separate Context and Cycle rows.  They describe
        // the same observations, so use only one family when reconstructing
        // RecalTable1 from a report that did not provide its own table.
        if (key.covariate == 'Y') continue;
        auto& destination = groups[key.read_group];
        if (key.quality >= 0 && key.quality < static_cast<int>(destination.bins.size())) {
            destination.bins[static_cast<std::size_t>(key.quality)].count += value.count;
            destination.bins[static_cast<std::size_t>(key.quality)].mismatches += value.mismatches;
        }
    }
    if (groups.empty()) groups.emplace(std::string{}, ReadGroupSummary{});
    std::uint64_t grouped = 0;
    for (const auto& [_, group] : groups)
        for (const auto& bin : group.bins) grouped += bin.count;
    if (grouped == 0) groups.begin()->second.bins = bins;
    return groups;
}

void write_gatk_report(const std::string& path,
                       const std::array<Bin, 94>& bins,
                       const CovariateTable& covariates,
                       const ReadGroupTable* supplied_groups = nullptr,
                       const ReportedDatumTable* supplied_reported = nullptr,
                       const std::vector<std::pair<std::string, std::string>>* supplied_arguments = nullptr) {
    const auto groups = summarize_read_groups(covariates, bins, supplied_groups);
    std::map<std::pair<std::string, char>, Bin> indel_groups;
    for (const auto& [key, value] : covariates) {
        if (key.event == 'M') continue;
        // GATK's standard Context and Cycle covariates both materialize the
        // same indel observations.  RecalTable1 has one aggregate row per
        // event, so a direct Java report must count only one family.  Native
        // sidecar rows use covariate='\0' and remain aggregate rows.
        if (key.covariate == 'Y') continue;
        auto& group = indel_groups[{key.read_group, key.event}];
        group.count += value.count;
        group.mismatches += value.mismatches;
    }
    std::size_t table1_rows = 0;
    std::size_t table0_rows = supplied_reported != nullptr && !supplied_reported->empty()
        ? supplied_reported->size() : 0;
    for (const auto& [read_group, group] : groups) {
        for (const auto& bin : group.bins) if (bin.count != 0) ++table1_rows;
        if (supplied_reported == nullptr || supplied_reported->empty()) ++table0_rows;
        for (const char event : {'I', 'D'}) {
            if (indel_groups.contains({read_group, event})) {
                ++table1_rows;
                if (supplied_reported == nullptr || supplied_reported->empty()) ++table0_rows;
            }
        }
    }
    std::size_t table2_rows = 0;
    for (const auto& [key, _] : covariates) {
        if (key.covariate == 'C' || key.covariate == 'Y') ++table2_rows;
        else if (key.event == 'M') table2_rows += 2;
        else table2_rows += (!key.context.empty() ? 1u : 0u) +
                            (key.cycle != std::numeric_limits<int>::min() ? 1u : 0u);
    }
    const auto& arguments = supplied_arguments != nullptr && !supplied_arguments->empty()
        ? *supplied_arguments : default_arguments();
    const auto levels = quantizing_levels(arguments);
    const auto quantized = quantize_gatk(bins, levels);
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open merged BQSR report: " + path);
    output << "#:GATKReport.v1.1:5\n"
           << "#:GATKTable:2:" << arguments.size() << ":%s:%s:;\n"
           << "#:GATKTable:Arguments:Recalibration argument collection values used in this run\n"
           << "Argument\tValue\n";
    for (const auto& [key, value] : arguments)
        output << key << '\t' << value << '\n';
    output << '\n'
           << "#:GATKTable:3:94:%d:%d:%d:;\n"
           << "#:GATKTable:Quantized:Quality quantization map\n"
           << "QualityScore\tCount\tQuantizedScore\n";
    for (std::size_t quality = 0; quality < bins.size(); ++quality) {
        output << quality << '\t' << quantized.counts[quality] << '\t'
               << quantized.values[quality] << '\n';
    }
    output << "\n#:GATKTable:6:" << table0_rows << ":%s:%s:%.4f:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable0:\n"
           << "ReadGroup\tEventType\tEmpiricalQuality\tEstimatedQReported\tObservations\tErrors\n";
    if (supplied_reported != nullptr && !supplied_reported->empty()) {
        for (const auto& [key, datum] : *supplied_reported) {
            const auto& [read_group, event] = key;
            const auto rg = display_read_group(read_group);
            // RecalUtils emits the insertion/deletion aggregate rows with
            // their default event quality (45), independent of the empirical
            // error count.  Preserve that Java writer boundary for direct
            // GATK reports while retaining the Bayesian estimate for M.
            const auto empirical = event == 'M'
                ? bayesian_empirical_quality(datum.bin, datum.reported_quality) : 45;
            output << rg << '\t' << event << '\t' << std::fixed << std::setprecision(4)
                   << static_cast<double>(empirical) << '\t' << datum.reported_quality
                   << '\t' << datum.bin.count << '\t' << std::setprecision(2)
                   << static_cast<double>(datum.bin.mismatches) << '\n';
        }
    } else for (const auto& [read_group, group] : groups) {
        std::uint64_t observations = 0, mismatches = 0;
        double weighted_quality = 0.0;
        for (std::size_t quality = 0; quality < group.bins.size(); ++quality) {
            observations += group.bins[quality].count;
            mismatches += group.bins[quality].mismatches;
            weighted_quality += static_cast<double>(quality) * group.bins[quality].count;
        }
        const double estimated = observations == 0 ? 0.0 : weighted_quality / observations;
        const auto rg = display_read_group(read_group);
        output << rg << "\tM\t" << std::fixed << std::setprecision(4)
               << empirical_quality(Bin{observations, mismatches}, static_cast<int>(std::lround(estimated)))
               << '\t' << estimated << '\t' << observations << '\t' << std::setprecision(2)
               << static_cast<double>(mismatches) << '\n';
        for (const char event : {'I', 'D'}) {
            const auto indel = indel_groups.find({read_group, event});
            if (indel == indel_groups.end()) continue;
            output << rg << '\t' << event << "\t45.0000\t45.0000\t"
                   << indel->second.count << '\t'
                   << static_cast<double>(indel->second.mismatches) << '\n';
        }
    }
    output << "\n#:GATKTable:6:" << table1_rows << ":%s:%s:%s:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable1:\n"
           << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations\tErrors\n";
    for (const auto& [read_group, group] : groups) {
        const auto rg = display_read_group(read_group);
        for (std::size_t quality = 0; quality < group.bins.size(); ++quality) {
            const auto& bin = group.bins[quality];
            if (bin.count == 0) continue;
            output << rg << '\t' << quality << "\tM\t" << std::fixed << std::setprecision(4)
                   // GATK's RecalTable1 empirical field is the covariate's
                   // reported quality for these rows; errors are carried in
                   // the separate Errors column.
                   << static_cast<double>(quality) << '\t' << bin.count
                   << '\t' << std::setprecision(2) << static_cast<double>(bin.mismatches) << '\n';
        }
        const auto insertion = indel_groups.find({read_group, 'I'});
        const auto deletion = indel_groups.find({read_group, 'D'});
        if (insertion != indel_groups.end())
            output << rg << "\t45\tI\t45.0000\t" << insertion->second.count
                   << "\t" << static_cast<double>(insertion->second.mismatches) << "\n";
        if (deletion != indel_groups.end())
            output << rg << "\t45\tD\t45.0000\t" << deletion->second.count
                   << "\t" << static_cast<double>(deletion->second.mismatches) << "\n";
    }
    output << "\n#:GATKTable:8:" << table2_rows << ":%s:%s:%s:%s:%s:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable2:\n"
           << "ReadGroup\tQualityScore\tCovariateValue\tCovariateName\tEventType\tEmpiricalQuality\tObservations\tErrors\n";
    for (const auto& [key, bin] : covariates) {
        const auto rg = display_read_group(key.read_group);
        const auto empirical = key.event == 'M' ? static_cast<double>(key.quality) : 45.0;
        const auto emit = [&](const char* name, const std::string& value) {
            output << rg << '\t' << key.quality << '\t' << value << '\t' << name << '\t'
                   << key.event << '\t' << std::fixed << std::setprecision(4) << empirical
                   << '\t' << bin.count << '\t' << std::setprecision(2)
                   << static_cast<double>(bin.mismatches) << '\n';
        };
        if (key.covariate == 'C') emit("Context", key.context);
        else if (key.covariate == 'Y') emit("Cycle", std::to_string(key.cycle));
        else {
            if (!key.context.empty()) emit("Context", key.context);
            if (key.cycle != std::numeric_limits<int>::min())
                emit("Cycle", std::to_string(key.cycle));
        }
    }
    if (!output) throw std::runtime_error("cannot write merged BQSR report: " + path);
}

void read_report_counts(const std::string& path, GatherInput& result) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open BQSR report: " + path);
    std::string line;
    enum class Section { None, Arguments, Table0, Table1, Table2 } section = Section::None;
    bool saw_gatk = false;
    while (std::getline(input, line)) {
        if (line.rfind("#:GATKTable:", 0) == 0) {
            saw_gatk = true;
            if (line.rfind("#:GATKTable:Arguments:", 0) == 0) section = Section::Arguments;
            else if (line.rfind("#:GATKTable:RecalTable0:", 0) == 0) section = Section::Table0;
            else if (line.rfind("#:GATKTable:RecalTable1:", 0) == 0) section = Section::Table1;
            else if (line.rfind("#:GATKTable:RecalTable2:", 0) == 0) section = Section::Table2;
            else section = Section::None;
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        if (section == Section::Arguments) {
            std::istringstream fields(line);
            std::string argument, value;
            if (fields >> argument >> value && argument != "Argument")
                result.arguments.emplace_back(std::move(argument), std::move(value));
            continue;
        }
        if (section == Section::Table1) {
            std::istringstream fields(line);
            std::string rg, event;
            int quality = 0;
            double empirical = 0.0, errors = 0.0;
            std::uint64_t observations = 0;
            if (fields >> rg >> quality >> event >> empirical >> observations >> errors &&
                quality >= 0 && quality < static_cast<int>(result.bins.size())) {
                result.has_recalibration_rows = true;
                const auto mismatches = static_cast<std::uint64_t>(std::llround(errors));
                if (event == "M") {
                    result.bins[static_cast<std::size_t>(quality)].count += observations;
                    result.bins[static_cast<std::size_t>(quality)].mismatches += mismatches;
                    auto& group = result.groups[rg].bins[static_cast<std::size_t>(quality)];
                    group.count += observations;
                    group.mismatches += mismatches;
                }
            }
            continue;
        }
        if (section == Section::Table0) {
            std::istringstream fields(line);
            std::string rg, event;
            double empirical = 0.0, reported_quality = 0.0, errors = 0.0;
            std::uint64_t observations = 0;
            if (!(fields >> rg >> event >> empirical >> reported_quality >> observations >> errors)) {
                if (line.rfind("ReadGroup", 0) == 0) continue;
                throw std::runtime_error("BAD_INPUT: malformed GATK RecalTable0 row: " + path);
            }
            if (event.empty()) throw std::runtime_error("BAD_INPUT: missing BQSR event type: " + path);
            const auto key = std::make_pair(rg == "UNKNOWN" ? std::string{} : rg, event.front());
            auto& datum = result.reported[key];
            const auto mismatch_count = static_cast<std::uint64_t>(std::llround(errors));
            if (!datum.initialized) {
                datum.reported_quality = reported_quality;
                datum.initialized = true;
            } else {
                const auto expected_errors = static_cast<double>(datum.bin.count) *
                    std::pow(10.0, -datum.reported_quality / 10.0) +
                    static_cast<double>(observations) * std::pow(10.0, -reported_quality / 10.0);
                const auto total = datum.bin.count + observations;
                datum.reported_quality = total == 0 ? 0.0 : -10.0 * std::log10(expected_errors / total);
            }
            datum.bin.count += observations;
            datum.bin.mismatches += mismatch_count;
            continue;
        }
        if (section == Section::Table2) {
            std::istringstream fields(line);
            std::string rg, covariate_name, event, value;
            int quality = 0;
            double empirical = 0.0, errors = 0.0;
            std::uint64_t observations = 0;
            if (!(fields >> rg >> quality >> value >> covariate_name >> event >> empirical
                  >> observations >> errors)) {
                if (line.rfind("ReadGroup", 0) == 0) continue;
                throw std::runtime_error("BAD_INPUT: malformed GATK RecalTable2 row: " + path);
            }
            result.has_recalibration_rows = true;
            CovariateKey key;
            key.read_group = rg == "UNKNOWN" ? std::string{} : rg;
            key.quality = quality;
            key.event = event.empty() ? 'M' : event.front();
            if (covariate_name == "Context") {
                key.covariate = 'C';
                key.context = value;
                key.cycle = std::numeric_limits<int>::min();
            } else if (covariate_name == "Cycle") {
                key.covariate = 'Y';
                try { key.cycle = std::stoi(value); }
                catch (...) { throw std::runtime_error("BAD_INPUT: malformed Cycle covariate: " + path); }
            } else {
                throw std::runtime_error("BAD_INPUT: unsupported GATK covariate: " + covariate_name);
            }
            auto& destination = result.covariates[key];
            destination.count += observations;
            destination.mismatches += static_cast<std::uint64_t>(std::llround(errors));
            continue;
        }
        std::istringstream fields(line);
        int quality = 0;
        std::uint64_t count = 0, mismatches = 0;
        double empirical = 0.0;
        int delta = 0;
        if (!saw_gatk && fields >> quality >> count >> mismatches >> empirical >> delta &&
            quality >= 0 && quality < static_cast<int>(result.bins.size())) {
            result.bins[static_cast<std::size_t>(quality)].count += count;
            result.bins[static_cast<std::size_t>(quality)].mismatches += mismatches;
        }
    }
}

void merge_covariates(const std::string& path, CovariateTable& table) {
    std::ifstream input(path);
    if (!input) return;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        CovariateKey key;
        Bin bin;
        double empirical = 0.0;
        int delta = 0;
        if (!(fields >> key.read_group >> key.cycle >> key.context >> key.quality >>
              bin.count >> bin.mismatches >> empirical >> delta))
            throw std::runtime_error("BAD_INPUT: malformed BQSR covariate row: " + path);
        std::string event;
        if (fields >> event && !event.empty()) key.event = event.front();
        if (key.read_group == "-") key.read_group.clear();
        key.covariate = '\0';
        auto& destination = table[key];
        destination.count += bin.count;
        destination.mismatches += bin.mismatches;
    }
}

int run(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto wall_begin = std::chrono::steady_clock::now();
    std::array<Bin, 94> bins{};
    CovariateTable covariates;
    ReadGroupTable groups;
    ReportedDatumTable reported;
    std::vector<std::pair<std::string, std::string>> arguments;
    std::size_t direct_gatk_report_inputs = 0;
    bool has_recalibration_rows = false;
    for (const auto& path : options.inputs) {
        GatherInput input;
        read_report_counts(path, input);
        has_recalibration_rows = has_recalibration_rows || input.has_recalibration_rows;
        if (input.arguments.empty()) input.arguments = default_arguments();
        if (arguments.empty()) arguments = input.arguments;
        else validate_structural_arguments(arguments, input.arguments, path);
        std::ifstream sidecar(path + ".covariates.tsv");
        if (sidecar.good()) {
            // Native reports have a richer sidecar representation.  Prefer it
            // over the expanded RecalTable2 rows to avoid counting each
            // context/cycle record twice.
            input.covariates.clear();
            merge_covariates(path + ".covariates.tsv", input.covariates);
        } else if (!input.covariates.empty()) ++direct_gatk_report_inputs;
        for (std::size_t quality = 0; quality < bins.size(); ++quality) {
            bins[quality].count += input.bins[quality].count;
            bins[quality].mismatches += input.bins[quality].mismatches;
        }
        for (const auto& [key, value] : input.covariates) {
            auto& destination = covariates[key];
            destination.count += value.count;
            destination.mismatches += value.mismatches;
        }
        for (const auto& [key, value] : input.reported) {
            auto& destination = reported[key];
            if (!destination.initialized) {
                destination = value;
            } else {
                const auto expected_errors = static_cast<double>(destination.bin.count) *
                    std::pow(10.0, -destination.reported_quality / 10.0) +
                    static_cast<double>(value.bin.count) * std::pow(10.0, -value.reported_quality / 10.0);
                const auto total = destination.bin.count + value.bin.count;
                destination.reported_quality = total == 0 ? 0.0 : -10.0 * std::log10(expected_errors / total);
                destination.bin.count = total;
                destination.bin.mismatches += value.bin.mismatches;
            }
        }
        for (const auto& [read_group, summary] : input.groups) {
            auto& destination = groups[read_group];
            for (std::size_t quality = 0; quality < summary.bins.size(); ++quality) {
                destination.bins[quality].count += summary.bins[quality].count;
                destination.bins[quality].mismatches += summary.bins[quality].mismatches;
            }
        }
    }
    // RecalibrationReport.gatherReports() filters empty reports and throws
    // "there is no usable data in any input file" when every report has only
    // headers/Arguments/Quantized tables.  Preserve that fail-closed writer
    // contract before creating any output or sidecar.
    if (!has_recalibration_rows)
        throw std::runtime_error("BAD_INPUT: there is no usable data in any input file");
    write_gatk_report(options.output, bins, covariates,
                      groups.empty() ? nullptr : &groups,
                      reported.empty() ? nullptr : &reported,
                      &arguments);
    const auto covariate_path = options.output + ".covariates.tsv";
    std::ofstream covariate_output(covariate_path);
    if (!covariate_output) throw std::runtime_error("cannot open merged BQSR covariate table: " + covariate_path);
    covariate_output << "# FASTGATK-BQSR-COVARIATES v2\n"
                     << "# read_group\tcycle\tcontext\traw_quality\tcount\tmismatches\tempirical_quality\tdelta\tevent_type\n"
                     << std::setprecision(12);
    for (const auto& [key, bin] : covariates) {
        const double error = bin.count == 0 ? 0.0 :
            static_cast<double>(bin.mismatches) / static_cast<double>(bin.count);
        const double empirical = bin.count == 0 ? static_cast<double>(key.quality) :
            -10.0 * std::log10(std::max(error, 1.0e-6));
        const auto delta = static_cast<int>(std::lround(std::clamp(
            empirical - static_cast<double>(key.quality), -20.0, 20.0)));
        covariate_output << (key.read_group.empty() ? "-" : key.read_group) << '\t'
                         << key.cycle << '\t' << key.context << '\t' << key.quality << '\t'
                         << bin.count << '\t' << bin.mismatches << '\t' << empirical << '\t'
                         << delta << '\t' << key.event << '\n';
    }
    covariate_output.close();
    std::error_code output_error;
    const auto output_bytes = std::filesystem::file_size(options.output, output_error);
    if (output_error || output_bytes == 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: BQSR merged report is missing or empty");
    const auto covariate_bytes = std::filesystem::file_size(covariate_path, output_error);
    if (output_error || covariate_bytes == 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: BQSR covariate table is missing or empty");
    const auto wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_begin).count();
    const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("cannot write BQSR gather manifest: " + manifest_path);
    std::uint64_t total_count = 0, total_mismatches = 0;
    for (const auto& bin : bins) { total_count += bin.count; total_mismatches += bin.mismatches; }
    manifest << "{\"schema_version\":1,\"tool\":\"GatherBQSRReports\"," 
             << "\"implementation\":\"fastgatk-gather-bqsr-reports\",\"status\":\"prototype\","
             << "\"execution_space\":\"Host\",\"determinism\":\"strict\","
             << "\"primary_output\":\"" << json_escape(options.output) << "\","
             << "\"primary_output_kind\":\"recalibration-report\","
             << "\"compatibility\":{\"integer_table\":true,\"gatk_report_v1_1\":true,\"covariate_table\":true,\"hierarchical_covariate_model\":true,\"deterministic_merge\":true,\"direct_java_gatk_report\":true,"
             << "\"quiet\":" << (options.quiet ? "true" : "false") << ",\"tmp_dir\":\""
             << json_escape(options.tmp_dir) << "\",\"verbosity\":\"" << json_escape(options.verbosity) << "\","
             << "\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
             << json_escape(options.output) << "\",\"kind\":\"recalibration-report\",\"complete\":true},{\"path\":\""
             << json_escape(covariate_path) << "\",\"kind\":\"bqsr-covariate-table\",\"complete\":true}],"
             << "\"telemetry\":{\"resources\":" << resources.to_json()
             << ",\"input_reports\":" << options.inputs.size()
             << ",\"direct_gatk_report_inputs\":" << direct_gatk_report_inputs
             << ",\"quantizing_levels\":" << quantizing_levels(arguments)
             << ",\"observations\":" << total_count << ",\"mismatches\":" << total_mismatches
             << ",\"output_bytes\":" << output_bytes
             << ",\"covariate_bytes\":" << covariate_bytes
             << ",\"quiet\":" << (options.quiet ? "true" : "false")
             << ",\"tmp_dir\":\"" << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
             << json_escape(options.verbosity) << "\""
             << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}}\n";
    manifest.close();
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize BQSR gather manifest");
    std::cout << "{\"tool\":\"GatherBQSRReports\",\"status\":\"prototype\",\"execution_space\":\"Host\",\"input_reports\":"
              << options.inputs.size() << ",\"observations\":" << total_count
              << ",\"covariate_rows\":" << covariates.size()
              << ",\"output_bytes\":" << output_bytes
              << ",\"covariate_bytes\":" << covariate_bytes
              << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse(argc, argv);
        return run(options, fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
