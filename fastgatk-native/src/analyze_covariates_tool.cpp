#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::string before;
    std::string after;
    std::string bqsr;
    std::string csv;
    std::string plots;
    std::string manifest;
    bool ignore_last_modification_times = false;
};

struct RowKey {
    std::string read_group;
    std::string covariate_value;
    std::string covariate_name;
    char event = 'M';

    static int covariate_rank(const std::string& name) {
        if (name == "Context") return 0;
        if (name == "Cycle") return 1;
        return 2; // QualityScore
    }

    static int base_rank(char base) {
        switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
        }
    }

    static int context_rank(const std::string& value) {
        int key = static_cast<int>(value.size());
        int shift = 4;
        for (const auto base : value) {
            const auto rank = base_rank(base);
            if (rank < 0) return std::numeric_limits<int>::max();
            key |= rank << shift;
            shift += 2;
        }
        return key;
    }

    static int cycle_rank(const std::string& value) {
        try {
            const auto cycle = std::stoi(value);
            return 2 * std::abs(cycle) + (cycle < 0 ? 1 : 0);
        } catch (...) {
            return std::numeric_limits<int>::max();
        }
    }

    static int value_rank(const std::string& value, const std::string& name) {
        if (name == "Context") return context_rank(value);
        if (name == "Cycle") return cycle_rank(value);
        try { return std::stoi(value); }
        catch (...) { return std::numeric_limits<int>::max(); }
    }

    static int event_rank(char value) {
        return value == 'M' ? 0 : value == 'I' ? 1 : value == 'D' ? 2 : 3;
    }

    bool operator<(const RowKey& other) const {
        if (read_group != other.read_group) return read_group < other.read_group;
        const auto left_covariate = covariate_rank(covariate_name);
        const auto right_covariate = covariate_rank(other.covariate_name);
        if (left_covariate != right_covariate) return left_covariate < right_covariate;
        const auto left_value = value_rank(covariate_value, covariate_name);
        const auto right_value = value_rank(other.covariate_value, other.covariate_name);
        if (left_value != right_value) return left_value < right_value;
        if (event_rank(event) != event_rank(other.event))
            return event_rank(event) < event_rank(other.event);
        return covariate_value < other.covariate_value;
    }
};

struct InputRow {
    RowKey key;
    std::uint64_t observations = 0;
    double errors = 0.0;
    double reported_quality = 0.0;
};

struct Report {
    std::string path;
    std::map<std::string, std::string> arguments;
    std::map<RowKey, InputRow> rows;
};

struct OutputRow {
    RowKey key;
    std::uint64_t observations = 0;
    double errors = 0.0;
    double empirical_quality = 0.0;
    double reported_quality = 0.0;
    double accuracy = 0.0;
    std::string mode;
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
    if ((argument == name || (short_name != nullptr && argument == short_name)) &&
        index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
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

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-analyze-covariates\n"
                         "  -before, --before-report-file FILE       first-pass BQSR report\n"
                         "  -after, --after-report-file FILE        second-pass BQSR report\n"
                         "  -bqsr, --bqsr-recal-file FILE           single BQSR report\n"
                         "  -csv, --intermediate-csv-file FILE      CSV output\n"
                         "  -plots, --plots-report-file FILE        PDF plot output\n"
                         "      --ignore-last-modification-times    suppress timestamp warning\n"
                         "      --output-manifest FILE              OutputManifest JSON\n";
            std::exit(0);
        } else if (has_option(argument, "--before-report-file") || argument == "-before") {
            options.before = require_value(index, argc, argv, argument,
                                           "--before-report-file", "-before");
        } else if (has_option(argument, "--after-report-file") || argument == "-after") {
            options.after = require_value(index, argc, argv, argument,
                                          "--after-report-file", "-after");
        } else if (argument.rfind("-bqsr=", 0) == 0) {
            // Barclay accepts the short alias only as a separate token; an
            // embedded '=' is rejected by GATK 4.6.2.0 before tool parsing.
            throw std::invalid_argument("unknown option: " + argument);
        } else if (has_option(argument, "--bqsr-recal-file") ||
                   has_option(argument, "-bqsr") || argument == "-BQSR") {
            const char* short_name = argument == "-bqsr" || has_option(argument, "-bqsr")
                ? "-bqsr" : "-BQSR";
            options.bqsr = require_value(index, argc, argv, argument,
                                         "--bqsr-recal-file", short_name);
        } else if (has_option(argument, "--intermediate-csv-file") || argument == "-csv") {
            options.csv = require_value(index, argc, argv, argument,
                                        "--intermediate-csv-file", "-csv");
        } else if (has_option(argument, "--plots-report-file") || argument == "-plots") {
            options.plots = require_value(index, argc, argv, argument,
                                          "--plots-report-file", "-plots");
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--ignore-last-modification-times") {
            options.ignore_last_modification_times = true;
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Launcher-compatible no-op flags.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.before.empty() && options.after.empty() && options.bqsr.empty())
        throw std::invalid_argument("at least one of -before, -after, or -BQSR is required");
    if (options.csv.empty() && options.plots.empty())
        throw std::invalid_argument("at least one output is required: -csv or -plots");
    return options;
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), ending) == 0;
}

template <typename Fn>
void for_each_line(const std::string& path, Fn&& callback) {
#if FASTGATK_HAS_ZLIB
    if (suffix(path, ".gz")) {
        gzFile stream = gzopen(path.c_str(), "rb");
        if (stream == nullptr) throw std::runtime_error("BAD_INPUT: cannot open report: " + path);
        char buffer[1 << 16];
        while (gzgets(stream, buffer, sizeof(buffer)) != nullptr) {
            std::string line(buffer);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            callback(line);
        }
        const auto error = gzclose(stream);
        if (error != Z_OK) throw std::runtime_error("BAD_INPUT: cannot read compressed report: " + path);
        return;
    }
#else
    if (suffix(path, ".gz"))
        throw std::runtime_error("BACKEND_UNAVAILABLE: compressed AnalyzeCovariates reports require zlib");
#endif
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open report: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        callback(line);
    }
    if (!input.eof()) throw std::runtime_error("BAD_INPUT: cannot read report: " + path);
}

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    for (std::string value; stream >> value;) result.push_back(std::move(value));
    return result;
}

double bayesian_quality(std::uint64_t observations, double errors, double prior) {
    if (observations == 0) return std::clamp(std::round(prior), 0.0, 40.0);
    const auto rounded_errors = std::floor(errors + 0.5) + 1.0;
    const auto smoothed_observations = static_cast<double>(observations) + 2.0;
    double best_score = -std::numeric_limits<double>::infinity();
    int best_quality = 0;
    for (int quality = 0; quality <= 40; ++quality) {
        const auto probability = std::pow(10.0, -static_cast<double>(quality) / 10.0);
        const auto difference = std::min(std::abs(static_cast<int>(quality - prior)), 40);
        const auto log_prior = -0.5 * static_cast<double>(difference * difference) / 0.25;
        const auto log_likelihood = std::lgamma(smoothed_observations + 1.0) -
            std::lgamma(rounded_errors + 1.0) -
            std::lgamma(smoothed_observations - rounded_errors + 1.0) +
            rounded_errors * std::log(probability) +
            (smoothed_observations - rounded_errors) * std::log1p(-probability);
        if (log_prior + log_likelihood > best_score) {
            best_score = log_prior + log_likelihood;
            best_quality = quality;
        }
    }
    return static_cast<double>(best_quality);
}

void add_row(Report& report, const RowKey& key, std::uint64_t observations,
             double errors, double reported_quality) {
    auto& destination = report.rows[key];
    if (destination.key.covariate_name.empty()) destination.key = key;
    const auto previous = static_cast<double>(destination.observations);
    const auto expected_errors = previous * std::pow(10.0, -destination.reported_quality / 10.0) +
        static_cast<double>(observations) * std::pow(10.0, -reported_quality / 10.0);
    destination.observations += observations;
    destination.errors += errors;
    const auto total = static_cast<double>(destination.observations);
    destination.reported_quality = total == 0.0 || expected_errors <= 0.0
        ? 0.0 : -10.0 * std::log10(expected_errors / total);
}

Report read_report(const std::string& path) {
    Report report;
    report.path = path;
    enum class Section { None, Arguments, Table1, Table2 } section = Section::None;
    bool saw_gatk = false;
    for_each_line(path, [&](const std::string& line) {
        if (line.rfind("#:GATKTable:", 0) == 0) {
            saw_gatk = true;
            if (line.rfind("#:GATKTable:Arguments:", 0) == 0) section = Section::Arguments;
            else if (line.rfind("#:GATKTable:RecalTable1:", 0) == 0) section = Section::Table1;
            else if (line.rfind("#:GATKTable:RecalTable2:", 0) == 0) section = Section::Table2;
            else section = Section::None;
            return;
        }
        if (line.empty() || line[0] == '#') return;
        const auto fields = split_fields(line);
        if (section == Section::Arguments) {
            if (fields.size() >= 2 && fields[0] != "Argument")
                report.arguments[fields[0]] = fields[1];
            return;
        }
        if (section == Section::Table1) {
            if (fields.size() < 6 || fields[0] == "ReadGroup") return;
            try {
                const auto quality = std::stoi(fields[1]);
                const auto event = fields[2].empty() ? 'M' : fields[2][0];
                const auto empirical = std::stod(fields[3]);
                const auto observations = static_cast<std::uint64_t>(std::stoull(fields[4]));
                const auto errors = std::stod(fields[5]);
                add_row(report, RowKey{fields[0], std::to_string(quality), "QualityScore", event},
                        observations, errors, static_cast<double>(quality));
                (void)empirical;
            } catch (...) {
                throw std::runtime_error("BAD_INPUT: malformed GATK RecalTable1 row: " + path);
            }
            return;
        }
        if (section == Section::Table2) {
            if (fields.size() < 8 || fields[0] == "ReadGroup") return;
            try {
                const auto quality = std::stod(fields[1]);
                const auto event = fields[4].empty() ? 'M' : fields[4][0];
                const auto empirical = std::stod(fields[5]);
                const auto observations = static_cast<std::uint64_t>(std::stoull(fields[6]));
                const auto errors = std::stod(fields[7]);
                add_row(report, RowKey{fields[0], fields[2], fields[3], event},
                        observations, errors, quality);
                (void)empirical;
            } catch (...) {
                throw std::runtime_error("BAD_INPUT: malformed GATK RecalTable2 row: " + path);
            }
            return;
        }
        if (!saw_gatk) {
            // Native sidecar-like TSV is accepted as a narrow compatibility
            // path: read_group cycle context quality count errors ... event.
            if (fields.size() >= 9) {
                try {
                    const auto event = fields[8].empty() ? 'M' : fields[8][0];
                    add_row(report, RowKey{fields[0], fields[2], "Context", event},
                            std::stoull(fields[4]), std::stod(fields[5]), std::stod(fields[3]));
                } catch (...) {
                    throw std::runtime_error("BAD_INPUT: malformed BQSR covariate row: " + path);
                }
            }
        }
    });
    if (report.rows.empty()) throw std::runtime_error("BAD_INPUT: report has no RecalTable1/2 rows: " + path);
    return report;
}

std::vector<OutputRow> expand_rows(const Report& report, const std::string& mode) {
    std::vector<OutputRow> rows;
    rows.reserve(report.rows.size());
    for (const auto& [key, row] : report.rows) {
        OutputRow value;
        value.key = key;
        value.observations = row.observations;
        value.errors = row.errors;
        value.reported_quality = row.reported_quality;
        value.empirical_quality = bayesian_quality(row.observations, row.errors,
                                                   row.reported_quality);
        value.accuracy = value.empirical_quality - value.reported_quality;
        value.mode = mode;
        rows.push_back(std::move(value));
    }
    return rows;
}

bool arguments_compatible(const Report& reference, const Report& candidate,
                          std::string& mismatch) {
    // Match RecalibrationArgumentCollection.compareReportArguments() in
    // GATK 4.6.2.0 exactly.  The Java implementation deliberately compares
    // this fixed set (and does not compare covariate/indels_context_size or
    // report bookkeeping); comparing every key would reject valid GATK
    // before/after pairs that Java accepts.
    static const std::vector<std::string> compared = {
        "no_standard_covs", "run_without_dbsnp", "solid_recal_mode",
        "solid_nocall_strategy", "mismatches_context_size",
        "mismatches_default_quality", "deletions_default_quality",
        "insertions_default_quality", "maximum_cycle_value",
        "low_quality_tail", "default_platform", "force_platform",
        "quantizing_levels", "binary_tag_name"
    };
    // Reports produced by current GATK contain all compared keys.  Defaults
    // preserve Java's constructor semantics for older reports that omit a
    // key, while keeping comparison fail-closed for non-default values.
    static const std::map<std::string, std::string> defaults = {
        {"no_standard_covs", "false"}, {"run_without_dbsnp", "false"},
        {"solid_recal_mode", "SET_Q_ZERO"},
        {"solid_nocall_strategy", "THROW_EXCEPTION"},
        {"mismatches_context_size", "2"},
        {"mismatches_default_quality", "-1"},
        {"deletions_default_quality", "45"},
        {"insertions_default_quality", "45"},
        {"maximum_cycle_value", "500"}, {"low_quality_tail", "2"},
        {"default_platform", "null"}, {"force_platform", "null"},
        {"quantizing_levels", "16"}, {"binary_tag_name", "null"}
    };
    const auto value_or_default = [&](const Report& report,
                                      const std::string& key) -> std::string {
        const auto found = report.arguments.find(key);
        if (found != report.arguments.end()) return found->second;
        return defaults.at(key);
    };
    for (const auto& key : compared) {
        const auto left = value_or_default(reference, key);
        const auto right = value_or_default(candidate, key);
        if (left != right) {
            mismatch = key + " (" + left + " vs " + right + ")";
            return false;
        }
    }
    return true;
}

const char* event_name(char event) {
    switch (event) {
    case 'I': return "Base Insertion";
    case 'D': return "Base Deletion";
    default: return "Base Substitution";
    }
}

void write_csv(const std::string& path, const std::vector<OutputRow>& rows) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open CSV: " + path);
    output << "ReadGroup,CovariateValue,CovariateName,EventType,Observations,Errors,"
              "EmpiricalQuality,AverageReportedQuality,Accuracy,Recalibration\n";
    output << std::fixed << std::setprecision(2);
    for (const auto& row : rows) {
        output << (row.key.read_group == "UNKNOWN" ? "" : row.key.read_group) << ','
               << row.key.covariate_value << ',' << row.key.covariate_name << ','
               << event_name(row.key.event) << ',' << row.observations << ',' << row.errors << ','
               << row.empirical_quality << ',' << row.reported_quality << ','
               << row.accuracy << ',' << row.mode << '\n';
    }
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write CSV: " + path);
}

void write_pdf(const std::string& path, const std::vector<OutputRow>& rows) {
    // A deterministic vector report: the CSV is the machine-readable output;
    // this PDF is intentionally dependency-free and replaces the R plotting
    // bridge with a bounded overview of empirical-vs-reported qualities.
    const double max_quality = std::max(40.0, std::accumulate(rows.begin(), rows.end(), 0.0,
        [](double current, const OutputRow& row) {
            return std::max({current, row.empirical_quality, row.reported_quality});
        }));
    std::ostringstream content;
    content << "BT /F1 12 Tf 36 760 Td (fastgatk AnalyzeCovariates) Tj ET\n"
            << "0.8 w 48 80 m 48 720 l 560 720 l S\n";
    const auto limit = std::min<std::size_t>(rows.size(), 240);
    for (std::size_t index = 0; index < limit; ++index) {
        const auto x = 52.0 + 504.0 * static_cast<double>(index) /
            std::max<std::size_t>(1, limit);
        const auto y_reported = 90.0 + 600.0 * rows[index].reported_quality / max_quality;
        const auto y_empirical = 90.0 + 600.0 * rows[index].empirical_quality / max_quality;
        content << std::fixed << std::setprecision(2)
                << x << ' ' << y_reported << " m " << x << ' ' << y_empirical << " l S\n";
    }
    content << "BT /F1 9 Tf 48 58 Td (reported-to-empirical quality segments; see CSV for exact rows) Tj ET\n";
    const std::string stream = content.str();
    std::vector<std::string> objects{
        "<< /Type /Catalog /Pages 2 0 R >>",
        "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>",
        "<< /Length " + std::to_string(stream.size()) + " >>\nstream\n" + stream + "endstream",
        "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    };
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open PDF: " + path);
    output << "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
    std::vector<std::streamoff> offsets(1, 0);
    for (std::size_t index = 0; index < objects.size(); ++index) {
        offsets.push_back(output.tellp());
        output << (index + 1) << " 0 obj\n" << objects[index] << "\nendobj\n";
    }
    const auto xref = output.tellp();
    output << "xref\n0 " << (objects.size() + 1) << "\n0000000000 65535 f \n";
    for (std::size_t index = 1; index < offsets.size(); ++index)
        output << std::setw(10) << std::setfill('0') << offsets[index]
               << " 00000 n \n";
    output << std::setfill(' ') << "trailer\n<< /Size " << (objects.size() + 1)
           << " /Root 1 0 R >>\nstartxref\n" << xref << "\n%%EOF\n";
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write PDF: " + path);
}

int run(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto wall_begin = std::chrono::steady_clock::now();
    std::vector<std::pair<std::string, std::string>> paths;
    if (!options.bqsr.empty()) paths.emplace_back("BQSR", options.bqsr);
    if (!options.before.empty()) paths.emplace_back("Before", options.before);
    if (!options.after.empty()) paths.emplace_back("After", options.after);
    std::vector<Report> reports;
    reports.reserve(paths.size());
    std::vector<OutputRow> output_rows;
    for (const auto& [role, path] : paths) {
        if (!std::filesystem::is_regular_file(path))
            throw std::runtime_error("BAD_INPUT: report is not a regular file: " + path);
        reports.push_back(read_report(path));
        const auto rows = expand_rows(reports.back(), role);
        output_rows.insert(output_rows.end(), rows.begin(), rows.end());
    }
    if (reports.size() > 1) {
        for (std::size_t index = 1; index < reports.size(); ++index) {
            std::string mismatch;
            if (!arguments_compatible(reports.front(), reports[index], mismatch))
                throw std::runtime_error("BAD_INPUT: incompatible AnalyzeCovariates report arguments: " + mismatch);
        }
    }
    if (!options.ignore_last_modification_times && !options.before.empty() && !options.after.empty()) {
        std::error_code error;
        const auto before_time = std::filesystem::last_write_time(options.before, error);
        const auto after_time = std::filesystem::last_write_time(options.after, error);
        if (!error && before_time > after_time)
            std::cerr << "warning: before report has a newer modification time than after report\n";
    }

    // Keep the numeric portion of report analysis on the selected Kokkos
    // execution space. Host owns report parsing and deterministic row order.
    Kokkos::View<double*> empirical("analyze_empirical", output_rows.size());
    Kokkos::View<double*> reported("analyze_reported", output_rows.size());
    Kokkos::View<double*> accuracy("analyze_accuracy", output_rows.size());
    auto empirical_host = Kokkos::create_mirror_view(empirical);
    auto reported_host = Kokkos::create_mirror_view(reported);
    for (std::size_t index = 0; index < output_rows.size(); ++index) {
        empirical_host(index) = output_rows[index].empirical_quality;
        reported_host(index) = output_rows[index].reported_quality;
    }
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch analysis_host("analyze-covariates-v1");
    analysis_host.records = output_rows.size();
    analysis_host.bytes = output_rows.size() * 3 * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> analysis_plan("analyze-covariates");
    analysis_plan.begin_prepare(analysis_host);
    Kokkos::deep_copy(empirical, empirical_host);
    Kokkos::deep_copy(reported, reported_host);
    fastgatk::core::DeviceBatch<ExecSpace> analysis_device(output_rows.size());
    analysis_device.bind("empirical", empirical);
    analysis_device.bind("reported", reported);
    analysis_device.bind("accuracy", accuracy);
    ExecSpace().fence();
    analysis_plan.end_prepare(analysis_device);
    analysis_plan.begin_execute();
    Kokkos::parallel_for("analyze_covariates_accuracy",
        Kokkos::RangePolicy<ExecSpace>(0, output_rows.size()),
        KOKKOS_LAMBDA(const std::size_t index) { accuracy(index) = empirical(index) - reported(index); });
    ExecSpace().fence();
    analysis_plan.end_execute();
    auto accuracy_host = Kokkos::create_mirror_view(accuracy);
    Kokkos::deep_copy(accuracy_host, accuracy);
    for (std::size_t index = 0; index < output_rows.size(); ++index)
        output_rows[index].accuracy = accuracy_host(index);
    const auto prepare_seconds = analysis_plan.telemetry().prepare_seconds;
    const auto execute_seconds = analysis_plan.telemetry().execute_seconds;

    if (!options.csv.empty()) write_csv(options.csv, output_rows);
    if (!options.plots.empty()) write_pdf(options.plots, output_rows);
    const auto file_bytes = [](const std::string& path) -> std::uintmax_t {
        if (path.empty()) return 0;
        std::error_code error;
        const auto bytes = std::filesystem::file_size(path, error);
        if (error || bytes == 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: AnalyzeCovariates output is missing or empty: " + path);
        return bytes;
    };
    const auto csv_bytes = file_bytes(options.csv);
    const auto pdf_bytes = file_bytes(options.plots);
    const auto wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_begin).count();
    const auto manifest_path = options.manifest.empty()
        ? (!options.csv.empty() ? options.csv : options.plots) + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + manifest_path);
    manifest << "{\"schema_version\":1,\"tool\":\"AnalyzeCovariates\","
             << "\"implementation\":\"fastgatk-analyze-covariates\",\"status\":\"prototype\","
             << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\","
             << "\"determinism\":\"strict\","
             << "\"compatibility\":{\"gatk_report_v1_1\":true,\"before_after\":"
             << (!options.before.empty() && !options.after.empty() ? "true" : "false")
             << ",\"csv\":" << (!options.csv.empty() ? "true" : "false")
             << ",\"pdf\":" << (!options.plots.empty() ? "true" : "false")
             << ",\"bit_identical_to_gatk\":false},\"outputs\":[";
    bool first = true;
    const auto output_entry = [&](const std::string& path, const char* kind) {
        if (!first) manifest << ',';
        first = false;
        manifest << "{\"path\":\"" << json_escape(path) << "\",\"kind\":\""
                 << kind << "\",\"complete\":true}";
    };
    if (!options.csv.empty()) output_entry(options.csv, "analyze-covariates-csv");
    if (!options.plots.empty()) output_entry(options.plots, "analyze-covariates-pdf");
    manifest << "],\"telemetry\":{\"resources\":" << resources.to_json()
             << ",\"input_reports\":" << reports.size()
             << ",\"rows\":" << output_rows.size()
             << ",\"prepare_seconds\":" << prepare_seconds
             << ",\"execute_seconds\":" << execute_seconds
             << ",\"csv_bytes\":" << csv_bytes
             << ",\"pdf_bytes\":" << pdf_bytes
             << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds
             << ",\"execution_space\":\""
             << json_escape(Kokkos::DefaultExecutionSpace::name()) << "\""
             << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
             << ",\"kernel_execution_policy\":\"RangePolicy\""
             << ",\"kernel_batches\":" << analysis_plan.telemetry().execute_calls
             << ",\"kernel_observations\":" << output_rows.size() << "}}\n";
    manifest.close();
    if (!manifest || !std::filesystem::is_regular_file(manifest_path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete AnalyzeCovariates manifest");
    std::cout << "{\"tool\":\"AnalyzeCovariates\",\"status\":\"prototype\",\"input_reports\":"
              << reports.size() << ",\"rows\":" << output_rows.size()
              << ",\"csv_bytes\":" << csv_bytes << ",\"pdf_bytes\":" << pdf_bytes
              << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds
              << ",\"kernel_batches\":" << analysis_plan.telemetry().execute_calls
              << ",\"kernel_observations\":" << output_rows.size()
              << ",\"kernel_execution_space\":\"" << ExecSpace::name() << "\"}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Kokkos::ScopeGuard guard(argc, argv);
    try {
        return run(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
