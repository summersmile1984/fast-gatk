#include <Kokkos_Core.hpp>

#include "fastgatk/core/plan.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::vector<std::string> denoised_copy_ratios;
    std::string allelic_counts;
    std::vector<std::string> allelic_count_files;
    std::string normal_allelic_counts;
    std::size_t min_total_allele_count_case = 0;
    std::size_t min_total_allele_count_normal = 30;
    double genotyping_homozygous_log_ratio_threshold = -10.0;
    double genotyping_base_error_rate = 0.05;
    std::string output;
    std::string output_directory;
    std::string output_prefix;
    std::string manifest;
    std::string sample;
    double change_point_threshold = 0.15;
    double allelic_change_point_threshold = 0.05;
    // GATK always uses KernelSegmenter when --segments is absent.  The
    // threshold options are native compatibility extensions; retain the
    // historical threshold path only when a caller explicitly supplies one.
    bool change_point_threshold_overridden = false;
    // GATK's advanced --segments input skips kernel segmentation and uses the
    // supplied Picard interval list as the modeling partition.
    std::string input_segments;
    std::size_t min_points = 1;
    std::size_t max_segments_per_contig = 1000;
    double kernel_variance_copy_ratio = 0.0;
    double kernel_variance_allele_fraction = 0.025;
    double kernel_scaling_allele_fraction = 1.0;
    std::size_t kernel_approximation_dimension = 100;
    std::vector<std::size_t> window_sizes{8, 16, 32, 64, 128, 256};
    double changepoint_penalty_factor = 1.0;
    bool kernel_requested = false;
    // The probabilistic ModelSegments controls use a bounded deterministic
    // posterior sampler.  The default threshold/KernelSegmenter path remains
    // byte-for-byte compatible when these options are not requested.
    std::string mode = "BOTH";
    std::size_t num_burn_in_iterations = 100;
    std::size_t num_samples = 1000;
    // GATK exposes separate copy-ratio/allele-fraction MCMC controls.  The
    // native bounded sampler shares one deterministic chain, so retain both
    // requested surfaces and use the largest effective chain dimensions.
    double minor_allele_fraction_prior_alpha = 25.0;
    std::size_t num_samples_copy_ratio = 100;
    std::size_t num_burn_in_samples_copy_ratio = 50;
    std::size_t num_samples_allele_fraction = 100;
    std::size_t num_burn_in_samples_allele_fraction = 50;
    bool gatk_mcmc_overridden = false;
    bool generic_mcmc_overridden = false;
    std::size_t maximum_number_of_smoothing_iterations = 25;
    std::size_t number_of_smoothing_iterations_per_fit = 0;
    double smoothing_credible_interval_threshold_copy_ratio = 2.0;
    double smoothing_credible_interval_threshold_allele_fraction = 2.0;
    bool smoothing_requested = false;
    bool probabilistic_requested = false;
    bool window_sizes_overridden = false;
    int threads = 1;
};

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

struct Point {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    double value = 0.0;
    // MultisampleMultidimensionalKernelSegmenter does not segment on minor
    // allele fraction.  In the combined copy-ratio/allelic-count mode GATK
    // takes the *first* heterozygous site in each copy-ratio interval and
    // uses its oriented ALT fraction; intervals without a site are imputed
    // to 0.5.  Keep that segmentation signal separate from the folded MAF
    // consumed by AlleleFractionModeller below.
    double segmentation_alt_fraction = std::numeric_limits<double>::quiet_NaN();
    double minor_allele_fraction = std::numeric_limits<double>::quiet_NaN();
    // Preserve the evidence consumed by AlleleFractionModeller.  The legacy
    // table path only needed the derived MAF, while the conditional model
    // also uses total/minor read counts when available.
    std::int64_t allele_ref_count = 0;
    std::int64_t allele_alt_count = 0;
    bool allele_counts_available = false;
    // In multisample segmentation GATK evaluates the kernel on the complete
    // sample vector (a sum of per-sample kernels).  Keep that vector on the
    // Host boundary instead of collapsing it to a mean, which is not
    // equivalent for linear kernels (and can erase anti-correlated events).
    std::vector<double> copy_ratio_samples;
    std::vector<double> segmentation_alt_fraction_samples;
    std::vector<double> minor_allele_fraction_samples;
};

struct Table {
    std::vector<std::string> header_lines;
    std::vector<std::string> contigs;
    std::vector<std::int64_t> contig_lengths;
    std::vector<Point> points;
    std::string sample = "UNKNOWN";
};

struct Segment {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::size_t points = 0;
    double mean = 0.0;
    double stddev = 0.0;
    std::size_t allelic_points = 0;
    double mean_minor_allele_fraction = std::numeric_limits<double>::quiet_NaN();
    double allelic_stddev = 0.0;
};

struct InputSegment {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

bool file_complete(const std::string& path);

std::uintmax_t file_bytes(const std::string& path) {
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : size;
}

struct PosteriorQuantiles {
    double copy_ratio_10 = std::numeric_limits<double>::quiet_NaN();
    double copy_ratio_50 = std::numeric_limits<double>::quiet_NaN();
    double copy_ratio_90 = std::numeric_limits<double>::quiet_NaN();
    double allele_fraction_10 = std::numeric_limits<double>::quiet_NaN();
    double allele_fraction_50 = std::numeric_limits<double>::quiet_NaN();
    double allele_fraction_90 = std::numeric_limits<double>::quiet_NaN();
};

// CopyRatioModeller in GATK uses a Gibbs state with two global parameters
// (variance and outlier probability) and one mean per segment.  The native
// implementation keeps the bounded, reproducible chain used by the existing
// contract, but now performs a deterministic conditional-model prepass over
// the actual copy-ratio points.  Keeping these diagnostics explicit avoids
// presenting the old segment-summary variance as a Java MCMC posterior.
struct CopyRatioConditionalModel {
    double variance = std::numeric_limits<double>::quiet_NaN();
    double outlier_probability = std::numeric_limits<double>::quiet_NaN();
    std::size_t iterations = 0;
    std::vector<double> segment_means;
};

// AlleleFractionModeller's exact Java likelihood marginalizes a Gamma
// allelic-bias variable and samples segment minor fractions with a bounded
// Beta prior.  This native slice keeps the same observable state variables
// while using deterministic responsibility-weighted binomial conditionals.
struct AlleleFractionConditionalModel {
    double mean_bias = std::numeric_limits<double>::quiet_NaN();
    double bias_variance = std::numeric_limits<double>::quiet_NaN();
    double outlier_probability = std::numeric_limits<double>::quiet_NaN();
    std::size_t iterations = 0;
    // GATK's initializer uses the integrated alt-minor/ref-minor
    // responsibility at 0.5 before the Gibbs chain starts.  Retain the
    // deterministic Kokkos-computed starting point for telemetry and the
    // pinned oracle; it is distinct from the legacy segment MAF mean.
    std::vector<double> initial_segment_means;
    std::vector<double> segment_means;
};

Segment summarize_range(const Table& table, const std::size_t begin, const std::size_t end) {
    double sum = 0.0;
    double sum_sq = 0.0;
    double allelic_sum = 0.0;
    double allelic_sum_sq = 0.0;
    std::size_t allelic_points = 0;
    for (std::size_t i = begin; i < end; ++i) {
        sum += table.points[i].value;
        sum_sq += table.points[i].value * table.points[i].value;
        if (std::isfinite(table.points[i].minor_allele_fraction)) {
            allelic_sum += table.points[i].minor_allele_fraction;
            allelic_sum_sq += table.points[i].minor_allele_fraction * table.points[i].minor_allele_fraction;
            ++allelic_points;
        }
    }
    const auto count = end - begin;
    const double mean = sum / static_cast<double>(count);
    const double variance = std::max(0.0, sum_sq / static_cast<double>(count) - mean * mean);
    const double allelic_mean = allelic_points == 0 ? std::numeric_limits<double>::quiet_NaN() :
        allelic_sum / static_cast<double>(allelic_points);
    const double allelic_variance = allelic_points == 0 ? 0.0 : std::max(0.0,
        allelic_sum_sq / static_cast<double>(allelic_points) - allelic_mean * allelic_mean);
    return Segment{table.points[begin].contig, table.points[begin].start,
                   table.points[end - 1].end, count, mean, std::sqrt(variance),
                   allelic_points, allelic_mean, std::sqrt(allelic_variance)};
}

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
            std::cout << "fastgatk-model-segments (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                    denoised copy-ratio TSV\n"
                         "      --denoised-copy-ratios FILE    alias for input; repeatable for bounded multisample segmentation\n"
                         "      --allelic-counts FILE           optional allelic-count TSV validation; repeatable per sample\n"
                         "      --normal-allelic-counts FILE    optional normal allelic-count TSV\n"
                         "      --minimum-total-allele-count-case N    default 0\n"
                         "      --minimum-total-allele-count-normal N  default 30\n"
                         "      --genotyping-homozygous-log-ratio-threshold F  default -10\n"
                         "      --genotyping-base-error-rate F       default 0.05\n"
                         "  -O, --output PATH                   output TSV, or output directory with --output-prefix\n"
                         "      --output-prefix PREFIX          derives PREFIX.modelFinal.segments.tsv\n"
                         "      --segments FILE                 Picard interval-list segmentation input\n"
                         "      --change-point-threshold F      default 0.15 log2 units\n"
                         "      --allelic-change-point-threshold F  default 0.05 minor-allele-fraction units\n"
                         "      --min-points N                  default 1\n"
                         "      --maximum-number-of-segments-per-chromosome N  kernel cap (default 1000)\n"
                         "      --kernel-variance-copy-ratio F  Gaussian variance; 0 selects linear kernel\n"
                         "      --kernel-variance-allele-fraction F  Gaussian variance (default 0.025)\n"
                         "      --kernel-scaling-allele-fraction F  AF kernel scale (default 1.0)\n"
                         "      --kernel-approximation-dimension N  deterministic anchor count (default 100)\n"
                         "      --window-size N                 repeatable kernel local window\n"
                         "      --number-of-changepoints-penalty-factor F  kernel penalty (default 1.0)\n"
                         "      --mode NAME                     BOTH, COPY_RATIO, or ALLELE_FRACTION\n"
                         "      --num-burn-in-iterations N      deterministic posterior burn-in\n"
                         "      --num-samples N                 deterministic posterior samples\n"
                         "      --minor-allele-fraction-prior-alpha F  GATK prior (default 25)\n"
                         "      --number-of-samples-copy-ratio N      GATK copy-ratio MCMC samples\n"
                         "      --number-of-burn-in-samples-copy-ratio N  GATK copy-ratio burn-in\n"
                         "      --number-of-samples-allele-fraction N    GATK allele-fraction MCMC samples\n"
                         "      --number-of-burn-in-samples-allele-fraction N  GATK allele-fraction burn-in\n"
                         "      --maximum-number-of-smoothing-iterations N  GATK smoothing bound\n"
                         "      --number-of-smoothing-iterations-per-fit N    GATK smoothing fit cadence\n"
                         "      --smoothing-credible-interval-threshold-copy-ratio F  GATK smoothing threshold\n"
                         "      --smoothing-credible-interval-threshold-allele-fraction F  GATK smoothing threshold\n"
                         "      --output-manifest FILE          OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input") ||
                   has_option(argument, "--denoised-copy-ratios")) {
            const char* name = argument.rfind("--denoised-copy-ratios", 0) == 0
                ? "--denoised-copy-ratios" : "--input";
            const auto value = require_value(index, argc, argv, argument, name,
                                             argument == "-I" ? "-I" : nullptr);
            if (options.input.empty()) options.input = value;
            options.denoised_copy_ratios.push_back(value);
        } else if (has_option(argument, "--allelic-counts")) {
            const auto value = require_value(index, argc, argv, argument, "--allelic-counts");
            if (options.allelic_counts.empty()) options.allelic_counts = value;
            options.allelic_count_files.push_back(value);
        } else if (has_option(argument, "--normal-allelic-counts")) {
            options.normal_allelic_counts = require_value(index, argc, argv, argument, "--normal-allelic-counts");
        } else if (has_option(argument, "--minimum-total-allele-count-case")) {
            options.min_total_allele_count_case = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--minimum-total-allele-count-case")));
        } else if (has_option(argument, "--minimum-total-allele-count-normal")) {
            options.min_total_allele_count_normal = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--minimum-total-allele-count-normal")));
        } else if (has_option(argument, "--genotyping-homozygous-log-ratio-threshold")) {
            options.genotyping_homozygous_log_ratio_threshold = std::stod(require_value(
                index, argc, argv, argument, "--genotyping-homozygous-log-ratio-threshold"));
        } else if (has_option(argument, "--genotyping-base-error-rate")) {
            options.genotyping_base_error_rate = std::stod(require_value(
                index, argc, argv, argument, "--genotyping-base-error-rate"));
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output_directory = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--output-prefix")) {
            options.output_prefix = require_value(index, argc, argv, argument, "--output-prefix");
        } else if (has_option(argument, "--segments")) {
            options.input_segments = require_value(index, argc, argv, argument, "--segments");
        } else if (has_option(argument, "--change-point-threshold")) {
            options.change_point_threshold = std::stod(require_value(index, argc, argv, argument,
                                                                       "--change-point-threshold"));
            options.change_point_threshold_overridden = true;
        } else if (has_option(argument, "--allelic-change-point-threshold")) {
            options.allelic_change_point_threshold = std::stod(require_value(
                index, argc, argv, argument, "--allelic-change-point-threshold"));
        } else if (has_option(argument, "--min-points")) {
            options.min_points = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--min-points")));
        } else if (has_option(argument, "--maximum-number-of-segments-per-chromosome") ||
                   has_option(argument, "--maximum-number-of-segments")) {
            const char* name = argument.rfind("--maximum-number-of-segments-per-chromosome", 0) == 0
                ? "--maximum-number-of-segments-per-chromosome" : "--maximum-number-of-segments";
            options.max_segments_per_contig = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, name)));
            options.kernel_requested = true;
        } else if (has_option(argument, "--kernel-variance-copy-ratio")) {
            options.kernel_variance_copy_ratio = std::stod(require_value(
                index, argc, argv, argument, "--kernel-variance-copy-ratio"));
            options.kernel_requested = true;
        } else if (has_option(argument, "--kernel-variance-allele-fraction")) {
            options.kernel_variance_allele_fraction = std::stod(require_value(
                index, argc, argv, argument, "--kernel-variance-allele-fraction"));
            options.kernel_requested = true;
        } else if (has_option(argument, "--kernel-scaling-allele-fraction")) {
            options.kernel_scaling_allele_fraction = std::stod(require_value(
                index, argc, argv, argument, "--kernel-scaling-allele-fraction"));
            options.kernel_requested = true;
        } else if (has_option(argument, "--kernel-approximation-dimension")) {
            options.kernel_approximation_dimension = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, "--kernel-approximation-dimension")));
            options.kernel_requested = true;
        } else if (has_option(argument, "--window-size")) {
            if (!options.window_sizes_overridden) {
                options.window_sizes.clear();
                options.window_sizes_overridden = true;
            }
            options.window_sizes.push_back(static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, "--window-size"))));
            options.kernel_requested = true;
        } else if (has_option(argument, "--number-of-changepoints-penalty-factor")) {
            options.changepoint_penalty_factor = std::stod(require_value(
                index, argc, argv, argument, "--number-of-changepoints-penalty-factor"));
            options.kernel_requested = true;
        } else if (has_option(argument, "--mode")) {
            options.mode = require_value(index, argc, argv, argument, "--mode");
            std::transform(options.mode.begin(), options.mode.end(), options.mode.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--num-burn-in-iterations")) {
            options.num_burn_in_iterations = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, "--num-burn-in-iterations")));
            options.generic_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--num-samples")) {
            options.num_samples = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, "--num-samples")));
            options.generic_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--minor-allele-fraction-prior-alpha")) {
            options.minor_allele_fraction_prior_alpha = std::stod(require_value(
                index, argc, argv, argument, "--minor-allele-fraction-prior-alpha"));
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--number-of-samples-copy-ratio")) {
            options.num_samples_copy_ratio = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--number-of-samples-copy-ratio")));
            options.gatk_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--number-of-burn-in-samples-copy-ratio")) {
            options.num_burn_in_samples_copy_ratio = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--number-of-burn-in-samples-copy-ratio")));
            options.gatk_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--number-of-samples-allele-fraction")) {
            options.num_samples_allele_fraction = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--number-of-samples-allele-fraction")));
            options.gatk_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--number-of-burn-in-samples-allele-fraction")) {
            options.num_burn_in_samples_allele_fraction = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--number-of-burn-in-samples-allele-fraction")));
            options.gatk_mcmc_overridden = true;
            options.probabilistic_requested = true;
        } else if (has_option(argument, "--maximum-number-of-smoothing-iterations")) {
            options.maximum_number_of_smoothing_iterations = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--maximum-number-of-smoothing-iterations")));
            options.smoothing_requested = true;
        } else if (has_option(argument, "--number-of-smoothing-iterations-per-fit")) {
            options.number_of_smoothing_iterations_per_fit = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--number-of-smoothing-iterations-per-fit")));
            options.smoothing_requested = true;
        } else if (has_option(argument, "--smoothing-credible-interval-threshold-copy-ratio")) {
            options.smoothing_credible_interval_threshold_copy_ratio = std::stod(require_value(
                index, argc, argv, argument, "--smoothing-credible-interval-threshold-copy-ratio"));
            options.smoothing_requested = true;
        } else if (has_option(argument, "--smoothing-credible-interval-threshold-allele-fraction")) {
            options.smoothing_credible_interval_threshold_allele_fraction = std::stod(require_value(
                index, argc, argv, argument, "--smoothing-credible-interval-threshold-allele-fraction"));
            options.smoothing_requested = true;
        } else if (has_option(argument, "--sample")) {
            options.sample = require_value(index, argc, argv, argument, "--sample");
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--disable-tool-default-read-filters") {
            // Compatibility switches for a materialized copy-number table.
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
    if (options.input.empty()) throw std::invalid_argument("-I/--input or --denoised-copy-ratios is required");
    if (options.output_prefix.empty() && options.output_directory.empty()) {
        throw std::invalid_argument("--denoised-copy-ratios and --output-prefix are required");
    }
    if (!options.output_prefix.empty()) {
        const auto prefix = std::filesystem::path(options.output_prefix);
        const auto leaf = prefix.filename().empty() ? prefix : prefix.filename();
        if (options.output_directory.empty()) {
            options.output = options.output_prefix + ".modelFinal.segments.tsv";
        } else {
            std::error_code error;
            std::filesystem::create_directories(options.output_directory, error);
            if (error) throw std::invalid_argument("cannot create ModelSegments output directory: " +
                                                   options.output_directory);
            options.output = (std::filesystem::path(options.output_directory) /
                              (leaf.string() + ".modelFinal.segments.tsv")).string();
        }
    } else {
        options.output = options.output_directory;
    }
    if (!std::isfinite(options.change_point_threshold) || options.change_point_threshold < 0.0 ||
        !std::isfinite(options.allelic_change_point_threshold) || options.allelic_change_point_threshold < 0.0 ||
        options.min_points == 0 || options.max_segments_per_contig == 0 ||
        !std::isfinite(options.kernel_variance_copy_ratio) || options.kernel_variance_copy_ratio < 0.0 ||
        !std::isfinite(options.kernel_variance_allele_fraction) || options.kernel_variance_allele_fraction < 0.0 ||
        !std::isfinite(options.kernel_scaling_allele_fraction) || options.kernel_scaling_allele_fraction < 0.0 ||
        options.kernel_approximation_dimension == 0 || options.window_sizes.empty() ||
        !std::isfinite(options.changepoint_penalty_factor) || options.changepoint_penalty_factor < 0.0 ||
        options.threads < 1 ||
        !std::isfinite(options.genotyping_homozygous_log_ratio_threshold) ||
        !std::isfinite(options.genotyping_base_error_rate) ||
        options.genotyping_base_error_rate < 0.0 || options.genotyping_base_error_rate > 1.0)
        throw std::invalid_argument("invalid ModelSegments prototype parameters");
    if (options.mode != "BOTH" && options.mode != "COPY_RATIO" &&
        options.mode != "ALLELE_FRACTION")
        throw std::invalid_argument("--mode must be BOTH, COPY_RATIO, or ALLELE_FRACTION");
    if (options.num_samples == 0 || options.num_samples > 1000000 ||
        options.num_burn_in_iterations > 1000000)
        throw std::invalid_argument("--num-samples must be in [1,1000000] and --num-burn-in-iterations in [0,1000000]");
    if (options.gatk_mcmc_overridden && !options.generic_mcmc_overridden) {
        options.num_samples = std::max(options.num_samples_copy_ratio, options.num_samples_allele_fraction);
        options.num_burn_in_iterations = std::max(options.num_burn_in_samples_copy_ratio,
                                                  options.num_burn_in_samples_allele_fraction);
    }
    if (!std::isfinite(options.minor_allele_fraction_prior_alpha) ||
        options.minor_allele_fraction_prior_alpha < 1.0 ||
        !std::isfinite(options.smoothing_credible_interval_threshold_copy_ratio) ||
        options.smoothing_credible_interval_threshold_copy_ratio < 0.0 ||
        !std::isfinite(options.smoothing_credible_interval_threshold_allele_fraction) ||
        options.smoothing_credible_interval_threshold_allele_fraction < 0.0 ||
        (options.gatk_mcmc_overridden &&
         (options.num_samples_copy_ratio == 0 || options.num_samples_allele_fraction == 0 ||
          options.num_burn_in_samples_copy_ratio >= options.num_samples_copy_ratio ||
          options.num_burn_in_samples_allele_fraction >= options.num_samples_allele_fraction)))
        throw std::invalid_argument("invalid GATK ModelSegments modeling parameters");
    if (options.num_samples == 0 || options.num_samples > 1000000 ||
        options.num_burn_in_iterations > 1000000)
        throw std::invalid_argument("--num-samples must be in [1,1000000] and --num-burn-in-iterations in [0,1000000]");
    if (options.num_burn_in_iterations >= options.num_samples)
        throw std::invalid_argument("number of MCMC samples must be greater than burn-in samples");
    for (const auto window : options.window_sizes)
        if (window == 0) throw std::invalid_argument("window sizes must be positive");
    std::sort(options.window_sizes.begin(), options.window_sizes.end());
    options.window_sizes.erase(std::unique(options.window_sizes.begin(), options.window_sizes.end()),
                               options.window_sizes.end());
    // In GATK ModelSegments the segmentation stage is KernelSegmenter by
    // default; --segments is the only normal way to skip it.  Native's
    // --change-point-threshold remains an explicit opt-in legacy extension so
    // existing threshold fixtures remain reproducible.
    if (!options.change_point_threshold_overridden)
        options.kernel_requested = true;
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

Table read_copy_ratios(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open denoised copy-ratio table: " + path);
    Table table;
    bool columns_seen = false;
    int contig_col = -1;
    int start_col = -1;
    int end_col = -1;
    int value_col = -1;
    for (std::string line; std::getline(input, line);) {
        if (line.empty()) continue;
        if (line.rfind("@RG", 0) == 0) {
            for (const auto& field : split_tab(line))
                if (field.rfind("SM:", 0) == 0) table.sample = field.substr(3);
            table.header_lines.push_back(line);
            continue;
        }
        if (line.rfind("@SQ", 0) == 0) {
            table.header_lines.push_back(line);
            const auto fields = split_tab(line);
            std::string contig;
            std::int64_t length = 0;
            for (const auto& field : fields) {
                if (field.rfind("SN:", 0) == 0) contig = field.substr(3);
                else if (field.rfind("LN:", 0) == 0) length = std::stoll(field.substr(3));
            }
            if (!contig.empty()) {
                table.contigs.push_back(contig);
                table.contig_lengths.push_back(length);
            }
            continue;
        }
        if (line[0] == '@') {
            table.header_lines.push_back(line);
            continue;
        }
        const auto fields = split_tab(line);
        if (!columns_seen) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "CONTIG") contig_col = static_cast<int>(i);
                else if (fields[i] == "START") start_col = static_cast<int>(i);
                else if (fields[i] == "END") end_col = static_cast<int>(i);
                else if (fields[i] == "LOG2_COPY_RATIO" || fields[i] == "MEAN_LOG2_COPY_RATIO") value_col = static_cast<int>(i);
            }
            if (contig_col < 0 || start_col < 0 || end_col < 0 || value_col < 0)
                throw std::runtime_error("BAD_INPUT: copy-ratio table requires CONTIG/START/END/LOG2_COPY_RATIO");
            columns_seen = true;
            continue;
        }
        const auto max_col = std::max({contig_col, start_col, end_col, value_col});
        if (fields.size() <= static_cast<std::size_t>(max_col))
            throw std::runtime_error("BAD_INPUT: malformed copy-ratio row");
        Point point;
        point.contig = fields[contig_col];
        point.start = std::stoll(fields[start_col]);
        point.end = std::stoll(fields[end_col]);
        point.value = std::stod(fields[value_col]);
        if (point.contig.empty() || point.start < 1 || point.end < point.start || !std::isfinite(point.value))
            throw std::runtime_error("BAD_INPUT: invalid copy-ratio row");
        table.points.push_back(std::move(point));
    }
    if (!columns_seen || table.points.empty()) throw std::runtime_error("BAD_INPUT: copy-ratio table has no records");
    return table;
}

std::vector<InputSegment> read_input_segments(const std::string& path, const Table& table) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open segments interval list: " + path);
    std::vector<InputSegment> segments;
    std::unordered_set<std::string> known_contigs(table.contigs.begin(), table.contigs.end());
    for (std::string line; std::getline(input, line);) {
        if (line.empty() || line[0] == '@') continue;
        const auto fields = split_tab(line);
        if (fields.size() < 3)
            throw std::runtime_error("BAD_INPUT: segments interval list row requires CONTIG/START/END");
        InputSegment segment;
        segment.contig = fields[0];
        try {
            segment.start = std::stoll(fields[1]);
            segment.end = std::stoll(fields[2]);
        } catch (const std::exception&) {
            throw std::runtime_error("BAD_INPUT: malformed segments interval list row: " + line);
        }
        if (segment.contig.empty() || segment.start < 1 || segment.end < segment.start)
            throw std::runtime_error("BAD_INPUT: invalid segments interval: " + line);
        if (!known_contigs.empty() && known_contigs.find(segment.contig) == known_contigs.end())
            throw std::runtime_error("BAD_INPUT: segments interval uses contig absent from copy-ratio dictionary: " + segment.contig);
        const auto contig_index = std::find(table.contigs.begin(), table.contigs.end(), segment.contig);
        if (contig_index != table.contigs.end()) {
            const auto index = static_cast<std::size_t>(contig_index - table.contigs.begin());
            if (index < table.contig_lengths.size() && table.contig_lengths[index] > 0 &&
                segment.end > table.contig_lengths[index])
                throw std::runtime_error("BAD_INPUT: segments interval exceeds contig length: " + line);
        }
        segments.push_back(std::move(segment));
    }
    if (segments.empty()) throw std::runtime_error("BAD_INPUT: segments interval list has no intervals");
    // Picard IntervalList is sequence-dictionary ordered and non-overlapping for
    // ModelSegments.  Validate rather than silently assigning a copy-ratio point
    // to multiple supplied segments, which would diverge from GATK's model.
    const auto contig_rank = [&](const std::string& contig) {
        const auto iterator = std::find(table.contigs.begin(), table.contigs.end(), contig);
        return iterator == table.contigs.end()
            ? table.contigs.size() : static_cast<std::size_t>(iterator - table.contigs.begin());
    };
    for (std::size_t index = 1; index < segments.size(); ++index) {
        const auto& previous = segments[index - 1];
        const auto& current = segments[index];
        const auto previous_rank = contig_rank(previous.contig);
        const auto current_rank = contig_rank(current.contig);
        if (current_rank < previous_rank ||
            (current_rank == previous_rank &&
             (current.start < previous.start || current.start <= previous.end)))
            throw std::runtime_error("BAD_INPUT: segments interval list must be dictionary-ordered and non-overlapping");
    }
    return segments;
}

std::vector<Segment> summarize_input_segments(const Table& table,
                                              const std::vector<InputSegment>& input_segments) {
    std::vector<Segment> segments;
    segments.reserve(input_segments.size());
    for (const auto& interval : input_segments) {
        std::vector<std::size_t> point_indices;
        for (std::size_t index = 0; index < table.points.size(); ++index) {
            const auto& point = table.points[index];
            if (point.contig != interval.contig) continue;
            // CopyRatioSegmentedData in GATK assigns a copy-ratio interval to a
            // segment by its midpoint, not by any-overlap.  Keep that detail at
            // the Host boundary so supplied --segments partitions are faithful.
            const auto midpoint = point.start + (point.end - point.start) / 2;
            if (midpoint >= interval.start && midpoint <= interval.end)
                point_indices.push_back(index);
        }
        Segment summary{interval.contig, interval.start, interval.end, point_indices.size(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(), 0,
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN()};
        if (!point_indices.empty()) {
            double sum = 0.0;
            double sum_sq = 0.0;
            double allelic_sum = 0.0;
            double allelic_sum_sq = 0.0;
            for (const auto index : point_indices) {
                const auto value = table.points[index].value;
                sum += value;
                sum_sq += value * value;
                if (std::isfinite(table.points[index].minor_allele_fraction)) {
                    allelic_sum += table.points[index].minor_allele_fraction;
                    allelic_sum_sq += table.points[index].minor_allele_fraction * table.points[index].minor_allele_fraction;
                    ++summary.allelic_points;
                }
            }
            summary.mean = sum / static_cast<double>(point_indices.size());
            summary.stddev = std::sqrt(std::max(0.0,
                sum_sq / static_cast<double>(point_indices.size()) - summary.mean * summary.mean));
            if (summary.allelic_points > 0) {
                summary.mean_minor_allele_fraction = allelic_sum / static_cast<double>(summary.allelic_points);
                summary.allelic_stddev = std::sqrt(std::max(0.0,
                    allelic_sum_sq / static_cast<double>(summary.allelic_points) -
                    summary.mean_minor_allele_fraction * summary.mean_minor_allele_fraction));
            }
        }
        segments.push_back(std::move(summary));
    }
    return segments;
}

struct AllelicCountRecord {
    std::int64_t ref = 0;
    std::int64_t alt = 0;
    std::string ref_nucleotide = "N";
    std::string alt_nucleotide = "N";
};

using AllelicCountMap = std::unordered_map<std::string, AllelicCountRecord>;

AllelicCountMap read_allelic_counts(const std::string& path, std::string* sample_name = nullptr) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open allelic-count table: " + path);
    AllelicCountMap counts;
    int contig_col = -1;
    int position_col = -1;
    int ref_count_col = -1;
    int alt_count_col = -1;
    int ref_nucleotide_col = -1;
    int alt_nucleotide_col = -1;
    bool header_seen = false;
    for (std::string line; std::getline(input, line);) {
        if (line.empty()) continue;
        if (line.rfind("@RG", 0) == 0) {
            if (sample_name != nullptr) {
                for (const auto& field : split_tab(line))
                    if (field.rfind("SM:", 0) == 0) *sample_name = field.substr(3);
            }
            continue;
        }
        if (line[0] == '@') continue;
        const auto fields = split_tab(line);
        if (!header_seen) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "CONTIG") contig_col = static_cast<int>(i);
                if (fields[i] == "POSITION") position_col = static_cast<int>(i);
                if (fields[i] == "REF_COUNT") ref_count_col = static_cast<int>(i);
                if (fields[i] == "ALT_COUNT") alt_count_col = static_cast<int>(i);
                if (fields[i] == "REF_NUCLEOTIDE") ref_nucleotide_col = static_cast<int>(i);
                if (fields[i] == "ALT_NUCLEOTIDE") alt_nucleotide_col = static_cast<int>(i);
            }
            if (contig_col < 0 || position_col < 0 || ref_count_col < 0 || alt_count_col < 0)
                throw std::runtime_error("BAD_INPUT: allelic-count header requires CONTIG/POSITION/REF_COUNT/ALT_COUNT");
            header_seen = true;
            continue;
        }
        if (fields.size() <= static_cast<std::size_t>(std::max({contig_col, position_col, ref_count_col, alt_count_col})))
            throw std::runtime_error("BAD_INPUT: malformed allelic-count row");
        const auto ref = std::stoll(fields[ref_count_col]);
        const auto alt = std::stoll(fields[alt_count_col]);
        if (ref < 0 || alt < 0)
            throw std::runtime_error("BAD_INPUT: allelic counts must be non-negative");
        const auto key = fields[contig_col] + ":" + fields[position_col];
        AllelicCountRecord record;
        record.ref = ref;
        record.alt = alt;
        if (ref_nucleotide_col >= 0 && static_cast<std::size_t>(ref_nucleotide_col) < fields.size())
            record.ref_nucleotide = fields[ref_nucleotide_col];
        if (alt_nucleotide_col >= 0 && static_cast<std::size_t>(alt_nucleotide_col) < fields.size())
            record.alt_nucleotide = fields[alt_nucleotide_col];
        const auto inserted = counts.emplace(key, std::move(record));
        if (!inserted.second && (inserted.first->second.ref != ref || inserted.first->second.alt != alt))
            throw std::runtime_error("BAD_INPUT: duplicate allelic-count locus with conflicting counts: " + key);
    }
    if (!header_seen) throw std::runtime_error("BAD_INPUT: allelic-count table has no header");
    if (counts.empty()) throw std::runtime_error("BAD_INPUT: allelic-count table has no records");
    return counts;
}

KOKKOS_INLINE_FUNCTION double model_segments_beta_continued_fraction(
    const double aa, const double bb, const double xx) {
    constexpr int max_iterations = 10000;
    constexpr double epsilon = 3.0e-14;
    constexpr double tiny = 1.0e-300;
    const double qab = aa + bb;
    const double qap = aa + 1.0;
    const double qam = aa - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * xx / qap;
    if (Kokkos::abs(d) < tiny) d = tiny;
    d = 1.0 / d;
    double h = d;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        const double m = static_cast<double>(iteration);
        const double m2 = 2.0 * m;
        double term = m * (bb - m) * xx / ((qam + m2) * (aa + m2));
        d = 1.0 + term * d;
        if (Kokkos::abs(d) < tiny) d = tiny;
        c = 1.0 + term / c;
        if (Kokkos::abs(c) < tiny) c = tiny;
        d = 1.0 / d;
        h *= d * c;
        term = -(aa + m) * (qab + m) * xx / ((aa + m2) * (qap + m2));
        d = 1.0 + term * d;
        if (Kokkos::abs(d) < tiny) d = tiny;
        c = 1.0 + term / c;
        if (Kokkos::abs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (Kokkos::abs(delta - 1.0) < epsilon) break;
    }
    return h;
}

KOKKOS_INLINE_FUNCTION double model_segments_regularized_beta(
    const double x, const double a, const double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double log_bt = Kokkos::lgamma(a + b) - Kokkos::lgamma(a) - Kokkos::lgamma(b) +
        a * Kokkos::log(x) + b * Kokkos::log1p(-x);
    const double bt = Kokkos::exp(log_bt);
    const double value = x < (a + 1.0) / (a + b + 2.0)
        ? bt * model_segments_beta_continued_fraction(a, b, x) / a
        : 1.0 - bt * model_segments_beta_continued_fraction(b, a, 1.0 - x) / b;
    return Kokkos::max(0.0, Kokkos::min(1.0, value));
}

KOKKOS_INLINE_FUNCTION double model_segments_homozygous_log_ratio(
    const std::int64_t ref, const std::int64_t alt, const double base_error_rate) {
    const double r = static_cast<double>(ref);
    const double n = static_cast<double>(ref + alt);
    const double beta_all = model_segments_regularized_beta(1.0, r + 1.0, n - r + 1.0);
    const double beta_error = model_segments_regularized_beta(base_error_rate, r + 1.0, n - r + 1.0);
    const double beta_one_minus_error = model_segments_regularized_beta(
        1.0 - base_error_rate, r + 1.0, n - r + 1.0);
    const double beta_hom = beta_error + beta_all - beta_one_minus_error;
    const double beta_het = beta_one_minus_error - beta_error;
    if (!(beta_hom > 0.0) || !(beta_het > 0.0))
        return beta_hom > beta_het ? 1.0e300 : -1.0e300;
    return Kokkos::log(beta_hom) - Kokkos::log(beta_het);
}

AllelicCountMap filter_heterozygous_counts(const AllelicCountMap& counts,
                                           const Table* copy_ratios,
                                           const std::size_t minimum_total_count,
                                           const double homozygous_log_ratio_threshold,
                                           const double base_error_rate,
                                           KernelStats* kernel_stats = nullptr) {
    AllelicCountMap result;
    std::unordered_map<std::string, std::vector<std::pair<std::int64_t, std::int64_t>>> intervals;
    if (copy_ratios != nullptr) {
        for (const auto& point : copy_ratios->points)
            intervals[point.contig].emplace_back(point.start, point.end);
        for (auto& entry : intervals)
            std::sort(entry.second.begin(), entry.second.end());
    }
    const auto is_in_interval = [&](const std::string& contig, const std::int64_t position) {
        if (copy_ratios == nullptr) return true;
        const auto iterator = intervals.find(contig);
        if (iterator == intervals.end()) return false;
        const auto& values = iterator->second;
        const auto upper = std::upper_bound(
            values.begin(), values.end(), position,
            [](const std::int64_t value, const auto& interval) { return value < interval.first; });
        if (upper == values.begin()) return false;
        const auto& candidate = *(upper - 1);
        return candidate.first <= position && position <= candidate.second;
    };
    std::vector<std::pair<std::string, AllelicCountRecord>> candidates;
    candidates.reserve(counts.size());
    for (const auto& entry : counts) {
        const auto& record = entry.second;
        if (record.ref + record.alt < static_cast<std::int64_t>(minimum_total_count)) continue;
        const auto separator = entry.first.find(':');
        const auto contig = separator == std::string::npos ? entry.first : entry.first.substr(0, separator);
        const auto position = separator == std::string::npos ? 0 : std::stoll(entry.first.substr(separator + 1));
        if (!is_in_interval(contig, position)) continue;
        candidates.emplace_back(entry.first, record);
    }
    if (candidates.empty()) return result;
    Kokkos::View<std::int64_t*> refs("model_segments_heterozygous_refs", candidates.size());
    Kokkos::View<std::int64_t*> alts("model_segments_heterozygous_alts", candidates.size());
    Kokkos::View<std::uint8_t*> retained("model_segments_heterozygous_retained", candidates.size());
    auto host_refs = Kokkos::create_mirror_view(refs);
    auto host_alts = Kokkos::create_mirror_view(alts);
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        host_refs(index) = candidates[index].second.ref;
        host_alts(index) = candidates[index].second.alt;
    }
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("model-segments-heterozygous-v1");
    host_batch.records = candidates.size();
    host_batch.bytes = candidates.size() * (2 * sizeof(std::int64_t) + sizeof(std::uint8_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-heterozygous");
    plan.begin_prepare(host_batch);
    Kokkos::deep_copy(refs, host_refs);
    Kokkos::deep_copy(alts, host_alts);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(candidates.size());
    device_batch.bind("refs", refs);
    device_batch.bind("alts", alts);
    device_batch.bind("retained", retained);
    ExecSpace().fence();
    plan.end_prepare(device_batch);
    plan.begin_execute();
    Kokkos::parallel_for("model_segments_heterozygous_genotyping",
                         Kokkos::RangePolicy<ExecSpace>(0, candidates.size()),
                         KOKKOS_LAMBDA(const std::size_t index) {
                             const auto ratio = model_segments_homozygous_log_ratio(
                                 refs(index), alts(index), base_error_rate);
                             retained(index) = ratio < homozygous_log_ratio_threshold ? 1 : 0;
                         });
    ExecSpace().fence();
    plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += candidates.size();
        kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
    }
    auto host_retained = Kokkos::create_mirror_view(retained);
    Kokkos::deep_copy(host_retained, retained);
    for (std::size_t index = 0; index < candidates.size(); ++index)
        if (host_retained(index) != 0) result.emplace(candidates[index].first, candidates[index].second);
    return result;
}

AllelicCountMap intersect_counts(const AllelicCountMap& counts, const AllelicCountMap& allowed) {
    AllelicCountMap result;
    for (const auto& entry : counts)
        if (allowed.find(entry.first) != allowed.end()) result.emplace(entry.first, entry.second);
    return result;
}

void write_allelic_count_table(const std::string& path, const Table& table,
                               const AllelicCountMap& counts, const std::string& sample) {
    std::vector<std::pair<std::string, AllelicCountRecord>> ordered(counts.begin(), counts.end());
    std::sort(ordered.begin(), ordered.end(), [](const auto& left, const auto& right) {
        const auto left_separator = left.first.find(':');
        const auto right_separator = right.first.find(':');
        const auto left_contig = left.first.substr(0, left_separator);
        const auto right_contig = right.first.substr(0, right_separator);
        if (left_contig != right_contig) return left_contig < right_contig;
        return std::stoll(left.first.substr(left_separator + 1)) < std::stoll(right.first.substr(right_separator + 1));
    });
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write heterozygous allelic counts: " + path);
    output << "@HD\tVN:1.6\n";
    for (std::size_t index = 0; index < table.contigs.size(); ++index)
        output << "@SQ\tSN:" << table.contigs[index] << "\tLN:" << table.contig_lengths[index] << '\n';
    output << "@RG\tID:GATKCopyNumber\tSM:" << (sample.empty() ? "UNKNOWN" : sample) << '\n';
    output << "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n";
    for (const auto& entry : ordered) {
        const auto separator = entry.first.find(':');
        output << entry.first.substr(0, separator) << '\t'
               << entry.first.substr(separator + 1) << '\t'
               << entry.second.ref << '\t' << entry.second.alt << '\t'
               << entry.second.ref_nucleotide << '\t' << entry.second.alt_nucleotide << '\n';
    }
    output.close();
    if (!output || !file_complete(path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete heterozygous allelic counts: " + path);
}

std::filesystem::path model_segments_prefix_path(const Options& options) {
    if (options.output_prefix.empty()) return {};
    return options.output_directory.empty()
        ? std::filesystem::path(options.output_prefix)
        : std::filesystem::path(options.output_directory) /
          std::filesystem::path(options.output_prefix).filename();
}

void apply_allelic_counts_to_table(Table& table, const AllelicCountMap& counts) {
    struct MAFObservation {
        std::int64_t position;
        double maf;
        double alt_fraction;
        std::int64_t ref;
        std::int64_t alt;
    };
    std::unordered_map<std::string, std::vector<MAFObservation>> by_contig;
    by_contig.reserve(counts.size());
    for (const auto& entry : counts) {
        const auto separator = entry.first.find(':');
        if (separator == std::string::npos) continue;
        const auto contig = entry.first.substr(0, separator);
        const auto position = std::stoll(entry.first.substr(separator + 1));
        const auto total = entry.second.ref + entry.second.alt;
        if (total > 0)
            by_contig[contig].push_back({position, static_cast<double>(std::min(
                entry.second.ref, entry.second.alt)) / static_cast<double>(total),
                static_cast<double>(entry.second.alt) / static_cast<double>(total),
                entry.second.ref, entry.second.alt});
    }
    for (auto& entry : by_contig)
        std::sort(entry.second.begin(), entry.second.end(),
                  [](const auto& left, const auto& right) { return left.position < right.position; });
    for (auto& point : table.points) {
        const auto iterator = by_contig.find(point.contig);
        if (iterator == by_contig.end()) continue;
        const auto& observations = iterator->second;
        const auto first = std::lower_bound(
            observations.begin(), observations.end(), point.start,
            [](const MAFObservation& observation, const std::int64_t position) {
                return observation.position < position;
            });
        const auto last = std::upper_bound(
            observations.begin(), observations.end(), point.end,
            [](const std::int64_t position, const MAFObservation& observation) {
                return position < observation.position;
            });
        double sum = 0.0;
        std::int64_t ref_sum = 0;
        std::int64_t alt_sum = 0;
        const auto observation_count = static_cast<std::size_t>(last - first);
        for (auto cursor = first; cursor != last; ++cursor) {
            sum += cursor->maf;
            ref_sum += cursor->ref;
            alt_sum += cursor->alt;
        }
        const auto observations_count = observation_count;
        if (observations_count > 0) {
            // GATK's combined-mode segmentation uses the first site in
            // dictionary/coordinate order, not an average of all sites in
            // the copy-ratio interval.  Modeling still consumes the complete
            // count evidence through the folded MAF/count summaries.
            point.segmentation_alt_fraction = first->alt_fraction;
            point.minor_allele_fraction = sum / static_cast<double>(observations_count);
            point.allele_ref_count = ref_sum;
            point.allele_alt_count = alt_sum;
            point.allele_counts_available = ref_sum + alt_sum > 0;
        }
    }
}

AllelicCountMap prepare_case_heterozygous_counts(
    const Table& copy_ratios, const AllelicCountMap& raw_counts,
    const AllelicCountMap* matched_normal_hets, const Options& options,
    KernelStats* kernel_stats = nullptr) {
    auto filtered = filter_heterozygous_counts(
        raw_counts, &copy_ratios, options.min_total_allele_count_case,
        options.genotyping_homozygous_log_ratio_threshold,
        options.genotyping_base_error_rate, kernel_stats);
    if (matched_normal_hets != nullptr) filtered = intersect_counts(filtered, *matched_normal_hets);
    return filtered;
}

void write_heterozygous_sidecars(const Options& options, const Table& table,
                                 const AllelicCountMap& case_hets,
                                 const AllelicCountMap* normal_hets,
                                 const std::string& normal_sample) {
    if (options.output_prefix.empty()) return;
    const auto prefix = model_segments_prefix_path(options);
    const auto sample = options.sample.empty() ? table.sample : options.sample;
    write_allelic_count_table(prefix.string() + ".hets.tsv", table, case_hets, sample);
    if (normal_hets != nullptr)
        write_allelic_count_table(prefix.string() + ".hets.normal.tsv", table, *normal_hets,
                                  normal_sample.empty() ? "NORMAL" : normal_sample);
}

std::vector<Segment> segment_points(const Table& table, double threshold, std::size_t min_points,
                                    double allelic_threshold, int requested_threads,
                                    KernelStats* kernel_stats = nullptr) {
    Kokkos::View<double*> values("model_segments_values", table.points.size());
    auto host_values = Kokkos::create_mirror_view(values);
    for (std::size_t i = 0; i < table.points.size(); ++i) host_values(i) = table.points[i].value;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<double*> allelic_values("model_segments_allelic_values", table.points.size());
    Kokkos::View<std::uint8_t*> allelic_valid("model_segments_allelic_valid", table.points.size());
    auto host_allelic_values = Kokkos::create_mirror_view(allelic_values);
    auto host_allelic_valid = Kokkos::create_mirror_view(allelic_valid);
    for (std::size_t i = 0; i < table.points.size(); ++i) {
        const bool valid = std::isfinite(table.points[i].segmentation_alt_fraction);
        host_allelic_values(i) = valid ? table.points[i].segmentation_alt_fraction : 0.0;
        host_allelic_valid(i) = valid ? 1 : 0;
    }
    fastgatk::core::HostBatch host_batch("model-segments-boundary-v1");
    host_batch.records = table.points.size();
    host_batch.bytes = 2 * table.points.size() * sizeof(double) +
                       3 * table.points.size() * sizeof(std::uint8_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-boundary");
    plan.begin_prepare(host_batch);
    Kokkos::deep_copy(values, host_values);
    Kokkos::deep_copy(allelic_values, host_allelic_values);
    Kokkos::deep_copy(allelic_valid, host_allelic_valid);
    // String/coordinate comparisons stay on Host, while the numeric change
    // point predicate is evaluated by the selected Kokkos execution space.
    // This keeps the segmentation decision identical for Serial/OpenMP/CUDA/
    // HIP/SYCL builds instead of silently falling back to a host-only loop.
    Kokkos::View<std::uint8_t*> contiguous("model_segments_contiguous", table.points.size());
    auto host_contiguous = Kokkos::create_mirror_view(contiguous);
    host_contiguous(0) = 0;
    for (std::size_t i = 1; i < table.points.size(); ++i) {
        const auto& previous = table.points[i - 1];
        const auto& current = table.points[i];
        host_contiguous(i) = previous.contig == current.contig &&
            current.start <= previous.end + 1 ? 1 : 0;
    }
    Kokkos::deep_copy(contiguous, host_contiguous);
    Kokkos::View<std::uint8_t*> boundaries("model_segments_boundaries", table.points.size());
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(table.points.size());
    device_batch.bind("values", values);
    device_batch.bind("allelic_values", allelic_values);
    device_batch.bind("allelic_valid", allelic_valid);
    device_batch.bind("contiguous", contiguous);
    device_batch.bind("boundaries", boundaries);
    ExecSpace().fence();
    plan.end_prepare(device_batch);
    plan.begin_execute();
    Kokkos::parallel_for("model_segments_boundary_kernel",
                         Kokkos::RangePolicy<ExecSpace>(0, table.points.size()),
                         KOKKOS_LAMBDA(const std::size_t i) {
                             if (i == 0 || contiguous(i) == 0) {
                                 boundaries(i) = 1;
                             } else {
                                const bool copy_ratio_boundary =
                                    Kokkos::abs(values(i) - values(i - 1)) > threshold;
                                const bool allelic_boundary = allelic_valid(i) != 0 &&
                                    allelic_valid(i - 1) != 0 &&
                                    Kokkos::abs(allelic_values(i) - allelic_values(i - 1)) > allelic_threshold;
                                boundaries(i) = copy_ratio_boundary || allelic_boundary ? 1 : 0;
                             }
                         });
    ExecSpace().fence();
    plan.end_execute();
    if (kernel_stats != nullptr) {
        ++kernel_stats->batches;
        kernel_stats->observations += table.points.size();
        kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
        kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
    }
    auto host_boundaries = Kokkos::create_mirror_view(boundaries);
    Kokkos::deep_copy(host_boundaries, boundaries);
    std::vector<Segment> segments;
    std::size_t begin = 0;
    while (begin < table.points.size()) {
        std::size_t end = begin + 1;
        double sum = table.points[begin].value;
        double sum_sq = table.points[begin].value * table.points[begin].value;
        double allelic_sum = 0.0;
        double allelic_sum_sq = 0.0;
        std::size_t allelic_points = 0;
        double segmentation_allele_sum = 0.0;
        std::size_t segmentation_allele_points = 0;
        if (std::isfinite(table.points[begin].segmentation_alt_fraction)) {
            segmentation_allele_sum = table.points[begin].segmentation_alt_fraction;
            segmentation_allele_points = 1;
        }
        if (std::isfinite(table.points[begin].minor_allele_fraction)) {
            allelic_sum = table.points[begin].minor_allele_fraction;
            allelic_sum_sq = allelic_sum * allelic_sum;
            allelic_points = 1;
        }
        while (end < table.points.size()) {
            const auto& previous = table.points[end - 1];
            const auto& current = table.points[end];
            const double mean = sum / static_cast<double>(end - begin);
            const bool contiguous = previous.contig == current.contig && current.start <= previous.end + 1;
            // The device predicate is the fast-path boundary.  Retain the
            // running-mean check as a conservative Host refinement so a
            // gradual drift cannot be merged differently from GATK's
            // segment-level mean rule.
            const bool current_segmentation_allelic = std::isfinite(current.segmentation_alt_fraction);
            const double segmentation_allelic_mean = segmentation_allele_points == 0 ? 0.0 :
                segmentation_allele_sum / static_cast<double>(segmentation_allele_points);
            const bool allelic_break = current_segmentation_allelic && segmentation_allele_points > 0 &&
                std::abs(current.segmentation_alt_fraction - segmentation_allelic_mean) > allelic_threshold;
            if (host_boundaries(end) != 0 || !contiguous || std::abs(current.value - mean) > threshold || allelic_break) break;
            sum += current.value;
            sum_sq += current.value * current.value;
            if (std::isfinite(current.minor_allele_fraction)) {
                allelic_sum += current.minor_allele_fraction;
                allelic_sum_sq += current.minor_allele_fraction * current.minor_allele_fraction;
                ++allelic_points;
            }
            if (current_segmentation_allelic) {
                segmentation_allele_sum += current.segmentation_alt_fraction;
                ++segmentation_allele_points;
            }
            ++end;
        }
        if (end - begin < min_points && end < table.points.size()) {
            sum += table.points[end].value;
            sum_sq += table.points[end].value * table.points[end].value;
            if (std::isfinite(table.points[end].minor_allele_fraction)) {
                allelic_sum += table.points[end].minor_allele_fraction;
                allelic_sum_sq += table.points[end].minor_allele_fraction * table.points[end].minor_allele_fraction;
                ++allelic_points;
            }
            ++end;
        }
        segments.push_back(summarize_range(table, begin, end));
        begin = end;
    }
    (void)requested_threads;
    return segments;
}

struct KernelCandidate {
    std::size_t boundary = 0;  // first point in the right-hand segment
    double score = 0.0;        // lower is a more persistent changepoint
};

// Deterministic Kokkos implementation of the bounded kernel-segmentation path.
// The Java implementation uses a randomized low-rank SVD approximation.  Here
// the same kernel/cost formulation is evaluated against a fixed anchor sample
// (seed 1216) so Serial/OpenMP/CUDA/HIP/SYCL produce the same boundaries while
// keeping memory O(N * approximation-dimension), rather than O(N^2).
std::vector<Segment> kernel_segment_points(const Table& table, const Options& options,
                                           KernelStats* kernel_stats = nullptr) {
    std::vector<Segment> segments;
    std::size_t chromosome_begin = 0;
    while (chromosome_begin < table.points.size()) {
        std::size_t chromosome_end = chromosome_begin + 1;
        while (chromosome_end < table.points.size() &&
               table.points[chromosome_end].contig == table.points[chromosome_begin].contig &&
               table.points[chromosome_end].start <= table.points[chromosome_end - 1].end + 1) {
            ++chromosome_end;
        }
        const std::size_t n = chromosome_end - chromosome_begin;
        if (n < 2 || options.max_segments_per_contig < 2) {
            segments.push_back(summarize_range(table, chromosome_begin, chromosome_end));
            chromosome_begin = chromosome_end;
            continue;
        }

        bool has_allelic_signal = false;
        for (std::size_t i = chromosome_begin; i < chromosome_end; ++i)
            has_allelic_signal = has_allelic_signal || std::isfinite(table.points[i].segmentation_alt_fraction);
        // A linear copy-ratio kernel has an exact one-dimensional feature map
        // (x -> x).  The Java KernelSegmenter obtains the same map through an
        // SVD of the sampled outer-product kernel.  Materialising x_i*x_j for
        // every sampled anchor is not equivalent: it scales the segmentation
        // cost by the sampled norm and can change the changepoint penalty.
        // Keep the exact rank-one representation for the common copy-ratio-
        // only path; multidimensional/Gaussian paths retain the bounded
        // anchor embedding below.
        const bool exact_linear_copy = options.kernel_variance_copy_ratio == 0.0 && !has_allelic_signal;
        const std::size_t dimension = exact_linear_copy
            ? 1 : std::min(n, options.kernel_approximation_dimension);
        std::vector<std::size_t> anchors(dimension);
        std::uint64_t state = 1216;  // GATK KernelSegmenter.RANDOM_SEED
        for (std::size_t j = 0; j < dimension; ++j) {
            state = state * 1664525ULL + 1013904223ULL;
            anchors[j] = static_cast<std::size_t>(state % n);
        }

        Kokkos::View<double*> copy_values("model_segments_kernel_copy", n);
        Kokkos::View<double*> maf_values("model_segments_kernel_maf", n);
        Kokkos::View<double*> anchor_copy("model_segments_kernel_anchor_copy", dimension);
        Kokkos::View<double*> anchor_maf("model_segments_kernel_anchor_maf", dimension);
        auto host_copy = Kokkos::create_mirror_view(copy_values);
        auto host_maf = Kokkos::create_mirror_view(maf_values);
        auto host_anchor_copy = Kokkos::create_mirror_view(anchor_copy);
        auto host_anchor_maf = Kokkos::create_mirror_view(anchor_maf);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& point = table.points[chromosome_begin + i];
            host_copy(i) = point.value;
            host_maf(i) = std::isfinite(point.segmentation_alt_fraction)
                ? point.segmentation_alt_fraction : 0.5;
        }
        for (std::size_t j = 0; j < dimension; ++j) {
            host_anchor_copy(j) = host_copy(anchors[j]);
            host_anchor_maf(j) = host_maf(anchors[j]);
        }
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host_batch("model-segments-embedding-v1");
        host_batch.records = n * dimension;
        host_batch.bytes = (2 * n + 2 * dimension + n * dimension) * sizeof(double);
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-embedding");
        plan.begin_prepare(host_batch);
        Kokkos::deep_copy(copy_values, host_copy);
        Kokkos::deep_copy(maf_values, host_maf);
        Kokkos::deep_copy(anchor_copy, host_anchor_copy);
        Kokkos::deep_copy(anchor_maf, host_anchor_maf);

        Kokkos::View<double**> embedding("model_segments_kernel_embedding", n, dimension);
        const double copy_variance = options.kernel_variance_copy_ratio;
        const double maf_variance = options.kernel_variance_allele_fraction;
        const double maf_scale = options.kernel_scaling_allele_fraction;
        const double inverse_sqrt_two_pi = 1.0 / std::sqrt(2.0 * 3.14159265358979323846);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(n * dimension);
        device_batch.bind("copy_values", copy_values);
        device_batch.bind("maf_values", maf_values);
        device_batch.bind("anchor_copy", anchor_copy);
        device_batch.bind("anchor_maf", anchor_maf);
        device_batch.bind("embedding", embedding);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        Kokkos::parallel_for("model_segments_kernel_embedding",
            Kokkos::RangePolicy<ExecSpace>(0, n * dimension), KOKKOS_LAMBDA(const std::size_t linear) {
                const std::size_t i = linear / dimension;
                const std::size_t j = linear % dimension;
                const double copy_delta = copy_values(i) - anchor_copy(j);
                double copy_kernel;
                if (exact_linear_copy) {
                    // SVD(K=x*x^T) projects the linear kernel to x (up to a
                    // numerically negligible positive scale).  Use this
                    // exact feature directly so costs and penalties are
                    // independent of the random anchor sample.
                    copy_kernel = copy_values(i);
                } else if (copy_variance == 0.0) {
                    copy_kernel = copy_values(i) * anchor_copy(j);
                } else {
                    const double sigma = Kokkos::sqrt(copy_variance);
                    const double z = copy_delta / sigma;
                    copy_kernel = inverse_sqrt_two_pi / sigma * Kokkos::exp(-0.5 * z * z);
                }
                double allele_kernel = 0.0;
                if (has_allelic_signal) {
                    const double allele_delta = maf_values(i) - anchor_maf(j);
                    if (maf_variance == 0.0) {
                        allele_kernel = maf_values(i) * anchor_maf(j);
                    } else {
                        const double sigma = Kokkos::sqrt(maf_variance);
                        const double z = allele_delta / sigma;
                        allele_kernel = inverse_sqrt_two_pi / sigma * Kokkos::exp(-0.5 * z * z);
                    }
                }
                embedding(i, j) = copy_kernel + maf_scale * allele_kernel;
            });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += n * dimension;
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        auto host_embedding = Kokkos::create_mirror_view(embedding);
        Kokkos::deep_copy(host_embedding, embedding);

        std::vector<double> diagonal_prefix(n + 1, 0.0);
        std::vector<double> sums((n + 1) * dimension, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double diagonal = 0.0;
            for (std::size_t j = 0; j < dimension; ++j) {
                const double value = host_embedding(i, j);
                diagonal += value * value;
                sums[(i + 1) * dimension + j] = sums[i * dimension + j] + value;
            }
            diagonal_prefix[i + 1] = diagonal_prefix[i] + diagonal;
        }
        const auto segment_cost = [&](const std::size_t begin, const std::size_t end) {
            if (begin > end) return 0.0;
            const double length = static_cast<double>(end - begin + 1);
            double norm = 0.0;
            for (std::size_t j = 0; j < dimension; ++j) {
                const double total = sums[(end + 1) * dimension + j] - sums[begin * dimension + j];
                norm += total * total;
            }
            return (diagonal_prefix[end + 1] - diagonal_prefix[begin]) - norm / length;
        };

        std::vector<KernelCandidate> candidates;
        for (const auto window : options.window_sizes) {
            if (2 * window > n) continue;
            std::vector<KernelCandidate> local;
            for (std::size_t center = window - 1; center + window < n; ++center) {
                const auto left_begin = center + 1 - window;
                const auto right_end = center + window;
                const double score = segment_cost(left_begin, center) +
                    segment_cost(center + 1, right_end) - segment_cost(left_begin, right_end);
                local.push_back(KernelCandidate{center + 1, score});
            }
            std::sort(local.begin(), local.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.score == rhs.score ? lhs.boundary < rhs.boundary : lhs.score < rhs.score;
            });
            const auto limit = std::min(options.max_segments_per_contig - 1, local.size());
            candidates.insert(candidates.end(), local.begin(), local.begin() + limit);
        }
        std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.score == rhs.score ? lhs.boundary < rhs.boundary : lhs.score < rhs.score;
        });
        std::vector<std::size_t> boundaries;
        for (const auto& candidate : candidates) {
            if (boundaries.size() >= options.max_segments_per_contig - 1) break;
            const double changepoints = static_cast<double>(boundaries.size() + 1);
            const double penalty = options.changepoint_penalty_factor * changepoints *
                (1.0 + std::log(static_cast<double>(n) / (changepoints + 1.0e-10)));
            if (candidate.score + penalty >= -1.0e-12) continue;
            const auto distance_from_left = candidate.boundary;
            const auto distance_from_right = n - candidate.boundary;
            if (distance_from_left < options.min_points || distance_from_right < options.min_points) continue;
            bool separated = true;
            for (const auto boundary : boundaries)
                if (boundary > candidate.boundary ? boundary - candidate.boundary < options.min_points
                                                   : candidate.boundary - boundary < options.min_points)
                    separated = false;
            if (separated) boundaries.push_back(candidate.boundary);
        }
        std::sort(boundaries.begin(), boundaries.end());
        std::size_t begin = chromosome_begin;
        for (const auto boundary : boundaries) {
            segments.push_back(summarize_range(table, begin, chromosome_begin + boundary));
            begin = chromosome_begin + boundary;
        }
        segments.push_back(summarize_range(table, begin, chromosome_end));
        chromosome_begin = chromosome_end;
    }
    return segments;
}

// Multisample KernelSegmenter must consume the complete vector at each
// interval.  Averaging samples before calling the single-sample segmenter is
// not a valid reduction: for a linear kernel GATK computes
//   K(x_i,x_j) = sum_s x_i,s * x_j,s,
// whereas K(mean(x_i),mean(x_j)) can cancel anti-correlated events entirely.
// This path keeps the same bounded local-cost implementation as the native
// single-sample path but materializes a concatenated Kokkos feature vector for
// linear copy-ratio/allele-fraction kernels.  Gaussian kernels use the same
// fixed-anchor approximation and sum the per-sample kernels, matching the
// Java kernel definition at the feature-map boundary.
std::vector<Segment> kernel_segment_points_multisample(
    const Table& table, const Options& options, const std::size_t sample_count,
    KernelStats* kernel_stats = nullptr) {
    std::vector<Segment> segments;
    std::size_t chromosome_begin = 0;
    while (chromosome_begin < table.points.size()) {
        std::size_t chromosome_end = chromosome_begin + 1;
        while (chromosome_end < table.points.size() &&
               table.points[chromosome_end].contig == table.points[chromosome_begin].contig &&
               table.points[chromosome_end].start <= table.points[chromosome_end - 1].end + 1)
            ++chromosome_end;
        const std::size_t n = chromosome_end - chromosome_begin;
        // MultisampleMultidimensionalKernelSegmenter deliberately skips
        // segmentation below ten points per chromosome.
        if (n < 10 || options.max_segments_per_contig < 2) {
            segments.push_back(summarize_range(table, chromosome_begin, chromosome_end));
            chromosome_begin = chromosome_end;
            continue;
        }

        bool has_allelic_signal = false;
        for (std::size_t i = chromosome_begin; i < chromosome_end; ++i)
            for (const auto value : table.points[i].segmentation_alt_fraction_samples)
                has_allelic_signal = has_allelic_signal || std::isfinite(value);
        const bool linear_copy = options.kernel_variance_copy_ratio == 0.0;
        const bool linear_allele = options.kernel_variance_allele_fraction == 0.0;
        const bool exact_linear = linear_copy && (!has_allelic_signal || linear_allele);
        const std::size_t dimension = exact_linear
            ? sample_count * (has_allelic_signal ? 2 : 1)
            : std::min(n, options.kernel_approximation_dimension);
        if (dimension == 0) {
            segments.push_back(summarize_range(table, chromosome_begin, chromosome_end));
            chromosome_begin = chromosome_end;
            continue;
        }

        std::vector<std::size_t> anchors;
        if (!exact_linear) {
            anchors.resize(dimension);
            std::uint64_t state = 1216;
            for (std::size_t j = 0; j < dimension; ++j) {
                state = state * 1664525ULL + 1013904223ULL;
                anchors[j] = static_cast<std::size_t>(state % n);
            }
        }
        const double copy_variance = options.kernel_variance_copy_ratio;
        const double allele_variance = options.kernel_variance_allele_fraction;
        const double allele_scale = options.kernel_scaling_allele_fraction;
        const double inverse_sqrt_two_pi = 1.0 / std::sqrt(2.0 * 3.14159265358979323846);
        auto copy_value = [&](const std::size_t point, const std::size_t sample) {
            const auto& values = table.points[chromosome_begin + point].copy_ratio_samples;
            return sample < values.size() ? values[sample] : table.points[chromosome_begin + point].value;
        };
        auto allele_value = [&](const std::size_t point, const std::size_t sample) {
            const auto& values = table.points[chromosome_begin + point].segmentation_alt_fraction_samples;
            if (sample < values.size() && std::isfinite(values[sample])) return values[sample];
            return 0.5;
        };
        std::vector<double> host_copy_values(n * sample_count, 0.0);
        std::vector<double> host_allele_values(n * sample_count, 0.5);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                host_copy_values[i * sample_count + sample] = copy_value(i, sample);
                host_allele_values[i * sample_count + sample] = allele_value(i, sample);
            }
        }

        using ExecSpace = Kokkos::DefaultExecutionSpace;
        Kokkos::View<double**> features("model_segments_multisample_kernel_features", n, dimension);
        Kokkos::View<double**> copy_values("model_segments_multisample_copy_values", n, sample_count);
        Kokkos::View<double**> allele_values("model_segments_multisample_allele_values", n, sample_count);
        auto host_copy_view = Kokkos::create_mirror_view(copy_values);
        auto host_allele_view = Kokkos::create_mirror_view(allele_values);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                host_copy_view(i, sample) = host_copy_values[i * sample_count + sample];
                host_allele_view(i, sample) = host_allele_values[i * sample_count + sample];
            }
        Kokkos::View<std::size_t*> anchor_view("model_segments_multisample_kernel_anchors", anchors.size());
        auto host_anchor_view = Kokkos::create_mirror_view(anchor_view);
        for (std::size_t index = 0; index < anchors.size(); ++index) host_anchor_view(index) = anchors[index];
        fastgatk::core::HostBatch host_batch("model-segments-multisample-embedding-v1");
        host_batch.records = n * dimension;
        host_batch.bytes = n * dimension * sizeof(double) +
                           2 * n * sample_count * sizeof(double) + anchors.size() * sizeof(std::size_t);
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-multisample-embedding");
        plan.begin_prepare(host_batch);
        Kokkos::deep_copy(copy_values, host_copy_view);
        Kokkos::deep_copy(allele_values, host_allele_view);
        if (!anchors.empty()) Kokkos::deep_copy(anchor_view, host_anchor_view);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(n * dimension);
        device_batch.bind("features", features);
        device_batch.bind("copy_values", copy_values);
        device_batch.bind("allele_values", allele_values);
        if (!anchors.empty()) device_batch.bind("anchors", anchor_view);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        Kokkos::parallel_for("model_segments_multisample_feature_kernel",
                             Kokkos::RangePolicy<ExecSpace>(0, n * dimension),
                             KOKKOS_LAMBDA(const std::size_t linear) {
                                 const std::size_t i = linear / dimension;
                                 const std::size_t j = linear % dimension;
                                 double feature = 0.0;
                                 if (exact_linear) {
                                     feature = j < sample_count
                                         ? copy_values(i, j)
                                         : allele_scale * allele_values(i, j - sample_count);
                                 } else {
                                     const auto anchor = anchor_view(j);
                                     for (std::size_t sample = 0; sample < sample_count; ++sample) {
                                         const auto copy_delta = copy_values(i, sample) - copy_values(anchor, sample);
                                         if (copy_variance == 0.0) feature += copy_values(i, sample) * copy_values(anchor, sample);
                                         else {
                                             const auto sigma = Kokkos::sqrt(copy_variance);
                                             const auto z = copy_delta / sigma;
                                             feature += inverse_sqrt_two_pi / sigma * Kokkos::exp(-0.5 * z * z);
                                         }
                                         if (has_allelic_signal) {
                                             const auto allele_delta = allele_values(i, sample) - allele_values(anchor, sample);
                                             if (allele_variance == 0.0)
                                                 feature += allele_scale * allele_values(i, sample) * allele_values(anchor, sample);
                                             else {
                                                 const auto sigma = Kokkos::sqrt(allele_variance);
                                                 const auto z = allele_delta / sigma;
                                                 feature += allele_scale * inverse_sqrt_two_pi / sigma * Kokkos::exp(-0.5 * z * z);
                                             }
                                         }
                                     }
                                 }
                                 features(i, j) = feature;
                             });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += n * dimension;
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        auto host_embedding = Kokkos::create_mirror_view(features);
        Kokkos::deep_copy(host_embedding, features);
        std::vector<double> diagonal_prefix(n + 1, 0.0);
        std::vector<double> sums((n + 1) * dimension, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double diagonal = 0.0;
            for (std::size_t j = 0; j < dimension; ++j) {
                const auto value = host_embedding(i, j);
                diagonal += value * value;
                sums[(i + 1) * dimension + j] = sums[i * dimension + j] + value;
            }
            diagonal_prefix[i + 1] = diagonal_prefix[i] + diagonal;
        }
        const auto segment_cost = [&](const std::size_t begin, const std::size_t end) {
            if (begin > end) return 0.0;
            const auto length = static_cast<double>(end - begin + 1);
            double norm = 0.0;
            for (std::size_t j = 0; j < dimension; ++j) {
                const auto total = sums[(end + 1) * dimension + j] - sums[begin * dimension + j];
                norm += total * total;
            }
            return (diagonal_prefix[end + 1] - diagonal_prefix[begin]) - norm / length;
        };
        std::vector<KernelCandidate> candidates;
        for (const auto window : options.window_sizes) {
            if (2 * window > n) continue;
            std::vector<KernelCandidate> local;
            for (std::size_t center = window - 1; center + window < n; ++center) {
                const auto left_begin = center + 1 - window;
                const auto right_end = center + window;
                const auto score = segment_cost(left_begin, center) +
                    segment_cost(center + 1, right_end) - segment_cost(left_begin, right_end);
                local.push_back(KernelCandidate{center + 1, score});
            }
            std::sort(local.begin(), local.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.score == rhs.score ? lhs.boundary < rhs.boundary : lhs.score < rhs.score;
            });
            const auto limit = std::min(options.max_segments_per_contig - 1, local.size());
            candidates.insert(candidates.end(), local.begin(), local.begin() + limit);
        }
        std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.score == rhs.score ? lhs.boundary < rhs.boundary : lhs.score < rhs.score;
        });
        std::vector<std::size_t> boundaries;
        for (const auto& candidate : candidates) {
            if (boundaries.size() >= options.max_segments_per_contig - 1) break;
            const auto changepoints = static_cast<double>(boundaries.size() + 1);
            const auto penalty = options.changepoint_penalty_factor * changepoints *
                (1.0 + std::log(static_cast<double>(n) / (changepoints + 1.0e-10)));
            if (candidate.score + penalty >= -1.0e-12) continue;
            if (candidate.boundary < options.min_points || n - candidate.boundary < options.min_points)
                continue;
            bool separated = true;
            for (const auto boundary : boundaries)
                if (boundary > candidate.boundary ? boundary - candidate.boundary < options.min_points
                                                   : candidate.boundary - boundary < options.min_points)
                    separated = false;
            if (separated) boundaries.push_back(candidate.boundary);
        }
        std::sort(boundaries.begin(), boundaries.end());
        std::size_t begin = chromosome_begin;
        for (const auto boundary : boundaries) {
            segments.push_back(summarize_range(table, begin, chromosome_begin + boundary));
            begin = chromosome_begin + boundary;
        }
        segments.push_back(summarize_range(table, begin, chromosome_end));
        chromosome_begin = chromosome_end;
    }
    return segments;
}

KOKKOS_INLINE_FUNCTION std::uint64_t model_segments_mix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

KOKKOS_INLINE_FUNCTION double model_segments_uniform(std::uint64_t seed) {
    // Use the top 53 bits, matching the precision of a host double while
    // remaining independent of the selected Kokkos execution space.
    const auto bits = model_segments_mix64(seed) >> 11;
    return (static_cast<double>(bits) + 0.5) / 9007199254740992.0;
}

KOKKOS_INLINE_FUNCTION double model_segments_normal(std::uint64_t seed) {
    constexpr double two_pi = 6.283185307179586476925286766559;
    const auto first = model_segments_uniform(seed);
    const auto second = model_segments_uniform(seed ^ 0xd1b54a32d192ed03ULL);
    return Kokkos::sqrt(-2.0 * Kokkos::log(first)) * Kokkos::cos(two_pi * second);
}

KOKKOS_INLINE_FUNCTION double model_segments_log_density(const double value,
                                                          const double mean,
                                                          const double sigma) {
    const auto z = (value - mean) / sigma;
    return -0.5 * z * z;
}

// These helpers are a direct transcription of the release-pinned GATK
// AlleleFractionLikelihoods Laplace approximation.  The previous native
// conditional pass treated allelic fractions as a Gaussian observation and
// estimated bias from (major + 1) / (minor + 1).  That is a useful fallback
// for MAF-only tables, but it is not the count model used by
// AlleleFractionModeller: GATK integrates the site-level Gamma bias and then
// mixes alt-minor, ref-minor and uniform-outlier components.  Keep this
// implementation Kokkos-callable so the same formula can be used by host
// coordinate control and a future device sampler.
KOKKOS_INLINE_FUNCTION double model_segments_af_log(const double value) {
    return Kokkos::log(Kokkos::fmax(1.0e-10, value));
}

KOKKOS_INLINE_FUNCTION double model_segments_af_prior_log_density(
    const double fraction, const double prior_alpha) {
    // AlleleFractionSamplers uses Beta(alpha, 1) evaluated at 2*f.  The
    // normalizing constant is independent of f and can be omitted from the
    // Metropolis/slice comparison.  Clamp only the log argument so the
    // physical boundary remains representable without a NaN on device.
    return (prior_alpha - 1.0) * model_segments_af_log(
        2.0 * Kokkos::fmax(1.0e-10, Kokkos::fmin(0.5, fraction)));
}

KOKKOS_INLINE_FUNCTION double model_segments_af_log_sum3(const double first,
                                                           const double second,
                                                           const double third) {
    const double maximum = Kokkos::fmax(first, Kokkos::fmax(second, third));
    return maximum + Kokkos::log(Kokkos::exp(first - maximum) +
                                 Kokkos::exp(second - maximum) +
                                 Kokkos::exp(third - maximum));
}

KOKKOS_INLINE_FUNCTION double model_segments_af_bias_mode(const double alpha,
                                                           const double beta,
                                                           const double fraction,
                                                           const double alt,
                                                           const double ref) {
    const double one_minus_fraction = 1.0 - fraction;
    const double w = one_minus_fraction * (alt - alpha + 1.0) + beta * fraction;
    const double discriminant = w * w + 4.0 * beta * fraction * one_minus_fraction *
        (ref + alpha - 1.0);
    return Kokkos::fmax((Kokkos::sqrt(Kokkos::fmax(0.0, discriminant)) - w) /
                        (2.0 * beta * one_minus_fraction), 1.0e-10);
}

KOKKOS_INLINE_FUNCTION double model_segments_af_het_log_likelihood(
    const double mean_bias, const double bias_variance,
    const double outlier_probability, const double fraction,
    const std::int64_t alt_count, const std::int64_t ref_count) {
    const double safe_variance = Kokkos::fmax(1.0e-10, bias_variance);
    const double beta = Kokkos::fmax(1.0e-10, mean_bias / safe_variance);
    const double alpha = Kokkos::fmax(1.0e-10, mean_bias * beta);
    const double minor = Kokkos::fmax(1.0e-10, Kokkos::fmin(0.5, fraction));
    const double major = 1.0 - minor;
    const double alt = Kokkos::fmax(0.0, static_cast<double>(alt_count));
    const double ref = Kokkos::fmax(0.0, static_cast<double>(ref_count));
    const double total = alt + ref;
    const double log_common = alpha * model_segments_af_log(beta) - Kokkos::lgamma(alpha);

    const auto component = [&](const double site_fraction, const bool alt_minor) {
        const double lambda = model_segments_af_bias_mode(alpha, beta, site_fraction, alt, ref);
        const double y = (1.0 - site_fraction) /
            (site_fraction + (1.0 - site_fraction) * lambda);
        const double curvature = total * y * y -
            (ref + alpha - 1.0) / (lambda * lambda);
        const double effective_alpha = Kokkos::fmax(1.0e-10,
            1.0 - curvature * lambda * lambda);
        const double effective_beta = Kokkos::fmax(1.0e-10,
            -curvature * lambda);
        // ``site_fraction`` is f for the ALT-minor component and 1-f for the
        // REF-minor component (the caller passes major as the second
        // component's site fraction).  In both cases the Java expression is
        // therefore a*log(site_fraction) + r*log(1-site_fraction).
        const double alt_fraction = site_fraction;
        const double ref_fraction = 1.0 - site_fraction;
        const double log_constant = log_common + alt * model_segments_af_log(alt_fraction) +
            ref * model_segments_af_log(ref_fraction) +
            (ref + alpha - effective_alpha) * model_segments_af_log(lambda) +
            (effective_beta - beta) * lambda - total *
            model_segments_af_log(alt_fraction + ref_fraction * lambda);
        return log_constant + Kokkos::lgamma(effective_alpha) -
            effective_alpha * model_segments_af_log(effective_beta);
    };
    const double log_not_outlier = model_segments_af_log(
        Kokkos::fmax(1.0e-10, (1.0 - outlier_probability) / 2.0));
    const double alt_minor_log = log_not_outlier + component(minor, true);
    const double ref_minor_log = log_not_outlier + component(major, false);
    // log(a! r! / (n + 1)!) is the uniform component used by GATK.
    const double outlier_log = model_segments_af_log(Kokkos::fmax(1.0e-10,
        outlier_probability)) - model_segments_af_log(total + 1.0) -
        Kokkos::lgamma(total + 1.0) + Kokkos::lgamma(alt + 1.0) + Kokkos::lgamma(ref + 1.0);
    return model_segments_af_log_sum3(alt_minor_log, ref_minor_log, outlier_log);
}

template <typename Function>
double model_segments_af_argmax(Function&& function, const double minimum,
                                 const double maximum, const double guess) {
    // Commons-Math Brent in GATK uses 1e-3 relative and absolute tolerances.
    // A fixed grid followed by golden-section refinement is deterministic,
    // bounded, and remains stable when the log-likelihood has a flat plateau
    // (a common occurrence for low-depth allele counts).
    (void)guess;
    constexpr std::size_t grid_points = 64;
    std::size_t best_index = 0;
    double best_value = function(minimum);
    for (std::size_t index = 1; index <= grid_points; ++index) {
        const double value = minimum + (maximum - minimum) *
            static_cast<double>(index) / static_cast<double>(grid_points);
        const double score = function(value);
        if (score > best_value) {
            best_value = score;
            best_index = index;
        }
    }
    if (best_index == 0 || best_index == grid_points)
        return best_index == 0 ? minimum : maximum;
    double left = minimum + (maximum - minimum) *
        static_cast<double>(best_index - 1) / static_cast<double>(grid_points);
    double right = minimum + (maximum - minimum) *
        static_cast<double>(best_index + 1) / static_cast<double>(grid_points);
    constexpr double golden = 0.6180339887498948482;
    double x1 = right - golden * (right - left);
    double x2 = left + golden * (right - left);
    double y1 = function(x1);
    double y2 = function(x2);
    for (std::size_t iteration = 0; iteration < 64; ++iteration) {
        if (y1 < y2) {
            left = x1;
            x1 = x2;
            y1 = y2;
            x2 = left + golden * (right - left);
            y2 = function(x2);
        } else {
            right = x2;
            x2 = x1;
            y2 = y1;
            x1 = right - golden * (right - left);
            y1 = function(x1);
        }
        if (right - left <= std::max(1.0e-3,
                                    1.0e-3 * std::abs(left + right) * 0.5)) break;
    }
    const double answer = (left + right) * 0.5;
    return std::max(minimum, std::min(maximum, answer));
}

double model_segments_quantile(std::vector<double> values, const double probability) {
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](const double value) { return !std::isfinite(value); }),
                 values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const auto position = probability * static_cast<double>(values.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(position));
    const auto upper = std::min(values.size() - 1, lower + 1);
    const auto fraction = position - static_cast<double>(lower);
    return values[lower] + fraction * (values[upper] - values[lower]);
}

// Fit the copy-ratio part of CopyRatioModeller's Gibbs state using a bounded
// deterministic conditional approximation.  GATK's model assigns each point
// by interval midpoint and alternates: point-level outlier indicators,
// segment means, global variance, and global outlier probability.  The native
// chain remains deliberately reproducible across Kokkos backends, so this
// prepass uses the conditional responsibilities (rather than backend-specific
// random draws) and then initializes the existing bounded sampler with the
// resulting state.  Raw points are carried in Kokkos Views; only the global
// scalar reduction is collected in a fixed host order, preserving Serial /
// OpenMP byte identity.
CopyRatioConditionalModel fit_copy_ratio_conditional_model(
    const Table& table, const std::vector<Segment>& segments,
    const Options& options, KernelStats* kernel_stats = nullptr) {
    CopyRatioConditionalModel model;
    model.segment_means.resize(segments.size(), 0.0);
    if (segments.empty() || !options.probabilistic_requested ||
        (options.mode != "BOTH" && options.mode != "COPY_RATIO"))
        return model;

    // Construct the same midpoint partition consumed by GATK's
    // CopyRatioSegmentedData.  Finite observations are the only values that
    // contribute to the numerical model; malformed/non-finite rows retain the
    // existing output contract but do not poison the conditional state.
    std::vector<double> host_values;
    std::vector<std::uint64_t> host_offsets(segments.size() + 1, 0);
    host_values.reserve(table.points.size());
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        host_offsets[segment] = static_cast<std::uint64_t>(host_values.size());
        for (const auto& point : table.points) {
            if (point.contig != segments[segment].contig || !std::isfinite(point.value)) continue;
            const auto midpoint = point.start + (point.end - point.start) / 2;
            if (midpoint < segments[segment].start || midpoint > segments[segment].end) continue;
            host_values.push_back(point.value);
            minimum = std::min(minimum, point.value);
            maximum = std::max(maximum, point.value);
        }
        const auto& segment_summary = segments[segment];
        model.segment_means[segment] = std::isfinite(segment_summary.mean)
            ? std::max(-50.0, std::min(10.0, segment_summary.mean)) : 0.0;
    }
    host_offsets[segments.size()] = static_cast<std::uint64_t>(host_values.size());
    if (host_values.empty()) return model;

    // GATK's uniform outlier component is normalized over the observed data
    // range.  Use the documented model bounds for a constant range in a
    // degenerate fixture, matching the Java fallback for a zero data range.
    const double observed_range = maximum > minimum ? maximum - minimum : 60.0;
    const double outlier_uniform_log_likelihood = -std::log(std::max(observed_range, 1.0e-12));
    constexpr double variance_min = 1.0e-6;
    double variance = 0.0;
    std::size_t variance_count = 0;
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        const auto begin = static_cast<std::size_t>(host_offsets[segment]);
        const auto end = static_cast<std::size_t>(host_offsets[segment + 1]);
        for (std::size_t index = begin; index < end; ++index) {
            const double residual = host_values[index] - model.segment_means[segment];
            variance += residual * residual;
            ++variance_count;
        }
    }
    variance = std::max(variance_min,
                        variance_count == 0 ? variance_min : variance / static_cast<double>(variance_count));
    double outlier_probability = 0.05;  // GATK OUTLIER_PROBABILITY_INITIAL.

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<double*> values("model_segments_copy_conditional_values", host_values.size());
    Kokkos::View<std::uint64_t*> offsets("model_segments_copy_conditional_offsets", host_offsets.size());
    Kokkos::View<double*> means("model_segments_copy_conditional_means", segments.size());
    Kokkos::View<double*> updated_means("model_segments_copy_conditional_updated_means", segments.size());
    Kokkos::View<double*> inlier_weights("model_segments_copy_conditional_inlier_weights", segments.size());
    Kokkos::View<double*> residual_sums("model_segments_copy_conditional_residual_sums", segments.size());
    Kokkos::View<double*> outlier_sums("model_segments_copy_conditional_outlier_sums", segments.size());
    auto host_mean_view = Kokkos::create_mirror_view(means);
    auto host_updated_means = Kokkos::create_mirror_view(updated_means);
    auto host_inlier_weights = Kokkos::create_mirror_view(inlier_weights);
    auto host_residual_sums = Kokkos::create_mirror_view(residual_sums);
    auto host_outlier_sums = Kokkos::create_mirror_view(outlier_sums);
    auto host_value_view = Kokkos::create_mirror_view(values);
    auto host_offset_view = Kokkos::create_mirror_view(offsets);
    for (std::size_t index = 0; index < host_values.size(); ++index) host_value_view(index) = host_values[index];
    for (std::size_t index = 0; index < host_offsets.size(); ++index) host_offset_view(index) = host_offsets[index];
    for (std::size_t index = 0; index < model.segment_means.size(); ++index) host_mean_view(index) = model.segment_means[index];
    Kokkos::deep_copy(values, host_value_view);
    Kokkos::deep_copy(offsets, host_offset_view);
    Kokkos::deep_copy(means, host_mean_view);

    // A short fixed-point pass is intentional: the exact Java implementation
    // uses slice samplers and release-specific RNGs, whereas this approximation
    // must have stable output across Serial/OpenMP/CUDA/HIP/SYCL.  Scaling the
    // bound with the requested chain still lets larger production runs tighten
    // the conditional state without creating an unbounded startup cost.
    const auto requested_iterations = options.num_burn_in_iterations + options.num_samples;
    const auto conditional_iterations = std::min<std::size_t>(32,
        std::max<std::size_t>(4, requested_iterations / 16 + 4));
    constexpr double two_pi = 6.283185307179586476925286766559;
    for (std::size_t iteration = 0; iteration < conditional_iterations; ++iteration) {
        fastgatk::core::HostBatch host_batch("model-segments-copy-ratio-conditional-v1");
        host_batch.records = segments.size();
        host_batch.bytes = (host_values.size() * sizeof(double)) +
                           (host_offsets.size() * sizeof(std::uint64_t)) +
                           (5 * segments.size() * sizeof(double));
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-copy-ratio-conditional");
        plan.begin_prepare(host_batch);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(segments.size());
        device_batch.bind("values", values);
        device_batch.bind("offsets", offsets);
        device_batch.bind("means", means);
        device_batch.bind("updated_means", updated_means);
        device_batch.bind("inlier_weights", inlier_weights);
        device_batch.bind("residual_sums", residual_sums);
        device_batch.bind("outlier_sums", outlier_sums);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        const double current_variance = variance;
        const double current_outlier_probability = outlier_probability;
        Kokkos::parallel_for("model_segments_copy_ratio_conditionals",
            Kokkos::RangePolicy<ExecSpace>(0, segments.size()),
            KOKKOS_LAMBDA(const std::size_t segment) {
                const auto begin = offsets(segment);
                const auto end = offsets(segment + 1);
                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                double residual_sum = 0.0;
                double outlier_sum = 0.0;
                for (auto index = begin; index < end; ++index) {
                    const double residual = values(index) - means(segment);
                    const double normal_log = Kokkos::log(Kokkos::fmax(1.0e-12,
                            1.0 - current_outlier_probability)) -
                        0.5 * Kokkos::log(two_pi * current_variance) -
                        0.5 * residual * residual / current_variance;
                    const double outlier_log = Kokkos::log(Kokkos::fmax(1.0e-12,
                            current_outlier_probability)) + outlier_uniform_log_likelihood;
                    const double max_log = Kokkos::fmax(normal_log, outlier_log);
                    const double normal_component = Kokkos::exp(normal_log - max_log);
                    const double outlier_component = Kokkos::exp(outlier_log - max_log);
                    const double outlier_probability_given_point = outlier_component /
                        Kokkos::fmax(1.0e-12, normal_component + outlier_component);
                    const double weight = 1.0 - outlier_probability_given_point;
                    weighted_sum += weight * values(index);
                    weight_sum += weight;
                    residual_sum += weight * residual * residual;
                    outlier_sum += outlier_probability_given_point;
                }
                updated_means(segment) = end > begin
                    ? Kokkos::fmax(-50.0, Kokkos::fmin(10.0,
                        weighted_sum / Kokkos::fmax(1.0e-12, weight_sum)))
                    : means(segment);
                inlier_weights(segment) = weight_sum;
                residual_sums(segment) = residual_sum;
                outlier_sums(segment) = outlier_sum;
            });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += segments.size();
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        Kokkos::deep_copy(host_updated_means, updated_means);
        Kokkos::deep_copy(host_inlier_weights, inlier_weights);
        Kokkos::deep_copy(host_residual_sums, residual_sums);
        Kokkos::deep_copy(host_outlier_sums, outlier_sums);
        double total_weight = 0.0;
        double total_residual = 0.0;
        double total_outliers = 0.0;
        for (std::size_t segment = 0; segment < segments.size(); ++segment) {
            model.segment_means[segment] = host_updated_means(segment);
            total_weight += host_inlier_weights(segment);
            total_residual += host_residual_sums(segment);
            total_outliers += host_outlier_sums(segment);
        }
        variance = std::max(variance_min,
            total_weight > 0.0 ? total_residual / total_weight : variance_min);
        // Beta(5,95), the same prior used by GATK's OutlierProbabilitySampler.
        outlier_probability = std::min(0.999, std::max(0.001,
            (5.0 + total_outliers) /
            (100.0 + static_cast<double>(host_values.size()))));
        Kokkos::deep_copy(host_mean_view, host_updated_means);
        Kokkos::deep_copy(means, host_mean_view);
    }
    model.variance = variance;
    model.outlier_probability = outlier_probability;
    model.iterations = conditional_iterations;
    return model;
}

AlleleFractionConditionalModel fit_allele_fraction_conditional_model(
    const Table& table, const std::vector<Segment>& segments,
    const Options& options, KernelStats* kernel_stats = nullptr) {
    AlleleFractionConditionalModel model;
    model.segment_means.resize(segments.size(), std::numeric_limits<double>::quiet_NaN());
    if (segments.empty() || !options.probabilistic_requested ||
        (options.mode != "BOTH" && options.mode != "ALLELE_FRACTION"))
        return model;

    // The production Java model consumes locus-level ref/alt counts.  The
    // native table has one copy-ratio interval per row, so retain the exact
    // counts aggregated into that interval when available and otherwise use
    // the already materialized MAF as a unit-weight observation.  This keeps
    // the conditional path useful for legacy MAF-only tables while making the
    // count-backed path faithful to the GATK evidence model.
    std::vector<double> host_values;
    std::vector<double> host_weights;
    std::vector<std::int64_t> host_refs;
    std::vector<std::int64_t> host_alts;
    std::vector<std::uint8_t> host_count_backed;
    std::vector<std::uint64_t> host_offsets(segments.size() + 1, 0);
    host_values.reserve(table.points.size());
    host_weights.reserve(table.points.size());
    host_refs.reserve(table.points.size());
    host_alts.reserve(table.points.size());
    host_count_backed.reserve(table.points.size());
    double ratio_sum = 0.0;
    double ratio_sq_sum = 0.0;
    double ratio_weight = 0.0;
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        host_offsets[segment] = static_cast<std::uint64_t>(host_values.size());
        for (const auto& point : table.points) {
            if (point.contig != segments[segment].contig ||
                !std::isfinite(point.minor_allele_fraction)) continue;
            const auto midpoint = point.start + (point.end - point.start) / 2;
            if (midpoint < segments[segment].start || midpoint > segments[segment].end) continue;
            double value = std::min(0.5, std::max(0.0, point.minor_allele_fraction));
            double weight = 1.0;
            if (point.allele_counts_available) {
                const auto total = point.allele_ref_count + point.allele_alt_count;
                if (total > 0) {
                    const auto minor = std::min(point.allele_ref_count, point.allele_alt_count);
                    value = std::min(0.5, std::max(0.0,
                        static_cast<double>(minor) / static_cast<double>(total)));
                    weight = static_cast<double>(total);
                    const auto major = std::max(point.allele_ref_count, point.allele_alt_count);
                    const double ratio = static_cast<double>(major + 1) /
                        static_cast<double>(minor + 1);
                    ratio_sum += ratio * weight;
                    ratio_sq_sum += ratio * ratio * weight;
                    ratio_weight += weight;
                }
            }
            host_values.push_back(value);
            host_weights.push_back(std::max(1.0, weight));
            host_refs.push_back(point.allele_ref_count);
            host_alts.push_back(point.allele_alt_count);
            host_count_backed.push_back(point.allele_counts_available ? 1 : 0);
        }
        model.segment_means[segment] = std::isfinite(segments[segment].mean_minor_allele_fraction)
            ? std::min(0.5, std::max(0.0, segments[segment].mean_minor_allele_fraction)) : 0.25;
    }
    host_offsets[segments.size()] = static_cast<std::uint64_t>(host_values.size());
    if (host_values.empty()) return model;

    double mean_bias = ratio_weight > 0.0 ? ratio_sum / ratio_weight : 1.0;
    mean_bias = std::min(5.0, std::max(0.0, mean_bias));
    double bias_variance = ratio_weight > 0.0
        ? std::max(1.0e-10, ratio_sq_sum / ratio_weight - mean_bias * mean_bias) : 0.05;
    bias_variance = std::min(0.5, bias_variance);
    double outlier_probability = 0.01;  // AlleleFractionInitializer default.

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<double*> values("model_segments_allele_conditional_values", host_values.size());
    Kokkos::View<double*> weights("model_segments_allele_conditional_weights", host_weights.size());
    Kokkos::View<std::int64_t*> refs("model_segments_allele_conditional_refs", host_refs.size());
    Kokkos::View<std::int64_t*> alts("model_segments_allele_conditional_alts", host_alts.size());
    Kokkos::View<std::uint8_t*> count_backed("model_segments_allele_conditional_count_backed", host_count_backed.size());
    Kokkos::View<std::uint64_t*> offsets("model_segments_allele_conditional_offsets", host_offsets.size());
    Kokkos::View<double*> means("model_segments_allele_conditional_means", segments.size());
    Kokkos::View<double*> initial_means("model_segments_allele_conditional_initial_means", segments.size());
    Kokkos::View<double*> updated_means("model_segments_allele_conditional_updated_means", segments.size());
    Kokkos::View<double*> inlier_weights("model_segments_allele_conditional_inlier_weights", segments.size());
    Kokkos::View<double*> residual_sums("model_segments_allele_conditional_residual_sums", segments.size());
    Kokkos::View<double*> outlier_sums("model_segments_allele_conditional_outlier_sums", segments.size());
    auto host_values_view = Kokkos::create_mirror_view(values);
    auto host_weights_view = Kokkos::create_mirror_view(weights);
    auto host_refs_view = Kokkos::create_mirror_view(refs);
    auto host_alts_view = Kokkos::create_mirror_view(alts);
    auto host_count_backed_view = Kokkos::create_mirror_view(count_backed);
    auto host_offsets_view = Kokkos::create_mirror_view(offsets);
    auto host_means_view = Kokkos::create_mirror_view(means);
    auto host_initial_means = Kokkos::create_mirror_view(initial_means);
    auto host_updated_means = Kokkos::create_mirror_view(updated_means);
    auto host_inlier_weights = Kokkos::create_mirror_view(inlier_weights);
    auto host_residual_sums = Kokkos::create_mirror_view(residual_sums);
    auto host_outlier_sums = Kokkos::create_mirror_view(outlier_sums);
    for (std::size_t index = 0; index < host_values.size(); ++index) {
        host_values_view(index) = host_values[index];
        host_weights_view(index) = host_weights[index];
        host_refs_view(index) = host_refs[index];
        host_alts_view(index) = host_alts[index];
        host_count_backed_view(index) = host_count_backed[index];
    }
    for (std::size_t index = 0; index < host_offsets.size(); ++index) host_offsets_view(index) = host_offsets[index];
    for (std::size_t index = 0; index < model.segment_means.size(); ++index) host_means_view(index) = model.segment_means[index];
    Kokkos::deep_copy(values, host_values_view);
    Kokkos::deep_copy(weights, host_weights_view);
    Kokkos::deep_copy(refs, host_refs_view);
    Kokkos::deep_copy(alts, host_alts_view);
    Kokkos::deep_copy(count_backed, host_count_backed_view);
    Kokkos::deep_copy(offsets, host_offsets_view);
    Kokkos::deep_copy(means, host_means_view);

    // Match AlleleFractionInitializer.calculateInitialMinorFractions.  The
    // incomplete-beta value is the posterior responsibility that ALT is the
    // minor allele at f=1/2.  Aggregate responsibility-weighted minor reads
    // per segment and add GATK's single flat-prior pseudocount.  This is a
    // real count-backed model initialization, rather than the old weighted
    // arithmetic MAF, and runs in the selected Kokkos execution space.
    {
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host_batch("model-segments-allele-initialization-v1");
        host_batch.records = segments.size();
        host_batch.bytes = (2 * host_values.size() * sizeof(std::int64_t)) +
                           host_count_backed.size() * sizeof(std::uint8_t) +
                           host_offsets.size() * sizeof(std::uint64_t) +
                           2 * segments.size() * sizeof(double);
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-allele-initialization");
        plan.begin_prepare(host_batch);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(segments.size());
        device_batch.bind("refs", refs);
        device_batch.bind("alts", alts);
        device_batch.bind("count_backed", count_backed);
        device_batch.bind("offsets", offsets);
        device_batch.bind("means", means);
        device_batch.bind("initial_means", initial_means);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        Kokkos::parallel_for("model_segments_allele_fraction_initialization",
            Kokkos::RangePolicy<ExecSpace>(0, segments.size()),
            KOKKOS_LAMBDA(const std::size_t segment) {
                const auto begin = offsets(segment);
                const auto end = offsets(segment + 1);
                double weighted_minor_reads = 0.0;
                double total_reads = 0.0;
                bool has_count_backed = false;
                for (auto index = begin; index < end; ++index) {
                    if (count_backed(index) == 0) continue;
                    const double alt = static_cast<double>(alts(index));
                    const double ref = static_cast<double>(refs(index));
                    const double total = alt + ref;
                    if (!(total > 0.0)) continue;
                    const double alt_minor_responsibility =
                        model_segments_regularized_beta(0.5, alt + 1.0, ref + 1.0);
                    weighted_minor_reads += alt_minor_responsibility * alt +
                        (1.0 - alt_minor_responsibility) * ref;
                    total_reads += total;
                    has_count_backed = true;
                }
                initial_means(segment) = has_count_backed && total_reads > 0.0
                    ? (weighted_minor_reads + 1.0) / (total_reads + 2.0)
                    : means(segment);
            });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += segments.size();
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        Kokkos::deep_copy(host_initial_means, initial_means);
        model.initial_segment_means.resize(segments.size(), std::numeric_limits<double>::quiet_NaN());
        for (std::size_t segment = 0; segment < segments.size(); ++segment) {
            model.initial_segment_means[segment] = host_initial_means(segment);
            model.segment_means[segment] = host_initial_means(segment);
        }
        Kokkos::deep_copy(host_means_view, initial_means);
        Kokkos::deep_copy(means, host_means_view);
    }

    const auto requested_iterations = options.num_burn_in_iterations + options.num_samples;
    const auto conditional_iterations = std::min<std::size_t>(32,
        std::max<std::size_t>(4, requested_iterations / 16 + 4));
    constexpr double two_pi = 6.283185307179586476925286766559;
    for (std::size_t iteration = 0; iteration < conditional_iterations; ++iteration) {
        fastgatk::core::HostBatch host_batch("model-segments-allele-conditional-v1");
        host_batch.records = segments.size();
        host_batch.bytes = (host_values.size() * 2 * sizeof(double)) +
                           (host_offsets.size() * sizeof(std::uint64_t)) +
                           (5 * segments.size() * sizeof(double)) +
                           (2 * host_values.size() * sizeof(std::int64_t)) +
                           host_count_backed.size() * sizeof(std::uint8_t);
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-allele-conditional");
        plan.begin_prepare(host_batch);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(segments.size());
        device_batch.bind("values", values);
        device_batch.bind("weights", weights);
        device_batch.bind("refs", refs);
        device_batch.bind("alts", alts);
        device_batch.bind("count_backed", count_backed);
        device_batch.bind("offsets", offsets);
        device_batch.bind("means", means);
        device_batch.bind("updated_means", updated_means);
        device_batch.bind("inlier_weights", inlier_weights);
        device_batch.bind("residual_sums", residual_sums);
        device_batch.bind("outlier_sums", outlier_sums);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        const double current_bias_variance = bias_variance;
        const double current_outlier_probability = outlier_probability;
        const double prior_alpha = std::max(1.0, options.minor_allele_fraction_prior_alpha);
        const double prior_mean = prior_alpha / (2.0 * (prior_alpha + 1.0));
        Kokkos::parallel_for("model_segments_allele_fraction_conditionals",
            Kokkos::RangePolicy<ExecSpace>(0, segments.size()),
            KOKKOS_LAMBDA(const std::size_t segment) {
                const auto begin = offsets(segment);
                const auto end = offsets(segment + 1);
                const double fraction_variance = Kokkos::fmax(1.0e-5,
                    0.25 * current_bias_variance);
                double weighted_sum = prior_alpha * prior_mean;
                double weight_sum = prior_alpha;
                double residual_sum = 0.0;
                double outlier_sum = 0.0;
                for (auto index = begin; index < end; ++index) {
                    const double residual = values(index) - means(segment);
                    const double normal_log = Kokkos::log(Kokkos::fmax(1.0e-12,
                            1.0 - current_outlier_probability)) -
                        0.5 * Kokkos::log(two_pi * fraction_variance) -
                        0.5 * residual * residual / fraction_variance;
                    const double outlier_log = Kokkos::log(Kokkos::fmax(1.0e-12,
                            current_outlier_probability)) - Kokkos::log(0.5);
                    const double max_log = Kokkos::fmax(normal_log, outlier_log);
                    const double normal_component = Kokkos::exp(normal_log - max_log);
                    const double outlier_component = Kokkos::exp(outlier_log - max_log);
                    const double point_outlier_probability = outlier_component /
                        Kokkos::fmax(1.0e-12, normal_component + outlier_component);
                    const double inlier_weight = (1.0 - point_outlier_probability) * weights(index);
                    weighted_sum += inlier_weight * values(index);
                    weight_sum += inlier_weight;
                    residual_sum += inlier_weight * residual * residual;
                    outlier_sum += point_outlier_probability;
                }
                updated_means(segment) = Kokkos::fmax(0.0, Kokkos::fmin(0.5,
                    weighted_sum / Kokkos::fmax(1.0e-12, weight_sum)));
                inlier_weights(segment) = weight_sum;
                residual_sums(segment) = residual_sum;
                outlier_sums(segment) = outlier_sum;
            });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += segments.size();
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        Kokkos::deep_copy(host_updated_means, updated_means);
        Kokkos::deep_copy(host_inlier_weights, inlier_weights);
        Kokkos::deep_copy(host_residual_sums, residual_sums);
        Kokkos::deep_copy(host_outlier_sums, outlier_sums);
        double total_weight = 0.0;
        double total_residual = 0.0;
        double total_outliers = 0.0;
        for (std::size_t segment = 0; segment < segments.size(); ++segment) {
            model.segment_means[segment] = host_updated_means(segment);
            total_weight += host_inlier_weights(segment);
            total_residual += host_residual_sums(segment);
            total_outliers += host_outlier_sums(segment);
        }
        const double fraction_variance = total_weight > 0.0
            ? total_residual / total_weight : 1.0e-5;
        bias_variance = std::min(0.5, std::max(1.0e-10, 4.0 * fraction_variance));
        outlier_probability = std::min(0.15, std::max(0.0,
            (1.0 + total_outliers) /
            (100.0 + static_cast<double>(host_values.size()))));
        Kokkos::deep_copy(host_means_view, host_updated_means);
        Kokkos::deep_copy(means, host_means_view);
    }

    // Refine count-backed data with the same partially-collapsed likelihood
    // used by GATK's AlleleFractionInitializer.  The responsibility pass
    // above remains useful for MAF-only legacy input and supplies telemetry;
    // for real AllelicCount evidence, however, publishing its Gaussian
    // surrogate as the final state would systematically overestimate bias
    // (especially for asymmetric low-depth loci).  Coordinate ascent mirrors
    // the Java initializer's order and convergence criterion, while the
    // bounded one-dimensional optimizer keeps runtime predictable for large
    // cohorts.
    bool has_count_backed = false;
    std::vector<std::vector<std::size_t>> exact_indices(segments.size());
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        const auto begin = static_cast<std::size_t>(host_offsets[segment]);
        const auto end = static_cast<std::size_t>(host_offsets[segment + 1]);
        for (std::size_t index = begin; index < end; ++index) {
            if (host_count_backed[index] != 0) {
                exact_indices[segment].push_back(index);
                has_count_backed = true;
            }
        }
    }
    if (has_count_backed) {
        std::vector<double> exact_fractions(segments.size(), 0.25);
        for (std::size_t segment = 0; segment < segments.size(); ++segment) {
            if (segment < model.initial_segment_means.size() &&
                std::isfinite(model.initial_segment_means[segment]))
                exact_fractions[segment] = std::min(0.5, std::max(0.0,
                    model.initial_segment_means[segment]));
        }
        double exact_mean_bias = 1.0;
        double exact_bias_variance = 0.05;
        double exact_outlier_probability = 0.01;
        const auto exact_log_likelihood = [&](const double mean_bias,
                                              const double bias_variance,
                                              const double outlier_probability,
                                              const std::vector<double>& fractions) {
            double total = 0.0;
            for (std::size_t segment = 0; segment < exact_indices.size(); ++segment) {
                for (const auto index : exact_indices[segment]) {
                    total += model_segments_af_het_log_likelihood(
                        mean_bias, bias_variance, outlier_probability,
                        fractions[segment], host_alts[index], host_refs[index]);
                }
            }
            return total;
        };
        std::size_t exact_iterations = 0;
        double previous_log_likelihood = -std::numeric_limits<double>::infinity();
        for (; exact_iterations < 50; ++exact_iterations) {
            const double current_log_likelihood = exact_log_likelihood(
                exact_mean_bias, exact_bias_variance, exact_outlier_probability,
                exact_fractions);
            exact_mean_bias = model_segments_af_argmax(
                [&](const double value) {
                    return exact_log_likelihood(value, exact_bias_variance,
                                                exact_outlier_probability, exact_fractions);
                }, 0.0, 5.0, exact_mean_bias);
            exact_bias_variance = model_segments_af_argmax(
                [&](const double value) {
                    return exact_log_likelihood(exact_mean_bias, value,
                                                exact_outlier_probability, exact_fractions);
                }, 1.0e-10, 0.5, exact_bias_variance);
            exact_outlier_probability = model_segments_af_argmax(
                [&](const double value) {
                    return exact_log_likelihood(exact_mean_bias, exact_bias_variance,
                                                value, exact_fractions);
                }, 0.0, 0.15, exact_outlier_probability);
            for (std::size_t segment = 0; segment < exact_indices.size(); ++segment) {
                if (exact_indices[segment].empty()) continue;
                exact_fractions[segment] = model_segments_af_argmax(
                    [&](const double value) {
                        auto candidate = exact_fractions;
                        candidate[segment] = value;
                        return exact_log_likelihood(exact_mean_bias, exact_bias_variance,
                                                    exact_outlier_probability, candidate);
                    }, 0.0, 0.5, exact_fractions[segment]);
            }
            const double next_log_likelihood = exact_log_likelihood(
                exact_mean_bias, exact_bias_variance, exact_outlier_probability,
                exact_fractions);
            if (std::isfinite(previous_log_likelihood) &&
                next_log_likelihood - previous_log_likelihood <= 0.5) {
                ++exact_iterations;
                break;
            }
            previous_log_likelihood = current_log_likelihood;
        }
        for (std::size_t segment = 0; segment < segments.size(); ++segment)
            model.segment_means[segment] = exact_fractions[segment];
        model.mean_bias = exact_mean_bias;
        model.bias_variance = exact_bias_variance;
        model.outlier_probability = exact_outlier_probability;
        // Include both the device responsibility pass and the bounded
        // Laplace-coordinate refinement in the reported iteration count.
        // This keeps telemetry monotonic when the Java-style optimizer
        // converges in fewer than four coordinate rounds.
        model.iterations = conditional_iterations + exact_iterations;
    }
    if (!has_count_backed) {
        model.mean_bias = mean_bias;
        model.bias_variance = bias_variance;
        model.outlier_probability = outlier_probability;
        model.iterations = conditional_iterations;
    }
    return model;
}

std::vector<PosteriorQuantiles> sample_model_segments_posterior_bounded(
    const Table& table, const std::vector<Segment>& segments, const Options& options,
    KernelStats* kernel_stats = nullptr,
    const CopyRatioConditionalModel* conditional_model = nullptr,
    const AlleleFractionConditionalModel* allele_conditional_model = nullptr) {
    std::vector<PosteriorQuantiles> result(segments.size());
    if (!options.probabilistic_requested || segments.empty()) return result;
    const bool sample_copy_ratio = options.mode == "BOTH" || options.mode == "COPY_RATIO";
    const bool sample_allele_fraction = options.mode == "BOTH" || options.mode == "ALLELE_FRACTION";
    const std::size_t sample_count = options.num_samples;
    constexpr std::size_t sample_block_size = 4096;
    constexpr std::size_t reservoir_capacity = 2048;
    Kokkos::View<double*> means("model_segments_posterior_copy_means_bounded", segments.size());
    Kokkos::View<double*> deviations("model_segments_posterior_copy_deviations_bounded", segments.size());
    Kokkos::View<double*> allele_means("model_segments_posterior_allele_means_bounded", segments.size());
    Kokkos::View<double*> allele_deviations("model_segments_posterior_allele_deviations_bounded", segments.size());
    Kokkos::View<std::uint64_t*> point_counts("model_segments_posterior_points_bounded", segments.size());
    Kokkos::View<std::uint64_t*> allele_counts("model_segments_posterior_allele_points_bounded", segments.size());
    // Count-backed allele evidence is carried into the posterior kernel so
    // proposals target AlleleFractionLikelihoods rather than only the
    // segment-level MAF Gaussian surrogate.  MAF-only inputs leave these
    // views empty and retain the deterministic legacy fallback.
    std::vector<std::int64_t> posterior_refs;
    std::vector<std::int64_t> posterior_alts;
    std::vector<std::uint64_t> posterior_offsets(segments.size() + 1, 0);
    posterior_refs.reserve(table.points.size());
    posterior_alts.reserve(table.points.size());
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        posterior_offsets[segment] = static_cast<std::uint64_t>(posterior_refs.size());
        for (const auto& point : table.points) {
            if (point.contig != segments[segment].contig || !point.allele_counts_available)
                continue;
            const auto midpoint = point.start + (point.end - point.start) / 2;
            if (midpoint < segments[segment].start || midpoint > segments[segment].end)
                continue;
            posterior_refs.push_back(point.allele_ref_count);
            posterior_alts.push_back(point.allele_alt_count);
        }
    }
    posterior_offsets[segments.size()] = static_cast<std::uint64_t>(posterior_refs.size());
    Kokkos::View<std::int64_t*> posterior_ref_view(
        "model_segments_posterior_allele_refs", posterior_refs.size());
    Kokkos::View<std::int64_t*> posterior_alt_view(
        "model_segments_posterior_allele_alts", posterior_alts.size());
    Kokkos::View<std::uint64_t*> posterior_offset_view(
        "model_segments_posterior_allele_offsets", posterior_offsets.size());
    auto host_posterior_refs = Kokkos::create_mirror_view(posterior_ref_view);
    auto host_posterior_alts = Kokkos::create_mirror_view(posterior_alt_view);
    auto host_posterior_offsets = Kokkos::create_mirror_view(posterior_offset_view);
    for (std::size_t index = 0; index < posterior_refs.size(); ++index) {
        host_posterior_refs(index) = posterior_refs[index];
        host_posterior_alts(index) = posterior_alts[index];
    }
    for (std::size_t index = 0; index < posterior_offsets.size(); ++index)
        host_posterior_offsets(index) = posterior_offsets[index];
    Kokkos::deep_copy(posterior_ref_view, host_posterior_refs);
    Kokkos::deep_copy(posterior_alt_view, host_posterior_alts);
    Kokkos::deep_copy(posterior_offset_view, host_posterior_offsets);
    auto host_means = Kokkos::create_mirror_view(means);
    auto host_deviations = Kokkos::create_mirror_view(deviations);
    auto host_allele_means = Kokkos::create_mirror_view(allele_means);
    auto host_allele_deviations = Kokkos::create_mirror_view(allele_deviations);
    auto host_point_counts = Kokkos::create_mirror_view(point_counts);
    auto host_allele_counts = Kokkos::create_mirror_view(allele_counts);
    for (std::size_t index = 0; index < segments.size(); ++index) {
        const bool has_conditional_mean = conditional_model != nullptr &&
            index < conditional_model->segment_means.size() &&
            std::isfinite(conditional_model->segment_means[index]);
        host_means(index) = has_conditional_mean
            ? conditional_model->segment_means[index] : segments[index].mean;
        host_deviations(index) = segments[index].stddev;
        const bool has_conditional_allele_mean = allele_conditional_model != nullptr &&
            index < allele_conditional_model->segment_means.size() &&
            std::isfinite(allele_conditional_model->segment_means[index]);
        host_allele_means(index) = has_conditional_allele_mean
            ? allele_conditional_model->segment_means[index]
            : segments[index].mean_minor_allele_fraction;
        host_allele_deviations(index) = segments[index].allelic_stddev;
        host_point_counts(index) = segments[index].points;
        host_allele_counts(index) = segments[index].allelic_points;
    }
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::deep_copy(means, host_means);
    Kokkos::deep_copy(deviations, host_deviations);
    Kokkos::deep_copy(allele_means, host_allele_means);
    Kokkos::deep_copy(allele_deviations, host_allele_deviations);
    Kokkos::deep_copy(point_counts, host_point_counts);
    Kokkos::deep_copy(allele_counts, host_allele_counts);

    Kokkos::View<double*> copy_states("model_segments_posterior_copy_states", segments.size());
    Kokkos::View<double*> allele_states("model_segments_posterior_allele_states", segments.size());
    auto host_copy_states = Kokkos::create_mirror_view(copy_states);
    auto host_allele_states = Kokkos::create_mirror_view(allele_states);
    const auto nan_value = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t index = 0; index < segments.size(); ++index) {
        host_copy_states(index) = sample_copy_ratio ? host_means(index) : nan_value;
        const bool valid_allele = sample_allele_fraction &&
            std::isfinite(host_allele_means(index));
        host_allele_states(index) = valid_allele
            ? std::min(0.5, std::max(0.0, host_allele_means(index))) : nan_value;
    }
    Kokkos::deep_copy(copy_states, host_copy_states);
    Kokkos::deep_copy(allele_states, host_allele_states);
    Kokkos::View<double**> copy_samples("model_segments_posterior_copy_block", segments.size(), sample_block_size);
    Kokkos::View<double**> allele_samples("model_segments_posterior_allele_block", segments.size(), sample_block_size);
    auto host_copy_samples = Kokkos::create_mirror_view(copy_samples);
    auto host_allele_samples = Kokkos::create_mirror_view(allele_samples);
    std::vector<std::vector<double>> copy_reservoirs(segments.size());
    std::vector<std::vector<double>> allele_reservoirs(segments.size());
    std::vector<std::size_t> copy_seen(segments.size(), 0);
    std::vector<std::size_t> allele_seen(segments.size(), 0);
    const auto append_reservoir = [](std::vector<double>& reservoir, std::size_t& seen,
                                     const double value, const std::uint64_t seed) {
        if (!std::isfinite(value)) return;
        ++seen;
        if (reservoir.size() < reservoir_capacity) {
            reservoir.push_back(value);
            return;
        }
        const auto slot = model_segments_mix64(seed) % seen;
        if (slot < reservoir_capacity) reservoir[static_cast<std::size_t>(slot)] = value;
    };
    for (std::size_t block_begin = 0; block_begin < sample_count; block_begin += sample_block_size) {
        const auto block_length = std::min(sample_block_size, sample_count - block_begin);
        const auto burn_in = options.num_burn_in_iterations;
        const auto first_block = block_begin == 0;
        const auto block_iterations = first_block ? burn_in + block_length : block_length;
        fastgatk::core::HostBatch host_batch("model-segments-posterior-v1");
        host_batch.records = segments.size() * block_iterations;
        host_batch.bytes = (6 * segments.size() + 2 * segments.size() * block_length) * sizeof(double) +
                           2 * segments.size() * sizeof(std::uint64_t);
        fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-posterior");
        plan.begin_prepare(host_batch);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(segments.size() * block_length);
        device_batch.bind("means", means);
        device_batch.bind("deviations", deviations);
        device_batch.bind("allele_means", allele_means);
        device_batch.bind("allele_deviations", allele_deviations);
        device_batch.bind("point_counts", point_counts);
        device_batch.bind("allele_counts", allele_counts);
        device_batch.bind("posterior_allele_refs", posterior_ref_view);
        device_batch.bind("posterior_allele_alts", posterior_alt_view);
        device_batch.bind("posterior_allele_offsets", posterior_offset_view);
        device_batch.bind("copy_states", copy_states);
        device_batch.bind("allele_states", allele_states);
        device_batch.bind("copy_samples", copy_samples);
        device_batch.bind("allele_samples", allele_samples);
        ExecSpace().fence();
        plan.end_prepare(device_batch);
        plan.begin_execute();
        const double conditional_variance = conditional_model != nullptr &&
            std::isfinite(conditional_model->variance) ? conditional_model->variance : 0.0;
        const double conditional_outlier_probability = conditional_model != nullptr &&
            std::isfinite(conditional_model->outlier_probability)
            ? std::min(0.999, std::max(0.001, conditional_model->outlier_probability)) : 0.0;
        const double conditional_allele_variance = allele_conditional_model != nullptr &&
            std::isfinite(allele_conditional_model->bias_variance)
            ? std::min(0.5, std::max(1.0e-10, 0.25 * allele_conditional_model->bias_variance)) : 0.0;
        const double allele_prior_alpha = std::max(1.0,
            options.minor_allele_fraction_prior_alpha);
        const bool use_exact_allele_likelihood = allele_conditional_model != nullptr &&
            std::isfinite(allele_conditional_model->mean_bias) &&
            std::isfinite(allele_conditional_model->bias_variance) &&
            posterior_refs.size() > 0;
        const double exact_mean_bias = use_exact_allele_likelihood
            ? allele_conditional_model->mean_bias : 1.0;
        const double exact_bias_variance = use_exact_allele_likelihood
            ? allele_conditional_model->bias_variance : 0.05;
        const double exact_outlier_probability = use_exact_allele_likelihood
            ? allele_conditional_model->outlier_probability : 0.01;
        Kokkos::parallel_for("model_segments_posterior_sampler_block",
            Kokkos::RangePolicy<ExecSpace>(0, segments.size()), KOKKOS_LAMBDA(const std::size_t segment) {
                const auto point_count = point_counts(segment) == 0 ? 1ULL : point_counts(segment);
                const auto allele_count = allele_counts(segment);
                const auto conditional_sigma = conditional_variance > 0.0
                    ? Kokkos::sqrt(conditional_variance /
                        (static_cast<double>(point_count) * (1.0 - conditional_outlier_probability))) : 0.0;
                const auto copy_sigma = Kokkos::fmax(
                    Kokkos::fmax(deviations(segment) / Kokkos::sqrt(static_cast<double>(point_count)),
                                 conditional_sigma), 0.005);
                const auto allele_sigma = Kokkos::fmax(
                    Kokkos::fmax(allele_deviations(segment) /
                        Kokkos::sqrt(static_cast<double>(allele_count == 0 ? 1ULL : allele_count)),
                        conditional_allele_variance > 0.0
                            ? Kokkos::sqrt(conditional_allele_variance /
                                static_cast<double>(allele_count == 0 ? 1ULL : allele_count)) : 0.0),
                    0.005);
                double copy_state = copy_states(segment);
                double allele_state = allele_states(segment);
                for (std::size_t local_iteration = 0; local_iteration < block_iterations; ++local_iteration) {
                    const auto sample_index = first_block
                        ? (local_iteration >= burn_in ? local_iteration - burn_in : sample_count)
                        : block_begin + local_iteration;
                    const auto global_iteration = first_block
                        ? local_iteration : burn_in + block_begin + local_iteration;
                    const auto seed = 1216ULL ^ (static_cast<std::uint64_t>(segment) * 0x9e3779b97f4a7c15ULL) ^
                                      (static_cast<std::uint64_t>(global_iteration) * 0xd1b54a32d192ed03ULL);
                    if (sample_copy_ratio) {
                        const auto proposal = copy_state + 0.75 * copy_sigma * model_segments_normal(seed);
                        const auto old_density = model_segments_log_density(copy_state, means(segment), copy_sigma);
                        const auto new_density = model_segments_log_density(proposal, means(segment), copy_sigma);
                        const auto log_acceptance = Kokkos::fmin(0.0, new_density - old_density);
                        if (Kokkos::log(model_segments_uniform(seed ^ 0x94d049bb133111ebULL)) < log_acceptance)
                            copy_state = proposal;
                    }
                    if (sample_allele_fraction && allele_count > 0 && allele_state == allele_state) {
                        auto proposal = allele_state + 0.75 * allele_sigma * model_segments_normal(seed ^ 0xbf58476d1ce4e5b9ULL);
                        if (proposal < 0.0) proposal = -proposal;
                        if (proposal > 0.5) proposal = 1.0 - proposal;
                        proposal = Kokkos::fmin(0.5, Kokkos::fmax(0.0, proposal));
                        double old_density = model_segments_log_density(
                            allele_state, allele_means(segment), allele_sigma) +
                            model_segments_af_prior_log_density(allele_state, allele_prior_alpha);
                        double new_density = model_segments_log_density(
                            proposal, allele_means(segment), allele_sigma) +
                            model_segments_af_prior_log_density(proposal, allele_prior_alpha);
                        if (use_exact_allele_likelihood) {
                            old_density = model_segments_af_prior_log_density(
                                allele_state, allele_prior_alpha);
                            new_density = model_segments_af_prior_log_density(
                                proposal, allele_prior_alpha);
                            const auto evidence_begin = posterior_offset_view(segment);
                            const auto evidence_end = posterior_offset_view(segment + 1);
                            for (auto evidence = evidence_begin; evidence < evidence_end; ++evidence) {
                                old_density += model_segments_af_het_log_likelihood(
                                    exact_mean_bias, exact_bias_variance,
                                    exact_outlier_probability, allele_state,
                                    posterior_alt_view(evidence), posterior_ref_view(evidence));
                                new_density += model_segments_af_het_log_likelihood(
                                    exact_mean_bias, exact_bias_variance,
                                    exact_outlier_probability, proposal,
                                    posterior_alt_view(evidence), posterior_ref_view(evidence));
                            }
                        }
                        const auto log_acceptance = Kokkos::fmin(0.0, new_density - old_density);
                        if (Kokkos::log(model_segments_uniform(seed ^ 0x94d049bb133111ebULL)) < log_acceptance)
                            allele_state = proposal;
                    }
                    if (sample_index < sample_count) {
                        copy_samples(segment, sample_index - block_begin) = copy_state;
                        allele_samples(segment, sample_index - block_begin) = allele_state;
                    }
                }
                copy_states(segment) = copy_state;
                allele_states(segment) = allele_state;
            });
        ExecSpace().fence();
        plan.end_execute();
        if (kernel_stats != nullptr) {
            ++kernel_stats->batches;
            kernel_stats->observations += segments.size() * block_iterations;
            kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
            kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
        }
        Kokkos::deep_copy(host_copy_samples, copy_samples);
        Kokkos::deep_copy(host_allele_samples, allele_samples);
        for (std::size_t segment = 0; segment < segments.size(); ++segment) {
            for (std::size_t local_sample = 0; local_sample < block_length; ++local_sample) {
                const auto global_sample = block_begin + local_sample;
                append_reservoir(copy_reservoirs[segment], copy_seen[segment], host_copy_samples(segment, local_sample),
                    0x6a09e667f3bcc909ULL ^ (static_cast<std::uint64_t>(segment) << 32) ^ global_sample);
                append_reservoir(allele_reservoirs[segment], allele_seen[segment], host_allele_samples(segment, local_sample),
                    0xbb67ae8584caa73bULL ^ (static_cast<std::uint64_t>(segment) << 32) ^ global_sample);
            }
        }
    }
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        result[segment].copy_ratio_10 = model_segments_quantile(copy_reservoirs[segment], 0.10);
        result[segment].copy_ratio_50 = model_segments_quantile(copy_reservoirs[segment], 0.50);
        result[segment].copy_ratio_90 = model_segments_quantile(copy_reservoirs[segment], 0.90);
        result[segment].allele_fraction_10 = model_segments_quantile(allele_reservoirs[segment], 0.10);
        result[segment].allele_fraction_50 = model_segments_quantile(allele_reservoirs[segment], 0.50);
        result[segment].allele_fraction_90 = model_segments_quantile(allele_reservoirs[segment], 0.90);
    }
    return result;
}

// Kept only as an historical reference; the bounded implementation below is
// the sole production posterior path and owns the Kokkos lifecycle.
#if 0
[[maybe_unused]] std::vector<PosteriorQuantiles> sample_model_segments_posterior_legacy(
    const std::vector<Segment>& segments, const Options& options) {
    std::vector<PosteriorQuantiles> result(segments.size());
    if (!options.probabilistic_requested || segments.empty()) return result;

    const bool sample_copy_ratio = options.mode == "BOTH" || options.mode == "COPY_RATIO";
    const bool sample_allele_fraction = options.mode == "BOTH" || options.mode == "ALLELE_FRACTION";
    const std::size_t sample_count = options.num_samples;
    const std::size_t total_iterations = options.num_burn_in_iterations + sample_count;
    Kokkos::View<double*> means("model_segments_posterior_copy_means", segments.size());
    Kokkos::View<double*> deviations("model_segments_posterior_copy_deviations", segments.size());
    Kokkos::View<double*> allele_means("model_segments_posterior_allele_means", segments.size());
    Kokkos::View<double*> allele_deviations("model_segments_posterior_allele_deviations", segments.size());
    Kokkos::View<std::uint64_t*> point_counts("model_segments_posterior_points", segments.size());
    Kokkos::View<std::uint64_t*> allele_counts("model_segments_posterior_allele_points", segments.size());
    auto host_means = Kokkos::create_mirror_view(means);
    auto host_deviations = Kokkos::create_mirror_view(deviations);
    auto host_allele_means = Kokkos::create_mirror_view(allele_means);
    auto host_allele_deviations = Kokkos::create_mirror_view(allele_deviations);
    auto host_point_counts = Kokkos::create_mirror_view(point_counts);
    auto host_allele_counts = Kokkos::create_mirror_view(allele_counts);
    for (std::size_t index = 0; index < segments.size(); ++index) {
        host_means(index) = segments[index].mean;
        host_deviations(index) = segments[index].stddev;
        host_allele_means(index) = segments[index].mean_minor_allele_fraction;
        host_allele_deviations(index) = segments[index].allelic_stddev;
        host_point_counts(index) = segments[index].points;
        host_allele_counts(index) = segments[index].allelic_points;
    }
    Kokkos::deep_copy(means, host_means);
    Kokkos::deep_copy(deviations, host_deviations);
    Kokkos::deep_copy(allele_means, host_allele_means);
    Kokkos::deep_copy(allele_deviations, host_allele_deviations);
    Kokkos::deep_copy(point_counts, host_point_counts);
    Kokkos::deep_copy(allele_counts, host_allele_counts);

    Kokkos::View<double***> samples("model_segments_posterior_samples",
                                    segments.size(), sample_count, 2);
    const auto nan_value = std::numeric_limits<double>::quiet_NaN();
    Kokkos::parallel_for(
        "model_segments_posterior_sampler",
        Kokkos::RangePolicy<>(0, segments.size()),
        KOKKOS_LAMBDA(const std::size_t segment) {
            const auto point_count = point_counts(segment) == 0 ? 1ULL : point_counts(segment);
            const auto allele_count = allele_counts(segment);
            const auto copy_sigma = Kokkos::fmax(
                deviations(segment) / Kokkos::sqrt(static_cast<double>(point_count)), 0.005);
            const auto allele_sigma = Kokkos::fmax(
                allele_deviations(segment) / Kokkos::sqrt(static_cast<double>(allele_count == 0 ? 1ULL : allele_count)),
                0.005);
            double copy_state = means(segment);
            double allele_state = Kokkos::fmin(0.5, Kokkos::fmax(0.0, allele_means(segment)));
            if (!sample_copy_ratio) copy_state = nan_value;
            if (!sample_allele_fraction || allele_count == 0 || allele_means(segment) != allele_means(segment))
                allele_state = nan_value;
            for (std::size_t iteration = 0; iteration < total_iterations; ++iteration) {
                const auto seed = 1216ULL ^ (static_cast<std::uint64_t>(segment) * 0x9e3779b97f4a7c15ULL) ^
                                  (static_cast<std::uint64_t>(iteration) * 0xd1b54a32d192ed03ULL);
                if (sample_copy_ratio) {
                    const auto proposal = copy_state + 0.75 * copy_sigma * model_segments_normal(seed);
                    const auto old_density = model_segments_log_density(copy_state, means(segment), copy_sigma);
                    const auto new_density = model_segments_log_density(proposal, means(segment), copy_sigma);
                    const auto log_acceptance = Kokkos::fmin(0.0, new_density - old_density);
                    if (Kokkos::log(model_segments_uniform(seed ^ 0x94d049bb133111ebULL)) < log_acceptance)
                        copy_state = proposal;
                }
                if (sample_allele_fraction && allele_count > 0 && allele_state == allele_state) {
                    auto proposal = allele_state + 0.75 * allele_sigma *
                        model_segments_normal(seed ^ 0xbf58476d1ce4e5b9ULL);
                    // Reflect proposals at the physical MAF boundaries rather
                    // than truncating them, preserving a symmetric proposal.
                    if (proposal < 0.0) proposal = -proposal;
                    if (proposal > 0.5) proposal = 1.0 - proposal;
                    proposal = Kokkos::fmin(0.5, Kokkos::fmax(0.0, proposal));
                    const auto old_density = model_segments_log_density(allele_state,
                        allele_means(segment), allele_sigma);
                    const auto new_density = model_segments_log_density(proposal,
                        allele_means(segment), allele_sigma);
                    const auto log_acceptance = Kokkos::fmin(0.0, new_density - old_density);
                    if (Kokkos::log(model_segments_uniform(seed ^ 0x94d049bb133111ebULL)) < log_acceptance)
                        allele_state = proposal;
                }
                if (iteration >= options.num_burn_in_iterations) {
                    const auto sample = iteration - options.num_burn_in_iterations;
                    samples(segment, sample, 0) = copy_state;
                    samples(segment, sample, 1) = allele_state;
                }
            }
        });
    Kokkos::fence();
    auto host_samples = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, samples);
    for (std::size_t segment = 0; segment < segments.size(); ++segment) {
        std::vector<double> copy_values(sample_count);
        std::vector<double> allele_values(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            copy_values[sample] = host_samples(segment, sample, 0);
            allele_values[sample] = host_samples(segment, sample, 1);
        }
        result[segment].copy_ratio_10 = model_segments_quantile(std::move(copy_values), 0.10);
        // Recreate the vector for each quantile so the helper owns sorting and
        // filtering; sample counts are bounded by the CLI validation above.
        copy_values.resize(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample)
            copy_values[sample] = host_samples(segment, sample, 0);
        result[segment].copy_ratio_50 = model_segments_quantile(copy_values, 0.50);
        result[segment].copy_ratio_90 = model_segments_quantile(std::move(copy_values), 0.90);
        result[segment].allele_fraction_10 = model_segments_quantile(std::move(allele_values), 0.10);
        allele_values.resize(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample)
            allele_values[sample] = host_samples(segment, sample, 1);
        result[segment].allele_fraction_50 = model_segments_quantile(allele_values, 0.50);
        result[segment].allele_fraction_90 = model_segments_quantile(std::move(allele_values), 0.90);
    }
    return result;
}
#endif

std::vector<PosteriorQuantiles> sample_model_segments_posterior(
    const Table& table, const std::vector<Segment>& segments, const Options& options,
    KernelStats* kernel_stats = nullptr,
    CopyRatioConditionalModel* conditional_model = nullptr,
    AlleleFractionConditionalModel* allele_conditional_model = nullptr) {
    CopyRatioConditionalModel fitted = fit_copy_ratio_conditional_model(
        table, segments, options, kernel_stats);
    AlleleFractionConditionalModel fitted_allele = fit_allele_fraction_conditional_model(
        table, segments, options, kernel_stats);
    if (conditional_model != nullptr) *conditional_model = fitted;
    if (allele_conditional_model != nullptr) *allele_conditional_model = fitted_allele;
    return sample_model_segments_posterior_bounded(
        table, segments, options, kernel_stats,
        fitted.iterations == 0 ? nullptr : &fitted,
        fitted_allele.iterations == 0 ? nullptr : &fitted_allele);
}

struct SmoothingResult {
    bool applied = false;
    std::size_t iterations = 0;
};

// This is the deterministic equivalent of GATK's
// MultidimensionalModeller.SimilarSegmentUtils.  The Java implementation
// merges adjacent segments when either posterior's 10--90% width contains the
// median difference, then combines the two posteriors as normal distributions
// using inverse-variance weighting.  Keep the candidate calculation in a
// Kokkos RangePolicy so the same decision kernel is available to Serial,
// OpenMP, CUDA, HIP, and SYCL; only the ordered list mutation is collected on
// Host after the kernel completes.
SmoothingResult smooth_model_segments(
    std::vector<Segment>& segments,
    std::vector<PosteriorQuantiles>& posterior,
    const Options& options,
    KernelStats* kernel_stats = nullptr) {
    SmoothingResult result;
    if (!options.probabilistic_requested || !options.smoothing_requested ||
        options.maximum_number_of_smoothing_iterations == 0 || segments.size() < 2)
        return result;

    auto merge_summary = [](const double low1, const double median1, const double high1,
                            const double low2, const double median2, const double high2) {
        if (!(median1 == median1)) {
            if (median2 == median2) return std::array<double, 3>{low2, median2, high2};
            return std::array<double, 3>{low1, median1, high1};
        }
        if (!(median2 == median2)) return std::array<double, 3>{low1, median1, high1};
        const double sd1 = 0.5 * (high1 - low1);
        const double sd2 = 0.5 * (high2 - low2);
        if (std::isfinite(sd1) && std::isfinite(sd2) && sd1 > 0.0 && sd2 > 0.0) {
            const double variance = 1.0 / (1.0 / (sd1 * sd1) + 1.0 / (sd2 * sd2));
            const double mean = (median1 / (sd1 * sd1) + median2 / (sd2 * sd2)) * variance;
            const double sd = std::sqrt(variance);
            return std::array<double, 3>{mean - sd, mean, mean + sd};
        }
        // The Java MCMC normally supplies positive credible widths.  Keep a
        // finite deterministic fallback for a degenerate bounded chain.
        const double mean = 0.5 * (median1 + median2);
        const double width = std::max(std::abs(high1 - low1), std::abs(high2 - low2)) * 0.25;
        return std::array<double, 3>{mean - width, mean, mean + width};
    };

    auto merge_segment = [](const Segment& left, const Segment& right) {
        Segment merged = left;
        merged.start = std::min(left.start, right.start);
        merged.end = std::max(left.end, right.end);
        const std::size_t points = left.points + right.points;
        if (points > 0) {
            merged.mean = (left.mean * static_cast<double>(left.points) +
                           right.mean * static_cast<double>(right.points)) /
                          static_cast<double>(points);
            const double left_second = std::isfinite(left.stddev)
                ? left.stddev * left.stddev + left.mean * left.mean : left.mean * left.mean;
            const double right_second = std::isfinite(right.stddev)
                ? right.stddev * right.stddev + right.mean * right.mean : right.mean * right.mean;
            merged.stddev = std::sqrt(std::max(0.0,
                (left_second * static_cast<double>(left.points) +
                 right_second * static_cast<double>(right.points)) /
                    static_cast<double>(points) - merged.mean * merged.mean));
        }
        merged.points = points;
        const std::size_t allele_points = left.allelic_points + right.allelic_points;
        if (allele_points > 0 && std::isfinite(left.mean_minor_allele_fraction) &&
            std::isfinite(right.mean_minor_allele_fraction)) {
            merged.mean_minor_allele_fraction =
                (left.mean_minor_allele_fraction * static_cast<double>(left.allelic_points) +
                 right.mean_minor_allele_fraction * static_cast<double>(right.allelic_points)) /
                static_cast<double>(allele_points);
            const double left_second = left.allelic_stddev * left.allelic_stddev +
                                       left.mean_minor_allele_fraction * left.mean_minor_allele_fraction;
            const double right_second = right.allelic_stddev * right.allelic_stddev +
                                        right.mean_minor_allele_fraction * right.mean_minor_allele_fraction;
            merged.allelic_stddev = std::sqrt(std::max(0.0,
                (left_second * static_cast<double>(left.allelic_points) +
                 right_second * static_cast<double>(right.allelic_points)) /
                    static_cast<double>(allele_points) -
                merged.mean_minor_allele_fraction * merged.mean_minor_allele_fraction));
        } else if (allele_points > 0) {
            merged.mean_minor_allele_fraction = std::isfinite(left.mean_minor_allele_fraction)
                ? left.mean_minor_allele_fraction : right.mean_minor_allele_fraction;
        } else {
            merged.mean_minor_allele_fraction = std::numeric_limits<double>::quiet_NaN();
            merged.allelic_stddev = 0.0;
        }
        merged.allelic_points = allele_points;
        return merged;
    };

    for (std::size_t iteration = 0;
         iteration < options.maximum_number_of_smoothing_iterations;
         ++iteration) {
        bool merged_in_iteration = false;
        // Java stays on the merged segment after a merge, allowing a chain of
        // similar neighbors to collapse in the same smoothing iteration.
        while (segments.size() > 1) {
            const std::size_t candidates = segments.size() - 1;
            Kokkos::View<double**> metrics("model_segments_smoothing_metrics", candidates, 7);
            Kokkos::View<std::uint8_t*> merge_flags("model_segments_smoothing_flags", candidates);
            auto host_metrics = Kokkos::create_mirror_view(metrics);
            for (std::size_t index = 0; index < candidates; ++index) {
                const auto& left = posterior[index];
                const auto& right = posterior[index + 1];
                host_metrics(index, 0) = segments[index].contig == segments[index + 1].contig ? 1.0 : 0.0;
                host_metrics(index, 1) = (left.copy_ratio_50 == left.copy_ratio_50 &&
                                         right.copy_ratio_50 == right.copy_ratio_50)
                    ? std::abs(left.copy_ratio_50 - right.copy_ratio_50) : -1.0;
                host_metrics(index, 2) = left.copy_ratio_90 - left.copy_ratio_10;
                host_metrics(index, 3) = right.copy_ratio_90 - right.copy_ratio_10;
                host_metrics(index, 4) = (left.allele_fraction_50 == left.allele_fraction_50 &&
                                         right.allele_fraction_50 == right.allele_fraction_50)
                    ? std::abs(left.allele_fraction_50 - right.allele_fraction_50) : -1.0;
                host_metrics(index, 5) = left.allele_fraction_90 - left.allele_fraction_10;
                host_metrics(index, 6) = right.allele_fraction_90 - right.allele_fraction_10;
            }
            using ExecSpace = Kokkos::DefaultExecutionSpace;
            fastgatk::core::HostBatch host_batch("model-segments-smoothing-v1");
            host_batch.records = candidates;
            host_batch.bytes = candidates * (7 * sizeof(double) + sizeof(std::uint8_t));
            fastgatk::core::KernelPlan<ExecSpace> plan("model-segments-smoothing");
            plan.begin_prepare(host_batch);
            Kokkos::deep_copy(metrics, host_metrics);
            fastgatk::core::DeviceBatch<ExecSpace> device_batch(candidates);
            device_batch.bind("metrics", metrics);
            device_batch.bind("merge_flags", merge_flags);
            ExecSpace().fence();
            plan.end_prepare(device_batch);
            const double copy_threshold = options.smoothing_credible_interval_threshold_copy_ratio;
            const double allele_threshold = options.smoothing_credible_interval_threshold_allele_fraction;
            plan.begin_execute();
            Kokkos::parallel_for("model_segments_smoothing_candidates",
                Kokkos::RangePolicy<ExecSpace>(0, candidates), KOKKOS_LAMBDA(const std::size_t index) {
                    const bool copy_similar = metrics(index, 1) < 0.0 ||
                        metrics(index, 1) < copy_threshold * metrics(index, 2) ||
                        metrics(index, 1) < copy_threshold * metrics(index, 3);
                    const bool allele_similar = metrics(index, 4) < 0.0 ||
                        metrics(index, 4) < allele_threshold * metrics(index, 5) ||
                        metrics(index, 4) < allele_threshold * metrics(index, 6);
                    merge_flags(index) = metrics(index, 0) > 0.5 && copy_similar && allele_similar;
                });
            ExecSpace().fence();
            plan.end_execute();
            if (kernel_stats != nullptr) {
                ++kernel_stats->batches;
                kernel_stats->observations += candidates;
                kernel_stats->prepare_seconds += plan.telemetry().prepare_seconds;
                kernel_stats->execute_seconds += plan.telemetry().execute_seconds;
            }
            auto host_flags = Kokkos::create_mirror_view(merge_flags);
            Kokkos::deep_copy(host_flags, merge_flags);
            std::size_t merge_index = candidates;
            for (std::size_t index = 0; index < candidates; ++index) {
                if (host_flags(index) != 0) { merge_index = index; break; }
            }
            if (merge_index == candidates) break;
            segments[merge_index] = merge_segment(segments[merge_index], segments[merge_index + 1]);
            const auto& left = posterior[merge_index];
            const auto& right = posterior[merge_index + 1];
            PosteriorQuantiles merged;
            const auto copy = merge_summary(left.copy_ratio_10, left.copy_ratio_50, left.copy_ratio_90,
                                            right.copy_ratio_10, right.copy_ratio_50, right.copy_ratio_90);
            merged.copy_ratio_10 = copy[0];
            merged.copy_ratio_50 = copy[1];
            merged.copy_ratio_90 = copy[2];
            const auto allele = merge_summary(left.allele_fraction_10, left.allele_fraction_50,
                                              left.allele_fraction_90, right.allele_fraction_10,
                                              right.allele_fraction_50, right.allele_fraction_90);
            merged.allele_fraction_10 = allele[0];
            merged.allele_fraction_50 = allele[1];
            merged.allele_fraction_90 = allele[2];
            posterior[merge_index] = merged;
            segments.erase(segments.begin() + static_cast<std::ptrdiff_t>(merge_index + 1));
            posterior.erase(posterior.begin() + static_cast<std::ptrdiff_t>(merge_index + 1));
            merged_in_iteration = true;
        }
        if (!merged_in_iteration) break;
        result.applied = true;
        ++result.iterations;
    }
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

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

struct ModelSegmentsSidecars {
    std::vector<std::string> paths;
};

ModelSegmentsSidecars write_model_segments_sidecars(
    const Options& options, const Table& table, const std::vector<Segment>& segments,
    const std::vector<PosteriorQuantiles>& posterior, const bool allelic_signal_used,
    const CopyRatioConditionalModel* conditional_model = nullptr,
    const AlleleFractionConditionalModel* allele_conditional_model = nullptr) {
    ModelSegmentsSidecars sidecars;
    if (options.output_prefix.empty()) return sidecars;
    const auto prefix = options.output_directory.empty()
        ? std::filesystem::path(options.output_prefix)
        : std::filesystem::path(options.output_directory) /
          std::filesystem::path(options.output_prefix).filename();
    const auto modeled_path = prefix.string() + ".modelFinal.seg";
    const auto begin_path = prefix.string() + ".modelBegin.seg";
    const auto copy_ratio_path = prefix.string() + ".cr.seg";
    const auto copy_ratio_igv_path = prefix.string() + ".cr.igv.seg";
    const auto allele_fraction_igv_path = prefix.string() + ".af.igv.seg";
    const auto begin_copy_parameter_path = prefix.string() + ".modelBegin.cr.param";
    const auto begin_allele_parameter_path = prefix.string() + ".modelBegin.af.param";
    const auto final_copy_parameter_path = prefix.string() + ".modelFinal.cr.param";
    const auto final_allele_parameter_path = prefix.string() + ".modelFinal.af.param";
    const auto sample = options.sample.empty() ? table.sample : options.sample;

    auto write_modeled = [&](const std::string& path) {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ModelSegments sidecar: " + path);
        output << "@HD\tVN:1.6\n";
        for (std::size_t i = 0; i < table.contigs.size(); ++i)
            output << "@SQ\tSN:" << table.contigs[i] << "\tLN:" << table.contig_lengths[i] << '\n';
        output << "@RG\tID:GATKCopyNumber\tSM:" << sample << '\n'
               << "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tNUM_POINTS_ALLELE_FRACTION"
               << "\tLOG2_COPY_RATIO_POSTERIOR_10\tLOG2_COPY_RATIO_POSTERIOR_50\tLOG2_COPY_RATIO_POSTERIOR_90"
               << "\tMINOR_ALLELE_FRACTION_POSTERIOR_10\tMINOR_ALLELE_FRACTION_POSTERIOR_50\tMINOR_ALLELE_FRACTION_POSTERIOR_90\n";
        output << std::fixed << std::setprecision(6);
        for (std::size_t index = 0; index < segments.size(); ++index) {
            const auto& segment = segments[index];
            const auto& quantiles = posterior[index];
            const double copy10 = options.probabilistic_requested ? quantiles.copy_ratio_10 : segment.mean - segment.stddev;
            const double copy50 = options.probabilistic_requested ? quantiles.copy_ratio_50 : segment.mean;
            const double copy90 = options.probabilistic_requested ? quantiles.copy_ratio_90 : segment.mean + segment.stddev;
            output << segment.contig << '\t' << segment.start << '\t' << segment.end << '\t'
                   << segment.points << '\t' << segment.allelic_points << '\t'
                   << copy10 << '\t' << copy50 << '\t' << copy90;
            if (!allelic_signal_used || segment.allelic_points == 0 ||
                !std::isfinite(segment.mean_minor_allele_fraction)) {
                output << "\tNaN\tNaN\tNaN";
            } else {
                const double allele10 = options.probabilistic_requested
                    ? quantiles.allele_fraction_10 : segment.mean_minor_allele_fraction - segment.allelic_stddev;
                const double allele50 = options.probabilistic_requested
                    ? quantiles.allele_fraction_50 : segment.mean_minor_allele_fraction;
                const double allele90 = options.probabilistic_requested
                    ? quantiles.allele_fraction_90 : segment.mean_minor_allele_fraction + segment.allelic_stddev;
                output << '\t' << allele10 << '\t' << allele50 << '\t' << allele90;
            }
            output << '\n';
        }
        output.close();
        if (!output || !file_complete(path))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete ModelSegments sidecar: " + path);
    };
    write_modeled(modeled_path);
    write_modeled(begin_path);

    // GATK always materializes global parameter reports, even for a
    // copy-ratio-only run (the allele-fraction report then contains the
    // initialized prior parameters).  The native bounded model does not yet
    // expose the full Gibbs state, but emitting the stable GATK report schema
    // makes the output directly consumable by downstream Picard/Nextflow
    // stages and records conservative deterministic summaries.
    const auto empirical_copy_variance = [&]() {
        double sum = 0.0;
        std::size_t count = 0;
        for (const auto& segment : segments) {
            if (std::isfinite(segment.stddev)) {
                sum += segment.stddev * segment.stddev;
                ++count;
            }
        }
        return count == 0 ? 1.0e-6 : std::max(1.0e-6, sum / static_cast<double>(count));
    }();
    const auto copy_variance = conditional_model != nullptr &&
        std::isfinite(conditional_model->variance)
        ? std::max(1.0e-6, conditional_model->variance) : empirical_copy_variance;
    const auto copy_outlier_probability = conditional_model != nullptr &&
        std::isfinite(conditional_model->outlier_probability)
        ? std::min(0.999, std::max(0.001, conditional_model->outlier_probability)) : 0.05;
    const auto write_parameter_table = [&](const std::string& path, const bool allele_fraction) {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ModelSegments parameter sidecar: " + path);
        output << "@HD\tVN:1.6\n@RG\tID:GATKCopyNumber\tSM:" << sample << '\n'
               << "PARAMETER_NAME\tPOSTERIOR_10\tPOSTERIOR_20\tPOSTERIOR_30\tPOSTERIOR_40\tPOSTERIOR_50\tPOSTERIOR_60\tPOSTERIOR_70\tPOSTERIOR_80\tPOSTERIOR_90\n";
        output << std::fixed << std::setprecision(6);
        const auto emit = [&](const std::string& name, const double value) {
            output << name;
            for (int decile = 0; decile < 9; ++decile) output << '\t' << value;
            output << '\n';
        };
        if (allele_fraction) {
            const double mean_bias = allele_conditional_model != nullptr &&
                std::isfinite(allele_conditional_model->mean_bias)
                ? allele_conditional_model->mean_bias : 1.0;
            const double bias_variance = allele_conditional_model != nullptr &&
                std::isfinite(allele_conditional_model->bias_variance)
                ? allele_conditional_model->bias_variance : 0.05;
            const double outlier_probability = allele_conditional_model != nullptr &&
                std::isfinite(allele_conditional_model->outlier_probability)
                ? allele_conditional_model->outlier_probability : 0.01;
            emit("MEAN_BIAS", mean_bias);
            emit("BIAS_VARIANCE", bias_variance);
            emit("OUTLIER_PROBABILITY", outlier_probability);
        } else {
            emit("VARIANCE", copy_variance);
            emit("OUTLIER_PROBABILITY", copy_outlier_probability);
        }
        output.close();
        if (!output || !file_complete(path))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete ModelSegments parameter sidecar: " + path);
    };
    write_parameter_table(begin_copy_parameter_path, false);
    write_parameter_table(begin_allele_parameter_path, true);
    write_parameter_table(final_copy_parameter_path, false);
    write_parameter_table(final_allele_parameter_path, true);

    const auto copy_ratio_median = [&](std::size_t index) {
        return options.probabilistic_requested ? posterior[index].copy_ratio_50 : segments[index].mean;
    };
    const auto allele_fraction_median = [&](std::size_t index) {
        return options.probabilistic_requested ? posterior[index].allele_fraction_50 :
            segments[index].mean_minor_allele_fraction;
    };
    {
        std::ofstream output(copy_ratio_path);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ModelSegments sidecar: " + copy_ratio_path);
        output << "@HD\tVN:1.6\n@RG\tID:GATKCopyNumber\tSM:" << sample << '\n'
               << "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n";
        output << std::fixed << std::setprecision(6);
        for (std::size_t index = 0; index < segments.size(); ++index)
            output << segments[index].contig << '\t' << segments[index].start << '\t' << segments[index].end << '\t'
                   << segments[index].points << '\t' << copy_ratio_median(index) << '\n';
    }
    {
        std::ofstream output(copy_ratio_igv_path);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ModelSegments sidecar: " + copy_ratio_igv_path);
        output << "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n";
        output << std::fixed << std::setprecision(6);
        for (std::size_t index = 0; index < segments.size(); ++index)
            output << sample << '\t' << segments[index].contig << '\t' << segments[index].start << '\t'
                   << segments[index].end << '\t' << segments[index].points << '\t' << copy_ratio_median(index) << '\n';
    }
    {
        std::ofstream output(allele_fraction_igv_path);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ModelSegments sidecar: " + allele_fraction_igv_path);
        output << "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n";
        output << std::fixed << std::setprecision(6);
        for (std::size_t index = 0; index < segments.size(); ++index)
            output << sample << '\t' << segments[index].contig << '\t' << segments[index].start << '\t'
                   << segments[index].end << '\t' << segments[index].allelic_points << '\t'
                   << allele_fraction_median(index) << '\n';
    }
    for (const auto& path : {modeled_path, begin_path, copy_ratio_path,
                             copy_ratio_igv_path, allele_fraction_igv_path,
                             begin_copy_parameter_path, begin_allele_parameter_path,
                             final_copy_parameter_path, final_allele_parameter_path}) {
        if (!file_complete(path))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete ModelSegments sidecar: " + path);
        sidecars.paths.push_back(path);
    }
    return sidecars;
}

std::string point_key(const Point& point) {
    return point.contig + ":" + std::to_string(point.start) + ":" + std::to_string(point.end);
}

std::string multisample_interval_output(const Options& options) {
    if (options.output_prefix.empty())
        throw std::invalid_argument("--output-prefix is required for multi-sample ModelSegments");
    const auto prefix = options.output_directory.empty()
        ? std::filesystem::path(options.output_prefix)
        : std::filesystem::path(options.output_directory) /
          std::filesystem::path(options.output_prefix).filename();
    return prefix.string() + ".interval_list";
}

int run_multisample_model_segments(const Options& options, KernelStats* kernel_stats = nullptr) {
    const auto wall_start = std::chrono::steady_clock::now();
    if (!options.input_segments.empty())
        throw std::invalid_argument("--segments cannot be specified in multisample mode");
    if (options.denoised_copy_ratios.size() < 2)
        throw std::invalid_argument("multi-sample ModelSegments requires at least two denoised-copy-ratios files");
    if (!options.allelic_count_files.empty() &&
        options.allelic_count_files.size() != options.denoised_copy_ratios.size())
        throw std::invalid_argument("multi-sample ModelSegments requires one --allelic-counts file per denoised-copy-ratios file");

    std::vector<Table> tables;
    tables.reserve(options.denoised_copy_ratios.size());
    for (const auto& path : options.denoised_copy_ratios) tables.push_back(read_copy_ratios(path));
    std::vector<AllelicCountMap> prepared_case_counts;
    const AllelicCountMap* matched_normal_hets = nullptr;
    AllelicCountMap normal_hets;
    if (!options.normal_allelic_counts.empty()) {
        const auto normal_counts = read_allelic_counts(options.normal_allelic_counts);
        normal_hets = filter_heterozygous_counts(
            normal_counts, &tables.front(), options.min_total_allele_count_normal,
            options.genotyping_homozygous_log_ratio_threshold,
            options.genotyping_base_error_rate, kernel_stats);
        matched_normal_hets = &normal_hets;
    }
    if (!options.allelic_count_files.empty()) {
        prepared_case_counts.reserve(tables.size());
        for (std::size_t sample = 0; sample < tables.size(); ++sample) {
            const auto raw_counts = read_allelic_counts(options.allelic_count_files[sample]);
            prepared_case_counts.push_back(prepare_case_heterozygous_counts(
                tables[sample], raw_counts, matched_normal_hets, options, kernel_stats));
        }
        if (matched_normal_hets == nullptr && prepared_case_counts.size() > 1) {
            AllelicCountMap common = prepared_case_counts.front();
            for (std::size_t sample = 1; sample < prepared_case_counts.size(); ++sample)
                common = intersect_counts(common, prepared_case_counts[sample]);
            for (auto& counts : prepared_case_counts) counts = intersect_counts(counts, common);
        }
        for (std::size_t sample = 0; sample < tables.size(); ++sample)
            apply_allelic_counts_to_table(tables[sample], prepared_case_counts[sample]);
    }

    // GATK's multidimensional segmenter uses common loci across samples.  The
    // bounded native path materializes a deterministic joint evidence table:
    // copy ratio and MAF are averaged over the samples that contain the exact
    // interval key, while Host retains the original dictionary/coordinates.
    std::vector<std::unordered_map<std::string, const Point*>> indexes(tables.size());
    for (std::size_t sample = 0; sample < tables.size(); ++sample)
        for (const auto& point : tables[sample].points)
            indexes[sample].emplace(point_key(point), &point);
    Table joint = tables.front();
    joint.sample = "MULTISAMPLE";
    joint.points.clear();
    for (const auto& anchor : tables.front().points) {
        const auto key = point_key(anchor);
        double copy_sum = 0.0;
        std::size_t copy_count = 0;
        double allele_sum = 0.0;
        std::size_t allele_count = 0;
        double segmentation_allele_sum = 0.0;
        std::size_t segmentation_allele_count = 0;
        for (const auto& index : indexes) {
            const auto iterator = index.find(key);
            if (iterator == index.end()) continue;
            copy_sum += iterator->second->value;
            ++copy_count;
            if (std::isfinite(iterator->second->minor_allele_fraction)) {
                allele_sum += iterator->second->minor_allele_fraction;
                ++allele_count;
            }
            if (std::isfinite(iterator->second->segmentation_alt_fraction)) {
                segmentation_allele_sum += iterator->second->segmentation_alt_fraction;
                ++segmentation_allele_count;
            }
        }
        if (copy_count != tables.size()) continue;
        Point point = anchor;
        point.value = copy_sum / static_cast<double>(copy_count);
        point.minor_allele_fraction = allele_count == 0
            ? std::numeric_limits<double>::quiet_NaN()
            : allele_sum / static_cast<double>(allele_count);
        point.segmentation_alt_fraction = segmentation_allele_count == 0
            ? std::numeric_limits<double>::quiet_NaN()
            : segmentation_allele_sum / static_cast<double>(segmentation_allele_count);
        point.copy_ratio_samples.clear();
        point.copy_ratio_samples.reserve(tables.size());
        point.minor_allele_fraction_samples.clear();
        point.minor_allele_fraction_samples.reserve(tables.size());
        point.segmentation_alt_fraction_samples.clear();
        point.segmentation_alt_fraction_samples.reserve(tables.size());
        for (const auto& index : indexes) {
            const auto iterator = index.find(key);
            if (iterator == index.end()) continue;
            point.copy_ratio_samples.push_back(iterator->second->value);
            point.minor_allele_fraction_samples.push_back(iterator->second->minor_allele_fraction);
            point.segmentation_alt_fraction_samples.push_back(
                iterator->second->segmentation_alt_fraction);
        }
        joint.points.push_back(std::move(point));
    }
    if (joint.points.empty())
        throw std::runtime_error("BAD_INPUT: multi-sample ModelSegments has no common interval rows");
    const bool supplied_segmentation = !options.input_segments.empty();
    const bool kernel_segmentation = options.kernel_requested && !supplied_segmentation;
    const auto segments = supplied_segmentation
        ? summarize_input_segments(joint, read_input_segments(options.input_segments, joint))
        : (kernel_segmentation
            ? kernel_segment_points_multisample(joint, options, tables.size(), kernel_stats)
            : segment_points(joint, options.change_point_threshold, options.min_points,
                             options.allelic_change_point_threshold, options.threads, kernel_stats));
    const auto output_path = multisample_interval_output(options);
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write multi-sample interval list: " + output_path);
    output << "@HD\tVN:1.6\n";
    for (std::size_t index = 0; index < joint.contigs.size(); ++index)
        output << "@SQ\tSN:" << joint.contigs[index] << "\tLN:" << joint.contig_lengths[index] << '\n';
    output << "@CO\tfastgatk multisample bounded joint segmentation\n";
    for (const auto& segment : segments)
        output << segment.contig << '\t' << segment.start << '\t' << segment.end
               << "\t+\tfastgatk_model_segments\n";
    output.close();
    if (!output || !file_complete(output_path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete multi-sample interval list: " + output_path);
    const auto manifest_path = options.manifest.empty() ? output_path + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write multi-sample manifest: " + manifest_path);
    const auto wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    manifest << "{\"schema_version\":1,\"tool\":\"ModelSegments\",\"implementation\":\"fastgatk-model-segments\","
             << "\"status\":\"prototype\",\"determinism\":\"strict\",\"multisample\":true,\"sample_count\":" << tables.size()
             << ",\"common_points\":" << joint.points.size() << ",\"segments\":" << segments.size()
             << ",\"segmentation_method\":\"" << (kernel_segmentation ? "KernelSegmenter" : "threshold") << "\""
             << ",\"segmentation_allele_fraction_semantics\":\"first-oriented-alt-per-copy-ratio-interval\""
             << ",\"execution_policy\":\"RangePolicy\",\"execution_space\":\""
             << Kokkos::DefaultExecutionSpace::name() << "\",\"bit_identical_to_gatk\":false,\"inputs\":[";
    for (std::size_t index = 0; index < options.denoised_copy_ratios.size(); ++index) {
        if (index != 0) manifest << ',';
        manifest << "{\"denoised_copy_ratios\":\"" << json_escape(options.denoised_copy_ratios[index]) << "\"";
        if (!options.allelic_count_files.empty())
            manifest << ",\"allelic_counts\":\"" << json_escape(options.allelic_count_files[index]) << "\"";
        manifest << '}';
    }
    manifest << "],\"outputs\":[{\"path\":\"" << json_escape(output_path)
             << "\",\"kind\":\"interval-list\",\"complete\":true,\"bytes\":"
             << file_bytes(output_path) << "}],\"telemetry\":{\"output_bytes\":"
             << file_bytes(output_path) << ",\"wall_seconds\":" << wall_seconds
             << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
             << "\",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
             << ",\"kernel_execution_policy\":\"RangePolicy\""
             << ",\"kernel_batches\":" << (kernel_stats == nullptr ? 0 : kernel_stats->batches)
             << ",\"kernel_observations\":" << (kernel_stats == nullptr ? 0 : kernel_stats->observations)
             << ",\"kernel_prepare_seconds\":" << (kernel_stats == nullptr ? 0.0 : kernel_stats->prepare_seconds)
             << ",\"kernel_execute_seconds\":" << (kernel_stats == nullptr ? 0.0 : kernel_stats->execute_seconds)
             << "},\"outputs_manifest_path\":\"" << json_escape(manifest_path) << "\"}\n";
    manifest.close();
    if (!manifest || !file_complete(manifest_path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete multi-sample manifest: " + manifest_path);
    std::cout << "{\"tool\":\"ModelSegments\",\"status\":\"prototype\",\"determinism\":\"strict\",\"multisample\":true,\"sample_count\":"
              << tables.size() << ",\"common_points\":" << joint.points.size() << ",\"segments\":"
              << segments.size() << ",\"output\":\"" << json_escape(output_path)
              << "\",\"output_bytes\":" << file_bytes(output_path)
              << ",\"wall_seconds\":" << wall_seconds
              << ",\"kernel_batches\":" << (kernel_stats == nullptr ? 0 : kernel_stats->batches)
              << ",\"kernel_observations\":" << (kernel_stats == nullptr ? 0 : kernel_stats->observations)
              << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
              << "\"}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto options = parse(argc, argv);
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(static_cast<unsigned>(std::max(1, options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        const auto wall_start = std::chrono::steady_clock::now();
        KernelStats kernel_stats;
        if (options.denoised_copy_ratios.size() > 1 || options.allelic_count_files.size() > 1) {
            const auto status = run_multisample_model_segments(options, &kernel_stats);
            Kokkos::finalize();
            initialized = false;
            return status;
        }
        auto table = read_copy_ratios(options.input);
        std::size_t allelic_loci = 0;
        std::size_t heterozygous_allelic_loci = 0;
        std::size_t missing_allelic_loci = 0;
        bool allelic_signal_used = false;
        AllelicCountMap case_heterozygous_counts;
        AllelicCountMap normal_heterozygous_counts;
        std::string normal_sample;
        const AllelicCountMap* matched_normal_hets = nullptr;
        if (!options.normal_allelic_counts.empty() && options.allelic_counts.empty())
            throw std::invalid_argument("--normal-allelic-counts requires --allelic-counts in single-sample mode");
        if (!options.normal_allelic_counts.empty()) {
            const auto normal_counts = read_allelic_counts(options.normal_allelic_counts, &normal_sample);
            normal_heterozygous_counts = filter_heterozygous_counts(
                normal_counts, &table, options.min_total_allele_count_normal,
                options.genotyping_homozygous_log_ratio_threshold,
                options.genotyping_base_error_rate, &kernel_stats);
            matched_normal_hets = &normal_heterozygous_counts;
        }
        if (!options.allelic_counts.empty()) {
            const auto counts = read_allelic_counts(options.allelic_counts);
            allelic_loci = counts.size();
            case_heterozygous_counts = prepare_case_heterozygous_counts(
                table, counts, matched_normal_hets, options, &kernel_stats);
            heterozygous_allelic_loci = case_heterozygous_counts.size();
            apply_allelic_counts_to_table(table, case_heterozygous_counts);
            for (const auto& point : table.points)
                if (std::isfinite(point.minor_allele_fraction)) allelic_signal_used = true;
            for (const auto& entry : counts) {
                const auto separator = entry.first.find(':');
                const auto contig = separator == std::string::npos ? std::string{} : entry.first.substr(0, separator);
                const auto position = separator == std::string::npos ? 0 : std::stoll(entry.first.substr(separator + 1));
                (void)contig;
                (void)position;
                if (case_heterozygous_counts.find(entry.first) == case_heterozygous_counts.end())
                    ++missing_allelic_loci;
            }
            write_heterozygous_sidecars(options, table, case_heterozygous_counts,
                                        matched_normal_hets, normal_sample);
        }
        std::size_t normal_loci = 0;
        if (!options.normal_allelic_counts.empty())
            normal_loci = read_allelic_counts(options.normal_allelic_counts).size();
        const bool supplied_segmentation = !options.input_segments.empty();
        const bool kernel_segmentation = options.kernel_requested && !supplied_segmentation;
        auto segments = supplied_segmentation
            ? summarize_input_segments(table, read_input_segments(options.input_segments, table))
            : (kernel_segmentation
                ? kernel_segment_points(table, options, &kernel_stats)
                : segment_points(table, options.change_point_threshold,
                                 options.min_points,
                                 options.allelic_change_point_threshold,
                                 options.threads, &kernel_stats));
        CopyRatioConditionalModel conditional_model;
        AlleleFractionConditionalModel allele_conditional_model;
        auto posterior = sample_model_segments_posterior(
            table, segments, options, &kernel_stats, &conditional_model,
            &allele_conditional_model);
        const auto segments_before_smoothing = segments.size();
        const auto smoothing = smooth_model_segments(
            segments, posterior, options, &kernel_stats);
        // GATK refits after a no-refit smoothing pass before writing the final
        // result.  Re-run the bounded native model on the merged partition so
        // the final posterior is conditioned on the same segment set.  The
        // default native contract remains unchanged because smoothing is only
        // activated when a smoothing option is explicitly supplied.
        if (smoothing.applied) {
            conditional_model = {};
            allele_conditional_model = {};
            posterior = sample_model_segments_posterior(
                table, segments, options, &kernel_stats, &conditional_model,
                &allele_conditional_model);
        }
        const std::string sample = options.sample.empty() ? table.sample : options.sample;
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write modeled segments: " + options.output);
        output << "@HD\tVN:1.6\n";
        for (std::size_t i = 0; i < table.contigs.size(); ++i)
            output << "@SQ\tSN:" << table.contigs[i] << "\tLN:" << table.contig_lengths[i] << '\n';
        output << "@RG\tID:GATKCopyNumber\tSM:" << sample << '\n';
        // GATK's modeled-segment schema carries both posterior families.  A
        // copy-ratio-only probabilistic run still emits that schema so a
        // downstream reader can distinguish an unavailable MAF posterior from
        // a silently omitted model.  Preserve the historical five-column
        // copy-ratio-only output when no probabilistic controls were requested.
        const bool emit_allelic_columns = allelic_signal_used || options.probabilistic_requested;
        output << "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO";
        if (emit_allelic_columns) {
            output << "\tNUM_POINTS_ALLELE_FRACTION"
                   << "\tLOG2_COPY_RATIO_POSTERIOR_10\tLOG2_COPY_RATIO_POSTERIOR_50\tLOG2_COPY_RATIO_POSTERIOR_90"
                   << "\tMINOR_ALLELE_FRACTION_POSTERIOR_10\tMINOR_ALLELE_FRACTION_POSTERIOR_50\tMINOR_ALLELE_FRACTION_POSTERIOR_90";
        }
        output << '\n';
        output << std::fixed << std::setprecision(6);
        for (std::size_t segment_index = 0; segment_index < segments.size(); ++segment_index) {
            const auto& segment = segments[segment_index];
            output << segment.contig << '\t' << segment.start << '\t' << segment.end << '\t'
                   << segment.points << '\t' << segment.mean;
            if (emit_allelic_columns) {
                const auto& quantiles = posterior[segment_index];
                const double copy_ratio_10 = options.probabilistic_requested
                    ? quantiles.copy_ratio_10 : segment.mean - segment.stddev;
                const double copy_ratio_50 = options.probabilistic_requested
                    ? quantiles.copy_ratio_50 : segment.mean;
                const double copy_ratio_90 = options.probabilistic_requested
                    ? quantiles.copy_ratio_90 : segment.mean + segment.stddev;
                output << '\t' << segment.allelic_points
                       << '\t' << copy_ratio_10
                       << '\t' << copy_ratio_50
                       << '\t' << copy_ratio_90;
                if (segment.allelic_points == 0 || !std::isfinite(segment.mean_minor_allele_fraction)) {
                    output << "\tNaN\tNaN\tNaN";
                } else {
                    const double allele_fraction_10 = options.probabilistic_requested
                        ? quantiles.allele_fraction_10
                        : segment.mean_minor_allele_fraction - segment.allelic_stddev;
                    const double allele_fraction_50 = options.probabilistic_requested
                        ? quantiles.allele_fraction_50
                        : segment.mean_minor_allele_fraction;
                    const double allele_fraction_90 = options.probabilistic_requested
                        ? quantiles.allele_fraction_90
                        : segment.mean_minor_allele_fraction + segment.allelic_stddev;
                    output << '\t' << allele_fraction_10
                           << '\t' << allele_fraction_50
                           << '\t' << allele_fraction_90;
                }
            }
            output << '\n';
        }
        output.close();
        if (!output || !file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete modeled segment output: " + options.output);
        const auto sidecars = write_model_segments_sidecars(
            options, table, segments, posterior, allelic_signal_used, &conditional_model,
            &allele_conditional_model);
        if (!options.manifest.empty()) {
            const auto wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"ModelSegments\","
                     << "\"implementation\":\"fastgatk-model-segments\",\"status\":\"prototype\","
                     << "\"determinism\":\"strict\","
                     << "\"input\":\"" << json_escape(options.input) << "\",\"output\":\""
                     << json_escape(options.output) << "\",\"sample\":\"" << json_escape(sample)
                     << "\",\"input_points\":" << table.points.size() << ",\"segments\":" << segments.size()
                     << ",\"change_point_threshold\":" << options.change_point_threshold
                     << ",\"allelic_change_point_threshold\":" << options.allelic_change_point_threshold
                     << ",\"allelic_loci\":" << allelic_loci
                     << ",\"heterozygous_allelic_loci\":" << heterozygous_allelic_loci
                     << ",\"normal_allelic_loci\":" << normal_loci
                     << ",\"missing_allelic_loci\":" << missing_allelic_loci
                     << ",\"matched_normal_mode\":" << (matched_normal_hets != nullptr ? "true" : "false")
                     << ",\"minimum_total_allele_count_case\":" << options.min_total_allele_count_case
                     << ",\"minimum_total_allele_count_normal\":" << options.min_total_allele_count_normal
                     << ",\"genotyping_homozygous_log_ratio_threshold\":"
                     << options.genotyping_homozygous_log_ratio_threshold
                     << ",\"genotyping_base_error_rate\":" << options.genotyping_base_error_rate
                     << ",\"genotyping_execution_policy\":\"RangePolicy\""
                     << ",\"genotyping_execution_space\":\""
                     << Kokkos::DefaultExecutionSpace::name() << "\""
                     << ",\"allelic_signal_used\":" << (allelic_signal_used ? "true" : "false")
                     << ",\"segmentation_allele_fraction_semantics\":\"first-oriented-alt-per-copy-ratio-interval\""
                     << ",\"segmentation_method\":\""
                     << (supplied_segmentation ? "provided" : (kernel_segmentation ? "KernelSegmenter" : "threshold")) << "\""
                     << ",\"segments_input\":\"" << json_escape(options.input_segments) << "\""
                     << ",\"maximum_number_of_segments_per_chromosome\":" << options.max_segments_per_contig
                     << ",\"kernel_variance_copy_ratio\":" << options.kernel_variance_copy_ratio
                     << ",\"kernel_variance_allele_fraction\":" << options.kernel_variance_allele_fraction
                     << ",\"kernel_scaling_allele_fraction\":" << options.kernel_scaling_allele_fraction
                     << ",\"kernel_approximation_dimension\":" << options.kernel_approximation_dimension
                     << ",\"window_sizes\":[";
            for (std::size_t i = 0; i < options.window_sizes.size(); ++i) {
                if (i != 0) manifest << ',';
                manifest << options.window_sizes[i];
            }
            manifest << "]"
                     << ",\"number_of_changepoints_penalty_factor\":" << options.changepoint_penalty_factor
                     << ",\"mode\":\"" << json_escape(options.mode) << "\""
                     << ",\"num_burn_in_iterations\":" << options.num_burn_in_iterations
                     << ",\"num_samples\":" << options.num_samples
                     << ",\"minor_allele_fraction_prior_alpha\":" << options.minor_allele_fraction_prior_alpha
                     << ",\"number_of_samples_copy_ratio\":" << options.num_samples_copy_ratio
                     << ",\"number_of_burn_in_samples_copy_ratio\":" << options.num_burn_in_samples_copy_ratio
                     << ",\"number_of_samples_allele_fraction\":" << options.num_samples_allele_fraction
                     << ",\"number_of_burn_in_samples_allele_fraction\":" << options.num_burn_in_samples_allele_fraction
                     // Keep requested GATK totals separate from the native
                     // shared-chain projection.  GATK's number-of-samples
                     // values include burn-in, whereas the native compact
                     // --num-samples control is the number retained after an
                     // additional burn-in.  Conflating these two conventions
                     // made manifests impossible to audit and could make a
                     // bounded run appear to be a release-equivalent chain.
                     << ",\"mcmc_controls_provenance\":\""
                     << (options.generic_mcmc_overridden ? "native-generic" :
                         (options.gatk_mcmc_overridden ? "gatk-family-max-projection" :
                          (options.probabilistic_requested ? "native-defaults" : "not-requested")))
                     << "\""
                     << ",\"copy_ratio_requested_total_samples\":" << options.num_samples_copy_ratio
                     << ",\"copy_ratio_requested_burn_in_samples\":" << options.num_burn_in_samples_copy_ratio
                     << ",\"copy_ratio_requested_retained_samples\":"
                     << (options.num_samples_copy_ratio - options.num_burn_in_samples_copy_ratio)
                     << ",\"allele_fraction_requested_total_samples\":" << options.num_samples_allele_fraction
                     << ",\"allele_fraction_requested_burn_in_samples\":" << options.num_burn_in_samples_allele_fraction
                     << ",\"allele_fraction_requested_retained_samples\":"
                     << (options.num_samples_allele_fraction - options.num_burn_in_samples_allele_fraction)
                     << ",\"native_chain_burn_in_iterations\":" << options.num_burn_in_iterations
                     << ",\"native_chain_retained_samples\":" << options.num_samples
                     << ",\"native_chain_total_iterations\":"
                     << (options.num_burn_in_iterations + options.num_samples)
                     << ",\"maximum_number_of_smoothing_iterations\":" << options.maximum_number_of_smoothing_iterations
                     << ",\"number_of_smoothing_iterations_per_fit\":" << options.number_of_smoothing_iterations_per_fit
                     << ",\"smoothing_credible_interval_threshold_copy_ratio\":"
                     << options.smoothing_credible_interval_threshold_copy_ratio
                     << ",\"smoothing_credible_interval_threshold_allele_fraction\":"
                     << options.smoothing_credible_interval_threshold_allele_fraction
                     << ",\"mcmc_parameter_surface\":\"gatk-compatible-bounded-shared-chain\""
                     << ",\"smoothing_requested\":" << (options.smoothing_requested ? "true" : "false")
                     << ",\"smoothing_applied\":" << (smoothing.applied ? "true" : "false")
                     << ",\"smoothing_iterations\":" << smoothing.iterations
                     << ",\"segments_before_smoothing\":" << segments_before_smoothing
                     << ",\"copy_ratio_conditional_model\":\"deterministic-gibbs-responsibility-v1\""
                     << ",\"copy_ratio_conditional_iterations\":" << conditional_model.iterations
                     << ",\"copy_ratio_variance_conditional\":"
                     << (std::isfinite(conditional_model.variance) ? conditional_model.variance : 0.0)
                     << ",\"copy_ratio_outlier_probability_conditional\":"
                     << (std::isfinite(conditional_model.outlier_probability)
                         ? conditional_model.outlier_probability : 0.0)
                     << ",\"allele_fraction_conditional_model\":\"deterministic-binomial-responsibility-v1\""
                     << ",\"allele_fraction_conditional_iterations\":" << allele_conditional_model.iterations
                     << ",\"allele_fraction_mean_bias_conditional\":"
                     << (std::isfinite(allele_conditional_model.mean_bias)
                         ? allele_conditional_model.mean_bias : 0.0)
                     << ",\"allele_fraction_bias_variance_conditional\":"
                     << (std::isfinite(allele_conditional_model.bias_variance)
                         ? allele_conditional_model.bias_variance : 0.0)
                     << ",\"allele_fraction_outlier_probability_conditional\":"
                     << (std::isfinite(allele_conditional_model.outlier_probability)
                         ? allele_conditional_model.outlier_probability : 0.0)
                     << ",\"allele_fraction_initial_segment_means\":[";
            for (std::size_t index = 0; index < allele_conditional_model.initial_segment_means.size(); ++index) {
                if (index != 0) manifest << ',';
                const auto value = allele_conditional_model.initial_segment_means[index];
                if (std::isfinite(value)) manifest << value;
                else manifest << "null";
            }
            manifest << "]"
                     << ",\"posterior_sample_block_size\":4096"
                     << ",\"posterior_reservoir_capacity\":2048"
                     << ",\"posterior_sampler\":\""
                     << (options.probabilistic_requested ? "deterministic-gibbs-conditional-random-walk-v2" : "none") << "\""
                     << ",\"posterior_rng\":\"splitmix64-box-muller-fixed-seed-1216\""
                     << ",\"posterior_report_semantics\":\"bounded-native-random-walk-quantiles\""
                     << ",\"parameter_report_semantics\":\"conditional-point-estimate-repeated-deciles\""
                     << ",\"parameter_report_distinct_draws\":1"
                     << ",\"model_begin_final_segment_reports_distinct\":false"
                     << ",\"model_begin_final_parameter_reports_distinct\":false"
                     << ",\"gatk_full_gibbs_slice_equivalent\":false"
                     << ",\"java_rng_bit_identity\":false"
                     << ",\"execution_policy\":\"RangePolicy\""
                     << ",\"poN\":false,\"gc_model\":false,\"mcmc\":"
                     << (options.probabilistic_requested ? "true" : "false") << ",\"execution_space\":\""
                     << Kokkos::DefaultExecutionSpace::name() << "\",\"sidecars\":[";
            for (std::size_t i = 0; i < sidecars.paths.size(); ++i) {
                if (i != 0) manifest << ',';
                manifest << "\"" << json_escape(sidecars.paths[i]) << "\"";
            }
            manifest << "],\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"modeled-segments\",\"complete\":true,\"bytes\":"
                     << file_bytes(options.output) << "}";
            std::uintmax_t sidecar_bytes = 0;
            for (const auto& path : sidecars.paths) {
                manifest << ",{\"path\":\"" << json_escape(path)
                         << "\",\"kind\":\"segment-sidecar\",\"complete\":true,\"bytes\":"
                         << file_bytes(path) << "}";
                sidecar_bytes += file_bytes(path);
            }
            manifest << "],\"telemetry\":{\"output_bytes\":" << file_bytes(options.output)
                     << ",\"sidecar_bytes\":" << sidecar_bytes
                     << ",\"wall_seconds\":" << wall_seconds
                     << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                     << "\",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                     << ",\"kernel_execution_policy\":\"RangePolicy\""
                     << ",\"kernel_batches\":" << kernel_stats.batches
                     << ",\"kernel_observations\":" << kernel_stats.observations
                     << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds
                     << "}}\n";
            manifest.close();
            if (!manifest || !file_complete(options.manifest))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest: " + options.manifest);
        }
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"ModelSegments\",\"status\":\"prototype\",\"segmentation_method\":\""
                  << (supplied_segmentation ? "provided" : (kernel_segmentation ? "KernelSegmenter" : "threshold"))
                  << "\",\"segments_input\":\"" << json_escape(options.input_segments)
                  << "\",\"input_points\":"
                  << table.points.size() << ",\"segments\":" << segments.size()
                  << ",\"smoothing_applied\":" << (smoothing.applied ? "true" : "false")
                  << ",\"smoothing_iterations\":" << smoothing.iterations
                  << ",\"mcmc\":" << (options.probabilistic_requested ? "true" : "false")
                  << ",\"mode\":\"" << json_escape(options.mode) << "\""
                  << ",\"output\":\"" << json_escape(options.output) << "\",\"output_bytes\":"
                  << file_bytes(options.output)
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
