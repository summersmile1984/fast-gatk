#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <zlib.h>

#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include "fmc_java_model.hpp"
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string reference;
    std::string stats;
    std::string filtering_stats;
    std::string manifest;
    // GATK's sites-only writer option is applied only after all filtering
    // calculations.  The input header/record retains FORMAT data until that
    // final boundary so AD/F1R2 and sample selection remain correct.
    bool sites_only_vcf_output = false;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> orientation_priors;
    float min_tlod = 0.0F;
    float max_contamination = -1.0F;
    std::string contamination_table;
    std::string tumor_segmentation;
    std::string tumor_sample;
    float min_orientation_balance = -1.0F;
    float min_allele_fraction = 0.0F;
    int min_reads_per_strand = 0;
    int unique_alt_read_count = 0;
    int max_alt_allele_count = 1;
    // Match M2FiltersArgumentCollection defaults.  A negative value is
    // reserved for the explicit "disable this hard filter" override; it is
    // not the Java default.  Keeping these defaults here is important for a
    // direct FilterMutectCalls replacement because GATK applies the
    // thresholds whenever the corresponding annotation is present.
    int min_median_base_quality = 20;
    // GATK's argument field is initialized to -1, but
    // M2FiltersArgumentCollection.getMinMedianMappingQuality() resolves that
    // sentinel to 30 for normal mode (and 20 for microbial mode) before the
    // MappingQualityFilter is built.  Store the effective normal-mode default
    // here so native matches the runtime behavior, not the raw help value.
    int min_median_mapping_quality = 30;
    int min_median_read_position = 1;
    int max_median_fragment_length_difference = 10000;
    bool min_median_mapping_quality_explicit = false;
    bool log_snv_prior_explicit = false;
    bool log_indel_prior_explicit = false;
    bool mitochondria_mode = false;
    bool microbial_mode = false;
    // GATK M2FiltersArgumentCollection defaults this to +Infinity.  Keep
    // -1 as an explicit native disable override, but do not conflate the
    // Java default with disabled parsing/manifest state.
    double max_n_ratio = std::numeric_limits<double>::infinity();
    double max_contamination_probability = 0.1;
    double normal_artifact_lod = 0.0;
    double normal_pileup_p_value_threshold = 0.001;
    double filter_error_probability_threshold = 0.1;
    int min_slippage_length = 8;
    double slippage_rate = 0.1;
    int long_indel_length = 5;
    double log_snv_prior = -6.0 * std::log(10.0);
    double log_indel_prior = -7.0 * std::log(10.0);
    double log_artifact_prior = -1.0 * std::log(10.0);
    float min_somatic_probability = 0.0F;
    float max_germline_probability = -1.0F;
    float max_artifact_probability = -1.0F;
    float max_orientation_artifact_probability = -1.0F;
    std::string orientation_threshold_strategy = "CONSTANT";
    std::string joint_threshold_strategy = "OPTIMAL_F_SCORE";
    double initial_threshold = 0.1;
    double max_false_discovery_rate = 0.05;
    double f_score_beta = 1.0;
    bool threshold_strategy_explicit = false;
    int max_events_in_region = 3;
    int max_events_in_haplotype = 2;
    int max_intra_haplotype_distance = 100;
    bool create_index = true;
};

struct SomaticModelParameters {
    static constexpr int kMaxClusters = 7;
    bool learned = false;
    double log_snv_prior = -6.0 * std::log(10.0);
    double log_indel_prior = -7.0 * std::log(10.0);
    double log_variant_vs_artifact_prior = -1.0 * std::log(10.0);
    std::array<double, 21> log_variant_priors{};
    double callable_sites = std::numeric_limits<double>::quiet_NaN();
    double background_alpha = 1.0;
    double background_beta = 1.0;
    double high_af_alpha = 10.0;
    double high_af_beta = 1.0;
    double background_weight = 0.99;
    double high_af_weight = 0.01;
    int cluster_count = 2;
    std::array<double, kMaxClusters> cluster_log_weights{};
    std::array<double, kMaxClusters> cluster_alpha{};
    std::array<double, kMaxClusters> cluster_beta{};
    std::array<double, kMaxClusters> cluster_mean{};
    std::array<int, kMaxClusters> cluster_binomial{};
};

struct SomaticModelDatum {
    int total_count = 0;
    int alt_count = 0;
    int indel_length = 0;
    double somatic_weight = 1.0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = [&]() {
        const auto primary = option_value(argument, name);
        if (!primary.empty() || short_name == nullptr) return primary;
        return option_value(argument, short_name);
    }();
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

bool looks_like_gatk_mutect_stats(const std::string& path) {
    if (path.empty()) return false;
    std::ifstream input(path);
    if (!input) return false;
    std::string line;
    for (int count = 0; count < 64 && std::getline(input, line); ++count) {
        if (line.find("#<METADATA>") != std::string::npos ||
            line.rfind("statistic\tvalue", 0) == 0 ||
            line.rfind("callable\t", 0) == 0)
            return true;
    }
    return false;
}

bool header_declares_mutect2(const bcf_hdr_t* header) {
    if (header == nullptr) return false;
    for (int index = 0; index < header->nhrec; ++index) {
        const auto* record = header->hrec[index];
        if (record != nullptr && record->type == BCF_HL_GEN &&
            record->key != nullptr && record->value != nullptr &&
            std::string_view(record->key) == "source" &&
            std::string_view(record->value) == "Mutect2")
            return true;
    }
    return false;
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

// JSON has no non-finite numeric literals. Preserve Java/GATK's +Infinity
// threshold as a stable string in telemetry rather than emitting invalid
// JSON (`inf`), while finite values remain JSON numbers.
std::string json_double(double value) {
    if (std::isinf(value)) return "\"Infinity\"";
    if (std::isnan(value)) return "null";
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

// QualityUtils.errorProbToQual(errorRate) with SAMUtils.MAX_PHRED_SCORE=93.
// GATK rounds the phred value, then bounds it to [1, 93]; keep this conversion
// on Host because it is an output encoding boundary, not posterior math.
int gatk_error_probability_to_quality(const double probability) {
    if (!std::isfinite(probability) || probability <= 0.0) return 93;
    if (probability >= 1.0) return 1;
    const auto raw = static_cast<long long>(std::lround(-10.0 * std::log10(probability)));
    return static_cast<int>(std::clamp(raw, 1LL, 93LL));
}

struct NativeFilterStatsEntry {
    double false_positive_count = 0.0;
    double false_negative_count = 0.0;
};

// FilterMutectCalls still has a few scalar posterior boundaries because the
// surrounding VCF/INFO/FORMAT representation is record-oriented.  Keep those
// calls on the same lifecycle as every bulk kernel so backend choice and the
// host/device hand-off remain observable and testable.
struct FilterKernelTelemetry {
    std::size_t batches = 0;
    std::size_t observations = 0;
    std::size_t batch_calls = 0;
    std::size_t batch_observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space = Kokkos::DefaultExecutionSpace::name();

    void record(const fastgatk::core::PlanTelemetry& plan, std::size_t count) {
        ++batches;
        observations += count;
        prepare_seconds += plan.prepare_seconds;
        execute_seconds += plan.execute_seconds;
        execution_space = Kokkos::DefaultExecutionSpace::name();
    }

    void record_batch(const fastgatk::core::PlanTelemetry& plan, std::size_t count) {
        record(plan, count);
        ++batch_calls;
        batch_observations += count;
    }
};

// Execute a one-result posterior through the common plan/lifetime boundary.
// The expression is expanded inside KOKKOS_LAMBDA, so all math remains in the
// Kokkos portability layer instead of selecting a host-only implementation.
#define FASTGATK_FILTER_SCALAR_KERNEL(label_literal, expression, telemetry_ptr) \
    do { \
        using FilterExecSpace = Kokkos::DefaultExecutionSpace; \
        fastgatk::core::HostBatch filter_host_batch(label_literal); \
        filter_host_batch.records = 1; \
        filter_host_batch.bytes = sizeof(double); \
        fastgatk::core::KernelPlan<FilterExecSpace> filter_plan(label_literal); \
        filter_plan.begin_prepare(filter_host_batch); \
        Kokkos::View<double*> filter_result_device(label_literal, 1); \
        fastgatk::core::DeviceBatch<FilterExecSpace> filter_device_batch(1); \
        filter_device_batch.bind("result", filter_result_device); \
        FilterExecSpace().fence(); \
        filter_plan.end_prepare(filter_device_batch); \
        filter_plan.begin_execute(); \
        Kokkos::parallel_for(label_literal, \
            Kokkos::RangePolicy<FilterExecSpace>(0, 1), \
            KOKKOS_LAMBDA(const int) { filter_result_device(0) = (expression); }); \
        FilterExecSpace().fence(); \
        filter_plan.end_execute(); \
        const auto filter_result_host = Kokkos::create_mirror_view_and_copy( \
            Kokkos::HostSpace{}, filter_result_device); \
        if ((telemetry_ptr) != nullptr) \
            (telemetry_ptr)->record(filter_plan.telemetry(), 1); \
        return filter_result_host(0); \
    } while (false)

// Host-side mirror of GATK FilteringOutputStats.  The probabilities are
// produced by the same per-ALT boundary used to write AS_FilterStatus; the
// accumulator intentionally stays off device because it is a small output
// accounting reduction, not a filtering kernel.
class NativeFilteringOutputStats {
public:
    NativeFilteringOutputStats() {
        static constexpr std::array<const char*, 21> names{
            "slippage", "panel_of_normals", "clustered_events", "weak_evidence",
            "germline", "multiallelic", "strand_bias", "normal_artifact",
            "base_qual", "map_qual", "fragment", "position", "contamination",
            "duplicate", "orientation", "haplotype", "strict_strand", "n_ratio",
            "low_allele_frac", "error_probability", "low_tlod"};
        for (const auto* name : names) entries_.emplace(name, NativeFilterStatsEntry{});
    }

    void record(const std::vector<double>& combined_error_probabilities,
                const std::vector<std::pair<std::string, std::vector<double>>>& probabilities_by_filter,
                const double threshold) {
        const auto count = combined_error_probabilities.size();
        std::vector<bool> filtered(count, false);
        for (std::size_t index = 0; index < count; ++index) {
            const auto probability = std::clamp(combined_error_probabilities[index], 0.0, 1.0);
            filtered[index] = probability > threshold;
            if (filtered[index]) {
                false_negative_total_ += 1.0 - probability;
            } else {
                ++pass_count_;
                false_positive_count_ += probability;
                true_positive_count_ += 1.0 - probability;
            }
        }
        for (const auto& [name, probabilities] : probabilities_by_filter) {
            auto& entry = entries_[name];
            const auto limit = std::min(count, probabilities.size());
            for (std::size_t index = 0; index < limit; ++index) {
                const auto probability = std::clamp(probabilities[index], 0.0, 1.0);
                if (probability > 1.0e-10 && probability > threshold - 1.0e-10) {
                    entry.false_negative_count += 1.0 -
                        std::clamp(combined_error_probabilities[index], 0.0, 1.0);
                } else if (!filtered[index]) {
                    entry.false_positive_count += probability;
                }
            }
        }
    }

    std::uint64_t pass_count() const { return pass_count_; }
    double true_positive_count() const { return true_positive_count_; }
    double false_positive_count() const { return false_positive_count_; }
    double false_negative_count() const { return false_negative_total_; }
    std::size_t filter_count() const { return entries_.size(); }

    std::string to_json() const {
        std::ostringstream output;
        output << std::setprecision(17)
               << "{\"pass\":" << pass_count_
               << ",\"TP\":" << true_positive_count_
               << ",\"FP\":" << false_positive_count_
               << ",\"FN\":" << false_negative_total_
               << ",\"filters\":{";
        bool first = true;
        const auto total_true = true_positive_count_ + false_negative_total_;
        for (const auto& [name, entry] : entries_) {
            if (!first) output << ',';
            first = false;
            const auto fdr = pass_count_ > 0 ? entry.false_positive_count /
                                              static_cast<double>(pass_count_) : 0.0;
            const auto fnr = total_true > 0.0 ? entry.false_negative_count / total_true : 0.0;
            output << "\"" << json_escape(name) << "\":{\"FP\":"
                   << entry.false_positive_count << ",\"FDR\":" << fdr
                   << ",\"FN\":" << entry.false_negative_count
                   << ",\"FNR\":" << fnr << '}';
        }
        output << "}}";
        return output.str();
    }

    void write_gatk_table(const std::string& path, const double threshold,
                          const std::vector<std::pair<std::string, std::string>>& clustering_metadata) const {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot write GATK filtering stats: " + path);
        output << std::setprecision(17);
        for (const auto& [key, value] : clustering_metadata)
            output << "#<METADATA>" << key << '=' << value << '\n';
        const auto total_true = true_positive_count_ + false_negative_total_;
        const auto fdr = pass_count_ > 0 ? false_positive_count_ / static_cast<double>(pass_count_) : 0.0;
        const auto sensitivity = total_true > 0.0 ? true_positive_count_ / total_true : 0.0;
        output << "#<METADATA>threshold=" << std::fixed << std::setprecision(3)
               << std::round(threshold * 1000.0) / 1000.0 << '\n'
               << "#<METADATA>fdr=" << std::round(fdr * 1000.0) / 1000.0 << '\n'
               << "#<METADATA>sensitivity=" << std::round(sensitivity * 1000.0) / 1000.0 << '\n'
               << "filter\tFP\tFDR\tFN\tFNR\n";
        output << std::setprecision(2);
        for (const auto& [name, entry] : entries_) {
            if (!(entry.false_positive_count > 0.0 || entry.false_negative_count > 0.0)) continue;
            const auto filter_fdr = pass_count_ > 0 ? entry.false_positive_count /
                                                    static_cast<double>(pass_count_) : 0.0;
            const auto filter_fnr = total_true > 0.0 ? entry.false_negative_count / total_true : 0.0;
            output << name << '\t' << entry.false_positive_count << '\t' << filter_fdr
                   << '\t' << entry.false_negative_count << '\t' << filter_fnr << '\n';
        }
        output.flush();
    }

private:
    std::map<std::string, NativeFilterStatsEntry> entries_;
    std::uint64_t pass_count_ = 0;
    double true_positive_count_ = 0.0;
    double false_positive_count_ = 0.0;
    double false_negative_total_ = 0.0;
};

std::string normalize_threshold_strategy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
    if (value == "FDR") value = "FALSE_DISCOVERY_RATE";
    if (value == "F1" || value == "F_SCORE") value = "OPTIMAL_F_SCORE";
    if (value != "CONSTANT" && value != "FALSE_DISCOVERY_RATE" && value != "OPTIMAL_F_SCORE")
        throw std::invalid_argument("--threshold-strategy must be CONSTANT, FALSE_DISCOVERY_RATE, or OPTIMAL_F_SCORE");
    return value;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-filter-mutect-calls (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE             reference FASTA (+ .fai) for orientation priors\n"
                         "  -V, --variant FILE              Mutect2 VCF/VCF.GZ\n"
                         "  -O, --output FILE               filtered VCF/VCF.GZ\n"
                         "  -L, --intervals REGION          interval selector (repeatable)\n"
                         "      --interval-set-rule RULE   UNION (default) or INTERSECTION\n"
                         "      --stats FILE                Mutect2 input stats table (existing GATK table is preserved) or native stats output\n"
                         "      --filtering-stats FILE      native filtering-stats JSON output\n"
                         "      --create-output-variant-index[=BOOL], -OVI BOOL  write VCF index (default true)\n"
                         "      --sites-only-vcf-output B   omit FORMAT/sample columns (default false)\n"
                         "      --min-tlod FLOAT            approximate TLOD threshold (default 0)\n"
                         "      --contamination-fraction F  contamination estimate (POPAF posterior; AF fallback)\n"
                         "      --contamination-estimate F  GATK alias for contamination fraction\n"
                         "      --max-contamination-probability F maximum contamination posterior\n"
                         "      --contamination-table FILE  GATK contamination table\n"
                         "      --tumor-segmentation FILE  GATK tumor segmentation table\n"
                         "      --tumor-sample NAME         sample whose AD evidence is filtered\n"
                         "                                   (with orientation priors, explicitly limits the prior filter to this sample)\n"
                         "      --min-orientation-balance F minimum minor F1R2/R1F2 fraction\n"
                         "      --min-allele-fraction F    minimum tumor allele fraction\n"
                         "      --min-reads-per-strand N   minimum alternate reads on each strand\n"
                         "      --unique-alt-read-count N minimum alternate read support\n"
                         "      --max-alt-allele-count N  maximum concrete ALT alleles\n"
                         "      --min-median-base-quality N minimum median ALT base quality (default 20)\n"
                         "      --min-median-mapping-quality N minimum median ALT mapping quality (default 30)\n"
                         "      --min-median-read-position N minimum median ALT read position (default 1)\n"
                         "      --max-median-fragment-length-difference N maximum MFRL ref/ALT difference (default 10000)\n"
                         "      --max-n-ratio F             maximum NCount/ALT-depth ratio (default +Infinity)\n"
                         "      --mitochondria-mode[=BOOL]  use GATK mitochondrial filter defaults\n"
                         "      --microbial-mode[=BOOL]     use GATK microbial filter defaults\n"
                         "      --normal-pileup-p-value-threshold F normal pileup p-value threshold\n"
                         "      --filter-error-probability-threshold F weak/slippage posterior threshold (default 0.1)\n"
                         "      --min-slippage-length N  minimum STR length for polymerase slippage (default 8)\n"
                         "      --slippage-rate F         polymerase slippage beta prior (default 0.1)\n"
                         "      --long-indel-length N     long-indel mapping-quality rescue length (default 5)\n"
                         "      --log-snv-prior F         natural-log prior for somatic SNVs (default log10(-6))\n"
                         "      --log-indel-prior F       natural-log prior for somatic indels (default log10(-7))\n"
                         "      --log-artifact-prior F   natural-log prior of variant versus artifact (default log10(-1))\n"
                         "      --min-somatic-probability F minimum PSOMATIC posterior\n"
                         "      --max-germline-probability F maximum PGERMLINE posterior\n"
                         "      --max-germline-posterior F  GATK alias for max germline probability\n"
                         "      --max-artifact-probability F maximum PARTIFACT posterior\n"
                         "      --max-strand-artifact-probability F GATK orientation artifact threshold\n"
                         "      --normal-artifact-lod F   GATK normal artifact LOD provenance\n"
                         "      --orientation-bias-artifact-priors FILE\n"
                         "                                   GATK orientation prior tar.gz (repeatable)\n"
                         "      --ob-priors FILE             alias for orientation prior tar.gz\n"
                         "      --max-orientation-artifact-probability F\n"
                         "                                   maximum orientation artifact posterior (default 0.5 when priors are supplied)\n"
                         "      --threshold-strategy NAME   CONSTANT, FALSE_DISCOVERY_RATE, or OPTIMAL_F_SCORE\n"
                         "      --initial-threshold F       initial posterior threshold for learned strategy\n"
                         "      --false-discovery-rate F    maximum learned false discovery rate\n"
                         "      --f-score-beta F            beta for learned optimal F-score\n"
                         "      --max-events-in-region N    maximum ECNT before clustered_events\n"
                         "      --max-events-in-haplotype N maximum ECNTH before clustered_events\n"
                         "      --max-intra-haplotype-distance N  phased haplotype filter radius (default 100)\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        } else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--variant") || argument == "-V")
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        else if (is_option(argument, "--input") || argument == "-I")
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                is_option(argument, "--region") ? "--region" :
                is_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (is_option(argument, "--interval-set-rule") || is_option(argument, "-isr")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (is_option(argument, "--stats"))
            options.stats = require_value(index, argc, argv, argument, "--stats");
        else if (is_option(argument, "--filtering-stats"))
            options.filtering_stats = require_value(index, argc, argv, argument, "--filtering-stats");
        else if (is_option(argument, "--sites-only-vcf-output")) {
            const auto value = require_value(index, argc, argv, argument,
                                             "--sites-only-vcf-output");
            if (value == "true" || value == "TRUE" || value == "1")
                options.sites_only_vcf_output = true;
            else if (value == "false" || value == "FALSE" || value == "0")
                options.sites_only_vcf_output = false;
            else
                throw std::invalid_argument(
                    "--sites-only-vcf-output must be true or false");
        }
        else if (is_option(argument, "--min-tlod"))
            options.min_tlod = std::stof(require_value(index, argc, argv, argument, "--min-tlod"));
        else if (is_option(argument, "--contamination-fraction") ||
                 is_option(argument, "--contamination-estimate"))
            options.max_contamination = std::stof(require_value(
                index, argc, argv, argument,
                argument.rfind("--contamination-estimate", 0) == 0
                    ? "--contamination-estimate" : "--contamination-fraction"));
        else if (is_option(argument, "--max-contamination-probability"))
            options.max_contamination_probability = std::stod(require_value(
                index, argc, argv, argument, "--max-contamination-probability"));
        else if (is_option(argument, "--contamination-table"))
            options.contamination_table = require_value(index, argc, argv, argument, "--contamination-table");
        else if (is_option(argument, "--tumor-segmentation"))
            options.tumor_segmentation = require_value(index, argc, argv, argument, "--tumor-segmentation");
        else if (is_option(argument, "--tumor-sample"))
            options.tumor_sample = require_value(index, argc, argv, argument, "--tumor-sample");
        else if (is_option(argument, "--min-orientation-balance"))
            options.min_orientation_balance = std::stof(require_value(index, argc, argv, argument, "--min-orientation-balance"));
        else if (is_option(argument, "--min-allele-fraction") || is_option(argument, "--min-af"))
            options.min_allele_fraction = std::stof(require_value(
                index, argc, argv, argument,
                argument.rfind("--min-af", 0) == 0 ? "--min-af" : "--min-allele-fraction"));
        else if (is_option(argument, "--min-reads-per-strand") ||
                 is_option(argument, "--min-reads-on-each-strand"))
            options.min_reads_per_strand = std::stoi(require_value(
                index, argc, argv, argument,
                argument.rfind("--min-reads-on-each-strand", 0) == 0
                    ? "--min-reads-on-each-strand" : "--min-reads-per-strand"));
        else if (is_option(argument, "--unique-alt-read-count") || argument == "--unique")
            options.unique_alt_read_count = std::stoi(require_value(
                index, argc, argv, argument, "--unique-alt-read-count", "--unique"));
        else if (is_option(argument, "--max-alt-allele-count"))
            options.max_alt_allele_count = std::stoi(require_value(
                index, argc, argv, argument, "--max-alt-allele-count"));
        else if (is_option(argument, "--min-median-base-quality"))
            options.min_median_base_quality = std::stoi(require_value(
                index, argc, argv, argument, "--min-median-base-quality"));
        else if (is_option(argument, "--min-median-mapping-quality")) {
            options.min_median_mapping_quality = std::stoi(require_value(
                index, argc, argv, argument, "--min-median-mapping-quality"));
            options.min_median_mapping_quality_explicit = true;
        }
        else if (is_option(argument, "--min-median-read-position"))
            options.min_median_read_position = std::stoi(require_value(
                index, argc, argv, argument, "--min-median-read-position"));
        else if (is_option(argument, "--max-median-fragment-length-difference"))
            options.max_median_fragment_length_difference = std::stoi(require_value(
                index, argc, argv, argument, "--max-median-fragment-length-difference"));
        else if (is_option(argument, "--max-n-ratio"))
            options.max_n_ratio = std::stod(require_value(index, argc, argv, argument, "--max-n-ratio"));
        else if (is_option(argument, "--normal-pileup-p-value-threshold") ||
                 is_option(argument, "--normal-p-value-threshold"))
            options.normal_pileup_p_value_threshold = std::stod(require_value(
                index, argc, argv, argument,
                argument.rfind("--normal-p-value-threshold", 0) == 0
                    ? "--normal-p-value-threshold" : "--normal-pileup-p-value-threshold"));
        else if (is_option(argument, "--filter-error-probability-threshold"))
            options.filter_error_probability_threshold = std::stod(require_value(
                index, argc, argv, argument, "--filter-error-probability-threshold"));
        else if (is_option(argument, "--min-slippage-length") ||
                 is_option(argument, "--min-pcr-slippage-size"))
            options.min_slippage_length = std::stoi(require_value(
                index, argc, argv, argument,
                argument.rfind("--min-pcr-slippage-size", 0) == 0
                    ? "--min-pcr-slippage-size" : "--min-slippage-length"));
        else if (is_option(argument, "--slippage-rate") ||
                 is_option(argument, "--pcr-slippage-rate"))
            options.slippage_rate = std::stod(require_value(
                index, argc, argv, argument,
                argument.rfind("--pcr-slippage-rate", 0) == 0
                    ? "--pcr-slippage-rate" : "--slippage-rate"));
        else if (is_option(argument, "--long-indel-length"))
            options.long_indel_length = std::stoi(require_value(
                index, argc, argv, argument, "--long-indel-length"));
        else if (is_option(argument, "--log-snv-prior")) {
            options.log_snv_prior = std::stod(require_value(index, argc, argv, argument, "--log-snv-prior"));
            options.log_snv_prior_explicit = true;
        }
        else if (is_option(argument, "--log-indel-prior")) {
            options.log_indel_prior = std::stod(require_value(index, argc, argv, argument, "--log-indel-prior"));
            options.log_indel_prior_explicit = true;
        }
        else if (is_option(argument, "--log-artifact-prior"))
            options.log_artifact_prior = std::stod(require_value(index, argc, argv, argument, "--log-artifact-prior"));
        else if (is_option(argument, "--log-somatic-prior")) {
            options.log_snv_prior = std::stod(require_value(index, argc, argv, argument, "--log-somatic-prior"));
            options.log_snv_prior_explicit = true;
        }
        else if (is_option(argument, "--min-somatic-probability"))
            options.min_somatic_probability = std::stof(require_value(
                index, argc, argv, argument, "--min-somatic-probability"));
        else if (is_option(argument, "--max-germline-probability") ||
                 is_option(argument, "--max-germline-posterior"))
            options.max_germline_probability = std::stof(require_value(
                index, argc, argv, argument,
                argument.rfind("--max-germline-posterior", 0) == 0
                    ? "--max-germline-posterior" : "--max-germline-probability"));
        else if (is_option(argument, "--max-artifact-probability"))
            options.max_artifact_probability = std::stof(require_value(
                index, argc, argv, argument, "--max-artifact-probability"));
        else if (is_option(argument, "--orientation-bias-artifact-priors") ||
                 is_option(argument, "--ob-priors"))
            options.orientation_priors.push_back(require_value(
                index, argc, argv, argument,
                argument.rfind("--ob-priors", 0) == 0
                    ? "--ob-priors" : "--orientation-bias-artifact-priors"));
        else if (is_option(argument, "--max-strand-artifact-probability") ||
                 is_option(argument, "--max-orientation-artifact-probability"))
            options.max_orientation_artifact_probability = std::stof(require_value(
                index, argc, argv, argument,
                argument.rfind("--max-strand-artifact-probability", 0) == 0
                    ? "--max-strand-artifact-probability" : "--max-orientation-artifact-probability"));
        else if (is_option(argument, "--normal-artifact-lod"))
            options.normal_artifact_lod = std::stod(require_value(
                index, argc, argv, argument, "--normal-artifact-lod"));
        else if (is_option(argument, "--threshold-strategy")) {
            const auto strategy = normalize_threshold_strategy(
                require_value(index, argc, argv, argument, "--threshold-strategy"));
            options.orientation_threshold_strategy = strategy;
            options.joint_threshold_strategy = strategy;
            options.threshold_strategy_explicit = true;
        } else if (is_option(argument, "--initial-threshold"))
            options.initial_threshold = std::stod(require_value(index, argc, argv, argument, "--initial-threshold"));
        else if (is_option(argument, "--false-discovery-rate"))
            options.max_false_discovery_rate = std::stod(require_value(index, argc, argv, argument, "--false-discovery-rate"));
        else if (is_option(argument, "--f-score-beta"))
            options.f_score_beta = std::stod(require_value(index, argc, argv, argument, "--f-score-beta"));
        else if (is_option(argument, "--max-events-in-region"))
            options.max_events_in_region = std::stoi(require_value(
                index, argc, argv, argument, "--max-events-in-region"));
        else if (is_option(argument, "--max-events-in-haplotype"))
            options.max_events_in_haplotype = std::stoi(require_value(
                index, argc, argv, argument, "--max-events-in-haplotype"));
        else if (is_option(argument, "--max-intra-haplotype-distance") ||
                 is_option(argument, "--distance-on-haplotype"))
            options.max_intra_haplotype_distance = std::stoi(require_value(
                index, argc, argv, argument,
                argument.rfind("--distance-on-haplotype", 0) == 0
                    ? "--distance-on-haplotype" : "--max-intra-haplotype-distance"));
        else if (argument == "--mitochondria-mode" ||
                 argument.rfind("--mitochondria-mode=", 0) == 0)
            options.mitochondria_mode = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--mitochondria-mode");
        else if (argument == "--microbial-mode" ||
                 argument.rfind("--microbial-mode=", 0) == 0)
            options.microbial_mode = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--microbial-mode");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0 ||
                 argument == "-OVI" || argument.rfind("-OVI=", 0) == 0) {
            const char* option_name = argument == "-OVI" ||
                    argument.rfind("-OVI=", 0) == 0 ? "-OVI" :
                    "--create-output-variant-index";
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, option_name);
        } else if (argument == "--quiet" || is_option(argument, "--java-options") ||
                   is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos &&
                (is_option(argument, "--java-options") || is_option(argument, "--verbosity")))
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.microbial_mode && !options.min_median_mapping_quality_explicit)
        options.min_median_mapping_quality = 20;
    if (options.mitochondria_mode) {
        if (!options.log_snv_prior_explicit)
            options.log_snv_prior = -2.5 * std::log(10.0);
        if (!options.log_indel_prior_explicit)
            options.log_indel_prior = -3.75 * std::log(10.0);
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.max_contamination < -1.0F || options.max_contamination > 1.0F)
        throw std::invalid_argument("--contamination-fraction must be in [0,1]");
    if (options.min_orientation_balance < -1.0F || options.min_orientation_balance > 0.5F)
        throw std::invalid_argument("--min-orientation-balance must be in [0,0.5]");
    if (options.max_orientation_artifact_probability < -1.0F ||
        options.max_orientation_artifact_probability > 1.0F ||
        !std::isfinite(options.max_orientation_artifact_probability))
        throw std::invalid_argument("--max-orientation-artifact-probability must be in [0,1] or -1");
    if (!std::isfinite(options.initial_threshold) || options.initial_threshold < 0.0 || options.initial_threshold > 1.0)
        throw std::invalid_argument("--initial-threshold must be in [0,1]");
    if (!std::isfinite(options.max_false_discovery_rate) || options.max_false_discovery_rate < 0.0 ||
        options.max_false_discovery_rate > 1.0)
        throw std::invalid_argument("--false-discovery-rate must be in [0,1]");
    if (!std::isfinite(options.max_contamination_probability) ||
        options.max_contamination_probability < 0.0 || options.max_contamination_probability > 1.0)
        throw std::invalid_argument("--max-contamination-probability must be in [0,1]");
    if (!std::isfinite(options.normal_artifact_lod))
        throw std::invalid_argument("--normal-artifact-lod must be finite");
    if (!std::isfinite(options.f_score_beta) || options.f_score_beta < 0.0)
        throw std::invalid_argument("--f-score-beta must be non-negative");
    if (!std::isfinite(options.min_allele_fraction) || options.min_allele_fraction < 0.0F ||
        options.min_allele_fraction > 1.0F)
        throw std::invalid_argument("--min-allele-fraction must be in [0,1]");
    if (options.min_reads_per_strand < 0 || options.unique_alt_read_count < 0 ||
        options.max_alt_allele_count < 0 || options.min_median_base_quality < -1 ||
        options.min_median_mapping_quality < -1 || options.min_median_read_position < -1 ||
        options.max_median_fragment_length_difference < -1 ||
        options.max_events_in_region < 0 || options.max_events_in_haplotype < 0 ||
        options.max_intra_haplotype_distance < 0 || options.long_indel_length < 0)
        throw std::invalid_argument("Mutect2 hard-filter counts must be non-negative");
    if (!std::isfinite(options.max_n_ratio) && options.max_n_ratio != std::numeric_limits<double>::infinity())
        throw std::invalid_argument("--max-n-ratio must be finite or +inf");
    if (options.max_n_ratio < -1.0)
        throw std::invalid_argument("--max-n-ratio must be non-negative or -1");
    if (!std::isfinite(options.normal_pileup_p_value_threshold) ||
        options.normal_pileup_p_value_threshold < 0.0 ||
        options.normal_pileup_p_value_threshold > 1.0)
        throw std::invalid_argument("--normal-pileup-p-value-threshold must be in [0,1]");
    if (!std::isfinite(options.filter_error_probability_threshold) ||
        options.filter_error_probability_threshold < 0.0 ||
        options.filter_error_probability_threshold > 1.0)
        throw std::invalid_argument("--filter-error-probability-threshold must be in [0,1]");
    if (options.min_slippage_length < 0 || !std::isfinite(options.slippage_rate) ||
        options.slippage_rate < 0.0 || options.slippage_rate > 1.0)
        throw std::invalid_argument("--min-slippage-length must be non-negative and --slippage-rate must be in [0,1]");
    if (!std::isfinite(options.log_snv_prior) || !std::isfinite(options.log_indel_prior) ||
        !std::isfinite(options.log_artifact_prior))
        throw std::invalid_argument("--log-snv-prior, --log-indel-prior, and --log-artifact-prior must be finite natural logs");
    const auto valid_probability = [](const float value) {
        return std::isfinite(value) && value >= -1.0F && value <= 1.0F;
    };
    if (!valid_probability(options.min_somatic_probability) || options.min_somatic_probability < 0.0F ||
        !valid_probability(options.max_germline_probability) ||
        !valid_probability(options.max_artifact_probability))
        throw std::invalid_argument("posterior filter probabilities must be in [0,1] or -1");
    return options;
}

std::optional<float> parse_numeric_token(const std::string& token) {
    try {
        std::size_t consumed = 0;
        const auto value = std::stof(token, &consumed);
        if (consumed != token.size() || !std::isfinite(value)) return std::nullopt;
        return value;
    } catch (...) {
        return std::nullopt;
    }
}

struct MutectStatsSummary {
    bool parsed = false;
    double callable_sites = std::numeric_limits<double>::quiet_NaN();
    std::unordered_map<std::string, double> values;
};

MutectStatsSummary read_mutect_stats_summary(const std::string& path) {
    MutectStatsSummary result;
    if (path.empty()) return result;
    std::ifstream input(path);
    if (!input) return result;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string statistic;
        std::string value;
        if (!(fields >> statistic >> value) || statistic == "statistic") continue;
        const auto parsed = parse_numeric_token(value);
        if (!parsed) continue;
        result.parsed = true;
        result.values[statistic] = static_cast<double>(*parsed);
        if (statistic == "callable") result.callable_sites = static_cast<double>(*parsed);
    }
    return result;
}

float read_contamination_estimate(
    const std::string& path, bool segmentation,
    std::unordered_map<std::string, float>* by_sample = nullptr) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open contamination resource: " + path);
    if (by_sample != nullptr) by_sample->clear();
    // ContaminationFilter has two distinct inputs.  A tumor-segmentation
    // resource supplies the global estimate used as the default for every
    // tumor.  A ContaminationRecord table, however, is a per-sample lookup;
    // GATK's filter uses the configured default contamination (0.0) when a
    // tumor sample is absent from that lookup.  Do not leak a value from an
    // unrelated sample into the global default.
    float estimate = segmentation ? -1.0F : 0.0F;
    bool saw_candidate = false;
    int value_column = -1;
    int sample_column = -1;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::vector<std::string> tokens;
        for (std::string token; fields >> token;) tokens.push_back(std::move(token));
        if (tokens.empty()) continue;
        bool header = false;
        for (std::size_t index = 0; index < tokens.size(); ++index) {
            std::string lower = tokens[index];
            std::transform(lower.begin(), lower.end(), lower.begin(),
                           [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
            if (lower.find("contamination") != std::string::npos) {
                value_column = static_cast<int>(index);
                header = true;
            }
            if (lower == "sample" || lower == "sample_name" || lower == "tumor_sample" ||
                lower == "tumor") sample_column = static_cast<int>(index);
        }
        if (header) continue;
        std::optional<float> candidate;
        if (value_column >= 0 && value_column < static_cast<int>(tokens.size()))
            candidate = parse_numeric_token(tokens[static_cast<std::size_t>(value_column)]);
        if (!candidate) {
            if (segmentation) {
                for (auto it = tokens.rbegin(); it != tokens.rend(); ++it) {
                    const auto parsed = parse_numeric_token(*it);
                    if (parsed && *parsed >= 0.0F && *parsed <= 1.0F) {
                        candidate = parsed;
                        break;
                    }
                }
            } else {
                for (const auto& token : tokens) {
                    const auto parsed = parse_numeric_token(token);
                    if (parsed && *parsed >= 0.0F && *parsed <= 1.0F) {
                        candidate = parsed;
                        break;
                    }
                }
            }
        }
        if (candidate && *candidate >= 0.0F && *candidate <= 1.0F) {
            saw_candidate = true;
            if (segmentation) estimate = std::max(estimate, *candidate);
        }
        if (candidate && by_sample != nullptr && sample_column >= 0 &&
            sample_column < static_cast<int>(tokens.size()) &&
            !tokens[static_cast<std::size_t>(sample_column)].empty() &&
            tokens[static_cast<std::size_t>(sample_column)] != ".") {
            (*by_sample)[tokens[static_cast<std::size_t>(sample_column)]] = *candidate;
        }
    }
    if (!saw_candidate || estimate < 0.0F)
        throw std::runtime_error("BAD_INPUT: contamination resource has no usable estimate: " + path);
    return estimate;
}

using PriorVector = std::array<double, 12>;

struct OrientationPriorCollection {
    std::string sample;
    std::unordered_map<std::string, PriorVector> contexts;
};

std::string trim_ascii(std::string value) {
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) value.pop_back();
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) ++begin;
    return value.substr(begin);
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

std::string tar_field(const unsigned char* field, std::size_t length) {
    std::size_t end = 0;
    while (end < length && field[end] != '\0') ++end;
    return std::string(reinterpret_cast<const char*>(field), end);
}

std::uint64_t parse_tar_octal(const unsigned char* field, std::size_t length) {
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < length; ++index) {
        const auto character = field[index];
        if (character >= '0' && character <= '7') value = value * 8 + (character - '0');
    }
    return value;
}

std::vector<unsigned char> read_gzip_bytes(const std::string& path) {
    gzFile stream = gzopen(path.c_str(), "rb");
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open orientation prior archive: " + path);
    std::vector<unsigned char> bytes;
    std::array<unsigned char, 64 * 1024> buffer{};
    while (true) {
        const auto read = gzread(stream, buffer.data(), static_cast<unsigned int>(buffer.size()));
        if (read < 0) {
            int error_number = 0;
            const auto* message = gzerror(stream, &error_number);
            gzclose(stream);
            throw std::runtime_error(std::string("BAD_INPUT: cannot decompress orientation prior archive: ") +
                                     (message == nullptr ? path : message));
        }
        if (read == 0) break;
        bytes.insert(bytes.end(), buffer.begin(), buffer.begin() + read);
    }
    gzclose(stream);
    return bytes;
}

double parse_finite_double(const std::string& token, const std::string& path) {
    try {
        std::size_t consumed = 0;
        const auto value = std::stod(token, &consumed);
        if (consumed != token.size() || !std::isfinite(value) || value < 0.0)
            throw std::runtime_error("invalid probability");
        return value;
    } catch (...) {
        throw std::runtime_error("BAD_INPUT: invalid orientation prior value in " + path + ": " + token);
    }
}

OrientationPriorCollection parse_orientation_prior_table(const std::string& member,
                                                           const std::string& content,
                                                           const std::string& path) {
    OrientationPriorCollection collection;
    std::vector<std::string> columns;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        line = trim_ascii(line);
        if (line.empty()) continue;
        if (line.rfind("#<METADATA>SAMPLE=", 0) == 0) {
            collection.sample = line.substr(std::string("#<METADATA>SAMPLE=").size());
            continue;
        }
        if (line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (columns.empty()) {
            columns = fields;
            continue;
        }
        if (fields.size() < columns.size())
            throw std::runtime_error("BAD_INPUT: malformed orientation prior row in " + path);
        std::unordered_map<std::string, std::size_t> position;
        for (std::size_t index = 0; index < columns.size(); ++index) position.emplace(columns[index], index);
        const auto context_it = position.find("context");
        if (context_it == position.end() || context_it->second >= fields.size())
            throw std::runtime_error("BAD_INPUT: orientation prior table has no context column: " + path);
        const auto context = fields[context_it->second];
        if (context.size() != 3) throw std::runtime_error("BAD_INPUT: orientation prior context is not a 3-mer: " + context);
        static constexpr std::array<const char*, 12> state_names{
            "f1r2_a", "f1r2_c", "f1r2_g", "f1r2_t",
            "f2r1_a", "f2r1_c", "f2r1_g", "f2r1_t",
            "hom_ref", "germline_het", "somatic_het", "hom_var"};
        PriorVector prior{};
        double sum = 0.0;
        for (std::size_t state = 0; state < state_names.size(); ++state) {
            const auto state_it = position.find(state_names[state]);
            if (state_it == position.end() || state_it->second >= fields.size())
                throw std::runtime_error("BAD_INPUT: orientation prior table is missing " +
                                         std::string(state_names[state]) + ": " + path);
            prior[state] = parse_finite_double(fields[state_it->second], path);
            sum += prior[state];
        }
        if (!(sum > 0.0)) throw std::runtime_error("BAD_INPUT: orientation prior row sums to zero: " + context);
        for (auto& value : prior) value /= sum;
        collection.contexts[context] = prior;
    }
    if (collection.sample.empty()) {
        const auto suffix = std::string(".orientation_priors");
        if (member.size() > suffix.size() && member.compare(member.size() - suffix.size(), suffix.size(), suffix) == 0)
            collection.sample = member.substr(0, member.size() - suffix.size());
    }
    if (collection.sample.empty() || collection.contexts.empty())
        throw std::runtime_error("BAD_INPUT: orientation prior archive member is empty or has no sample/context rows: " + path);
    return collection;
}

std::vector<OrientationPriorCollection> read_orientation_prior_archives(
    const std::vector<std::string>& paths) {
    std::vector<OrientationPriorCollection> collections;
    for (const auto& path : paths) {
        const auto bytes = read_gzip_bytes(path);
        std::size_t offset = 0;
        bool saw_member = false;
        while (offset + 512 <= bytes.size()) {
            const auto* header = bytes.data() + offset;
            bool all_zero = true;
            for (std::size_t index = 0; index < 512; ++index) {
                if (header[index] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero) break;
            const auto member = tar_field(header, 100);
            const auto size = parse_tar_octal(header + 124, 12);
            offset += 512;
            if (size > bytes.size() - offset)
                throw std::runtime_error("BAD_INPUT: truncated orientation prior tar member: " + path);
            static constexpr std::string_view prior_suffix = ".orientation_priors";
            if (member.size() >= prior_suffix.size() &&
                member.compare(member.size() - prior_suffix.size(), prior_suffix.size(), prior_suffix) == 0) {
                const std::string content(reinterpret_cast<const char*>(bytes.data() + offset),
                                          static_cast<std::size_t>(size));
                collections.push_back(parse_orientation_prior_table(member, content, path));
                saw_member = true;
            }
            offset += static_cast<std::size_t>(size);
            const auto padding = (512 - (size % 512)) % 512;
            if (padding > bytes.size() - offset)
                throw std::runtime_error("BAD_INPUT: truncated orientation prior tar padding: " + path);
            offset += static_cast<std::size_t>(padding);
        }
        if (!saw_member)
            throw std::runtime_error("BAD_INPUT: orientation prior archive has no .orientation_priors member: " + path);
    }
    return collections;
}

#if FASTGATK_HAS_HTSLIB

int filter_record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, static_cast<int>(values[0]));
    free(values);
    return end;
}

const char* filter_interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_filter_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<fastgatk::io::IndexedInterval> merged;
    merged.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (interval.end <= interval.begin) continue;
        if (!merged.empty() && merged.back().rid == interval.rid && interval.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    intervals.swap(merged);
}

void append_filter_interval_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<fastgatk::io::IndexedInterval>& intervals,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<fastgatk::io::IndexedInterval> incoming;
    fastgatk::io::append_interval_selector(selector, header, incoming, stats);
    normalize_filter_intervals(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        intervals.insert(intervals.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<fastgatk::io::IndexedInterval> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < intervals.size() && right < incoming.size()) {
        if (intervals[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < intervals[left].rid) { ++right; continue; }
        const int begin = std::max(intervals[left].begin, incoming[right].begin);
        const int end = std::min(intervals[left].end, incoming[right].end);
        if (begin < end) {
            auto overlap = intervals[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (intervals[left].end < incoming[right].end) ++left;
        else ++right;
    }
    intervals.swap(intersection);
}

bool filter_record_in_intervals(const bcf_hdr_t* header, bcf1_t* record,
                                const std::vector<fastgatk::io::IndexedInterval>& intervals) {
    // An empty vector is meaningful when selectors were supplied and their
    // INTERSECTION is disjoint; callers bypass this predicate when no selector
    // was requested, so fail closed here rather than silently unfiltering.
    if (intervals.empty()) return false;
    const int end = filter_record_end_exclusive(header, record);
    return std::any_of(intervals.begin(), intervals.end(), [&](const auto& interval) {
        return interval.rid == record->rid && end > interval.begin && record->pos < interval.end;
    });
}

float info_float(const bcf_hdr_t* header, bcf1_t* record, const char* name, float fallback) {
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, name, &values, &count);
    const auto result = length > 0 && count > 0 ? values[0] : fallback;
    free(values);
    return result;
}

int info_int(const bcf_hdr_t* header, bcf1_t* record, const char* name, int fallback) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, name, &values, &count);
    const auto result = length > 0 && count > 0 && values[0] != bcf_int32_missing ? values[0] : fallback;
    free(values);
    return result;
}

std::vector<float> info_float_values(const bcf_hdr_t* header, bcf1_t* record, const char* name) {
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, name, &values, &count);
    std::vector<float> result;
    if (length > 0 && count > 0) {
        result.reserve(static_cast<std::size_t>(count));
        for (int index = 0; index < count; ++index) {
            if (std::isfinite(values[index])) result.push_back(values[index]);
        }
    }
    free(values);
    return result;
}

std::vector<int> info_int_values(const bcf_hdr_t* header, bcf1_t* record, const char* name) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, name, &values, &count);
    std::vector<int> result;
    if (length > 0 && count > 0) {
        result.reserve(static_cast<std::size_t>(count));
        for (int index = 0; index < count; ++index) {
            const auto value = values[index];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end)
                result.push_back(static_cast<int>(value));
        }
    }
    free(values);
    return result;
}

// Unlike info_float_values(), this helper preserves the cardinality and
// position of missing INFO vector elements.  HTSlib represents missing and
// vector-end values as sentinels; expose them as NaN so Number=A annotations
// can be indexed without shifting a later ALT toward ALT1.
std::vector<double> info_float_values_indexed(const bcf_hdr_t* header,
                                              bcf1_t* record, const char* name) {
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, name, &values, &count);
    std::vector<double> result;
    if (length > 0 && count > 0) {
        result.reserve(static_cast<std::size_t>(count));
        for (int index = 0; index < count; ++index) {
            const auto value = values[index];
            if (bcf_float_is_vector_end(value)) break;
            result.push_back((bcf_float_is_missing(value) || !std::isfinite(value))
                                 ? std::numeric_limits<double>::quiet_NaN()
                                 : static_cast<double>(value));
        }
    }
    free(values);
    return result;
}

std::optional<std::string> info_string_value(const bcf_hdr_t* header, bcf1_t* record,
                                             const char* name) {
    if (header == nullptr || record == nullptr) return std::nullopt;
    char* encoded = nullptr;
    int count = 0;
    const auto length = bcf_get_info_string(header, record, name, &encoded, &count);
    std::optional<std::string> result;
    if (length > 0 && encoded != nullptr && count > 0) {
        const std::string value(encoded, static_cast<std::size_t>(length));
        if (!value.empty() && value != ".") result = value;
    }
    free(encoded);
    return result;
}

std::optional<float> info_float_value_at(const bcf_hdr_t* header, bcf1_t* record,
                                         const char* name, int index) {
    if (index < 0) return std::nullopt;
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, name, &values, &count);
    std::optional<float> result;
    if (length > 0 && count > index && std::isfinite(values[index]))
        result = values[index];
    free(values);
    return result;
}

bool info_present(const bcf_hdr_t* header, const bcf1_t* record, const char* name) {
    if (header == nullptr || record == nullptr) return false;
    const auto id = bcf_hdr_id2int(header, BCF_DT_ID, name);
    if (id < 0) return false;
    // Iterate the unpacked INFO vector instead of relying on the declared VCF
    // type.  PON is a Flag in GATK callsets, but callers in the wild also emit
    // a scalar/string marker; the filter is presence-based in both cases.
    for (int index = 0; index < record->n_info; ++index)
        if (record->d.info[index].key == id) return true;
    return false;
}

bool format_present(const bcf_hdr_t* header, const bcf1_t* record, const char* name) {
    if (header == nullptr || record == nullptr) return false;
    const auto id = bcf_hdr_id2int(header, BCF_DT_ID, name);
    if (id < 0) return false;
    for (int index = 0; index < record->n_fmt; ++index)
        if (record->d.fmt[index].id == id) return true;
    return false;
}

std::optional<int> info_alt_int(const bcf_hdr_t* header, bcf1_t* record,
                                const char* name, int alt_index) {
    if (alt_index < 0) return std::nullopt;
    const auto values = info_int_values(header, record, name);
    if (alt_index < static_cast<int>(values.size()))
        return values[static_cast<std::size_t>(alt_index)];
    // GATK declares AS_UNIQ_ALT_READ_COUNT as Number=A/Integer, but the
    // annotation writer may serialize allele-specific values with the raw
    // `|` delimiter.  HTSlib cannot decode that representation as int32;
    // accept it explicitly at this VCF boundary.
    char* encoded = nullptr;
    int count = 0;
    const auto length = bcf_get_info_string(header, record, name, &encoded, &count);
    if (length <= 0 || encoded == nullptr || count <= 0) {
        free(encoded);
        return std::nullopt;
    }
    std::string value(encoded, static_cast<std::size_t>(length));
    free(encoded);
    std::size_t begin = 0;
    for (int index = 0; index <= alt_index; ++index) {
        const auto end = value.find('|', begin);
        if (index == alt_index) {
            const auto token = value.substr(begin, end == std::string::npos
                                                       ? std::string::npos : end - begin);
            try {
                std::size_t consumed = 0;
                const auto parsed = std::stoi(token, &consumed);
                if (consumed == token.size() && parsed >= 0) return parsed;
            } catch (...) {
            }
            return std::nullopt;
        }
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return std::nullopt;
}

std::optional<std::pair<int, int>> info_alt_strand_counts(
    const bcf_hdr_t* header, bcf1_t* record, const char* name, int alt_index) {
    if (alt_index < 0) return std::nullopt;
    char* encoded = nullptr;
    int count = 0;
    const auto length = bcf_get_info_string(header, record, name, &encoded, &count);
    if (length <= 0 || encoded == nullptr || count <= 0) {
        free(encoded);
        return std::nullopt;
    }
    std::string value(encoded, static_cast<std::size_t>(length));
    free(encoded);
    std::vector<std::string> groups;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find('|', begin);
        groups.push_back(value.substr(begin, end == std::string::npos
                                                 ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    // The canonical representation has one reference pair followed by one
    // alternate pair per ALT.  A few producers collapse a single ALT into a
    // single four-value group; accept that form too.
    const auto group_index = groups.size() > 1
        ? std::min<std::size_t>(groups.size() - 1, static_cast<std::size_t>(alt_index + 1))
        : 0;
    if (group_index < groups.size()) {
            const auto& token = groups[group_index];
            std::vector<int> fields;
            std::size_t field_begin = 0;
            while (field_begin <= token.size()) {
                const auto field_end = token.find(',', field_begin);
                const auto field = token.substr(field_begin,
                                                field_end == std::string::npos
                                                    ? std::string::npos : field_end - field_begin);
                try {
                    std::size_t consumed = 0;
                    const auto parsed = std::stoi(field, &consumed);
                    if (consumed != field.size() || parsed < 0) return std::nullopt;
                    fields.push_back(parsed);
                } catch (...) {
                    return std::nullopt;
                }
                if (field_end == std::string::npos) break;
                field_begin = field_end + 1;
            }
            // GATK emits ref-forward,ref-reverse|alt-forward,alt-reverse
            // for the single-ALT case; multi-ALT rows may contain a four-value
            // token.  Normalize both forms to the alternate pair.
            if (fields.size() == 2) return std::make_pair(fields[0], fields[1]);
            if (fields.size() >= 4) return std::make_pair(fields[fields.size() - 2], fields.back());
            return std::nullopt;
    }
    return std::nullopt;
}

// Decode the canonical AS_SB_TABLE representation once when a filter needs
// both the selected ALT counts and the all-allele forward/reverse totals.
// StrandArtifactFilter sums reference plus every ALT before evaluating the
// beta-binomial model; using only the selected ALT depth changes both the
// null likelihood and its learned responsibilities.
std::optional<std::vector<std::pair<int, int>>> info_allele_strand_counts(
    const bcf_hdr_t* header, bcf1_t* record, const char* name) {
    char* encoded = nullptr;
    int count = 0;
    const auto length = bcf_get_info_string(header, record, name, &encoded, &count);
    if (length <= 0 || encoded == nullptr || count <= 0) {
        free(encoded);
        return std::nullopt;
    }
    std::string value(encoded, static_cast<std::size_t>(length));
    free(encoded);
    std::vector<std::pair<int, int>> result;
    std::size_t group_begin = 0;
    while (group_begin <= value.size()) {
        const auto group_end = value.find('|', group_begin);
        const auto group = value.substr(group_begin, group_end == std::string::npos
            ? std::string::npos : group_end - group_begin);
        std::vector<int> fields;
        std::size_t field_begin = 0;
        while (field_begin <= group.size()) {
            const auto field_end = group.find(',', field_begin);
            const auto field = group.substr(field_begin, field_end == std::string::npos
                ? std::string::npos : field_end - field_begin);
            try {
                std::size_t consumed = 0;
                const int parsed = std::stoi(field, &consumed);
                if (consumed != field.size() || parsed < 0) return std::nullopt;
                fields.push_back(parsed);
            } catch (...) {
                return std::nullopt;
            }
            if (field_end == std::string::npos) break;
            field_begin = field_end + 1;
        }
        if (fields.size() != 2) return std::nullopt;
        result.emplace_back(fields[0], fields[1]);
        if (group_end == std::string::npos) break;
        group_begin = group_end + 1;
    }
    return result.empty() ? std::nullopt
                          : std::optional<std::vector<std::pair<int, int>>>{std::move(result)};
}

std::optional<int> info_int_max(const bcf_hdr_t* header, bcf1_t* record, const char* name) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, name, &values, &count);
    std::optional<int> maximum;
    if (length > 0 && count > 0) {
        for (int index = 0; index < count; ++index) {
            const auto value = values[index];
            if (value == bcf_int32_missing || value == bcf_int32_vector_end) continue;
            maximum = maximum ? std::max(*maximum, static_cast<int>(value)) : static_cast<int>(value);
        }
    }
    free(values);
    return maximum;
}

float info_numeric(const bcf_hdr_t* header, bcf1_t* record, const char* name, float fallback) {
    const auto floating = info_float(header, record, name, std::numeric_limits<float>::quiet_NaN());
    if (std::isfinite(floating)) return floating;
    const auto integer = info_int(header, record, name, std::numeric_limits<int>::min());
    return integer == std::numeric_limits<int>::min() ? fallback : static_cast<float>(integer);
}

double info_probability_at(const bcf_hdr_t* header, bcf1_t* record,
                           const char* name, const int alt_index) {
    const auto values = info_float_values_indexed(header, record, name);
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (values.size() == 1) return values.front();
    if (alt_index <= 0 || alt_index > static_cast<int>(values.size()))
        return std::numeric_limits<double>::quiet_NaN();
    return values[static_cast<std::size_t>(alt_index - 1)];
}

// Return one probability per concrete ALT, preserving the ALT position even
// when an annotation is missing.  Number=1 annotations are broadcast to all
// ALTs (the same compatibility rule used by info_probability_at), while
// Number=A annotations stay allele-specific.  Keeping NaN placeholders is
// important: filtering a multi-ALT site must not silently shift ALT #2 into
// ALT #1 when one vector element is missing.
std::vector<double> info_probability_values_by_alt(
    const bcf_hdr_t* header, bcf1_t* record, const char* name) {
    std::vector<double> result;
    if (header == nullptr || record == nullptr) return result;
    for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
        const auto* allele = record->d.allele[alt_index];
        if (allele == nullptr || allele[0] == '<') continue;
        result.push_back(info_probability_at(header, record, name, alt_index));
    }
    return result;
}

double aggregate_probability_values(const std::vector<double>& values,
                                    const bool take_max,
                                    const bool require_complete,
                                    const double fallback) {
    if (values.empty()) return fallback;
    double aggregate = take_max ? -std::numeric_limits<double>::infinity()
                                : std::numeric_limits<double>::infinity();
    bool observed = false;
    for (const auto value : values) {
        if (!std::isfinite(value)) {
            if (require_complete) return fallback;
            continue;
        }
        aggregate = take_max ? std::max(aggregate, value) : std::min(aggregate, value);
        observed = true;
    }
    return observed ? aggregate : fallback;
}

int concrete_alt_count(const bcf1_t* record) {
    int count = 0;
    for (int allele = 1; allele < record->n_allele; ++allele) {
        const auto* value = record->d.allele[allele];
        if (value != nullptr && value[0] != '<') ++count;
    }
    return count;
}

std::optional<int> info_allele_int(const bcf_hdr_t* header, bcf1_t* record,
                                   const char* name, int allele_index,
                                   bool reference_first) {
    if (allele_index <= 0) return std::nullopt;
    const auto values = info_int_values(header, record, name);
    if (values.empty()) return std::nullopt;
    // MBQ/MMQ are Number=R (reference first), while MPOS is Number=A.  A
    // few compatibility VCFs use Number=1 for a single ALT; accepting the
    // compact representation keeps the scalar contract intact.
    int offset = allele_index - 1;
    if (reference_first && static_cast<int>(values.size()) >= record->n_allele)
        offset = allele_index;
    if (offset < 0 || offset >= static_cast<int>(values.size())) return std::nullopt;
    return values[static_cast<std::size_t>(offset)];
}

std::optional<int> info_first_int_optional(const bcf_hdr_t* header, bcf1_t* record,
                                           const char* name) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, name, &values, &count);
    std::optional<int> result;
    if (length > 0 && count > 0 && values[0] != bcf_int32_missing &&
        values[0] != bcf_int32_vector_end) result = static_cast<int>(values[0]);
    free(values);
    return result;
}

int tumor_alt_depth(const bcf_hdr_t* header, bcf1_t* record, int sample_index) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
    if (length <= 0 || count <= 0 || record->n_sample <= 0 || count % record->n_sample != 0 ||
        sample_index < 0 || sample_index >= record->n_sample) {
        free(values);
        return -1;
    }
    const int width = count / record->n_sample;
    int depth = 0;
    for (int index = 1; index < width && index < record->n_allele; ++index) {
        const auto value = values[sample_index * width + index];
        if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0 &&
            record->d.allele[index] != nullptr && record->d.allele[index][0] != '<')
            depth += value;
    }
    free(values);
    return depth;
}

int format_int_value(const bcf_hdr_t* header, bcf1_t* record, const char* name,
                     int sample_index, int allele_index, int fallback = -1) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, name, &values, &count);
    if (length <= 0 || count <= 0 || record->n_sample <= 0 || count % record->n_sample != 0 ||
        sample_index < 0 || sample_index >= record->n_sample) {
        free(values);
        return fallback;
    }
    const int width = count / record->n_sample;
    if (allele_index < 0 || allele_index >= width) {
        free(values);
        return fallback;
    }
    const auto value = values[sample_index * width + allele_index];
    free(values);
    return value == bcf_int32_missing || value == bcf_int32_vector_end ? fallback : static_cast<int>(value);
}

// Return one FORMAT integer vector for a sample without imposing Number=A/R
// semantics.  Orientation annotations have been emitted as Number=R by
// Mutect2 since the read-orientation model was introduced, but accepting a
// legacy Number=A/scalar fixture is useful for old native callsets.  The
// caller decides how to map the vector to an allele.
std::vector<int> format_int_values(const bcf_hdr_t* header, bcf1_t* record,
                                   const char* name, int sample_index) {
    std::vector<int> result;
    if (header == nullptr || record == nullptr || sample_index < 0 ||
        sample_index >= bcf_hdr_nsamples(header)) return result;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, name, &values, &count);
    if (length <= 0 || count <= 0 || record->n_sample <= 0 ||
        count % record->n_sample != 0) {
        free(values);
        return result;
    }
    const int width = count / record->n_sample;
    result.reserve(static_cast<std::size_t>(width));
    for (int index = 0; index < width; ++index) {
        const auto value = values[sample_index * width + index];
        if (value == bcf_int32_missing || value == bcf_int32_vector_end) {
            result.push_back(-1);
        } else {
            result.push_back(static_cast<int>(value));
        }
    }
    free(values);
    return result;
}

std::optional<float> format_float_value_at(const bcf_hdr_t* header, bcf1_t* record,
                                           const char* name, int sample_index,
                                           int allele_index) {
    if (header == nullptr || record == nullptr || sample_index < 0 || allele_index < 0 ||
        sample_index >= bcf_hdr_nsamples(header)) return std::nullopt;
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_float(header, record, name, &values, &count);
    std::optional<float> result;
    const auto samples = bcf_hdr_nsamples(header);
    if (length > 0 && count > 0 && samples > 0 && count % samples == 0) {
        const int width = count / samples;
        if (allele_index < width) {
            const auto value = values[sample_index * width + allele_index];
            if (!bcf_float_is_missing(value) && !bcf_float_is_vector_end(value) &&
                std::isfinite(value)) result = value;
        }
    }
    free(values);
    return result;
}

int format_sum_alt_depth(const bcf_hdr_t* header, bcf1_t* record, int sample_index) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
    if (length <= 0 || count <= 0 || record->n_sample <= 0 || count % record->n_sample != 0 ||
        sample_index < 0 || sample_index >= record->n_sample) {
        free(values);
        return -1;
    }
    const int width = count / record->n_sample;
    int depth = 0;
    for (int index = 1; index < width && index < record->n_allele; ++index) {
        const auto value = values[sample_index * width + index];
        if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0 &&
            record->d.allele[index] != nullptr && record->d.allele[index][0] != '<') depth += value;
    }
    free(values);
    return depth;
}

// GATK's ContaminationFilter uses Genotype.getAD() as totalAD, i.e. the sum
// of reference and alternate observations for this tumour genotype.  Keep a
// separate helper from format_sum_alt_depth(), whose historical native
// callers intentionally exclude the reference allele.
int format_total_depth(const bcf_hdr_t* header, bcf1_t* record, int sample_index) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
    if (length <= 0 || count <= 0 || record->n_sample <= 0 || count % record->n_sample != 0 ||
        sample_index < 0 || sample_index >= record->n_sample) {
        free(values);
        return -1;
    }
    const int width = count / record->n_sample;
    int depth = 0;
    for (int index = 0; index < width; ++index) {
        const auto value = values[sample_index * width + index];
        if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0)
            depth += value;
    }
    free(values);
    return depth;
}

int total_alt_depth_all_samples(const bcf_hdr_t* header, bcf1_t* record) {
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_int32(header, record, "AD", &values, &count);
    const auto samples = bcf_hdr_nsamples(header);
    if (length <= 0 || count <= 0 || samples <= 0 || count % samples != 0) {
        free(values);
        return -1;
    }
    const int width = count / samples;
    int depth = 0;
    for (int sample = 0; sample < samples; ++sample) {
        for (int index = 1; index < width && index < record->n_allele; ++index) {
            const auto value = values[sample * width + index];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0 &&
                record->d.allele[index] != nullptr && record->d.allele[index][0] != '<')
                depth += value;
        }
    }
    free(values);
    return depth;
}

// Rust/GATK's filtering engine sums AD over all samples classified as tumour.
// Keep the full allele vector here because both TumorEvidenceFilter and
// PolymeraseSlippageFilter distinguish the reference depth from each ALT.
std::vector<int> tumor_allele_depths(const bcf_hdr_t* header, bcf1_t* record,
                                     const std::vector<int>& sample_indices) {
    if (header == nullptr || record == nullptr || record->n_allele <= 0 ||
        sample_indices.empty()) return {};
    std::vector<int> result(static_cast<std::size_t>(record->n_allele), 0);
    bool observed = false;
    for (const auto sample_index : sample_indices) {
        for (int allele_index = 0; allele_index < record->n_allele; ++allele_index) {
            const int depth = format_int_value(header, record, "AD", sample_index, allele_index);
            if (depth < 0) return {};
            result[static_cast<std::size_t>(allele_index)] += depth;
            observed = true;
        }
    }
    return observed ? result : std::vector<int>{};
}

double tumor_allele_fraction_for_alt(const bcf_hdr_t* header, bcf1_t* record,
                                     int sample_index, int allele_index) {
    if (header == nullptr || record == nullptr || allele_index <= 0) return 1.0;
    const bool has_format_af = format_present(header, record, "AF");
    if (has_format_af) {
        // Rust/GATK treats a missing genotype AF as one for this hard filter;
        // this is deliberately different from a numeric zero.
        const auto value = format_float_value_at(header, record, "AF", sample_index,
                                                 allele_index - 1);
        return value ? std::clamp(static_cast<double>(*value), 0.0, 1.0) : 1.0;
    }
    const auto info_af = info_float_value_at(header, record, "AF", allele_index - 1);
    if (info_af) return std::clamp(static_cast<double>(*info_af), 0.0, 1.0);
    const int reference_depth = format_int_value(header, record, "AD", sample_index, 0);
    const int alternate_depth = format_int_value(header, record, "AD", sample_index, allele_index);
    const int total_depth = format_sum_alt_depth(header, record, sample_index);
    const int denominator = reference_depth >= 0 && total_depth >= 0
        ? reference_depth + total_depth : 0;
    if (alternate_depth >= 0 && denominator > 0)
        return std::clamp(static_cast<double>(alternate_depth) / denominator, 0.0, 1.0);
    return 1.0;
}

int format_alt_depth(const bcf_hdr_t* header, bcf1_t* record, int sample_index,
                     int allele_index) {
    return format_int_value(header, record, "AD", sample_index, allele_index, -1);
}

// FORMAT strings are allocated as one contiguous block by HTSlib.  Keep the
// ownership rule in one small helper so PGT/PID parsing cannot leak per-record
// buffers while FilterMutectCalls performs its learning pre-pass.
std::optional<std::string> format_string_value(const bcf_hdr_t* header, bcf1_t* record,
                                                const char* name, int sample_index) {
    if (header == nullptr || record == nullptr || sample_index < 0 ||
        sample_index >= bcf_hdr_nsamples(header)) return std::nullopt;
    char** values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_string(header, record, name, &values, &count);
    std::optional<std::string> result;
    if (length > sample_index && values != nullptr && values[sample_index] != nullptr) {
        const std::string value(values[sample_index]);
        if (!value.empty() && value != ".") result = value;
    }
    if (values != nullptr) {
        if (count > 0 && values[0] != nullptr) free(values[0]);
        free(values);
    }
    return result;
}

std::vector<double> format_float_values(const bcf_hdr_t* header, bcf1_t* record,
                                        const char* name, int sample_index) {
    std::vector<double> result;
    if (header == nullptr || record == nullptr || sample_index < 0 ||
        sample_index >= bcf_hdr_nsamples(header)) return result;
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_format_float(header, record, name, &values, &count);
    const auto samples = bcf_hdr_nsamples(header);
    if (length > 0 && count > 0 && samples > 0 && count % samples == 0) {
        const int width = count / samples;
        result.reserve(static_cast<std::size_t>(width));
        for (int index = 0; index < width; ++index) {
            const auto value = values[sample_index * width + index];
            if (bcf_float_is_missing(value) || bcf_float_is_vector_end(value) || !std::isfinite(value))
                continue;
            result.push_back(static_cast<double>(value));
        }
    }
    free(values);
    return result;
}

std::optional<std::string> phased_haplotype_key(const bcf_hdr_t* header, bcf1_t* record,
                                                int sample_index) {
    const auto pgt = format_string_value(header, record, "PGT", sample_index);
    const auto pid = format_string_value(header, record, "PID", sample_index);
    if (!pgt || !pid) return std::nullopt;
    // GATK's makePhasingString() deliberately concatenates PGT and PID with
    // no separator.  Matching this wire-level key keeps records in one local
    // haplotype even when the PID itself contains an underscore.
    return *pgt + *pid;
}

struct PhasedHaplotypeObservation {
    std::string key;
    double allele_fraction = 0.0;
};

double genotype_allele_fraction(const bcf_hdr_t* header, bcf1_t* record, int sample_index) {
    const auto values = format_float_values(header, record, "AF", sample_index);
    double maximum = 0.0;
    for (const auto value : values)
        if (std::isfinite(value)) maximum = std::max(maximum, value);
    if (!values.empty()) return std::clamp(maximum, 0.0, 1.0);
    // Native Mutect2 historically omitted FORMAT/AF on some records.  AD is
    // a deterministic compatibility fallback; when AF exists this branch is
    // never taken and preserves GATK's genotype-annotation semantics.
    const int reference_depth = format_int_value(header, record, "AD", sample_index, 0);
    const int alternate_depth = format_sum_alt_depth(header, record, sample_index);
    const int depth = reference_depth >= 0 && alternate_depth >= 0
        ? reference_depth + alternate_depth : 0;
    return depth > 0 ? static_cast<double>(alternate_depth) / static_cast<double>(depth) : 0.0;
}

std::vector<PhasedHaplotypeObservation> phased_haplotype_observations(
    const bcf_hdr_t* header, bcf1_t* record, const std::vector<int>& sample_indices) {
    std::vector<PhasedHaplotypeObservation> result;
    result.reserve(sample_indices.size());
    for (const auto sample_index : sample_indices) {
        const auto key = phased_haplotype_key(header, record, sample_index);
        if (key) result.push_back(PhasedHaplotypeObservation{
            *key, genotype_allele_fraction(header, record, sample_index)});
    }
    return result;
}

// A compact native representation of GATK's FilteredHaplotypeFilter state.
// The Java filter receives one probability per named filter from the Mutect2
// filtering engine.  The VCF contract exposes the already-combined PARTIFACT
// posterior, so this adapter uses that value as the artifact probability and
// preserves the exact two-pass state machine, inclusive distance test, and
// max-over-loci reduction from the reference implementation.
class LearnedHaplotypeFilter {
public:
    explicit LearnedHaplotypeFilter(const int max_distance)
        : max_distance_(max_distance) {}

    void accumulate(const std::string& key, const std::int64_t locus,
                    const double probability) {
        auto& values = accumulating_[key];
        values.emplace_back(locus, probability);
        ++accumulated_records_;
    }

    void learn_parameters_and_clear() {
        learned_ = accumulating_;
        accumulating_.clear();
    }

    double error_probability(const std::string& key, const std::int64_t locus) const {
        const auto iterator = learned_.find(key);
        if (iterator == learned_.end()) return 0.0;
        double maximum = 0.0;
        bool observed = false;
        for (const auto& [observed_locus, probability] : iterator->second) {
            const auto distance = observed_locus >= locus
                ? observed_locus - locus : locus - observed_locus;
            if (distance <= static_cast<std::int64_t>(max_distance_)) {
                maximum = observed ? std::max(maximum, probability) : probability;
                observed = true;
            }
        }
        return observed ? maximum : 0.0;
    }

    std::uint64_t accumulated_records() const { return accumulated_records_; }
    std::size_t learned_haplotypes() const { return learned_.size(); }

private:
    int max_distance_ = 100;
    std::unordered_map<std::string, std::vector<std::pair<std::int64_t, double>>> accumulating_;
    std::unordered_map<std::string, std::vector<std::pair<std::int64_t, double>>> learned_;
    std::uint64_t accumulated_records_ = 0;
};

struct OrientationCounts {
    int f1r2 = -1;
    int f2r1 = -1;
    int ref_f1r2 = -1;
    int ref_f2r1 = -1;
    int ref = -1;
    int alt = -1;
};

OrientationCounts read_orientation_counts(const bcf_hdr_t* header, bcf1_t* record,
                                          int sample_index, int alt_index) {
    OrientationCounts counts;
    // GATK stores F1R2/F2R1 as Number=R FORMAT arrays.  Older native
    // compatibility VCFs used Number=A or scalar INFO annotations; accept
    // those representations explicitly while preferring standard FORMAT.
    const auto f1r2_values = format_int_values(header, record, "F1R2", sample_index);
    const auto f2r1_values = format_int_values(header, record, "F2R1", sample_index);
    // The production Mutect2 header is Number=R (REF followed by each ALT),
    // while older native fixtures used Number=A.  Keep the source shape
    // explicit so the model's depth is built from the same vectors as GATK.
    const auto read_orientation_vector = [&](const std::vector<int>& values,
                                              int& ref_count, int& alt_count) {
        if (values.size() >= static_cast<std::size_t>(record->n_allele)) {
            ref_count = values[0];
            alt_count = values[static_cast<std::size_t>(alt_index)];
        } else if (values.size() >= static_cast<std::size_t>(record->n_allele - 1)) {
            // Number=A has no orientation-specific REF count.  Use AD's REF
            // count as a compatibility fallback for the denominator, while
            // preserving the orientation ALT observation.
            alt_count = values[static_cast<std::size_t>(alt_index - 1)];
        } else if (values.size() == 1 && alt_index == 1) {
            alt_count = values.front();
        }
    };
    read_orientation_vector(f1r2_values, counts.ref_f1r2, counts.f1r2);
    read_orientation_vector(f2r1_values, counts.ref_f2r1, counts.f2r1);
    if (counts.f1r2 < 0) counts.f1r2 = info_int(header, record, "F1R2", -1);
    if (counts.f2r1 < 0) counts.f2r1 = info_int(header, record, "R1F2", -1);
    if (counts.f1r2 < 0 || counts.f2r1 < 0) {
        const auto allele_specific = info_alt_strand_counts(header, record, "AS_SB_TABLE", alt_index);
        if (allele_specific) {
            if (counts.f1r2 < 0) counts.f1r2 = allele_specific->first;
            if (counts.f2r1 < 0) counts.f2r1 = allele_specific->second;
        }
    }
    counts.ref = format_int_value(header, record, "AD", sample_index, 0);
    counts.alt = format_int_value(header, record, "AD", sample_index, alt_index);
    if (counts.alt < 0) counts.alt = format_sum_alt_depth(header, record, sample_index);
    if (counts.ref_f1r2 < 0 && counts.ref_f2r1 < 0 && counts.ref >= 0) {
        // Legacy scalar/Number=A records do not carry REF orientation
        // counts.  Their total depth still follows GATK's AD-compatible
        // boundary, but all orientation ALT evidence remains independent.
        counts.ref_f1r2 = counts.ref;
        counts.ref_f2r1 = 0;
    } else {
        if (counts.ref_f1r2 < 0) counts.ref_f1r2 = 0;
        if (counts.ref_f2r1 < 0) counts.ref_f2r1 = 0;
    }
    return counts;
}

std::unordered_set<std::string> normal_sample_names(const bcf_hdr_t* header) {
    std::unordered_set<std::string> names;
    if (header == nullptr) return names;
    kstring_t formatted = {0, 0, nullptr};
    if (bcf_hdr_format(header, 0, &formatted) != 0 || formatted.s == nullptr || formatted.l == 0) {
        free(formatted.s);
        return names;
    }
    const std::string header_text(formatted.s, formatted.l);
    free(formatted.s);
    static constexpr std::string_view prefix = "##normal_sample=";
    std::size_t begin = 0;
    while (begin < header_text.size()) {
        const auto end = header_text.find('\n', begin);
        const auto line = header_text.substr(begin, end == std::string::npos
                                                       ? std::string::npos : end - begin);
        if (line.rfind(prefix, 0) == 0) {
            auto value = line.substr(prefix.size());
            if (!value.empty() && value.back() == '\r') value.pop_back();
            if (!value.empty()) names.insert(std::move(value));
        }
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return names;
}

bool sample_is_hom_ref(const bcf_hdr_t* header, bcf1_t* record, int sample_index) {
    if (header == nullptr || record == nullptr || record->n_sample <= 0 ||
        sample_index < 0 || sample_index >= record->n_sample) return false;
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &values, &count);
    if (length <= 0 || count <= 0 || count % record->n_sample != 0) {
        free(values);
        return false;
    }
    const int width = count / record->n_sample;
    bool saw_allele = false;
    bool hom_ref = true;
    for (int index = 0; index < width; ++index) {
        const auto value = values[sample_index * width + index];
        if (bcf_gt_is_missing(value) || value == bcf_int32_vector_end) continue;
        saw_allele = true;
        if (bcf_gt_allele(value) != 0) {
            hom_ref = false;
            break;
        }
    }
    free(values);
    return saw_allele && hom_ref;
}

double weighted_median_posterior(std::vector<std::pair<int, double>> evidence) {
    if (evidence.empty()) return 0.0;
    std::stable_sort(evidence.begin(), evidence.end(),
                     [](const auto& left, const auto& right) { return left.second < right.second; });
    std::int64_t total_depth = 0;
    for (const auto& item : evidence) total_depth += std::max(item.first, 0);
    if (total_depth <= 0) return 0.0;
    std::int64_t cumulative_depth = 0;
    for (const auto& item : evidence) {
        cumulative_depth += std::max(item.first, 0);
        if (cumulative_depth * 2 >= total_depth) return item.second;
    }
    return 0.0;
}

// Forward declaration: the record adapter below is kept next to the
// host-side weighted-median reduction, while the batch implementation lives
// with the other Kokkos posterior kernels above.
std::vector<double> contamination_probability_batch_kokkos(
    const std::vector<int>& total_counts, const std::vector<int>& alt_counts,
    const std::vector<double>& contaminations,
    const std::vector<double>& negative_log10_population_afs,
    const std::vector<double>& log_somatic_priors,
    const SomaticModelParameters& model, FilterKernelTelemetry* telemetry);

// Reconstruct the per-ALT ContaminationFilter posterior once per record.  The
// VCF/FORMAT decode and weighted-median reduction remain on Host, while the
// actual SomaticClusteringModel/binomial evaluation is delegated to the
// shared Kokkos batch above.  Keeping this helper shared by the empirical
// ThresholdCalculator pass and the final filtering pass is important: Java
// ErrorProbabilities sees the same ContaminationFilter evidence in both
// phases, and a duplicated adapter would otherwise silently diverge.
std::map<int, double> contamination_probabilities_for_record_kokkos(
    const bcf_hdr_t* header, bcf1_t* record, const float effective_contamination,
    const std::unordered_map<std::string, float>& contamination_by_sample,
    const std::vector<int>& tumor_sample_indices, const SomaticModelParameters& model,
    FilterKernelTelemetry* telemetry = nullptr, bool* evaluated = nullptr) {
    if (evaluated != nullptr) *evaluated = false;
    std::map<int, std::vector<std::pair<int, double>>> evidence_by_alt;
    if (header == nullptr || record == nullptr || effective_contamination < 0.0F ||
        bcf_hdr_id2int(header, BCF_DT_ID, "POPAF") < 0 || tumor_sample_indices.empty())
        return {};

    std::vector<int> totals;
    std::vector<int> alts;
    std::vector<double> contaminations;
    std::vector<double> population_afs;
    std::vector<double> somatic_priors;
    std::vector<int> alt_indices;
    for (const auto sample_index : tumor_sample_indices) {
        const auto total_depth = format_total_depth(header, record, sample_index);
        if (total_depth < 0) continue;
        double sample_contamination = static_cast<double>(effective_contamination);
        if (sample_index >= 0 && sample_index < bcf_hdr_nsamples(header)) {
            const auto* sample_name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample_index);
            if (sample_name != nullptr) {
                const auto estimate = contamination_by_sample.find(sample_name);
                if (estimate != contamination_by_sample.end())
                    sample_contamination = estimate->second;
            }
        }
        if (!std::isfinite(sample_contamination) || sample_contamination < 0.0) continue;
        for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
            const auto* allele = record->d.allele[alt_index];
            if (allele == nullptr || allele[0] == '<') continue;
            const auto alt_depth = format_alt_depth(header, record, sample_index, alt_index);
            if (alt_depth < 0) continue;
            const auto population_af = info_probability_at(header, record, "POPAF", alt_index);
            // GATK's POPAF accessor maps a missing Number=A element to
            // +infinity, and 10^-infinity is therefore the zero AF limit.
            const auto population_af_for_kernel = std::isfinite(population_af)
                ? population_af : 1.0e6;
            const auto reference_length = record->d.allele[0] == nullptr
                ? 0 : static_cast<int>(std::strlen(record->d.allele[0]));
            const auto alternate_length = static_cast<int>(std::strlen(allele));
            const auto indel_length = alternate_length - reference_length;
            double log_somatic_prior = alternate_length == reference_length
                ? model.log_snv_prior - std::log(3.0) : model.log_indel_prior;
            if (indel_length >= -10 && indel_length <= 10 &&
                std::isfinite(model.log_variant_priors[
                    static_cast<std::size_t>(indel_length + 10)])) {
                log_somatic_prior = model.log_variant_priors[
                    static_cast<std::size_t>(indel_length + 10)] +
                    (indel_length == 0 ? -std::log(3.0) : 0.0);
            }
            totals.push_back(total_depth);
            alts.push_back(alt_depth);
            contaminations.push_back(sample_contamination);
            population_afs.push_back(population_af_for_kernel);
            somatic_priors.push_back(log_somatic_prior);
            alt_indices.push_back(alt_index);
        }
    }
    if (totals.empty()) return {};
    const auto probabilities = contamination_probability_batch_kokkos(
        totals, alts, contaminations, population_afs, somatic_priors, model, telemetry);
    for (std::size_t offset = 0; offset < probabilities.size(); ++offset) {
        if (!std::isfinite(probabilities[offset])) continue;
        evidence_by_alt[alt_indices[offset]].emplace_back(alts[offset], probabilities[offset]);
    }
    std::map<int, double> result;
    for (auto& entry : evidence_by_alt)
        result[entry.first] = weighted_median_posterior(std::move(entry.second));
    if (evaluated != nullptr) *evaluated = !result.empty();
    return result;
}

double calculate_optimal_f_score_threshold(std::vector<double> posteriors, double beta) {
    if (posteriors.empty()) return 0.0;
    std::sort(posteriors.begin(), posteriors.end());
    double expected_true_positives = 0.0;
    for (const auto posterior : posteriors) expected_true_positives += 1.0 - posterior;
    double true_positives = 0.0;
    double false_positives = 0.0;
    double false_negatives = expected_true_positives;
    int optimal_index = -1;
    double optimal_f_score = 0.0;
    const double beta_squared = beta * beta;
    for (int index = 0; index < static_cast<int>(posteriors.size()); ++index) {
        const auto posterior = posteriors[static_cast<std::size_t>(index)];
        true_positives += 1.0 - posterior;
        false_positives += posterior;
        false_negatives -= 1.0 - posterior;
        const auto denominator = (1.0 + beta_squared) * true_positives +
                                 beta_squared * false_negatives + false_positives;
        const auto score = denominator > 0.0
            ? (1.0 + beta_squared) * true_positives / denominator : 0.0;
        if (score >= optimal_f_score) {
            optimal_index = index;
            optimal_f_score = score;
        }
    }
    if (optimal_index < 0) return 0.0;
    if (optimal_index == static_cast<int>(posteriors.size()) - 1) return 1.0;
    return posteriors[static_cast<std::size_t>(optimal_index)];
}

double calculate_false_discovery_rate_threshold(std::vector<double> posteriors,
                                                double requested_fdr) {
    if (posteriors.empty()) return 1.0;
    std::sort(posteriors.begin(), posteriors.end());
    double cumulative_expected_false_positives = 0.0;
    for (std::size_t index = 0; index < posteriors.size(); ++index) {
        const auto posterior = posteriors[index];
        const auto expected_fdr = (cumulative_expected_false_positives + posterior) /
                                  static_cast<double>(index + 1);
        if (expected_fdr > requested_fdr)
            return index > 0 ? posteriors[index - 1] : 0.0;
        cumulative_expected_false_positives += posterior;
    }
    return 1.0;
}

double learn_orientation_threshold(const Options& options, const std::vector<double>& posteriors) {
    if (!options.threshold_strategy_explicit)
        return options.max_orientation_artifact_probability >= 0.0F
            ? options.max_orientation_artifact_probability : 0.5;
    if (options.orientation_threshold_strategy == "CONSTANT")
        return options.max_orientation_artifact_probability >= 0.0F
            ? options.max_orientation_artifact_probability : options.initial_threshold;
    if (options.orientation_threshold_strategy == "FALSE_DISCOVERY_RATE")
        return calculate_false_discovery_rate_threshold(posteriors, options.max_false_discovery_rate);
    return calculate_optimal_f_score_threshold(posteriors, options.f_score_beta);
}

// FilterMutectCalls combines independent error *types* multiplicatively:
// 1 - prod(1 - p_type).  Filters within one ErrorType are correlated and are
// reduced with max() before the cross-type product.  The native VCF exposes
// the already-calibrated posterior fields produced by Mutect2 (and, for older
// records, TLOD/AD), so this helper reconstructs that ErrorProbabilities
// boundary without inventing a Java-only filter map.  A negative/non-finite
// argument means that the corresponding type was not present in the input
// record.
KOKKOS_INLINE_FUNCTION double kokkos_combine_error_probabilities(
    const double sequencing_error, const double non_somatic_error,
    const double artifact_error, const double orientation_error) {
    double probability_of_no_error = 1.0;
    bool observed = false;
    const auto consume = [&](const double value) {
        if (Kokkos::isfinite(value) && value >= 0.0) {
            probability_of_no_error *= 1.0 - Kokkos::min(1.0, Kokkos::max(0.0, value));
            observed = true;
        }
    };
    consume(sequencing_error);
    consume(non_somatic_error);
    consume(artifact_error);
    consume(orientation_error);
    return observed ? 1.0 - probability_of_no_error : -1.0;
}

double combine_error_probabilities_kokkos(
    const double sequencing_error, const double non_somatic_error,
    const double artifact_error, const double orientation_error,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_combined_error",
        kokkos_combine_error_probabilities(
            sequencing_error, non_somatic_error, artifact_error, orientation_error),
        telemetry);
}

std::vector<double> combine_error_probabilities_batch_kokkos(
    const std::vector<double>& sequencing_errors,
    const std::vector<double>& non_somatic_errors,
    const std::vector<double>& artifact_errors,
    const std::vector<double>& orientation_errors,
    FilterKernelTelemetry* telemetry = nullptr) {
    const auto count = sequencing_errors.size();
    if (non_somatic_errors.size() != count || artifact_errors.size() != count ||
        orientation_errors.size() != count)
        throw std::invalid_argument("internal FilterMutectCalls error batch size mismatch");
    if (count == 0) return {};
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("filter-mutect-error-probability-batch-v1");
    host_batch.records = count;
    host_batch.bytes = count * sizeof(double) * 5;
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("filter-mutect-error-probability-batch");
    kernel_plan.begin_prepare(host_batch);
    Kokkos::View<double*> sequencing_view("filter_mutect_batch_sequencing", count);
    Kokkos::View<double*> non_somatic_view("filter_mutect_batch_non_somatic", count);
    Kokkos::View<double*> artifact_view("filter_mutect_batch_artifact", count);
    Kokkos::View<double*> orientation_view("filter_mutect_batch_orientation", count);
    Kokkos::View<double*> result_view("filter_mutect_batch_error", count);
    auto host_sequencing = Kokkos::create_mirror_view(sequencing_view);
    auto host_non_somatic = Kokkos::create_mirror_view(non_somatic_view);
    auto host_artifact = Kokkos::create_mirror_view(artifact_view);
    auto host_orientation = Kokkos::create_mirror_view(orientation_view);
    for (std::size_t index = 0; index < count; ++index) {
        host_sequencing(index) = sequencing_errors[index];
        host_non_somatic(index) = non_somatic_errors[index];
        host_artifact(index) = artifact_errors[index];
        host_orientation(index) = orientation_errors[index];
    }
    Kokkos::deep_copy(sequencing_view, host_sequencing);
    Kokkos::deep_copy(non_somatic_view, host_non_somatic);
    Kokkos::deep_copy(artifact_view, host_artifact);
    Kokkos::deep_copy(orientation_view, host_orientation);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(count);
    device_batch.bind("sequencing", sequencing_view);
    device_batch.bind("non_somatic", non_somatic_view);
    device_batch.bind("artifact", artifact_view);
    device_batch.bind("orientation", orientation_view);
    device_batch.bind("result", result_view);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "filter_mutect_error_probability_batch",
        Kokkos::RangePolicy<ExecSpace>(0, count),
        KOKKOS_LAMBDA(const std::size_t index) {
            result_view(index) = kokkos_combine_error_probabilities(
                sequencing_view(index), non_somatic_view(index), artifact_view(index),
                orientation_view(index));
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) telemetry->record_batch(kernel_plan.telemetry(), count);
    const auto host_result = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, result_view);
    std::vector<double> result(count);
    for (std::size_t index = 0; index < count; ++index) result[index] = host_result(index);
    return result;
}

std::string reference_context(faidx_t* reference, const bcf_hdr_t* header, const bcf1_t* record,
                              int offset = 0) {
    if (reference == nullptr || record->rid < 0) return {};
    const auto* contig = bcf_hdr_id2name(header, record->rid);
    if (contig == nullptr || record->pos + offset < 1) return {};
    int length = 0;
    char* bases = faidx_fetch_seq(reference, contig, record->pos + offset - 1,
                                  record->pos + offset + 1, &length);
    if (bases == nullptr || length != 3) {
        free(bases);
        return {};
    }
    std::string context(bases, static_cast<std::size_t>(length));
    free(bases);
    for (auto& base : context) base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    for (const auto base : context) {
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') return {};
    }
    return context;
}

KOKKOS_INLINE_FUNCTION double kokkos_log_beta_binomial_probability(
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

KOKKOS_INLINE_FUNCTION double kokkos_log_binomial_mass(
    const int trials, const int successes, const double probability) {
    if (trials < 0 || successes < 0 || successes > trials ||
        probability < 0.0 || probability > 1.0)
        return -1.0e300;
    if (probability == 0.0) return successes == 0 ? 0.0 : -1.0e300;
    if (probability == 1.0) return successes == trials ? 0.0 : -1.0e300;
    return Kokkos::lgamma(static_cast<double>(trials) + 1.0) -
           Kokkos::lgamma(static_cast<double>(successes) + 1.0) -
           Kokkos::lgamma(static_cast<double>(trials - successes) + 1.0) +
           static_cast<double>(successes) * Kokkos::log(probability) +
           static_cast<double>(trials - successes) * Kokkos::log(1.0 - probability);
}

KOKKOS_INLINE_FUNCTION double kokkos_log_beta_probability(
    const double alpha, const double beta) {
    if (!(alpha > 0.0) || !(beta > 0.0)) return -1.0e300;
    return Kokkos::lgamma(alpha) + Kokkos::lgamma(beta) - Kokkos::lgamma(alpha + beta);
}

KOKKOS_INLINE_FUNCTION double kokkos_log_sum_exp_pair(
    const double first, const double second) {
    const auto maximum = first > second ? first : second;
    if (maximum < -1.0e299) return -1.0e300;
    return maximum + Kokkos::log(Kokkos::exp(first - maximum) +
                                 Kokkos::exp(second - maximum));
}

KOKKOS_INLINE_FUNCTION double kokkos_regularized_beta_continued_fraction(
    const double a, const double b, const double x) {
    // Numerical Recipes' modified Lentz recurrence.  The iteration cap is
    // deterministic across Kokkos backends; normal Mutect AD depths converge
    // in a few dozen steps, while the cap keeps pathological records finite.
    constexpr int max_iterations = 256;
    constexpr double epsilon = 3.0e-14;
    constexpr double floor = 1.0e-300;
    double c = 1.0;
    double d = 1.0 - (a + b) * x / (a + 1.0);
    if (Kokkos::abs(d) < floor) d = floor;
    d = 1.0 / d;
    double h = d;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        const double m = static_cast<double>(iteration);
        const double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((a - 1.0 + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (Kokkos::abs(d) < floor) d = floor;
        c = 1.0 + aa / c;
        if (Kokkos::abs(c) < floor) c = floor;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1.0 + m2));
        d = 1.0 + aa * d;
        if (Kokkos::abs(d) < floor) d = floor;
        c = 1.0 + aa / c;
        if (Kokkos::abs(c) < floor) c = floor;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (Kokkos::abs(delta - 1.0) <= epsilon) break;
    }
    return h;
}

KOKKOS_INLINE_FUNCTION double kokkos_regularized_beta(
    const double x, const double a, const double b) {
    if (!(a > 0.0) || !(b > 0.0) || x < 0.0 || x > 1.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x == 1.0) return 1.0;
    const double log_bt = kokkos_log_beta_probability(a, b) +
                          a * Kokkos::log(x) + b * Kokkos::log(1.0 - x);
    const double bt = Kokkos::exp(log_bt);
    const double cutoff = (a + 1.0) / (a + b + 2.0);
    if (x < cutoff)
        return bt * kokkos_regularized_beta_continued_fraction(a, b, x) / a;
    return 1.0 - bt * kokkos_regularized_beta_continued_fraction(b, a, 1.0 - x) / b;
}

KOKKOS_INLINE_FUNCTION double kokkos_posterior_error(
    const double log_odds_of_real, const double log_prior_of_real) {
    const double log_error_prior = Kokkos::log1p(-Kokkos::exp(log_prior_of_real));
    const double maximum = log_odds_of_real + log_prior_of_real > log_error_prior
        ? log_odds_of_real + log_prior_of_real : log_error_prior;
    if (!(maximum > -1.0e299)) return 0.0;
    const double denominator = Kokkos::exp(log_odds_of_real + log_prior_of_real - maximum) +
                               Kokkos::exp(log_error_prior - maximum);
    return Kokkos::exp(log_error_prior - maximum) / denominator;
}

KOKKOS_INLINE_FUNCTION double kokkos_tumor_evidence_probability(
    const double tumor_log_odds, const int total_count, const int alt_count) {
    if (total_count < 0 || alt_count < 0 || alt_count > total_count ||
        !Kokkos::isfinite(tumor_log_odds)) return 0.0;
    const int ref_count = total_count - alt_count;
    const double flat_correction =
        -kokkos_log_beta_probability(1.0, 1.0) +
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count)) +
        kokkos_log_beta_probability(1.0, 1.0) -
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count));
    const double high_correction =
        -kokkos_log_beta_probability(10.0, 1.0) +
        kokkos_log_beta_probability(10.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count)) +
        kokkos_log_beta_probability(1.0, 1.0) -
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count));
    const double log_likelihood = kokkos_log_sum_exp_pair(
        Kokkos::log1p(0.01) + tumor_log_odds + flat_correction,
        Kokkos::log(0.01) + tumor_log_odds + high_correction);
    // DEFAULT_INITIAL_LOG_PRIOR_OF_VARIANT_VERSUS_ARTIFACT = log10(-1).
    return kokkos_posterior_error(log_likelihood, -Kokkos::log(10.0));
}

KOKKOS_INLINE_FUNCTION double kokkos_tumor_evidence_probability_with_prior(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const double log_somatic_prior) {
    if (total_count < 0 || alt_count < 0 || alt_count > total_count ||
        !Kokkos::isfinite(tumor_log_odds) || !Kokkos::isfinite(log_somatic_prior)) return 0.0;
    const int ref_count = total_count - alt_count;
    const double flat_correction =
        -kokkos_log_beta_probability(1.0, 1.0) +
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count)) +
        kokkos_log_beta_probability(1.0, 1.0) -
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count));
    const double high_correction =
        -kokkos_log_beta_probability(10.0, 1.0) +
        kokkos_log_beta_probability(10.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count)) +
        kokkos_log_beta_probability(1.0, 1.0) -
        kokkos_log_beta_probability(1.0 + static_cast<double>(alt_count),
                                    1.0 + static_cast<double>(ref_count));
    const double log_likelihood = kokkos_log_sum_exp_pair(
        Kokkos::log1p(0.01) + tumor_log_odds + flat_correction,
        Kokkos::log(0.01) + tumor_log_odds + high_correction);
    return kokkos_posterior_error(log_likelihood, log_somatic_prior);
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_cluster_log_likelihood(
    const int total_count, const int alt_count, const double tumor_log_odds,
    const int cluster, const SomaticModelParameters& model) {
    if (cluster < 0 || cluster >= model.cluster_count || total_count < 0 || alt_count < 0 ||
        alt_count > total_count || !Kokkos::isfinite(tumor_log_odds)) return -1.0e300;
    const double alpha = model.cluster_alpha[static_cast<std::size_t>(cluster)];
    const double beta = model.cluster_beta[static_cast<std::size_t>(cluster)];
    if (!(alpha > 0.0) || !(beta > 0.0)) return -1.0e300;
    const double ref_count = static_cast<double>(total_count - alt_count);
    const double flat_log_beta = kokkos_log_beta_probability(1.0, 1.0);
    const double observed_flat_log_beta = kokkos_log_beta_probability(
        1.0 + static_cast<double>(alt_count), 1.0 + ref_count);
    const double correction = -kokkos_log_beta_probability(alpha, beta) +
        kokkos_log_beta_probability(alpha + static_cast<double>(alt_count), beta + ref_count) +
        flat_log_beta - observed_flat_log_beta;
    return tumor_log_odds + correction;
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_cluster_mixture_log_likelihood(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const SomaticModelParameters& model) {
    if (!model.learned || model.cluster_count < 1 || model.cluster_count > SomaticModelParameters::kMaxClusters)
        return -1.0e300;
    double maximum = -1.0e300;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const double candidate = model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
            kokkos_empirical_cluster_log_likelihood(total_count, alt_count, tumor_log_odds, cluster, model);
        if (candidate > maximum) maximum = candidate;
    }
    if (!(maximum > -1.0e299)) return -1.0e300;
    double sum = 0.0;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const double candidate = model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
            kokkos_empirical_cluster_log_likelihood(total_count, alt_count, tumor_log_odds, cluster, model);
        sum += Kokkos::exp(candidate - maximum);
    }
    return maximum + Kokkos::log(sum);
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_somatic_log_likelihood(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const SomaticModelParameters& model) {
    if (!model.learned || total_count < 0 || alt_count < 0 || alt_count > total_count ||
        !Kokkos::isfinite(tumor_log_odds)) return -1.0e300;
    return kokkos_empirical_cluster_mixture_log_likelihood(
        tumor_log_odds, total_count, alt_count, model);
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_somatic_count_log_likelihood(
    const int total_count, const int alt_count, const SomaticModelParameters& model) {
    if (!model.learned || total_count < 0 || alt_count < 0 || alt_count > total_count)
        return -1.0e300;
    if (model.cluster_count < 1 || model.cluster_count > SomaticModelParameters::kMaxClusters)
        return -1.0e300;
    double maximum = -1.0e300;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const double candidate = model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
            kokkos_log_beta_binomial_probability(
                alt_count, total_count,
                model.cluster_alpha[static_cast<std::size_t>(cluster)],
                model.cluster_beta[static_cast<std::size_t>(cluster)]);
        if (candidate > maximum) maximum = candidate;
    }
    if (!(maximum > -1.0e299)) return -1.0e300;
    double sum = 0.0;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const double candidate = model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
            kokkos_log_beta_binomial_probability(
                alt_count, total_count,
                model.cluster_alpha[static_cast<std::size_t>(cluster)],
                model.cluster_beta[static_cast<std::size_t>(cluster)]);
        sum += Kokkos::exp(candidate - maximum);
    }
    return maximum + Kokkos::log(sum);
}

// SomaticClusteringModel::logLikelihoodGivenSomatic() is available even
// before the optional empirical-learning pass.  Keep that distinction here:
// the existing empirical helper intentionally returns -inf when the model is
// not learned, while contamination filtering must still use GATK's initial
// two-cluster beta-binomial model.
KOKKOS_INLINE_FUNCTION double kokkos_model_somatic_count_log_likelihood(
    const int total_count, const int alt_count, const SomaticModelParameters& model) {
    if (total_count < 0 || alt_count < 0 || alt_count > total_count ||
        model.cluster_count < 1 || model.cluster_count > SomaticModelParameters::kMaxClusters)
        return -1.0e300;
    double maximum = -1.0e300;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const auto index = static_cast<std::size_t>(cluster);
        const auto mean = Kokkos::min(1.0, Kokkos::max(0.0, model.cluster_mean[index]));
        const auto cluster_log_likelihood = model.cluster_binomial[index] != 0
            ? kokkos_log_binomial_mass(total_count, alt_count, mean)
            : kokkos_log_beta_binomial_probability(
                alt_count, total_count, model.cluster_alpha[index], model.cluster_beta[index]);
        maximum = Kokkos::max(maximum,
            model.cluster_log_weights[index] + cluster_log_likelihood);
    }
    if (!(maximum > -1.0e299)) return -1.0e300;
    double sum = 0.0;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster) {
        const auto index = static_cast<std::size_t>(cluster);
        const auto mean = Kokkos::min(1.0, Kokkos::max(0.0, model.cluster_mean[index]));
        const auto cluster_log_likelihood = model.cluster_binomial[index] != 0
            ? kokkos_log_binomial_mass(total_count, alt_count, mean)
            : kokkos_log_beta_binomial_probability(
                alt_count, total_count, model.cluster_alpha[index], model.cluster_beta[index]);
        const auto candidate = model.cluster_log_weights[index] + cluster_log_likelihood;
        if (candidate > -1.0e299) sum += Kokkos::exp(candidate - maximum);
    }
    return sum > 0.0 ? maximum + Kokkos::log(sum) : -1.0e300;
}

// GATK GermlineFilter.germlineProbability(), kept as a device-callable
// boundary so the caller can materialize the same posterior that is exposed
// as INFO/GERMQ.  The two germline hypotheses (het and hom-alt) share the
// normal odds, while the somatic hypothesis carries the complementary
// population prior.  This is deliberately separate from PGERMLINE: that
// annotation is an optional upstream Mutect2 posterior and is not the
// GATK GermlineFilter calculation.
KOKKOS_INLINE_FUNCTION double kokkos_germline_probability(
    const double normal_log_odds,
    const double log_odds_of_germline_het_vs_somatic,
    const double log_odds_of_germline_hom_alt_vs_somatic,
    const double population_af,
    const double log_prior_somatic) {
    if (!Kokkos::isfinite(normal_log_odds) ||
        !Kokkos::isfinite(log_odds_of_germline_het_vs_somatic) ||
        !Kokkos::isfinite(population_af) ||
        population_af < 0.0 || population_af > 1.0)
        return 0.0;
    if (population_af < 1.0e-10) return 0.0;
    if (population_af > 1.0 - 1.0e-10) return 1.0;
    const double log_prior_not_somatic = Kokkos::log1p(-Kokkos::exp(log_prior_somatic));
    const double log_prior_germline_het = Kokkos::log(
        2.0 * population_af * (1.0 - population_af));
    const double log_prior_germline_hom_alt = Kokkos::log(
        population_af * population_af);
    const double log_prior_not_germline = Kokkos::log(
        (1.0 - population_af) * (1.0 - population_af));
    const double log_prob_germline_het = log_prior_germline_het +
        log_odds_of_germline_het_vs_somatic + normal_log_odds + log_prior_not_somatic;
    const double log_prob_germline_hom_alt = log_prior_germline_hom_alt +
        log_odds_of_germline_hom_alt_vs_somatic + normal_log_odds + log_prior_not_somatic;
    const double log_prob_germline = kokkos_log_sum_exp_pair(
        log_prob_germline_het, log_prob_germline_hom_alt);
    const double log_prob_somatic = log_prior_not_germline + log_prior_somatic;
    const double maximum = log_prob_germline > log_prob_somatic
        ? log_prob_germline : log_prob_somatic;
    if (!(maximum > -1.0e299)) return 0.0;
    const double germline_weight = Kokkos::exp(log_prob_germline - maximum);
    const double somatic_weight = Kokkos::exp(log_prob_somatic - maximum);
    const double normalizer = germline_weight + somatic_weight;
    return normalizer > 0.0 ? germline_weight / normalizer : 0.0;
}

double germline_probability_kokkos(
    const double normal_log_odds,
    const double log_odds_of_germline_het_vs_somatic,
    const double log_odds_of_germline_hom_alt_vs_somatic,
    const double population_af,
    const double log_prior_somatic,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_germline_probability",
        kokkos_germline_probability(
            normal_log_odds, log_odds_of_germline_het_vs_somatic,
            log_odds_of_germline_hom_alt_vs_somatic, population_af,
            log_prior_somatic),
        telemetry);
}

// Complete the bounded no-segmentation GermlineFilter path from the counts
// already decoded by Host.  The allele with the largest TLOD is selected by
// the caller; Kokkos evaluates the germline-het binomial likelihood, the
// initial/learned somatic beta-binomial mixture and the final posterior.
KOKKOS_INLINE_FUNCTION double kokkos_germline_probability_from_counts(
    const double normal_log10_odds,
    const int total_count,
    const int alt_count,
    const double weighted_tumor_af,
    const double minor_allele_fraction,
    const double population_af,
    const double log_prior_somatic,
    const SomaticModelParameters& model) {
    if (total_count < 0 || alt_count < 0 || alt_count > total_count ||
        !Kokkos::isfinite(normal_log10_odds) ||
        !Kokkos::isfinite(weighted_tumor_af) ||
        !Kokkos::isfinite(minor_allele_fraction) ||
        !Kokkos::isfinite(population_af)) return 0.0;
    const double log_germline_likelihood = -Kokkos::log(2.0) +
        kokkos_log_sum_exp_pair(
            kokkos_log_binomial_mass(
                total_count, alt_count, minor_allele_fraction),
            kokkos_log_binomial_mass(
                total_count, alt_count, 1.0 - minor_allele_fraction));
    const double log_somatic_likelihood =
        kokkos_model_somatic_count_log_likelihood(total_count, alt_count, model);
    if (!(log_somatic_likelihood > -1.0e299)) return 0.0;
    const double log_odds_het = log_germline_likelihood - log_somatic_likelihood;
    const double log_odds_hom = weighted_tumor_af < 0.9
        ? -std::numeric_limits<double>::infinity() : 0.0;
    return kokkos_germline_probability(
        -normal_log10_odds * Kokkos::log(10.0), log_odds_het,
        log_odds_hom, population_af, log_prior_somatic);
}

double germline_probability_from_counts_kokkos(
    const double normal_log10_odds,
    const int total_count,
    const int alt_count,
    const double weighted_tumor_af,
    const double minor_allele_fraction,
    const double population_af,
    const double log_prior_somatic,
    const SomaticModelParameters& model,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_germline_from_counts",
        kokkos_germline_probability_from_counts(
            normal_log10_odds, total_count, alt_count, weighted_tumor_af,
            minor_allele_fraction, population_af, log_prior_somatic, model),
        telemetry);
}

KOKKOS_INLINE_FUNCTION double kokkos_contamination_probability(
    const int total_count, const int alt_count, const double contamination,
    const double negative_log10_population_af, const double log_somatic_prior,
    const SomaticModelParameters& model) {
    if (total_count < 0 || alt_count < 0 || alt_count > total_count ||
        !Kokkos::isfinite(contamination) || !Kokkos::isfinite(negative_log10_population_af) ||
        !Kokkos::isfinite(log_somatic_prior))
        return std::numeric_limits<double>::quiet_NaN();
    const auto bounded_contamination = Kokkos::min(1.0 - 1.0e-10,
                                                    Kokkos::max(0.0, contamination));
    const auto population_af = Kokkos::min(1.0,
                                             Kokkos::max(0.0,
                                                         Kokkos::pow(10.0, -negative_log10_population_af)));
    const auto somatic = kokkos_model_somatic_count_log_likelihood(total_count, alt_count, model);
    if (!(somatic > -1.0e299)) return std::numeric_limits<double>::quiet_NaN();
    const auto log_binomial_half = kokkos_log_binomial_mass(
        total_count, alt_count, bounded_contamination * 0.5);
    const auto log_binomial_full = kokkos_log_binomial_mass(
        total_count, alt_count, bounded_contamination);
    const auto single_weight = 2.0 * population_af * (1.0 - population_af);
    const auto homozygous_weight = population_af * population_af;
    const auto log_single_first = single_weight > 0.0
        ? Kokkos::log(single_weight) + log_binomial_half : -1.0e300;
    const auto log_single_second = homozygous_weight > 0.0
        ? Kokkos::log(homozygous_weight) + log_binomial_full : -1.0e300;
    const auto log_single = kokkos_log_sum_exp_pair(log_single_first, log_single_second);
    const auto mixed_probability = bounded_contamination * population_af;
    const auto log_many = kokkos_log_binomial_mass(total_count, alt_count, mixed_probability);
    const auto log_contamination = Kokkos::max(log_single, log_many);
    if (!(log_contamination > -1.0e299)) return 0.0;
    return kokkos_posterior_error(somatic - log_contamination, log_somatic_prior);
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_tumor_error_probability(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const double log_somatic_prior, const SomaticModelParameters& model) {
    return kokkos_posterior_error(
        kokkos_empirical_somatic_log_likelihood(tumor_log_odds, total_count, alt_count, model),
        log_somatic_prior);
}

// The E-step is a small POD so the complete posterior remains device-owned
// while the M-step can retain GATK's Host-side, cross-pass optimizer.  Keeping
// the directional responsibilities (rather than only their sum) is necessary
// for StrandArtifactFilter.learnParameters().
struct StrandArtifactEStep {
    double forward_artifact_responsibility = 0.0;
    double reverse_artifact_responsibility = 0.0;
    int forward_count = 0;
    int reverse_count = 0;
    int forward_alt_count = 0;
    int reverse_alt_count = 0;
};

// StrandArtifactFilter evaluates three mutually exclusive hypotheses: forward
// strand artifact, reverse strand artifact, and non-artifact.  This is the
// direct GATK 4.6.2 calculation, including the combinatorial correction in
// the non-artifact likelihood.  The E-step is a Kokkos kernel just like the
// other FilterMutect error boundaries; its learned state is Host-owned.
KOKKOS_INLINE_FUNCTION StrandArtifactEStep kokkos_strand_artifact_e_step(
    const double artifact_prior, const double alpha_strand, const double beta_strand,
    const int forward_count, const int reverse_count,
    const int forward_alt_count, const int reverse_alt_count,
    const int indel_size) {
    StrandArtifactEStep result;
    result.forward_count = forward_count;
    result.reverse_count = reverse_count;
    result.forward_alt_count = forward_alt_count;
    result.reverse_alt_count = reverse_alt_count;
    if (!(artifact_prior >= 0.0 && artifact_prior < 1.0) ||
        !(alpha_strand > 0.0) || !(beta_strand > 0.0) ||
        forward_count < 0 || reverse_count < 0 || forward_alt_count < 0 ||
        reverse_alt_count < 0 || forward_alt_count > forward_count ||
        reverse_alt_count > reverse_count || indel_size > 4 ||
        forward_alt_count + reverse_alt_count == 0) return result;
    const double beta_sequencing = indel_size == 0 ? 1000.0
        : (indel_size < 3 ? 5000.0 : 50000.0);
    const double forward_log_likelihood =
        kokkos_log_beta_binomial_probability(
            forward_alt_count, forward_count, alpha_strand, beta_strand) +
        kokkos_log_beta_binomial_probability(
            reverse_alt_count, reverse_count, 1.0, beta_sequencing);
    const double reverse_log_likelihood =
        kokkos_log_beta_binomial_probability(
            reverse_alt_count, reverse_count, alpha_strand, beta_strand) +
        kokkos_log_beta_binomial_probability(
            forward_alt_count, forward_count, 1.0, beta_sequencing);
    const auto log_choose = [](const int trials, const int successes) {
        return Kokkos::lgamma(static_cast<double>(trials) + 1.0) -
            Kokkos::lgamma(static_cast<double>(successes) + 1.0) -
            Kokkos::lgamma(static_cast<double>(trials - successes) + 1.0);
    };
    const double none_log_likelihood =
        log_choose(forward_count, forward_alt_count) +
        log_choose(reverse_count, reverse_alt_count) -
        log_choose(forward_count + reverse_count,
                   forward_alt_count + reverse_alt_count) +
        kokkos_log_beta_binomial_probability(
            forward_alt_count + reverse_alt_count, forward_count + reverse_count,
            1.0, 1.0);
    const double forward = forward_log_likelihood + Kokkos::log(artifact_prior / 2.0);
    const double reverse = reverse_log_likelihood + Kokkos::log(artifact_prior / 2.0);
    const double none = none_log_likelihood + Kokkos::log1p(-artifact_prior);
    const double maximum = Kokkos::max(forward, Kokkos::max(reverse, none));
    const double forward_weight = Kokkos::exp(forward - maximum);
    const double reverse_weight = Kokkos::exp(reverse - maximum);
    const double none_weight = Kokkos::exp(none - maximum);
    const double normalizer = forward_weight + reverse_weight + none_weight;
    if (normalizer > 0.0) {
        result.forward_artifact_responsibility = forward_weight / normalizer;
        result.reverse_artifact_responsibility = reverse_weight / normalizer;
    }
    return result;
}

KOKKOS_INLINE_FUNCTION double kokkos_strand_artifact_probability(
    const double artifact_prior, const double alpha_strand, const double beta_strand,
    const int forward_count, const int reverse_count,
    const int forward_alt_count, const int reverse_alt_count,
    const int indel_size) {
    const auto e_step = kokkos_strand_artifact_e_step(
        artifact_prior, alpha_strand, beta_strand, forward_count, reverse_count,
        forward_alt_count, reverse_alt_count, indel_size);
    return e_step.forward_artifact_responsibility +
        e_step.reverse_artifact_responsibility;
}

StrandArtifactEStep strand_artifact_e_step_kokkos(
    const double artifact_prior, const double alpha_strand, const double beta_strand,
    const int forward_count, const int reverse_count,
    const int forward_alt_count, const int reverse_alt_count,
    const int indel_size, FilterKernelTelemetry* telemetry = nullptr) {
    using FilterExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("filter_mutect_strand_artifact_estep");
    host_batch.records = 1;
    host_batch.bytes = sizeof(StrandArtifactEStep);
    fastgatk::core::KernelPlan<FilterExecSpace> plan("filter_mutect_strand_artifact_estep");
    plan.begin_prepare(host_batch);
    Kokkos::View<StrandArtifactEStep*> result_device(
        "filter_mutect_strand_artifact_estep", 1);
    fastgatk::core::DeviceBatch<FilterExecSpace> device_batch(1);
    device_batch.bind("result", result_device);
    FilterExecSpace().fence();
    plan.end_prepare(device_batch);
    plan.begin_execute();
    Kokkos::parallel_for(
        "filter_mutect_strand_artifact_estep",
        Kokkos::RangePolicy<FilterExecSpace>(0, 1),
        KOKKOS_LAMBDA(const int) {
            result_device(0) = kokkos_strand_artifact_e_step(
                artifact_prior, alpha_strand, beta_strand, forward_count,
                reverse_count, forward_alt_count, reverse_alt_count, indel_size);
        });
    FilterExecSpace().fence();
    plan.end_execute();
    const auto result_host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, result_device);
    if (telemetry != nullptr) telemetry->record(plan.telemetry(), 1);
    return result_host(0);
}

double strand_artifact_probability_kokkos(
    const double artifact_prior, const double alpha_strand, const double beta_strand,
    const int forward_count, const int reverse_count,
    const int forward_alt_count, const int reverse_alt_count,
    const int indel_size, FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_strand_artifact",
        kokkos_strand_artifact_probability(
            artifact_prior, alpha_strand, beta_strand, forward_count, reverse_count,
            forward_alt_count, reverse_alt_count, indel_size), telemetry);
}

// Apache Commons Math's BrentOptimizer is what GATK calls from
// StrandArtifactFilter.  This compact Host transcription keeps the same
// bracket, tolerances, starting point and 100-evaluation ceiling.  The
// objective itself is a sum of beta-binomial log probabilities; its posterior
// responsibilities were produced by the Kokkos E-step above.
template <typename Objective>
double brent_maximize(Objective&& objective, const double lower, const double upper,
                      const double initial, const double relative_tolerance,
                      const double absolute_tolerance, const int max_evaluations) {
    double a = lower;
    double b = upper;
    double x = std::clamp(initial, lower, upper);
    double w = x;
    double v = x;
    double fx = objective(x);
    double fw = fx;
    double fv = fx;
    double d = 0.0;
    double e = 0.0;
    int evaluations = 1;
    constexpr double golden_section = 0.5 * (3.0 - 2.2360679774997896964);
    while (evaluations < max_evaluations) {
        const double midpoint = 0.5 * (a + b);
        const double tolerance1 = relative_tolerance * std::abs(x) + absolute_tolerance;
        const double tolerance2 = 2.0 * tolerance1;
        if (std::abs(x - midpoint) <= tolerance2 - 0.5 * (b - a)) break;

        bool accept_parabolic_step = false;
        if (std::abs(e) > tolerance1) {
            const double r = (x - w) * (fx - fv);
            const double q0 = (x - v) * (fx - fw);
            double p = (x - v) * q0 - (x - w) * r;
            double q = 2.0 * (q0 - r);
            if (q > 0.0) p = -p;
            q = std::abs(q);
            const double previous_e = e;
            e = d;
            if (q > 0.0 && p > q * (a - x) && p < q * (b - x) &&
                std::abs(p) < std::abs(0.5 * q * previous_e)) {
                d = p / q;
                const double candidate = x + d;
                if (candidate - a < tolerance2 || b - candidate < tolerance2)
                    d = std::copysign(tolerance1, midpoint - x);
                accept_parabolic_step = true;
            }
        }
        if (!accept_parabolic_step) {
            e = x < midpoint ? b - x : a - x;
            d = golden_section * e;
        }
        const double u = x + (std::abs(d) >= tolerance1
            ? d : std::copysign(tolerance1, d == 0.0 ? midpoint - x : d));
        const double fu = objective(u);
        ++evaluations;
        if (fu >= fx) {
            if (u < x) b = x; else a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu >= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu >= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return x;
}

class LearnedStrandArtifactFilter {
public:
    StrandArtifactEStep calculate_e_step(const int forward_count, const int reverse_count,
                                         const int forward_alt_count,
                                         const int reverse_alt_count,
                                         const int indel_size,
                                         FilterKernelTelemetry* telemetry) const {
        return strand_artifact_e_step_kokkos(
            prior_, alpha_, beta_, forward_count, reverse_count,
            forward_alt_count, reverse_alt_count, indel_size, telemetry);
    }

    double error_probability(const int forward_count, const int reverse_count,
                             const int forward_alt_count, const int reverse_alt_count,
                             const int indel_size,
                             FilterKernelTelemetry* telemetry) const {
        return strand_artifact_probability_kokkos(
            prior_, alpha_, beta_, forward_count, reverse_count,
            forward_alt_count, reverse_alt_count, indel_size, telemetry);
    }

    void accumulate(const StrandArtifactEStep& e_step) { e_steps_.push_back(e_step); }

    void learn_parameters_and_clear() {
        std::vector<const StrandArtifactEStep*> potential_artifacts;
        potential_artifacts.reserve(e_steps_.size());
        double total_artifacts = 0.0;
        double total_non_artifacts = 0.0;
        for (const auto& e_step : e_steps_) {
            const double artifact_probability =
                e_step.forward_artifact_responsibility +
                e_step.reverse_artifact_responsibility;
            total_non_artifacts += 1.0 - artifact_probability;
            if (artifact_probability > 0.1) {
                potential_artifacts.push_back(&e_step);
                total_artifacts += artifact_probability;
            }
        }
        prior_ = (total_artifacts + kArtifactPseudocount) /
            (total_artifacts + kArtifactPseudocount + total_non_artifacts +
             kNonArtifactPseudocount);

        double artifact_alt_count = 0.0;
        double artifact_depth = 0.0;
        for (const auto* e_step : potential_artifacts) {
            artifact_alt_count +=
                e_step->forward_artifact_responsibility * e_step->forward_alt_count +
                e_step->reverse_artifact_responsibility * e_step->reverse_alt_count;
            artifact_depth +=
                e_step->forward_artifact_responsibility * e_step->forward_count +
                e_step->reverse_artifact_responsibility * e_step->reverse_count;
        }
        const double mean = (artifact_alt_count + kInitialAlpha) /
            (artifact_depth + kInitialAlpha + kInitialBeta);
        const double safe_mean = std::clamp(mean, 1.0e-12, 1.0 - 1.0e-12);
        const auto objective = [&](const double alpha) {
            const double beta = (1.0 / safe_mean - 1.0) * alpha;
            double result = 0.0;
            for (const auto* e_step : potential_artifacts) {
                result += e_step->forward_artifact_responsibility *
                    kokkos_log_beta_binomial_probability(
                        e_step->forward_alt_count, e_step->forward_count, alpha, beta);
                result += e_step->reverse_artifact_responsibility *
                    kokkos_log_beta_binomial_probability(
                        e_step->reverse_alt_count, e_step->reverse_count, alpha, beta);
            }
            return result;
        };
        alpha_ = brent_maximize(objective, 0.01, 100.0, kInitialAlpha,
                                0.01, 0.01, 100);
        beta_ = (1.0 / safe_mean - 1.0) * alpha_;
        e_steps_.clear();
    }

    double prior() const { return prior_; }
    double alpha() const { return alpha_; }
    double beta() const { return beta_; }
    std::size_t accumulated_e_steps() const { return e_steps_.size(); }

private:
    static constexpr double kInitialAlpha = 1.0;
    static constexpr double kInitialBeta = 20.0;
    static constexpr double kArtifactPseudocount = 1.0;
    static constexpr double kNonArtifactPseudocount = 1000.0;
    double alpha_ = kInitialAlpha;
    double beta_ = kInitialBeta;
    double prior_ = 0.001;
    std::vector<StrandArtifactEStep> e_steps_;
};

KOKKOS_INLINE_FUNCTION double kokkos_slippage_probability(
    const double slippage_rate, const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length) {
    if (reference_repeats < 0 || alternate_repeats < 0 || repeat_unit_length < 0 ||
        total_count < 0 || alternate_count < 0 || alternate_count > total_count ||
        repeat_unit_length * reference_repeats < min_slippage_length ||
        Kokkos::abs(reference_repeats - alternate_repeats) != 1 ||
        !(slippage_rate >= 0.0 && slippage_rate <= 1.0)) return 0.0;
    const double somatic = kokkos_log_sum_exp_pair(
        Kokkos::log1p(0.01) + kokkos_log_beta_binomial_probability(
            alternate_count, total_count, 1.0, 1.0),
        Kokkos::log(0.01) + kokkos_log_beta_binomial_probability(
            alternate_count, total_count, 10.0, 1.0));
    const double artifact_probability = kokkos_regularized_beta(
        slippage_rate, static_cast<double>(alternate_count) + 1.0,
        static_cast<double>(total_count - alternate_count) + 1.0);
    if (!(artifact_probability > 0.0) || !Kokkos::isfinite(artifact_probability)) return 0.0;
    const int indel_length = alternate_length - reference_length;
    const double prior = indel_length == 0
        ? -6.0 * Kokkos::log(10.0) - Kokkos::log(3.0)
        : -7.0 * Kokkos::log(10.0);
    return kokkos_posterior_error(somatic - Kokkos::log(artifact_probability), prior);
}

KOKKOS_INLINE_FUNCTION double kokkos_slippage_probability_with_prior(
    const double slippage_rate, const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length, const double log_somatic_prior) {
    if (reference_repeats < 0 || alternate_repeats < 0 || repeat_unit_length < 0 ||
        total_count < 0 || alternate_count < 0 || alternate_count > total_count ||
        repeat_unit_length * reference_repeats < min_slippage_length ||
        Kokkos::abs(reference_repeats - alternate_repeats) != 1 ||
        !(slippage_rate >= 0.0 && slippage_rate <= 1.0)) return 0.0;
    const double somatic = kokkos_log_sum_exp_pair(
        Kokkos::log1p(0.01) + kokkos_log_beta_binomial_probability(
            alternate_count, total_count, 1.0, 1.0),
        Kokkos::log(0.01) + kokkos_log_beta_binomial_probability(
            alternate_count, total_count, 10.0, 1.0));
    const double artifact_probability = kokkos_regularized_beta(
        slippage_rate, static_cast<double>(alternate_count) + 1.0,
        static_cast<double>(total_count - alternate_count) + 1.0);
    if (!(artifact_probability > 0.0) || !Kokkos::isfinite(artifact_probability) ||
        !Kokkos::isfinite(log_somatic_prior)) return 0.0;
    return kokkos_posterior_error(somatic - Kokkos::log(artifact_probability), log_somatic_prior);
}

KOKKOS_INLINE_FUNCTION double kokkos_empirical_slippage_probability(
    const SomaticModelParameters& model, const double slippage_rate,
    const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length, const double log_somatic_prior) {
    if (!model.learned || reference_repeats < 0 || alternate_repeats < 0 ||
        repeat_unit_length < 0 || total_count < 0 || alternate_count < 0 ||
        alternate_count > total_count || repeat_unit_length * reference_repeats < min_slippage_length ||
        Kokkos::abs(reference_repeats - alternate_repeats) != 1 ||
        !(slippage_rate >= 0.0 && slippage_rate <= 1.0)) return 0.0;
    const auto somatic = kokkos_empirical_somatic_count_log_likelihood(
        total_count, alternate_count, model);
    const auto artifact_probability = kokkos_regularized_beta(
        slippage_rate, static_cast<double>(alternate_count) + 1.0,
        static_cast<double>(total_count - alternate_count) + 1.0);
    if (!(artifact_probability > 0.0) || !Kokkos::isfinite(artifact_probability)) return 0.0;
    const int indel_length = alternate_length - reference_length;
    (void)indel_length;
    return kokkos_posterior_error(somatic - Kokkos::log(artifact_probability), log_somatic_prior);
}

double tumor_evidence_probability_kokkos(const double tumor_log_odds,
                                         const int total_count, const int alt_count,
                                         FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_tumor_evidence",
        kokkos_tumor_evidence_probability(tumor_log_odds, total_count, alt_count),
        telemetry);
}

double tumor_evidence_probability_with_prior_kokkos(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const double log_somatic_prior, FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_tumor_evidence_prior",
        kokkos_tumor_evidence_probability_with_prior(
            tumor_log_odds, total_count, alt_count, log_somatic_prior),
        telemetry);
}

double empirical_tumor_evidence_probability_kokkos(
    const double tumor_log_odds, const int total_count, const int alt_count,
    const double log_somatic_prior,
    const SomaticModelParameters& model,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_empirical_tumor_evidence",
        kokkos_empirical_tumor_error_probability(
            tumor_log_odds, total_count, alt_count, log_somatic_prior, model),
        telemetry);
}

// Batch the per-ALT tumor-evidence boundary.  FilterMutectCalls keeps VCF
// interpretation and allele ordering on Host, but all ALT-local posterior
// math is independent and can be launched once for the record (or a larger
// shard in a future persistent plan).
std::vector<double> tumor_evidence_probability_batch_kokkos(
    const std::vector<double>& tumor_log_odds,
    const std::vector<int>& total_counts,
    const std::vector<int>& alt_counts,
    const std::vector<double>& log_somatic_priors,
    const SomaticModelParameters* empirical_model,
    FilterKernelTelemetry* telemetry = nullptr) {
    const auto count = tumor_log_odds.size();
    if (total_counts.size() != count || alt_counts.size() != count ||
        log_somatic_priors.size() != count)
        throw std::invalid_argument("internal FilterMutectCalls batch size mismatch");
    if (count == 0) return {};
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("filter-mutect-tumor-evidence-batch-v1");
    host_batch.records = count;
    host_batch.bytes = count * (sizeof(double) * 3 + sizeof(int) * 2);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("filter-mutect-tumor-evidence-batch");
    kernel_plan.begin_prepare(host_batch);
    Kokkos::View<double*> log_odds_view("filter_mutect_batch_log_odds", count);
    Kokkos::View<int*> total_view("filter_mutect_batch_total", count);
    Kokkos::View<int*> alt_view("filter_mutect_batch_alt", count);
    Kokkos::View<double*> prior_view("filter_mutect_batch_prior", count);
    Kokkos::View<double*> result_view("filter_mutect_batch_result", count);
    auto host_log_odds = Kokkos::create_mirror_view(log_odds_view);
    auto host_total = Kokkos::create_mirror_view(total_view);
    auto host_alt = Kokkos::create_mirror_view(alt_view);
    auto host_prior = Kokkos::create_mirror_view(prior_view);
    for (std::size_t index = 0; index < count; ++index) {
        host_log_odds(index) = tumor_log_odds[index];
        host_total(index) = total_counts[index];
        host_alt(index) = alt_counts[index];
        host_prior(index) = log_somatic_priors[index];
    }
    Kokkos::deep_copy(log_odds_view, host_log_odds);
    Kokkos::deep_copy(total_view, host_total);
    Kokkos::deep_copy(alt_view, host_alt);
    Kokkos::deep_copy(prior_view, host_prior);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(count);
    device_batch.bind("log_odds", log_odds_view);
    device_batch.bind("total", total_view);
    device_batch.bind("alt", alt_view);
    device_batch.bind("prior", prior_view);
    device_batch.bind("result", result_view);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    const bool use_empirical = empirical_model != nullptr && empirical_model->learned;
    const auto model = use_empirical ? *empirical_model : SomaticModelParameters{};
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "filter_mutect_tumor_evidence_batch",
        Kokkos::RangePolicy<ExecSpace>(0, count),
        KOKKOS_LAMBDA(const std::size_t index) {
            result_view(index) = use_empirical
                ? kokkos_empirical_tumor_error_probability(
                    log_odds_view(index), total_view(index), alt_view(index),
                    prior_view(index), model)
                : kokkos_tumor_evidence_probability_with_prior(
                    log_odds_view(index), total_view(index), alt_view(index),
                    prior_view(index));
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) telemetry->record_batch(kernel_plan.telemetry(), count);
    const auto host_result = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, result_view);
    std::vector<double> result(count);
    for (std::size_t index = 0; index < count; ++index) result[index] = host_result(index);
    return result;
}

// Per-tumour/per-ALT contamination posteriors are independent.  Keep VCF
// decoding and the weighted-median reduction on Host, but execute the
// SomaticClusteringModel + binomial likelihoods through one Kokkos launch.
std::vector<double> contamination_probability_batch_kokkos(
    const std::vector<int>& total_counts,
    const std::vector<int>& alt_counts,
    const std::vector<double>& contaminations,
    const std::vector<double>& negative_log10_population_afs,
    const std::vector<double>& log_somatic_priors,
    const SomaticModelParameters& model,
    FilterKernelTelemetry* telemetry = nullptr) {
    const auto count = total_counts.size();
    if (alt_counts.size() != count || contaminations.size() != count ||
        negative_log10_population_afs.size() != count || log_somatic_priors.size() != count)
        throw std::invalid_argument("internal FilterMutectCalls contamination batch size mismatch");
    if (count == 0) return {};
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("filter-mutect-contamination-batch-v1");
    host_batch.records = count;
    host_batch.bytes = count * (sizeof(double) * 3 + sizeof(int) * 2);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("filter-mutect-contamination-batch");
    kernel_plan.begin_prepare(host_batch);
    Kokkos::View<int*> total_view("filter_mutect_contamination_total", count);
    Kokkos::View<int*> alt_view("filter_mutect_contamination_alt", count);
    Kokkos::View<double*> contamination_view("filter_mutect_contamination_fraction", count);
    Kokkos::View<double*> popaf_view("filter_mutect_contamination_popaf", count);
    Kokkos::View<double*> prior_view("filter_mutect_contamination_prior", count);
    Kokkos::View<double*> result_view("filter_mutect_contamination_result", count);
    auto host_total = Kokkos::create_mirror_view(total_view);
    auto host_alt = Kokkos::create_mirror_view(alt_view);
    auto host_contamination = Kokkos::create_mirror_view(contamination_view);
    auto host_popaf = Kokkos::create_mirror_view(popaf_view);
    auto host_prior = Kokkos::create_mirror_view(prior_view);
    for (std::size_t index = 0; index < count; ++index) {
        host_total(index) = total_counts[index];
        host_alt(index) = alt_counts[index];
        host_contamination(index) = contaminations[index];
        host_popaf(index) = negative_log10_population_afs[index];
        host_prior(index) = log_somatic_priors[index];
    }
    Kokkos::deep_copy(total_view, host_total);
    Kokkos::deep_copy(alt_view, host_alt);
    Kokkos::deep_copy(contamination_view, host_contamination);
    Kokkos::deep_copy(popaf_view, host_popaf);
    Kokkos::deep_copy(prior_view, host_prior);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(count);
    device_batch.bind("total", total_view);
    device_batch.bind("alt", alt_view);
    device_batch.bind("contamination", contamination_view);
    device_batch.bind("popaf", popaf_view);
    device_batch.bind("prior", prior_view);
    device_batch.bind("result", result_view);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    const auto device_model = model;
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "filter_mutect_contamination_probability_batch",
        Kokkos::RangePolicy<ExecSpace>(0, count),
        KOKKOS_LAMBDA(const std::size_t index) {
            result_view(index) = kokkos_contamination_probability(
                total_view(index), alt_view(index), contamination_view(index),
                popaf_view(index), prior_view(index), device_model);
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) telemetry->record_batch(kernel_plan.telemetry(), count);
    const auto host_result = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, result_view);
    std::vector<double> result(count);
    for (std::size_t index = 0; index < count; ++index) result[index] = host_result(index);
    return result;
}

void initialize_somatic_model_clusters(SomaticModelParameters& model) {
    model.cluster_count = 2;
    model.cluster_log_weights.fill(-1.0e300);
    model.cluster_alpha.fill(1.0);
    model.cluster_beta.fill(1.0);
    model.cluster_mean.fill(0.5);
    model.cluster_binomial.fill(0);
    model.cluster_log_weights[0] = std::log1p(0.01);
    model.cluster_log_weights[1] = std::log(0.01);
    model.cluster_alpha[0] = model.background_alpha = 1.0;
    model.cluster_beta[0] = model.background_beta = 1.0;
    model.cluster_alpha[1] = model.high_af_alpha = 10.0;
    model.cluster_beta[1] = model.high_af_beta = 1.0;
    model.cluster_mean[0] = 0.5;
    model.cluster_mean[1] = 10.0 / 11.0;
    model.background_weight = std::exp(std::log1p(0.01));
    model.high_af_weight = 0.01;
}

double host_cluster_count_log_likelihood(const int total_count, const int alt_count,
                                         const int cluster, const SomaticModelParameters& model) {
    if (cluster < 0 || cluster >= model.cluster_count || total_count < 0 || alt_count < 0 ||
        alt_count > total_count) return -1.0e300;
    return kokkos_log_beta_binomial_probability(
        alt_count, total_count,
        model.cluster_alpha[static_cast<std::size_t>(cluster)],
        model.cluster_beta[static_cast<std::size_t>(cluster)]);
}

double host_model_count_log_likelihood(const int total_count, const int alt_count,
                                       const SomaticModelParameters& model) {
    double maximum = -1.0e300;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster)
        maximum = std::max(maximum, model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
            host_cluster_count_log_likelihood(total_count, alt_count, cluster, model));
    if (!(maximum > -1.0e299)) return -1.0e300;
    double sum = 0.0;
    for (int cluster = 0; cluster < model.cluster_count; ++cluster)
        sum += std::exp(model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
                        host_cluster_count_log_likelihood(total_count, alt_count, cluster, model) - maximum);
    return maximum + std::log(std::max(sum, 1.0e-300));
}

// Apache Commons Math's positive-argument digamma approximation used by the
// GATK beta-binomial learner.  Recurrence to x>=8 followed by the Bernoulli
// expansion is deterministic and device-callable if the learner is moved to
// a Kokkos TeamPolicy in a later backend implementation.
KOKKOS_INLINE_FUNCTION double filter_gatk_digamma(double value) {
    if (!(value > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    double result = 0.0;
    while (value < 8.0) {
        result -= 1.0 / value;
        value += 1.0;
    }
    const double inverse = 1.0 / value;
    const double inverse_squared = inverse * inverse;
    result += Kokkos::log(value) - 0.5 * inverse - inverse_squared * (
        1.0 / 12.0 - inverse_squared * (
            1.0 / 120.0 - inverse_squared * (
                1.0 / 252.0 - inverse_squared * (1.0 / 240.0))));
    return result;
}

void perform_somatic_model_em_iteration(
    SomaticModelParameters& model, const std::vector<SomaticModelDatum>& data) {
    const int cluster_count = model.cluster_count;
    std::array<double, SomaticModelParameters::kMaxClusters> responsibility_sum{};
    std::array<double, SomaticModelParameters::kMaxClusters> weighted_af_sum{};
    std::array<double, SomaticModelParameters::kMaxClusters> weighted_af2_sum{};
    std::vector<std::pair<std::size_t, std::array<double, SomaticModelParameters::kMaxClusters>>> datum_responsibilities;
    datum_responsibilities.reserve(data.size());
    for (const auto& datum : data) {
        if (datum.total_count <= 0 || datum.alt_count < 0 || datum.alt_count > datum.total_count) continue;
        std::array<double, SomaticModelParameters::kMaxClusters> log_responsibilities{};
        double maximum = -1.0e300;
        for (int cluster = 0; cluster < cluster_count; ++cluster) {
            log_responsibilities[static_cast<std::size_t>(cluster)] =
                model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
                host_cluster_count_log_likelihood(datum.total_count, datum.alt_count, cluster, model);
            maximum = std::max(maximum, log_responsibilities[static_cast<std::size_t>(cluster)]);
        }
        double normalization = 0.0;
        for (int cluster = 0; cluster < cluster_count; ++cluster)
            normalization += std::exp(log_responsibilities[static_cast<std::size_t>(cluster)] - maximum);
        if (!(normalization > 0.0) || !std::isfinite(normalization)) continue;
        const auto somatic_weight = std::clamp(datum.somatic_weight, 0.0, 1.0);
        const auto allele_fraction = static_cast<double>(datum.alt_count) /
                                     static_cast<double>(datum.total_count);
        std::array<double, SomaticModelParameters::kMaxClusters> current_responsibilities{};
        for (int cluster = 0; cluster < cluster_count; ++cluster) {
            const auto responsibility = somatic_weight *
                std::exp(log_responsibilities[static_cast<std::size_t>(cluster)] - maximum) / normalization;
            current_responsibilities[static_cast<std::size_t>(cluster)] = responsibility;
            responsibility_sum[static_cast<std::size_t>(cluster)] += responsibility;
            weighted_af_sum[static_cast<std::size_t>(cluster)] += responsibility * allele_fraction;
            weighted_af2_sum[static_cast<std::size_t>(cluster)] += responsibility * allele_fraction * allele_fraction;
        }
        datum_responsibilities.emplace_back(
            static_cast<std::size_t>(&datum - data.data()), current_responsibilities);
    }
    double total_weight = 0.0;
    for (int cluster = 0; cluster < cluster_count; ++cluster)
        total_weight += responsibility_sum[static_cast<std::size_t>(cluster)] + 1.0;
    if (!(total_weight > 0.0)) return;
    for (int cluster = 0; cluster < cluster_count; ++cluster) {
        const auto index = static_cast<std::size_t>(cluster);
        model.cluster_log_weights[index] = std::log(
            (responsibility_sum[index] + 1.0) / total_weight);
        if (!(responsibility_sum[index] > 1.0e-9)) continue;
        if (cluster < 2) {
            // Match BetaBinomialCluster.learn: fixed responsibilities, ten
            // online epochs, rate 0.01, and the Java lower bounds.
            double alpha = model.cluster_alpha[index];
            double beta = model.cluster_beta[index];
            for (int epoch = 0; epoch < 10; ++epoch) {
                for (const auto& entry : datum_responsibilities) {
                    const auto& datum = data[entry.first];
                    const auto responsibility = entry.second[index];
                    if (!(responsibility > 0.0) || datum.total_count < 0 ||
                        datum.alt_count < 0 || datum.alt_count > datum.total_count) continue;
                    const auto total_plus = static_cast<double>(datum.total_count) + alpha + beta;
                    const auto alpha_gradient = filter_gatk_digamma(
                        alpha + static_cast<double>(datum.alt_count)) - filter_gatk_digamma(total_plus) -
                        filter_gatk_digamma(alpha) + filter_gatk_digamma(alpha + beta);
                    const auto beta_gradient = filter_gatk_digamma(
                        beta + static_cast<double>(datum.total_count - datum.alt_count)) -
                        filter_gatk_digamma(total_plus) - filter_gatk_digamma(beta) +
                        filter_gatk_digamma(alpha + beta);
                    alpha = std::max(alpha + 0.01 * alpha_gradient * responsibility, 1.0);
                    beta = std::max(beta + 0.01 * beta_gradient * responsibility, 0.5);
                }
            }
            model.cluster_alpha[index] = alpha;
            model.cluster_beta[index] = beta;
            model.cluster_mean[index] = alpha / (alpha + beta);
        } else {
            const auto mean = std::clamp(weighted_af_sum[index] / responsibility_sum[index], 1.0e-4, 1.0 - 1.0e-4);
            // The Java BinomialCluster is represented as a very concentrated
            // beta distribution around its weighted AF mean.
            const auto concentration = std::clamp(
                (1.0 - mean) / (mean * 1.0e-4) - 1.0, 10.0, 1.0e6);
            model.cluster_mean[index] = mean;
            model.cluster_alpha[index] = std::max(mean * concentration, 0.1);
            model.cluster_beta[index] = std::max((1.0 - mean) * concentration, 0.1);
        }
    }
    model.background_alpha = model.cluster_alpha[0];
    model.background_beta = model.cluster_beta[0];
    model.high_af_alpha = model.cluster_alpha[1];
    model.high_af_beta = model.cluster_beta[1];
    model.background_weight = std::exp(model.cluster_log_weights[0]);
    model.high_af_weight = std::exp(model.cluster_log_weights[1]);
}

SomaticModelParameters learn_empirical_somatic_model(
    const std::vector<SomaticModelDatum>& data, const double callable_sites) {
    SomaticModelParameters model;
    initialize_somatic_model_clusters(model);
    if (data.size() < 2) return model;
    model.learned = true;
    model.callable_sites = callable_sites;

    // GATK initializes additional binomial clusters from probability-weighted
    // AF peaks, accepting a split only when it improves the BIC score.
    double previous_bic = -std::numeric_limits<double>::infinity();
    for (int split = 0; split < 5 && model.cluster_count < SomaticModelParameters::kMaxClusters; ++split) {
        std::vector<std::pair<double, double>> af_and_somatic;
        af_and_somatic.reserve(data.size());
        double total_somatic = 0.0;
        for (const auto& datum : data) {
            if (datum.total_count <= 0 || datum.alt_count < 0 || datum.alt_count > datum.total_count) continue;
            const auto weight = std::clamp(datum.somatic_weight, 0.0, 1.0);
            af_and_somatic.emplace_back(
                static_cast<double>(datum.alt_count) / static_cast<double>(datum.total_count), weight);
            total_somatic += weight;
        }
        if (af_and_somatic.empty() || !(total_somatic > 0.0)) break;
        std::sort(af_and_somatic.begin(), af_and_somatic.end());
        std::vector<double> quantiles;
        const auto step_size = total_somatic / 50.0;
        double cumulative = 0.0;
        double next_quantile = step_size;
        for (const auto& item : af_and_somatic) {
            cumulative += item.second;
            while (step_size > 0.0 && cumulative > next_quantile) {
                if (quantiles.empty() || quantiles.back() != item.first) quantiles.push_back(item.first);
                next_quantile += step_size;
            }
        }
        if (quantiles.empty()) break;
        std::vector<double> responsibilities(quantiles.size(), 0.0);
        for (const auto& datum : data) {
            if (datum.total_count <= 0 || datum.alt_count < 0 || datum.alt_count > datum.total_count) continue;
            const auto background_log = model.cluster_log_weights[0] +
                host_cluster_count_log_likelihood(datum.total_count, datum.alt_count, 0, model);
            double maximum = -std::numeric_limits<double>::infinity();
            for (int cluster = 0; cluster < model.cluster_count; ++cluster)
                maximum = std::max(maximum, model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
                    host_cluster_count_log_likelihood(datum.total_count, datum.alt_count, cluster, model));
            double mixture = 0.0;
            for (int cluster = 0; cluster < model.cluster_count; ++cluster)
                mixture += std::exp(model.cluster_log_weights[static_cast<std::size_t>(cluster)] +
                                    host_cluster_count_log_likelihood(datum.total_count, datum.alt_count, cluster, model) - maximum);
            const auto background_probability = (mixture > 0.0 && std::isfinite(mixture))
                ? std::exp(background_log - maximum) / mixture : 0.0;
            const auto somatic_weight = std::clamp(datum.somatic_weight, 0.0, 1.0);
            for (std::size_t q = 0; q < quantiles.size(); ++q) {
                const auto f = quantiles[q];
                const auto log_mass = kokkos_log_binomial_mass(
                    datum.total_count, datum.alt_count, std::clamp(f, 0.0, 1.0));
                const auto density = log_mass < -700.0 ? 0.0 : std::exp(log_mass);
                responsibilities[q] += somatic_weight * background_probability * density * (datum.total_count + 1.0);
            }
        }
        std::vector<std::pair<double, double>> peaks;
        double peak_mass = 0.0;
        double peak_af = 0.0;
        double peak_height = 0.0;
        for (std::size_t q = 0; q < quantiles.size(); ++q) {
            const auto left_height = q == 0 ? 0.0 : responsibilities[q - 1];
            const auto right_height = q + 1 == quantiles.size() ? 0.0 : responsibilities[q + 1];
            const auto left_af = q == 0 ? 0.0 : quantiles[q - 1];
            peak_mass += (quantiles[q] - left_af) * (left_height + responsibilities[q]) / 2.0;
            if (responsibilities[q] > peak_height) {
                peak_height = responsibilities[q];
                peak_af = quantiles[q];
            }
            const bool local_min = (responsibilities[q] < left_height && responsibilities[q] <= right_height) ||
                (responsibilities[q] <= left_height && responsibilities[q] < right_height);
            if ((local_min && q > 0) || q + 1 == quantiles.size()) {
                peaks.emplace_back(peak_af, peak_mass);
                peak_mass = 0.0;
                peak_af = quantiles[q];
                peak_height = responsibilities[q];
            }
        }
        if (peaks.empty()) break;
        const auto biggest = std::max_element(peaks.begin(), peaks.end(),
            [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        const auto minimum_index = std::min<std::size_t>(5, quantiles.size() - 1);
        if (biggest->first < quantiles[minimum_index]) break;
        double total_mass = 0.0;
        for (const auto& peak : peaks) total_mass += peak.second;
        if (!(total_mass > 0.0)) break;
        const auto fraction = std::min(0.9, biggest->second / total_mass);
        const auto old_model = model;
        const auto new_index = model.cluster_count++;
        model.cluster_log_weights[0] = std::log1p(fraction) + old_model.cluster_log_weights[0];
        model.cluster_log_weights[static_cast<std::size_t>(new_index)] =
            std::log(std::max(fraction, 1.0e-12)) + old_model.cluster_log_weights[0];
        model.cluster_alpha[static_cast<std::size_t>(new_index)] =
            std::max(biggest->first * 10000.0, 0.1);
        model.cluster_beta[static_cast<std::size_t>(new_index)] =
            std::max((1.0 - biggest->first) * 10000.0, 0.1);
        model.cluster_mean[static_cast<std::size_t>(new_index)] = biggest->first;
        model.cluster_binomial[static_cast<std::size_t>(new_index)] = 1;
        for (int iteration = 0; iteration < 5; ++iteration)
            perform_somatic_model_em_iteration(model, data);
        double weighted_log_likelihood = 0.0;
        for (const auto& datum : data)
            weighted_log_likelihood += std::clamp(datum.somatic_weight, 0.0, 1.0) *
                host_model_count_log_likelihood(datum.total_count, datum.alt_count, model);
        const auto effective_count = std::max(total_somatic, 1.0);
        const auto current_bic = weighted_log_likelihood -
            2.0 * static_cast<double>(model.cluster_count) * std::log(effective_count);
        if (current_bic < previous_bic) {
            model = old_model;
            break;
        }
        previous_bic = current_bic;
    }
    for (int iteration = 0; iteration < 5; ++iteration)
        perform_somatic_model_em_iteration(model, data);
    double variant_count = 0.0;
    for (const auto& datum : data) variant_count += std::clamp(datum.somatic_weight, 0.0, 1.0);
    const auto artifact_count = std::max(0.0, static_cast<double>(data.size()) - variant_count);
    model.log_variant_vs_artifact_prior = std::log(
        (variant_count + 1.0) / (variant_count + artifact_count + 2.0));
    for (auto& prior : model.log_variant_priors) prior = model.log_indel_prior;
    model.log_variant_priors[10] = model.log_snv_prior;
    if (std::isfinite(callable_sites) && callable_sites >= 1.0) {
        std::array<double, 21> variant_counts{};
        for (const auto& datum : data) {
            if (datum.indel_length < -10 || datum.indel_length > 10) continue;
            variant_counts[static_cast<std::size_t>(datum.indel_length + 10)] +=
                std::clamp(datum.somatic_weight, 0.0, 1.0);
        }
        for (int indel_length = -10; indel_length <= 10; ++indel_length) {
            const auto count = variant_counts[static_cast<std::size_t>(indel_length + 10)];
            const auto minimum = indel_length == 0 ? 1.0e-8 : 1.0e-9;
            model.log_variant_priors[static_cast<std::size_t>(indel_length + 10)] =
                std::log(std::max(count / callable_sites, minimum));
        }
        model.log_snv_prior = model.log_variant_priors[10];
        model.log_indel_prior = model.log_variant_priors[9];
    }
    return model;
}

double slippage_probability_kokkos(
    const double slippage_rate, const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length, FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_slippage",
        kokkos_slippage_probability(
            slippage_rate, reference_repeats, alternate_repeats, repeat_unit_length,
            min_slippage_length, total_count, alternate_count, reference_length,
            alternate_length),
        telemetry);
}

double slippage_probability_with_prior_kokkos(
    const double slippage_rate, const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length, const double log_somatic_prior,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_slippage_prior",
        kokkos_slippage_probability_with_prior(
            slippage_rate, reference_repeats, alternate_repeats, repeat_unit_length,
            min_slippage_length, total_count, alternate_count, reference_length,
            alternate_length, log_somatic_prior),
        telemetry);
}

double empirical_slippage_probability_kokkos(
    const SomaticModelParameters& model, const double slippage_rate,
    const int reference_repeats, const int alternate_repeats,
    const int repeat_unit_length, const int min_slippage_length,
    const int total_count, const int alternate_count, const int reference_length,
    const int alternate_length, const double log_somatic_prior,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_empirical_slippage",
        kokkos_empirical_slippage_probability(
            model, slippage_rate, reference_repeats, alternate_repeats,
            repeat_unit_length, min_slippage_length, total_count, alternate_count,
            reference_length, alternate_length, log_somatic_prior),
        telemetry);
}

double combined_error_probability_for_alt_kokkos(
    const bcf_hdr_t* header, bcf1_t* record, const int alt_index,
    const std::vector<int>& tumor_depths, const std::vector<float>& tlod_values,
    const SomaticModelParameters* empirical_model = nullptr,
    FilterKernelTelemetry* telemetry = nullptr,
    const std::map<int, double>* contamination_by_alt = nullptr,
    const double orientation_posterior = std::numeric_limits<double>::quiet_NaN()) {
    if (header == nullptr || record == nullptr || alt_index <= 0 ||
        alt_index >= record->n_allele) return std::numeric_limits<double>::quiet_NaN();
    double sequencing_error = info_probability_at(header, record, "PSOMATIC", alt_index);
    if (std::isfinite(sequencing_error)) sequencing_error = 1.0 - sequencing_error;
    if (!std::isfinite(sequencing_error) &&
        alt_index <= static_cast<int>(tlod_values.size()) &&
        tumor_depths.size() == static_cast<std::size_t>(record->n_allele)) {
        int total_depth = 0;
        for (const auto depth : tumor_depths) if (depth > 0) total_depth += depth;
        if (total_depth > 0 && tumor_depths[static_cast<std::size_t>(alt_index)] >= 0) {
            const auto log_odds = static_cast<double>(
                tlod_values[static_cast<std::size_t>(alt_index - 1)]) * std::log(10.0);
            const auto reference_length = record->d.allele[0] == nullptr
                ? 0 : static_cast<int>(std::strlen(record->d.allele[0]));
            const auto alternate_length = record->d.allele[alt_index] == nullptr
                ? reference_length : static_cast<int>(std::strlen(record->d.allele[alt_index]));
            const auto log_snv_prior = empirical_model != nullptr
                ? empirical_model->log_snv_prior : -6.0 * std::log(10.0);
            const auto log_indel_prior = empirical_model != nullptr
                ? empirical_model->log_indel_prior : -7.0 * std::log(10.0);
            const auto indel_length = alternate_length - reference_length;
            double log_somatic_prior = alternate_length == reference_length
                ? log_snv_prior - std::log(3.0) : log_indel_prior;
            if (empirical_model != nullptr && indel_length >= -10 && indel_length <= 10 &&
                std::isfinite(empirical_model->log_variant_priors[
                    static_cast<std::size_t>(indel_length + 10)]))
                log_somatic_prior = empirical_model->log_variant_priors[
                    static_cast<std::size_t>(indel_length + 10)] +
                    (indel_length == 0 ? -std::log(3.0) : 0.0);
            sequencing_error = empirical_model != nullptr && empirical_model->learned
                ? empirical_tumor_evidence_probability_kokkos(
                    log_odds, total_depth, tumor_depths[static_cast<std::size_t>(alt_index)],
                    log_somatic_prior, *empirical_model, telemetry)
                : tumor_evidence_probability_with_prior_kokkos(
                    log_odds, total_depth, tumor_depths[static_cast<std::size_t>(alt_index)],
                    log_somatic_prior, telemetry);
        }
    }
    const auto germline = info_probability_at(header, record, "PGERMLINE", alt_index);
    // ErrorProbabilities groups GermlineFilter and ContaminationFilter under
    // NON_SOMATIC and takes the maximum probability within that type.  The
    // contamination map is optional because legacy records may not have the
    // POPAF/AD evidence required by ContaminationFilter.
    double non_somatic = germline;
    if (contamination_by_alt != nullptr) {
        const auto iterator = contamination_by_alt->find(alt_index);
        if (iterator != contamination_by_alt->end() && std::isfinite(iterator->second))
            non_somatic = std::isfinite(non_somatic)
                ? std::max(non_somatic, iterator->second) : iterator->second;
    }
    const auto artifact = info_probability_at(header, record, "PARTIFACT", alt_index);
    const auto orientation_annotation = info_probability_at(header, record, "OBP", alt_index);
    // ReadOrientationFilter is a Mutect2VariantFilter in GATK, so its
    // weighted-median posterior is a site-level value copied to every ALT.
    // Prefer the model result when orientation priors are supplied, while
    // retaining OBP as a compatibility fallback for pre-annotated VCFs.
    const auto orientation = std::isfinite(orientation_posterior)
        ? orientation_posterior : orientation_annotation;
    // ErrorProbabilities reduces all technical-artifact filters (including
    // PARTIFACT and orientation-bias) with max(), then combines that ARTIFACT
    // type independently with NON_SOMATIC and SEQUENCING.
    const auto artifact_type = std::isfinite(artifact) && std::isfinite(orientation)
        ? std::max(artifact, orientation)
        : (std::isfinite(artifact) ? artifact : orientation);
    return combine_error_probabilities_kokkos(
        sequencing_error, non_somatic, artifact_type, -1.0, telemetry);
}

KOKKOS_INLINE_FUNCTION double kokkos_normal_artifact_probability(
    const double normal_artifact_log10_odds,
    const int tumor_depth, const int tumor_alt_depth,
    const int normal_depth, const int normal_alt_depth,
    const int median_reference_base_quality,
    const double normal_pileup_p_value_threshold) {
    // This is the deterministic, pre-EM NormalArtifactFilter boundary from
    // GATK/Rust.  A clean normal whose AF is <10% of the tumour AF cannot be
    // called an artifact in normal; otherwise NALOD and the normal pileup
    // binomial tail determine the error probability.
    const double tumor_af = tumor_depth == 0
        ? 0.0 : static_cast<double>(tumor_alt_depth) / static_cast<double>(tumor_depth);
    const double normal_af = normal_depth == 0
        ? 0.0
        : static_cast<double>(normal_alt_depth) / static_cast<double>(normal_depth);
    if (tumor_depth != 0 && normal_af < 0.1 * tumor_af) return 0.0;
    const double error_probability = Kokkos::pow(
        10.0, -static_cast<double>(median_reference_base_quality) / 10.0);
    double tail = 0.0;
    for (int successes = normal_alt_depth; successes <= normal_depth; ++successes) {
        const auto log_mass = kokkos_log_binomial_mass(normal_depth, successes, error_probability);
        if (log_mass > -1.0e299) tail += Kokkos::exp(log_mass);
    }
    if (tail < normal_pileup_p_value_threshold) return 1.0;
    const double log_prior_real = -Kokkos::log(10.0);
    const double log_real = -normal_artifact_log10_odds * Kokkos::log(10.0) + log_prior_real;
    const double log_error_prior = Kokkos::log1p(-Kokkos::exp(log_prior_real));
    const double maximum = log_real > log_error_prior ? log_real : log_error_prior;
    const double normalizer = Kokkos::exp(log_real - maximum) +
                              Kokkos::exp(log_error_prior - maximum);
    return Kokkos::exp(log_error_prior - maximum) / normalizer;
}

KOKKOS_INLINE_FUNCTION double kokkos_normal_artifact_probability_with_prior(
    const double normal_artifact_log10_odds,
    const int tumor_depth, const int tumor_alt_depth,
    const int normal_depth, const int normal_alt_depth,
    const int median_reference_base_quality,
    const double normal_pileup_p_value_threshold,
    const double log_prior_real) {
    const double tumor_af = tumor_depth == 0
        ? 0.0 : static_cast<double>(tumor_alt_depth) / static_cast<double>(tumor_depth);
    const double normal_af = normal_depth == 0
        ? 0.0 : static_cast<double>(normal_alt_depth) / static_cast<double>(normal_depth);
    if (tumor_depth != 0 && normal_af < 0.1 * tumor_af) return 0.0;
    const double error_probability = Kokkos::pow(
        10.0, -static_cast<double>(median_reference_base_quality) / 10.0);
    double tail = 0.0;
    for (int successes = normal_alt_depth; successes <= normal_depth; ++successes) {
        const auto log_mass = kokkos_log_binomial_mass(normal_depth, successes, error_probability);
        if (log_mass > -1.0e299) tail += Kokkos::exp(log_mass);
    }
    if (tail < normal_pileup_p_value_threshold) return 1.0;
    const double log_real = -normal_artifact_log10_odds * Kokkos::log(10.0) + log_prior_real;
    const double log_error_prior = Kokkos::log1p(-Kokkos::exp(log_prior_real));
    const double maximum = log_real > log_error_prior ? log_real : log_error_prior;
    const double normalizer = Kokkos::exp(log_real - maximum) +
                              Kokkos::exp(log_error_prior - maximum);
    return Kokkos::exp(log_error_prior - maximum) / normalizer;
}

double normal_artifact_probability_kokkos(
    const double normal_artifact_log10_odds,
    const int tumor_depth, const int tumor_alt_depth,
    const int normal_depth, const int normal_alt_depth,
    const int median_reference_base_quality,
    const double normal_pileup_p_value_threshold,
    FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_normal_artifact",
        kokkos_normal_artifact_probability(
            normal_artifact_log10_odds, tumor_depth, tumor_alt_depth,
            normal_depth, normal_alt_depth, median_reference_base_quality,
            normal_pileup_p_value_threshold),
        telemetry);
}

double normal_artifact_probability_with_prior_kokkos(
    const double normal_artifact_log10_odds,
    const int tumor_depth, const int tumor_alt_depth,
    const int normal_depth, const int normal_alt_depth,
    const int median_reference_base_quality,
    const double normal_pileup_p_value_threshold,
    const double log_prior_real, FilterKernelTelemetry* telemetry = nullptr) {
    FASTGATK_FILTER_SCALAR_KERNEL(
        "filter_mutect_normal_artifact_prior",
        kokkos_normal_artifact_probability_with_prior(
            normal_artifact_log10_odds, tumor_depth, tumor_alt_depth,
            normal_depth, normal_alt_depth, median_reference_base_quality,
            normal_pileup_p_value_threshold, log_prior_real),
        telemetry);
}

double orientation_artifact_probability_kokkos(
    const PriorVector& prior, char ref_base, char alt_base,
    int alt_depth, int alt_f1r2, int alt_f2r1, int depth,
    Kokkos::View<double*>& prior_device, Kokkos::View<double*>& log_device,
    FilterKernelTelemetry* telemetry = nullptr) {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    int ref_index = -1;
    int alt_index = -1;
    for (int index = 0; index < 4; ++index) {
        if (bases[index] == ref_base) ref_index = index;
        if (bases[index] == alt_base) alt_index = index;
    }
    // AD is the allele-depth denominator for the allele-fraction component,
    // while F1R2/F2R1 are the independent orientation observations.  GATK's
    // ReadOrientationModel does not infer one from the other: a caller may
    // drop a read from AD while retaining an orientation count (or vice
    // versa).  Keep the two trial counts separate and fail closed only when
    // either observed vector is internally invalid.
    if (ref_index < 0 || alt_index < 0 || ref_index == alt_index || depth <= 0 ||
        alt_depth < 0 || alt_depth > depth || alt_f1r2 < 0 || alt_f2r1 < 0) return 0.0;
    const auto orientation_depth = alt_f1r2 + alt_f2r1;
    if (orientation_depth <= 0) return 0.0;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("filter-mutect-orientation-v1");
    host_batch.records = 12;
    host_batch.bytes = 24 * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("filter-mutect-orientation");
    kernel_plan.begin_prepare(host_batch);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(12);
    device_batch.bind("prior", prior_device);
    device_batch.bind("log_posterior", log_device);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    auto prior_host = Kokkos::create_mirror_view(prior_device);
    for (int state = 0; state < 12; ++state) prior_host(state) = prior[static_cast<std::size_t>(state)];
    Kokkos::deep_copy(prior_device, prior_host);
    kernel_plan.begin_execute();
    Kokkos::parallel_for("filter_mutect_orientation_posterior", Kokkos::RangePolicy<ExecSpace>(0, 12),
        KOKKOS_LAMBDA(const int state) {
            // FilterMutectCalls calls the responsibility routine with
            // givenNotHomRef=true: HOM_REF is not a candidate once an ALT
            // genotype is being evaluated.
            if (state == 8 || (state < 8 && ((state % 4) == ref_index || (state % 4) != alt_index))) {
                log_device(state) = -1.0e300;
                return;
            }
            const double alpha_af = state < 8 ? 1.0 : state == 8 ? 3.0 : state == 9 ? 5.0 : state == 10 ? 2.0 : 10000.0;
            const double beta_af = state < 8 ? 9.0 : state == 8 ? 10000.0 : state == 9 ? 5.0 : state == 10 ? 5.0 : 3.0;
            const double alpha_orientation = state < 4 ? 100.0 : state < 8 ? 1.0 : 10.0;
            const double beta_orientation = state < 4 ? 1.0 : state < 8 ? 100.0 : 10.0;
            const int orientation_successes = state < 4 ? alt_f1r2 : alt_f2r1;
            const auto prior_value = Kokkos::fmax(prior_device(state), 1.0e-300);
            log_device(state) = Kokkos::log(prior_value) +
                kokkos_log_beta_binomial_probability(alt_depth, depth, alpha_af, beta_af) +
                kokkos_log_beta_binomial_probability(orientation_successes, orientation_depth,
                                                     alpha_orientation, beta_orientation);
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) telemetry->record(kernel_plan.telemetry(), 12);
    auto log_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, log_device);
    double maximum = -std::numeric_limits<double>::infinity();
    for (int state = 0; state < 12; ++state) maximum = std::max(maximum, log_host(state));
    if (!std::isfinite(maximum)) return 0.0;
    std::array<double, 12> posterior{};
    double normalizer = 0.0;
    for (int state = 0; state < 12; ++state) {
        posterior[static_cast<std::size_t>(state)] = log_host(state) < -1.0e299
            ? 0.0 : std::exp(log_host(state) - maximum);
        normalizer += posterior[static_cast<std::size_t>(state)];
    }
    if (!(normalizer > 0.0)) return 0.0;
    return std::max(posterior[static_cast<std::size_t>(alt_index)],
                    posterior[static_cast<std::size_t>(4 + alt_index)]) / normalizer;
}

#undef FASTGATK_FILTER_SCALAR_KERNEL

const PriorVector* find_orientation_prior(
    const std::vector<OrientationPriorCollection>& collections,
    const std::string& sample, const std::string& context) {
    for (const auto& collection : collections) {
        if (collection.sample != sample) continue;
        const auto iterator = collection.contexts.find(context);
        if (iterator != collection.contexts.end()) return &iterator->second;
    }
    return nullptr;
}

int max_tlod_alt_index(const bcf_hdr_t* header, bcf1_t* record) {
    const auto values = info_float_values(header, record, "TLOD");
    int best = 1;
    float maximum = -std::numeric_limits<float>::infinity();
    for (int index = 0; index < static_cast<int>(values.size()) && index + 1 < record->n_allele; ++index) {
        if (values[static_cast<std::size_t>(index)] > maximum) {
            maximum = values[static_cast<std::size_t>(index)];
            best = index + 1;
        }
    }
    return best;
}

struct OrientationEvaluation {
    double probability = 0.0;
    bool has_prior = false;
    std::uint64_t sample_evaluations = 0;
};

std::string join_strings(const std::vector<std::string>& values, const char* separator) {
    std::ostringstream result;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) result << separator;
        result << values[index];
    }
    return result.str();
}

int indel_length_for_alt(const bcf1_t* record, int allele_index) {
    if (record == nullptr || allele_index <= 0 || allele_index >= record->n_allele ||
        record->d.allele[0] == nullptr || record->d.allele[allele_index] == nullptr)
        return 0;
    const auto* ref = record->d.allele[0];
    const auto* alt = record->d.allele[allele_index];
    if (ref[0] == '<' || alt[0] == '<') return 0;
    return static_cast<int>(std::strlen(alt)) - static_cast<int>(std::strlen(ref));
}

OrientationEvaluation evaluate_orientation_posterior(
    const bcf_hdr_t* header, bcf1_t* record, faidx_t* reference,
    const std::vector<OrientationPriorCollection>& collections,
    const std::vector<int>& sample_indices,
    Kokkos::View<double*>& prior_device, Kokkos::View<double*>& log_device,
    FilterKernelTelemetry* telemetry = nullptr) {
    OrientationEvaluation result;
    const int alt_index = max_tlod_alt_index(header, record);
    if (collections.empty() || record == nullptr || record->n_allele <= alt_index ||
        record->d.allele[0] == nullptr || record->d.allele[alt_index] == nullptr ||
        std::strlen(record->d.allele[0]) == 0 || std::strlen(record->d.allele[alt_index]) == 0 ||
        std::strlen(record->d.allele[0]) != std::strlen(record->d.allele[alt_index]))
        return result;
    const auto* ref_allele = record->d.allele[0];
    const auto* alt_allele = record->d.allele[alt_index];
    const auto allele_length = std::strlen(ref_allele);
    std::vector<std::pair<int, double>> evidence;
    evidence.reserve(sample_indices.size());
    for (const auto sample_index : sample_indices) {
        const auto counts = read_orientation_counts(header, record, sample_index, alt_index);
        if (counts.ref < 0 || counts.alt < 0 || counts.f1r2 < 0 || counts.f2r1 < 0) continue;
        const auto alt_weight = format_sum_alt_depth(header, record, sample_index);
        if (alt_weight < 0) continue;
        double sample_probability = 0.0;
        bool sample_has_prior = false;
        if (!sample_is_hom_ref(header, record, sample_index)) {
            const auto* sample_name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample_index);
            // ReadOrientationFilter derives both AD and the model depth from
            // the F1R2/F2R1 vectors: REF is f1r2[0]+f2r1[0], ALT is the
            // selected allele in each vector.  Do not substitute FORMAT/AD
            // when these standard Number=R annotations are present; AD can
            // legitimately differ after filtering or duplicate removal.
            const int orientation_alt_count = counts.f1r2 + counts.f2r1;
            const int orientation_ref_count = counts.ref_f1r2 + counts.ref_f2r1;
            const int depth = orientation_ref_count + orientation_alt_count;
            for (std::size_t offset = 0; offset < allele_length; ++offset) {
                const auto context = reference_context(reference, header, record, static_cast<int>(offset));
                if (context.empty()) continue;
                const auto* prior = sample_name == nullptr ? nullptr : find_orientation_prior(
                    collections, sample_name, context);
                if (prior == nullptr) continue;
                const auto ref_base = static_cast<char>(std::toupper(
                    static_cast<unsigned char>(ref_allele[offset])));
                const auto alt_base = static_cast<char>(std::toupper(
                    static_cast<unsigned char>(alt_allele[offset])));
                sample_probability = std::max(sample_probability, orientation_artifact_probability_kokkos(
                    *prior, ref_base, alt_base, orientation_alt_count,
                    counts.f1r2, counts.f2r1, depth,
                    prior_device, log_device, telemetry));
                sample_has_prior = true;
            }
        }
        if (sample_has_prior) {
            result.has_prior = true;
            ++result.sample_evaluations;
        }
        evidence.emplace_back(alt_weight, sample_probability);
    }
    if (result.has_prior) result.probability = weighted_median_posterior(std::move(evidence));
    return result;
}


namespace fmc_j = fastgatk_fmc_java;

// Java-parity learning pipeline for FilterMutectCalls: three accumulate
// passes over the input VCF; the first two learn the
// SomaticClusteringModel transcription (fmc_java_model.hpp) from Java
// record() data; the third supplies pass-2 errors for the learned threshold.
// Results are transferred into the native SomaticModelParameters so the
// existing final filtering loop consumes the Java-learned model.
static void run_java_learning_pipeline(
    const Options& options, SomaticModelParameters& native_model,
    double& effective_threshold, std::size_t& observations_out,
    std::vector<double>* threshold_errors_out,
    const std::string& mutect_stats_path,
    faidx_t* reference_index,
    const std::vector<OrientationPriorCollection>& orientation_collections,
    const std::vector<int>& orientation_sample_indices,
    const std::vector<int>& tumor_sample_indices,
    const std::vector<int>& haplotype_sample_indices,
    LearnedHaplotypeFilter* haplotype_filter,
    LearnedStrandArtifactFilter& strand_artifact_filter,
    Kokkos::View<double*>& orientation_prior_device,
    Kokkos::View<double*>& orientation_log_device,
    const std::unordered_map<std::string, float>& contamination_by_sample,
    const float effective_contamination,
    FilterKernelTelemetry* telemetry,
    std::uint64_t* empirical_records_out) {
    fmc_j::SomaticModel model;
    model.clusters.push_back(fmc_j::Cluster{0, 1.0, 1.0, 0.5});
    model.clusters.push_back(fmc_j::Cluster{0, 10.0, 1.0, 10.0 / 11.0});
    model.log_cluster_weights = {std::log1p(0.01), std::log(0.01)};
    for (int i = 0; i < 21; ++i)
        model.log_variant_priors[static_cast<std::size_t>(i)] = -7.0 * fmc_j::kLn10;
    model.log_variant_priors[10] = -6.0 * fmc_j::kLn10;
    // Callable sites come from the Mutect stats summary when available.
    const auto stats = read_mutect_stats_summary(mutect_stats_path);
    model.callable_sites = stats.callable_sites;
    model.has_callable_sites = stats.callable_sites >= 1.0;
    const double contamination =
        contamination_by_sample.empty()
            ? static_cast<double>(effective_contamination)
            : static_cast<double>(contamination_by_sample.begin()->second);

    auto accumulate_pass = [&](const bool learn_after,
                                std::vector<double>* errors_ptr = nullptr) {
        model.data.clear();
        // Phase 1: stream the VCF collecting the guarded per-allele evidence
        // (sequencing/germline/contamination posteriors are a pure function of
        // the current model + per-datum features, so batching preserves the
        // exact values and Java record order while the numerics run in the
        // shared kernels module).
        std::vector<int> ev_alt;
        std::vector<int> ev_total;
        std::vector<double> ev_tlo;
        std::vector<double> ev_pop_af;
        std::vector<double> ev_nlod;
        std::vector<double> ev_art;
        // Keep the artifact sources separate until ErrorProbabilities has
        // supplied the germline posterior.  FilteredHaplotypeFilter excludes
        // itself and conditionally excludes NormalArtifactFilter when the
        // record is likely germline; reducing them early loses that source
        // ordering.
        std::vector<double> ev_non_normal_non_haplotype_art;
        std::vector<double> ev_normal_art;
        std::vector<int> ev_indel_length;
        std::vector<std::uint8_t> ev_germ_gate;
        struct HaplotypePassRecord {
            std::size_t evidence_begin = 0;
            std::size_t evidence_end = 0;
            std::int64_t locus = 0;
            std::vector<PhasedHaplotypeObservation> observations;
        };
        std::vector<HaplotypePassRecord> haplotype_pass_records;
        {
            htsFile* input = bcf_open(options.input.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen Mutect2 VCF for Java learning");
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) throw std::runtime_error("BAD_INPUT: cannot reread Mutect2 VCF header");
            const auto normal_names = normal_sample_names(header);
            std::vector<int> normal_sample_indices;
            for (int sample_index = 0; sample_index < bcf_hdr_nsamples(header); ++sample_index) {
                const auto* sample_name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample_index);
                if (sample_name != nullptr && normal_names.count(sample_name) != 0)
                    normal_sample_indices.push_back(sample_index);
            }
            bcf1_t* record = bcf_init();
            while (fastgatk::io::read_variant_record(input, header, record,
                                                   options.input) == 0) {
                bcf_unpack(record, BCF_UN_ALL);
                const auto evidence_begin = ev_alt.size();
                const auto phased_observations = haplotype_filter != nullptr
                    ? phased_haplotype_observations(
                        header, record, haplotype_sample_indices)
                    : std::vector<PhasedHaplotypeObservation>{};
                double learned_haplotype_artifact = 0.0;
                if (haplotype_filter != nullptr && !phased_observations.empty()) {
                    const auto best_haplotype = std::max_element(
                        phased_observations.begin(), phased_observations.end(),
                        [](const auto& left, const auto& right) {
                            // Java Stream.max keeps the first genotype on an
                            // AF tie, which is std::max_element's behavior.
                            return left.allele_fraction < right.allele_fraction;
                        });
                    learned_haplotype_artifact = haplotype_filter->error_probability(
                        best_haplotype->key, static_cast<std::int64_t>(record->pos) + 1);
                }
                const auto tlod_values = info_float_values(header, record, "TLOD");
                const auto popaf_values = info_float_values(header, record, "POPAF");
                const auto nlod_values = info_float_values(header, record, "NLOD");
                const auto nalod_values = info_float_values(header, record, "NALOD");
                const auto mbq_values = info_int_values(header, record, "MBQ");
                // ReadOrientationFilter is a variant-level ARTIFACT contributor
                // to Java ErrorProbabilities: evaluate its posterior for the
                // max-TLOD ALT exactly like the final filtering loop, and let
                // that ALT carry it into the Java record()/threshold errors.
                OrientationEvaluation orientation_evaluation;
                if (!orientation_collections.empty()) {
                    orientation_evaluation = evaluate_orientation_posterior(
                        header, record, reference_index, orientation_collections,
                        orientation_sample_indices, orientation_prior_device,
                        orientation_log_device);
                }
                const int orientation_alt = max_tlod_alt_index(header, record);
                // NormalArtifactFilter is a variant-level ARTIFACT source in
                // ErrorProbabilities.  It must participate in the two model
                // learning passes (not just final VCF writing), otherwise the
                // learned F-score threshold incorrectly admits low-evidence
                // calls that GATK filters.
                double normal_artifact_probability = 0.0;
                if (!normal_sample_indices.empty() && !tlod_values.empty() &&
                    orientation_alt > 0 &&
                    orientation_alt <= static_cast<int>(nalod_values.size())) {
                    const auto tumor_depths = tumor_allele_depths(
                        header, record, tumor_sample_indices);
                    const auto normal_depths = tumor_allele_depths(
                        header, record, normal_sample_indices);
                    if (tumor_depths.size() == static_cast<std::size_t>(record->n_allele) &&
                        normal_depths.size() == static_cast<std::size_t>(record->n_allele)) {
                        int tumor_depth = 0;
                        int normal_depth = 0;
                        for (const int depth : tumor_depths) if (depth > 0) tumor_depth += depth;
                        for (const int depth : normal_depths) if (depth > 0) normal_depth += depth;
                        const int tumor_alt = tumor_depths[static_cast<std::size_t>(orientation_alt)];
                        const int normal_alt = normal_depths[static_cast<std::size_t>(orientation_alt)];
                        if (tumor_depth > 0 && normal_depth >= 0 && tumor_alt >= 0 && normal_alt >= 0) {
                            normal_artifact_probability = normal_artifact_probability_with_prior_kokkos(
                                nalod_values[static_cast<std::size_t>(orientation_alt - 1)],
                                tumor_depth, tumor_alt, normal_depth, normal_alt,
                                mbq_values.empty() ? 30 : mbq_values.front(),
                                options.normal_pileup_p_value_threshold,
                                model.log_variant_vs_artifact_prior, telemetry);
                        }
                    }
                }
                const int alts = record->n_allele;
                // ErrorProbabilities reduces every ARTIFACT filter with a
                // per-ALT max before SomaticClusteringModel.record().  The
                // learning path used to carry only NormalArtifact and
                // ReadOrientation, leaving deterministic hard artifacts (in
                // particular clustered calls) in the clustering data.  That
                // changes the learned priors on dense real callsets.
                const auto strand_counts_by_allele = info_allele_strand_counts(
                    header, record, "AS_SB_TABLE");
                int strand_total_forward = 0;
                int strand_total_reverse = 0;
                if (strand_counts_by_allele &&
                    strand_counts_by_allele->size() == static_cast<std::size_t>(alts)) {
                    for (const auto& counts : *strand_counts_by_allele) {
                        strand_total_forward += counts.first;
                        strand_total_reverse += counts.second;
                    }
                }
                std::map<int, double> strand_artifact_by_alt;
                if (strand_counts_by_allele &&
                    strand_counts_by_allele->size() == static_cast<std::size_t>(alts)) {
                    for (int alt_index = 1; alt_index < alts; ++alt_index) {
                        if (record->d.allele[alt_index] == nullptr ||
                            record->d.allele[alt_index][0] == '<') continue;
                        const auto& counts = (*strand_counts_by_allele)[
                            static_cast<std::size_t>(alt_index)];
                        const auto e_step = strand_artifact_filter.calculate_e_step(
                            strand_total_forward, strand_total_reverse,
                            counts.first, counts.second,
                            std::abs(indel_length_for_alt(record, alt_index)), telemetry);
                        strand_artifact_filter.accumulate(e_step);
                        strand_artifact_by_alt.emplace(alt_index,
                            e_step.forward_artifact_responsibility +
                            e_step.reverse_artifact_responsibility);
                    }
                }
                const auto mmq_values = info_int_values(header, record, "MMQ");
                const auto mfrl_values = info_int_values(header, record, "MFRL");
                const auto n_count = info_first_int_optional(header, record, "NCount");
                const auto region_event_count = info_int(header, record, "ECNT", -1);
                const auto haplotype_event_count = info_int_max(header, record, "ECNTH");
                const bool panel_of_normals = info_present(header, record, "PON");
                int concrete_alt_count = 0;
                int high_lod_alt_count = 0;
                for (int candidate = 1; candidate < alts; ++candidate) {
                    const auto* allele = record->d.allele[candidate];
                    if (allele == nullptr || allele[0] == '<') continue;
                    ++concrete_alt_count;
                    if (candidate - 1 < static_cast<int>(tlod_values.size()) &&
                        static_cast<double>(tlod_values[static_cast<std::size_t>(candidate - 1)]) *
                            fmc_j::kLn10 > 5.0)
                        ++high_lod_alt_count;
                }
                const bool clustered_artifact = !options.mitochondria_mode &&
                    !options.microbial_mode &&
                    ((region_event_count >= 0 && region_event_count > options.max_events_in_region) ||
                     (haplotype_event_count && *haplotype_event_count > options.max_events_in_haplotype));
                const bool multiallelic_artifact = !options.mitochondria_mode &&
                    !options.microbial_mode && high_lod_alt_count > options.max_alt_allele_count;
                const bool fragment_artifact = !options.mitochondria_mode &&
                    !options.microbial_mode && options.max_median_fragment_length_difference >= 0 &&
                    mfrl_values.size() >= 2 &&
                    std::abs(mfrl_values[1] - mfrl_values[0]) >
                        options.max_median_fragment_length_difference;
                for (int alt_index = 1; alt_index < alts; ++alt_index) {
                    if (record->d.allele[alt_index] == nullptr ||
                        record->d.allele[alt_index][0] == '<')
                        continue;
                    // Mutect2's VCF does not require the tumour to be the
                    // first sample column.  FilterMutectCalls identifies
                    // normal samples from ##normal_sample and sums every
                    // remaining genotype.  In particular, the common
                    // normal,tumour column order must not train the somatic
                    // model on the normal AD vector.
                    const auto depths = tumor_allele_depths(
                        header, record, tumor_sample_indices);
                    int total = 0;
                    for (const int depth : depths)
                        if (depth > 0) total += depth;
                    const int alt =
                        static_cast<int>(depths.size()) == alts
                            ? depths[static_cast<std::size_t>(alt_index)]
                            : -1;
                    if (total <= 0 || alt < 0) continue;
                    double tlo = 8.0 * fmc_j::kLn10;
                    if (alt_index - 1 < static_cast<int>(tlod_values.size()))
                        tlo = static_cast<double>(tlod_values[
                            static_cast<std::size_t>(alt_index - 1)]) *
                            fmc_j::kLn10;
                    double pop_af = 0.01;
                    bool has_pop_af = false;
                    if (!popaf_values.empty() &&
                        alt_index - 1 < static_cast<int>(popaf_values.size()) &&
                        std::isfinite(static_cast<double>(
                            popaf_values[static_cast<std::size_t>(alt_index - 1)]))) {
                        pop_af = std::pow(10.0, -static_cast<double>(
                            popaf_values[static_cast<std::size_t>(alt_index - 1)]));
                        has_pop_af = true;
                    }
                    double nlod = 0.0;
                    if (alt_index - 1 < static_cast<int>(nlod_values.size()) &&
                        std::isfinite(static_cast<double>(nlod_values[
                            static_cast<std::size_t>(alt_index - 1)])))
                        nlod = static_cast<double>(nlod_values[
                            static_cast<std::size_t>(alt_index - 1)]);
                    const double orientation_artifact = orientation_evaluation.has_prior &&
                            alt_index == orientation_alt
                        ? orientation_evaluation.probability : 0.0;
                    double hard_artifact = 0.0;
                    const auto alt_mbq = info_allele_int(header, record, "MBQ", alt_index, true);
                    if (alt_mbq && options.min_median_base_quality >= 0 &&
                        *alt_mbq < options.min_median_base_quality)
                        hard_artifact = 1.0;
                    auto alt_mmq = info_allele_int(header, record, "MMQ", alt_index, true);
                    const auto indel_length = indel_length_for_alt(record, alt_index);
                    if (indel_length >= options.long_indel_length && !mmq_values.empty())
                        alt_mmq = mmq_values.front();
                    if (alt_mmq && options.min_median_mapping_quality >= 0 &&
                        *alt_mmq < options.min_median_mapping_quality)
                        hard_artifact = 1.0;
                    const auto alt_mpos = info_allele_int(header, record, "MPOS", alt_index, false);
                    if (alt_mpos && options.min_median_read_position >= 0 && *alt_mpos > -1 &&
                        *alt_mpos < options.min_median_read_position)
                        hard_artifact = 1.0;
                    if (options.unique_alt_read_count > 0) {
                        const auto unique_count = info_alt_int(
                            header, record, "AS_UNIQ_ALT_READ_COUNT", alt_index - 1);
                        if (unique_count && *unique_count <= options.unique_alt_read_count)
                            hard_artifact = 1.0;
                    }
                    if (options.min_reads_per_strand > 0 && strand_counts_by_allele &&
                        static_cast<std::size_t>(alt_index) < strand_counts_by_allele->size()) {
                        const auto& counts = (*strand_counts_by_allele)[static_cast<std::size_t>(alt_index)];
                        if (counts.first == 0 || counts.second == 0) hard_artifact = 1.0;
                    }
                    if (options.min_allele_fraction > 0.0F &&
                        tumor_allele_fraction_for_alt(header, record, tumor_sample_indices.front(), alt_index) <
                            options.min_allele_fraction)
                        hard_artifact = 1.0;
                    if (std::isfinite(options.max_n_ratio) && n_count) {
                        const auto all_depths = tumor_allele_depths(header, record, tumor_sample_indices);
                        int all_alt_depth = 0;
                        for (std::size_t index = 1; index < all_depths.size(); ++index)
                            if (all_depths[index] > 0) all_alt_depth += all_depths[index];
                        if (all_alt_depth > 0 && static_cast<double>(*n_count) /
                            static_cast<double>(all_alt_depth) >= options.max_n_ratio)
                            hard_artifact = 1.0;
                    }
                    if (panel_of_normals || clustered_artifact || multiallelic_artifact ||
                        fragment_artifact)
                        hard_artifact = 1.0;
                    const auto strand_iterator = strand_artifact_by_alt.find(alt_index);
                    const double strand_artifact = strand_iterator == strand_artifact_by_alt.end()
                        ? 0.0 : strand_iterator->second;
                    const double non_normal_non_haplotype_artifact = std::max(
                        {orientation_artifact, strand_artifact, hard_artifact});
                    const double art = std::max(
                        {non_normal_non_haplotype_artifact,
                         normal_artifact_probability, learned_haplotype_artifact});
                    ev_alt.push_back(alt);
                    ev_total.push_back(total);
                    ev_tlo.push_back(tlo);
                    ev_pop_af.push_back(pop_af);
                    ev_nlod.push_back(nlod);
                    ev_art.push_back(art);
                    ev_non_normal_non_haplotype_art.push_back(
                        non_normal_non_haplotype_artifact);
                    ev_normal_art.push_back(normal_artifact_probability);
                    ev_indel_length.push_back(indel_length);
                    ev_germ_gate.push_back(has_pop_af ? 1 : 0);
                }
                if (haplotype_filter != nullptr) {
                    haplotype_pass_records.push_back(HaplotypePassRecord{
                        evidence_begin, ev_alt.size(),
                        static_cast<std::int64_t>(record->pos) + 1,
                        phased_observations});
                }
                bcf_clear(record);
            }
            bcf_destroy(record);
            bcf_hdr_destroy(header);
            hts_close(input);
        }
        if (!ev_alt.empty()) {
            // Phase 2: batch sequencing/germline/contamination posteriors in
            // the kernels module under the current (pre-pass) model state.
            fastgatk::kernels::FmcClusterParams params;
            params.log_cluster_weights = model.log_cluster_weights;
            params.cluster_alpha.assign(
                fastgatk::kernels::kFmcMaxClusters, 1.0);
            params.cluster_beta.assign(
                fastgatk::kernels::kFmcMaxClusters, 1.0);
            params.cluster_kind.assign(
                fastgatk::kernels::kFmcMaxClusters, 0);
            for (std::size_t i = 0; i < model.clusters.size(); ++i) {
                params.cluster_alpha[i] = model.clusters[i].alpha;
                params.cluster_beta[i] = model.clusters[i].beta;
                params.cluster_kind[i] = model.clusters[i].kind;
            }
            fastgatk::kernels::FmcEvidence evidence;
            evidence.total_count = ev_total;
            evidence.alt_count = ev_alt;
            evidence.tumor_log_odds = ev_tlo;
            evidence.artifact_prob = ev_art;
            evidence.non_somatic_prob.assign(ev_alt.size(), 0.0);
            evidence.germline_evaluated = ev_germ_gate;
            evidence.population_af = ev_pop_af;
            evidence.normal_log10_odds = ev_nlod;
            evidence.contamination.assign(ev_alt.size(), contamination);
            evidence.indel_length = ev_indel_length;
            std::vector<double> priors;
            priors.reserve(21);
            for (const double value : model.log_variant_priors)
                priors.push_back(value);
            const auto ep = fastgatk::kernels::fmc_error_probabilities_kokkos(
                params, evidence, priors);
            if (!ep.valid)
                throw std::runtime_error(
                    "BAD_INPUT: Java-model error-probability kernel rejected "
                    "the accumulated FilterMutectCalls evidence");
            if (ev_non_normal_non_haplotype_art.size() != ev_alt.size() ||
                ev_normal_art.size() != ev_alt.size()) {
                throw std::runtime_error(
                    "internal FilterMutectCalls artifact-source accounting mismatch");
            }
            // Phase 3: Java record() exclusions and threshold errors in VCF
            // order (record() drops non-somatic > 0.9 but the pass-2 error is
            // computed before that exclusion, as in GATK).
            for (std::size_t i = 0; i < ev_alt.size(); ++i) {
                const double seq = ep.sequencing_error[i];
                const double ns_value = ep.non_somatic_error[i];
                const double art_value = ev_art[i];
                if (errors_ptr != nullptr)
                    errors_ptr->push_back(
                        1.0 - (1.0 - seq) * (1.0 - ns_value) *
                                  (1.0 - art_value));
                model.record(ev_tlo[i], art_value, ns_value, ev_alt[i],
                             ev_total[i], ev_indel_length[i]);
            }
            if (haplotype_filter != nullptr) {
                // FilteredHaplotypeFilter sees all filter values for a VCF
                // record, not just the ALT that later survives model.record.
                // Its exclusion of normal_artifact depends only on the
                // maximum *GermlineFilter* posterior (not contamination).
                for (const auto& record_state : haplotype_pass_records) {
                    double maximum_germline = 0.0;
                    for (std::size_t index = record_state.evidence_begin;
                         index < record_state.evidence_end; ++index) {
                        if (index < ep.germline_error.size() &&
                            std::isfinite(ep.germline_error[index])) {
                            maximum_germline = std::max(
                                maximum_germline, ep.germline_error[index]);
                        }
                    }
                    const bool ignore_normal_artifact = maximum_germline > 0.25;
                    double artifact_probability = 0.0;
                    for (std::size_t index = record_state.evidence_begin;
                         index < record_state.evidence_end; ++index) {
                        artifact_probability = std::max(
                            artifact_probability,
                            ev_non_normal_non_haplotype_art[index]);
                        if (!ignore_normal_artifact) {
                            artifact_probability = std::max(
                                artifact_probability, ev_normal_art[index]);
                        }
                    }
                    for (const auto& observation : record_state.observations) {
                        haplotype_filter->accumulate(
                            observation.key, record_state.locus, artifact_probability);
                    }
                }
            }
        }
        if (learn_after) {
            // Mutect2FilteringEngine learns every filter before it advances
            // SomaticClusteringModel.  This makes pass 2 consume the exact
            // StrandArtifact/FilteredHaplotype state learned from pass 1.
            strand_artifact_filter.learn_parameters_and_clear();
            if (haplotype_filter != nullptr)
                haplotype_filter->learn_parameters_and_clear();
            model.learn_and_clear();
        }
    };

    auto dump_model = [&](const char* tag) {
        if (std::getenv("FASTGATK_DEBUG_FMC") == nullptr) return;
        std::fprintf(stderr, "[FMCSTATE] %s k=%zu w=", tag,
                     model.clusters.size());
        for (std::size_t i = 0; i < model.log_cluster_weights.size(); ++i)
            std::fprintf(stderr, "%.17g,", model.log_cluster_weights[i]);
        std::fprintf(stderr, " ab=");
        for (std::size_t i = 0; i < model.clusters.size(); ++i)
            std::fprintf(stderr, "%.17g/%.17g/%d,",
                         model.clusters[i].alpha, model.clusters[i].beta,
                         model.clusters[i].kind);
        std::fprintf(stderr, " p10=%.17g vv=%.17g obs=%d\n",
                     model.log_variant_priors[10],
                     model.log_variant_vs_artifact_prior,
                     model.obvious_artifact_count);
    };
    accumulate_pass(true);
    dump_model("pass1");
    accumulate_pass(true);
    dump_model("pass2");
    std::vector<double> errors;
    accumulate_pass(false, &errors);
    dump_model("pass3");
    // The final (non-learning) pass leaves model.data populated; report its
    // size as the empirical somatic model record count so the stats JSON no
    // longer reads 0 while the model is in fact learned.
    if (empirical_records_out != nullptr)
        *empirical_records_out = static_cast<std::uint64_t>(model.data.size());

    if (std::getenv("FASTGATK_DEBUG_FMC") != nullptr) {
        std::fprintf(stderr, "[FMCERRORS]");
        for (const double value : errors) std::fprintf(stderr, " %.17g", value);
        std::fprintf(stderr, "\n");
    }

    effective_threshold = fmc_j::optimal_f_score_threshold(errors);
    observations_out = errors.size();
    if (threshold_errors_out != nullptr)
        *threshold_errors_out = std::move(errors);

    // Transfer into native model parameters used by the final filtering loop.
    native_model.cluster_count = model.clusters.size();
    for (std::size_t i = 0; i < model.clusters.size() &&
                            i < SomaticModelParameters::kMaxClusters; ++i) {
        native_model.cluster_log_weights[i] = model.log_cluster_weights[i];
        native_model.cluster_alpha[i] = model.clusters[i].alpha;
        native_model.cluster_beta[i] = model.clusters[i].beta;
        native_model.cluster_mean[i] = model.clusters[i].alpha /
            (model.clusters[i].alpha + model.clusters[i].beta);
        native_model.cluster_binomial[i] = 0;
    }
    for (int i = 0; i < 21; ++i)
        native_model.log_variant_priors[static_cast<std::size_t>(i)] =
            model.log_variant_priors[static_cast<std::size_t>(i)];
    native_model.log_snv_prior = model.log_variant_priors[10];
    native_model.log_indel_prior = model.log_variant_priors[9];
    native_model.log_variant_vs_artifact_prior = model.log_variant_vs_artifact_prior;
    native_model.learned = true;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    htsFile* input = nullptr;
    htsFile* output = nullptr;
    bcf_hdr_t* input_header = nullptr;
    bcf_hdr_t* output_header = nullptr;
    bcf1_t* record = nullptr;
    std::uint64_t input_records = 0;
    std::uint64_t pass_records = 0;
    std::uint64_t interval_skipped_records = 0;
    std::uint64_t low_tlod_records = 0;
    std::uint64_t germline_records = 0;
    std::uint64_t contamination_records = 0;
    std::uint64_t contamination_posterior_records = 0;
    std::uint64_t contamination_posterior_alleles = 0;
    std::uint64_t orientation_records = 0;
    std::uint64_t low_af_records = 0;
    std::uint64_t strand_records = 0;
    std::uint64_t low_alt_read_records = 0;
    std::uint64_t multiallelic_records = 0;
    std::uint64_t low_median_base_quality_records = 0;
    std::uint64_t low_median_mapping_quality_records = 0;
    std::uint64_t low_median_read_position_records = 0;
    std::uint64_t fragment_length_records = 0;
    std::uint64_t n_ratio_records = 0;
    std::uint64_t allele_specific_strand_records = 0;
    std::uint64_t allele_specific_unique_records = 0;
    std::uint64_t strict_strand_records = 0;
    std::uint64_t duplicate_records = 0;
    std::uint64_t panel_of_normals_records = 0;
    std::uint64_t as_filter_records = 0;
    std::uint64_t normal_artifact_records = 0;
    std::uint64_t weak_evidence_records = 0;
    std::uint64_t slippage_records = 0;
    std::uint64_t low_somatic_probability_records = 0;
    std::uint64_t high_germline_probability_records = 0;
    std::uint64_t high_artifact_probability_records = 0;
    std::uint64_t orientation_prior_records = 0;
    std::uint64_t orientation_prior_sample_evaluations = 0;
    std::uint64_t orientation_threshold_observations = 0;
    std::uint64_t orientation_depth_disagreements = 0;
    std::uint64_t joint_error_threshold_observations = 0;
    std::uint64_t joint_error_records = 0;
    std::uint64_t empirical_somatic_model_records = 0;
    std::uint64_t clustered_events_records = 0;
    std::uint64_t haplotype_records = 0;
    FilterKernelTelemetry kernel_telemetry;
    float effective_contamination = options.max_contamination;
    std::string contamination_source = options.max_contamination >= 0.0F ? "argument" : "none";
    std::unordered_map<std::string, float> contamination_by_sample;
    fastgatk::io::IntervalFileStats interval_stats;
    std::vector<fastgatk::io::IndexedInterval> intervals;
    float effective_orientation_threshold = options.max_orientation_artifact_probability;
    double effective_joint_error_threshold = -1.0;
    SomaticModelParameters empirical_somatic_model;
    initialize_somatic_model_clusters(empirical_somatic_model);
    NativeFilteringOutputStats filtering_output_stats;
    empirical_somatic_model.log_snv_prior = options.log_snv_prior;
    empirical_somatic_model.log_indel_prior = options.log_indel_prior;
    empirical_somatic_model.log_variant_vs_artifact_prior = options.log_artifact_prior;
    for (auto& prior : empirical_somatic_model.log_variant_priors) prior = options.log_indel_prior;
    empirical_somatic_model.log_variant_priors[10] = options.log_snv_prior;
    // FilterMutectCalls derives the Mutect2 summary path from the driving VCF
    // when --stats is omitted.  The source implementation uses exactly
    // ``drivingVariantFile + ".stats"``; using only an explicit --stats here
    // silently skipped all four learning/calling passes for ordinary GATK
    // Mutect2 output.
    const std::string mutect_stats_path = options.stats.empty()
        ? options.input + ".stats" : options.stats;
    const bool input_stats_is_gatk_table = looks_like_gatk_mutect_stats(mutect_stats_path);
    // A real GATK FilterMutectCalls invocation normally receives Mutect2's
    // stats table and therefore runs the default OPTIMAL_F_SCORE empirical
    // pass even when --threshold-strategy was omitted.  Synthetic/legacy VCF
    // calls without that table keep the explicit adapter behavior.
    const bool empirical_learning_requested = options.threshold_strategy_explicit || input_stats_is_gatk_table;
    const auto mutect_stats = input_stats_is_gatk_table
        ? read_mutect_stats_summary(mutect_stats_path) : MutectStatsSummary{};
    if (std::isfinite(mutect_stats.callable_sites) && mutect_stats.callable_sites >= 1.0)
        empirical_somatic_model.callable_sites = mutect_stats.callable_sites;
    const std::string stats_input_path = input_stats_is_gatk_table ? mutect_stats_path : std::string{};
    const std::string stats_output_path = !options.filtering_stats.empty()
        ? options.filtering_stats
        : (input_stats_is_gatk_table || options.stats.empty()
            ? options.output + ".stats.json" : options.stats);
    std::vector<OrientationPriorCollection> orientation_collections;
    faidx_t* reference_index = nullptr;
    try {
        if (!options.orientation_priors.empty()) {
            if (options.reference.empty())
                throw std::invalid_argument("-R/--reference is required with --orientation-bias-artifact-priors");
            orientation_collections = read_orientation_prior_archives(options.orientation_priors);
            if (effective_orientation_threshold < 0.0F) effective_orientation_threshold = 0.5F;
            reference_index = fai_load(options.reference.c_str());
            if (reference_index == nullptr)
                throw std::runtime_error("BAD_INPUT: cannot load reference FASTA index: " + options.reference);
        }
        if (!options.contamination_table.empty()) {
            // GATK always reads the per-sample table, even when an explicit
            // contamination estimate supplies the fallback for samples that
            // are absent from it.  A matching table row therefore overrides
            // the global estimate in contamination_probabilities_for_record.
            const float table_default = read_contamination_estimate(
                options.contamination_table, false, &contamination_by_sample);
            if (effective_contamination < 0.0F) effective_contamination = table_default;
            contamination_source = "contamination-table";
        }
        if (effective_contamination < 0.0F && !options.tumor_segmentation.empty()) {
            effective_contamination = read_contamination_estimate(options.tumor_segmentation, true);
            contamination_source = "tumor-segmentation";
        }
        input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open Mutect2 VCF: " + options.input);
        input_header = bcf_hdr_read(input);
        if (!input_header) throw std::runtime_error("BAD_INPUT: cannot read Mutect2 VCF header");
        // GATK rejects a genuine Mutect2 callset when its companion stats
        // table is absent.  Retain the adapter's legacy synthetic-VCF mode
        // (which has no Mutect2 provenance) for focused hard-filter tests,
        // while making production Mutect2 input follow the source contract.
        if (!input_stats_is_gatk_table && header_declares_mutect2(input_header)) {
            throw std::runtime_error(
                "BAD_INPUT: Mutect stats table " + mutect_stats_path +
                " not found or is not a GATK Mutect2 stats table");
        }
        bool first_interval_selector = true;
        for (const auto& selector : options.regions) {
            append_filter_interval_selector_with_rule(
                selector, input_header, intervals, interval_stats,
                options.interval_set_rule, first_interval_selector);
            first_interval_selector = false;
        }
        normalize_filter_intervals(intervals);
        output_header = bcf_hdr_dup(input_header);
        if (!output_header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate VCF header");
        if (options.sites_only_vcf_output) {
            // NULL means exclude all samples.  HTSlib keeps a mask so
            // bcf_subset_format() can strip each record's FORMAT payload at
            // the writer boundary without affecting input-side calculations.
            if (bcf_hdr_set_samples(output_header, nullptr, 0) != 0)
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: cannot configure sites-only VCF header");
        }
        // FilterMutectCalls replaces Mutect2's provisional status and emits
        // the canonical GATK filter vocabulary.  Do not leave native-only
        // filter names or descriptions in the public VCF schema: downstream
        // tools use these exact IDs to interpret AS_FilterStatus.
        bcf_hdr_remove(output_header, BCF_HL_GEN, "filtering_status");
        bcf_hdr_remove(output_header, BCF_HL_GEN, "source");
        bcf_hdr_remove(output_header, BCF_HL_INFO, "AS_FilterStatus");
        bcf_hdr_remove(output_header, BCF_HL_INFO, "CONTQ");
        bcf_hdr_remove(output_header, BCF_HL_INFO, "GERMQ");
        bcf_hdr_remove(output_header, BCF_HL_INFO, "ROQ");
        bcf_hdr_remove(output_header, BCF_HL_INFO, "STRQ");
        static constexpr std::array<const char*, 21> kMutectFilterNames{{
            "FAIL", "PASS", "base_qual", "clustered_events", "contamination",
            "duplicate", "fragment", "germline", "haplotype", "low_allele_frac",
            "map_qual", "multiallelic", "n_ratio", "normal_artifact", "orientation",
            "panel_of_normals", "position", "possible_numt", "slippage", "strand_bias",
            "strict_strand"}};
        for (const auto* name : kMutectFilterNames)
            bcf_hdr_remove(output_header, BCF_HL_FLT, name);
        // weak_evidence is listed separately to keep the array above aligned
        // with the fixed GATK filter-header ordering used by the 4.6 source.
        bcf_hdr_remove(output_header, BCF_HL_FLT, "weak_evidence");
        static constexpr std::array<const char*, 22> kMutectFilterHeaders{{
            "##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">",
            "##FILTER=<ID=PASS,Description=\"Site contains at least one allele that passes filters\">",
            "##FILTER=<ID=base_qual,Description=\"alt median base quality\">",
            "##FILTER=<ID=clustered_events,Description=\"Clustered events observed in the tumor\">",
            "##FILTER=<ID=contamination,Description=\"contamination\">",
            "##FILTER=<ID=duplicate,Description=\"evidence for alt allele is overrepresented by apparent duplicates\">",
            "##FILTER=<ID=fragment,Description=\"abs(ref - alt) median fragment length\">",
            "##FILTER=<ID=germline,Description=\"Evidence indicates this site is germline, not somatic\">",
            "##FILTER=<ID=haplotype,Description=\"Variant near filtered variant on same haplotype.\">",
            "##FILTER=<ID=low_allele_frac,Description=\"Allele fraction is below specified threshold\">",
            "##FILTER=<ID=map_qual,Description=\"ref - alt median mapping quality\">",
            "##FILTER=<ID=multiallelic,Description=\"Site filtered because too many alt alleles pass tumor LOD\">",
            "##FILTER=<ID=n_ratio,Description=\"Ratio of N to alt exceeds specified ratio\">",
            "##FILTER=<ID=normal_artifact,Description=\"artifact_in_normal\">",
            "##FILTER=<ID=orientation,Description=\"orientation bias detected by the orientation bias mixture model\">",
            "##FILTER=<ID=panel_of_normals,Description=\"Blacklisted site in panel of normals\">",
            "##FILTER=<ID=position,Description=\"median distance of alt variants from end of reads\">",
            "##FILTER=<ID=possible_numt,Description=\"Allele depth is below expected coverage of NuMT in autosome\">",
            "##FILTER=<ID=slippage,Description=\"Site filtered due to contraction of short tandem repeat region\">",
            "##FILTER=<ID=strand_bias,Description=\"Evidence for alt allele comes from one read direction only\">",
            "##FILTER=<ID=strict_strand,Description=\"Evidence for alt allele is not represented in both directions\">",
            "##FILTER=<ID=weak_evidence,Description=\"Mutation does not meet likelihood threshold\">"}};
        for (const auto* line : kMutectFilterHeaders)
            bcf_hdr_append(output_header, line);
        bcf_hdr_append(output_header,
            "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=\"Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles.\">");
        bcf_hdr_append(output_header,
            "##INFO=<ID=CONTQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to contamination\">");
        bcf_hdr_append(output_header,
            "##INFO=<ID=GERMQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not germline variants\">");
        bcf_hdr_append(output_header,
            "##INFO=<ID=ROQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to read orientation artifact\">");
        bcf_hdr_append(output_header,
            "##INFO=<ID=STRQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors\">");
        bcf_hdr_append(output_header,
            "##filtering_status=These calls have been filtered by FilterMutectCalls to label false positives with a list of failed filters and true positives with PASS.");
        bcf_hdr_append(output_header, "##source=FilterMutectCalls");
        bcf_hdr_append(output_header, "##source=Mutect2");
        if (bcf_hdr_sync(output_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync filtered VCF header");
        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open filtered VCF: " + options.output);
        if (bcf_hdr_write(output, output_header) != 0)
            throw std::runtime_error("cannot write filtered VCF header");
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
        const auto normal_samples = normal_sample_names(input_header);
        std::vector<int> normal_sample_indices;
        std::vector<int> tumor_sample_indices;
        for (int index = 0; index < bcf_hdr_nsamples(input_header); ++index) {
            const auto* sample_name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, index);
            if (sample_name != nullptr && normal_samples.count(sample_name) != 0)
                normal_sample_indices.push_back(index);
            else
                tumor_sample_indices.push_back(index);
        }
        int tumor_sample_index = tumor_sample_indices.empty() ? 0 : tumor_sample_indices.front();
        const std::string requested_tumor = options.tumor_sample.empty() ? "TUMOR" : options.tumor_sample;
        if (!options.tumor_sample.empty()) {
            const auto requested_index = bcf_hdr_id2int(
                input_header, BCF_DT_SAMPLE, requested_tumor.c_str());
            if (requested_index < 0)
                throw std::runtime_error("BAD_INPUT: tumor sample not found in VCF header: " + options.tumor_sample);
            tumor_sample_index = requested_index;
        }
        const std::string selected_tumor_name =
            tumor_sample_index >= 0 && tumor_sample_index < bcf_hdr_nsamples(input_header)
                ? bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, tumor_sample_index) : requested_tumor;
        if (!options.tumor_sample.empty()) {
            tumor_sample_indices.clear();
            tumor_sample_indices.push_back(tumor_sample_index);
        }
        if (tumor_sample_indices.empty()) tumor_sample_indices.push_back(tumor_sample_index);
        std::vector<int> orientation_sample_indices;
        if (!options.orientation_priors.empty()) {
            if (!options.tumor_sample.empty()) {
                orientation_sample_indices.push_back(tumor_sample_index);
            } else {
                for (int index = 0; index < bcf_hdr_nsamples(input_header); ++index) {
                    const auto* sample_name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, index);
                    if (sample_name != nullptr && normal_samples.count(sample_name) == 0)
                        orientation_sample_indices.push_back(index);
                }
            }
        }
        Kokkos::View<double*> orientation_prior_device("filter_orientation_prior", 12);
        Kokkos::View<double*> orientation_log_device("filter_orientation_log_posterior", 12);
        if (!orientation_collections.empty() && options.threshold_strategy_explicit &&
            options.orientation_threshold_strategy == "CONSTANT" &&
            options.max_orientation_artifact_probability < 0.0F)
            effective_orientation_threshold = static_cast<float>(options.initial_threshold);
        if (!orientation_collections.empty() && options.threshold_strategy_explicit &&
            options.orientation_threshold_strategy != "CONSTANT" &&
            options.max_orientation_artifact_probability < 0.0F) {
            std::vector<double> learned_posteriors;
            while (fastgatk::io::read_variant_record(input, input_header, record,
                                                   options.input) == 0) {
                if (!options.regions.empty() &&
                    !filter_record_in_intervals(input_header, record, intervals)) {
                    bcf_clear(record);
                    continue;
                }
                bcf_unpack(record, BCF_UN_ALL);
                const auto evaluation = evaluate_orientation_posterior(
                    input_header, record, reference_index, orientation_collections,
                    orientation_sample_indices, orientation_prior_device, orientation_log_device,
                    &kernel_telemetry);
                if (evaluation.has_prior) learned_posteriors.push_back(evaluation.probability);
                bcf_clear(record);
            }
            orientation_threshold_observations = learned_posteriors.size();
            effective_orientation_threshold = static_cast<float>(
                learn_orientation_threshold(options, learned_posteriors));
            bcf_close(input);
            input = nullptr;
            bcf_hdr_destroy(input_header);
            input_header = nullptr;
            input = bcf_open(options.input.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen Mutect2 VCF for threshold learning");
            input_header = bcf_hdr_read(input);
            if (!input_header) throw std::runtime_error("BAD_INPUT: cannot reread Mutect2 VCF header");
        }
        LearnedHaplotypeFilter haplotype_filter(options.max_intra_haplotype_distance);
        LearnedStrandArtifactFilter strand_artifact_filter;
        const bool has_phasing_fields = bcf_hdr_id2int(input_header, BCF_DT_ID, "PGT") >= 0 &&
                                        bcf_hdr_id2int(input_header, BCF_DT_ID, "PID") >= 0;
        std::vector<int> haplotype_sample_indices;
        if (has_phasing_fields) {
            if (!options.tumor_sample.empty()) {
                haplotype_sample_indices.push_back(tumor_sample_index);
            } else {
                // FilteredHaplotypeFilter sees every tumour genotype, not a
                // hard-coded first sample.  ##normal_sample metadata is the
                // same boundary used by the orientation-prior adapter.
                for (int index = 0; index < bcf_hdr_nsamples(input_header); ++index) {
                    const auto* sample_name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, index);
                    if (sample_name != nullptr && normal_samples.count(sample_name) == 0)
                        haplotype_sample_indices.push_back(index);
                }
                if (haplotype_sample_indices.empty())
                    haplotype_sample_indices.push_back(tumor_sample_index);
            }
        }
        if (has_phasing_fields && !empirical_learning_requested) {
            // Seed the first FilteredHaplotypeFilter state from deterministic
            // artifact evidence.  In particular, ErrorProbabilities carries
            // ClusteredEventsFilter as ARTIFACT=1; using only the legacy
            // PARTIFACT INFO annotation here lost that signal and made every
            // clustered phased call evade GATK's haplotype filter.
            //
            // This remains a streaming pass: only per-haplotype
            // (locus, probability) vectors remain resident.
            while (fastgatk::io::read_variant_record(input, input_header, record,
                                                   options.input) == 0) {
                if (!options.regions.empty() &&
                    !filter_record_in_intervals(input_header, record, intervals)) {
                    bcf_clear(record);
                    continue;
                }
                bcf_unpack(record, BCF_UN_ALL);
                const auto observations = phased_haplotype_observations(
                    input_header, record, haplotype_sample_indices);
                const auto region_event_count = info_int(input_header, record, "ECNT", -1);
                const auto haplotype_event_count = info_int_max(input_header, record, "ECNTH");
                const bool clustered_artifact = !options.mitochondria_mode &&
                    !options.microbial_mode &&
                    ((region_event_count >= 0 && region_event_count > options.max_events_in_region) ||
                     (haplotype_event_count && *haplotype_event_count > options.max_events_in_haplotype));
                const auto probability = info_numeric(input_header, record, "PARTIFACT", 0.0F);
                const auto finite_probability = clustered_artifact ? 1.0 :
                    (std::isfinite(probability)
                        ? std::clamp(static_cast<double>(probability), 0.0, 1.0) : 0.0);
                for (const auto& observation : observations)
                    haplotype_filter.accumulate(
                        observation.key, static_cast<std::int64_t>(record->pos) + 1,
                        finite_probability);
                bcf_clear(record);
            }
            haplotype_filter.learn_parameters_and_clear();
            bcf_close(input);
            input = nullptr;
            bcf_hdr_destroy(input_header);
            input_header = nullptr;
            input = bcf_open(options.input.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen Mutect2 VCF for haplotype learning");
            input_header = bcf_hdr_read(input);
            if (!input_header) throw std::runtime_error("BAD_INPUT: cannot reread Mutect2 VCF header");
        }
        if (empirical_learning_requested) {
            // GATK's ThresholdCalculator learns from the combined per-ALT
            // ErrorProbabilities before the final pass.  Reconstruct the
            // exposed native error types (PSOMATIC/TLOD, PGERMLINE,
            // PARTIFACT/OBP) in a streaming pass so the threshold strategy
            // does not retain an entire VCF in memory.
            // GATK always makes two parameter-learning passes followed by a
            // threshold pass before writing calls.  This is the production
            // path, not an opt-in compatibility mode.
            const bool java_passes = true;
            if (java_passes) {
                std::size_t java_observations = 0;
                std::vector<double> java_threshold_errors;
                run_java_learning_pipeline(
                    options, empirical_somatic_model,
                    effective_joint_error_threshold, java_observations,
                    &java_threshold_errors, mutect_stats_path, reference_index,
                    orientation_collections, orientation_sample_indices, tumor_sample_indices,
                    haplotype_sample_indices,
                    has_phasing_fields ? &haplotype_filter : nullptr,
                    strand_artifact_filter,
                    orientation_prior_device, orientation_log_device,
                    contamination_by_sample, effective_contamination,
                    &kernel_telemetry, &empirical_somatic_model_records);
                joint_error_threshold_observations = java_observations;
                // GATK ThresholdCalculator only learns the threshold for the
                // OPTIMAL_F_SCORE strategy; CONSTANT and FALSE_DISCOVERY_RATE
                // resolve from the same per-ALT errors, exactly as the native
                // learning branch below does.
                if (options.joint_threshold_strategy == "CONSTANT") {
                    effective_joint_error_threshold =
                        options.initial_threshold;
                } else if (options.joint_threshold_strategy ==
                           "FALSE_DISCOVERY_RATE") {
                    effective_joint_error_threshold =
                        calculate_false_discovery_rate_threshold(
                            java_threshold_errors,
                            options.max_false_discovery_rate);
                }
            } else {
            std::vector<double> learned_joint_errors;
            std::vector<SomaticModelDatum> empirical_model_data;
            std::vector<double> threshold_errors;
            const int learning_total = java_passes ? 3 : 1;
            for (int learning_pass = 0; learning_pass < learning_total; ++learning_pass) {
                learned_joint_errors.clear();
                empirical_model_data.clear();
                if (learning_pass > 0) {
                    bcf_close(input);
                    input = nullptr;
                    bcf_hdr_destroy(input_header);
                    input_header = nullptr;
                    input = bcf_open(options.input.c_str(), "r");
                    if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen Mutect2 VCF for ErrorProbabilities learning");
                    input_header = bcf_hdr_read(input);
                    if (!input_header) throw std::runtime_error("BAD_INPUT: cannot reread Mutect2 VCF header");
                }
                while (fastgatk::io::read_variant_record(input, input_header,
                                                       record, options.input) == 0) {
                if (!options.regions.empty() &&
                    !filter_record_in_intervals(input_header, record, intervals)) {
                    bcf_clear(record);
                    continue;
                }
                bcf_unpack(record, BCF_UN_ALL);
                const auto tumor_depths = tumor_allele_depths(
                    input_header, record, tumor_sample_indices);
                const auto tlod_values = info_float_values(input_header, record, "TLOD");
                // The Java ThresholdCalculator consumes ErrorProbabilities
                // after all enabled filters have contributed.  In particular,
                // ContaminationFilter is a NON_SOMATIC contributor; include
                // its per-ALT posterior in this learning pass as well as the
                // final filtering pass.
                const auto contamination_probabilities =
                    contamination_probabilities_for_record_kokkos(
                        input_header, record, effective_contamination,
                        contamination_by_sample, tumor_sample_indices,
                        empirical_somatic_model, &kernel_telemetry);
                // ReadOrientationFilter is a variant-level ARTIFACT
                // contributor to ErrorProbabilities.  Evaluate the same
                // weighted-median F1R2/F2R1 posterior during threshold
                // learning instead of waiting until final output.
                const auto orientation_learning_evaluation =
                    evaluate_orientation_posterior(
                        input_header, record, reference_index, orientation_collections,
                        orientation_sample_indices, orientation_prior_device,
                        orientation_log_device, &kernel_telemetry);
                const auto orientation_learning_posterior =
                    orientation_learning_evaluation.has_prior
                    ? orientation_learning_evaluation.probability
                    : std::numeric_limits<double>::quiet_NaN();
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    if (record->d.allele[alt_index] == nullptr || record->d.allele[alt_index][0] == '<')
                        continue;
                    if (tumor_depths.size() == static_cast<std::size_t>(record->n_allele)) {
                        int total_depth = 0;
                        for (const auto depth : tumor_depths) if (depth > 0) total_depth += depth;
                        const auto alt_depth = tumor_depths[static_cast<std::size_t>(alt_index)];
                        // Legacy note: java_passes is constant false here (the
                        // Java path returns via run_java_learning_pipeline
                        // above), so this branch is dead control flow kept for
                        // historical reference; the else below is the active
                        // INFO-driven native-default semantics.
                        if (java_passes) {
                        // Java record() semantics (only in the isolated
                        // JAVA_PASSES path): recompute the group errors from
                        // counts when INFO is absent, exclude ns/art > 0.9 and
                        // weight by (1-art)(1-ns)(1-seq).
                        const auto artifact_info = info_probability_at(
                            input_header, record, "PARTIFACT", alt_index);
                        const auto orientation_info = info_probability_at(
                            input_header, record, "OBP", alt_index);
                        const auto artifact_probability = std::isfinite(artifact_info)
                                && std::isfinite(orientation_info)
                            ? std::max(artifact_info, orientation_info)
                            : (std::isfinite(artifact_info) ? artifact_info
                                                           : orientation_info);
                        double sequencing_error = std::numeric_limits<double>::quiet_NaN();
                        const auto posterior_info = info_probability_at(
                            input_header, record, "PSOMATIC", alt_index);
                        if (std::isfinite(posterior_info))
                            sequencing_error = 1.0 - posterior_info;
                        const auto tlo_alt = alt_index - 1 < static_cast<int>(tlod_values.size())
                            ? static_cast<double>(tlod_values[
                                static_cast<std::size_t>(alt_index - 1)])
                            : std::numeric_limits<double>::quiet_NaN();
                        const auto reference_length = record->d.allele[0] == nullptr
                            ? 0 : static_cast<int>(std::strlen(record->d.allele[0]));
                        const auto alternate_length = record->d.allele[alt_index] == nullptr
                            ? reference_length
                            : static_cast<int>(std::strlen(record->d.allele[alt_index]));
                        const auto indel_length = alternate_length - reference_length;
                        double log_somatic_prior = -6.0 * std::log(10.0) - std::log(3.0);
                        if (indel_length >= -10 && indel_length <= 10 &&
                            std::isfinite(empirical_somatic_model.log_variant_priors[
                                static_cast<std::size_t>(indel_length + 10)])) {
                            log_somatic_prior = empirical_somatic_model.log_variant_priors[
                                static_cast<std::size_t>(indel_length + 10)] +
                                (indel_length == 0 ? -std::log(3.0) : 0.0);
                        }
                        if (!std::isfinite(sequencing_error) && std::isfinite(tlo_alt)) {
                            sequencing_error = empirical_tumor_evidence_probability_kokkos(
                                tlo_alt * std::log(10.0), total_depth, alt_depth,
                                log_somatic_prior, empirical_somatic_model,
                                &kernel_telemetry);
                        }
                        const auto pger_info = info_probability_at(
                            input_header, record, "PGERMLINE", alt_index);
                        double germline_error = std::numeric_limits<double>::quiet_NaN();
                        const auto popaf_values = info_float_values(
                            input_header, record, "POPAF");
                        if (!popaf_values.empty() &&
                            alt_index - 1 < static_cast<int>(popaf_values.size()) &&
                            std::isfinite(static_cast<double>(popaf_values[alt_index - 1]))) {
                            const auto population_af = std::pow(10.0,
                                -static_cast<double>(popaf_values[alt_index - 1]));
                            const auto weighted_tumor_af = tumor_allele_fraction_for_alt(
                                input_header, record, 0, alt_index);
                            germline_error = germline_probability_from_counts_kokkos(
                                0.0, total_depth, alt_depth, weighted_tumor_af, 0.5,
                                population_af, log_somatic_prior,
                                empirical_somatic_model, &kernel_telemetry);
                        }
                        if (!std::isfinite(germline_error) && std::isfinite(pger_info))
                            germline_error = pger_info;
                        const auto contamination_iterator =
                            contamination_probabilities.find(alt_index);
                        const auto contamination_value =
                            contamination_iterator != contamination_probabilities.end()
                            ? contamination_iterator->second
                            : std::numeric_limits<double>::quiet_NaN();
                        double non_somatic_error = std::numeric_limits<double>::quiet_NaN();
                        if (std::isfinite(germline_error))
                            non_somatic_error = germline_error;
                        if (std::isfinite(contamination_value))
                            non_somatic_error = std::isfinite(non_somatic_error)
                                ? std::max(non_somatic_error, contamination_value)
                                : contamination_value;
                        const bool keep_datum = total_depth > 0 && alt_depth >= 0 &&
                            !(std::isfinite(non_somatic_error) && non_somatic_error > 0.9) &&
                            !(std::isfinite(artifact_probability) && artifact_probability > 0.9);
                        if (keep_datum) {
                            const auto artifact_keep = std::isfinite(artifact_probability)
                                ? artifact_probability : 0.0;
                            const auto non_somatic_keep = std::isfinite(non_somatic_error)
                                ? non_somatic_error : 0.0;
                            const auto sequencing_keep = std::isfinite(sequencing_error)
                                ? sequencing_error : 0.0;
                            const auto somatic_weight = std::clamp(
                                (1.0 - artifact_keep) * (1.0 - non_somatic_keep)
                                * (1.0 - sequencing_keep), 0.0, 1.0);
                            empirical_model_data.push_back(SomaticModelDatum{
                                total_depth, alt_depth,
                                indel_length_for_alt(record, alt_index),
                                somatic_weight});
                        }
                        } else {
                        const auto posterior = info_probability_at(
                            input_header, record, "PSOMATIC", alt_index);
                        const auto artifact = info_probability_at(
                            input_header, record, "PARTIFACT", alt_index);
                        const auto orientation = info_probability_at(
                            input_header, record, "OBP", alt_index);
                        const auto non_somatic = info_probability_at(
                            input_header, record, "PGERMLINE", alt_index);
                        const auto artifact_probability = std::isfinite(artifact) && std::isfinite(orientation)
                            ? std::max(artifact, orientation)
                            : (std::isfinite(artifact) ? artifact : orientation);
                        const auto somatic_weight = std::isfinite(posterior)
                            ? posterior
                            : (std::isfinite(non_somatic) || std::isfinite(artifact_probability)
                                ? (1.0 - std::max(
                                    std::isfinite(non_somatic) ? non_somatic : 0.0,
                                    std::isfinite(artifact_probability) ? artifact_probability : 0.0))
                                : 1.0);
                        if (total_depth > 0 && alt_depth >= 0)
                            empirical_model_data.push_back(SomaticModelDatum{
                                total_depth, alt_depth,
                                indel_length_for_alt(record, alt_index),
                                std::clamp(somatic_weight, 0.0, 1.0)});
                    }
                        }
                    if (const char* debug_gate = std::getenv("FASTGATK_DEBUG_FMC")) {
                        if (std::string(debug_gate) == "1") {
                            int total_depth = 0;
                            if (tumor_depths.size() ==
                                static_cast<std::size_t>(record->n_allele))
                                for (const auto depth : tumor_depths)
                                    if (depth > 0) total_depth += depth;
                            const int alt_depth =
                                (tumor_depths.size() ==
                                 static_cast<std::size_t>(record->n_allele))
                                ? tumor_depths[static_cast<std::size_t>(alt_index)]
                                : -1;
                            const auto contam_it = contamination_probabilities.find(alt_index);
                            const auto contam_value = contam_it != contamination_probabilities.end()
                                ? contam_it->second : std::numeric_limits<double>::quiet_NaN();
                            double germ_candidate = std::numeric_limits<double>::quiet_NaN();
                            const auto popaf_values = info_float_values(input_header, record, "POPAF");
                            if (!popaf_values.empty() &&
                                alt_index - 1 < static_cast<int>(popaf_values.size()) &&
                                std::isfinite(static_cast<double>(popaf_values[alt_index - 1]))) {
                                const auto population_af = std::pow(10.0,
                                    -static_cast<double>(popaf_values[alt_index - 1]));
                                const auto weighted_tumor_af = tumor_allele_fraction_for_alt(
                                    input_header, record, 0, alt_index);
                                const auto reference_length = record->d.allele[0] == nullptr
                                    ? 0 : static_cast<int>(std::strlen(record->d.allele[0]));
                                const auto alternate_length = record->d.allele[alt_index] == nullptr
                                    ? reference_length
                                    : static_cast<int>(std::strlen(record->d.allele[alt_index]));
                                const auto indel_length = alternate_length - reference_length;
                                double log_somatic_prior = -6.0 * std::log(10.0) - std::log(3.0);
                                if (indel_length >= -10 && indel_length <= 10 &&
                                    std::isfinite(empirical_somatic_model.log_variant_priors[
                                        static_cast<std::size_t>(indel_length + 10)])) {
                                    log_somatic_prior = empirical_somatic_model.log_variant_priors[
                                        static_cast<std::size_t>(indel_length + 10)] +
                                        (indel_length == 0 ? -std::log(3.0) : 0.0);
                                }
                                germ_candidate = germline_probability_from_counts_kokkos(
                                    0.0, total_depth, alt_depth, weighted_tumor_af, 0.5,
                                    population_af, log_somatic_prior,
                                    empirical_somatic_model, &kernel_telemetry);
                            }
                            std::cerr << "[FASTGATK_DEBUG_FMC] pos=" << record->pos + 1
                                << " alt=" << alt_index << " depth=" << total_depth
                                << " alt_depth=" << alt_depth
                                << " contam=" << contam_value
                                << " germ_recompute=" << germ_candidate
                                << " pger=" << info_probability_at(
                                    input_header, record, "PGERMLINE", alt_index)
                                << " psom=" << info_probability_at(
                                    input_header, record, "PSOMATIC", alt_index) << '\n';
                        }
                    }
                    const auto probability = combined_error_probability_for_alt_kokkos(
                        input_header, record, alt_index, tumor_depths, tlod_values,
                        &empirical_somatic_model, &kernel_telemetry,
                        &contamination_probabilities, orientation_learning_posterior);
                    if (std::isfinite(probability) && probability >= 0.0) {
                        learned_joint_errors.push_back(std::clamp(probability, 0.0, 1.0));
                    }
                }
                bcf_clear(record);
            }
                const bool learn_after = java_passes ? (learning_pass < 2) : true;
                if (learn_after) {
                    const auto log_snv_prior = empirical_somatic_model.log_snv_prior;
            const auto log_indel_prior = empirical_somatic_model.log_indel_prior;
            const auto log_artifact_prior = empirical_somatic_model.log_variant_vs_artifact_prior;
            const auto initial_prior_map = empirical_somatic_model.log_variant_priors;
            const auto callable_sites = empirical_somatic_model.callable_sites;
            empirical_somatic_model = learn_empirical_somatic_model(
                empirical_model_data, callable_sites);
            empirical_somatic_model.log_snv_prior = log_snv_prior;
            empirical_somatic_model.log_indel_prior = log_indel_prior;
            empirical_somatic_model.log_variant_vs_artifact_prior = log_artifact_prior;
            empirical_somatic_model.callable_sites = callable_sites;
            if (!std::isfinite(callable_sites) || callable_sites < 1.0) {
                empirical_somatic_model.log_variant_priors = initial_prior_map;
            } else {
                empirical_somatic_model.log_snv_prior = empirical_somatic_model.log_variant_priors[10];
                empirical_somatic_model.log_indel_prior = empirical_somatic_model.log_variant_priors[9];
            }
                    empirical_somatic_model_records = empirical_model_data.size();
                }
                if (learning_pass == learning_total - 1)
                    threshold_errors = learned_joint_errors;
                joint_error_threshold_observations = threshold_errors.size();
            if (!threshold_errors.empty()) {
                if (options.joint_threshold_strategy == "CONSTANT")
                    effective_joint_error_threshold = options.initial_threshold;
                else if (options.joint_threshold_strategy == "FALSE_DISCOVERY_RATE")
                    effective_joint_error_threshold = calculate_false_discovery_rate_threshold(
                        threshold_errors, options.max_false_discovery_rate);
                else
                    effective_joint_error_threshold = calculate_optimal_f_score_threshold(
                        threshold_errors, options.f_score_beta);
            }
            bcf_close(input);
            input = nullptr;
            bcf_hdr_destroy(input_header);
            input_header = nullptr;
            input = bcf_open(options.input.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen Mutect2 VCF for ErrorProbabilities learning");
            input_header = bcf_hdr_read(input);
            if (!input_header) throw std::runtime_error("BAD_INPUT: cannot reread Mutect2 VCF header");
            }
            }
        }
        while (fastgatk::io::read_variant_record(input, input_header, record,
                                               options.input) == 0) {
            if (!options.regions.empty() &&
                !filter_record_in_intervals(input_header, record, intervals)) {
                ++interval_skipped_records;
                bcf_clear(record);
                continue;
            }
            ++input_records;
            bcf_unpack(record, BCF_UN_ALL);
            const auto tlod_values = info_float_values(input_header, record, "TLOD");
            const auto nalod_values = info_float_values(input_header, record, "NALOD");
            const int concrete_alts = concrete_alt_count(record);
            // Mutect2 emits most evidence annotations as Number=A.  Read all
            // concrete ALTs before deriving site-level filters so a strong
            // second ALT cannot be hidden by the first vector element.  For
            // legacy scalar annotations these vectors simply broadcast the
            // scalar to every concrete ALT.
            const auto tlod_by_alt = info_probability_values_by_alt(
                input_header, record, "TLOD");
            const auto somatic_probability_by_alt = info_probability_values_by_alt(
                input_header, record, "PSOMATIC");
            const auto germline_probability_by_alt = info_probability_values_by_alt(
                input_header, record, "PGERMLINE");
            const auto artifact_probability_by_alt = info_probability_values_by_alt(
                input_header, record, "PARTIFACT");
            const auto tlod = static_cast<float>(aggregate_probability_values(
                tlod_by_alt, true, false, info_float(input_header, record, "TLOD", 0.0F)));
            const auto somatic_probability = static_cast<float>(aggregate_probability_values(
                somatic_probability_by_alt, false, true,
                std::numeric_limits<float>::quiet_NaN()));
            const auto germline_probability = static_cast<float>(aggregate_probability_values(
                germline_probability_by_alt, true, true,
                std::numeric_limits<float>::quiet_NaN()));
            const auto artifact_probability = static_cast<float>(aggregate_probability_values(
                artifact_probability_by_alt, true, true,
                std::numeric_limits<float>::quiet_NaN()));
            const int orientation_alt_index = max_tlod_alt_index(input_header, record);
            const auto scalar_f1r2 = format_int_value(
                input_header, record, "F1R2", tumor_sample_index, orientation_alt_index - 1);
            const auto scalar_f2r1 = format_int_value(
                input_header, record, "F2R1", tumor_sample_index, orientation_alt_index - 1);
            const auto info_f1r2 = info_int(input_header, record, "F1R2", -1);
            const auto info_f2r1 = info_int(input_header, record, "R1F2", -1);
            const auto allele_specific_strand = info_alt_strand_counts(
                input_header, record, "AS_SB_TABLE", orientation_alt_index - 1);
            const auto orientation_counts = read_orientation_counts(
                input_header, record, tumor_sample_index, orientation_alt_index);
            const auto f1r2 = orientation_counts.f1r2;
            const auto r1f2 = orientation_counts.f2r1;
            const auto median_base_quality = info_numeric(input_header, record, "MBQ", -1.0F);
            const auto median_mapping_quality = info_numeric(input_header, record, "MMQ", -1.0F);
            const auto median_read_position = info_numeric(input_header, record, "MPOS", -1.0F);
            const auto median_fragment_lengths = info_int_values(input_header, record, "MFRL");
            const auto n_count_value = info_first_int_optional(input_header, record, "NCount");
            const auto n_count = n_count_value.value_or(0);
            const auto region_event_count = info_int(input_header, record, "ECNT", -1);
            const auto haplotype_event_count = info_int_max(input_header, record, "ECNTH");
            const bool germline = bcf_has_filter(input_header, record, const_cast<char*>("germline")) > 0;
            const auto tumor_af_by_alt = [&]() {
                std::vector<double> values;
                values.reserve(static_cast<std::size_t>(concrete_alts));
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    const auto* allele = record->d.allele[alt_index];
                    if (allele == nullptr || allele[0] == '<') continue;
                    values.push_back(tumor_allele_fraction_for_alt(
                        input_header, record, tumor_sample_index, alt_index));
                }
                return values;
            }();
            const auto has_failing_or_missing = [](const std::vector<double>& values,
                                                   const double threshold) {
                return !values.empty() && std::all_of(values.begin(), values.end(),
                    [&](const double value) {
                        return !std::isfinite(value) || value < threshold;
                    });
            };
            const bool low_tlod = has_failing_or_missing(tlod_by_alt, options.min_tlod) ||
                (tlod_by_alt.empty() && tlod < options.min_tlod);
            std::map<int, double> contamination_probability_by_alt;
            bool contamination_probability_evaluated = false;
            if (effective_contamination >= 0.0F) {
                contamination_probability_by_alt =
                    contamination_probabilities_for_record_kokkos(
                        input_header, record, effective_contamination,
                        contamination_by_sample, tumor_sample_indices,
                        empirical_somatic_model, &kernel_telemetry,
                        &contamination_probability_evaluated);
                if (contamination_probability_evaluated) {
                    ++contamination_posterior_records;
                    contamination_posterior_alleles += contamination_probability_by_alt.size();
                }
            }
            bool contamination_site = false;
            if (contamination_probability_evaluated) {
                int observed_alts = 0;
                int failing_alts = 0;
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    const auto* allele = record->d.allele[alt_index];
                    if (allele == nullptr || allele[0] == '<') continue;
                    ++observed_alts;
                    const auto iterator = contamination_probability_by_alt.find(alt_index);
                    const auto probability = iterator == contamination_probability_by_alt.end()
                        ? std::numeric_limits<double>::quiet_NaN() : iterator->second;
                    if (std::isfinite(probability) &&
                        probability > options.max_contamination_probability) {
                        ++failing_alts;
                    }
                }
                contamination_site = observed_alts > 0 && failing_alts == observed_alts;
            }
            const bool contamination = contamination_probability_evaluated
                ? contamination_site
                : (effective_contamination >= 0.0F &&
                   has_failing_or_missing(tumor_af_by_alt, effective_contamination));
            const int orientation_total = f1r2 >= 0 && r1f2 >= 0 ? f1r2 + r1f2 : 0;
            const int orientation_minor = orientation_total > 0 ? std::min(f1r2, r1f2) : 0;
            const bool orientation = options.min_orientation_balance >= 0.0F && orientation_total > 0 &&
                                     static_cast<float>(orientation_minor) / static_cast<float>(orientation_total) <
                                     options.min_orientation_balance;
            const bool low_af = options.min_allele_fraction > 0.0F &&
                has_failing_or_missing(tumor_af_by_alt, options.min_allele_fraction);
            const bool strand = options.min_reads_per_strand > 0 &&
                                (f1r2 < options.min_reads_per_strand ||
                                 r1f2 < options.min_reads_per_strand);
            const auto alt_depth = tumor_alt_depth(input_header, record, tumor_sample_index);
            if (orientation_total > 0 && alt_depth >= 0 && orientation_total != alt_depth)
                ++orientation_depth_disagreements;
            const auto allele_specific_unique_alt = info_alt_int(
                input_header, record, "AS_UNIQ_ALT_READ_COUNT", orientation_alt_index - 1);
            const bool used_allele_specific_strand = allele_specific_strand &&
                ((scalar_f1r2 < 0 && info_f1r2 < 0) || (scalar_f2r1 < 0 && info_f2r1 < 0));
            if (used_allele_specific_strand) ++allele_specific_strand_records;
            if (allele_specific_unique_alt) ++allele_specific_unique_records;
            const int unique_alt_depth = allele_specific_unique_alt.value_or(alt_depth);
            const bool low_alt_reads = options.unique_alt_read_count > 0 &&
                                       (unique_alt_depth < 0 || unique_alt_depth <= options.unique_alt_read_count);
            int tlod_passing_alts = 0;
            if (static_cast<int>(tlod_values.size()) >= concrete_alts) {
                for (int index = 0; index < concrete_alts; ++index)
                    if (static_cast<double>(tlod_values[static_cast<std::size_t>(index)]) * std::log(10.0) > 5.0)
                        ++tlod_passing_alts;
            }
            // Older native Mutect2 records sometimes carry a scalar TLOD on a
            // multi-ALT record.  Preserve the historical concrete-ALT guard
            // for that malformed-but-common representation; canonical
            // Number=A TLOD uses Rust/GATK's natural-log threshold above.
            const bool mode_excludes_genomic_filters =
                options.mitochondria_mode || options.microbial_mode;
            const bool multiallelic = !mode_excludes_genomic_filters &&
                concrete_alts > options.max_alt_allele_count &&
                (static_cast<int>(tlod_values.size()) < concrete_alts ||
                 tlod_passing_alts > options.max_alt_allele_count);
            const bool low_median_base_quality =
                options.min_median_base_quality >= 0 && median_base_quality >= 0.0F &&
                median_base_quality < options.min_median_base_quality;
            const bool low_median_mapping_quality =
                options.min_median_mapping_quality >= 0 && median_mapping_quality >= 0.0F &&
                median_mapping_quality < options.min_median_mapping_quality;
            const bool low_median_read_position =
                options.min_median_read_position >= 0 && median_read_position >= 0.0F &&
                median_read_position < options.min_median_read_position;
            const bool fragment_length = !mode_excludes_genomic_filters &&
                options.max_median_fragment_length_difference >= 0 &&
                median_fragment_lengths.size() >= 2 &&
                std::abs(median_fragment_lengths[1] - median_fragment_lengths[0]) >
                    options.max_median_fragment_length_difference;
            const int alt_depth_total = total_alt_depth_all_samples(input_header, record);
            // Java NRatioFilter skips records with no NCount annotation or
            // zero alternate depth; only a positive ALT denominator reaches
            // the ratio comparison. This also makes the GATK default
            // +Infinity a true no-op without manufacturing a filter on an
            // annotation-only record.
            const bool n_ratio = options.max_n_ratio >= 0.0 && n_count_value &&
                alt_depth_total > 0 &&
                static_cast<double>(n_count) / static_cast<double>(alt_depth_total) >=
                    options.max_n_ratio;
            const bool panel_of_normals = info_present(input_header, record, "PON");
            // The first deterministic TumorEvidenceFilter boundary uses the
            // unlearned SomaticClusteringModel.  Keep one probability per
            // concrete ALT so AS_FilterStatus remains allele-specific.
            const auto tumor_depths_for_filters = tumor_allele_depths(
                input_header, record, tumor_sample_indices);
            // GATK's GermlineFilter is a variant-level posterior.  Reproduce
            // its common case (no tumor-segmentation MAF override) from the
            // max-TLOD ALT, tumour AD, POPAF and optional NLOD.  PGERMLINE is
            // not substituted here: it is an upstream Mutect2 annotation,
            // whereas GATK computes GERMQ from the model below.
            double germline_filter_probability = std::numeric_limits<double>::quiet_NaN();
            if (!tlod_values.empty() &&
                bcf_hdr_id2int(input_header, BCF_DT_ID, "POPAF") >= 0 &&
                tumor_depths_for_filters.size() == static_cast<std::size_t>(record->n_allele)) {
                const int germline_alt_index = max_tlod_alt_index(input_header, record);
                const auto population_negative_log10_af = info_probability_at(
                    input_header, record, "POPAF", germline_alt_index);
                const auto population_af = std::isfinite(population_negative_log10_af)
                    ? std::pow(10.0, -population_negative_log10_af)
                    : std::numeric_limits<double>::quiet_NaN();
                int total_tumor_count = 0;
                for (const auto depth : tumor_depths_for_filters)
                    if (depth > 0) total_tumor_count += depth;
                const int alt_tumor_count = germline_alt_index > 0 &&
                    germline_alt_index < static_cast<int>(tumor_depths_for_filters.size())
                    ? tumor_depths_for_filters[static_cast<std::size_t>(germline_alt_index)] : -1;
                const auto weighted_tumor_af = tumor_allele_fraction_for_alt(
                    input_header, record, tumor_sample_index, germline_alt_index);
                double normal_log10_odds = 0.0;
                const auto normal_lod_values = info_float_values_indexed(
                    input_header, record, "NLOD");
                if (germline_alt_index > 0 &&
                    germline_alt_index <= static_cast<int>(normal_lod_values.size()) &&
                    std::isfinite(normal_lod_values[static_cast<std::size_t>(germline_alt_index - 1)]))
                    normal_log10_odds = normal_lod_values[
                        static_cast<std::size_t>(germline_alt_index - 1)];
                const auto delta = indel_length_for_alt(record, germline_alt_index);
                double germline_log_prior = delta == 0
                    ? empirical_somatic_model.log_snv_prior - std::log(3.0)
                    : empirical_somatic_model.log_indel_prior;
                if (empirical_somatic_model.learned && delta >= -10 && delta <= 10)
                    germline_log_prior = empirical_somatic_model.log_variant_priors[
                        static_cast<std::size_t>(delta + 10)] +
                        (delta == 0 ? -std::log(3.0) : 0.0);
                if (total_tumor_count > 0 && alt_tumor_count >= 0 &&
                    std::isfinite(population_af) && population_af >= 0.0 && population_af <= 1.0) {
                    germline_filter_probability = germline_probability_from_counts_kokkos(
                        normal_log10_odds, total_tumor_count, alt_tumor_count,
                        weighted_tumor_af, 0.5, population_af, germline_log_prior,
                        empirical_somatic_model, &kernel_telemetry);
                }
            }
            // GermlineFilter is a variant-level ErrorProbabilities source in
            // GATK, not merely a GERMQ writer.  Once its posterior exceeds
            // the effective threshold it must contribute the site-level
            // germline FILTER as well.  The previous adapter materialized
            // GERMQ but silently omitted this filtering decision, which was
            // observable for high-population-AF calls in the contamination
            // table oracle.
            const auto germline_filter_threshold = effective_joint_error_threshold >= 0.0
                ? effective_joint_error_threshold : options.filter_error_probability_threshold;
            const bool germline_filter = std::isfinite(germline_filter_probability) &&
                germline_filter_probability > germline_filter_threshold;
            std::vector<double> weak_evidence_probabilities;
            weak_evidence_probabilities.reserve(static_cast<std::size_t>(concrete_alts));
            // Keep the probability keyed by original ALT index as well as in
            // the compact vector used for site-level reduction.  Symbolic
            // alleles (for example <NON_REF>) can be interleaved with concrete
            // ALTs; using the compact offset for FilteringOutputStats would
            // then attach the wrong weak-evidence value to a concrete ALT.
            std::map<int, double> weak_evidence_by_alt;
            std::vector<double> weak_batch_log_odds;
            std::vector<int> weak_batch_total_counts;
            std::vector<int> weak_batch_alt_counts;
            std::vector<double> weak_batch_priors;
            std::vector<int> weak_batch_alt_indices;
            std::vector<std::size_t> weak_batch_offsets;
            if (tumor_depths_for_filters.size() == static_cast<std::size_t>(record->n_allele) &&
                !tlod_values.empty()) {
                int total_tumor_depth = 0;
                for (const auto depth : tumor_depths_for_filters)
                    if (depth > 0) total_tumor_depth += depth;
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    const auto* allele = record->d.allele[alt_index];
                    if (allele == nullptr || allele[0] == '<' ||
                        alt_index > static_cast<int>(tlod_values.size())) {
                        if (allele != nullptr && allele[0] != '<') {
                            weak_evidence_probabilities.push_back(0.0);
                            weak_evidence_by_alt[alt_index] = 0.0;
                        }
                        continue;
                    }
                    const int alt_depth_for_filter = tumor_depths_for_filters[
                        static_cast<std::size_t>(alt_index)];
                    const auto log_odds = static_cast<double>(
                        tlod_values[static_cast<std::size_t>(alt_index - 1)]) * std::log(10.0);
                    const auto reference_length = record->d.allele[0] == nullptr
                        ? 0 : static_cast<int>(std::strlen(record->d.allele[0]));
                    const auto alternate_length = record->d.allele[alt_index] == nullptr
                        ? reference_length : static_cast<int>(std::strlen(record->d.allele[alt_index]));
                    // Preserve the historical native initial-model boundary
                    // when no empirical pass was requested; once the model
                    // is learned, use GATK's SNV/indel prior map.
                    double log_somatic_prior = empirical_somatic_model.learned
                        ? (alternate_length == reference_length
                            ? empirical_somatic_model.log_snv_prior - std::log(3.0)
                            : empirical_somatic_model.log_indel_prior)
                        : -std::log(10.0);
                    const auto indel_length = alternate_length - reference_length;
                    if (empirical_somatic_model.learned && indel_length >= -10 && indel_length <= 10)
                        log_somatic_prior = empirical_somatic_model.log_variant_priors[
                            static_cast<std::size_t>(indel_length + 10)] +
                            (indel_length == 0 ? -std::log(3.0) : 0.0);
                    weak_batch_alt_indices.push_back(alt_index);
                    weak_batch_offsets.push_back(weak_evidence_probabilities.size());
                    weak_evidence_probabilities.push_back(0.0);
                    weak_batch_log_odds.push_back(log_odds);
                    weak_batch_total_counts.push_back(total_tumor_depth);
                    weak_batch_alt_counts.push_back(alt_depth_for_filter);
                    weak_batch_priors.push_back(log_somatic_prior);
                }
            }
            if (!weak_batch_log_odds.empty()) {
                std::vector<double> probabilities;
                if (weak_batch_log_odds.size() == 1) {
                    // A scalar ALT is the overwhelmingly common case in
                    // ordinary Mutect2 VCFs. Avoid paying batch packing cost
                    // when there is no parallelism to expose; multi-ALT
                    // records use the one-launch path below.
                    probabilities.push_back(empirical_somatic_model.learned
                        ? empirical_tumor_evidence_probability_kokkos(
                            weak_batch_log_odds[0], weak_batch_total_counts[0],
                            weak_batch_alt_counts[0], weak_batch_priors[0],
                            empirical_somatic_model, &kernel_telemetry)
                        : tumor_evidence_probability_with_prior_kokkos(
                            weak_batch_log_odds[0], weak_batch_total_counts[0],
                            weak_batch_alt_counts[0], weak_batch_priors[0],
                            &kernel_telemetry));
                } else {
                    probabilities = tumor_evidence_probability_batch_kokkos(
                        weak_batch_log_odds, weak_batch_total_counts, weak_batch_alt_counts,
                        weak_batch_priors,
                        empirical_somatic_model.learned ? &empirical_somatic_model : nullptr,
                        &kernel_telemetry);
                }
                for (std::size_t index = 0; index < probabilities.size(); ++index) {
                    const auto offset = weak_batch_offsets[index];
                    weak_evidence_probabilities[offset] = probabilities[index];
                    weak_evidence_by_alt[weak_batch_alt_indices[index]] = probabilities[index];
                }
            }
            // Match Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats:
            // the learned threshold is bounded only to avoid finite-precision
            // 0/1 decisions.  It is shared by every default posterior filter,
            // including weak evidence and strand artifact.
            constexpr double kFilterEpsilon = 1.0e-10;
            const double gatk_error_threshold = std::min(1.0 - kFilterEpsilon,
                std::max(kFilterEpsilon, effective_joint_error_threshold >= 0.0
                    ? effective_joint_error_threshold
                    : options.filter_error_probability_threshold));
            const bool weak_evidence = !weak_evidence_probabilities.empty() &&
                std::any_of(weak_evidence_probabilities.begin(), weak_evidence_probabilities.end(),
                            [&](const double probability) {
                                return std::isfinite(probability) &&
                                       probability > gatk_error_threshold;
                            });
            const bool weak_evidence_site = !weak_evidence_probabilities.empty() &&
                std::all_of(weak_evidence_probabilities.begin(), weak_evidence_probabilities.end(),
                            [&](const double probability) {
                                return std::isfinite(probability) &&
                                       probability > gatk_error_threshold;
                            });
            std::map<int, double> strand_artifact_probability_by_alt;
            const auto strand_counts_by_allele = info_allele_strand_counts(
                input_header, record, "AS_SB_TABLE");
            if (strand_counts_by_allele &&
                strand_counts_by_allele->size() == static_cast<std::size_t>(record->n_allele)) {
                int total_forward = 0;
                int total_reverse = 0;
                for (const auto& counts : *strand_counts_by_allele) {
                    total_forward += counts.first;
                    total_reverse += counts.second;
                }
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    if (record->d.allele[alt_index] == nullptr ||
                        record->d.allele[alt_index][0] == '<') continue;
                    const auto& counts = (*strand_counts_by_allele)[
                        static_cast<std::size_t>(alt_index)];
                    strand_artifact_probability_by_alt.emplace(alt_index,
                        strand_artifact_filter.error_probability(
                            total_forward, total_reverse,
                            counts.first, counts.second,
                            std::abs(indel_length_for_alt(record, alt_index)),
                            &kernel_telemetry));
                }
            }
            const auto strand_artifact_applies = [&](const int alt_index) {
                const auto iterator = strand_artifact_probability_by_alt.find(alt_index);
                return iterator != strand_artifact_probability_by_alt.end() &&
                    std::isfinite(iterator->second) && iterator->second > gatk_error_threshold;
            };
            const bool strand_artifact_site = !strand_artifact_probability_by_alt.empty() &&
                std::all_of(strand_artifact_probability_by_alt.begin(),
                            strand_artifact_probability_by_alt.end(),
                            [&](const auto& entry) { return strand_artifact_applies(entry.first); });
            // ReadOrientationFilter is a variant-level ARTIFACT contributor
            // to ErrorProbabilities.  Evaluate it before the joint posterior
            // so both final filtering and the learned threshold see it.
            const auto orientation_evaluation = evaluate_orientation_posterior(
                input_header, record, reference_index, orientation_collections,
                orientation_sample_indices, orientation_prior_device, orientation_log_device,
                &kernel_telemetry);
            const auto orientation_artifact_probability_value = orientation_evaluation.probability;
            const auto has_orientation_prior = orientation_evaluation.has_prior;
            if (has_orientation_prior) {
                ++orientation_prior_records;
                orientation_prior_sample_evaluations += orientation_evaluation.sample_evaluations;
            }
            const bool orientation_prior_filter = has_orientation_prior &&
                effective_orientation_threshold >= 0.0F &&
                orientation_artifact_probability_value > effective_orientation_threshold;
            std::vector<double> joint_error_probabilities;
            std::map<int, double> joint_error_by_alt;
            if (effective_joint_error_threshold >= 0.0) {
                joint_error_probabilities.reserve(static_cast<std::size_t>(concrete_alts));
                for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                    if (record->d.allele[alt_index] == nullptr || record->d.allele[alt_index][0] == '<')
                        continue;
                    const auto probability = combined_error_probability_for_alt_kokkos(
                        input_header, record, alt_index, tumor_depths_for_filters, tlod_values,
                        &empirical_somatic_model, &kernel_telemetry,
                        &contamination_probability_by_alt,
                        has_orientation_prior ? orientation_artifact_probability_value
                                              : std::numeric_limits<double>::quiet_NaN());
                    joint_error_probabilities.push_back(probability);
                    joint_error_by_alt[alt_index] = probability;
                }
            }
            bool joint_error = !joint_error_probabilities.empty() &&
                std::any_of(joint_error_probabilities.begin(), joint_error_probabilities.end(),
                            [&](const double probability) {
                                return std::isfinite(probability) &&
                                       probability > effective_joint_error_threshold;
                            });
            bool joint_error_site = !joint_error_probabilities.empty() &&
                std::all_of(joint_error_probabilities.begin(), joint_error_probabilities.end(),
                            [&](const double probability) {
                                return std::isfinite(probability) &&
                                       probability > effective_joint_error_threshold;
                            });
            double slippage_probability = 0.0;
            bool slippage = false;
            const auto repeats_per_allele = info_int_values(input_header, record, "RPA");
            const auto repeat_unit = info_string_value(input_header, record, "RU");
            if (!options.mitochondria_mode && repeats_per_allele.size() >= 2 && repeat_unit &&
                tumor_depths_for_filters.size() == static_cast<std::size_t>(record->n_allele) &&
                record->n_allele > 1 && record->d.allele[0] != nullptr &&
                record->d.allele[1] != nullptr) {
                int total_tumor_depth = 0;
                for (const auto depth : tumor_depths_for_filters)
                    if (depth > 0) total_tumor_depth += depth;
                const auto reference_length = static_cast<int>(std::strlen(record->d.allele[0]));
                const auto alternate_length = static_cast<int>(std::strlen(record->d.allele[1]));
                double log_somatic_prior = alternate_length == reference_length
                    ? empirical_somatic_model.log_snv_prior - std::log(3.0)
                    : empirical_somatic_model.log_indel_prior;
                const auto indel_length = alternate_length - reference_length;
                if (empirical_somatic_model.learned && indel_length >= -10 && indel_length <= 10)
                    log_somatic_prior = empirical_somatic_model.log_variant_priors[
                        static_cast<std::size_t>(indel_length + 10)] +
                        (indel_length == 0 ? -std::log(3.0) : 0.0);
                if (empirical_somatic_model.learned)
                    slippage_probability = empirical_slippage_probability_kokkos(
                        empirical_somatic_model, options.slippage_rate,
                        repeats_per_allele[0], repeats_per_allele[1],
                        static_cast<int>(repeat_unit->size()), options.min_slippage_length,
                        total_tumor_depth, tumor_depths_for_filters[1],
                        reference_length, alternate_length, log_somatic_prior, &kernel_telemetry);
                else
                    slippage_probability = slippage_probability_with_prior_kokkos(
                        options.slippage_rate, repeats_per_allele[0], repeats_per_allele[1],
                        static_cast<int>(repeat_unit->size()), options.min_slippage_length,
                        total_tumor_depth, tumor_depths_for_filters[1],
                        reference_length, alternate_length, log_somatic_prior, &kernel_telemetry);
                slippage = std::isfinite(slippage_probability) &&
                           slippage_probability > options.filter_error_probability_threshold;
            }
            double normal_artifact_probability = 0.0;
            bool normal_artifact = false;
            const auto mbq_values_for_normal = info_int_values(input_header, record, "MBQ");
            if (!normal_sample_indices.empty() && !tlod_values.empty() &&
                !nalod_values.empty() && !mbq_values_for_normal.empty()) {
                const int normal_alt_index = max_tlod_alt_index(input_header, record);
                if (normal_alt_index > 0 &&
                    normal_alt_index <= static_cast<int>(nalod_values.size())) {
                    int tumor_depth_for_normal_artifact = 0;
                    int tumor_alt_for_normal_artifact = 0;
                    for (const auto sample_index : tumor_sample_indices) {
                        const auto reference_depth = format_int_value(
                            input_header, record, "AD", sample_index, 0);
                        const auto alternate_depth = format_int_value(
                            input_header, record, "AD", sample_index, normal_alt_index);
                        const auto alternate_total = format_sum_alt_depth(
                            input_header, record, sample_index);
                        if (reference_depth >= 0 && alternate_total >= 0)
                            tumor_depth_for_normal_artifact += reference_depth + alternate_total;
                        if (alternate_depth >= 0)
                            tumor_alt_for_normal_artifact += alternate_depth;
                    }
                    int normal_depth_for_normal_artifact = 0;
                    int normal_alt_for_normal_artifact = 0;
                    for (const auto sample_index : normal_sample_indices) {
                        const auto reference_depth = format_int_value(
                            input_header, record, "AD", sample_index, 0);
                        const auto alternate_depth = format_int_value(
                            input_header, record, "AD", sample_index, normal_alt_index);
                        const auto alternate_total = format_sum_alt_depth(
                            input_header, record, sample_index);
                        if (reference_depth >= 0 && alternate_total >= 0)
                            normal_depth_for_normal_artifact += reference_depth + alternate_total;
                        if (alternate_depth >= 0)
                            normal_alt_for_normal_artifact += alternate_depth;
                    }
                    normal_artifact_probability = normal_artifact_probability_with_prior_kokkos(
                        nalod_values[static_cast<std::size_t>(normal_alt_index - 1)],
                        tumor_depth_for_normal_artifact, tumor_alt_for_normal_artifact,
                        normal_depth_for_normal_artifact, normal_alt_for_normal_artifact,
                        mbq_values_for_normal.front(), options.normal_pileup_p_value_threshold,
                        empirical_somatic_model.learned
                            ? empirical_somatic_model.log_variant_vs_artifact_prior
                            : options.log_artifact_prior,
                        &kernel_telemetry);
                    // NormalArtifactFilter returns an ARTIFACT error
                    // probability.  GATK compares it with the engine's
                    // effective error threshold (the constant/learned joint
                    // threshold when configured, otherwise the native
                    // compatibility threshold); a merely positive posterior
                    // is not sufficient to add normal_artifact.
                    const auto normal_artifact_threshold =
                        effective_joint_error_threshold >= 0.0
                        ? effective_joint_error_threshold
                        : options.filter_error_probability_threshold;
                    normal_artifact = std::isfinite(normal_artifact_probability) &&
                                      normal_artifact_probability > normal_artifact_threshold;
                }
            }
            const bool low_somatic_probability = options.min_somatic_probability > 0.0F &&
                (!std::isfinite(somatic_probability) || somatic_probability < options.min_somatic_probability);
            const bool high_germline_probability = options.max_germline_probability >= 0.0F &&
                (!std::isfinite(germline_probability) || germline_probability > options.max_germline_probability);
            const bool high_artifact_probability = options.max_artifact_probability >= 0.0F &&
                (!std::isfinite(artifact_probability) || artifact_probability > options.max_artifact_probability);
            const bool clustered_events = !mode_excludes_genomic_filters &&
                ((region_event_count >= 0 && region_event_count > options.max_events_in_region) ||
                 (haplotype_event_count && *haplotype_event_count > options.max_events_in_haplotype));
            const auto phased_observations = has_phasing_fields
                ? phased_haplotype_observations(input_header, record, haplotype_sample_indices)
                : std::vector<PhasedHaplotypeObservation>{};
            const auto best_haplotype = std::max_element(
                phased_observations.begin(), phased_observations.end(),
                [](const auto& left, const auto& right) {
                    // std::max_element keeps the first element on ties,
                    // matching Java Stream.max(maxBy)'s left-tie behavior.
                    return left.allele_fraction < right.allele_fraction;
                });
            const bool haplotype = !mode_excludes_genomic_filters &&
                best_haplotype != phased_observations.end() &&
                haplotype_filter.error_probability(
                    best_haplotype->key, static_cast<std::int64_t>(record->pos) + 1) > 0.0;

            // Reconstruct the observable subset of GATK ErrorProbabilities
            // for FilteringOutputStats.  Each vector is in concrete-ALT
            // order, and repeated sources for one filter are combined by max,
            // matching ErrorProbabilities' per-filter/per-allele reduction.
            std::vector<int> stats_alt_indices;
            for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                if (record->d.allele[alt_index] != nullptr && record->d.allele[alt_index][0] != '<')
                    stats_alt_indices.push_back(alt_index);
            }
            std::map<std::string, std::vector<double>> stats_probabilities;
            const auto stats_mmq_values = info_int_values(input_header, record, "MMQ");
            const int stats_ref_mmq = stats_mmq_values.empty() ? -1 : stats_mmq_values.front();
            auto add_stats_allele_probability = [&](const std::string& name, const int alt_index,
                                                    const double probability) {
                if (!std::isfinite(probability) || alt_index <= 0) return;
                const auto iterator = std::find(stats_alt_indices.begin(), stats_alt_indices.end(), alt_index);
                if (iterator == stats_alt_indices.end()) return;
                auto& values = stats_probabilities[name];
                if (values.empty()) values.assign(stats_alt_indices.size(), 0.0);
                const auto offset = static_cast<std::size_t>(iterator - stats_alt_indices.begin());
                values[offset] = std::max(values[offset], std::clamp(probability, 0.0, 1.0));
            };
            auto add_stats_site_probability = [&](const std::string& name, const double probability) {
                if (!std::isfinite(probability) || stats_alt_indices.empty()) return;
                auto& values = stats_probabilities[name];
                if (values.empty()) values.assign(stats_alt_indices.size(), 0.0);
                const auto bounded = std::clamp(probability, 0.0, 1.0);
                for (auto& value : values) value = std::max(value, bounded);
            };
            for (std::size_t offset = 0; offset < stats_alt_indices.size(); ++offset) {
                const auto alt_index = stats_alt_indices[offset];
                const auto alt_low_tlod = alt_index <= static_cast<int>(tlod_values.size()) &&
                    static_cast<double>(tlod_values[static_cast<std::size_t>(alt_index - 1)]) < options.min_tlod;
                add_stats_allele_probability("low_tlod", alt_index, alt_low_tlod ? 1.0 : 0.0);
                const auto alt_af = tumor_allele_fraction_for_alt(
                    input_header, record, tumor_sample_index, alt_index);
                add_stats_allele_probability("low_allele_frac", alt_index, alt_af >= 0.0 &&
                    alt_af < options.min_allele_fraction ? 1.0 : 0.0);
                const auto alt_mbq = info_allele_int(input_header, record, "MBQ", alt_index, true);
                add_stats_allele_probability("base_qual", alt_index,
                    options.min_median_base_quality >= 0 && alt_mbq &&
                    *alt_mbq < options.min_median_base_quality ? 1.0 : 0.0);
                auto alt_mmq = info_allele_int(input_header, record, "MMQ", alt_index, true);
                const auto alt_indel_length = indel_length_for_alt(record, alt_index);
                if (alt_indel_length >= options.long_indel_length && stats_ref_mmq >= 0) alt_mmq = stats_ref_mmq;
                add_stats_allele_probability("map_qual", alt_index,
                    options.min_median_mapping_quality >= 0 && alt_mmq &&
                    *alt_mmq < options.min_median_mapping_quality ? 1.0 : 0.0);
                const auto alt_mpos = info_allele_int(input_header, record, "MPOS", alt_index, false);
                add_stats_allele_probability("position", alt_index,
                    options.min_median_read_position >= 0 && alt_mpos &&
                    *alt_mpos > -1 && *alt_mpos < options.min_median_read_position ? 1.0 : 0.0);
                const auto strand_counts = info_alt_strand_counts(
                    input_header, record, "AS_SB_TABLE", alt_index - 1);
                const auto orientation_alt_counts = read_orientation_counts(
                    input_header, record, tumor_sample_index, alt_index);
                const auto forward = strand_counts ? strand_counts->first : orientation_alt_counts.f1r2;
                const auto reverse = strand_counts ? strand_counts->second : orientation_alt_counts.f2r1;
                const auto strand_total = forward >= 0 && reverse >= 0 ? forward + reverse : 0;
                const auto strand_error = strand_total > 0
                    ? 1.0 - 2.0 * static_cast<double>(std::min(forward, reverse)) /
                        static_cast<double>(strand_total) : 0.0;
                add_stats_allele_probability("strand_bias", alt_index, strand_error);
                add_stats_allele_probability("strict_strand", alt_index, strand_counts &&
                    (strand_counts->first == 0 || strand_counts->second == 0) ? 1.0 : 0.0);
                const auto unique_count = info_alt_int(
                    input_header, record, "AS_UNIQ_ALT_READ_COUNT", alt_index - 1);
                const auto alt_unique = unique_count.value_or(
                    format_alt_depth(input_header, record, tumor_sample_index, alt_index));
                add_stats_allele_probability("duplicate", alt_index,
                    options.unique_alt_read_count > 0 &&
                    (alt_unique < 0 || alt_unique <= options.unique_alt_read_count) ? 1.0 : 0.0);
                if (const auto iterator = weak_evidence_by_alt.find(alt_index);
                    iterator != weak_evidence_by_alt.end())
                    add_stats_allele_probability("weak_evidence", alt_index, iterator->second);
                if (const auto iterator = contamination_probability_by_alt.find(alt_index);
                    iterator != contamination_probability_by_alt.end())
                    add_stats_allele_probability("contamination", alt_index, iterator->second);
                if (const auto iterator = joint_error_by_alt.find(alt_index);
                    iterator != joint_error_by_alt.end())
                    add_stats_allele_probability("error_probability", alt_index, iterator->second);
            }
            add_stats_site_probability("germline", std::isfinite(germline_probability)
                ? germline_probability : (germline ? 1.0 : 0.0));
            if (!contamination_probability_evaluated)
                add_stats_site_probability("contamination", contamination ? 1.0 : 0.0);
            add_stats_site_probability("normal_artifact", normal_artifact_probability);
            add_stats_site_probability("orientation", has_orientation_prior
                ? orientation_artifact_probability_value : (orientation ? 1.0 : 0.0));
            add_stats_site_probability("multiallelic", multiallelic ? 1.0 : 0.0);
            add_stats_site_probability("fragment", fragment_length ? 1.0 : 0.0);
            add_stats_site_probability("n_ratio", n_ratio ? 1.0 : 0.0);
            add_stats_site_probability("panel_of_normals", panel_of_normals ? 1.0 : 0.0);
            add_stats_site_probability("slippage", slippage_probability);
            add_stats_site_probability("haplotype", haplotype ? 1.0 : 0.0);
            add_stats_site_probability("clustered_events", clustered_events ? 1.0 : 0.0);
            add_stats_site_probability("low_allele_frac", low_af ? 1.0 : 0.0);
            add_stats_site_probability("low_tlod", low_tlod ? 1.0 : 0.0);
            add_stats_site_probability("weak_evidence", weak_evidence_site ? 1.0 : 0.0);
            add_stats_site_probability("error_probability", joint_error_site ? 1.0 : 0.0);
            // FilteringOutputStats uses the Java ErrorProbabilities reduction:
            // max within each correlated ErrorType, then 1-prod(1-p) across
            // SEQUENCING/NON_SOMATIC/ARTIFACT.  Keep this reduction separate
            // from the VCF FILTER adapter above so the stats contract can use
            // every deterministic filter probability available in this pass.
            const auto is_sequencing_stats_filter = [](const std::string& name) {
                return name == "weak_evidence" || name == "low_tlod";
            };
            const auto is_non_somatic_stats_filter = [](const std::string& name) {
                return name == "germline" || name == "contamination";
            };
            std::vector<double> stats_combined_errors;
            stats_combined_errors.reserve(stats_alt_indices.size());
            std::vector<double> stats_sequencing_errors;
            std::vector<double> stats_non_somatic_errors;
            std::vector<double> stats_artifact_errors;
            std::vector<double> stats_orientation_errors;
            stats_sequencing_errors.reserve(stats_alt_indices.size());
            stats_non_somatic_errors.reserve(stats_alt_indices.size());
            stats_artifact_errors.reserve(stats_alt_indices.size());
            stats_orientation_errors.reserve(stats_alt_indices.size());
            for (const auto alt_index : stats_alt_indices) {
                double sequencing_error = std::numeric_limits<double>::quiet_NaN();
                double non_somatic_error = std::numeric_limits<double>::quiet_NaN();
                double artifact_error = std::numeric_limits<double>::quiet_NaN();
                const auto update_max = [](double& target, const double candidate) {
                    if (!std::isfinite(candidate)) return;
                    target = std::isfinite(target) ? std::max(target, candidate) : candidate;
                };
                const auto posterior_somatic = info_probability_at(
                    input_header, record, "PSOMATIC", alt_index);
                if (std::isfinite(posterior_somatic))
                    update_max(sequencing_error, 1.0 - posterior_somatic);
                const auto posterior_germline = info_probability_at(
                    input_header, record, "PGERMLINE", alt_index);
                if (std::isfinite(posterior_germline))
                    update_max(non_somatic_error, posterior_germline);
                const auto posterior_artifact = info_probability_at(
                    input_header, record, "PARTIFACT", alt_index);
                if (std::isfinite(posterior_artifact))
                    update_max(artifact_error, posterior_artifact);
                const auto posterior_orientation = info_probability_at(
                    input_header, record, "OBP", alt_index);
                if (std::isfinite(posterior_orientation))
                    update_max(artifact_error, posterior_orientation);
                for (const auto& entry : stats_probabilities) {
                    // error_probability is already the joint value and must
                    // not be fed back into its own ErrorType reduction.
                    if (entry.first == "error_probability") continue;
                    const auto iterator = std::find(
                        stats_alt_indices.begin(), stats_alt_indices.end(), alt_index);
                    if (iterator == stats_alt_indices.end()) continue;
                    const auto offset = static_cast<std::size_t>(
                        iterator - stats_alt_indices.begin());
                    if (offset >= entry.second.size()) continue;
                    const auto probability = entry.second[offset];
                    if (is_sequencing_stats_filter(entry.first))
                        update_max(sequencing_error, probability);
                    else if (is_non_somatic_stats_filter(entry.first))
                        update_max(non_somatic_error, probability);
                    else
                        update_max(artifact_error, probability);
                }
                stats_sequencing_errors.push_back(sequencing_error);
                stats_non_somatic_errors.push_back(non_somatic_error);
                stats_artifact_errors.push_back(artifact_error);
                stats_orientation_errors.push_back(-1.0);
            }
            if (!stats_sequencing_errors.empty()) {
                const auto probabilities = stats_sequencing_errors.size() == 1
                    ? std::vector<double>{combine_error_probabilities_kokkos(
                        stats_sequencing_errors[0], stats_non_somatic_errors[0],
                        stats_artifact_errors[0], stats_orientation_errors[0],
                        &kernel_telemetry)}
                    : combine_error_probabilities_batch_kokkos(
                        stats_sequencing_errors, stats_non_somatic_errors,
                        stats_artifact_errors, stats_orientation_errors,
                        &kernel_telemetry);
                for (auto probability : probabilities) {
                    if (!std::isfinite(probability) || probability < 0.0) probability = 0.0;
                    stats_combined_errors.push_back(std::clamp(probability, 0.0, 1.0));
                }
            }
            if (effective_joint_error_threshold >= 0.0 && !stats_combined_errors.empty()) {
                // The stats reduction is the authoritative joint posterior
                // once all filter-specific probabilities are available.  Use
                // the same per-ALT values for AS_FilterStatus and the site
                // FILTER so explicit threshold strategies do not silently
                // discard hard-filter artifact evidence.
                joint_error_probabilities = stats_combined_errors;
                joint_error_by_alt.clear();
                for (std::size_t offset = 0; offset < stats_alt_indices.size(); ++offset)
                    joint_error_by_alt[stats_alt_indices[offset]] = stats_combined_errors[offset];
                joint_error = std::any_of(
                    stats_combined_errors.begin(), stats_combined_errors.end(),
                    [&](const double probability) {
                        return std::isfinite(probability) &&
                               probability > effective_joint_error_threshold;
                    });
                joint_error_site = std::all_of(
                    stats_combined_errors.begin(), stats_combined_errors.end(),
                    [&](const double probability) {
                        return std::isfinite(probability) &&
                               probability > effective_joint_error_threshold;
                    });
                auto& error_probability_values = stats_probabilities["error_probability"];
                error_probability_values = stats_combined_errors;
            }
            std::vector<std::pair<std::string, std::vector<double>>> stats_probability_list;
            stats_probability_list.reserve(stats_probabilities.size());
            for (auto& entry : stats_probabilities)
                stats_probability_list.emplace_back(entry.first, std::move(entry.second));
            const auto stats_threshold = effective_joint_error_threshold >= 0.0
                ? effective_joint_error_threshold : options.filter_error_probability_threshold;
            filtering_output_stats.record(stats_combined_errors, stats_probability_list,
                                          std::max(0.0, stats_threshold - 1.0e-10));
            // Build the Number=A AS_FilterStatus contract from the same
            // annotations consumed by the Rust/GATK hard filters.  Site-level
            // filters are copied to every concrete ALT; allele-specific
            // annotations remain isolated to their own ALT.
            const int as_alt_count = std::max(0, record->n_allele - 1);
            std::vector<std::vector<std::string>> allele_filters(
                static_cast<std::size_t>(as_alt_count));
            const auto mmq_values = info_int_values(input_header, record, "MMQ");
            const int ref_mmq = mmq_values.empty() ? -1 : mmq_values.front();
            bool strict_strand = false;
            bool duplicate = false;
            auto add_allele_filter = [&](int alt_index, const char* name) {
                if (alt_index <= 0 || alt_index > as_alt_count) return;
                auto& values = allele_filters[static_cast<std::size_t>(alt_index - 1)];
                if (std::find(values.begin(), values.end(), name) == values.end()) values.emplace_back(name);
            };
            for (int alt_index = 1; alt_index <= as_alt_count; ++alt_index) {
                if (record->d.allele[alt_index] != nullptr && record->d.allele[alt_index][0] != '<') {
                    const auto iterator = contamination_probability_by_alt.find(alt_index);
                    if (iterator != contamination_probability_by_alt.end() &&
                        std::isfinite(iterator->second) &&
                        iterator->second > options.max_contamination_probability)
                        add_allele_filter(alt_index, "contamination");
                }
                if (alt_index <= static_cast<int>(weak_evidence_probabilities.size()) &&
                    std::isfinite(weak_evidence_probabilities[static_cast<std::size_t>(alt_index - 1)]) &&
                    weak_evidence_probabilities[static_cast<std::size_t>(alt_index - 1)] >
                        gatk_error_threshold)
                    add_allele_filter(alt_index, "weak_evidence");
                const auto alt_mbq = info_allele_int(input_header, record, "MBQ", alt_index, true);
                if (options.min_median_base_quality >= 0 && alt_mbq &&
                    *alt_mbq < options.min_median_base_quality) {
                    add_allele_filter(alt_index, "base_qual");
                }
                auto alt_mmq = info_allele_int(input_header, record, "MMQ", alt_index, true);
                const auto indel_length = indel_length_for_alt(record, alt_index);
                if (indel_length >= options.long_indel_length && ref_mmq >= 0) alt_mmq = ref_mmq;
                if (options.min_median_mapping_quality >= 0 && alt_mmq &&
                    *alt_mmq < options.min_median_mapping_quality) {
                    add_allele_filter(alt_index, "map_qual");
                }
                const auto alt_mpos = info_allele_int(input_header, record, "MPOS", alt_index, false);
                if (options.min_median_read_position >= 0 && alt_mpos && *alt_mpos > -1 &&
                    *alt_mpos < options.min_median_read_position) {
                    add_allele_filter(alt_index, "position");
                }
                const auto strand_counts = info_alt_strand_counts(
                    input_header, record, "AS_SB_TABLE", alt_index - 1);
                if (options.min_reads_per_strand > 0 && strand_counts &&
                    (strand_counts->first == 0 || strand_counts->second == 0)) {
                    add_allele_filter(alt_index, "strict_strand");
                    strict_strand = true;
                }
                const auto unique_count = info_alt_int(
                    input_header, record, "AS_UNIQ_ALT_READ_COUNT", alt_index - 1);
                const int alt_unique = unique_count.value_or(
                    format_alt_depth(input_header, record, tumor_sample_index, alt_index));
                if (options.unique_alt_read_count > 0 &&
                    (alt_unique < 0 || alt_unique <= options.unique_alt_read_count)) {
                    add_allele_filter(alt_index, "duplicate");
                    duplicate = true;
                }
                if (strand_artifact_applies(alt_index))
                    add_allele_filter(alt_index, "strand_bias");
                const auto alt_af = tumor_allele_fraction_for_alt(
                    input_header, record, tumor_sample_index, alt_index);
                if (options.min_allele_fraction > 0.0F && alt_af < options.min_allele_fraction) {
                    add_allele_filter(alt_index, "low_allele_frac");
                }
                // Keep low-TLOD status aligned with the same per-ALT vector
                // used for the site-level reduction.  Scalar legacy TLOD is
                // broadcast, while missing Number=A elements remain
                // unfiltered here and are handled by the site-level
                // fail-closed boundary below.
                const auto alt_tlod = info_probability_at(
                    input_header, record, "TLOD", alt_index);
                if (std::isfinite(alt_tlod) && alt_tlod < options.min_tlod)
                    add_allele_filter(alt_index, "low_tlod");
            }
            bool has_allele_status = false;
            std::vector<std::string> allele_status;
            allele_status.reserve(static_cast<std::size_t>(as_alt_count));
            for (const auto& values : allele_filters) {
                if (values.empty()) allele_status.emplace_back("SITE");
                else {
                    has_allele_status = true;
                    allele_status.push_back(join_strings(values, ","));
                }
            }
            if (as_alt_count > 0) {
                const auto encoded_status = join_strings(allele_status, "|");
                if (bcf_update_info_string(output_header, record, "AS_FilterStatus",
                                           encoded_status.c_str()) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AS_FilterStatus");
                if (has_allele_status) ++as_filter_records;
            }
            if (germline || germline_filter) ++germline_records;
            if (low_tlod) ++low_tlod_records;
            if (contamination) ++contamination_records;
            if (orientation) ++orientation_records;
            if (low_af) ++low_af_records;
            if (strand || strand_artifact_site) ++strand_records;
            if (low_alt_reads) ++low_alt_read_records;
            if (multiallelic) ++multiallelic_records;
            if (low_median_base_quality) ++low_median_base_quality_records;
            if (low_median_mapping_quality) ++low_median_mapping_quality_records;
            if (low_median_read_position) ++low_median_read_position_records;
            if (fragment_length) ++fragment_length_records;
            if (n_ratio) ++n_ratio_records;
            if (strict_strand) ++strict_strand_records;
            if (duplicate) ++duplicate_records;
            if (panel_of_normals) ++panel_of_normals_records;
            if (normal_artifact) ++normal_artifact_records;
            if (weak_evidence) ++weak_evidence_records;
            if (slippage) ++slippage_records;
            if (joint_error) ++joint_error_records;
            if (low_somatic_probability) ++low_somatic_probability_records;
            if (high_germline_probability) ++high_germline_probability_records;
            if (high_artifact_probability) ++high_artifact_probability_records;
            if (orientation_prior_filter) ++orientation_records;
            if (clustered_events) ++clustered_events_records;
            if (haplotype) ++haplotype_records;
            // Mutect2FilteringEngine first turns an allele-specific filter
            // into a site FILTER only when *every* concrete ALT failed that
            // filter.  Variant-level filters retain their posterior.  It then
            // reports only filters within the largest posterior (with the
            // GATK 0.1 reporting floor), rather than dumping every positive
            // signal into FILTER.  This is why a weak-evidence call must not
            // acquire an unrelated low-probability germline label.
            std::map<std::string, double> site_filter_probabilities;
            const auto add_site_probability = [&](const char* name, const double probability) {
                if (!std::isfinite(probability) || probability <= gatk_error_threshold) return;
                auto iterator = site_filter_probabilities.find(name);
                if (iterator == site_filter_probabilities.end())
                    site_filter_probabilities.emplace(name, probability);
                else
                    iterator->second = std::max(iterator->second, probability);
            };
            const auto every_concrete_alt_has = [&](const char* name) {
                bool observed = false;
                for (int alt_index = 1; alt_index <= as_alt_count; ++alt_index) {
                    if (record->d.allele[alt_index] == nullptr ||
                        record->d.allele[alt_index][0] == '<') continue;
                    observed = true;
                    const auto& values = allele_filters[static_cast<std::size_t>(alt_index - 1)];
                    if (std::find(values.begin(), values.end(), name) == values.end()) return false;
                }
                return observed;
            };
            static constexpr std::array<const char*, 9> kAlleleSpecificFilters{{
                "weak_evidence", "base_qual", "map_qual", "duplicate", "strand_bias",
                "strict_strand", "position", "low_allele_frac", "contamination"}};
            for (const auto* name : kAlleleSpecificFilters)
                if (every_concrete_alt_has(name)) add_site_probability(name, 1.0);
            if (germline_filter)
                add_site_probability("germline", germline_filter_probability);
            if (normal_artifact)
                add_site_probability("normal_artifact", normal_artifact_probability);
            if (slippage)
                add_site_probability("slippage", slippage_probability);
            if (orientation_prior_filter)
                add_site_probability("orientation", orientation_artifact_probability_value);
            if (multiallelic) add_site_probability("multiallelic", 1.0);
            if (fragment_length) add_site_probability("fragment", 1.0);
            if (n_ratio) add_site_probability("n_ratio", 1.0);
            if (panel_of_normals) add_site_probability("panel_of_normals", 1.0);
            if (clustered_events) add_site_probability("clustered_events", 1.0);
            if (haplotype) add_site_probability("haplotype", 1.0);
            // GATK emits FAIL only when every concrete ALT has an
            // allele-specific failure and no variant/site filter was already
            // selected.  A record such as `low_tlod|SITE` is therefore PASS:
            // the second ALT remains callable even though the first does not.
            bool any_concrete_alt_passes = false;
            for (int alt_index = 1; alt_index <= as_alt_count; ++alt_index) {
                if (record->d.allele[alt_index] == nullptr ||
                    record->d.allele[alt_index][0] == '<') continue;
                if (allele_filters[static_cast<std::size_t>(alt_index - 1)].empty()) {
                    any_concrete_alt_passes = true;
                    break;
                }
            }
            if (site_filter_probabilities.empty() && has_allele_status &&
                !any_concrete_alt_passes)
                add_site_probability("FAIL", 1.0);

            double maximum_site_error = 0.0;
            for (const auto& entry : site_filter_probabilities)
                maximum_site_error = std::max(maximum_site_error, entry.second);
            const double report_cutoff = std::min(maximum_site_error, 0.1);
            std::vector<int> filters;
            filters.reserve(site_filter_probabilities.size());
            for (const auto& entry : site_filter_probabilities) {
                if (entry.second + kFilterEpsilon < report_cutoff) continue;
                const auto id = bcf_hdr_id2int(output_header, BCF_DT_ID, entry.first.c_str());
                if (id >= 0) filters.push_back(id);
            }
            if (has_orientation_prior) {
                // FilterMutectCalls serializes the orientation posterior with
                // QualityUtils.errorProbToQual(), the same bounded phred
                // conversion used by GERMQ/CONTQ.  The previous native
                // adapter returned 1000 for a zero artifact posterior and 0
                // for a certain artifact; both values are outside GATK's
                // [1,93] SAMUtils quality range and made otherwise compatible
                // records differ at the writer boundary.
                const auto quality = static_cast<float>(gatk_error_probability_to_quality(
                    std::clamp(static_cast<double>(orientation_artifact_probability_value), 0.0, 1.0)));
                if (bcf_update_info_float(output_header, record, "ROQ", &quality, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write orientation ROQ annotation");
            }
            if (std::isfinite(germline_filter_probability)) {
                const int quality = gatk_error_probability_to_quality(
                    germline_filter_probability);
                if (bcf_update_info_int32(output_header, record, "GERMQ", &quality, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write germline GERMQ annotation");
            }
            if (contamination_probability_evaluated && !contamination_probability_by_alt.empty()) {
                double maximum_contamination_probability = 0.0;
                for (const auto& entry : contamination_probability_by_alt)
                    if (std::isfinite(entry.second))
                        maximum_contamination_probability = std::max(
                            maximum_contamination_probability, std::clamp(entry.second, 0.0, 1.0));
                const auto quality = maximum_contamination_probability <= 0.0
                    ? 93.0F
                    : static_cast<float>(std::clamp(static_cast<int>(std::lround(
                        -10.0 * std::log10(maximum_contamination_probability))), 1, 93));
                if (bcf_update_info_float(output_header, record, "CONTQ", &quality, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write contamination CONTQ annotation");
            }
            if (!options.mitochondria_mode && repeat_unit && repeats_per_allele.size() >= 2 &&
                std::isfinite(slippage_probability)) {
                int quality = 93;
                if (slippage_probability >= 1.0) quality = 1;
                else if (slippage_probability > 0.0)
                    quality = std::clamp(static_cast<int>(std::lround(
                        -10.0 * std::log10(slippage_probability))), 1, 93);
                if (bcf_update_info_int32(output_header, record, "STRQ", &quality, 1) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write slippage STRQ annotation");
            }
            const auto filter_count = static_cast<int>(filters.size());
            if (filter_count == 0) {
                ++pass_records;
                bcf_update_filter(output_header, record, nullptr, 0);
                bcf_add_filter(output_header, record, bcf_hdr_id2int(output_header, BCF_DT_ID, "PASS"));
            } else if (bcf_update_filter(output_header, record, filters.data(), filter_count) != 0) {
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update Mutect2 FILTER");
            }
            if (options.sites_only_vcf_output && bcf_subset_format(output_header, record) != 0)
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: cannot write sites-only VCF record");
            if (bcf_write(output, output_header, record) != 0)
                throw std::runtime_error("cannot write filtered Mutect2 record");
            bcf_clear(record);
        }
        bcf_close(output); output = nullptr;
        bcf_close(input); input = nullptr;
        std::string index_path;
        if (options.create_index) {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build filtered VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const auto& stats_path = stats_output_path;
        const bool filtering_stats_table = suffix(stats_path, ".table") || suffix(stats_path, ".tsv");
        const std::string json_stats_path = filtering_stats_table ? stats_path + ".json" : stats_path;
        std::ofstream stats(json_stats_path);
        if (!stats) throw std::runtime_error("cannot write FilterMutectCalls stats: " + json_stats_path);
        stats << std::setprecision(17)
              << "{\"schema_version\":1,\"tool\":\"FilterMutectCalls\",\"status\":\"prototype\"," 
              << "\"input_records\":" << input_records << ",\"pass_records\":" << pass_records
              << ",\"interval_set_rule\":\"" << filter_interval_set_rule_name(options.interval_set_rule) << "\""
              << ",\"interval_skipped_records\":" << interval_skipped_records
              << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
              << ",\"microbial_mode\":" << (options.microbial_mode ? "true" : "false")
              << ",\"interval_list_inputs\":" << interval_stats.files
              << ",\"interval_list_records\":" << interval_stats.records
              << ",\"low_tlod_records\":" << low_tlod_records << ",\"germline_records\":"
              << germline_records << ",\"contamination_records\":" << contamination_records
              << ",\"contamination_posterior_records\":" << contamination_posterior_records
              << ",\"contamination_posterior_alleles\":" << contamination_posterior_alleles
              << ",\"orientation_records\":" << orientation_records
              << ",\"low_af_records\":" << low_af_records
              << ",\"strand_records\":" << strand_records
              << ",\"low_alt_read_records\":" << low_alt_read_records
              << ",\"multiallelic_records\":" << multiallelic_records
              << ",\"low_median_base_quality_records\":" << low_median_base_quality_records
              << ",\"low_median_mapping_quality_records\":" << low_median_mapping_quality_records
              << ",\"low_median_read_position_records\":" << low_median_read_position_records
              << ",\"fragment_length_records\":" << fragment_length_records
              << ",\"n_ratio_records\":" << n_ratio_records
              << ",\"allele_specific_strand_records\":" << allele_specific_strand_records
              << ",\"allele_specific_unique_records\":" << allele_specific_unique_records
              << ",\"strict_strand_records\":" << strict_strand_records
              << ",\"duplicate_records\":" << duplicate_records
              << ",\"panel_of_normals_records\":" << panel_of_normals_records
              << ",\"as_filter_records\":" << as_filter_records
              << ",\"normal_artifact_records\":" << normal_artifact_records
              << ",\"normal_pileup_p_value_threshold\":" << options.normal_pileup_p_value_threshold
              << ",\"weak_evidence_records\":" << weak_evidence_records
              << ",\"slippage_records\":" << slippage_records
              << ",\"joint_error_records\":" << joint_error_records
              << ",\"joint_error_threshold\":" << effective_joint_error_threshold
              << ",\"joint_error_threshold_observations\":" << joint_error_threshold_observations
              << ",\"joint_error_threshold_learned\":"
              << (effective_joint_error_threshold >= 0.0 ? "true" : "false")
              << ",\"joint_threshold_strategy\":\"" << options.joint_threshold_strategy << "\""
              << ",\"empirical_somatic_model_learned\":"
              << (empirical_somatic_model.learned ? "true" : "false")
              << ",\"empirical_somatic_model_records\":"
              << empirical_somatic_model_records
              << ",\"empirical_dynamic_peak_clustering\":"
              << (empirical_somatic_model.learned ? "true" : "false")
              << ",\"empirical_background_alpha\":" << empirical_somatic_model.background_alpha
              << ",\"empirical_background_beta\":" << empirical_somatic_model.background_beta
              << ",\"empirical_high_af_alpha\":" << empirical_somatic_model.high_af_alpha
              << ",\"empirical_high_af_beta\":" << empirical_somatic_model.high_af_beta
              << ",\"callable_sites\":";
        if (std::isfinite(empirical_somatic_model.callable_sites))
            stats << empirical_somatic_model.callable_sites;
        else
            stats << "null";
        stats << ",\"empirical_variant_priors\":[";
        for (std::size_t index = 0; index < empirical_somatic_model.log_variant_priors.size(); ++index) {
            if (index != 0) stats << ',';
            stats << empirical_somatic_model.log_variant_priors[index];
        }
        stats << ']'
              << ",\"filtering_output_stats\":" << filtering_output_stats.to_json()
              << ",\"empirical_cluster_count\":" << empirical_somatic_model.cluster_count
              << ",\"empirical_cluster_weights\":[";
        for (int cluster = 0; cluster < empirical_somatic_model.cluster_count; ++cluster) {
            if (cluster != 0) stats << ',';
            stats << std::exp(empirical_somatic_model.cluster_log_weights[static_cast<std::size_t>(cluster)]);
        }
        stats << "]"
              << ",\"empirical_cluster_means\":[";
        for (int cluster = 0; cluster < empirical_somatic_model.cluster_count; ++cluster) {
            if (cluster != 0) stats << ',';
            stats << empirical_somatic_model.cluster_mean[static_cast<std::size_t>(cluster)];
        }
        stats << ']'
              << ",\"filter_error_probability_threshold\":" << options.filter_error_probability_threshold
              << ",\"log_snv_prior\":" << options.log_snv_prior
              << ",\"log_indel_prior\":" << options.log_indel_prior
              << ",\"log_artifact_prior\":" << options.log_artifact_prior
              << ",\"min_slippage_length\":" << options.min_slippage_length
              << ",\"slippage_rate\":" << options.slippage_rate
              << ",\"input_stats\":\"" << json_escape(stats_input_path) << "\""
              << ",\"clustered_events_records\":" << clustered_events_records
              << ",\"min_tlod\":" << options.min_tlod
              << ",\"max_contamination\":" << effective_contamination
              << ",\"max_contamination_probability\":" << options.max_contamination_probability
              << ",\"contamination_source\":\"" << contamination_source << "\""
              << ",\"contamination_sample_estimates\":" << contamination_by_sample.size()
              << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
              << ",\"normal_artifact_lod\":" << options.normal_artifact_lod
              << ",\"min_orientation_balance\":" << options.min_orientation_balance
              << ",\"min_allele_fraction\":" << options.min_allele_fraction
              << ",\"min_reads_per_strand\":" << options.min_reads_per_strand
              << ",\"unique_alt_read_count\":" << options.unique_alt_read_count
              << ",\"max_alt_allele_count\":" << options.max_alt_allele_count
              << ",\"min_median_base_quality\":" << options.min_median_base_quality
              << ",\"min_median_mapping_quality\":" << options.min_median_mapping_quality
              << ",\"long_indel_length\":" << options.long_indel_length
              << ",\"min_median_read_position\":" << options.min_median_read_position
              << ",\"max_median_fragment_length_difference\":"
              << options.max_median_fragment_length_difference
              << ",\"max_n_ratio\":" << json_double(options.max_n_ratio)
              << ",\"genomic_hard_filters_disabled\":"
              << ((options.mitochondria_mode || options.microbial_mode) ? "true" : "false")
              << ",\"low_somatic_probability_records\":" << low_somatic_probability_records
              << ",\"high_germline_probability_records\":" << high_germline_probability_records
              << ",\"high_artifact_probability_records\":" << high_artifact_probability_records
              << ",\"min_somatic_probability\":" << options.min_somatic_probability
              << ",\"max_germline_probability\":" << options.max_germline_probability
              << ",\"max_artifact_probability\":" << options.max_artifact_probability
              << ",\"orientation_prior_records\":" << orientation_prior_records
              << ",\"orientation_prior_sample_evaluations\":" << orientation_prior_sample_evaluations
              << ",\"orientation_depth_disagreements\":" << orientation_depth_disagreements
              << ",\"orientation_threshold_strategy\":\"" << options.orientation_threshold_strategy << "\""
              << ",\"orientation_threshold_observations\":" << orientation_threshold_observations
              << ",\"orientation_threshold_learned\":"
              << (options.threshold_strategy_explicit &&
                  options.orientation_threshold_strategy != "CONSTANT" &&
                  options.max_orientation_artifact_probability < 0.0F ? "true" : "false")
              << ",\"max_orientation_artifact_probability\":" << effective_orientation_threshold
              << ",\"tumor_sample\":\"" << json_escape(selected_tumor_name) << "\""
              << ",\"max_events_in_region\":" << options.max_events_in_region
              << ",\"max_events_in_haplotype\":" << options.max_events_in_haplotype
              << ",\"haplotype_records\":" << haplotype_records
              << ",\"haplotype_learning_records\":" << haplotype_filter.accumulated_records()
              << ",\"learned_haplotypes\":" << haplotype_filter.learned_haplotypes()
              << ",\"max_intra_haplotype_distance\":" << options.max_intra_haplotype_distance << "}\n";
        stats.flush();
        if (filtering_stats_table) {
            std::vector<std::pair<std::string, std::string>> clustering_metadata{
                {"Ln prior of SNV", std::to_string(empirical_somatic_model.log_snv_prior)},
                {"Ln prior of indel", std::to_string(empirical_somatic_model.log_indel_prior)},
                {"Background beta-binomial cluster", "weight = " +
                    std::to_string(empirical_somatic_model.background_weight)},
                {"High-AF beta-binomial cluster", "weight = " +
                    std::to_string(empirical_somatic_model.high_af_weight)}};
            const auto table_threshold = effective_joint_error_threshold >= 0.0
                ? effective_joint_error_threshold : options.filter_error_probability_threshold;
            filtering_output_stats.write_gatk_table(stats_path, table_threshold, clustering_metadata);
        }
        if (!file_complete(options.output) || !file_complete(json_stats_path) ||
            (filtering_stats_table && !file_complete(stats_path)) ||
            (!index_path.empty() && !file_complete(index_path)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: FilterMutectCalls output or sidecar is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write FilterMutectCalls manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"FilterMutectCalls\",\"implementation\":\"fastgatk-filter-mutect-calls\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"tlod_filter\":true,\"germline_filter\":true,\"contamination_filter\":"
                 << (effective_contamination >= 0.0F ? "true" : "false")
                 << ",\"contamination_posterior_filter\":"
                 << (effective_contamination >= 0.0F &&
                     bcf_hdr_id2int(input_header, BCF_DT_ID, "POPAF") >= 0 ? "true" : "false")
                 << ",\"contamination_quality_annotation\":true"
                 << ",\"contamination_resource\":"
                 << ((!options.contamination_table.empty() || !options.tumor_segmentation.empty()) ? "true" : "false")
                 << ",\"interval_subset\":" << (!intervals.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << filter_interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
                 << ",\"microbial_mode\":" << (options.microbial_mode ? "true" : "false")
                 << ",\"genomic_hard_filters_disabled\":"
                 << ((options.mitochondria_mode || options.microbial_mode) ? "true" : "false")
                 << ",\"interval_list_inputs\":" << interval_stats.files
                 << ",\"interval_list_records\":" << interval_stats.records
                 << ",\"orientation_filter\":" << (options.min_orientation_balance >= 0.0F ? "true" : "false")
                 << ",\"min_allele_fraction_filter\":" << (options.min_allele_fraction > 0.0F ? "true" : "false")
                 << ",\"strand_support_filter\":" << (options.min_reads_per_strand > 0 ? "true" : "false")
                 << ",\"unique_alt_read_filter\":" << (options.unique_alt_read_count > 0 ? "true" : "false")
                 << ",\"multiallelic_filter\":"
                 << ((options.mitochondria_mode || options.microbial_mode) ? "false" : "true")
                 << ",\"median_base_quality_filter\":"
                 << (options.min_median_base_quality >= 0 ? "true" : "false")
                 << ",\"median_mapping_quality_filter\":"
                 << (options.min_median_mapping_quality >= 0 ? "true" : "false")
                 << ",\"median_read_position_filter\":"
                 << (options.min_median_read_position >= 0 ? "true" : "false")
                 << ",\"fragment_length_filter\":"
                 << ((options.max_median_fragment_length_difference >= 0 &&
                      !options.mitochondria_mode && !options.microbial_mode) ? "true" : "false")
                 << ",\"n_ratio_filter\":"
                 << (options.max_n_ratio >= 0.0 ? "true" : "false")
                 << ",\"strict_strand_filter\":"
                 << (options.min_reads_per_strand > 0 ? "true" : "false")
                 << ",\"panel_of_normals_filter\":true"
                 << ",\"normal_artifact_filter\":true"
                 << ",\"weak_evidence_filter\":true"
                 << ",\"slippage_filter\":true"
                 << ",\"slippage_initial_model\":true"
                 << ",\"joint_error_probability_filter\":"
                 << (effective_joint_error_threshold >= 0.0 ? "true" : "false")
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"joint_error_threshold_learning\":"
                 << (effective_joint_error_threshold >= 0.0 ? "true" : "false")
                 << ",\"empirical_somatic_model\":"
                 << (empirical_somatic_model.learned ? "true" : "false")
                 << ",\"empirical_dynamic_peak_clustering\":"
                 << (empirical_somatic_model.learned ? "true" : "false")
                 << ",\"empirical_cluster_count\":" << empirical_somatic_model.cluster_count
                 << ",\"callable_site_prior_learning\":"
                 << (empirical_somatic_model.learned && std::isfinite(empirical_somatic_model.callable_sites) ? "true" : "false")
                 << ",\"gatk_prior_flags\":true"
                 << ",\"long_indel_length\":" << options.long_indel_length
                 << ",\"full_empirical_threshold_learning\":"
                 << (empirical_learning_requested ? "true" : "false")
                 << ",\"as_filter_status\":true"
                 << ",\"allele_specific_strand_filter\":true"
                 << ",\"allele_specific_unique_alt_filter\":true"
                 << ",\"somatic_probability_filter\":"
                 << (options.min_somatic_probability > 0.0F ? "true" : "false")
                 << ",\"germline_probability_filter\":"
                 << (options.max_germline_probability >= 0.0F ? "true" : "false")
                 << ",\"artifact_probability_filter\":"
                 << (options.max_artifact_probability >= 0.0F ? "true" : "false")
                 << ",\"orientation_prior_filter\":"
                 << (!options.orientation_priors.empty() && effective_orientation_threshold >= 0.0F ? "true" : "false")
                 << ",\"orientation_prior_table\":"
                 << (!options.orientation_priors.empty() ? "true" : "false")
                 << ",\"orientation_weighted_median\":"
                 << (!options.orientation_priors.empty() ? "true" : "false")
                 << ",\"orientation_threshold_learning\":"
                 << (options.threshold_strategy_explicit && !options.orientation_priors.empty() ? "true" : "false")
                 << ",\"kokkos_filter_kernel\":true"
                 << ",\"input_stats_preserved\":" << (input_stats_is_gatk_table ? "true" : "false")
                 << ",\"tumor_sample_selection\":true"
                 << ",\"clustered_events_filter\":"
                 << ((options.mitochondria_mode || options.microbial_mode) ? "false" : "true")
                 << ",\"haplotype_filter\":"
                 << ((options.mitochondria_mode || options.microbial_mode) ? "false" : "true")
                 << ",\"haplotype_learning\":" << (has_phasing_fields ? "true" : "false")
                 << ",\"filtering_stats_format\":\""
                 << (filtering_stats_table ? "gatk-filter-stats-table" : "native-json") << "\""
                 << ",\"contamination_source\":\"" << contamination_source << "\""
                 << ",\"contamination_sample_estimates\":" << contamination_by_sample.size()
                 << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"vcf_index\":"
                 << (index_path.empty() ? "false" : "true") << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":"
                 << (file_complete(options.output) ? "true" : "false") << "},{\"path\":\""
                 << json_escape(stats_path) << "\",\"kind\":\"filter-stats\",\"complete\":"
                 << (file_complete(stats_path) ? "true" : "false") << "}"
                 << (filtering_stats_table ? ",{\"path\":\"" + json_escape(json_stats_path) +
                     "\",\"kind\":\"filter-stats-json\",\"complete\":" +
                     (file_complete(json_stats_path) ? "true}" : "false}") : "")
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":" +
                     (file_complete(index_path) ? "true}" : "false}"))
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"kernel_execution_space\":\"" << json_escape(kernel_telemetry.execution_space) << "\""
                 << ",\"kernel_execution_policy\":\"RangePolicy\""
                 << ",\"kernel_batches\":" << kernel_telemetry.batches
                 << ",\"kernel_observations\":" << kernel_telemetry.observations
                 << ",\"kernel_batch_calls\":" << kernel_telemetry.batch_calls
                 << ",\"kernel_batch_observations\":" << kernel_telemetry.batch_observations
                 << ",\"kernel_prepare_seconds\":" << kernel_telemetry.prepare_seconds
                 << ",\"kernel_execute_seconds\":" << kernel_telemetry.execute_seconds
                 << ",\"input_records\":" << input_records << ",\"pass_records\":" << pass_records
                 << ",\"interval_skipped_records\":" << interval_skipped_records
                 << ",\"interval_set_rule\":\"" << filter_interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
                 << ",\"microbial_mode\":" << (options.microbial_mode ? "true" : "false")
                 << ",\"orientation_depth_disagreements\":" << orientation_depth_disagreements
                 << ",\"allele_specific_strand_records\":" << allele_specific_strand_records
                 << ",\"allele_specific_unique_records\":" << allele_specific_unique_records
                 << ",\"strict_strand_records\":" << strict_strand_records
                 << ",\"duplicate_records\":" << duplicate_records
                 << ",\"panel_of_normals_records\":" << panel_of_normals_records
                 << ",\"as_filter_records\":" << as_filter_records
                 << ",\"normal_artifact_records\":" << normal_artifact_records
                 << ",\"contamination_posterior_records\":" << contamination_posterior_records
                 << ",\"contamination_posterior_alleles\":" << contamination_posterior_alleles
                 << ",\"weak_evidence_records\":" << weak_evidence_records
                 << ",\"slippage_records\":" << slippage_records
                 << ",\"haplotype_records\":" << haplotype_records
                 << ",\"learned_haplotypes\":" << haplotype_filter.learned_haplotypes()
                 << ",\"max_n_ratio\":" << json_double(options.max_n_ratio)
                 << ",\"filtering_output_pass\":" << filtering_output_stats.pass_count()
                 << ",\"filtering_output_TP\":" << filtering_output_stats.true_positive_count()
                 << ",\"filtering_output_FP\":" << filtering_output_stats.false_positive_count()
                 << ",\"filtering_output_FN\":" << filtering_output_stats.false_negative_count()
                 << ",\"filtering_output_filter_count\":" << filtering_output_stats.filter_count()
                 << "}}\n";
        std::cout << "{\"tool\":\"FilterMutectCalls\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"pass_records\":" << pass_records
                  << ",\"interval_skipped_records\":" << interval_skipped_records
                  << ",\"interval_set_rule\":\"" << filter_interval_set_rule_name(options.interval_set_rule)
                  << "\"}\n";
        bcf_destroy(record);
        bcf_hdr_destroy(output_header);
        bcf_hdr_destroy(input_header);
        if (reference_index != nullptr) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        bcf_destroy(record);
        if (output) bcf_close(output);
        if (input) bcf_close(input);
        if (output_header) bcf_hdr_destroy(output_header);
        if (input_header) bcf_hdr_destroy(input_header);
        if (reference_index != nullptr) fai_destroy(reference_index);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for FilterMutectCalls");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto options = parse(argc, argv);
        Kokkos::initialize();
        initialized = true;
        const auto status = run_tool(options, fastgatk::runtime::ResourceSnapshot::probe());
        Kokkos::finalize();
        initialized = false;
        return status;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << "fastgatk-filter-mutect-calls: " << error.what() << '\n';
        return 2;
    }
}
