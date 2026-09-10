// Kokkos kernels for the FilterMutectCalls Java SomaticClusteringModel path.
// Host code (filter_mutect_tool.cpp) owns VCF passes, record() exclusion and
// the serial EM control flow; the per-datum numeric work (beta-binomial /
// binomial masses, corrected log likelihoods, responsibilities, cluster
// totals) runs here on Views so CPU and device callers share one numeric
// implementation built on gamma_math.hpp.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

constexpr std::size_t kFmcMaxClusters = 7;

struct FmcClusterParams {
    std::vector<double> log_cluster_weights;
    std::vector<double> cluster_alpha;
    std::vector<double> cluster_beta;
    // kind 0 = beta-binomial (GATK BetaBinomialCluster), 1 = fuzzy binomial.
    std::vector<int> cluster_kind;
};

struct FmcEvidence {
    std::vector<int> total_count;
    std::vector<int> alt_count;
    std::vector<double> tumor_log_odds;
    std::vector<double> artifact_prob;
    std::vector<double> non_somatic_prob;
    // Indel length per datum (SNV = 0).  Empty means all SNVs.
    std::vector<int> indel_length;
    // Optional per-datum population AF and contamination estimate used by the
    // germline/contamination posteriors.  Empty entries default to 0.01/0.1.
    std::vector<double> population_af;
    std::vector<double> contamination;
    // Mutect2 NLOD is the log10 odds that an ALT is absent from the matched
    // normal.  GermlineFilter negates it before combining the normal evidence
    // with the tumour germline-vs-somatic likelihood.  Empty means no normal
    // evidence, matching tumour-only mode.
    std::vector<double> normal_log10_odds;
    // Optional per-datum gate: germline posterior is evaluated only when the
    // upstream VCF carried a POPAF annotation (GATK GermlineFilter skips the
    // site otherwise and the NON_SOMATIC group then only sees contamination).
    // Empty = evaluate for every datum (backwards compatible).
    std::vector<std::uint8_t> germline_evaluated;
};

struct FmcErrorProbabilitiesResult {
    std::vector<double> sequencing_error;   // per datum
    std::vector<double> non_somatic_error;  // max(germline, contamination)
    std::vector<double> germline_error;     // per datum (diagnostics/tests)
    std::vector<double> contamination_error;  // per datum (diagnostics/tests)
    bool valid = false;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Batch sequencing/germline/contamination posteriors for a pass of evidence
// under the current model.  Host accumulates records and calls this once per
// pass so every per-datum numeric in the learning pipeline is kernel-owned.
FmcErrorProbabilitiesResult fmc_error_probabilities_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in);

// Per-split-iteration quantities the Java SomaticClusteringModel
// initialization needs, computed under the current (iteration-start) model:
// the sequencing posterior and Java probabilityOfSomaticVariant per datum,
// the count-based somatic mixture log likelihood and background-cluster
// posterior (softmax[0]) per datum, and, when a quantile grid is supplied,
// the Java quantile background responsibilities
// sum_i exp(log_binom(total_i, alt_i, q)) * bg_prob_i * (total_i + 1) with
// bg_prob_i = fixed_somatic_i * background_given_i.  The Host still owns the
// quantile grid bookkeeping, peak/trapezoid scan and BIC control flow, but
// every lgamma/binomial/logsumexp value is kernel-owned.
struct FmcSplitQuantitiesResult {
    std::vector<double> sequencing_error;           // per datum
    std::vector<double> somatic_probability;        // per datum (1-art)(1-ns)(1-seq)
    std::vector<double> log_somatic_mixture;        // per datum
    std::vector<double> background_given_somatic;   // per datum softmax[0]
    std::vector<double> background_probability;     // per datum
    std::vector<double> quantile_responsibility;    // per quantile
    bool valid = false;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

FmcSplitQuantitiesResult fmc_model_split_quantities_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in,
    const std::vector<double>& somatic_external,
    const std::vector<double>& quantiles_in);

// Fixed per-SNV variant-prior table indices: the Java model keeps a map over
// indel lengths -10..10 stored at [length + 10], with index 10 the SNV prior.
constexpr std::size_t kFmcPriorTable = 21;
constexpr std::size_t kFmcSnvPriorIndex = 10;

struct FmcEmIterationResult {
    std::vector<double> log_cluster_weights;
    std::vector<double> cluster_alpha;
    std::vector<double> cluster_beta;
    std::vector<double> log_variant_priors;  // size kFmcPriorTable
    double log_variant_vs_artifact_prior = 0.0;
    double somatic_variant_count = 0.0;
    double technical_artifact_count = 0.0;
    bool valid = false;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// One Java SomaticClusteringModel EM iteration, executed as a Kokkos work
// item over Views.  The GATK update is inherently sequential (responsibilities
// -> weight normalization -> optional prior update -> per-cluster gradient
// epochs over data), so a single device work item reproduces the exact
// operation order while every numeric primitive runs as a device kernel built
// on gamma_math.hpp.  `params`, `evidence` and `log_variant_priors_in` are
// inputs; updated model state and per-length variant counts are returned.
FmcEmIterationResult fmc_model_em_iteration_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in,
    const double log_variant_vs_artifact_prior_in,
    const double callable_sites,
    const bool update_somatic_priors);

}  // namespace fastgatk::kernels
