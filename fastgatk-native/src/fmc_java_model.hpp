// Host orchestration for the GATK 4.6.2.0 SomaticClusteringModel (from the
// machine-validated Python reference a5_p2_reference.py).  Used by the
// production FilterMutectCalls learning path.  Architecture split (1:1 goal): all learning
// numerics run as Kokkos kernels in fastgatk-kernels (filter_model_kokkos:
// error probabilities, split quantities for the quantile peak-split / BIC,
// responsibilities, EM iteration); this Host header only keeps the Java
// record()/pass orchestration, the quantile-grid / peak-scan control flow and
// scalar model bookkeeping.  Per-datum std-lib math was removed once each
// kernel covered it; remaining std::log uses are scalar weight arithmetic.
#pragma once

#include <algorithm>
#include <array>

#include <utility>
#include <vector>

#include "fastgatk/kernels/filter_model_kokkos.hpp"

namespace fastgatk_fmc_java {

constexpr double kLn10 = 2.3025850929940456840;
inline std::array<double, 2> fuzzy_binomial_shape(const double mean) {
    constexpr double std_over_mean = 0.01;
    const double bounded = std::min(mean, 1.0 - std_over_mean);
    const double alpha_plus_beta =
        (1.0 - bounded) / (bounded * std_over_mean * std_over_mean) - 1.0;
    const double alpha = bounded * alpha_plus_beta;
    return {alpha, alpha_plus_beta - alpha};
}

struct JavaDatum {
    double tumor_log_odds = 0.0;
    double artifact_prob = 0.0;
    double non_somatic_prob = 0.0;
    int alt = 0;
    int total = 0;
    int indel_length = 0;
};

struct Cluster {
    // kind 0 = beta-binomial, 1 = fuzzy binomial (beta shape with fixed
    // std/mean from the current mean).
    int kind = 0;
    double alpha = 1.0;
    double beta = 1.0;
    double mean = 0.5;  // used only for kind 1 refresh
};

struct SomaticModel {
    std::vector<Cluster> clusters;
    std::vector<double> log_cluster_weights;
    std::array<double, 21> log_variant_priors{};
    double log_variant_vs_artifact_prior = -1.0 * kLn10;
    bool initialized = false;
    int obvious_artifact_count = 0;
    // Transient accumulators owned by the caller.
    std::vector<JavaDatum> data;
    double callable_sites = 0.0;
    bool has_callable_sites = false;
    void record(const double tumor_log_odds, const double artifact_prob,
                const double non_somatic_prob, const int alt, const int total,
                const int indel_length) {
        if (artifact_prob > 0.9) {
            ++obvious_artifact_count;
            return;
        }
        if (non_somatic_prob > 0.9) return;
        data.push_back(JavaDatum{tumor_log_odds, artifact_prob,
                                 non_somatic_prob, alt, total, indel_length});
    }

    // Marshalling helpers shared by the kernel adapters below.  Numerics run
    // in fastgatk-kernels; these only copy model state into View inputs.
    fastgatk::kernels::FmcClusterParams kernel_params() const {
        fastgatk::kernels::FmcClusterParams params;
        params.log_cluster_weights.assign(clusters.size(), 0.0);
        for (std::size_t i = 0; i < clusters.size(); ++i)
            params.log_cluster_weights[i] = log_cluster_weights[i];
        params.cluster_alpha.assign(fastgatk::kernels::kFmcMaxClusters, 1.0);
        params.cluster_beta.assign(fastgatk::kernels::kFmcMaxClusters, 1.0);
        params.cluster_kind.assign(fastgatk::kernels::kFmcMaxClusters, 0);
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            params.cluster_alpha[i] = clusters[i].alpha;
            params.cluster_beta[i] = clusters[i].beta;
            params.cluster_kind[i] = clusters[i].kind;
        }
        return params;
    }

    fastgatk::kernels::FmcEvidence kernel_evidence() const {
        fastgatk::kernels::FmcEvidence evidence;
        evidence.alt_count.reserve(data.size());
        evidence.total_count.reserve(data.size());
        evidence.tumor_log_odds.reserve(data.size());
        evidence.artifact_prob.reserve(data.size());
        evidence.non_somatic_prob.reserve(data.size());
        for (const JavaDatum& datum : data) {
            evidence.alt_count.push_back(datum.alt);
            evidence.total_count.push_back(datum.total);
            evidence.tumor_log_odds.push_back(datum.tumor_log_odds);
            evidence.artifact_prob.push_back(datum.artifact_prob);
            evidence.non_somatic_prob.push_back(datum.non_somatic_prob);
            evidence.indel_length.push_back(datum.indel_length);
        }
        return evidence;
    }

    void initialize_clusters() {
        if (initialized) return;
        constexpr int kMaxBinomial = 5;
        constexpr int kNumQuantiles = 50;
        constexpr int kMinQuantileIndex = 5;
        constexpr double kMaxSplitFraction = 0.9;

        if (data.empty()) {
            // Mirrors the Java flow on an empty pass: no quantiles, no splits.
            initialized = true;
            return;
        }
        const fastgatk::kernels::FmcEvidence evidence = kernel_evidence();
        const std::vector<double> priors(log_variant_priors.begin(),
                                         log_variant_priors.end());
        // probabilityOfSomaticVariant snapshot at the pristine model, kept
        // fixed through the split/BIC flow exactly as Java computes it once.
        const auto pristine =
            fastgatk::kernels::fmc_model_split_quantities_kokkos(
                kernel_params(), evidence, priors, {}, {});
        if (!pristine.valid) return;
        const std::vector<double> somatic_probs =
            pristine.somatic_probability;

        double previous_bic = -std::numeric_limits<double>::infinity();
        for (int split = 0; split < kMaxBinomial; ++split) {
            const auto old_weights = log_cluster_weights;
            // Probability-weighted background responsibilities at quantiles.
            // The per-datum probabilityOfSomaticVariant is recomputed under
            // the current (iteration-start) model as Java does; every numeric
            // value comes from the shared kernels module.
            const auto sq =
                fastgatk::kernels::fmc_model_split_quantities_kokkos(
                    kernel_params(), evidence, priors, somatic_probs, {});
            if (!sq.valid) break;
            std::vector<std::pair<double, double>> af_somatic;
            af_somatic.reserve(data.size());
            for (std::size_t i = 0; i < data.size(); ++i)
                af_somatic.emplace_back(
                    static_cast<double>(data[i].alt) /
                        static_cast<double>(data[i].total),
                    sq.somatic_probability[i]);
            std::sort(af_somatic.begin(), af_somatic.end());
            double total_somatic = 0.0;
            for (const auto& entry : af_somatic) total_somatic += entry.second;
            const double step = total_somatic / kNumQuantiles;
            double cumulative = 0.0;
            double quantile_prob = step;
            std::vector<double> quantiles;
            for (const auto& entry : af_somatic) {
                cumulative += entry.second;
                if (cumulative > quantile_prob) {
                    quantiles.push_back(entry.first);
                    while (cumulative > quantile_prob) quantile_prob += step;
                }
            }
            if (quantiles.empty()) break;
            // Kernel-owned background responsibilities (softmax over the
            // count-based cluster likelihoods, weighted by the fixed somatic
            // probabilities) and the Java quantile background responsibilities.
            const auto bg =
                fastgatk::kernels::fmc_model_split_quantities_kokkos(
                    kernel_params(), evidence, priors, somatic_probs,
                    quantiles);
            if (!bg.valid) break;
            (void)bg.background_probability;  // fed into the qresp above
            const std::vector<double>& quantile_responsibilities =
                bg.quantile_responsibility;
            // Peaks (trapezoid quadrature between local minima).
            struct Peak {
                double fraction;
                double mass;
            };
            std::vector<Peak> peaks;
            double current_mass = 0.0;
            double current_peak = 0.0;
            double current_responsibility = 0.0;
            const std::size_t nq = quantiles.size();
            for (std::size_t q = 0; q < nq; ++q) {
                const double left_resp =
                    q == 0 ? 0.0 : quantile_responsibilities[q - 1];
                const double responsibility = quantile_responsibilities[q];
                const double right_resp = q == nq - 1
                    ? 0.0 : quantile_responsibilities[q + 1];
                const double left_fraction =
                    q == 0 ? 0.0 : quantiles[q - 1];
                current_mass += (quantiles[q] - left_fraction) *
                    (left_resp + responsibility) / 2.0;
                if (responsibility > current_responsibility) {
                    current_peak = quantiles[q];
                    current_responsibility = responsibility;
                }
                const bool local_min =
                    (responsibility < left_resp && responsibility <= right_resp) ||
                    (responsibility <= left_resp && responsibility < right_resp);
                if ((local_min && q > 0) || q == nq - 1) {
                    peaks.push_back({current_peak, current_mass});
                    current_mass = 0.0;
                    current_peak = quantiles[q];
                    current_responsibility = responsibility;
                }
            }
            if (peaks.empty()) break;
            Peak* biggest = &peaks[0];
            for (Peak& peak : peaks)
                if (peak.mass > biggest->mass) biggest = &peak;
            if (biggest->fraction <
                quantiles[std::min<std::size_t>(
                    kMinQuantileIndex, quantiles.size() - 1)])
                break;
            double total_mass = 0.0;
            for (const Peak& peak : peaks) total_mass += peak.mass;
            const double fraction = std::min(
                kMaxSplitFraction,
                total_mass > 0.0 ? biggest->mass / total_mass : 0.0);
            if (fraction <= 0.0) break;
            const auto shape = fuzzy_binomial_shape(biggest->fraction);
            clusters.push_back(Cluster{1, shape[0], shape[1],
                                       biggest->fraction});
            log_cluster_weights.push_back(std::log(fraction) +
                                          log_cluster_weights[0]);
            log_cluster_weights[0] = std::log1p(fraction) +
                                     log_cluster_weights[0];
            // Java runs NUM_ITERATIONS = 5 EM(false) iterations per accepted
            // split before the BIC comparison; this drives the background
            // beta to its GATK value (a single iteration lagged bg beta and
            // produced the old ~0.00035 threshold drift).
            for (int init_em = 0; init_em < 5; ++init_em)
                perform_em_iteration(false);
            // BIC under the refreshed (EM-updated) model: kernel-owned mixture
            // log likelihood per datum weighted by the fixed somatic probs.
            const auto bic_sq =
                fastgatk::kernels::fmc_model_split_quantities_kokkos(
                    kernel_params(), evidence, priors, somatic_probs, {});
            if (!bic_sq.valid) break;
            double weighted = 0.0;
            for (std::size_t i = 0; i < data.size(); ++i)
                weighted += somatic_probs[i] *
                            bic_sq.log_somatic_mixture[i];
            double effective = 0.0;
            for (const double prob : somatic_probs) effective += prob;
            const double num_params = 2.0 * clusters.size();
            const double bic = weighted - num_params * std::log(effective);
            if (bic < previous_bic) {
                clusters.pop_back();
                log_cluster_weights = old_weights;
                break;
            }
            previous_bic = bic;
        }
        initialized = true;
    }

    void perform_em_iteration(const bool update_somatic_priors) {
        // The EM numerics run as a Kokkos kernel in fastgatk-kernels; this
        // Host method only marshals state into Views and applies the result.
        const std::vector<double> priors(log_variant_priors.begin(),
                                         log_variant_priors.end());
        const double callable = has_callable_sites ? callable_sites : 0.0;
        const auto result = fastgatk::kernels::fmc_model_em_iteration_kokkos(
            kernel_params(), kernel_evidence(), priors,
            log_variant_vs_artifact_prior, callable,
            update_somatic_priors && has_callable_sites);
        // The kernel is a no-op for an empty pass (no recorded datums) and
        // returns empty result vectors; keep the pristine model untouched.
        if (!result.valid || data.empty()) return;

        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (i < result.log_cluster_weights.size())
                log_cluster_weights[i] = result.log_cluster_weights[i];
            clusters[i].alpha = result.cluster_alpha[i];
            clusters[i].beta = result.cluster_beta[i];
        }
        for (std::size_t i = 0; i < log_variant_priors.size(); ++i)
            log_variant_priors[i] = result.log_variant_priors[i];
        log_variant_vs_artifact_prior =
            result.log_variant_vs_artifact_prior;
        if (std::getenv("FASTGATK_DEBUG_FMC") != nullptr) {
            std::fprintf(stderr,
                "[CXXEM] upd=%d k=%zu w0=%.17g bg=%.17g/%.17g p10=%.17g\n",
                update_somatic_priors ? 1 : 0, clusters.size(),
                log_cluster_weights.empty() ? 0.0 : log_cluster_weights[0],
                clusters[0].alpha, clusters[0].beta,
                log_variant_priors[10]);
        }
    }

    void learn_and_clear() {
        if (!initialized) initialize_clusters();
        for (int iteration = 0; iteration < 5; ++iteration)
            perform_em_iteration(true);
        data.clear();
        obvious_artifact_count = 0;
    }
};

inline double optimal_f_score_threshold(std::vector<double> posteriors,
                                        const double beta = 1.0) {
    std::sort(posteriors.begin(), posteriors.end());
    double expected_tp = 0.0;
    for (const double posterior : posteriors) expected_tp += 1.0 - posterior;
    double tp = 0.0;
    double fp = 0.0;
    double fn = expected_tp;
    int optimal_index = -1;
    double optimal_f = 0.0;
    for (std::size_t n = 0; n < posteriors.size(); ++n) {
        tp += 1.0 - posteriors[n];
        fp += posteriors[n];
        fn -= 1.0 - posteriors[n];
        const double f = (1.0 + beta * beta) * tp /
            ((1.0 + beta * beta) * tp + beta * beta * fn + fp);
        if (f >= optimal_f) {
            optimal_index = static_cast<int>(n);
            optimal_f = f;
        }
    }
    if (optimal_index == -1) return 0.0;
    if (optimal_index == static_cast<int>(posteriors.size()) - 1) return 1.0;
    return posteriors[static_cast<std::size_t>(optimal_index)];
}

}  // namespace fastgatk_fmc_java
