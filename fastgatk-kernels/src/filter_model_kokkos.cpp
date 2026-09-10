#include "fastgatk/kernels/filter_model_kokkos.hpp"
#include "fastgatk/kernels/gamma_math.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

// Commons/GATK log beta-binomial PMF (log choose + logBeta terms) and log
// binomial PMF built on the shared Commons-mirror log-gamma kernels.  The
// shared `log_choose` cancels in most odds ratios but is kept for parity.
KOKKOS_INLINE_FUNCTION double log_choose(const int n, const int k) {
    return detail::log_gamma_positive(n + 1.0) -
           detail::log_gamma_positive(k + 1.0) -
           detail::log_gamma_positive(n - k + 1.0);
}

KOKKOS_INLINE_FUNCTION double log_beta_binom_pmf(const int n, const int k,
                                                 const double alpha,
                                                 const double beta) {
    return log_choose(n, k) +
           detail::log_dirichlet_normalization_two(alpha, beta) -
           detail::log_dirichlet_normalization_two(alpha + k, beta + n - k);
}

KOKKOS_INLINE_FUNCTION double log_binom_pmf(const int n, const int k,
                                            const double p) {
    if (!(p > 0.0)) return k == 0 ? 0.0 : -std::numeric_limits<double>::infinity();
    if (p >= 1.0)
        return k == n ? 0.0 : -std::numeric_limits<double>::infinity();
    return log_choose(n, k) + k * Kokkos::log(p) +
           (n - k) * Kokkos::log1p(-p);
}

KOKKOS_INLINE_FUNCTION double log_sum_exp_pair(const double a,
                                               const double b) {
    const double high = a > b ? a : b;
    if (high == -std::numeric_limits<double>::infinity()) return high;
    return high + Kokkos::log(Kokkos::exp(a - high) + Kokkos::exp(b - high));
}

KOKKOS_INLINE_FUNCTION double log_sum_exp(const double* values,
                                          const int count) {
    double high = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < count; ++i) high = high > values[i] ? high : values[i];
    if (high == -std::numeric_limits<double>::infinity()) return high;
    double total = 0.0;
    for (int i = 0; i < count; ++i) total += Kokkos::exp(values[i] - high);
    return high + Kokkos::log(total);
}

// Java correctedLogLikelihood: tumor log odds + log-odds correction measured
// against the flat Beta(1,1) reference.
KOKKOS_INLINE_FUNCTION double corrected_log_likelihood(
    const double tumor_log_odds, const int alt, const int ref,
    const double alpha, const double beta) {
    const double g_new = detail::log_dirichlet_normalization_two(alpha, beta);
    const double g_new_counts =
        detail::log_dirichlet_normalization_two(alpha + alt, beta + ref);
    const double g_flat_counts = detail::log_dirichlet_normalization_two(
        1.0 + alt, 1.0 + ref);
    return tumor_log_odds + g_new - g_new_counts + g_flat_counts;
}

}  // namespace

KOKKOS_INLINE_FUNCTION double log1mexp_fmc(const double x) {
    if (x >= 0.0) return -std::numeric_limits<double>::infinity();
    return Kokkos::log(-Kokkos::expm1(x));
}

// Softmax over two natural-log odds: posterior error probability for
// logOdds-of-real vs error with logPriorReal.
KOKKOS_INLINE_FUNCTION double posterior_error_fmc(
    const double log_odds_real_vs_error, const double log_prior_real) {
    const double real = log_odds_real_vs_error + log_prior_real;
    const double error = log1mexp_fmc(log_prior_real);
    const double high = real > error ? real : error;
    const double total = Kokkos::exp(real - high) + Kokkos::exp(error - high);
    return Kokkos::exp(error - high) / total;
}

// Java probabilityOfSequencingError under the current cluster model: softmax
// posterior error of the corrected-log-odds mixture with the per-length
// variant prior (SNV adds log(1/3)).  Shared by the error-probabilities, split
// and EM kernels so the per-datum sequencing posterior has one numeric source.
KOKKOS_INLINE_FUNCTION double fmc_seq_posterior_device(
    const double* logw, const double* alpha, const double* beta,
    const int ccount, const double tumor_log_odds, const int alt,
    const int total, const int indel, const double* priors) {
    double corrected[kFmcMaxClusters];
    for (int c = 0; c < ccount; ++c)
        corrected[c] = logw[c] + corrected_log_likelihood(
            tumor_log_odds, alt, total - alt, alpha[c], beta[c]);
    const std::size_t prior_index = static_cast<std::size_t>(
        (indel + 10 < 0) ? 0 : (indel + 10 > 20 ? 20 : indel + 10));
    double seq_prior = priors[prior_index];
    if (indel == 0) seq_prior += Kokkos::log(1.0 / 3.0);
    return posterior_error_fmc(log_sum_exp(corrected, ccount), seq_prior);
}

// Java logLikelihoodGivenSomatic: logsumexp over the weighted count-based
// beta-binomial cluster likelihoods.  Shared numeric source for the
// background-given-somatic softmax, contamination posterior and BIC mixture.
KOKKOS_INLINE_FUNCTION double fmc_count_mixture_device(
    const double* logw, const double* alpha, const double* beta,
    const int ccount, const int total, const int alt) {
    double values[kFmcMaxClusters];
    for (int c = 0; c < ccount; ++c)
        values[c] = logw[c] + log_beta_binom_pmf(total, alt, alpha[c], beta[c]);
    return log_sum_exp(values, ccount);
}

FmcEmIterationResult fmc_model_em_iteration_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in,
    const double log_variant_vs_artifact_prior_in,
    const double callable_sites,
    const bool update_somatic_priors) {
    (void)log_variant_vs_artifact_prior_in;  // recomputed from counts when
                                             // update_somatic_priors is set.
    FmcEmIterationResult result;
    const std::size_t cluster_count = params.log_cluster_weights.size();
    const std::size_t data_count = evidence.alt_count.size();
    if (cluster_count == 0 || cluster_count > kFmcMaxClusters ||
        data_count == 0 || evidence.total_count.size() != data_count ||
        evidence.tumor_log_odds.size() != data_count ||
        evidence.artifact_prob.size() != data_count ||
        evidence.non_somatic_prob.size() != data_count ||
        (!evidence.indel_length.empty() &&
         evidence.indel_length.size() != data_count) ||
        log_variant_priors_in.size() != kFmcPriorTable ||
        params.cluster_alpha.size() != kFmcMaxClusters ||
        params.cluster_beta.size() != kFmcMaxClusters ||
        params.cluster_kind.size() != kFmcMaxClusters)
        return result;
    if (!Kokkos::is_initialized()) return result;

    const auto prepare_begin = std::chrono::steady_clock::now();
    using DView = Kokkos::View<double*, MemorySpace>;
    using KView = Kokkos::View<int*, MemorySpace>;
    DView device_logw("fmc_em_logw", kFmcMaxClusters);
    DView device_alpha("fmc_em_alpha", kFmcMaxClusters);
    DView device_beta("fmc_em_beta", kFmcMaxClusters);
    KView device_kind("fmc_em_kind", kFmcMaxClusters);
    DView device_tlo("fmc_em_tlo", data_count);
    DView device_art("fmc_em_art", data_count);
    DView device_ns("fmc_em_ns", data_count);
    KView device_alt("fmc_em_alt", data_count);
    KView device_total("fmc_em_total", data_count);
    KView device_indel("fmc_em_indel", data_count);
    DView device_priors("fmc_em_priors", kFmcPriorTable);
    DView device_resp("fmc_em_resp", data_count * kFmcMaxClusters);
    DView device_somatic("fmc_em_somatic", data_count);
    DView device_out_logw("fmc_em_out_logw", kFmcMaxClusters);
    DView device_out_alpha("fmc_em_out_alpha", kFmcMaxClusters);
    DView device_out_beta("fmc_em_out_beta", kFmcMaxClusters);
    DView device_out_priors("fmc_em_out_priors", kFmcPriorTable);
    DView device_counts("fmc_em_counts", kFmcPriorTable);
    DView device_scalar("fmc_em_scalar", 3);
    auto host_logw = Kokkos::create_mirror_view(device_logw);
    auto host_alpha = Kokkos::create_mirror_view(device_alpha);
    auto host_beta = Kokkos::create_mirror_view(device_beta);
    auto host_tlo = Kokkos::create_mirror_view(device_tlo);
    auto host_art = Kokkos::create_mirror_view(device_art);
    auto host_ns = Kokkos::create_mirror_view(device_ns);
    auto host_alt = Kokkos::create_mirror_view(device_alt);
    auto host_total = Kokkos::create_mirror_view(device_total);
    auto host_indel = Kokkos::create_mirror_view(device_indel);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    auto host_out_logw = Kokkos::create_mirror_view(device_out_logw);
    auto host_out_alpha = Kokkos::create_mirror_view(device_out_alpha);
    auto host_out_beta = Kokkos::create_mirror_view(device_out_beta);
    auto host_out_priors = Kokkos::create_mirror_view(device_out_priors);
    auto host_counts = Kokkos::create_mirror_view(device_counts);
    auto host_scalar = Kokkos::create_mirror_view(device_scalar);
    auto host_kind = Kokkos::create_mirror_view(device_kind);
    for (std::size_t i = 0; i < kFmcMaxClusters; ++i) {
        host_logw(i) = i < cluster_count ? params.log_cluster_weights[i]
                                         : -std::numeric_limits<double>::infinity();
        host_alpha(i) = params.cluster_alpha[i];
        host_beta(i) = params.cluster_beta[i];
        host_kind(i) = params.cluster_kind[i];
    }
    const bool has_indel = !evidence.indel_length.empty();
    for (std::size_t i = 0; i < data_count; ++i) {
        host_alt(i) = evidence.alt_count[i];
        host_total(i) = evidence.total_count[i];
        host_tlo(i) = evidence.tumor_log_odds[i];
        host_art(i) = evidence.artifact_prob[i];
        host_ns(i) = evidence.non_somatic_prob[i];
        host_indel(i) = has_indel ? evidence.indel_length[i] : 0;
    }
    for (std::size_t i = 0; i < kFmcPriorTable; ++i)
        host_priors(i) = log_variant_priors_in[i];
    Kokkos::deep_copy(device_logw, host_logw);
    Kokkos::deep_copy(device_alpha, host_alpha);
    Kokkos::deep_copy(device_beta, host_beta);
    Kokkos::deep_copy(device_tlo, host_tlo);
    Kokkos::deep_copy(device_art, host_art);
    Kokkos::deep_copy(device_ns, host_ns);
    Kokkos::deep_copy(device_alt, host_alt);
    Kokkos::deep_copy(device_total, host_total);
    Kokkos::deep_copy(device_indel, host_indel);
    Kokkos::deep_copy(device_priors, host_priors);
    Kokkos::deep_copy(device_kind, host_kind);
    Kokkos::deep_copy(device_resp, 0.0);
    Kokkos::deep_copy(device_somatic, 0.0);
    Kokkos::deep_copy(device_counts, 0.0);
    Kokkos::deep_copy(device_scalar, 0.0);
    const double callable = callable_sites;
    const bool update_priors = update_somatic_priors;
    const auto prepare_end = std::chrono::steady_clock::now();
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_fmc_em_iteration",
        Kokkos::RangePolicy<ExecSpace>(0, 1),
        KOKKOS_LAMBDA(const std::size_t) {
            const int ccount = static_cast<int>(cluster_count);
            // Pass 1: somatic weight + scaled responsibilities snapshot.
            for (std::size_t datum = 0; datum < data_count; ++datum) {
                const int alt = device_alt(datum);
                const int total = device_total(datum);
                if (total <= 0 || alt < 0 || alt > total) continue;
                const double sequencing_error = fmc_seq_posterior_device(
                    device_logw.data(), device_alpha.data(),
                    device_beta.data(), ccount, device_tlo(datum), alt,
                    total, device_indel(datum), device_priors.data());
                const double somatic = (1.0 - device_art(datum)) *
                    (1.0 - device_ns(datum)) * (1.0 - sequencing_error);
                device_somatic(datum) = somatic;
                const int indel2 = device_indel(datum);
                const std::size_t count_index = static_cast<std::size_t>(
                    (indel2 + 10 < 0) ? 0
                    : (indel2 + 10 > 20 ? 20 : indel2 + 10));
                device_counts(count_index) += somatic;
                device_scalar(0) += somatic;
                device_scalar(1) += device_art(datum);
                double values[kFmcMaxClusters];
                for (int c = 0; c < ccount; ++c) {
                    values[c] = device_logw(c) + log_beta_binom_pmf(
                        total, alt, device_alpha(c), device_beta(c));
                }
                const double high = log_sum_exp(values, ccount);
                for (int c = 0; c < ccount; ++c) {
                    const double responsibility = somatic *
                        Kokkos::exp(values[c] - high);
                    device_resp(datum * kFmcMaxClusters +
                                static_cast<std::size_t>(c)) = responsibility;
                }
            }
            // Weight normalization with REGULARIZING_PSEUDOCOUNT = 1.
            for (int c = 0; c < ccount; ++c) {
                double cluster_total = 1.0;
                for (std::size_t datum = 0; datum < data_count; ++datum)
                    cluster_total += device_resp(
                        datum * kFmcMaxClusters + static_cast<std::size_t>(c));
                device_out_logw(c) = Kokkos::log(cluster_total);
            }
            double weight_denom = 0.0;
            for (int c = 0; c < ccount; ++c)
                weight_denom += Kokkos::exp(device_out_logw(c));
            for (int c = 0; c < ccount; ++c)
                device_out_logw(c) -= Kokkos::log(weight_denom);
            if (update_priors) {
                const double variant_count = device_scalar(0);
                const double artifact_count = device_scalar(1);
                device_scalar(2) = Kokkos::log(
                    (variant_count + 1.0) /
                    (variant_count + artifact_count + 2.0));
                if (callable >= 1.0) {
                    for (int length = -10; length <= 10; ++length) {
                        const std::size_t index =
                            static_cast<std::size_t>(length + 10);
                        const double ratio =
                            device_counts(index) / callable;
                        const double floor =
                            length == 0 ? 1.0e-8 : 1.0e-9;
                        device_out_priors(index) =
                            Kokkos::log(ratio > floor ? ratio : floor);
                    }
                } else {
                    for (std::size_t i = 0; i < kFmcPriorTable; ++i)
                        device_out_priors(i) = device_priors(i);
                }
            } else {
                for (std::size_t i = 0; i < kFmcPriorTable; ++i)
                    device_out_priors(i) = device_priors(i);
            }
            // Cluster learning uses the pass-1 responsibility snapshot.
            for (int c = 0; c < ccount; ++c) {
                if (device_kind(c) == 1) {
                    double alt_sum = 0.0;
                    double total_sum = 0.0;
                    for (std::size_t datum = 0; datum < data_count; ++datum) {
                        const int alt = device_alt(datum);
                        const int total = device_total(datum);
                        if (total <= 0 || alt < 0 || alt > total) continue;
                        const double responsibility = device_resp(
                            datum * kFmcMaxClusters +
                            static_cast<std::size_t>(c));
                        alt_sum += alt * responsibility;
                        total_sum += total * responsibility;
                    }
                    const double mean = (alt_sum + 1.0e-4) /
                                        (total_sum + 1.0e-4);
                    const double bounded =
                        mean < 1.0 - 0.01 ? mean : 1.0 - 0.01;
                    const double alpha_plus_beta =
                        (1.0 - bounded) / (bounded * 0.01 * 0.01) - 1.0;
                    const double alpha = bounded * alpha_plus_beta;
                    device_out_alpha(c) = alpha;
                    device_out_beta(c) = alpha_plus_beta - alpha;
                } else {
                    double alpha = device_alpha(c);
                    double beta = device_beta(c);
                    constexpr double rate = 0.01;
                    for (int epoch = 0; epoch < 10; ++epoch) {
                        for (std::size_t datum = 0;
                             datum < data_count; ++datum) {
                            const int alt = device_alt(datum);
                            const int total = device_total(datum);
                            if (total <= 0 || alt < 0 || alt > total) continue;
                            const double responsibility = device_resp(
                                datum * kFmcMaxClusters +
                                static_cast<std::size_t>(c));
                            if (!(responsibility > 0.0)) continue;
                            const int ref = total - alt;
                            const double digamma_total =
                                detail::digamma_positive(total + alpha + beta);
                            const double digamma_sum =
                                detail::digamma_positive(alpha + beta);
                            const double alpha_gradient =
                                detail::digamma_positive(alpha + alt) -
                                digamma_total -
                                detail::digamma_positive(alpha) + digamma_sum;
                            const double beta_gradient =
                                detail::digamma_positive(beta + ref) -
                                digamma_total -
                                detail::digamma_positive(beta) + digamma_sum;
                            alpha += rate * alpha_gradient * responsibility;
                            beta += rate * beta_gradient * responsibility;
                            if (alpha < 1.0) alpha = 1.0;
                            if (beta < 0.5) beta = 0.5;
                        }
                    }
                    device_out_alpha(c) = alpha;
                    device_out_beta(c) = beta;
                }
            }
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();

    Kokkos::deep_copy(host_out_logw, device_out_logw);
    Kokkos::deep_copy(host_out_alpha, device_out_alpha);
    Kokkos::deep_copy(host_out_beta, device_out_beta);
    Kokkos::deep_copy(host_out_priors, device_out_priors);
    Kokkos::deep_copy(host_counts, device_counts);
    Kokkos::deep_copy(host_scalar, device_scalar);
    result.log_cluster_weights.assign(cluster_count, 0.0);
    result.cluster_alpha.assign(kFmcMaxClusters, 0.0);
    result.cluster_beta.assign(kFmcMaxClusters, 0.0);
    result.log_variant_priors.assign(kFmcPriorTable, 0.0);
    for (std::size_t i = 0; i < cluster_count; ++i)
        result.log_cluster_weights[i] = host_out_logw(i);
    for (std::size_t i = 0; i < kFmcMaxClusters; ++i) {
        result.cluster_alpha[i] = host_out_alpha(i);
        result.cluster_beta[i] = host_out_beta(i);
    }
    for (std::size_t i = 0; i < kFmcPriorTable; ++i)
        result.log_variant_priors[i] = host_out_priors(i);
    result.somatic_variant_count = host_scalar(0);
    result.technical_artifact_count = host_scalar(1);
    result.log_variant_vs_artifact_prior = host_scalar(2);
    result.valid = true;
    return result;
}

FmcErrorProbabilitiesResult fmc_error_probabilities_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in) {
    FmcErrorProbabilitiesResult result;
    const std::size_t cluster_count = params.log_cluster_weights.size();
    const std::size_t data_count = evidence.alt_count.size();
    if (cluster_count == 0 || cluster_count > kFmcMaxClusters ||
        data_count == 0 || evidence.total_count.size() != data_count ||
        evidence.alt_count.size() != data_count ||
        evidence.tumor_log_odds.size() != data_count ||
        evidence.non_somatic_prob.size() != data_count ||
        (!evidence.indel_length.empty() &&
         evidence.indel_length.size() != data_count) ||
        (!evidence.normal_log10_odds.empty() &&
         evidence.normal_log10_odds.size() != data_count) ||
        (!evidence.germline_evaluated.empty() &&
         evidence.germline_evaluated.size() != data_count) ||
        log_variant_priors_in.size() != kFmcPriorTable ||
        params.cluster_alpha.size() != kFmcMaxClusters ||
        params.cluster_beta.size() != kFmcMaxClusters ||
        params.cluster_kind.size() != kFmcMaxClusters)
        return result;
    result.sequencing_error.assign(data_count, 0.0);
    result.non_somatic_error.assign(data_count, 0.0);
    result.germline_error.assign(data_count, 0.0);
    result.contamination_error.assign(data_count, 0.0);
    if (!Kokkos::is_initialized()) return result;

    const auto prepare_begin = std::chrono::steady_clock::now();
    using DView = Kokkos::View<double*, MemorySpace>;
    using KView = Kokkos::View<int*, MemorySpace>;
    DView device_logw("fmc_ep_logw", cluster_count);
    DView device_alpha("fmc_ep_alpha", kFmcMaxClusters);
    DView device_beta("fmc_ep_beta", kFmcMaxClusters);
    KView device_alt("fmc_ep_alt", data_count);
    KView device_total("fmc_ep_total", data_count);
    KView device_indel("fmc_ep_indel", data_count);
    KView device_gate("fmc_ep_gate", data_count);
    DView device_tlo("fmc_ep_tlo", data_count);
    DView device_af("fmc_ep_af", data_count);
    DView device_cont("fmc_ep_cont", data_count);
    DView device_nlod("fmc_ep_nlod", data_count);
    DView device_priors("fmc_ep_priors", kFmcPriorTable);
    DView device_seq("fmc_ep_seq", data_count);
    DView device_ns("fmc_ep_ns", data_count);
    DView device_germ("fmc_ep_germ", data_count);
    DView device_cont_out("fmc_ep_cont_out", data_count);
    auto host_logw = Kokkos::create_mirror_view(device_logw);
    auto host_alpha = Kokkos::create_mirror_view(device_alpha);
    auto host_beta = Kokkos::create_mirror_view(device_beta);
    auto host_alt = Kokkos::create_mirror_view(device_alt);
    auto host_total = Kokkos::create_mirror_view(device_total);
    auto host_indel = Kokkos::create_mirror_view(device_indel);
    auto host_gate = Kokkos::create_mirror_view(device_gate);
    auto host_tlo = Kokkos::create_mirror_view(device_tlo);
    auto host_af = Kokkos::create_mirror_view(device_af);
    auto host_cont = Kokkos::create_mirror_view(device_cont);
    auto host_nlod = Kokkos::create_mirror_view(device_nlod);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    auto host_seq = Kokkos::create_mirror_view(device_seq);
    auto host_ns = Kokkos::create_mirror_view(device_ns);
    auto host_germ = Kokkos::create_mirror_view(device_germ);
    auto host_cont_out = Kokkos::create_mirror_view(device_cont_out);
    for (std::size_t i = 0; i < cluster_count; ++i)
        host_logw(i) = params.log_cluster_weights[i];
    for (std::size_t i = 0; i < kFmcMaxClusters; ++i) {
        host_alpha(i) = params.cluster_alpha[i];
        host_beta(i) = params.cluster_beta[i];
    }
    const bool has_indel = !evidence.indel_length.empty();
    const bool has_af = !evidence.population_af.empty();
    const bool has_cont = !evidence.contamination.empty();
    const bool has_nlod = !evidence.normal_log10_odds.empty();
    const bool has_gate = !evidence.germline_evaluated.empty();
    for (std::size_t i = 0; i < data_count; ++i) {
        host_alt(i) = evidence.alt_count[i];
        host_total(i) = evidence.total_count[i];
        host_tlo(i) = evidence.tumor_log_odds[i];
        host_indel(i) = has_indel ? evidence.indel_length[i] : 0;
        host_gate(i) = has_gate ? evidence.germline_evaluated[i] : 1;
        // Host defaults: population AF 0.01, contamination estimate 0.1.
        host_af(i) = has_af ? evidence.population_af[i] : 0.01;
        host_cont(i) = has_cont ? evidence.contamination[i] : 0.1;
        host_nlod(i) = has_nlod ? evidence.normal_log10_odds[i] : 0.0;
    }
    for (std::size_t i = 0; i < kFmcPriorTable; ++i)
        host_priors(i) = log_variant_priors_in[i];
    Kokkos::deep_copy(device_logw, host_logw);
    Kokkos::deep_copy(device_alpha, host_alpha);
    Kokkos::deep_copy(device_beta, host_beta);
    Kokkos::deep_copy(device_alt, host_alt);
    Kokkos::deep_copy(device_total, host_total);
    Kokkos::deep_copy(device_indel, host_indel);
    Kokkos::deep_copy(device_gate, host_gate);
    Kokkos::deep_copy(device_tlo, host_tlo);
    Kokkos::deep_copy(device_af, host_af);
    Kokkos::deep_copy(device_cont, host_cont);
    Kokkos::deep_copy(device_nlod, host_nlod);
    Kokkos::deep_copy(device_priors, host_priors);
    const auto prepare_end = std::chrono::steady_clock::now();
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();

    const auto execute_begin = std::chrono::steady_clock::now();
    // One thread per datum; every per-datum numeric below is the device
    // transcription of the Java / host formulas (SomaticClusteringModel
    // sequencing posterior, GermlineFilter het germline posterior, and
    // ContaminationFilter contaminant posterior).
    Kokkos::parallel_for(
        "fastgatk_fmc_error_probabilities",
        Kokkos::RangePolicy<ExecSpace>(0, data_count),
        KOKKOS_LAMBDA(const std::size_t datum) {
            const int ccount = static_cast<int>(cluster_count);
            const int alt = device_alt(datum);
            const int total = device_total(datum);
            if (total <= 0 || alt < 0 || alt > total) return;
            const int indel = device_indel(datum);
            const std::size_t prior_index = static_cast<std::size_t>(
                (indel + 10 < 0) ? 0 : (indel + 10 > 20 ? 20 : indel + 10));
            double log_prior = device_priors(prior_index);
            if (indel == 0) log_prior += Kokkos::log(1.0 / 3.0);
            const double af = device_af(datum);
            const double cont = device_cont(datum);
            // NLOD is log10 odds of absence from the normal; GermlineFilter
            // needs log odds of presence, hence the sign reversal.
            const double normal_log_odds = -device_nlod(datum) * Kokkos::log(10.0);

            // Sequencing posterior via the shared device helper.
            const double seq = fmc_seq_posterior_device(
                device_logw.data(), device_alpha.data(), device_beta.data(),
                ccount, device_tlo(datum), alt, total, indel,
                device_priors.data());
            device_seq(datum) = seq;

            // Count-based somatic mixture (logLikelihoodGivenSomatic) via the
            // shared device helper, used by germline and contamination.
            const double log_somatic = fmc_count_mixture_device(
                device_logw.data(), device_alpha.data(), device_beta.data(),
                ccount, total, alt);

            // Germline posterior: diploid het at maf 0.5 plus the GATK
            // high-AF homozygous-alt branch.  Normal evidence is carried by
            // NLOD, whose sign convention is opposite that expected by
            // GermlineFilter.germlineProbability().
            // Skipped when the site carried no POPAF annotation, mirroring
            // GATK's GermlineFilter which has no hypothesis to evaluate then.
            double germ = 0.0;
            if (device_gate(datum) != 0) {
                const double pmf_het = log_binom_pmf(total, alt, 0.5);
                const double log_germ_lik =
                    -Kokkos::log(2.0) +
                    log_sum_exp_pair(pmf_het, pmf_het);
                const double lodds_het = log_germ_lik - log_somatic;
                const double lpns = log1mexp_fmc(log_prior);
                const double lpg_het =
                    Kokkos::log(2.0 * af * (1.0 - af)) + lodds_het +
                    normal_log_odds + lpns;
                const double lpg_hom = Kokkos::log(af * af) +
                    (static_cast<double>(alt) / static_cast<double>(total) < 0.9
                        ? -std::numeric_limits<double>::infinity() : 0.0) +
                    normal_log_odds + lpns;
                const double lpg = log_sum_exp_pair(lpg_het, lpg_hom);
                const double lps = Kokkos::log((1.0 - af) * (1.0 - af)) +
                                   log_prior;
                // softmax({lpg, lps})[0].
                const double high = lpg > lps ? lpg : lps;
                const double total2 =
                    Kokkos::exp(lpg - high) + Kokkos::exp(lps - high);
                germ = Kokkos::exp(lpg - high) / total2;
            }

            // Contamination posterior (contamination fraction `cont`,
            // population AF `af`, count-based somatic likelihood).
            const double single =
                2.0 * af * (1.0 - af) *
                    Kokkos::exp(log_binom_pmf(total, alt, cont / 2.0)) +
                af * af *
                    Kokkos::exp(log_binom_pmf(total, alt, cont));
            const double many = Kokkos::exp(
                log_binom_pmf(total, alt, cont * af));
            // Degenerate no-contamination case (estimate 0): no contaminant
            // mass at any hypothesis, so the contamination posterior is 0
            // (limit of posteriorError(+inf odds, prior)); -inf -> exp of
            // (inf-inf) would produce NaN in the softmax below.
            double cont_p = 0.0;
            if (single > 0.0 || many > 0.0) {
                const double log_contaminant =
                    Kokkos::log(single > many ? single : many);
                cont_p = posterior_error_fmc(
                    log_somatic - log_contaminant, log_prior);
            }
            device_germ(datum) = germ;
            device_cont_out(datum) = cont_p;
            device_ns(datum) = germ > cont_p ? germ : cont_p;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    Kokkos::deep_copy(host_seq, device_seq);
    Kokkos::deep_copy(host_ns, device_ns);
    Kokkos::deep_copy(host_germ, device_germ);
    Kokkos::deep_copy(host_cont_out, device_cont_out);
    for (std::size_t i = 0; i < data_count; ++i) {
        result.sequencing_error[i] = host_seq(i);
        result.non_somatic_error[i] = host_ns(i);
        result.germline_error[i] = host_germ(i);
        result.contamination_error[i] = host_cont_out(i);
    }
    result.valid = true;
    return result;
}

FmcSplitQuantitiesResult fmc_model_split_quantities_kokkos(
    const FmcClusterParams& params,
    const FmcEvidence& evidence,
    const std::vector<double>& log_variant_priors_in,
    const std::vector<double>& somatic_external,
    const std::vector<double>& quantiles_in) {
    FmcSplitQuantitiesResult result;
    const std::size_t cluster_count = params.log_cluster_weights.size();
    const std::size_t data_count = evidence.alt_count.size();
    const bool use_external = !somatic_external.empty();
    if (use_external && somatic_external.size() != data_count) return result;
    if (cluster_count == 0 || cluster_count > kFmcMaxClusters ||
        data_count == 0 || evidence.total_count.size() != data_count ||
        evidence.tumor_log_odds.size() != data_count ||
        evidence.artifact_prob.size() != data_count ||
        evidence.non_somatic_prob.size() != data_count ||
        (!evidence.indel_length.empty() &&
         evidence.indel_length.size() != data_count) ||
        log_variant_priors_in.size() != kFmcPriorTable ||
        params.cluster_alpha.size() != kFmcMaxClusters ||
        params.cluster_beta.size() != kFmcMaxClusters ||
        params.cluster_kind.size() != kFmcMaxClusters)
        return result;
    const std::size_t q_count = quantiles_in.size();
    result.sequencing_error.assign(data_count, 0.0);
    result.somatic_probability.assign(data_count, 0.0);
    result.log_somatic_mixture.assign(data_count, 0.0);
    result.background_given_somatic.assign(data_count, 0.0);
    result.background_probability.assign(data_count, 0.0);
    result.quantile_responsibility.assign(q_count, 0.0);
    if (!Kokkos::is_initialized()) return result;

    const auto prepare_begin = std::chrono::steady_clock::now();
    using DView = Kokkos::View<double*, MemorySpace>;
    using KView = Kokkos::View<int*, MemorySpace>;
    DView device_logw("fmc_sq_logw", cluster_count);
    DView device_alpha("fmc_sq_alpha", kFmcMaxClusters);
    DView device_beta("fmc_sq_beta", kFmcMaxClusters);
    KView device_kind("fmc_sq_kind", kFmcMaxClusters);
    KView device_alt("fmc_sq_alt", data_count);
    KView device_total("fmc_sq_total", data_count);
    KView device_indel("fmc_sq_indel", data_count);
    DView device_tlo("fmc_sq_tlo", data_count);
    DView device_art("fmc_sq_art", data_count);
    DView device_ns("fmc_sq_ns", data_count);
    DView device_ext("fmc_sq_ext", data_count);
    DView device_priors("fmc_sq_priors", kFmcPriorTable);
    DView device_q("fmc_sq_q", q_count > 0 ? q_count : 1);
    DView device_qresp("fmc_sq_qresp", q_count > 0 ? q_count : 1);
    DView device_seq("fmc_sq_seq", data_count);
    DView device_som("fmc_sq_som", data_count);
    DView device_mix("fmc_sq_mix", data_count);
    DView device_bg("fmc_sq_bg", data_count);
    DView device_bgp("fmc_sq_bgp", data_count);
    auto host_logw = Kokkos::create_mirror_view(device_logw);
    auto host_alpha = Kokkos::create_mirror_view(device_alpha);
    auto host_beta = Kokkos::create_mirror_view(device_beta);
    auto host_kind = Kokkos::create_mirror_view(device_kind);
    auto host_alt = Kokkos::create_mirror_view(device_alt);
    auto host_total = Kokkos::create_mirror_view(device_total);
    auto host_indel = Kokkos::create_mirror_view(device_indel);
    auto host_tlo = Kokkos::create_mirror_view(device_tlo);
    auto host_art = Kokkos::create_mirror_view(device_art);
    auto host_ns = Kokkos::create_mirror_view(device_ns);
    auto host_ext = Kokkos::create_mirror_view(device_ext);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    auto host_q = Kokkos::create_mirror_view(device_q);
    auto host_qresp = Kokkos::create_mirror_view(device_qresp);
    auto host_seq = Kokkos::create_mirror_view(device_seq);
    auto host_som = Kokkos::create_mirror_view(device_som);
    auto host_mix = Kokkos::create_mirror_view(device_mix);
    auto host_bg = Kokkos::create_mirror_view(device_bg);
    auto host_bgp = Kokkos::create_mirror_view(device_bgp);
    for (std::size_t i = 0; i < cluster_count; ++i)
        host_logw(i) = params.log_cluster_weights[i];
    for (std::size_t i = 0; i < kFmcMaxClusters; ++i) {
        host_alpha(i) = params.cluster_alpha[i];
        host_beta(i) = params.cluster_beta[i];
        host_kind(i) = params.cluster_kind[i];
    }
    const bool has_indel = !evidence.indel_length.empty();
    for (std::size_t i = 0; i < data_count; ++i) {
        host_alt(i) = evidence.alt_count[i];
        host_total(i) = evidence.total_count[i];
        host_tlo(i) = evidence.tumor_log_odds[i];
        host_art(i) = evidence.artifact_prob[i];
        host_ns(i) = evidence.non_somatic_prob[i];
        host_indel(i) = has_indel ? evidence.indel_length[i] : 0;
        host_ext(i) = use_external ? somatic_external[i] : 0.0;
    }
    for (std::size_t i = 0; i < kFmcPriorTable; ++i)
        host_priors(i) = log_variant_priors_in[i];
    for (std::size_t q = 0; q < q_count; ++q) host_q(q) = quantiles_in[q];
    Kokkos::deep_copy(device_logw, host_logw);
    Kokkos::deep_copy(device_alpha, host_alpha);
    Kokkos::deep_copy(device_beta, host_beta);
    Kokkos::deep_copy(device_kind, host_kind);
    Kokkos::deep_copy(device_alt, host_alt);
    Kokkos::deep_copy(device_total, host_total);
    Kokkos::deep_copy(device_indel, host_indel);
    Kokkos::deep_copy(device_tlo, host_tlo);
    Kokkos::deep_copy(device_art, host_art);
    Kokkos::deep_copy(device_ns, host_ns);
    Kokkos::deep_copy(device_ext, host_ext);
    Kokkos::deep_copy(device_priors, host_priors);
    Kokkos::deep_copy(device_q, host_q);
    Kokkos::deep_copy(device_qresp, 0.0);
    const auto prepare_end = std::chrono::steady_clock::now();
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();

    const auto execute_begin = std::chrono::steady_clock::now();
    const std::size_t qn = q_count;
    const bool external = use_external;
    Kokkos::parallel_for(
        "fastgatk_fmc_split_quantities",
        Kokkos::RangePolicy<ExecSpace>(0, data_count),
        KOKKOS_LAMBDA(const std::size_t datum) {
            const int ccount = static_cast<int>(cluster_count);
            const int alt = device_alt(datum);
            const int total = device_total(datum);
            if (total <= 0 || alt < 0 || alt > total) return;
            // Sequencing posterior under the current model via the shared
            // device helper (Java probabilityOfSequencingError with the
            // per-length variant prior).
            const double seq = fmc_seq_posterior_device(
                device_logw.data(), device_alpha.data(), device_beta.data(),
                ccount, device_tlo(datum), alt, total,
                device_indel(datum), device_priors.data());
            device_seq(datum) = seq;
            // probabilityOfSomaticVariant = (1-art)(1-ns)(1-seq).
            const double somatic =
                (1.0 - device_art(datum)) * (1.0 - device_ns(datum)) *
                (1.0 - seq);
            device_som(datum) = somatic;
            // Count-based mixture (Java logLikelihoodGivenSomatic) and the
            // background-cluster posterior softmax[0] (cluster 0 is the Java
            // background cluster), both from the shared device helper.
            const double high = fmc_count_mixture_device(
                device_logw.data(), device_alpha.data(), device_beta.data(),
                ccount, total, alt);
            device_mix(datum) = high;
            const double bg_given = Kokkos::exp(
                device_logw(0) +
                log_beta_binom_pmf(total, alt, device_alpha(0),
                                   device_beta(0)) -
                high);
            device_bg(datum) = bg_given;
            const double bg_prob =
                (external ? device_ext(datum) : somatic) * bg_given;
            device_bgp(datum) = bg_prob;
            if (qn > 0) {
                for (std::size_t q = 0; q < qn; ++q) {
                    Kokkos::atomic_add(
                        &device_qresp(q),
                        Kokkos::exp(log_binom_pmf(total, alt, device_q(q))) *
                            bg_prob * static_cast<double>(total + 1));
                }
            }
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    Kokkos::deep_copy(host_seq, device_seq);
    Kokkos::deep_copy(host_som, device_som);
    Kokkos::deep_copy(host_mix, device_mix);
    Kokkos::deep_copy(host_bg, device_bg);
    Kokkos::deep_copy(host_bgp, device_bgp);
    Kokkos::deep_copy(host_qresp, device_qresp);
    for (std::size_t i = 0; i < data_count; ++i) {
        result.sequencing_error[i] = host_seq(i);
        result.somatic_probability[i] = host_som(i);
        result.log_somatic_mixture[i] = host_mix(i);
        result.background_given_somatic[i] = host_bg(i);
        result.background_probability[i] = host_bgp(i);
    }
    for (std::size_t q = 0; q < q_count; ++q)
        result.quantile_responsibility[q] = host_qresp(q);
    result.valid = true;
    return result;
}

}  // namespace fastgatk::kernels
