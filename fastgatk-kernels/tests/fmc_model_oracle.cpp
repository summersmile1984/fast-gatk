// Kernel-level oracle for one Java SomaticClusteringModel EM iteration.
// Prints the updated cluster weights/shapes and prior table for the pinned
// single-datum scenario; the paired Python verifier compares against the
// machine-validated transcription (a5_p1_full.perform_em_iteration).
#include <Kokkos_Core.hpp>
#include "fastgatk/kernels/filter_model_kokkos.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    Kokkos::initialize();
    using fastgatk::kernels::FmcClusterParams;
    using fastgatk::kernels::FmcEvidence;
    using fastgatk::kernels::fmc_model_em_iteration_kokkos;
    using fastgatk::kernels::kFmcMaxClusters;
    using fastgatk::kernels::kFmcPriorTable;
    const double kLn10 = 2.3025850929940456840;
    FmcClusterParams params;
    params.log_cluster_weights = {std::log1p(0.01), std::log(0.01)};
    params.cluster_alpha.assign(kFmcMaxClusters, 1.0);
    params.cluster_beta.assign(kFmcMaxClusters, 1.0);
    params.cluster_kind.assign(kFmcMaxClusters, 0);
    params.cluster_alpha[0] = 1.0; params.cluster_beta[0] = 1.0;
    params.cluster_alpha[1] = 10.0; params.cluster_beta[1] = 1.0;
    FmcEvidence evidence;
    evidence.alt_count = {25};
    evidence.total_count = {100};
    evidence.tumor_log_odds = {8.0 * kLn10};
    evidence.artifact_prob = {0.0};
    evidence.non_somatic_prob = {0.5369};
    std::vector<double> priors(kFmcPriorTable, -7.0 * kLn10);
    priors[10] = -6.0 * kLn10;
    auto result = fmc_model_em_iteration_kokkos(
        params, evidence, priors, -1.0 * kLn10, 1e6, false);
    if (!result.valid) return 2;
    std::printf("weights %.17g %.17g\n", result.log_cluster_weights[0],
                result.log_cluster_weights[1]);
    std::printf("alpha %.17g %.17g\n", result.cluster_alpha[0],
                result.cluster_alpha[1]);
    std::printf("beta %.17g %.17g\n", result.cluster_beta[0],
                result.cluster_beta[1]);
    std::printf("somatic %.17g art %.17g prior0 %.17g\n",
                result.somatic_variant_count, result.technical_artifact_count,
                result.log_variant_priors[10]);

    // Pinned single-iteration expectations from the machine-validated Python
    // transcription (a5_p1_full.perform_em_iteration).
    const double ew0 = -0.52459762908212737;
    const double ew1 = -0.89598195067409214;
    const double ea0 = 1.0;
    const double ea1 = 9.9999999464796474;
    const double eb0 = 1.0310462572784453;
    const double eb1 = 1.0000001281686792;
    const double esomatic = 0.44974136498786943;
    bool ok = true;
    auto near = [&](double a, double b, double eps) {
        if (a < b - eps || a > b + eps) ok = false;
    };
    near(result.log_cluster_weights[0], ew0, 1e-12);
    near(result.log_cluster_weights[1], ew1, 1e-12);
    near(result.cluster_alpha[0], ea0, 1e-12);
    near(result.cluster_alpha[1], ea1, 1e-10);
    near(result.cluster_beta[0], eb0, 1e-10);
    near(result.cluster_beta[1], eb1, 1e-10);
    near(result.somatic_variant_count, esomatic, 1e-12);
    if (!ok) return 3;
    std::printf("scenario-1 PASS\n");

    // Scenario 2: three datums, prior empirical update enabled.
    FmcEvidence evidence2;
    evidence2.alt_count = {25, 30, 40};
    evidence2.total_count = {100, 100, 100};
    evidence2.tumor_log_odds = {8.0 * kLn10, 8.0 * kLn10, 8.0 * kLn10};
    evidence2.artifact_prob = {0.0, 0.0, 0.0};
    evidence2.non_somatic_prob = {0.5369, 0.0005, 0.0};
    auto result2 = fmc_model_em_iteration_kokkos(
        params, evidence2, priors, -1.0 * kLn10, 1e6, true);
    if (!result2.valid) return 2;
    std::printf("s2 weights %.17g %.17g\n", result2.log_cluster_weights[0], result2.log_cluster_weights[1]);
    std::printf("s2 alpha %.17g %.17g\n", result2.cluster_alpha[0], result2.cluster_alpha[1]);
    std::printf("s2 beta %.17g %.17g\n", result2.cluster_beta[0], result2.cluster_beta[1]);
    std::printf("s2 snv %.17g vv %.17g somatic %.17g art %.17g\n", result2.log_variant_priors[10], result2.log_variant_vs_artifact_prior, result2.somatic_variant_count, result2.technical_artifact_count);
    near(result2.log_cluster_weights[0], -0.25840818374048119, 1e-12);
    near(result2.log_cluster_weights[1], -1.4796381923111619, 1e-12);
    near(result2.cluster_alpha[0], 1.0015174529569879, 5e-9);
    near(result2.cluster_alpha[1], 9.999996548337343, 5e-9);
    near(result2.cluster_beta[0], 1.1262779927115272, 5e-9);
    near(result2.cluster_beta[1], 1.0000111134799075, 5e-9);
    near(result2.log_variant_priors[10], -12.943562656179168, 1e-11);
    near(result2.log_variant_vs_artifact_prior, -0.25839419869249314, 1e-11);
    if (!ok) return 4;
    std::printf("scenario-2 PASS\n");

    // Scenario 3: batch error probabilities (sequencing / germline /
    // contamination posteriors) under the pristine two-cluster model for the
    // three alt depths of the contamination-joint fixture.  Pinned values are
    // from the machine-validated Python transcription (a5_p2_reference.py
    // germline_posterior / contamination_posterior / sequencing_posterior at
    // pop AF 0.01, contamination 0.1, TLOD 8, SNV prior -6*ln10).
    FmcEvidence evidence3;
    evidence3.alt_count = {25, 30, 40};
    evidence3.total_count = {100, 100, 100};
    evidence3.tumor_log_odds = {8.0 * kLn10, 8.0 * kLn10, 8.0 * kLn10};
    evidence3.artifact_prob = {0.0, 0.0, 0.0};
    evidence3.non_somatic_prob = {0.0, 0.0, 0.0};
    evidence3.population_af = {0.01, 0.01, 0.01};
    evidence3.contamination = {0.1, 0.1, 0.1};
    auto result3 = fmc_error_probabilities_kokkos(
        params, evidence3, priors);
    if (!result3.valid) return 2;
    std::printf("s3 seq  %.17g %.17g %.17g\n",
                result3.sequencing_error[0], result3.sequencing_error[1],
                result3.sequencing_error[2]);
    std::printf("s3 germ %.17g %.17g %.17g\n",
                result3.germline_error[0], result3.germline_error[1],
                result3.germline_error[2]);
    std::printf("s3 cont %.17g %.17g %.17g\n",
                result3.contamination_error[0],
                result3.contamination_error[1],
                result3.contamination_error[2]);
    std::printf("s3 ns   %.17g %.17g %.17g\n",
                result3.non_somatic_error[0], result3.non_somatic_error[1],
                result3.non_somatic_error[2]);
    const double seq_pin[3] = {0.028846113176701632,
                               0.028846017922356672,
                               0.028844917342607534};
    const double germ_pin[3] = {0.5369248971183561,
                                0.992929251681918,
                                0.9999847835853424};
    const double cont_pin[3] = {0.2121513079953642,
                                5.518196476275156e-4,
                                7.41031024178629e-11};
    for (int i = 0; i < 3; ++i) {
        near(result3.sequencing_error[static_cast<std::size_t>(i)],
             seq_pin[i], 1e-10);
        near(result3.germline_error[static_cast<std::size_t>(i)],
             germ_pin[i], 1e-10);
        near(result3.contamination_error[static_cast<std::size_t>(i)],
             cont_pin[i], 1e-10);
        near(result3.non_somatic_error[static_cast<std::size_t>(i)],
             germ_pin[i] > cont_pin[i] ? germ_pin[i] : cont_pin[i], 1e-10);
    }
    if (!ok) return 5;
    std::printf("scenario-3 PASS\n");

    // Scenario 4: Java SomaticClusteringModel initialization split-quantity
    // kernel under the pristine model for the single recorded pass-0 datum
    // (alt25/100, TLOD 8, ns = max(germ,contam) = 0.5369).  Pinned from
    // a5_p1_full probabilityOfSomaticVariant / logLikelihoodGivenSomatic /
    // background softmax and the Java quantile background responsibility at
    // the pass-0 quantile grid value 0.25.
    FmcEvidence evidence4;
    evidence4.alt_count = {25};
    evidence4.total_count = {100};
    evidence4.tumor_log_odds = {8.0 * kLn10};
    evidence4.artifact_prob = {0.0};
    evidence4.non_somatic_prob = {0.5369};
    auto result4 = fmc_model_split_quantities_kokkos(
        params, evidence4, priors, {}, {0.25});
    if (!result4.valid) return 2;
    std::printf("s4 seq  %.17g\n", result4.sequencing_error[0]);
    std::printf("s4 som  %.17g\n", result4.somatic_probability[0]);
    std::printf("s4 mix  %.17g\n", result4.log_somatic_mixture[0]);
    std::printf("s4 bg   %.17g\n", result4.background_given_somatic[0]);
    std::printf("s4 bgp  %.17g\n", result4.background_probability[0]);
    std::printf("s4 qr   %.17g\n", result4.quantile_responsibility[0]);
    near(result4.sequencing_error[0], 0.028846113176701632, 1e-10);
    near(result4.somatic_probability[0], 0.44974136498786943, 1e-10);
    near(result4.log_somatic_mixture[0], -4.605169067568704, 1e-10);
    near(result4.background_given_somatic[0], 0.9999988815812243, 1e-10);
    near(result4.background_probability[0], 0.4497408619886826, 1e-10);
    near(result4.quantile_responsibility[0], 4.169893323056797, 1e-10);
    if (!ok) return 6;
    std::printf("scenario-4 PASS\n");

    // Scenario 5: contamination estimate 0 (degenerate no-contamination).
    // Regression for the NaN fix: with zero contaminant mass at both
    // hypotheses the contamination posterior must be exactly 0, and
    // sequencing/germline values must be unchanged vs scenario 3, so the
    // non-somatic group error stays the germline posterior.
    FmcEvidence evidence5 = evidence3;
    evidence5.contamination = {0.0, 0.0, 0.0};
    auto result5 = fmc_error_probabilities_kokkos(params, evidence5, priors);
    if (!result5.valid) return 2;
    std::printf("s5 seq  %.17g %.17g %.17g\n",
                result5.sequencing_error[0], result5.sequencing_error[1],
                result5.sequencing_error[2]);
    std::printf("s5 germ %.17g %.17g %.17g\n",
                result5.germline_error[0], result5.germline_error[1],
                result5.germline_error[2]);
    std::printf("s5 cont %.17g %.17g %.17g\n",
                result5.contamination_error[0], result5.contamination_error[1],
                result5.contamination_error[2]);
    for (int i = 0; i < 3; ++i) {
        near(result5.sequencing_error[static_cast<std::size_t>(i)],
             seq_pin[i], 1e-10);
        near(result5.germline_error[static_cast<std::size_t>(i)],
             germ_pin[i], 1e-10);
        near(result5.contamination_error[static_cast<std::size_t>(i)],
             0.0, 0.0);
        near(result5.non_somatic_error[static_cast<std::size_t>(i)],
             germ_pin[i], 1e-10);
    }
    if (!ok) return 7;
    std::printf("scenario-5 PASS\n");
    std::printf("fmc-model-oracle PASS\n");
    Kokkos::finalize();
    return 0;
}
