// Full-trajectory oracle for the Java-parity FilterMutectCalls learning path
// (FASTGATK_FMC_JAVA_PASSES semantics) without needing a GATK jar or VCF:
// replays the three accumulate passes over the contamination-joint fixture's
// ten-alt-depth callset directly through the shared kernels module and the
// fmc_java_model.hpp host orchestration, then pins the learned model state
// and threshold.  Guards the R17 fix (Java runs NUM_ITERATIONS = 5 EM(false)
// iterations per accepted split): a regression to 1 EM per split moves the
// background beta from ~1.0455 toward ~1.02 and the threshold by ~5e-5.
#include <Kokkos_Core.hpp>
#include "fmc_java_model.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    Kokkos::initialize();
    namespace fmc_j = fastgatk_fmc_java;
    using fastgatk::kernels::FmcClusterParams;
    using fastgatk::kernels::FmcEvidence;
    using fastgatk::kernels::fmc_error_probabilities_kokkos;

    fmc_j::SomaticModel model;
    model.clusters.push_back(fmc_j::Cluster{0, 1.0, 1.0, 0.5});
    model.clusters.push_back(fmc_j::Cluster{0, 10.0, 1.0, 10.0 / 11.0});
    model.log_cluster_weights = {std::log1p(0.01), std::log(0.01)};
    for (int i = 0; i < 21; ++i)
        model.log_variant_priors[static_cast<std::size_t>(i)] =
            -7.0 * fmc_j::kLn10;
    model.log_variant_priors[10] = -6.0 * fmc_j::kLn10;
    model.callable_sites = 1.0e6;
    model.has_callable_sites = true;

    const int alt_depths[] = {1, 2, 3, 5, 10, 15, 20, 25, 30, 40};
    const double kTlo = 8.0 * fmc_j::kLn10;   // TLOD 8
    const double kPopAf = 0.01;               // POPAF=2
    const double kContamination = 0.1;        // --contamination-estimate 0.1

    auto accumulate_pass = [&](const bool learn_after,
                               std::vector<double>* errors_out) {
        model.data.clear();
        FmcEvidence evidence;
        std::vector<double> tlo;
        for (const int alt : alt_depths) {
            evidence.alt_count.push_back(alt);
            evidence.total_count.push_back(100);
            tlo.push_back(kTlo);
            evidence.artifact_prob.push_back(0.0);
            evidence.non_somatic_prob.push_back(0.0);
            evidence.population_af.push_back(kPopAf);
            evidence.contamination.push_back(kContamination);
        }
        evidence.tumor_log_odds = tlo;
        std::vector<double> priors(model.log_variant_priors.begin(),
                                   model.log_variant_priors.end());
        const auto ep = fmc_error_probabilities_kokkos(
            model.kernel_params(), evidence, priors);
        if (!ep.valid) return 2;
        for (std::size_t i = 0; i < evidence.alt_count.size(); ++i) {
            const double seq = ep.sequencing_error[i];
            const double ns = ep.non_somatic_error[i];
            if (errors_out != nullptr)
                errors_out->push_back(1.0 - (1.0 - seq) * (1.0 - ns));
            model.record(tlo[i], 0.0, ns, evidence.alt_count[i],
                         evidence.total_count[i], 0);
        }
        if (learn_after) model.learn_and_clear();
        return 0;
    };

    if (accumulate_pass(true, nullptr) != 0) return 2;
    if (accumulate_pass(true, nullptr) != 0) return 2;
    std::vector<double> errors;
    if (accumulate_pass(false, &errors) != 0) return 2;
    const double threshold = fmc_j::optimal_f_score_threshold(errors);

    const auto near = [](double a, double b, double eps) {
        return a < b - eps || a > b + eps;
    };
    const auto print = [&]() {
        std::printf("trajectory k=%zu w0=%.17g bg=%.17g/%.17g p10=%.17g "
                    "vv=%.17g threshold=%.17g obs=%zu\n",
                    model.clusters.size(), model.log_cluster_weights[0],
                    model.clusters[0].alpha, model.clusters[0].beta,
                    model.log_variant_priors[10],
                    model.log_variant_vs_artifact_prior, threshold,
                    errors.size());
    };
    print();

    // Machine-validated state (a5_p2 Python reference / GATK EngineProbe,
    // reproduced bit-closely by the kernels trajectory; tolerances cover the
    // gamma_math-vs-libm last-ulp differences).
    bool ok = true;
    if (model.clusters.size() != 7) ok = false;
    if (near(model.clusters[0].alpha, 1.0, 1e-9)) ok = false;
    if (near(model.clusters[0].beta, 1.0454669114043891, 1e-6)) ok = false;
    if (near(model.log_cluster_weights[0], -2.0305792632544235, 1e-6))
        ok = false;
    if (near(model.log_variant_priors[10], -14.136134772191198, 1e-6))
        ok = false;
    if (near(model.log_variant_vs_artifact_prior, -0.4570933809721523, 1e-6))
        ok = false;
    if (near(threshold, 0.19259174653718336, 1e-9)) ok = false;
    if (errors.size() != 10) ok = false;
    Kokkos::finalize();
    if (!ok) return 3;
    std::printf("fmc-trajectory PASS\n");
    return 0;
}
