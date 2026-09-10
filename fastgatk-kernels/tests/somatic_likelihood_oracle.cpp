#include "fastgatk/kernels/somatic.hpp"

#include <Kokkos_Core.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

// This executable deliberately has no file or HTSlib dependency.  It emits
// fixed, allele-major log10 likelihood matrices that are small enough to
// audit by hand and are also used by verify_somatic_likelihood_oracle.py as an
// independent implementation of GATK's SomaticLikelihoodsEngine.
int main() {
    Kokkos::initialize();
    try {
        const auto biallelic =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01}, 1, 4, 101);
        const auto multiallelic =
            fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0,
                 -0.01, -0.01, -3.0, -3.0,
                 -3.0, -3.0, -0.01, -0.01}, 3, 4);
        const auto one_sided =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {0.0}, {-std::numeric_limits<double>::infinity()}, 1, 1, 101);
        const auto min_af =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01}, 1, 4, 101, 0.1);
        // GATK keeps tumor and normal AlleleLikelihoods evidence counts
        // independent.  This deliberately uses four tumor fragments and
        // two normal fragments: the normal evidence must be consumed rather
        // than silently discarded because the widths differ.
        const auto unequal_normal =
            fastgatk::kernels::calculate_somatic_posterior_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01},
                {-1.0, -1.0}, {-0.01, -0.01},
                {0}, {0}, 1, 4);
        const auto no_normal =
            fastgatk::kernels::calculate_somatic_posterior_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01},
                {}, {}, {0}, {0}, 1, 4);
        std::cout << std::setprecision(17)
                  << "{\"status\":\"pass\",\"execution_space\":\""
                  << biallelic.execution_space
                  << "\",\"biallelic_tlod\":" << biallelic.tlod.at(0)
                  << ",\"biallelic_nlod\":" << biallelic.normal_log10_odds.at(0)
                  << ",\"biallelic_af\":" << biallelic.best_allele_fraction.at(0)
                  << ",\"multiallelic_tlod\":["
                  << multiallelic.tlod.at(0) << "," << multiallelic.tlod.at(1)
                  << "],\"multiallelic_nlod\":["
                  << multiallelic.normal_log10_odds.at(0) << ","
                  << multiallelic.normal_log10_odds.at(1)
                  << "],\"one_sided_tlod\":" << one_sided.tlod.at(0)
                  << ",\"one_sided_nlod\":" << one_sided.normal_log10_odds.at(0)
                  << ",\"one_sided_af\":" << one_sided.best_allele_fraction.at(0)
                  << ",\"min_af_tlod\":" << min_af.tlod.at(0)
                  << ",\"min_af_af\":" << min_af.best_allele_fraction.at(0)
                  << ",\"unequal_normal_somatic_evidence\":"
                  << unequal_normal.somatic_log10_evidence.at(0)
                  << ",\"no_normal_somatic_evidence\":"
                  << no_normal.somatic_log10_evidence.at(0)
                  << ",\"unequal_normal_informative\":"
                  << unequal_normal.informative_reads.at(0)
                  << "}\n";
        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
