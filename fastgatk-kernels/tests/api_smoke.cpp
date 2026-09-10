#include "fastgatk/kernels/pairhmm.hpp"
#include "fastgatk/kernels/pairhmm_kokkos.hpp"
#include "fastgatk/kernels/smith_waterman.hpp"
#include "fastgatk/kernels/kmer_graph.hpp"
#include "fastgatk/kernels/activity_profile.hpp"
#include "fastgatk/kernels/read_filter.hpp"
#include "fastgatk/kernels/bqsr.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "fastgatk/kernels/reference_confidence.hpp"
#include "fastgatk/kernels/somatic.hpp"
#include "fastgatk/io/flow_codec.hpp"
#include "fastgatk/io/hts_reader.hpp"

#include <Kokkos_Core.hpp>

#include "fastgatk/kernels/read_error_correction.hpp"
#include "fastgatk/kernels/pileup_error_correction.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

int main() {
    Kokkos::initialize();
    try {
        fastgatk::pairhmm::PairInput input;
        input.haplotype = {'A', 'C', 'G', 'T'};
        input.read = {'A', 'C', 'G', 'T'};
        input.read_qual = {30, 30, 30, 30};
        input.insertion_gop = {40, 40, 40, 40};
        input.deletion_gop = {40, 40, 40, 40};
        input.gap_continuation = {10, 10, 10, 10};
        fastgatk::pairhmm::PairIndexBatch pairs{{0}, {0}};
        const auto result = fastgatk::pairhmm::compute_kokkos({input}, pairs, 1);
        if (result.likelihoods.size() != 1 || !std::isfinite(result.likelihoods[0]))
            throw std::runtime_error("PairHMM reusable API smoke failed");
        // PairHMMModel's Java table path clamps the sum of insertion/deletion
        // event probabilities at one and converts it through log1p.  Quality
        // zero is therefore a valid boundary input with a non-negative
        // no-event transition; the old direct `1 - 10^x` construction made
        // this transition negative and could produce NaNs downstream.
        auto boundary_input = input;
        boundary_input.insertion_gop.assign(boundary_input.read.size(), 0);
        boundary_input.deletion_gop.assign(boundary_input.read.size(), 0);
        const auto boundary_result = fastgatk::pairhmm::compute_kokkos(
            {boundary_input}, pairs, 1);
        if (boundary_result.scaled_sums.size() != 1 ||
            !std::isfinite(boundary_result.scaled_sums[0]) ||
            boundary_result.scaled_sums[0] < 0.0 ||
            !std::isfinite(boundary_result.likelihoods[0]))
            throw std::runtime_error("PairHMM quality-zero transition boundary failed");
        const auto compat_scalar = fastgatk::pairhmm::compute_scalar({input}, 1, 1);
        const auto compat_simd = fastgatk::pairhmm::compute_simd(
            {input}, 1, 1, fastgatk::pairhmm::SimdBackend::Auto);
        const auto compiled_backend = fastgatk::pairhmm::compiled_simd_backend();
        if (fastgatk::pairhmm::compiled_simd_width() == 0 ||
            fastgatk::pairhmm::compiled_simd_backend_name().empty() ||
            (compiled_backend == fastgatk::pairhmm::SimdBackend::Avx512 &&
             !fastgatk::pairhmm::cpu_supports_avx512()) ||
            (compiled_backend == fastgatk::pairhmm::SimdBackend::Avx2 &&
             !fastgatk::pairhmm::cpu_supports_avx2()))
            throw std::runtime_error("Kokkos SIMD ABI/runtime capability contract failed");
        const auto compat_selected = fastgatk::pairhmm::compute_simd(
            {input}, 1, 1, compiled_backend);
        if (compat_selected.likelihoods.size() != 1 ||
            compat_selected.likelihoods[0] != result.likelihoods[0])
            throw std::runtime_error("compiled Kokkos SIMD ABI selection failed");
        const auto incompatible_backend = compiled_backend == fastgatk::pairhmm::SimdBackend::Avx512
            ? fastgatk::pairhmm::SimdBackend::Avx2
            : fastgatk::pairhmm::SimdBackend::Avx512;
        bool rejected_incompatible = false;
        try {
            (void)fastgatk::pairhmm::compute_simd({input}, 1, 1, incompatible_backend);
        } catch (const std::runtime_error&) {
            rejected_incompatible = true;
        }
        if (!rejected_incompatible)
            throw std::runtime_error("incompatible Kokkos SIMD ABI request was not rejected");
        if (compat_scalar.likelihoods.size() != 1 || compat_simd.likelihoods.size() != 1 ||
            compat_scalar.likelihoods[0] != result.likelihoods[0] ||
            compat_simd.likelihoods[0] != result.likelihoods[0])
            throw std::runtime_error("Kokkos PairHMM compatibility forwarding API failed");
        const auto normalized = fastgatk::pairhmm::normalize_likelihoods_kokkos(
            {0.0, -10.0, -20.0, -3.0}, {0, 0, 0, 1}, {1, 1, 0, 0}, 2, 4.5);
        if (normalized.likelihoods.size() != 4 || normalized.best_by_read.size() != 2 ||
            normalized.best_by_read[0] != 0.0 || !std::isinf(normalized.best_by_read[1]) ||
            normalized.best_by_read[1] >= 0.0 ||
            normalized.likelihoods[1] != -4.5 || normalized.likelihoods[2] != -4.5 ||
            normalized.likelihoods[3] != -3.0 || normalized.execution_space.empty() ||
            normalized.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos GATK likelihood-normalization API smoke failed");
        const auto uncapped = fastgatk::pairhmm::normalize_likelihoods_kokkos(
            {0.0, -std::numeric_limits<double>::infinity(), -3.0},
            {0, 0, 0}, {1, 1, 1}, 1,
            -std::numeric_limits<double>::infinity());
        if (uncapped.likelihoods.size() != 3 || uncapped.best_by_read.at(0) != 0.0 ||
            !std::isinf(uncapped.likelihoods.at(1)) || uncapped.likelihoods.at(1) >= 0.0 ||
            uncapped.likelihoods.at(2) != -3.0)
            throw std::runtime_error("Kokkos GATK disabled likelihood cap changed -Infinity");
        const auto marginalized =
            fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos(
                {{0, 0, -2.0}, {0, 0, -1.0}, {0, 1, -4.0}, {0, 1, -0.5},
                 {1, 0, -3.0}},
                2);
        if (marginalized.best_by_row_allele.size() != 4 ||
            marginalized.allele_count != 2 ||
            marginalized.best_by_row_allele[0] != -1.0 ||
            marginalized.best_by_row_allele[1] != -0.5 ||
            marginalized.best_by_row_allele[2] != -3.0 ||
            !std::isinf(marginalized.best_by_row_allele[3]) ||
            marginalized.execution_space.empty() ||
            marginalized.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos allele max-marginalization API smoke failed");
        const auto multiallelic_marginalized =
            fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos(
                {{0, 0, -4.0}, {0, 1, -2.0}, {0, 2, -1.0},
                 {0, 2, -3.0}, {1, 0, -5.0}, {1, 2, -0.25}},
                2, 3);
        if (multiallelic_marginalized.allele_count != 3 ||
            multiallelic_marginalized.best_by_row_allele.size() != 6 ||
            multiallelic_marginalized.best_by_row_allele[0] != -4.0 ||
            multiallelic_marginalized.best_by_row_allele[1] != -2.0 ||
            multiallelic_marginalized.best_by_row_allele[2] != -1.0 ||
            multiallelic_marginalized.best_by_row_allele[3] != -5.0 ||
            !std::isinf(multiallelic_marginalized.best_by_row_allele[4]) ||
            multiallelic_marginalized.best_by_row_allele[5] != -0.25 ||
            multiallelic_marginalized.execution_space.empty() ||
            multiallelic_marginalized.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos multi-allelic max-marginalization API smoke failed");
        const auto fragment_aggregation =
            fastgatk::pairhmm::aggregate_fragment_haplotype_likelihoods_kokkos(
                {{0, -1.0},
                 {1, -2.0},
                 {0, -std::numeric_limits<double>::infinity()},
                 {1, -3.0}},
                3);
        if (fragment_aggregation.sums_by_cell.size() != 3 ||
            !std::isinf(fragment_aggregation.sums_by_cell[0]) ||
            fragment_aggregation.sums_by_cell[0] >= 0.0 ||
            fragment_aggregation.sums_by_cell[1] != -5.0 ||
            fragment_aggregation.sums_by_cell[2] != 0.0 ||
            fragment_aggregation.execution_space.empty() ||
            fragment_aggregation.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos fragment-haplotype aggregation API smoke failed");
        const auto uncertainty =
            fastgatk::pairhmm::reduce_read_allele_uncertainty_kokkos(
                {{0, -2.0}, {0, -1.0}, {0, -1.0}, {1, -3.0}}, 2);
        if (uncertainty.best_second_by_row.size() != 4 ||
            uncertainty.best_second_by_row[0] != -1.0 ||
            uncertainty.best_second_by_row[1] != -1.0 ||
            uncertainty.best_second_by_row[2] != -3.0 ||
            !std::isinf(uncertainty.best_second_by_row[3]) ||
            uncertainty.execution_space.empty() ||
            uncertainty.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos read-allele uncertainty API smoke failed");
        const auto best_alleles =
            fastgatk::pairhmm::reduce_read_allele_best_kokkos(
                {{0, 0, -0.5}, {0, 1, -0.5}, {0, 2, -2.0},
                 {1, 0, -3.0}, {1, 1, -1.0}}, 2);
        if (best_alleles.best_second_allele_by_row.size() != 4 ||
            best_alleles.best_second_likelihood_by_row.size() != 4 ||
            best_alleles.best_second_allele_by_row[0] != 0U ||
            best_alleles.best_second_allele_by_row[1] != 1U ||
            best_alleles.best_second_likelihood_by_row[0] != -0.5 ||
            best_alleles.best_second_likelihood_by_row[1] != -0.5 ||
            best_alleles.best_second_allele_by_row[2] != 1U ||
            best_alleles.best_second_allele_by_row[3] != 0U ||
            best_alleles.best_second_likelihood_by_row[2] != -1.0 ||
            best_alleles.best_second_likelihood_by_row[3] != -3.0 ||
            best_alleles.execution_space.empty() ||
            best_alleles.execution_policy != "RangePolicy")
            throw std::runtime_error("Kokkos BestAllele reduction API smoke failed");
        const auto somatic = fastgatk::kernels::calculate_somatic_likelihood_kokkos(
            {-1.0, -1.0, -1.0, -1.0},
            {-0.01, -0.01, -0.01, -0.01}, 1, 4, 101);
        if (somatic.tlod.size() != 1 || somatic.best_allele_fraction.size() != 1 ||
            somatic.tlod[0] <= 0.0 || somatic.best_allele_fraction[0] < 0.75 ||
            somatic.informative_reads[0] != 4 || somatic.execution_space.empty())
            throw std::runtime_error("Somatic mixture likelihood API smoke failed");
        // Pinned GATK 4.6.2.0 SomaticLikelihoodsEngine.logEvidence oracle for
        // this matrix: logEvidence(all)-logEvidence(REF-only) in log10 space.
        if (std::abs(somatic.tlod[0] - 3.2838961908386195) > 1.0e-11 ||
            std::abs(somatic.best_allele_fraction[0] - 0.8240236497641983) > 1.0e-11)
            throw std::runtime_error("Somatic biallelic GATK oracle drift");
        // A CSR row with all four retained fragment cells must be bitwise
        // equivalent to the established dense Kokkos path.  This is the
        // production representation for wide Mutect2 windows, where Host
        // membership is sparse but numerical evidence remains device-owned.
        const auto sparse_somatic =
            fastgatk::kernels::calculate_somatic_likelihood_sparse_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01}, {0, 4}, 1, 4, 101);
        if (sparse_somatic.tlod.size() != 1 ||
            sparse_somatic.tlod[0] != somatic.tlod[0] ||
            sparse_somatic.normal_log10_odds[0] != somatic.normal_log10_odds[0] ||
            sparse_somatic.best_allele_fraction[0] != somatic.best_allele_fraction[0] ||
            sparse_somatic.reference_log10_likelihood[0] != somatic.reference_log10_likelihood[0] ||
            sparse_somatic.best_log10_likelihood[0] != somatic.best_log10_likelihood[0] ||
            sparse_somatic.informative_reads[0] != somatic.informative_reads[0])
            throw std::runtime_error("Sparse/dense Kokkos somatic parity failed");
        // Omitted CSR cells are precisely dense both-missing cells. Exercise
        // unequal row widths and one-sided -Infinity values as used by real
        // retained fragment sets, rather than only the rectangular happy path.
        const std::vector<double> irregular_reference{
            -1.0, -1.0e300, -2.0, -3.0,
            -0.5, -1.0e300, -4.0, -1.0e300};
        const std::vector<double> irregular_alternate{
            -0.01, -std::numeric_limits<double>::infinity(), -5.0, -0.2,
            -4.0, -std::numeric_limits<double>::infinity(), -0.1,
            -std::numeric_limits<double>::infinity()};
        const auto irregular_dense = fastgatk::kernels::calculate_somatic_likelihood_kokkos(
            irregular_reference, irregular_alternate, 2, 4, 101);
        const auto irregular_sparse =
            fastgatk::kernels::calculate_somatic_likelihood_sparse_kokkos(
                {-1.0, -2.0, -3.0, -0.5, -4.0},
                {-0.01, -5.0, -0.2, -4.0, -0.1}, {0, 3, 5}, 2, 4, 101);
        for (std::size_t candidate = 0; candidate < 2; ++candidate) {
            if (irregular_sparse.tlod[candidate] != irregular_dense.tlod[candidate] ||
                irregular_sparse.normal_log10_odds[candidate] !=
                    irregular_dense.normal_log10_odds[candidate] ||
                irregular_sparse.best_allele_fraction[candidate] !=
                    irregular_dense.best_allele_fraction[candidate] ||
                irregular_sparse.reference_log10_likelihood[candidate] !=
                    irregular_dense.reference_log10_likelihood[candidate] ||
                irregular_sparse.best_log10_likelihood[candidate] !=
                    irregular_dense.best_log10_likelihood[candidate] ||
                irregular_sparse.informative_reads[candidate] !=
                    irregular_dense.informative_reads[candidate])
                throw std::runtime_error("Irregular sparse/dense Kokkos somatic parity failed");
        }
        // A positive GATK minimum-allele-fraction prior must reach the same
        // Kokkos evidence path (rather than being silently accepted and
        // ignored at the CLI boundary).  minAF=0.1 maps to an ALT
        // pseudocount of 1.301..., so both the posterior AF and evidence
        // should move while remaining finite.
        const auto min_af_somatic =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0},
                {-0.01, -0.01, -0.01, -0.01}, 1, 4, 101, 0.1);
        if (!std::isfinite(min_af_somatic.tlod[0]) ||
            !std::isfinite(min_af_somatic.best_allele_fraction[0]) ||
            std::abs(min_af_somatic.tlod[0] - somatic.tlod[0]) < 1.0e-6 ||
            std::abs(min_af_somatic.best_allele_fraction[0] -
                     somatic.best_allele_fraction[0]) < 1.0e-6)
            throw std::runtime_error("Somatic minimum-allele-fraction prior was ignored");
        // Stress Commons-Math Gamma.digamma's x >= 49 asymptotic branch.  A
        // previous native implementation omitted the 1/120 Bernoulli term,
        // which is small per call but shifts the converged posterior once a
        // locus has roughly one hundred reads.  The matrix is log10 here;
        // GATK's Java oracle uses the equivalent natural-log values.  The
        // expected values are from GATK 4.6.2.0's
        // SomaticLikelihoodsEngine.logEvidence/alleleFractionsPosterior.
        std::vector<double> digamma_stress_reference(100, 0.0);
        std::vector<double> digamma_stress_alternate(100, -0.01);
        const auto digamma_stress =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                digamma_stress_reference, digamma_stress_alternate, 1, 100, 101);
        if (digamma_stress.informative_reads[0] != 100 ||
            std::abs(digamma_stress.tlod[0] - (-1.229861580035300)) > 1.0e-11 ||
            std::abs(digamma_stress.best_allele_fraction[0] -
                     0.1965127942094241) > 1.0e-11)
            throw std::runtime_error("Somatic digamma asymptotic GATK oracle drift");
        const auto biallelic_general =
            fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0,
                 -0.01, -0.01, -0.01, -0.01}, 2, 4);
        if (biallelic_general.tlod.size() != 1 ||
            std::abs(biallelic_general.tlod[0] - somatic.tlod[0]) > 1.0e-10 ||
            std::abs(biallelic_general.best_allele_fraction[0] -
                     somatic.best_allele_fraction[0]) > 1.0e-10)
            throw std::runtime_error("Biallelic/multiallelic somatic parity failed");
        // GATK's NaturalLogUtils treats -Infinity as a zero-probability
        // allele, not as a missing evidence column.  With one read supporting
        // REF only, the [1,1] Dirichlet posterior is [2,1], giving AF=1/3 and
        // log10 evidence=-log10(2).  This also guards against 0 * -Infinity
        // in the likelihood contribution.
        const auto one_sided_missing =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {0.0}, {-std::numeric_limits<double>::infinity()}, 1, 1, 101);
        if (one_sided_missing.informative_reads[0] != 1 ||
            std::abs(one_sided_missing.best_allele_fraction[0] - 1.0 / 3.0) > 1.0e-10 ||
            std::abs(one_sided_missing.tlod[0] + std::log10(2.0)) > 1.0e-10)
            throw std::runtime_error("Somatic one-sided missing likelihood semantics failed");
        const auto all_missing_somatic =
            fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                {-1.0e300}, {-std::numeric_limits<double>::infinity()}, 1, 1, 101);
        if (all_missing_somatic.informative_reads[0] != 0 ||
            all_missing_somatic.tlod[0] != 0.0 ||
            all_missing_somatic.best_allele_fraction[0] != 0.0)
            throw std::runtime_error("Somatic all-missing likelihood semantics failed");
        const auto multiallelic_somatic =
            fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                {-1.0, -1.0, -1.0, -1.0,
                 -0.01, -0.01, -3.0, -3.0,
                 -3.0, -3.0, -0.01, -0.01}, 3, 4);
        if (multiallelic_somatic.tlod.size() != 2 ||
            multiallelic_somatic.best_allele_fraction.size() != 2 ||
            multiallelic_somatic.informative_reads[0] != 4 ||
            multiallelic_somatic.informative_reads[1] != 4 ||
            !std::isfinite(multiallelic_somatic.tlod[0]) ||
            !std::isfinite(multiallelic_somatic.tlod[1]) ||
            multiallelic_somatic.tlod[0] <= 0.0 ||
            multiallelic_somatic.tlod[1] <= 0.0 ||
            multiallelic_somatic.best_allele_fraction[0] <= 0.0 ||
            multiallelic_somatic.best_allele_fraction[1] <= 0.0 ||
            multiallelic_somatic.best_allele_fraction[0] +
                    multiallelic_somatic.best_allele_fraction[1] >= 1.0 ||
            multiallelic_somatic.execution_space.empty())
            throw std::runtime_error("Multiallelic somatic evidence API smoke failed");
        // The same pinned Java oracle with REF/ALT1/ALT2 rows (the API uses
        // log10 likelihoods; Java's RealMatrix path uses natural logs) yields
        // 1.448530966859346 for each ALT after excluding that ALT in turn.
        if (std::abs(multiallelic_somatic.tlod[0] - 1.448530966859346) > 1.0e-11 ||
            std::abs(multiallelic_somatic.tlod[1] - 1.448530966859346) > 1.0e-11)
            throw std::runtime_error("Somatic multiallelic GATK oracle drift");
        const auto somatic_posterior = fastgatk::kernels::calculate_somatic_posterior_kokkos(
            {-5.0, -5.0, -5.0, -5.0},
            {-0.01, -0.01, -0.01, -0.01},
            {-0.01, -0.01, -0.01, -0.01},
            {-5.0, -5.0, -5.0, -5.0},
            {2}, {2}, 1, 4, 0.10);
        const auto posterior_sum = somatic_posterior.somatic_probability[0] +
            somatic_posterior.germline_probability[0] + somatic_posterior.artifact_probability[0];
        if (somatic_posterior.somatic_probability.size() != 1 ||
            somatic_posterior.somatic_probability[0] <= somatic_posterior.germline_probability[0] ||
            somatic_posterior.somatic_probability[0] <= somatic_posterior.artifact_probability[0] ||
            somatic_posterior.orientation_bias_probability[0] > 0.5000001 ||
            std::abs(posterior_sum - 1.0) > 1.0e-9 || somatic_posterior.execution_space.empty())
            throw std::runtime_error("Somatic posterior API smoke failed");
        const auto biased_posterior = fastgatk::kernels::calculate_somatic_posterior_kokkos(
            {-5.0, -5.0, -5.0, -5.0},
            {-0.01, -0.01, -0.01, -0.01},
            {}, {}, {4}, {0}, 1, 4);
        if (biased_posterior.orientation_bias_probability[0] <= 0.5 ||
            biased_posterior.artifact_probability[0] <= somatic_posterior.artifact_probability[0])
            throw std::runtime_error("Somatic orientation posterior API smoke failed");
        fastgatk::pairhmm::PairHmmRead bucket_read;
        bucket_read.bases = {'A', 'C', 'G'};
        bucket_read.qualities = {30, 30, 30};
        bucket_read.insertion_gop = {40, 40, 40};
        bucket_read.deletion_gop = {40, 40, 40};
        bucket_read.gap_continuation = {10, 10, 10};
        fastgatk::pairhmm::PairHmmHaplotype bucket_haplotype;
        bucket_haplotype.bases = {'A', 'C', 'G', 'T', 'A'};
        const auto bucketed = fastgatk::pairhmm::compute_kokkos_bucketed(
            {bucket_read}, {bucket_haplotype}, {{0, 0}}, 1);
        if (bucketed.likelihoods.size() != 1 || !std::isfinite(bucketed.likelihoods[0]))
            throw std::runtime_error("PairHMM variable-length bucket API smoke failed");
        const auto float_bucketed = fastgatk::pairhmm::compute_kokkos_bucketed(
            {bucket_read}, {bucket_haplotype}, {{0, 0}}, 1,
            fastgatk::pairhmm::PairHmmPrecision::Float32);
        if (float_bucketed.precision != "float32" || float_bucketed.likelihoods.size() != 1 ||
            !std::isfinite(float_bucketed.likelihoods[0]) ||
            std::abs(float_bucketed.likelihoods[0] - bucketed.likelihoods[0]) > 1.0)
            throw std::runtime_error("PairHMM Kokkos float32 API smoke failed");
        fastgatk::pairhmm::PairHmmRead second_read = bucket_read;
        second_read.bases.push_back('T');
        second_read.qualities.push_back(30);
        second_read.insertion_gop.push_back(40);
        second_read.deletion_gop.push_back(40);
        second_read.gap_continuation.push_back(10);
        fastgatk::pairhmm::PairHmmHaplotype second_haplotype = bucket_haplotype;
        second_haplotype.bases.push_back('C');
        const auto full_matrix = fastgatk::pairhmm::compute_kokkos_full_matrix(
            {bucket_read, second_read}, {bucket_haplotype, second_haplotype}, 1);
        if (full_matrix.matrix_rows != 2 || full_matrix.matrix_columns != 2 ||
            full_matrix.likelihoods.size() != 4 || full_matrix.scaled_sums.size() != 4 ||
            !std::isfinite(full_matrix.likelihoods[0]) ||
            full_matrix.likelihoods[0] != bucketed.likelihoods[0])
            throw std::runtime_error("PairHMM ragged full-matrix API smoke failed");
        fastgatk::pairhmm::PersistentBucketPlan persistent_plan;
        const auto persistent_a = persistent_plan.execute(
            {bucket_read}, {bucket_haplotype}, {{0, 0}}, 1);
        const auto persistent_b = persistent_plan.execute(
            {bucket_read}, {bucket_haplotype}, {{0, 0}}, 1);
        if (persistent_a.likelihoods != persistent_b.likelihoods ||
            persistent_plan.cached_shapes() != 1 || persistent_plan.cache_hits() < 1 ||
            persistent_b.cache_hits < 1 || !std::isfinite(persistent_b.likelihoods[0]))
            throw std::runtime_error("Persistent PairHMM plan/cache API smoke failed");
        // Flow-space PairHMM is a separate recurrence (4-base flow stride),
        // not a regular-base approximation.  A calibrated probability table
        // makes the test independent of HTSlib flow-tag parsing while still
        // checking model dispatch, length bucketing and stable request order.
        fastgatk::pairhmm::FlowPairHmmRead flow_read;
        flow_read.key = {1, 1, 1, 1};
        flow_read.flow_order = {'A', 'C', 'G', 'T'};
        flow_read.insertion_gop = {40, 40, 40, 40};
        flow_read.deletion_gop = {40, 40, 40, 40};
        flow_read.gap_continuation = {10, 10, 10, 10};
        flow_read.probabilities.assign(4 * 256, 0.001);
        for (std::size_t i = 0; i < 4; ++i)
            flow_read.probabilities[i * 256 + 1] = 0.999;
        fastgatk::pairhmm::FlowPairHmmHaplotype flow_haplotype;
        flow_haplotype.key = {1, 1, 1, 1};
        flow_haplotype.flow_order = {'A', 'C', 'G', 'T'};
        fastgatk::pairhmm::FlowPairHmmHaplotype flow_alt = flow_haplotype;
        flow_alt.key[0] = 2;
        const auto flow_result = fastgatk::pairhmm::compute_kokkos_flow(
            {flow_read}, {flow_haplotype, flow_alt}, {{0, 0}, {0, 1}}, 1);
        if (flow_result.error_model != "flow" || flow_result.execution_policy != "TeamPolicy" ||
            flow_result.likelihoods.size() != 2 ||
            !std::isfinite(flow_result.likelihoods[0]) ||
            !std::isfinite(flow_result.likelihoods[1]) ||
            flow_result.likelihoods[0] <= flow_result.likelihoods[1] ||
            flow_result.execution_space.empty())
            throw std::runtime_error("Flow PairHMM Kokkos API smoke failed");
        // Host decoder contract: this mirrors FlowBasedKeyCodec plus the
        // production tp/t0 matrix rules before the flat table crosses into
        // Kokkos.  Keep boundary spreading disabled here so the asserted
        // cells are directly attributable to the tags.
        fastgatk::io::FlowReadTags flow_tags;
        flow_tags.bases = {'A', 'G'};
        flow_tags.qualities = {30, 30};
        flow_tags.tp = {0, 1};
        flow_tags.t0_phred = {20, 20};
        flow_tags.flow_order = "ACGT";
        flow_tags.max_hmer = 12;
        flow_tags.filling_value = 0.001;
        flow_tags.use_t0_tag = true;
        flow_tags.keep_boundary_flows = true;
        const auto decoded_flow = fastgatk::io::decode_flow_read(flow_tags);
        const auto encoded_flow = fastgatk::io::encode_flow_key(flow_tags.bases, flow_tags.flow_order);
        if (decoded_flow.key != std::vector<std::int32_t>({1, 0, 1}) ||
            decoded_flow.flow_order != std::vector<std::uint8_t>({'A', 'C', 'G'}) ||
            decoded_flow.probabilities.size() != 3 * 256 ||
            decoded_flow.probabilities[1 * 256 + 1] < 0.0099 ||
            decoded_flow.probabilities[0 * 256 + 1] <= 0.0 ||
            decoded_flow.probabilities[0 * 256 + 255] != decoded_flow.probabilities[0 * 256 + 12])
            throw std::runtime_error("FlowBasedRead Host tag decoder smoke failed");
        if (encoded_flow.key != decoded_flow.key ||
            encoded_flow.flow_to_base != std::vector<std::int32_t>({-1, 0, 0}) ||
            encoded_flow.reverse_key != std::vector<std::int32_t>({1, 0, 1}) ||
            fastgatk::io::find_left_flow_clipping(1, encoded_flow) != std::make_pair<std::size_t, std::size_t>(2, 0) ||
            fastgatk::io::find_right_flow_clipping(1, encoded_flow) != std::make_pair<std::size_t, std::size_t>(2, 0))
            throw std::runtime_error("FlowBasedHaplotype key/clipping Host mapping smoke failed");
        // Base-space AssemblyRegion clipping must update the flow key and
        // the two boundary probability columns exactly like
        // FlowBasedRead.applyBaseClipping.  Keep spreading disabled here so
        // the asserted values isolate the hmer shift itself.
        fastgatk::io::DecodedFlowRead clipped_flow;
        clipped_flow.key = {3, 0, 2};
        clipped_flow.flow_order = {'A', 'C', 'G'};
        clipped_flow.max_hmer = 4;
        clipped_flow.filling_value = 0.001;
        clipped_flow.probabilities.resize(3 * 256);
        for (std::size_t flow = 0; flow < 3; ++flow)
            for (std::size_t hmer = 0; hmer < 256; ++hmer)
                clipped_flow.probabilities[flow * 256 + hmer] =
                    static_cast<double>(100 * flow + hmer);
        fastgatk::io::clip_decoded_flow_read(clipped_flow, 1, 1, false);
        if (clipped_flow.key != std::vector<std::int32_t>({2, 0, 1}) ||
            clipped_flow.flow_order != std::vector<std::uint8_t>({'A', 'C', 'G'}) ||
            clipped_flow.probabilities[0 * 256 + 0] != 1.0 ||
            clipped_flow.probabilities[0 * 256 + 3] != 4.0 ||
            clipped_flow.probabilities[2 * 256 + 0] != 201.0 ||
            clipped_flow.probabilities[2 * 256 + 3] != 204.0 ||
            clipped_flow.probabilities[0 * 256 + 255] != 0.0 ||
            clipped_flow.probabilities[2 * 256 + 255] != 0.0)
            throw std::runtime_error("FlowBasedRead base clipping/key-matrix shift smoke failed");
        // GATK invokes applyClipping even when trimToHaplotype computes no
        // base removal.  Its zero-width operation still drops terminal zero
        // flow columns, which changes the phase of the direct flow scorer.
        fastgatk::io::DecodedFlowRead zero_clip_flow;
        zero_clip_flow.key = {0, 1, 0};
        zero_clip_flow.flow_order = {'A', 'C', 'G'};
        zero_clip_flow.max_hmer = 4;
        zero_clip_flow.filling_value = 0.001;
        zero_clip_flow.probabilities.resize(3 * 256);
        for (std::size_t flow = 0; flow < 3; ++flow)
            for (std::size_t hmer = 0; hmer < 256; ++hmer)
                zero_clip_flow.probabilities[flow * 256 + hmer] =
                    static_cast<double>(100 * flow + hmer);
        fastgatk::io::clip_decoded_flow_read(zero_clip_flow, 0, 0, false);
        if (zero_clip_flow.key != std::vector<std::int32_t>({1}) ||
            zero_clip_flow.flow_order != std::vector<std::uint8_t>({'C'}) ||
            zero_clip_flow.probabilities.size() != 256 ||
            zero_clip_flow.probabilities[0] != 100.0 ||
            zero_clip_flow.probabilities[255] != 355.0)
            throw std::runtime_error("FlowBasedRead zero-width clipping diverges from GATK");
        const auto first_hmer_flow = fastgatk::io::encode_flow_key(
            std::vector<std::uint8_t>{'A', 'A', 'A', 'C'}, "ACGT");
        if (fastgatk::io::find_left_flow_clipping(1, first_hmer_flow) !=
                std::make_pair<std::size_t, std::size_t>(0, 1) ||
            fastgatk::io::find_left_flow_clipping(2, first_hmer_flow) !=
                std::make_pair<std::size_t, std::size_t>(0, 2) ||
            fastgatk::io::find_right_flow_clipping(1, first_hmer_flow) !=
                std::make_pair<std::size_t, std::size_t>(1, 0))
            throw std::runtime_error("FlowBasedHaplotype first-hmer clipping diverges from GATK");
        const std::vector<std::uint8_t> gatk_haplotype_bases{
            'A','T','C','G','C','A','G','G','G','A','A','T','T','G','T','C','C','C','C','A','T','G','A','A','A','C','T','A','A','G'};
        const auto gatk_haplotype_flow = fastgatk::io::encode_flow_key(gatk_haplotype_bases, "TACG");
        const std::vector<std::int32_t> gatk_expected_key{
            0,1,0,0,1,0,1,1,0,0,1,0,0,1,0,3,0,2,0,0,2,0,0,1,1,0,4,0,0,1,0,0,1,0,0,1,0,3,1,0,1,2,0,1};
        if (gatk_haplotype_flow.key != gatk_expected_key)
            throw std::runtime_error("FlowBasedHaplotype key diverges from GATK FlowBasedKeyCodec");
        // Host-side equivalent of GATK's LongHomopolymerHaplotypeCollapsingEngine
        // collapseBases(): preserve the first homopolymer and cap subsequent
        // runs at the requested threshold.  This is deliberately tested before
        // flow-key encoding because the production path collapses haplotypes at
        // that same boundary.
        const auto collapsed_hmer = fastgatk::io::collapse_flow_homopolymers(
            std::vector<std::uint8_t>{'A','A','A','A','A','C','C','C','C','C','C','G'}, 4);
        const std::vector<std::uint8_t> expected_collapsed_hmer{
            'A','A','A','A','A','C','C','C','C','G'};
        if (collapsed_hmer != expected_collapsed_hmer ||
            fastgatk::io::collapse_flow_homopolymers(expected_collapsed_hmer, 0) !=
                expected_collapsed_hmer ||
            !fastgatk::io::flow_homopolymer_exceeds_threshold(
                std::vector<std::uint8_t>{'A','A','A'}, 2) ||
            fastgatk::io::flow_homopolymer_exceeds_threshold(
                std::vector<std::uint8_t>{'A','A'}, 2))
            throw std::runtime_error("Flow homopolymer collapse API diverges from GATK");
        fastgatk::io::ReadBatch flow_batch;
        flow_batch.offsets = {0, 2};
        flow_batch.bases = flow_tags.bases;
        flow_batch.qualities = flow_tags.qualities;
        flow_batch.positions = {0};
        flow_batch.flow_tp_offsets = {0, 2};
        flow_batch.flow_tp = flow_tags.tp;
        flow_batch.flow_t0_offsets = {0, 2};
        flow_batch.flow_t0_phred = flow_tags.t0_phred;
        flow_batch.flow_order_offsets = {0, 4};
        flow_batch.flow_orders = {'A', 'C', 'G', 'T'};
        flow_batch.flow_max_hmer = {12};
        const auto batch_tags = fastgatk::io::flow_tags_from_batch(flow_batch, 0, true, 0.001, true);
        const auto decoded_batch_flow = fastgatk::io::decode_flow_read(batch_tags);
        if (decoded_batch_flow.key != decoded_flow.key ||
            decoded_batch_flow.probabilities != decoded_flow.probabilities)
            throw std::runtime_error("FlowBatch-to-codec adapter changed decoded flow data");
        const auto score = fastgatk::kernels::smith_waterman_score_reference(
            input.read.data(), input.read.size(), input.haplotype.data(), input.haplotype.size());
        if (score <= 0) throw std::runtime_error("Smith-Waterman reference API smoke failed");
        const auto alignment = fastgatk::kernels::smith_waterman_align_reference(
            input.read.data(), input.read.size(), input.haplotype.data(), input.haplotype.size());
        if (alignment.score != score || alignment.cigar != "4M" ||
            alignment.read_start != 0 || alignment.read_end != 4 ||
            alignment.reference_start != 0 || alignment.reference_end != 4)
            throw std::runtime_error("Smith-Waterman traceback/CIGAR API smoke failed");
        const std::uint8_t deletion_read[] = {'A', 'C', 'G', 'T'};
        const std::uint8_t deletion_reference[] = {'A', 'C', 'G', 'G', 'T'};
        const auto deletion = fastgatk::kernels::smith_waterman_align_reference(
            deletion_read, 4, deletion_reference, 5,
            fastgatk::kernels::SmithWatermanParameters{10, -15, -2, -1});
        const auto deletion_score = fastgatk::kernels::smith_waterman_score_reference(
            deletion_read, 4, deletion_reference, 5,
            fastgatk::kernels::SmithWatermanParameters{10, -15, -2, -1});
        if (deletion.score != deletion_score || deletion.cigar != "2M1D2M" || deletion.read_start != 0 ||
            deletion.read_end != 4 || deletion.reference_start != 0 ||
            deletion.reference_end != 5)
            throw std::runtime_error("Smith-Waterman deletion CIGAR smoke failed: " +
                                     deletion.cigar + " score=" + std::to_string(deletion.score));
        // Golden values generated by GATK 4.6.2.0's
        // SmithWatermanJavaAligner (SWParameters 10,-15,-30,-5):
        // read=TTACGTAA, reference=ACGT.
        const std::uint8_t overhang_read[] = {'T', 'T', 'A', 'C', 'G', 'T', 'A', 'A'};
        const std::uint8_t overhang_reference[] = {'A', 'C', 'G', 'T'};
        const auto softclip = fastgatk::kernels::smith_waterman_align_reference(
            overhang_read, 8, overhang_reference, 4,
            fastgatk::kernels::SmithWatermanParameters{},
            fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip);
        const auto indel = fastgatk::kernels::smith_waterman_align_reference(
            overhang_read, 8, overhang_reference, 4,
            fastgatk::kernels::SmithWatermanParameters{},
            fastgatk::kernels::SmithWatermanOverhangStrategy::Indel);
        const auto leading_indel = fastgatk::kernels::smith_waterman_align_reference(
            overhang_read, 8, overhang_reference, 4,
            fastgatk::kernels::SmithWatermanParameters{},
            fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel);
        const auto ignore = fastgatk::kernels::smith_waterman_align_reference(
            overhang_read, 8, overhang_reference, 4,
            fastgatk::kernels::SmithWatermanParameters{},
            fastgatk::kernels::SmithWatermanOverhangStrategy::Ignore);
        if (softclip.cigar != "2S4M2S" || softclip.reference_start != 0 ||
            indel.cigar != "2I4M2I" || leading_indel.cigar != "2I4M2I" ||
            ignore.cigar != "8M" || ignore.alignment_offset != -2)
            throw std::runtime_error("Smith-Waterman GATK overhang strategy smoke failed: " +
                softclip.cigar + "," + indel.cigar + "," + leading_indel.cigar + "," + ignore.cigar);

        // These vectors are copied from GATK's SmithWatermanAlignerAbstractUnitTest
        // (the Java aligner is the compatibility oracle for the host traceback).
        // Keep the complete strategy/parameter matrix here so a change to the
        // Kokkos score path cannot silently alter the CIGAR tie-break contract.
        const auto assert_sw_golden = [](const std::string& reference,
                                         const std::string& read,
                                         const fastgatk::kernels::SmithWatermanParameters parameters,
                                         const fastgatk::kernels::SmithWatermanOverhangStrategy strategy,
                                         const std::size_t expected_offset,
                                         const std::string& expected_cigar) {
            const std::vector<std::uint8_t> ref(reference.begin(), reference.end());
            const std::vector<std::uint8_t> alt(read.begin(), read.end());
            const auto result = fastgatk::kernels::smith_waterman_align_reference(
                alt.data(), alt.size(), ref.data(), ref.size(), parameters, strategy);
            if (result.alignment_offset != static_cast<int>(expected_offset) ||
                result.cigar != expected_cigar)
                throw std::runtime_error("Smith-Waterman GATK golden mismatch: offset=" +
                    std::to_string(result.alignment_offset) + " cigar=" + result.cigar +
                    " expected_offset=" + std::to_string(expected_offset) +
                    " expected_cigar=" + expected_cigar);
        };
        const fastgatk::kernels::SmithWatermanParameters original_default{3, -1, -4, -3};
        assert_sw_golden("AAAGGACTGACTG", "ACTGACTGACTG", original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip, 1, "12M");
        assert_sw_golden("AAAGACTACTG", "AACGGACACTG", {50, -100, -220, -12},
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip, 1,
                         "2M2I3M1D4M");
        assert_sw_golden("AAAGACTACTG", "AACGGACACTG", {200, -50, -300, -22},
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip, 0, "11M");
        assert_sw_golden("AAACCCCC", "CCCCCGGG", original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip, 3, "5M3S");
        const std::string substring_reference = "AAACCCCC";
        const std::string substring_read = "CCCCC";
        assert_sw_golden(substring_reference, substring_read, original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip, 3, "5M");
        assert_sw_golden(substring_reference, substring_read, original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Indel, 0, "3D5M");
        assert_sw_golden(substring_reference, substring_read, original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel, 0, "3D5M");
        assert_sw_golden(substring_reference, substring_read, original_default,
                         fastgatk::kernels::SmithWatermanOverhangStrategy::Ignore, 3, "5M");
        const std::string long_reference =
            "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        for (const auto strategy : {fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip,
                                    fastgatk::kernels::SmithWatermanOverhangStrategy::Indel,
                                    fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel,
                                    fastgatk::kernels::SmithWatermanOverhangStrategy::Ignore}) {
            const std::size_t expected_offset =
                strategy == fastgatk::kernels::SmithWatermanOverhangStrategy::Indel ||
                strategy == fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel ? 0 : 359;
            const std::string expected_cigar =
                strategy == fastgatk::kernels::SmithWatermanOverhangStrategy::Indel
                    ? "1M358D6M29D"
                    : strategy == fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel
                        ? "1M1D6M" : "7M";
            assert_sw_golden(long_reference, "AAAAAAA", original_default, strategy, expected_offset, expected_cigar);
        }
        const auto overhang_scores = fastgatk::kernels::smith_waterman_score_kokkos(
            {{{'T', 'T', 'A', 'C', 'G', 'T', 'A', 'A'},
              {'A', 'C', 'G', 'T'}}},
            fastgatk::kernels::SmithWatermanParameters{},
            fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip);
        if (overhang_scores.scores.size() != 1 ||
            overhang_scores.scores[0] != softclip.score)
            throw std::runtime_error("Smith-Waterman Kokkos/GATK overhang score mismatch");
        const auto sw_batch = fastgatk::kernels::smith_waterman_score_kokkos({
            {std::vector<std::uint8_t>{'A', 'C', 'G', 'T'},
             std::vector<std::uint8_t>{'A', 'C', 'G', 'T'}},
            {{'A', 'C', 'G'}, {'A', 'C', 'G', 'T', 'A'}}});
        if (sw_batch.scores.size() != 2 || sw_batch.scores[0] != score ||
            sw_batch.scores[1] <= 0 || sw_batch.execution_space.empty())
            throw std::runtime_error("Smith-Waterman Kokkos batch API smoke failed");
        std::vector<fastgatk::kernels::SmithWatermanRequest> sw_uniform;
        for (std::size_t request = 0; request < 8; ++request) {
            fastgatk::kernels::SmithWatermanRequest item;
            item.read.assign({'A', 'C', 'G', 'T', 'A', 'C'});
            item.reference.assign({'A', 'C', 'G', 'T', 'A', 'C', 'G'});
            item.read[request % item.read.size()] = static_cast<std::uint8_t>('T');
            sw_uniform.push_back(std::move(item));
        }
        const auto sw_uniform_result = fastgatk::kernels::smith_waterman_score_kokkos(sw_uniform);
        if (sw_uniform_result.scores.size() != sw_uniform.size() ||
            sw_uniform_result.simd_width == 0 || sw_uniform_result.execution_space.empty())
            throw std::runtime_error("Smith-Waterman SIMD telemetry failed");
        for (std::size_t request = 0; request < sw_uniform.size(); ++request) {
            const auto expected = fastgatk::kernels::smith_waterman_score_reference(
                sw_uniform[request].read.data(), sw_uniform[request].read.size(),
                sw_uniform[request].reference.data(), sw_uniform[request].reference.size());
            if (sw_uniform_result.scores[request] != expected)
                throw std::runtime_error("Smith-Waterman SIMD/reference score mismatch");
        }
        // Exercise a final partial SIMD group (production batch sizes are not
        // required to be multiples of the selected Kokkos SIMD width).  This
        // is a memory-safety regression for padded lane loads.
        auto sw_partial = sw_uniform;
        sw_partial.push_back(sw_uniform.front());
        sw_partial.push_back(sw_uniform[1]);
        sw_partial.push_back(sw_uniform[2]);
        const auto sw_partial_result = fastgatk::kernels::smith_waterman_score_kokkos(sw_partial);
        if (sw_partial_result.scores.size() != sw_partial.size())
            throw std::runtime_error("Smith-Waterman partial SIMD result size mismatch");
        for (std::size_t request = 0; request < sw_partial.size(); ++request) {
            const auto expected = fastgatk::kernels::smith_waterman_score_reference(
                sw_partial[request].read.data(), sw_partial[request].read.size(),
                sw_partial[request].reference.data(), sw_partial[request].reference.size());
            if (sw_partial_result.scores[request] != expected)
                throw std::runtime_error("Smith-Waterman partial SIMD/reference score mismatch");
        }
        std::vector<fastgatk::kernels::SmithWatermanRequest> sw_repeat(16);
        for (auto& item : sw_repeat) {
            item.read.assign({'A','A','A','A','A','A','A','A','A','A','G','A','A','A','A','A','A','A','A','A','A'});
            item.reference.assign({'A','A','A','A','A','C','A','A','A','A','A','G','A','A','A','A','G','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A'});
        }
        const auto sw_repeat_result = fastgatk::kernels::smith_waterman_score_kokkos(sw_repeat);
        const auto sw_repeat_expected = fastgatk::kernels::smith_waterman_score_reference(
            sw_repeat.front().read.data(), sw_repeat.front().read.size(),
            sw_repeat.front().reference.data(), sw_repeat.front().reference.size());
        if ((sw_repeat_result.simd_width > 1 && sw_repeat_result.simd_groups == 0) ||
            sw_repeat_result.scores.front() != sw_repeat_expected)
            throw std::runtime_error("Smith-Waterman repetitive SIMD/reference score mismatch simd=" +
                std::to_string(sw_repeat_result.scores.front()) + " ref=" +
                std::to_string(sw_repeat_expected));
        std::vector<fastgatk::kernels::SmithWatermanRequest> sw_golden;
        for (std::size_t request = 0; request < 8; ++request) {
            fastgatk::kernels::SmithWatermanRequest item;
            item.read.resize(5 + request % 4);
            item.reference.resize(7 + request % 5);
            for (std::size_t i = 0; i < item.read.size(); ++i)
                item.read[i] = static_cast<std::uint8_t>((i * 3 + request) & 3U);
            for (std::size_t i = 0; i < item.reference.size(); ++i)
                item.reference[i] = static_cast<std::uint8_t>((i + request * 2) & 3U);
            sw_golden.push_back(std::move(item));
        }
        const auto sw_golden_result = fastgatk::kernels::smith_waterman_score_kokkos(sw_golden);
        for (std::size_t i = 0; i < sw_golden.size(); ++i) {
            const auto expected = fastgatk::kernels::smith_waterman_score_reference(
                sw_golden[i].read.data(), sw_golden[i].read.size(),
                sw_golden[i].reference.data(), sw_golden[i].reference.size());
            if (sw_golden_result.scores[i] != expected)
                throw std::runtime_error("Smith-Waterman Kokkos/reference score mismatch");
        }
        const std::array<fastgatk::kernels::SmithWatermanOverhangStrategy, 4> strategies{
            fastgatk::kernels::SmithWatermanOverhangStrategy::Softclip,
            fastgatk::kernels::SmithWatermanOverhangStrategy::Indel,
            fastgatk::kernels::SmithWatermanOverhangStrategy::LeadingIndel,
            fastgatk::kernels::SmithWatermanOverhangStrategy::Ignore};
        for (const auto strategy : strategies) {
            const auto strategy_result = fastgatk::kernels::smith_waterman_score_kokkos(
                sw_golden, {}, strategy);
            for (std::size_t i = 0; i < sw_golden.size(); ++i) {
                const auto expected = fastgatk::kernels::smith_waterman_score_reference(
                    sw_golden[i].read.data(), sw_golden[i].read.size(),
                    sw_golden[i].reference.data(), sw_golden[i].reference.size(), {}, strategy);
                if (strategy_result.scores[i] != expected)
                    throw std::runtime_error("Smith-Waterman strategy score mismatch");
            }
        }
        fastgatk::kernels::KmerGraphInput graph_input{
            {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'}, {0, 8}, {},
            {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'}, {0, 8}};
        // Mutect2 enables adaptive pruning with the ReadThreadingAssembler
        // defaults expressed in natural-log units: ln(10) and ln(10^4).
        // These values govern whether a low-support chain exists when
        // dangling-end recovery runs, so merely testing an explicitly-set
        // custom threshold would not protect the native CLI default.
        const auto default_pruning_options = fastgatk::kernels::KmerGraphOptions{};
        if (std::abs(default_pruning_options.pruning_log_odds_threshold - std::log(10.0)) > 1.0e-12 ||
            std::abs(default_pruning_options.pruning_seeding_log_odds_threshold -
                     4.0 * std::log(10.0)) > 1.0e-12)
            throw std::runtime_error("GATK adaptive-pruning default thresholds drifted");
        const auto graph = fastgatk::kernels::build_kmer_graph_kokkos(
            graph_input, fastgatk::kernels::KmerGraphOptions{3, 1});
        if (!graph.used || graph.input_reads != 1 || graph.input_kmers != 6 ||
            graph.nodes == 0 || graph.edges == 0 || graph.reference_nodes == 0 ||
            graph.reference_edges == 0 || graph.reference_connected_nodes == 0 ||
            graph.reference_path_count == 0 || graph.haplotype_path_count == 0 ||
            graph.haplotype_path_sequences.empty() ||
            graph.haplotype_path_tids.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_starts.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_ends.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_support.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_has_non_reference_edge.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_alt_read_starts.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_alt_read_ends.size() != graph.haplotype_path_sequences.size() ||
            graph.haplotype_path_sequences.front().size() < 3 || graph.execution_space.empty())
            throw std::runtime_error("K-mer graph Kokkos API smoke failed: reads=" +
                std::to_string(graph.input_reads) + " kmers=" + std::to_string(graph.input_kmers) +
                " nodes=" + std::to_string(graph.nodes) + " edges=" + std::to_string(graph.edges));
        // GATK SeqGraphUnitTest.testMergeNonVariation expects a linear
        // de-Bruijn chain (GGTTAACC) to zip to one sequence vertex for every
        // valid k-mer size. The raw topology remains available above while
        // the native SeqGraph counters expose the compressed representation.
        for (std::uint32_t k = 3; k <= 7; ++k) {
            fastgatk::kernels::KmerGraphInput linear_graph_input;
            const std::string linear_sequence = "GGTTAACC";
            linear_graph_input.bases.assign(linear_sequence.begin(), linear_sequence.end());
            linear_graph_input.offsets = {0, static_cast<std::uint32_t>(linear_sequence.size())};
            linear_graph_input.reference_bases = linear_graph_input.bases;
            linear_graph_input.reference_offsets = linear_graph_input.offsets;
            const auto linear_graph = fastgatk::kernels::build_kmer_graph_kokkos(
                linear_graph_input, fastgatk::kernels::KmerGraphOptions{k, 1});
            if (linear_graph.seqgraph_nodes != 1 || linear_graph.seqgraph_edges != 0 ||
                linear_graph.nodes == 0 || linear_graph.edges == 0) {
                throw std::runtime_error("SeqGraph linear-chain compression failed for k=" +
                    std::to_string(k) + " raw=" + std::to_string(linear_graph.nodes) + "/" +
                    std::to_string(linear_graph.edges) + " compressed=" +
                    std::to_string(linear_graph.seqgraph_nodes) + "/" +
                std::to_string(linear_graph.seqgraph_edges));
        }
        // SeqGraph's sequence-level rewrites are exercised independently of
        // k-mer overlap details.  The first corpus is a true diamond with a
        // shared prefix/suffix; the second is a MergeTails corpus with two
        // sink paths sharing seven terminal bases.  The returned language
        // must still contain every input haplotype after rewriting.
        const std::vector<fastgatk::kernels::SeqGraphPath> diamond_paths{
            {"AACCGGTT", 0, 100, 108, true},
            {"AATCGGTT", 0, 100, 108, false},
            {"AAGCGGTT", 0, 100, 108, false}};
        const auto diamond = fastgatk::kernels::simplify_seqgraph_paths(diamond_paths, 4);
        if (diamond.diamond_merges == 0 || diamond.suffix_splits == 0 ||
            diamond.suffix_merges == 0 || diamond.path_sequences.size() != diamond_paths.size())
            throw std::runtime_error("SeqGraph diamond/suffix rewrite failed");
        for (const auto& path : diamond_paths)
            if (std::find(diamond.path_sequences.begin(), diamond.path_sequences.end(), path.sequence) ==
                diamond.path_sequences.end())
                throw std::runtime_error("SeqGraph diamond changed path language");
        const std::vector<fastgatk::kernels::SeqGraphPath> tail_paths{
            {"TTTACCGGTT", 0, 200, -1, false},
            {"GGGACCGGTT", 0, 200, -1, false}};
        const auto tails = fastgatk::kernels::simplify_seqgraph_paths(tail_paths, 7);
        if (tails.tail_merges == 0 || tails.suffix_splits == 0 ||
            tails.path_sequences.size() != tail_paths.size())
            throw std::runtime_error("SeqGraph tail rewrite failed");
        for (const auto& path : tail_paths)
            if (std::find(tails.path_sequences.begin(), tails.path_sequences.end(), path.sequence) ==
                tails.path_sequences.end())
                throw std::runtime_error("SeqGraph tail changed path language");
        // AdaptiveChainPruner corpus mirroring GATK's adjacent-bad-edge test:
        // a high-weight reference chain and a real 50x variant survive, while
        // two 5x error chains are removed as one chain each.
        const std::vector<fastgatk::kernels::SeqGraphEdge> pruning_edges{
            {0, 1, 100, true, {100}}, {1, 2, 100, true, {100}},
            {2, 5, 100, true, {100}}, {1, 3, 5, false, {5}},
            {3, 2, 5, false, {5}}, {1, 4, 50, false, {50}},
            {4, 2, 50, false, {50}}};
        auto pruning_options = fastgatk::kernels::KmerGraphOptions{};
        pruning_options.min_pruning = 1;
        pruning_options.num_pruning_samples = 1;
        pruning_options.use_adaptive_pruning = true;
        pruning_options.initial_error_rate_for_pruning = 0.01;
        pruning_options.pruning_log_odds_threshold = 2.0;
        pruning_options.pruning_seeding_log_odds_threshold = std::log(4.0);
        pruning_options.max_unpruned_variants = 50;
        const auto adaptive_pruning = fastgatk::kernels::prune_seqgraph_chains(
            6, pruning_edges, pruning_options);
        if (adaptive_pruning.chain_count == 0 || adaptive_pruning.adaptive_pruned_chain_count == 0 ||
            adaptive_pruning.keep_edges[3] != 0 || adaptive_pruning.keep_edges[4] != 0 ||
            adaptive_pruning.keep_edges[5] == 0 || adaptive_pruning.keep_edges[6] == 0)
            throw std::runtime_error("AdaptiveChainPruner graph corpus failed");

        // Numerical oracle values from GATK 4.6.2.0's
        // Mutect2Engine.logLikelihoodRatio(ref, alt, errorProbability). The
        // cases straddle Commons Math's exact-coefficient and explicit
        // log-sum branches, including the digamma cutoff at depth 49.
        const std::vector<std::tuple<std::uint32_t, std::uint32_t, double, double>>
            adaptive_likelihood_oracle{
                {5, 1, 0.001, 3.1634726070818580},
                {10, 2, 0.001, 7.0457378619341780},
                {30, 10, 0.001, 44.751626277016540},
                {500, 100, 0.001, 416.27067099743640},
                {529, 500, 0.001, 2736.5256988795670},
                {530, 500, 0.001, 2735.8587961669637},
                {1023, 37, 0.001, 88.519960528930030},
                {1024, 37, 0.001, 88.480673867195260},
                {50000, 100, 0.001, -204.94155734433104},
                {5, 1, 0.5, -3.2169986460726040},
                {10, 2, 0.1, -2.5115694521232728}};
        for (const auto& [n_ref, n_alt, error_probability, expected] : adaptive_likelihood_oracle) {
            const auto observed = fastgatk::kernels::seqgraph_constant_error_log_likelihood_ratio(
                n_ref, n_alt, error_probability);
            if (!std::isfinite(observed) || std::abs(observed - expected) > 2.0e-8)
                throw std::runtime_error("AdaptiveChainPruner GATK likelihood oracle drift: observed=" +
                    std::to_string(observed) + " expected=" + std::to_string(expected));
        }

        // GATK seeds the max-weight chain by its largest single edge, not by
        // the sum across a long chain.  This discriminating corpus prevents a
        // regression to the old sum-based ordering: Java keeps the isolated
        // 100x chain (edge 2), while the 60x+60x chain is pruned.
        const std::vector<fastgatk::kernels::SeqGraphEdge> max_edge_seed_edges{
            {0, 1, 60, false, {60}}, {1, 2, 60, false, {60}},
            {3, 4, 100, false, {100}}};
        auto max_edge_options = pruning_options;
        max_edge_options.pruning_log_odds_threshold = 1.0;
        max_edge_options.pruning_seeding_log_odds_threshold = 1.0;
        const auto max_edge_seed = fastgatk::kernels::prune_seqgraph_chains(
            5, max_edge_seed_edges, max_edge_options);
        if (max_edge_seed.keep_edges.size() != max_edge_seed_edges.size() ||
            max_edge_seed.keep_edges[0] != 0 || max_edge_seed.keep_edges[1] != 0 ||
            max_edge_seed.keep_edges[2] == 0)
            throw std::runtime_error("AdaptiveChainPruner max-edge seed ordering failed");
        }
        // GATK ReadThreadingGraph consumes GATKRead.getBases() in the stored
        // BAM order. FLAG 0x10 is alignment metadata, not an instruction to
        // reverse-complement the graph payload.
        fastgatk::kernels::KmerGraphInput reverse_graph_input;
        const std::string reverse_reference = "ACGTTGCAATC";
        const std::string reverse_read = reverse_reference;
        reverse_graph_input.bases.assign(reverse_read.begin(), reverse_read.end());
        reverse_graph_input.offsets = {0, static_cast<std::uint32_t>(reverse_read.size())};
        reverse_graph_input.flags = {0x10U};
        reverse_graph_input.reference_bases.assign(reverse_reference.begin(), reverse_reference.end());
        reverse_graph_input.reference_offsets = {0, static_cast<std::uint32_t>(reverse_reference.size())};
        const auto reverse_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            reverse_graph_input, fastgatk::kernels::KmerGraphOptions{5, 1});
        if (reverse_graph.reference_connected_nodes == 0 || reverse_graph.dangling_nodes != 0 ||
            reverse_graph.haplotype_path_sequences.empty()) {
            throw std::runtime_error("reverse-strand k-mer orientation failed: connected=" +
                std::to_string(reverse_graph.reference_connected_nodes) + " dangling=" +
                std::to_string(reverse_graph.dangling_nodes));
        }
        // GraphBasedKBestHaplotypeFinder returns only complete source-to-sink
        // paths. The default must therefore traverse this 320-base scaffold
        // even though the historical native 256-step cutoff would have
        // emitted a truncated prefix.
        const std::string source_sink_reference =
            "AGAGTCGGCCCTGGCGCGGAATTGCATCGGCGGGATTAAATGACGGGAGCCTCTCGGTAGCGCAACAAGTACAAACGAACTCTTAACAACGTCCCTCCGAAAGTGCCTAGCATTACGGTCAAATCTGGTTATACTCGTTGCGTGTACTTCACGGATTATGTTCGTTGAAATGACCCTGACGCTTGAGCAAAGAGCACTGGTGAATATCGGTGACATTTACTTCTAGCGGTCTATGTTCCCCACGCTTGCCAGGCTAAAATCATGGACATCACAATGAGATGTAAGCAGAACTTGGTTCCGCATGTGTCGCTTACGATCGA";
        fastgatk::kernels::KmerGraphInput long_graph_input;
        long_graph_input.bases.assign(source_sink_reference.begin(), source_sink_reference.end());
        long_graph_input.offsets = {0, static_cast<std::uint32_t>(source_sink_reference.size())};
        long_graph_input.reference_bases = long_graph_input.bases;
        long_graph_input.reference_offsets = long_graph_input.offsets;
        const auto long_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            long_graph_input, fastgatk::kernels::KmerGraphOptions{10, 1});
        if (std::find(long_graph.haplotype_path_sequences.begin(),
                      long_graph.haplotype_path_sequences.end(), source_sink_reference) ==
            long_graph.haplotype_path_sequences.end()) {
            throw std::runtime_error("K-best source-to-sink traversal truncated a long reference path");
        }
        // ReadThreadingGraphUnitTest.testSimpleHaplotypeRethreading: a
        // single SNP-containing haplotype must add an alternate branch and
        // remain discoverable in deterministic path materialization.
        const std::string rethread_ref =
            "CATGCACTTTAAAACTTGCCTTTTTAACAAGACTTCCAGATG";
        const std::string rethread_alt =
            "CATGCACTTTAAAACTTGCCGTTTTAACAAGACTTCCAGATG";
        fastgatk::kernels::KmerGraphInput rethread_input;
        rethread_input.bases.insert(rethread_input.bases.end(),
                                    rethread_ref.begin(), rethread_ref.end());
        rethread_input.bases.insert(rethread_input.bases.end(),
                                    rethread_alt.begin(), rethread_alt.end());
        rethread_input.offsets = {0, static_cast<std::uint32_t>(rethread_ref.size()),
                                  static_cast<std::uint32_t>(rethread_ref.size() + rethread_alt.size())};
        rethread_input.reference_bases.assign(rethread_ref.begin(), rethread_ref.end());
        rethread_input.reference_offsets = {0, static_cast<std::uint32_t>(rethread_ref.size())};
        const auto rethread_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            rethread_input, fastgatk::kernels::KmerGraphOptions{11, 1, 1, 64, 256, 0, false, true});
        const auto has_rethread_alt = std::any_of(
            rethread_graph.haplotype_path_has_non_reference_edge.begin(),
            rethread_graph.haplotype_path_has_non_reference_edge.end(),
            [](const auto value) { return value != 0; });
        if (rethread_graph.nodes <= rethread_ref.size() - 11 + 1 || !has_rethread_alt ||
            rethread_graph.haplotype_path_sequences.empty()) {
            throw std::runtime_error("ReadThreadingGraph SNP rethreading corpus failed: nodes=" +
                std::to_string(rethread_graph.nodes) + " paths=" +
                std::to_string(rethread_graph.haplotype_path_count) + " alt=" +
                std::to_string(has_rethread_alt));
        }
        // Direct ports of the source-to-sink cases in
        // ReadThreadingAssemblerUnitTest.  Unlike that class's raw-graph
        // diagnostics, each alternate here remains connected to both the
        // reference source and sink, so it is also part of the default caller
        // assembly language after the connected-path cleanup.
        const auto check_source_assembler_paths =
            [](const std::string& label, const std::string& reference,
               const std::vector<std::string>& reads,
               const std::vector<std::string>& expected,
               const std::uint32_t kmer_size = 3) {
                fastgatk::kernels::KmerGraphInput input;
                input.offsets.push_back(0);
                for (const auto& read : reads) {
                    input.bases.insert(input.bases.end(), read.begin(), read.end());
                    input.offsets.push_back(static_cast<std::uint32_t>(input.bases.size()));
                }
                input.reference_bases.assign(reference.begin(), reference.end());
                input.reference_offsets = {0, static_cast<std::uint32_t>(reference.size())};
                auto options = fastgatk::kernels::KmerGraphOptions{};
                options.k = kmer_size;
                options.min_count = 1;
                // This is the source TestAssembler configuration:
                // minPruneFactor=0 means do not prune by read support.
                options.min_pruning = 0;
                options.max_paths = 64;
                options.min_dangling_branch_length = 0;
                options.recover_all_dangling_branches = false;
                options.allow_non_unique_kmers_in_ref = true;
                const auto graph = fastgatk::kernels::build_kmer_graph_kokkos(input, options);
                auto actual = graph.haplotype_path_sequences;
                std::sort(actual.begin(), actual.end());
                auto sorted_expected = expected;
                std::sort(sorted_expected.begin(), sorted_expected.end());
                if (graph.reference_kmer_rejected || actual != sorted_expected) {
                    throw std::runtime_error(
                        "ReadThreadingAssembler source-to-sink path mismatch: " + label +
                        " paths=" + std::to_string(graph.haplotype_path_count));
                }
            };
        check_source_assembler_paths(
            "single full bubble", "ACAACTGA", {"ACAGCTGA"},
            {"ACAACTGA", "ACAGCTGA"});
        check_source_assembler_paths(
            "partial reads complete a bubble", "ACAACTGA", {"ACAGCT", "GCTGA"},
            {"ACAACTGA", "ACAGCTGA"});
        check_source_assembler_paths(
            "middle-start branch reconnects", "CAAAATGGGG", {"AAATCGGG"},
            {"CAAAATGGGG", "CAAAATCGGG"});
        // ReadThreadingAssemblerUnitTest.testSingleIndelAsDoubleIndel3Reads:
        // a deletion crosses two repetitive structures. The graph must keep
        // the one complete alternate path rather than splitting it into
        // incompatible partial branches or emitting recombinants.
        const std::string source_indel_reference =
            "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCTCTCTGTGTGTGTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        std::string source_indel_alternate =
            "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCT----------GTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        source_indel_alternate.erase(
            std::remove(source_indel_alternate.begin(), source_indel_alternate.end(), '-'),
            source_indel_alternate.end());
        check_source_assembler_paths(
            "repetitive single indel", source_indel_reference,
            {source_indel_alternate, source_indel_alternate},
            {source_indel_reference, source_indel_alternate}, 25);
        // MultiSampleEdge-compatible pruning: requiring two independent
        // samples must remove a branch observed only in sample 0, while the
        // reference scaffold remains intact.
        auto single_sample_branch = rethread_input;
        single_sample_branch.sample_ids = {0, 0};
        const auto multi_sample_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            single_sample_branch,
            fastgatk::kernels::KmerGraphOptions{11, 1, 1, 64, 256, 0, false, true, 3, true, 2});
        const auto multi_sample_alt = std::any_of(
            multi_sample_graph.haplotype_path_has_non_reference_edge.begin(),
            multi_sample_graph.haplotype_path_has_non_reference_edge.end(),
            [](const auto value) { return value != 0; });
        if (multi_sample_alt || multi_sample_graph.pruned_nodes == 0)
            throw std::runtime_error("MultiSampleEdge pruning sample gate failed");
        // MultiSampleEdge-like ordering: two alternate branches share the
        // same prefix/suffix, and the branch observed twice must be emitted
        // before the branch observed once (independent of input order).
        const std::string branch_ref = "AAACCCGGGTTT";
        const std::string branch_high = "AAAGCCGGGTTT";
        const std::string branch_low = "AAATCCGGGTTT";
        fastgatk::kernels::KmerGraphInput branch_input;
        const std::array<std::string, 4> branch_reads{
            branch_ref, branch_low, branch_high, branch_high};
        for (const auto& sequence : branch_reads) {
            branch_input.bases.insert(branch_input.bases.end(), sequence.begin(), sequence.end());
            branch_input.offsets.push_back(static_cast<std::uint32_t>(branch_input.bases.size()));
        }
        branch_input.offsets.insert(branch_input.offsets.begin(), 0);
        branch_input.reference_bases.assign(branch_ref.begin(), branch_ref.end());
        branch_input.reference_offsets = {0, static_cast<std::uint32_t>(branch_ref.size())};
        const auto branch_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            branch_input, fastgatk::kernels::KmerGraphOptions{3, 1, 1, 64, 256, 0, false, true});
        std::vector<std::uint32_t> alternate_supports;
        std::vector<std::string> alternate_sequences;
        for (std::size_t path = 0; path < branch_graph.haplotype_path_sequences.size(); ++path) {
            if (path < branch_graph.haplotype_path_has_non_reference_edge.size() &&
                branch_graph.haplotype_path_has_non_reference_edge[path] != 0) {
                alternate_supports.push_back(branch_graph.haplotype_path_support[path]);
                alternate_sequences.push_back(branch_graph.haplotype_path_sequences[path]);
            }
        }
        if (alternate_supports.size() < 2 || alternate_supports[0] < alternate_supports[1] ||
            alternate_supports[0] < 2) {
            throw std::runtime_error("MultiSampleEdge support ordering failed: paths=" +
                std::to_string(alternate_supports.size()) + " first=" +
                (alternate_supports.empty() ? std::string("none") : std::to_string(alternate_supports[0])) +
                " second=" + (alternate_supports.size() < 2 ? std::string("none") :
                    std::to_string(alternate_supports[1])) + " seq=" +
                (alternate_sequences.empty() ? std::string("none") : alternate_sequences.front()));
        }
        // ReadThreadingGraphUnitTest.testNonUniqueMiddle: repeated kmers in
        // one source sequence are reported as non-unique, while merely
        // seeing a kmer once in several independent reads is not a repeat.
        fastgatk::kernels::KmerGraphInput non_unique_input;
        const std::array<std::string, 2> non_unique_reads{
            "GACACGTCA", "CACGTCA"};
        for (const auto& sequence : non_unique_reads) {
            non_unique_input.bases.insert(non_unique_input.bases.end(), sequence.begin(), sequence.end());
            non_unique_input.offsets.push_back(static_cast<std::uint32_t>(non_unique_input.bases.size()));
        }
        non_unique_input.offsets.insert(non_unique_input.offsets.begin(), 0);
        const std::string non_unique_ref = "GACACACAGTCA";
        non_unique_input.reference_bases.assign(non_unique_ref.begin(), non_unique_ref.end());
        non_unique_input.reference_offsets = {0, static_cast<std::uint32_t>(non_unique_ref.size())};
        const auto non_unique_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            non_unique_input, fastgatk::kernels::KmerGraphOptions{3, 1, 1, 64, 256, 0, false, true});
        const std::set<std::string> non_unique_set(non_unique_graph.non_unique_kmers.begin(),
                                                   non_unique_graph.non_unique_kmers.end());
        if (non_unique_set != std::set<std::string>{"ACA", "CAC"}) {
            throw std::runtime_error("ReadThreadingGraph non-unique kmer corpus failed");
        }
        const auto rejected_non_unique_reference = fastgatk::kernels::build_kmer_graph_kokkos(
            non_unique_input, fastgatk::kernels::KmerGraphOptions{
                3, 1, 1, 64, 256, 0, false, true, 3, false});
        if (!rejected_non_unique_reference.reference_kmer_rejected ||
            rejected_non_unique_reference.reference_non_unique_kmers == 0 ||
            rejected_non_unique_reference.nodes != 0)
            throw std::runtime_error(
                "ReadThreadingAssembler reference non-unique kmer rejection failed");
        // ReadThreadingGraphUnitTest.testReadsCreateNonUnique: a duplicate
        // k-mer inside one read becomes non-unique for later read threading;
        // a separate read starting at that k-mer must not erase the source
        // classification or turn a per-read duplicate into a global count.
        fastgatk::kernels::KmerGraphInput read_created_non_unique_input;
        const std::array<std::string, 2> read_created_non_unique_reads{
            "GCACACGTCA", "CACGTCA"};
        read_created_non_unique_input.offsets.push_back(0);
        for (const auto& sequence : read_created_non_unique_reads) {
            read_created_non_unique_input.bases.insert(read_created_non_unique_input.bases.end(),
                                                        sequence.begin(), sequence.end());
            read_created_non_unique_input.offsets.push_back(
                static_cast<std::uint32_t>(read_created_non_unique_input.bases.size()));
        }
        const std::string read_created_non_unique_reference = "GCACGTCA";
        read_created_non_unique_input.reference_bases.assign(
            read_created_non_unique_reference.begin(), read_created_non_unique_reference.end());
        read_created_non_unique_input.reference_offsets = {
            0, static_cast<std::uint32_t>(read_created_non_unique_reference.size())};
        const auto read_created_non_unique_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            read_created_non_unique_input,
            fastgatk::kernels::KmerGraphOptions{3, 1, 1, 64, 256, 0, false, true});
        if (std::set<std::string>(read_created_non_unique_graph.non_unique_kmers.begin(),
                                  read_created_non_unique_graph.non_unique_kmers.end()) !=
            std::set<std::string>{"CAC"})
            throw std::runtime_error("ReadThreadingGraph read-created non-unique kmer contract failed");
        // ReadThreadingGraphUnitTest.testNsInReadsAreNotUsedForGraph: a
        // single N anywhere in a read makes that read unavailable to the
        // source read-threading graph.  The reference path must remain the
        // sole reference-source-to-sink haplotype even after every possible
        // N position is supplied as a read.
        const std::string n_read_reference(100, 'A');
        fastgatk::kernels::KmerGraphInput n_read_input;
        n_read_input.offsets.push_back(0);
        for (std::size_t position = 0; position < n_read_reference.size(); ++position) {
            auto read = n_read_reference;
            read[position] = 'N';
            n_read_input.bases.insert(n_read_input.bases.end(), read.begin(), read.end());
            n_read_input.offsets.push_back(static_cast<std::uint32_t>(n_read_input.bases.size()));
        }
        n_read_input.reference_bases.assign(n_read_reference.begin(), n_read_reference.end());
        n_read_input.reference_offsets = {0, static_cast<std::uint32_t>(n_read_reference.size())};
        const auto n_read_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            n_read_input, fastgatk::kernels::KmerGraphOptions{25, 1, 1, 64, 0, 0, false, true});
        if (n_read_graph.haplotype_path_count != 1U ||
            n_read_graph.haplotype_path_sequences != std::vector<std::string>{n_read_reference})
            throw std::runtime_error("ReadThreadingGraph N-read exclusion contract failed: paths=" +
                std::to_string(n_read_graph.haplotype_path_count) + " nodes=" +
                std::to_string(n_read_graph.nodes));
        // A read consisting solely of repeated kmers must not be globally
        // collapsed into a de-Bruijn self-loop. ReadThreadingGraph marks its
        // kmers non-unique, cannot select a threading start, and therefore
        // leaves the reference graph acyclic. This guards the distinction
        // between a real alternate topology cycle and the old synthetic
        // string-key cycle.
        fastgatk::kernels::KmerGraphInput repeated_kmer_input{
            {'A','C','G','A','C','G','A','C','G','A','C','G','A','C','G','A','C','G',
             'A','C','G','A','C','G','A','C','G','A','C','G','A','C','G','A','C','G',
             'A','C','G','A','C','G','A','C','G','A','C','G'}, {0, 48}, {},
            {'A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T',
             'A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T',
             'A','C','G','T','A','C','G','T','A','C','G','T'}, {0, 52}};
        const auto repeated_kmer_graph = fastgatk::kernels::build_kmer_graph_kokkos(
            repeated_kmer_input, fastgatk::kernels::KmerGraphOptions{3, 1, 1, 64, 256, 0, false, false});
        if (!repeated_kmer_graph.used || repeated_kmer_graph.has_non_reference_cycles ||
            repeated_kmer_graph.kmer_iterations != 1)
            throw std::runtime_error("ReadThreadingGraph repeated-kmer false-cycle guard failed: iterations=" +
                std::to_string(repeated_kmer_graph.kmer_iterations) + "/" +
                std::to_string(repeated_kmer_graph.kmer_size) + "/" +
                std::to_string(repeated_kmer_graph.has_non_reference_cycles) + "/" +
                std::to_string(repeated_kmer_graph.nodes) + "/" +
                std::to_string(repeated_kmer_graph.edges));
        // Direct port of ReadThreadingGraphUnitTest.testCyclesInGraph.  The
        // b37 chr20 source fixture has one SNP whose overlapping 100-base
        // reads form a non-reference cycle at k=25 but not at k=75.  Disable
        // retries here so each assertion observes the exact GATK graph
        // attempt rather than the caller's later k-mer escalation policy.
        const std::string gatk_cycle_reference = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";
        const std::string gatk_cycle_alternate = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";
        fastgatk::kernels::KmerGraphInput gatk_cycle_input;
        gatk_cycle_input.offsets.push_back(0);
        for (std::size_t begin = 0; begin + 100U < gatk_cycle_alternate.size(); begin += 20U) {
            gatk_cycle_input.bases.insert(gatk_cycle_input.bases.end(),
                                          gatk_cycle_alternate.begin() +
                                              static_cast<std::ptrdiff_t>(begin),
                                          gatk_cycle_alternate.begin() +
                                              static_cast<std::ptrdiff_t>(begin + 100U));
            gatk_cycle_input.offsets.push_back(
                static_cast<std::uint32_t>(gatk_cycle_input.bases.size()));
        }
        gatk_cycle_input.reference_bases.assign(gatk_cycle_reference.begin(),
                                                gatk_cycle_reference.end());
        gatk_cycle_input.reference_offsets = {
            0, static_cast<std::uint32_t>(gatk_cycle_reference.size())};
        const auto gatk_cycle_k25 = fastgatk::kernels::build_kmer_graph_kokkos(
            gatk_cycle_input,
            fastgatk::kernels::KmerGraphOptions{25, 1, 1, 64, 0, 0, false, true});
        const auto gatk_cycle_k75 = fastgatk::kernels::build_kmer_graph_kokkos(
            gatk_cycle_input,
            fastgatk::kernels::KmerGraphOptions{75, 1, 1, 64, 0, 0, false, true});
        if (gatk_cycle_k25.reference_kmer_rejected ||
            !gatk_cycle_k25.has_non_reference_cycles ||
            gatk_cycle_k75.reference_kmer_rejected ||
            gatk_cycle_k75.has_non_reference_cycles)
            throw std::runtime_error("GATK ReadThreadingGraph cycle-size contract failed");
        fastgatk::kernels::KmerGraphInput dangling_input{
            {'C', 'C', 'C', 'G', 'G', 'A', 'T'}, {0, 7}, {},
            {'A', 'A', 'A', 'T', 'T', 'T', 'G'}, {0, 7}};
        const auto dangling = fastgatk::kernels::build_kmer_graph_kokkos(
            dangling_input, fastgatk::kernels::KmerGraphOptions{3, 1});
        if (dangling.reference_connected_nodes != 0 || dangling.dangling_nodes == 0 ||
            dangling.pruned_nodes == 0)
            throw std::runtime_error("K-mer graph dangling-branch classification failed");
        const auto recovered_dangling = fastgatk::kernels::build_kmer_graph_kokkos(
            dangling_input, fastgatk::kernels::KmerGraphOptions{
                3, 1, 1, 64, 256, 4, true});
        // GATK's recoverAllDanglingBranches only changes `findPath` inside
        // reference-connected graph components.  The subsequent
        // removePathsNotConnectedToRef() still removes this component, so it
        // must not become an unanchored or globally SW-rescued haplotype.
        if (recovered_dangling.dangling_nodes == 0 ||
            recovered_dangling.pruned_nodes == 0 ||
            recovered_dangling.haplotype_path_count != dangling.haplotype_path_count ||
            recovered_dangling.haplotype_path_sequences != dangling.haplotype_path_sequences)
            throw std::runtime_error("GATK recover-all disconnected-branch removal failed: " +
                std::to_string(dangling.haplotype_path_count) + "/" +
                std::to_string(recovered_dangling.haplotype_path_count) +
                " dangling=" + std::to_string(recovered_dangling.dangling_nodes) +
                " pruned=" + std::to_string(recovered_dangling.pruned_nodes));
        // A disconnected branch with no shared k-mer is not a dangling-end
        // recovery candidate in AbstractReadThreadingGraph.  GATK never
        // globally aligns it to the reference to manufacture a graph edge.
        fastgatk::kernels::KmerGraphInput recoverable_dangling_input{
            {'A','A','T','A','A','T'}, {0, 6}, {},
            {'A','A','A','A','A','A'}, {0, 6}, {0}, {100}, {106}};
        const auto recoverable_dangling = fastgatk::kernels::build_kmer_graph_kokkos(
            recoverable_dangling_input, fastgatk::kernels::KmerGraphOptions{
                3, 1, 1, 64, 256, 0, true, true});
        if (recoverable_dangling.dangling_nodes == 0 ||
            recoverable_dangling.dangling_recovered_paths != 0 ||
            std::any_of(recoverable_dangling.haplotype_path_sequences.begin(),
                        recoverable_dangling.haplotype_path_sequences.end(),
                        [](const auto& sequence) { return sequence.find("AAT") != std::string::npos; }))
            throw std::runtime_error("GATK disconnected branch was spuriously recovered: dangling=" +
                std::to_string(recoverable_dangling.dangling_nodes) +
                " recovered=" + std::to_string(recoverable_dangling.dangling_recovered_paths) +
                " paths=" + std::to_string(recoverable_dangling.haplotype_path_count));
        // AbstractReadThreadingGraph.mergeDanglingTail uses a LeadingIndel
        // CIGAR and reconnects at the exact matching suffix.  With k=5 this
        // read has only three reference-matching bases after its insertion,
        // so it cannot rejoin through a k-mer and must be recovered as a
        // dangling tail.  The recovered language contains the insertion and
        // the downstream reference suffix, rather than an artificial terminal
        // deletion.  This is the one-gap case that the old exact-offset-only
        // recovery could not represent.
        const std::string gapped_tail_reference = "ACGTTGCATAGCTACGATCTGCA";
        const std::string gapped_tail_read = "ACGTTGCATATGCT";
        const std::string gapped_tail_expected = "ACGTTGCATATGCTACGATCTGCA";
        fastgatk::kernels::KmerGraphInput gapped_tail_input{
            std::vector<std::uint8_t>(gapped_tail_read.begin(), gapped_tail_read.end()),
            {0, static_cast<std::uint32_t>(gapped_tail_read.size())}, {},
            std::vector<std::uint8_t>(gapped_tail_reference.begin(), gapped_tail_reference.end()),
            {0, static_cast<std::uint32_t>(gapped_tail_reference.size())}, {0}, {200},
            {200 + static_cast<std::int32_t>(gapped_tail_reference.size())}};
        const auto gapped_tail = fastgatk::kernels::build_kmer_graph_kokkos(
            gapped_tail_input, fastgatk::kernels::KmerGraphOptions{
                5, 1, 1, 64, 256, 0, false, false, 3, false});
        const auto gapped_tail_path = std::find(
            gapped_tail.haplotype_path_sequences.begin(),
            gapped_tail.haplotype_path_sequences.end(), gapped_tail_expected);
        if (gapped_tail_path == gapped_tail.haplotype_path_sequences.end() ||
            gapped_tail.dangling_recovered_paths != 1U) {
            throw std::runtime_error(
                "GATK gapped dangling-tail recovery failed: recovered=" +
                std::to_string(gapped_tail.dangling_recovered_paths) + " paths=" +
                std::to_string(gapped_tail.haplotype_path_count));
        }
        const auto gapped_tail_index = static_cast<std::size_t>(
            gapped_tail_path - gapped_tail.haplotype_path_sequences.begin());
        if (gapped_tail_index >= gapped_tail.haplotype_path_alt_read_ends.size() ||
            gapped_tail.haplotype_path_alt_read_ends[gapped_tail_index] >=
                gapped_tail_expected.size()) {
            throw std::runtime_error(
                "GATK gapped dangling-tail merge did not retain the branch boundary");
        }
        // ReadThreadingAssembler.assemble() builds every explicit k-mer graph
        // before it decides whether automatic k+10 expansion is necessary.
        // The first (k=5) graph already has an ALT here; a retry-style loop
        // would stop early and silently skip the explicit k=7 graph.  GATK
        // subsequently de-duplicates its paths by sequence, so the stable
        // result retains the k=5 coordinate metadata while reporting both
        // actual graph attempts.
        auto explicit_kmer_list_options = fastgatk::kernels::KmerGraphOptions{
            5, 1, 1, 64, 0, 0, false, true};
        explicit_kmer_list_options.requested_kmer_sizes = {5, 7};
        const auto explicit_kmer_list = fastgatk::kernels::build_kmer_graph_kokkos(
            gapped_tail_input, explicit_kmer_list_options);
        if (explicit_kmer_list.kmer_iterations != 2U ||
            explicit_kmer_list.kmer_size != 5U ||
            std::find(explicit_kmer_list.haplotype_path_sequences.begin(),
                      explicit_kmer_list.haplotype_path_sequences.end(),
                      gapped_tail_expected) == explicit_kmer_list.haplotype_path_sequences.end()) {
            throw std::runtime_error(
                "GATK explicit multi-kmer graph assembly contract failed: attempts=" +
                std::to_string(explicit_kmer_list.kmer_iterations) + " selected=" +
                std::to_string(explicit_kmer_list.kmer_size) + " paths=" +
                std::to_string(explicit_kmer_list.haplotype_path_count));
        }
        // `allowNonUniqueKmersInRef` allows this graph build, but does not
        // skip ReadThreadingAssembler.recoverDanglingTails().  Put a repeated
        // 5-mer outside the recovery join and assert that the same one-gap
        // alternate reconnects to the now longer reference suffix.
        const std::string non_unique_gapped_tail_reference =
            gapped_tail_reference + "ACGTT";
        const std::string non_unique_gapped_tail_expected =
            gapped_tail_expected + "ACGTT";
        fastgatk::kernels::KmerGraphInput non_unique_gapped_tail_input{
            std::vector<std::uint8_t>(gapped_tail_read.begin(), gapped_tail_read.end()),
            {0, static_cast<std::uint32_t>(gapped_tail_read.size())}, {},
            std::vector<std::uint8_t>(non_unique_gapped_tail_reference.begin(),
                                      non_unique_gapped_tail_reference.end()),
            {0, static_cast<std::uint32_t>(non_unique_gapped_tail_reference.size())}, {0}, {250},
            {250 + static_cast<std::int32_t>(non_unique_gapped_tail_reference.size())}};
        const auto non_unique_gapped_tail = fastgatk::kernels::build_kmer_graph_kokkos(
            non_unique_gapped_tail_input, fastgatk::kernels::KmerGraphOptions{
                5, 1, 1, 64, 256, 0, false, false, 3, true});
        const auto non_unique_gapped_tail_path = std::find(
            non_unique_gapped_tail.haplotype_path_sequences.begin(),
            non_unique_gapped_tail.haplotype_path_sequences.end(),
            non_unique_gapped_tail_expected);
        if (non_unique_gapped_tail.reference_non_unique_kmers == 0 ||
            non_unique_gapped_tail.reference_kmer_rejected ||
            non_unique_gapped_tail_path == non_unique_gapped_tail.haplotype_path_sequences.end() ||
            non_unique_gapped_tail.dangling_recovered_paths != 1U) {
            throw std::runtime_error(
                "GATK non-unique-reference dangling-tail recovery failed: nonunique=" +
                std::to_string(non_unique_gapped_tail.reference_non_unique_kmers) +
                " rejected=" + std::to_string(non_unique_gapped_tail.reference_kmer_rejected) +
                " recovered=" + std::to_string(non_unique_gapped_tail.dangling_recovered_paths) +
                " paths=" + std::to_string(non_unique_gapped_tail.haplotype_path_count));
        }
        // GATK's LowWeightChainPruner retains a heterogeneous 2/2/2/1 tail
        // because the chain has high-support edges. findPath() clears the
        // final one-read suffix, then continues upstream and reconnects the
        // remaining path at the reference. This is verified directly against
        // ReadThreadingGraph with k=5, min-pruning=2 (not inferred from the
        // chain-pruner comment alone).
        const std::string mixed_support_tail_read = "ACGTTGCATATGC";
        fastgatk::kernels::KmerGraphInput mixed_support_tail_input{
            std::vector<std::uint8_t>(mixed_support_tail_read.begin(),
                                      mixed_support_tail_read.end()),
            {0, static_cast<std::uint32_t>(mixed_support_tail_read.size()),
             static_cast<std::uint32_t>(mixed_support_tail_read.size() +
                                        gapped_tail_read.size())}, {},
            std::vector<std::uint8_t>(gapped_tail_reference.begin(),
                                      gapped_tail_reference.end()),
            {0, static_cast<std::uint32_t>(gapped_tail_reference.size())}, {0}, {275},
            {275 + static_cast<std::int32_t>(gapped_tail_reference.size())}};
        mixed_support_tail_input.bases.insert(mixed_support_tail_input.bases.end(),
                                              gapped_tail_read.begin(), gapped_tail_read.end());
        const auto mixed_support_tail = fastgatk::kernels::build_kmer_graph_kokkos(
            mixed_support_tail_input, fastgatk::kernels::KmerGraphOptions{
                5, 1, 2, 64, 256, 0, false, false, -1, false});
        if (mixed_support_tail.dangling_recovered_paths != 1U ||
            std::find(mixed_support_tail.haplotype_path_sequences.begin(),
                      mixed_support_tail.haplotype_path_sequences.end(),
                      gapped_tail_expected) ==
                mixed_support_tail.haplotype_path_sequences.end()) {
            throw std::runtime_error(
                "GATK heterogeneous dangling-tail recovery failed: recovered=" +
                std::to_string(mixed_support_tail.dangling_recovered_paths) +
                " paths=" + std::to_string(mixed_support_tail.haplotype_path_count));
        }
        // MutectReadThreadingAssemblerArgumentCollection uses adaptive
        // chain pruning, then passes pruneFactor=0 to dangling-end recovery.
        // The singleton suffix must therefore stay available to the CIGAR
        // merge even though the CLI min-pruning floor remains two. This is
        // deliberately distinct from the fixed-pruning corpus above: using
        // min-pruning here would truncate the tail before it can rejoin.
        auto adaptive_mixed_support_options = fastgatk::kernels::KmerGraphOptions{
            5, 1, 2, 64, 256, 0, false, false, -1, false};
        adaptive_mixed_support_options.use_adaptive_pruning = true;
        const auto adaptive_mixed_support_tail = fastgatk::kernels::build_kmer_graph_kokkos(
            mixed_support_tail_input, adaptive_mixed_support_options);
        if (adaptive_mixed_support_tail.dangling_recovered_paths != 1U ||
            std::find(adaptive_mixed_support_tail.haplotype_path_sequences.begin(),
                      adaptive_mixed_support_tail.haplotype_path_sequences.end(),
                      gapped_tail_expected) ==
                adaptive_mixed_support_tail.haplotype_path_sequences.end()) {
            throw std::runtime_error(
                "Mutect adaptive dangling-tail prune-factor contract failed: recovered=" +
                std::to_string(adaptive_mixed_support_tail.dangling_recovered_paths) +
                " paths=" + std::to_string(adaptive_mixed_support_tail.haplotype_path_count));
        }
        // Direct port of ReadThreadingGraphUnitTest.testForkedDanglingEnds:
        // two incomplete alternate tails share their first divergence, then
        // fork again. GATK's recoverAll walk follows the heaviest incoming
        // edge only while finding a valid chain; each terminal branch still
        // gets its own CIGAR-selected rejoin and survives SeqGraph/K-best
        // materialization.
        const std::string forked_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        const std::string forked_reference = forked_prefix + "GCTAGCTAATCG";
        const std::string forked_alt1 = forked_prefix + "ACTAGCTAATCG";
        const std::string forked_alt2 = forked_prefix + "ACTAGATAATCG";
        fastgatk::kernels::KmerGraphInput forked_tail_input;
        forked_tail_input.offsets.push_back(0);
        for (const auto* sequence : {&forked_alt1, &forked_alt2}) {
            forked_tail_input.bases.insert(forked_tail_input.bases.end(),
                                           sequence->begin(), sequence->end());
            forked_tail_input.offsets.push_back(
                static_cast<std::uint32_t>(forked_tail_input.bases.size()));
        }
        forked_tail_input.reference_bases.assign(forked_reference.begin(), forked_reference.end());
        forked_tail_input.reference_offsets = {
            0, static_cast<std::uint32_t>(forked_reference.size())};
        const auto forked_tails = fastgatk::kernels::build_kmer_graph_kokkos(
            forked_tail_input,
            fastgatk::kernels::KmerGraphOptions{15, 1, 1, 64, 0, 4, true, false, -1, false});
        for (const auto* expected : {&forked_reference, &forked_alt1, &forked_alt2}) {
            if (std::find(forked_tails.haplotype_path_sequences.begin(),
                          forked_tails.haplotype_path_sequences.end(), *expected) ==
                forked_tails.haplotype_path_sequences.end()) {
                throw std::runtime_error("GATK forked dangling-tail path missing: recovered=" +
                    std::to_string(forked_tails.dangling_recovered_paths) +
                    " paths=" + std::to_string(forked_tails.haplotype_path_count));
            }
        }
        // Direct port of ReadThreadingGraphUnitTest.
        // testForkedDanglingEndsWithSuffixCode.  The two tails fork from the
        // same k-mer.  GATK proves that the long `...GCCGATGGCT` tail has a
        // superficially legal SW CIGAR but *zero* usable suffix match, so
        // mergeDanglingTail() must reject it.  The source test intentionally
        // invokes the default `recoverAll=false` walk: allowing recover-all
        // is a different GATK mode, whose heaviest-edge traversal can select
        // a different reference path.  This guards against turning a CIGAR
        // score into an unconditional global-alignment rescue of the long
        // tail in the default mode.
        const std::string suffix_reference = forked_prefix + "GCTAGCTAATCGTTAAGCTTTAAC";
        const std::string suffix_alt1 = forked_prefix + "GCTAGCTAAGGCG";
        const std::string suffix_alt2 = forked_prefix + "GCTAGCTAAGCCGATGGCT";
        const std::string suffix_alt2_spurious_rejoined = suffix_alt2 + "TTAAGCTTTAAC";
        fastgatk::kernels::KmerGraphInput suffix_tail_input;
        suffix_tail_input.offsets.push_back(0);
        for (const auto* sequence : {&suffix_alt2, &suffix_alt1}) {
            suffix_tail_input.bases.insert(suffix_tail_input.bases.end(),
                                           sequence->begin(), sequence->end());
            suffix_tail_input.offsets.push_back(
                static_cast<std::uint32_t>(suffix_tail_input.bases.size()));
        }
        suffix_tail_input.reference_bases.assign(suffix_reference.begin(),
                                                 suffix_reference.end());
        suffix_tail_input.reference_offsets = {
            0, static_cast<std::uint32_t>(suffix_reference.size())};
        const auto suffix_tails = fastgatk::kernels::build_kmer_graph_kokkos(
            suffix_tail_input,
            fastgatk::kernels::KmerGraphOptions{15, 1, 1, 64, 0, 2, false, false, 1, false});
        const auto contains_suffix_path = [&](const std::string& sequence) {
            return std::find(suffix_tails.haplotype_path_sequences.begin(),
                             suffix_tails.haplotype_path_sequences.end(), sequence) !=
                suffix_tails.haplotype_path_sequences.end();
        };
        // GATK's complete default pass may still attempt the short sibling
        // against the forked path; the source assertion is specifically that
        // the long `GCCGATGGCT` tail does not acquire a reference rejoin.
        if (suffix_tails.dangling_recovered_paths != 1U ||
            !contains_suffix_path(suffix_reference) ||
            contains_suffix_path(suffix_alt2_spurious_rejoined)) {
            throw std::runtime_error(
                "GATK forked suffix-tail recovery contract failed: recovered=" +
                std::to_string(suffix_tails.dangling_recovered_paths) +
                " paths=" + std::to_string(suffix_tails.haplotype_path_count));
        }
        // Source-equivalent subset of
        // ReadThreadingGraphUnitTest.makeDanglingTailsData.  These run the
        // complete native recover-and-clean path rather than only a local SW
        // helper, so a positive case must create exactly one Host topology
        // edge and a rejected CIGAR must leave no alternate source-to-sink
        // path.  `CCCC C` is the valid-DNA counterpart of the source's
        // intentionally nonmatching XXXXX whole-insertion probe.
        const auto check_source_dangling_tail =
            [&](const std::string& label, const std::string& reference_end,
                const std::string& read_end, const bool should_recover) {
                const auto reference = forked_prefix + reference_end;
                const auto read = forked_prefix + read_end;
                fastgatk::kernels::KmerGraphInput input;
                input.bases.assign(read.begin(), read.end());
                input.offsets = {0, static_cast<std::uint32_t>(read.size())};
                input.reference_bases.assign(reference.begin(), reference.end());
                input.reference_offsets = {0, static_cast<std::uint32_t>(reference.size())};
                const auto graph = fastgatk::kernels::build_kmer_graph_kokkos(
                    input, fastgatk::kernels::KmerGraphOptions{
                        15, 1, 1, 64, 0, 4, false, false, -1, true});
                const auto contains_reference = std::find(graph.haplotype_path_sequences.begin(),
                    graph.haplotype_path_sequences.end(), reference) !=
                    graph.haplotype_path_sequences.end();
                if (!contains_reference ||
                    graph.dangling_recovered_paths != (should_recover ? 1U : 0U) ||
                    graph.haplotype_path_count != (should_recover ? 2U : 1U)) {
                    throw std::runtime_error(
                        "GATK dangling-tail source matrix failed: " + label +
                        " recovered=" + std::to_string(graph.dangling_recovered_paths) +
                        " paths=" + std::to_string(graph.haplotype_path_count));
                }
            };
        check_source_dangling_tail("incomplete haplotype", "AAAAAAAAAA", "CAAA", true);
        check_source_dangling_tail("insertion", "AAAAAAAAAA", "CAAAAAAAAAA", true);
        check_source_dangling_tail("deletion", "CCAAAAAAAAAA", "AAAAAAAAAA", true);
        check_source_dangling_tail("multiple SNPs", "AAAAAAAA", "CAAGATAA", true);
        check_source_dangling_tail("too little data", "AAAAA", "CA", false);
        // GATK accepts this 8M CIGAR but its final base mismatch produces a
        // zero usable suffix, so mergeDanglingTail() returns zero.
        check_source_dangling_tail("ends in mismatch", "AAAAAAA", "CAAAAAC", false);
        check_source_dangling_tail("complex CIGAR", "AAAAAA", "CGAAAACGAA", false);
        check_source_dangling_tail("whole-tail insertion", "AAAAA", "CCCCC", false);
        // AbstractReadThreadingGraph.mergeDanglingHead is a topology rewrite,
        // not merely a coordinate annotation: after aligning the reversed
        // head and reference prefix with LEADING_INDEL, it adds the upstream
        // reference-to-alternate edge.  This source-derived one-SNP case has
        // exactly two matching bases at the head merge boundary, so the
        // GATK default-like threshold of two recovers the full alternate
        // haplotype while preserving the reference path.
        const std::string gapped_head_reference = "TGTAGCAAACCGGTTACGT";
        const std::string gapped_head_read = "AATCGGTTACGT";
        // `extendDanglingPathAgainstReference` removes the original dangling
        // source before appending its synthetic k-mer chain.  The recovery
        // edge must enter the CIGAR-selected appended node (not its
        // predecessor), yielding the full source-derived alternate path.
        const std::string gapped_head_expected = "TGTAGCAAATCGGTTACGT";
        fastgatk::kernels::KmerGraphInput gapped_head_input{
            std::vector<std::uint8_t>(gapped_head_read.begin(), gapped_head_read.end()),
            {0, static_cast<std::uint32_t>(gapped_head_read.size())}, {},
            std::vector<std::uint8_t>(gapped_head_reference.begin(), gapped_head_reference.end()),
            {0, static_cast<std::uint32_t>(gapped_head_reference.size())}, {0}, {300},
            {300 + static_cast<std::int32_t>(gapped_head_reference.size())}};
        const auto gapped_head = fastgatk::kernels::build_kmer_graph_kokkos(
            gapped_head_input, fastgatk::kernels::KmerGraphOptions{
                5, 1, 1, 64, 256, 0, false, false, 2, false});
        const auto gapped_head_path = std::find(
            gapped_head.haplotype_path_sequences.begin(),
            gapped_head.haplotype_path_sequences.end(), gapped_head_expected);
        if (gapped_head_path == gapped_head.haplotype_path_sequences.end() ||
            gapped_head.dangling_recovered_paths == 0) {
            throw std::runtime_error(
                "GATK dangling-head graph-edge recovery failed: recovered=" +
                std::to_string(gapped_head.dangling_recovered_paths) + " paths=" +
                std::to_string(gapped_head.haplotype_path_count) + " seq0=" +
                (gapped_head.haplotype_path_sequences.empty() ? std::string("none") :
                    gapped_head.haplotype_path_sequences.front()) + " seq1=" +
                (gapped_head.haplotype_path_sequences.size() < 2 ? std::string("none") :
                    gapped_head.haplotype_path_sequences[1]));
        }
        // The graph-attempt opt-in for a repeated reference k-mer is not a
        // dangling-head opt-out.  ReadThreadingAssembler runs both recovery
        // passes after accepting this graph, so the one-SNP head above must
        // still create its recovery edge and retain the full path.
        const std::string non_unique_gapped_head_reference =
            gapped_head_reference + "TGTAG";
        const std::string non_unique_gapped_head_expected =
            gapped_head_expected + "TGTAG";
        fastgatk::kernels::KmerGraphInput non_unique_gapped_head_input{
            std::vector<std::uint8_t>(gapped_head_read.begin(), gapped_head_read.end()),
            {0, static_cast<std::uint32_t>(gapped_head_read.size())}, {},
            std::vector<std::uint8_t>(non_unique_gapped_head_reference.begin(),
                                      non_unique_gapped_head_reference.end()),
            {0, static_cast<std::uint32_t>(non_unique_gapped_head_reference.size())}, {0}, {350},
            {350 + static_cast<std::int32_t>(non_unique_gapped_head_reference.size())}};
        const auto non_unique_gapped_head = fastgatk::kernels::build_kmer_graph_kokkos(
            non_unique_gapped_head_input, fastgatk::kernels::KmerGraphOptions{
                5, 1, 1, 64, 256, 0, false, false, 2, true});
        if (non_unique_gapped_head.reference_non_unique_kmers == 0 ||
            non_unique_gapped_head.reference_kmer_rejected ||
            std::find(non_unique_gapped_head.haplotype_path_sequences.begin(),
                      non_unique_gapped_head.haplotype_path_sequences.end(),
                      non_unique_gapped_head_expected) ==
                non_unique_gapped_head.haplotype_path_sequences.end() ||
            non_unique_gapped_head.dangling_recovered_paths == 0U) {
            throw std::runtime_error(
                "GATK non-unique-reference dangling-head recovery failed: nonunique=" +
                std::to_string(non_unique_gapped_head.reference_non_unique_kmers) +
                " rejected=" + std::to_string(non_unique_gapped_head.reference_kmer_rejected) +
                " recovered=" + std::to_string(non_unique_gapped_head.dangling_recovered_paths) +
                " paths=" + std::to_string(non_unique_gapped_head.haplotype_path_count));
        }
        // Port the decision matrix in GATK's
        // ReadThreadingGraphUnitTest.makeDanglingHeadsData.  The source test
        // uses X/Y placeholder bytes; these valid-DNA equivalents preserve
        // the same CIGAR and merge outcomes.  A C homopolymer represents the
        // source's repeated X prefix, hence non-unique reference kmers must
        // be accepted just as ReadThreadingGraph itself accepts them.
        const auto check_source_dangling_head =
            [](const std::string& label, const std::string& reference,
               const std::string& read, const std::int32_t min_matching_bases,
               const bool should_merge, const std::string& expected_alternate) {
                fastgatk::kernels::KmerGraphInput input;
                input.bases.assign(read.begin(), read.end());
                input.offsets = {0, static_cast<std::uint32_t>(read.size())};
                input.reference_bases.assign(reference.begin(), reference.end());
                input.reference_offsets = {0, static_cast<std::uint32_t>(reference.size())};
                const auto graph = fastgatk::kernels::build_kmer_graph_kokkos(
                    input, fastgatk::kernels::KmerGraphOptions{
                        5, 1, 1, 64, 256, 0, false, false, min_matching_bases, true});
                const auto contains = [&](const std::string& sequence) {
                    return std::find(graph.haplotype_path_sequences.begin(),
                                     graph.haplotype_path_sequences.end(), sequence) !=
                        graph.haplotype_path_sequences.end();
                };
                if (graph.reference_kmer_rejected || !contains(reference) ||
                    graph.haplotype_path_count != (should_merge ? 2U : 1U) ||
                    (should_merge && (!contains(expected_alternate) ||
                                      graph.dangling_recovered_paths != 1U)) ||
                    (!should_merge && graph.dangling_recovered_paths != 0U)) {
                    throw std::runtime_error(
                        "GATK dangling-head source matrix failed: " + label +
                        " recovered=" + std::to_string(graph.dangling_recovered_paths) +
                        " paths=" + std::to_string(graph.haplotype_path_count));
                }
            };
        check_source_dangling_head("exact threshold rejects one SNP",
                                   "TGTAGCAAACCGGTTACGT", "AATCGGTTACGT", 3,
                                   false, "");
        check_source_dangling_head("near-reference-start rejection",
                                   "CCCAACCGGTTACGT", "CAAACCGGTTACGT", 0,
                                   false, "");
        check_source_dangling_head("legacy deletion rejection",
                                   "CCCCCCCAACCGGTTACGT", "CAACGGTTACGT", -1,
                                   false, "");
        check_source_dangling_head("exact deletion recovery",
                                   "CCCCCCCAACCGGTTACGT", "CAACGGTTACGT", 1,
                                   true, "CCCCCCCAACGGTTACGT");
        check_source_dangling_head("legacy insertion rejection",
                                   "CCCCCCCAACCGGTTACGT", "CAACTCGGTTACGT", -1,
                                   false, "");
        check_source_dangling_head("exact insertion recovery",
                                   "CCCCCCCAACCGGTTACGT", "CAACTCGGTTACGT", 1,
                                   true, "CCCCCCCAACTCGGTTACGT");
        check_source_dangling_head("legacy multi-SNP rejection",
                                   "TGTAGCAAACCGGTTACGT", "ATTCGGTTACGT", -1,
                                   false, "");
        check_source_dangling_head("exact multi-SNP recovery",
                                   "TGTAGCAAACCGGTTACGT", "ATTCGGTTACGT", 1,
                                   true, "TGTAGCAATTCGGTTACGT");
        fastgatk::kernels::KmerGraphInput pruning_input{
            {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'A'}, {0, 8}, {},
            {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'}, {0, 8}};
        const auto min_pruning = fastgatk::kernels::build_kmer_graph_kokkos(
            pruning_input, fastgatk::kernels::KmerGraphOptions{3, 1, 2});
        const auto retained_alternate = std::any_of(
            min_pruning.haplotype_path_has_non_reference_edge.begin(),
            min_pruning.haplotype_path_has_non_reference_edge.end(),
            [](const auto value) { return value != 0; });
        if (retained_alternate || min_pruning.pruned_nodes == 0)
            throw std::runtime_error("K-mer graph --min-pruning support floor failed");
        // Linked-de-Bruijn mode keeps the raw graph (no SeqGraph zipping) and
        // exposes the uncovered-edge artificial-haplotype recovery switch.
        fastgatk::kernels::KmerGraphInput linked_input;
        const std::string linked_reference = "ACGTTGCATAGCTACGATCTGCA";
        // 192 reads exercise a bounded multi-sample assembly corpus rather
        // than only a three-read toy graph; both alternate branches occur in
        // both samples so the sample gate cannot hide a real branch.
        std::vector<std::string> linked_reads;
        std::vector<std::int32_t> linked_sample_ids;
        linked_reads.reserve(192);
        linked_sample_ids.reserve(192);
        for (std::size_t read = 0; read < 192; ++read) {
            linked_reads.push_back(read % 6 == 0 ? linked_reference
                : ((read / 2U) & 1U) != 0 ? "ACGTTGCATATCTACGATCTGCA"
                                           : "ACGTTGCATAGCTCCGATCTGCA");
            linked_sample_ids.push_back(static_cast<std::int32_t>(read & 1U));
        }
        linked_input.reference_bases.assign(linked_reference.begin(), linked_reference.end());
        linked_input.reference_offsets = {0, static_cast<std::uint32_t>(linked_reference.size())};
        linked_input.offsets.push_back(0);
        for (std::size_t read = 0; read < linked_reads.size(); ++read) {
            const auto& sequence = linked_reads[read];
            linked_input.bases.insert(linked_input.bases.end(), sequence.begin(), sequence.end());
            linked_input.offsets.push_back(static_cast<std::uint32_t>(linked_input.bases.size()));
        }
        linked_input.sample_ids = linked_sample_ids;
        auto linked_options = fastgatk::kernels::KmerGraphOptions{};
        linked_options.k = 5;
        linked_options.min_pruning = 1;
        linked_options.num_pruning_samples = 2;
        linked_options.max_paths = 16;
        linked_options.linked_de_bruijn_graph = true;
        linked_options.disable_seqgraph_simplification = true;
        const auto linked_graph = fastgatk::kernels::build_kmer_graph_kokkos(linked_input, linked_options);
        auto linked_disabled_options = linked_options;
        linked_disabled_options.disable_artificial_haplotype_recovery = true;
        const auto linked_disabled = fastgatk::kernels::build_kmer_graph_kokkos(
            linked_input, linked_disabled_options);
        if (linked_graph.seqgraph_linear_chain_merges != 0 ||
            linked_graph.seqgraph_diamond_merges != 0 ||
            linked_graph.seqgraph_tail_merges != 0 ||
            linked_graph.seqgraph_suffix_splits != 0 ||
            linked_graph.seqgraph_suffix_merges != 0 ||
            linked_disabled.artificial_haplotype_recovery_paths != 0 ||
            linked_graph.haplotype_path_count == 0)
            throw std::runtime_error("linked-de-Bruijn raw/recovery contract failed: seqgraph=" +
                std::to_string(linked_graph.seqgraph_nodes) + "/" +
                std::to_string(linked_graph.seqgraph_edges) + " paths=" +
                std::to_string(linked_graph.haplotype_path_count) + "/" +
                std::to_string(linked_disabled.haplotype_path_count) + " artificial=" +
                std::to_string(linked_graph.artificial_haplotype_recovery_paths));
        fastgatk::kernels::ReadErrorCorrectionInput correction_input;
        correction_input.offsets.push_back(0);
        const std::string correction_reference = "ACGTACGTACGTACGT";
        for (std::size_t read = 0; read < 5; ++read) {
            auto sequence = correction_reference;
            if (read == 4) sequence[7] = sequence[7] == 'A' ? 'C' : 'A';
            correction_input.bases.insert(correction_input.bases.end(), sequence.begin(), sequence.end());
            correction_input.qualities.insert(correction_input.qualities.end(), sequence.size(), 30);
            correction_input.offsets.push_back(static_cast<std::uint32_t>(correction_input.bases.size()));
        }
        const auto correction = fastgatk::kernels::correct_read_errors_kokkos(
            correction_input, fastgatk::kernels::ReadErrorCorrectionOptions{5, 3, 1, 30});
        if (!correction.used || correction.corrected_bases == 0 || correction.corrected_reads == 0 ||
            correction.bases.size() != correction_input.bases.size() ||
            correction.qualities.size() != correction_input.qualities.size() ||
            correction.bases[4 * correction_reference.size() + 7] != correction_reference[7])
            throw std::runtime_error("Kokkos read error correction API smoke failed");
        // GATK NearbyKmerErrorCorrector maps a sparse k-mer to a solid
        // Hamming neighbor before collecting overlapping-base consensus. A
        // length==k probe makes the two-mismatch path observable without any
        // competing one-mismatch windows.
        fastgatk::kernels::ReadErrorCorrectionInput two_mismatch_input;
        two_mismatch_input.offsets.push_back(0);
        for (std::size_t read = 0; read < 6; ++read) {
            const std::string sequence = "ACGTA";
            two_mismatch_input.bases.insert(two_mismatch_input.bases.end(), sequence.begin(), sequence.end());
            two_mismatch_input.qualities.insert(two_mismatch_input.qualities.end(), sequence.size(), 30);
            two_mismatch_input.offsets.push_back(static_cast<std::uint32_t>(two_mismatch_input.bases.size()));
        }
        const std::string two_mismatch_bad = "ATCTA";
        two_mismatch_input.bases.insert(two_mismatch_input.bases.end(),
                                        two_mismatch_bad.begin(), two_mismatch_bad.end());
        two_mismatch_input.qualities.insert(two_mismatch_input.qualities.end(), two_mismatch_bad.size(), 30);
        two_mismatch_input.offsets.push_back(static_cast<std::uint32_t>(two_mismatch_input.bases.size()));
        const auto two_mismatch = fastgatk::kernels::correct_read_errors_kokkos(
            two_mismatch_input, fastgatk::kernels::ReadErrorCorrectionOptions{5, 3, 2, 30, 1});
        const auto bad_begin = 6U * 5U;
        if (two_mismatch.corrected_kmers == 0 || two_mismatch.corrected_bases != 2 ||
            two_mismatch.corrected_reads != 1 ||
            two_mismatch.bases[bad_begin] != 'A' || two_mismatch.bases[bad_begin + 1] != 'C' ||
            two_mismatch.bases[bad_begin + 2] != 'G' || two_mismatch.bases[bad_begin + 3] != 'T' ||
            two_mismatch.bases[bad_begin + 4] != 'A')
            throw std::runtime_error("GATK two-mismatch k-mer correction semantics failed");
        // PileupReadErrorCorrector uses the Mutect2 log-likelihood ratio and
        // corrects a low-frequency plurality outlier.  A threshold just above
        // the six-observation (5 ref, 1 alt, Q30) oracle value makes the test
        // independent of any caller-level defaults.
        fastgatk::kernels::PileupReadErrorCorrectionInput pileup_input{
            {0, 1, 2, 3, 4, 5, 6}, {'A', 'A', 'A', 'A', 'A', 'C'},
            {30, 30, 30, 30, 30, 30}, {0, 0, 0, 0, 0, 0},
            {100, 100, 100, 100, 100, 100}, {}, {}};
        const auto pileup = fastgatk::kernels::correct_reads_by_pileup_kokkos(
            pileup_input, fastgatk::kernels::PileupReadErrorCorrectionOptions{3.5, 30, 15, 3});
        if (!pileup.used || pileup.loci != 1 || pileup.corrected_loci != 1 ||
            pileup.corrected_reads != 1 || pileup.corrected_bases != 1 ||
            pileup.bases.back() != 'A' || pileup.qualities.back() != 30 ||
            pileup.execution_space.empty())
            throw std::runtime_error("GATK pileup error correction API smoke failed");
        const auto pileup_disabled = fastgatk::kernels::correct_reads_by_pileup_kokkos(
            pileup_input);
        if (pileup_disabled.used || pileup_disabled.corrected_bases != 0 ||
            pileup_disabled.bases != pileup_input.bases)
            throw std::runtime_error("pileup error correction disable policy failed");
        const auto activity = fastgatk::kernels::compute_activity_profile_kokkos(
            {{0, 0, 0}, {100, 101, 500}, {10, 0, 0, 0, 3, 1, 0, 0, 10, 0, 0, 0}},
            fastgatk::kernels::ActivityProfileOptions{4, 0.002, 2, 300});
        // The raw Kokkos seed is 101; the Host-side BandPassActivityProfile
        // constructs the AssemblyRegion from its filtered segment, which
        // begins at the adjacent 100 locus.  This is the span GATK pads and
        // gives to local graph/PairHMM work.
        if (!activity.used || activity.loci != 3 || activity.active_loci != 1 ||
            activity.regions.size() != 1 || activity.regions[0].active_start != 100 ||
            activity.filter_size != 50 || activity.effective_max_prob_propagation_distance != 350 ||
            activity.execution_space.empty())
            throw std::runtime_error("Activity profile Kokkos API smoke failed");
        fastgatk::kernels::ActivityProfileInput reference_activity{
            {0, 0}, {200, 201}, {8, 2, 0, 0, 1, 9, 0, 0}, {0, 0}};
        const auto reference_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            reference_activity, fastgatk::kernels::ActivityProfileOptions{4, 0.002, 2, 300});
        // A reference-aware profile must count non-reference observations
        // against the projected reference, even when ALT is the majority;
        // the old diversity signal would report only 0.1 for the second
        // locus instead of 0.9.
        if (reference_profile.active_loci != 2 ||
            reference_profile.activity.size() != 2 ||
            std::abs(reference_profile.activity[0] - 0.2) > 1e-12 ||
            std::abs(reference_profile.activity[1] - 0.9) > 1e-12)
            throw std::runtime_error("reference-aware activity profile semantics failed");
        // --alleles is an explicit HaplotypeCaller ActivityProfileState. It
        // force-activates only the matching locus; ordinary observation
        // counts and the rest of the compact batch stay unchanged.
        fastgatk::kernels::ActivityProfileInput forced_alleles_activity{
            {0, 0}, {202, 203}, {10, 0, 0, 0, 10, 0, 0, 0}, {0, 0}};
        forced_alleles_activity.forced_allele_active = {0, 1};
        const auto forced_alleles_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            forced_alleles_activity, fastgatk::kernels::ActivityProfileOptions{4, 0.002, 2, 300});
        if (forced_alleles_profile.active_loci != 1 || forced_alleles_profile.active[0] != 0 ||
            forced_alleles_profile.active[1] != 1 ||
            forced_alleles_profile.activity[0] != 0.0 ||
            forced_alleles_profile.activity[1] != 1.0)
            throw std::runtime_error("GenotypeGivenAlleles activity-mask Kokkos semantics failed");
        // AssemblyRegionIterator deliberately wraps its pileup source in
        // IntervalAlignmentContextIterator, so an otherwise uncovered -L
        // coordinate is still a zero ActivityProfileState.  A single active
        // seed at 100 must therefore retain the complete default +/-37bp
        // Gaussian interval inside a 0..200 traversal span, rather than be
        // truncated at the compact evidence row.  Kokkos still receives just
        // that one raw pileup; only the Host profile supplies sparse zeros.
        fastgatk::kernels::ActivityProfileInput traversal_span_activity{
            {0}, {100}, {10, 0, 0, 0}, {0}, {1}};
        fastgatk::kernels::ActivityProfileOptions traversal_span_options;
        traversal_span_options.min_depth = 1;
        traversal_span_options.halo = 1;
        traversal_span_options.min_region_size = 1;
        traversal_span_options.max_region_size = 300;
        traversal_span_options.max_prob_propagation_distance = 50;
        traversal_span_options.traversal_spans = {{0, 0, 201}};
        const auto traversal_span_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            traversal_span_activity, traversal_span_options);
        const auto source_active_span = std::find_if(
            traversal_span_profile.profile_regions.begin(),
            traversal_span_profile.profile_regions.end(),
            [](const auto& region) { return region.active; });
        if (traversal_span_profile.loci != 1 || traversal_span_profile.active_loci != 1 ||
            source_active_span == traversal_span_profile.profile_regions.end() ||
            source_active_span->start != 63 || source_active_span->end != 137 ||
            traversal_span_profile.regions.size() != 1 ||
            traversal_span_profile.regions.front().start != 62 ||
            traversal_span_profile.regions.front().end != 138)
            throw std::runtime_error("GATK empty-pileup traversal ActivityProfile boundary failed");
        // A long uncovered stretch is still traversed by GATK, but it cannot
        // carry activity farther than the finite filter shoulder. The Host
        // sparse path may force-convert that zero desert; it must preserve
        // both source active segments rather than merging the two seeds or
        // clipping either tail.
        fastgatk::kernels::ActivityProfileInput separated_span_activity{
            {0, 0}, {100, 1000},
            {10, 0, 0, 0, 10, 0, 0, 0}, {0, 0}, {1, 1}};
        auto separated_span_options = traversal_span_options;
        separated_span_options.traversal_spans = {{0, 0, 1101}};
        const auto separated_span_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            separated_span_activity, separated_span_options);
        std::vector<std::pair<std::int32_t, std::int32_t>> separated_active_spans;
        for (const auto& region : separated_span_profile.profile_regions)
            if (region.active)
                separated_active_spans.emplace_back(region.start, region.end);
        if (separated_active_spans !=
                std::vector<std::pair<std::int32_t, std::int32_t>>{{63, 137}, {963, 1037}} ||
            separated_span_profile.regions.size() != 2)
            throw std::runtime_error("GATK sparse empty-pileup ActivityProfile gap contract failed");
        // No compact pileup rows at all is still a real source traversal:
        // IntervalAlignmentContextIterator supplies zero states for the
        // complete -L range, and AssemblyRegionWalker applies --force-active
        // only after the profile has been cut at its hard 300bp boundary.
        // The Host must retain those regions without manufacturing an empty
        // Kokkos input row for every genomic coordinate.
        fastgatk::kernels::ActivityProfileInput zero_span_activity;
        auto zero_span_options = traversal_span_options;
        zero_span_options.force_active = true;
        zero_span_options.traversal_spans = {{0, 1000, 1701}};
        const auto zero_span_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            zero_span_activity, zero_span_options);
        const std::vector<std::pair<std::int32_t, std::int32_t>> zero_span_boundaries{
            {1000, 1299}, {1300, 1599}, {1600, 1700}};
        if (zero_span_profile.loci != 0 || zero_span_profile.active_loci != 0 ||
            zero_span_profile.profile_regions.size() != zero_span_boundaries.size() ||
            zero_span_profile.regions.size() != zero_span_boundaries.size())
            throw std::runtime_error("GATK force-active zero-coverage traversal size failed");
        for (std::size_t index = 0; index < zero_span_boundaries.size(); ++index) {
            const auto& profile_region = zero_span_profile.profile_regions[index];
            if (!profile_region.active || profile_region.start != zero_span_boundaries[index].first ||
                profile_region.end != zero_span_boundaries[index].second ||
                profile_region.active_loci != 0 ||
                !zero_span_profile.regions[index].active)
                throw std::runtime_error("GATK force-active zero-coverage traversal boundary failed");
        }
        // A whole-contig no-coverage traversal is still split at GATK's hard
        // 300 bp AssemblyRegion boundary.  This guards the Host profile's
        // FIFO drain: erasing a prefix from a vector for each region makes
        // this exact 1 Mb source traversal quadratic, even though its Kokkos
        // evidence payload remains empty.
        auto large_zero_span_options = traversal_span_options;
        large_zero_span_options.traversal_spans = {{0, 0, 1000000}};
        const auto large_zero_span_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            fastgatk::kernels::ActivityProfileInput{}, large_zero_span_options);
        if (large_zero_span_profile.loci != 0 ||
            !large_zero_span_profile.regions.empty() ||
            large_zero_span_profile.profile_regions.size() != 3334 ||
            large_zero_span_profile.profile_regions.front().start != 0 ||
            large_zero_span_profile.profile_regions.front().end != 299 ||
            large_zero_span_profile.profile_regions.back().start != 999900 ||
            large_zero_span_profile.profile_regions.back().end != 999999 ||
            large_zero_span_profile.profile_regions.front().active ||
            large_zero_span_profile.profile_regions.back().active)
            throw std::runtime_error("GATK large zero-coverage ActivityProfile traversal failed");
        // Mutect2's quality-aware active detector uses the same flat-Beta
        // pileup likelihood-ratio formula as GATK.  Two reads at Q30 (Q35
        // after the default multi-substitution correction) exceed the
        // log10-LOD=2 seed threshold, while Q10/Q15 does not.
        fastgatk::kernels::ActivityProfileInput quality_activity{
            {0, 0}, {300, 301}, {1, 1, 0, 0, 1, 1, 0, 0}, {0, 0},
            {0, 0}, {0, 2, 4}, {0, 1, 0, 1}, {30, 30, 10, 10}};
        fastgatk::kernels::ActivityProfileOptions quality_options;
        quality_options.min_depth = 1;
        quality_options.halo = 1;
        quality_options.min_region_size = 1;
        quality_options.max_region_size = 10;
        quality_options.max_prob_propagation_distance = 10;
        quality_options.quality_aware_somatic = true;
        quality_options.initial_tumor_log10_odds = 2.0;
        const auto quality_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            quality_activity, quality_options);
        if (quality_profile.active_loci != 1 || quality_profile.active[0] != 1 ||
            quality_profile.active[1] != 0 || quality_profile.activity[0] != 1.0 ||
            quality_profile.activity[1] != 0.0)
            throw std::runtime_error("Mutect2 quality-aware activity semantics failed");
        // A source PileupQualBuffer with no alternate qualities is not a
        // zero-LOD observation: logLikelihoodRatio(nRef, empty) retains the
        // flat-Beta entropy term -log(nRef + 1).  This matters exactly at
        // --initial-tumor-lod 0, where returning 0 would activate every
        // reference-only pileup and change AssemblyRegion ownership.
        fastgatk::kernels::ActivityProfileInput reference_only_somatic;
        reference_only_somatic.tids = {0};
        reference_only_somatic.positions = {302};
        reference_only_somatic.counts = {78, 0, 0, 0};
        reference_only_somatic.reference_bases = {0};
        reference_only_somatic.indel_counts = {0};
        reference_only_somatic.quality_offsets = {0, 78};
        reference_only_somatic.quality_bases.assign(78, 0);
        reference_only_somatic.quality_values.assign(78, 40);
        reference_only_somatic.quality_alt_flags.assign(78, 0);
        reference_only_somatic.indel_quality_offsets = {0, 0};
        fastgatk::kernels::ActivityProfileOptions zero_lod_options = quality_options;
        zero_lod_options.initial_tumor_log10_odds = 0.0;
        // The source all-zero ActivityProfile is split only by its hard
        // AssemblyRegion size.  The quality threshold probe above uses a
        // deliberately tiny 10bp region cap; restore Mutect2's 300bp cap
        // here before asserting the source traversal boundaries.
        zero_lod_options.max_region_size = 300;
        zero_lod_options.traversal_spans = {{0, 0, 600}};
        const auto reference_only_somatic_profile =
            fastgatk::kernels::compute_activity_profile_kokkos(
                reference_only_somatic, zero_lod_options);
        if (reference_only_somatic_profile.active_loci != 0 ||
            reference_only_somatic_profile.active[0] != 0 ||
            reference_only_somatic_profile.activity[0] != 0.0 ||
            reference_only_somatic_profile.profile_regions.size() != 2 ||
            reference_only_somatic_profile.profile_regions[0].start != 0 ||
            reference_only_somatic_profile.profile_regions[0].end != 299 ||
            reference_only_somatic_profile.profile_regions[1].start != 300 ||
            reference_only_somatic_profile.profile_regions[1].end != 599)
            throw std::runtime_error("Mutect2 zero-LOD reference activity semantics failed");
        // One Q30 CIGAR indel among 78 pileup elements is likewise below a
        // zero initial LOD in GATK's exact PileupQualBuffer likelihood.  Keep
        // its separate indel ALT bucket on the kernel path, rather than
        // masking it on Host or declaring all indels unconditionally active.
        auto sparse_indel_somatic = reference_only_somatic;
        sparse_indel_somatic.indel_counts = {1};
        sparse_indel_somatic.indel_quality_offsets = {0, 1};
        sparse_indel_somatic.indel_quality_values = {30};
        const auto sparse_indel_somatic_profile =
            fastgatk::kernels::compute_activity_profile_kokkos(
                sparse_indel_somatic, zero_lod_options);
        if (sparse_indel_somatic_profile.active_loci != 0 ||
            sparse_indel_somatic_profile.active[0] != 0 ||
            sparse_indel_somatic_profile.activity[0] != 0.0)
            throw std::runtime_error("Mutect2 sparse indel activity semantics failed");
        // max_region_size is a hard active-span bound, independent of the
        // propagation distance.  A padded interval may overlap another
        // interval, but the Host state machine must not merge two active
        // chunks back into an unbounded graph/PairHMM window.
        fastgatk::kernels::ActivityProfileInput bounded_activity{
            {0, 0, 0, 0, 0}, {0, 50, 100, 130, 180},
            {8, 2, 0, 0, 8, 2, 0, 0, 8, 2, 0, 0, 8, 2, 0, 0, 8, 2, 0, 0},
            {0, 0, 0, 0, 0}};
        const auto bounded_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            bounded_activity, fastgatk::kernels::ActivityProfileOptions{4, 0.002, 2, 120, 200});
        // The GATK ActivityProfile hard cap chooses the rightmost eligible
        // local minimum in the first max-size span, not the sparse raw
        // evidence coordinates.  The band-pass segment therefore splits at
        // 75/76 while retaining the 120-base hard bound.
        if (bounded_profile.regions.size() != 2 ||
            bounded_profile.regions[0].active_start != 0 ||
            bounded_profile.regions[0].active_end != 75 ||
            bounded_profile.regions[1].active_start != 76 ||
            bounded_profile.regions[1].active_end != 180 ||
            bounded_profile.filter_size != 50 ||
            bounded_profile.effective_max_prob_propagation_distance != 250)
            throw std::runtime_error("AssemblyRegion hard-size/propagation controls failed");
        // Direct port of ActivityProfileUnitTest.ActiveRegionCutTests'
        // "lowest of two minima" case.  With filtering effectively disabled,
        // an active six-base hard window has minima 0.5 at base 3 and 0.75
        // at base 5; GATK selects the lower first cut (bases 1..3), rather
        // than merely the rightmost eligible local minimum.  The Kokkos
        // counts encode those exact reference-aware probabilities.
        fastgatk::kernels::ActivityProfileInput two_minima_activity{
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
            {0, 100, 0, 0,
             0, 100, 0, 0,
             50, 50, 0, 0,
             0, 100, 0, 0,
             25, 75, 0, 0,
             0, 100, 0, 0,
             0, 100, 0, 0,
             0, 100, 0, 0,
             0, 100, 0, 0,
             0, 100, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        fastgatk::kernels::ActivityProfileOptions two_minima_options;
        two_minima_options.min_depth = 1;
        two_minima_options.halo = 1;
        two_minima_options.min_region_size = 1;
        two_minima_options.max_region_size = 6;
        two_minima_options.max_prob_propagation_distance = 1;
        two_minima_options.bandpass_max_filter_size = 1;
        two_minima_options.bandpass_sigma = 0.01;
        two_minima_options.bandpass_adaptive_filter = true;
        const auto two_minima_profile = fastgatk::kernels::compute_activity_profile_kokkos(
            two_minima_activity, two_minima_options);
        if (two_minima_profile.filter_size != 0 ||
            two_minima_profile.profile_regions.empty() ||
            two_minima_profile.profile_regions.front().start != 1 ||
            two_minima_profile.profile_regions.front().end != 3 ||
            !two_minima_profile.profile_regions.front().active)
            throw std::runtime_error("GATK ActivityProfile lowest-minimum cut contract failed");
        // Additional direct cases from ActivityProfileUnitTest's
        // ActiveRegionCutTests generator.  In particular, source accepts the
        // left edge of a flat minimum (`<=` on the right, `<` on the left),
        // and it must ignore an otherwise lower minimum that would create a
        // region shorter than minRegionSize.  These checks cover choices the
        // simple two-distinct-minima probe cannot observe.
        const auto profile_from_percentages = [](const std::vector<std::uint32_t>& alt_percent,
                                                 const std::uint32_t min_region_size) {
            fastgatk::kernels::ActivityProfileInput input;
            for (std::size_t index = 0; index < alt_percent.size(); ++index) {
                input.tids.push_back(0);
                input.positions.push_back(static_cast<std::int32_t>(index + 1U));
                const auto alt = std::min<std::uint32_t>(100U, alt_percent[index]);
                input.counts.insert(input.counts.end(), {100U - alt, alt, 0U, 0U});
                input.reference_bases.push_back(0U);
            }
            fastgatk::kernels::ActivityProfileOptions options;
            options.min_depth = 1;
            options.halo = 1;
            options.min_region_size = min_region_size;
            options.max_region_size = 6;
            options.max_prob_propagation_distance = 1;
            options.bandpass_max_filter_size = 1;
            options.bandpass_sigma = 0.01;
            return fastgatk::kernels::compute_activity_profile_kokkos(input, options);
        };
        const auto flat_cut_profile = profile_from_percentages(
            {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, 1);
        const auto plateau_cut_profile = profile_from_percentages(
            {100, 50, 50, 100, 100, 100, 100, 100, 100, 100}, 1);
        const auto min_size_cut_profile = profile_from_percentages(
            {100, 50, 100, 100, 75, 100, 100, 100, 100, 100}, 5);
        const auto min_size_fallback_profile = profile_from_percentages(
            {100, 50, 100, 75, 100, 100, 100, 100, 100, 100}, 5);
        const auto first_profile_span = [](const auto& profile) {
            return profile.profile_regions.empty()
                ? std::pair<std::int32_t, std::int32_t>{-1, -1}
                : std::pair<std::int32_t, std::int32_t>{
                    profile.profile_regions.front().start, profile.profile_regions.front().end};
        };
        if (first_profile_span(flat_cut_profile) != std::pair<std::int32_t, std::int32_t>{1, 6} ||
            first_profile_span(plateau_cut_profile) != std::pair<std::int32_t, std::int32_t>{1, 2} ||
            first_profile_span(min_size_cut_profile) != std::pair<std::int32_t, std::int32_t>{1, 5} ||
            first_profile_span(min_size_fallback_profile) != std::pair<std::int32_t, std::int32_t>{1, 6})
            throw std::runtime_error("GATK ActivityProfile cut/minimum-size source matrix failed");
        const auto filters = fastgatk::kernels::filter_reads_kokkos(
            {{0, 4, 8, 12}, {0, 0, 0}, {10, 10, -1}, {60, 10, 60}, {0, 0x100, 0}},
            fastgatk::kernels::ReadFilterOptions{20, true, true, true, false, true, 0, 7});
        if (!filters.used || filters.input_reads != 3 || filters.passed_reads != 1 ||
            filters.filtered_reads != 2 || filters.keep[0] != 1 || filters.keep[1] != 0 ||
            filters.keep[2] != 0 || filters.execution_space.empty())
            throw std::runtime_error("Read filter Kokkos API smoke failed");
        fastgatk::kernels::ReadFilterOptions mapq_available_options;
        mapq_available_options.min_mapq = 0;
        mapq_available_options.exclude_unmapped = false;
        mapq_available_options.exclude_secondary = false;
        mapq_available_options.exclude_supplementary = false;
        mapq_available_options.exclude_qcfail = false;
        mapq_available_options.exclude_mapping_quality_unavailable = true;
        const auto mapq_available = fastgatk::kernels::filter_reads_kokkos(
            {{0, 1, 2}, {0, 0}, {10, 20}, {255, 60}, {0, 0}},
            mapq_available_options);
        if (mapq_available.keep != std::vector<std::uint8_t>({0, 1}) ||
            mapq_available.passed_reads != 1)
            throw std::runtime_error("MappingQualityAvailableReadFilter Kokkos API smoke failed");
        fastgatk::kernels::ReadFilterOptions mapq_nonzero_options;
        mapq_nonzero_options.min_mapq = 0;
        mapq_nonzero_options.exclude_unmapped = false;
        mapq_nonzero_options.exclude_secondary = false;
        mapq_nonzero_options.exclude_supplementary = false;
        mapq_nonzero_options.exclude_qcfail = false;
        mapq_nonzero_options.exclude_mapping_quality_zero = true;
        const auto mapq_nonzero = fastgatk::kernels::filter_reads_kokkos(
            {{0, 1, 2}, {0, 0}, {10, 20}, {0, 60}, {0, 0}},
            mapq_nonzero_options);
        if (mapq_nonzero.keep != std::vector<std::uint8_t>({0, 1}) ||
            mapq_nonzero.passed_reads != 1)
            throw std::runtime_error("MappingQualityNotZeroReadFilter Kokkos API smoke failed");
        fastgatk::kernels::ReadFilterOptions length_options;
        length_options.min_mapq = 0;
        length_options.exclude_unmapped = false;
        length_options.exclude_secondary = false;
        length_options.exclude_supplementary = false;
        length_options.exclude_qcfail = false;
        length_options.require_read_length = true;
        length_options.min_read_length = 4;
        length_options.max_read_length = 4;
        const auto length_filters = fastgatk::kernels::filter_reads_kokkos(
            {{0, 3, 7}, {0, 0}, {10, 10}, {60, 60}, {0, 0}}, length_options);
        if (length_filters.passed_reads != 1 || length_filters.filtered_reads != 1 ||
            length_filters.keep != std::vector<std::uint8_t>({0, 1}) ||
            length_filters.execution_space.empty())
            throw std::runtime_error("ReadLengthReadFilter Kokkos API smoke failed");
        fastgatk::kernels::ReadFilterInput wellformed_input{
            {0, 10, 20, 30, 40}, {0, 0, 0, 0}, {10, 20, 30, 40},
            {60, 60, 60, 60}, {0, 0, 0, 0}};
        wellformed_input.cigar_offsets = {0, 1, 4, 5, 6};
        wellformed_input.cigar_ops = {
            (10U << 4) | 0U,
            (5U << 4) | 0U, (5U << 4) | 3U, (5U << 4) | 0U,
            (10U << 4) | 0U,
            (10U << 4) | 0U};
        wellformed_input.read_group_offsets = {0, 2, 4, 4, 6};
        wellformed_input.read_groups = {'r', 'g', 'r', 'g', 'r', 'g'};
        fastgatk::kernels::ReadFilterOptions wellformed_options;
        wellformed_options.min_mapq = 0;
        wellformed_options.exclude_unmapped = false;
        wellformed_options.exclude_secondary = false;
        wellformed_options.exclude_supplementary = false;
        wellformed_options.exclude_qcfail = false;
        wellformed_options.require_good_cigar = true;
        wellformed_options.require_nonzero_reference_span = true;
        wellformed_options.require_no_n_cigar = true;
        wellformed_options.require_read_group = true;
        const auto wellformed_filters = fastgatk::kernels::filter_reads_kokkos(
            wellformed_input, wellformed_options);
        if (wellformed_filters.keep != std::vector<std::uint8_t>({1, 0, 0, 1}) ||
            wellformed_filters.passed_reads != 2)
            throw std::runtime_error("WellformedReadFilter Kokkos API smoke failed");
        fastgatk::kernels::ReadFilterInput chimeric_input{
            {0, 1, 2, 3, 4}, {0, 0, 0, 0}, {10, 20, 30, 40},
            {60, 60, 60, 60}, {0, 0, 0, 0}};
        auto append_tag = [](const std::string& value, std::vector<std::uint8_t>& payload,
                             std::vector<std::uint32_t>& offsets) {
            payload.insert(payload.end(), value.begin(), value.end());
            offsets.push_back(static_cast<std::uint32_t>(payload.size()));
        };
        chimeric_input.original_alignment_offsets = {0};
        chimeric_input.mate_contig_offsets = {0};
        append_tag("chr1,1,+,1M,60,0", chimeric_input.original_alignments,
                   chimeric_input.original_alignment_offsets);
        append_tag("chr1,1,+,1M,60,0", chimeric_input.original_alignments,
                   chimeric_input.original_alignment_offsets);
        chimeric_input.original_alignment_offsets.push_back(
            static_cast<std::uint32_t>(chimeric_input.original_alignments.size()));
        append_tag("chr2,1,+,1M,60,0", chimeric_input.original_alignments,
                   chimeric_input.original_alignment_offsets);
        append_tag("chr1", chimeric_input.mate_contigs, chimeric_input.mate_contig_offsets);
        append_tag("chr2", chimeric_input.mate_contigs, chimeric_input.mate_contig_offsets);
        chimeric_input.mate_contig_offsets.push_back(
            static_cast<std::uint32_t>(chimeric_input.mate_contigs.size()));
        chimeric_input.mate_contig_offsets.push_back(
            static_cast<std::uint32_t>(chimeric_input.mate_contigs.size()));
        fastgatk::kernels::ReadFilterOptions chimeric_options;
        chimeric_options.min_mapq = 0;
        chimeric_options.exclude_unmapped = false;
        chimeric_options.exclude_secondary = false;
        chimeric_options.exclude_supplementary = false;
        chimeric_options.exclude_qcfail = false;
        chimeric_options.require_non_chimeric_original_alignment = true;
        const auto chimeric_filters = fastgatk::kernels::filter_reads_kokkos(
            chimeric_input, chimeric_options);
        if (chimeric_filters.keep != std::vector<std::uint8_t>({1, 0, 1, 1}) ||
            chimeric_filters.passed_reads != 3)
            throw std::runtime_error("NonChimericOriginalAlignmentReadFilter Kokkos API smoke failed");
        // Presence is separate from payload length: GATK evaluates OA/XM
        // when both tags exist even if one is an empty Z value, while a
        // missing tag must pass the NonChimeric filter.  Exercise both
        // cases explicitly so flat-span encoding cannot regress this edge.
        fastgatk::kernels::ReadFilterInput empty_tag_input{
            {0, 1, 2}, {0, 0}, {10, 20}, {60, 60}, {0, 0}};
        empty_tag_input.original_alignment_offsets = {0, 0, 0};
        empty_tag_input.mate_contig_offsets = {0, 1, 1};
        empty_tag_input.mate_contigs = {'c'};
        empty_tag_input.original_alignment_present = {1, 0};
        empty_tag_input.mate_contig_present = {1, 0};
        fastgatk::kernels::ReadFilterOptions empty_tag_options;
        empty_tag_options.min_mapq = 0;
        empty_tag_options.exclude_unmapped = false;
        empty_tag_options.exclude_secondary = false;
        empty_tag_options.exclude_supplementary = false;
        empty_tag_options.exclude_qcfail = false;
        empty_tag_options.require_non_chimeric_original_alignment = true;
        const auto empty_tag_filters = fastgatk::kernels::filter_reads_kokkos(
            empty_tag_input, empty_tag_options);
        if (empty_tag_filters.keep != std::vector<std::uint8_t>({0, 1}) ||
            empty_tag_filters.passed_reads != 1)
            throw std::runtime_error("NonChimeric empty/missing tag presence semantics failed");
        const fastgatk::kernels::ReadFilterInput downsample_input{
            {0, 1, 2, 3, 4}, {0, 0, 0, 0}, {100, 100, 100, 100},
            {60, 60, 60, 60}, {0, 0, 0, 0}};
        const auto downsample_a = fastgatk::kernels::filter_reads_kokkos(
            downsample_input, fastgatk::kernels::ReadFilterOptions{0, true, true, true, false, true, 2, 99});
        const auto downsample_b = fastgatk::kernels::filter_reads_kokkos(
            downsample_input, fastgatk::kernels::ReadFilterOptions{0, true, true, true, false, true, 2, 99});
        if (downsample_a.passed_reads != 2 || downsample_a.downsampled_reads != 2 ||
            downsample_a.keep != downsample_b.keep)
            throw std::runtime_error("deterministic read downsampling failed");
        // GATK Mutect2 replaces the ordinary positional sampler with
        // MutectDownsampler.  An overloaded pool discards MAPQ <= 50 before
        // its Java reservoir is updated; seed 99 chooses the third eligible
        // read for the second reservoir slot (Random.nextInt(3) == 1).
        const fastgatk::kernels::ReadFilterInput mutect_downsample_input{
            {0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0}, {100, 100, 100, 100, 100},
            {60, 50, 52, 51, 40}, {0, 0, 0, 0, 0}};
        fastgatk::kernels::ReadFilterOptions mutect_downsample_options;
        mutect_downsample_options.min_mapq = 0;
        mutect_downsample_options.exclude_unmapped = true;
        mutect_downsample_options.exclude_secondary = true;
        mutect_downsample_options.exclude_supplementary = true;
        mutect_downsample_options.exclude_qcfail = true;
        mutect_downsample_options.max_reads_per_locus = 2;
        mutect_downsample_options.seed = 99;
        mutect_downsample_options.mutect2_downsampling = true;
        const auto mutect_downsample = fastgatk::kernels::filter_reads_kokkos(
            mutect_downsample_input, mutect_downsample_options);
        if (mutect_downsample.passed_reads != 2 ||
            mutect_downsample.downsampled_reads != 3 ||
            mutect_downsample.keep != std::vector<std::uint8_t>({1, 0, 0, 1, 0}))
            throw std::runtime_error("MutectDownsampler Kokkos API smoke failed");
        // A Mutect stride starts at the first observed alignment start, not
        // an absolute coordinate-window boundary. With cap=1 and stride=2,
        // reads at 100/101/101 form one three-read pool (capacity two),
        // while the following 102/103 pool remains intact.
        const fastgatk::kernels::ReadFilterInput mutect_stride_input{
            {0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0}, {100, 101, 101, 102, 103},
            {60, 60, 60, 60, 60}, {0, 0, 0, 0, 0}};
        auto mutect_stride_options = mutect_downsample_options;
        mutect_stride_options.max_reads_per_locus = 1;
        mutect_stride_options.mutect2_downsampling_stride = 2;
        const auto mutect_stride = fastgatk::kernels::filter_reads_kokkos(
            mutect_stride_input, mutect_stride_options);
        if (mutect_stride.passed_reads != 4 || mutect_stride.downsampled_reads != 1 ||
            mutect_stride.keep != std::vector<std::uint8_t>({1, 0, 1, 1, 1}))
            throw std::runtime_error("MutectDownsampler stride semantics failed");
        // The suspicious-read limit operates even when coverage sampling is
        // disabled: once the count reaches limit*stride, Java discards the
        // entire pending stride rather than only the low-MAPQ reads.
        const fastgatk::kernels::ReadFilterInput mutect_suspicious_input{
            {0, 1, 2, 3, 4}, {0, 0, 0, 0}, {100, 100, 101, 101},
            {60, 50, 60, 40}, {0, 0, 0, 0}};
        auto mutect_suspicious_options = mutect_downsample_options;
        mutect_suspicious_options.max_reads_per_locus = 0;
        mutect_suspicious_options.mutect2_downsampling_stride = 2;
        mutect_suspicious_options.mutect2_max_suspicious_reads_per_alignment_start = 1;
        const auto mutect_suspicious = fastgatk::kernels::filter_reads_kokkos(
            mutect_suspicious_input, mutect_suspicious_options);
        if (mutect_suspicious.passed_reads != 0 || mutect_suspicious.downsampled_reads != 4 ||
            mutect_suspicious.keep != std::vector<std::uint8_t>({0, 0, 0, 0}))
            throw std::runtime_error("MutectDownsampler suspicious-read semantics failed");
        const auto bqsr = fastgatk::kernels::count_bqsr_quality_kokkos(
            {{30, 30, 20, 30, 93}, {0, 1, 1, 0, 1}});
        if (bqsr.observations != 5 || bqsr.bins.size() != 94 ||
            bqsr.bins[30].count != 3 || bqsr.bins[30].mismatches != 1 ||
            bqsr.bins[20].count != 1 || bqsr.bins[20].mismatches != 1 ||
            bqsr.bins[93].count != 1 || bqsr.bins[93].mismatches != 1 ||
            bqsr.execution_space.empty())
            throw std::runtime_error("BQSR Kokkos quality reduction failed");
        fastgatk::kernels::BqsrCovariateObservationBatch covariate_input;
        covariate_input.covariate_ids = {0, 1, 0, 2, 1};
        covariate_input.mismatches = {0, 1, 1, 0, 1};
        const auto covariates = fastgatk::kernels::count_bqsr_covariates_kokkos(
            covariate_input);
        if (covariates.observations != 5 || covariates.covariates.size() != 3 ||
            covariates.covariates[0].count != 2 ||
            covariates.covariates[0].mismatches != 1 ||
            covariates.covariates[1].count != 2 ||
            covariates.covariates[1].mismatches != 2 ||
            covariates.covariates[2].count != 1 ||
            covariates.covariates[2].mismatches != 0 ||
            covariates.execution_space.empty() ||
            covariates.execution_policy != "TeamPolicy" ||
            !covariates.team_local_histogram || covariates.workspace_bytes == 0)
            throw std::runtime_error("BQSR Kokkos covariate reduction failed");
        fastgatk::kernels::BqsrCovariateObservationBatch sparse_covariate_input;
        sparse_covariate_input.covariate_ids = {7};
        sparse_covariate_input.mismatches = {0};
        bool sparse_covariate_rejected = false;
        try {
            (void)fastgatk::kernels::count_bqsr_covariates_kokkos(
                sparse_covariate_input);
        } catch (const std::invalid_argument&) {
            sparse_covariate_rejected = true;
        }
        if (!sparse_covariate_rejected)
            throw std::runtime_error("BQSR sparse covariate ids were not rejected");
        fastgatk::kernels::BqsrCovariateObservationBatch wide_covariate_input;
        constexpr std::size_t wide_covariate_count = 32768;
        wide_covariate_input.covariate_ids.resize(wide_covariate_count);
        wide_covariate_input.mismatches.assign(wide_covariate_count, 0);
        for (std::size_t id = 0; id < wide_covariate_count; ++id)
            wide_covariate_input.covariate_ids[id] = static_cast<std::uint32_t>(id);
        const auto wide_covariates = fastgatk::kernels::count_bqsr_covariates_kokkos(
            wide_covariate_input);
        if (wide_covariates.execution_policy != "RangePolicy" ||
            wide_covariates.team_local_histogram || wide_covariates.workspace_bytes == 0 ||
            wide_covariates.covariates.size() != wide_covariate_count ||
            wide_covariates.covariates.front().count != 1 ||
            wide_covariates.covariates.back().count != 1)
            throw std::runtime_error("BQSR bounded RangePolicy fallback failed");
        const auto transformed = fastgatk::kernels::apply_bqsr_quality_kokkos(
            {0, 30, 93}, {5, -10, 5});
        if (transformed.adjusted != std::vector<std::uint8_t>({5, 20, 93}) ||
            transformed.execution_space.empty())
            throw std::runtime_error("BQSR Kokkos quality transform failed");
        const auto genotype = fastgatk::kernels::derive_diploid_gt_gq_kokkos(
            {50, 0, 80, 100, 80, 0}, 2, 2);
        if (genotype.first_allele != std::vector<std::int32_t>({0, 1}) ||
            genotype.second_allele != std::vector<std::int32_t>({1, 1}) ||
            genotype.gq != std::vector<std::int32_t>({50, 80}) ||
            genotype.execution_space.empty() || genotype.seconds <= 0.0)
            throw std::runtime_error("Genotype PL Kokkos API smoke failed");
        const auto remapped_pl = fastgatk::kernels::remap_genotype_pl_kokkos(
            {0, 10, 20, 30, 40, 50}, 1, 3, 2, 2, {0, 2});
        if (remapped_pl.pl != std::vector<std::int32_t>({0, 30, 50}) ||
            remapped_pl.ploidy != 2 || remapped_pl.execution_space.empty() ||
            remapped_pl.seconds <= 0.0) {
            throw std::runtime_error("Genotype PL remap Kokkos API smoke failed");
        }
        const auto reordered_pl = fastgatk::kernels::remap_genotype_pl_kokkos(
            {0, 10, 20, 30, 40, 50}, 1, 3, 3, 2, {0, 2, 1});
        if (reordered_pl.pl != std::vector<std::int32_t>({0, 30, 50, 10, 40, 20}))
            throw std::runtime_error("Genotype PL reordered-ALT remap failed");
        const auto missing_projection = fastgatk::kernels::remap_genotype_pl_kokkos(
            {0, 10, 20, 30, 40, 50}, 1, 3, 3, 2, {0, -1, 2});
        if (missing_projection.pl != std::vector<std::int32_t>({
                0, std::numeric_limits<std::int32_t>::min(),
                std::numeric_limits<std::int32_t>::min(), 30,
                std::numeric_limits<std::int32_t>::min(), 50}))
            throw std::runtime_error("Genotype PL missing-allele projection failed");
        const auto remapped_ad = fastgatk::kernels::remap_allele_field_kokkos(
            {10, 20, 30, 40, 50, 60}, 2, 3, 2, {0, 2});
        if (remapped_ad.values != std::vector<std::int32_t>({10, 30, 40, 60}) ||
            remapped_ad.execution_space.empty() || remapped_ad.seconds <= 0.0)
            throw std::runtime_error("Genotype allele-field remap Kokkos API smoke failed");
        const auto missing_ad = fastgatk::kernels::remap_allele_field_kokkos(
            {10, 20, 30}, 1, 3, 3, {0, -1, 2});
        if (missing_ad.values != std::vector<std::int32_t>({
                10, std::numeric_limits<std::int32_t>::min(), 30}))
            throw std::runtime_error("Genotype allele-field missing projection failed");
        const auto allele_counts = fastgatk::kernels::count_alleles_kokkos(
            {0, 1, 1, 1}, 2, 2, 2);
        if (allele_counts.counts != std::vector<std::int32_t>({1, 3}) ||
            allele_counts.an != 4 || allele_counts.execution_space.empty() ||
            allele_counts.seconds <= 0.0)
            throw std::runtime_error("Genotype allele-count Kokkos API smoke failed");
        const auto triploid = fastgatk::kernels::derive_genotype_gt_gq_kokkos(
            {50, 0, 80, 100, 100, 80, 0, 90}, 2, 2, 3);
        if (triploid.ploidy != 3 ||
            triploid.alleles != std::vector<std::int32_t>({0, 0, 1, 0, 1, 1}) ||
            triploid.gq != std::vector<std::int32_t>({50, 80}) ||
            triploid.execution_space.empty() || triploid.seconds <= 0.0)
            throw std::runtime_error("Genotype arbitrary-ploidy Kokkos API smoke failed");
        const auto joint_pl = fastgatk::kernels::calculate_joint_genotype_pl_kokkos(
            {-10.0, -1.0, -1.0, -10.0, -1.0, -1.0}, 2, 3, 2);
        if (joint_pl.ploidy != 2 || joint_pl.allele_count != 3 ||
            joint_pl.read_count != 2 || joint_pl.pl.size() != 6 ||
            *std::min_element(joint_pl.pl.begin(), joint_pl.pl.end()) != 0 ||
            joint_pl.execution_space.empty() || joint_pl.seconds <= 0.0)
            throw std::runtime_error("Joint multi-allelic PL Kokkos API smoke failed");
        const auto genotype_priors = fastgatk::kernels::calculate_genotype_priors_kokkos(
            {0.0, -3.9542425094393248, -4.0},
            {0.0, -6.4771212547196626, -8.0}, 2);
        if (genotype_priors.ploidy != 2 || genotype_priors.allele_count != 3 ||
            genotype_priors.log10_priors.size() != 6 ||
            std::abs(genotype_priors.log10_priors[0]) > 1.0e-15 ||
            std::abs(genotype_priors.log10_priors[1] + 3.9542425094393248) > 1.0e-12 ||
            std::abs(genotype_priors.log10_priors[2] + 6.4771212547196626) > 1.0e-12 ||
            std::abs(genotype_priors.log10_priors[4] + 7.954242509439325) > 1.0e-12 ||
            genotype_priors.execution_space.empty() || genotype_priors.seconds <= 0.0)
            throw std::runtime_error("Genotype prior Kokkos API smoke failed");
        const auto posterior_assignment =
            fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                {5, 0, 5}, 1, 2, 2, {0.0, -1.0, -2.0});
        if (posterior_assignment.alleles != std::vector<std::int32_t>({0, 0}) ||
            posterior_assignment.gq.empty() || posterior_assignment.gq[0] <= 0 ||
            posterior_assignment.posterior_phred.size() != 3 ||
            posterior_assignment.prior_phred.size() != 3 ||
            posterior_assignment.posterior_phred[1] < 0.0 ||
            posterior_assignment.execution_space.empty() ||
            posterior_assignment.seconds <= 0.0)
            throw std::runtime_error("Posterior genotype-assignment Kokkos API smoke failed");
        const auto posterior = fastgatk::kernels::calculate_site_posterior_kokkos(
            {50, 0, 80, 100, 80, 0}, 2, 2, 2,
            {0.0, -3.4771212547196626, -6.477121254719663});
        if (posterior.samples_with_likelihoods != 2 || posterior.qual <= 0.0 ||
            posterior.log10_no_variant >= 0.0 || posterior.execution_space.empty() ||
            posterior.seconds <= 0.0)
            throw std::runtime_error("Genotype posterior Kokkos API smoke failed");
        const auto spanning_posterior =
            fastgatk::kernels::calculate_site_posterior_kokkos(
                {50, 0, 80, 100, 80, 0}, 2, 2, 2,
                {0.0, -3.4771212547196626, -6.477121254719663}, 1);
        if (spanning_posterior.samples_with_likelihoods != 2 ||
            spanning_posterior.log10_no_variant < posterior.log10_no_variant ||
            spanning_posterior.qual > posterior.qual)
            throw std::runtime_error("Spanning-deletion posterior no-variant API smoke failed");
        const auto cross_sample_reference =
            fastgatk::kernels::calculate_cross_sample_reference_confidence_kokkos(
                {50, 0, 80, 100, 80, 0}, 2, 2, 2,
                {0.0, -3.4771212547196626, -6.477121254719663});
        if (cross_sample_reference.samples_with_likelihoods != 2 ||
            cross_sample_reference.sample_log10_p_reference.size() != 2 ||
            cross_sample_reference.joint_qual <= 0.0 ||
            cross_sample_reference.joint_log10_p_reference >= 0.0 ||
            cross_sample_reference.execution_space.empty() ||
            cross_sample_reference.seconds <= 0.0)
            throw std::runtime_error("Cross-sample reference-confidence Kokkos API smoke failed");
        const auto cohort = fastgatk::kernels::calculate_allele_frequency_kokkos(
            {50, 0, 80, 100, 80, 0}, 2, 2, 2, {10.0, 0.01});
        const auto cohort_count_sum = cohort.effective_allele_counts[0] +
                                      cohort.effective_allele_counts[1];
        if (cohort.integer_allele_counts.size() != 2 ||
            cohort.log10_p_allele_absent.size() != 2 ||
            cohort.log10_allele_frequencies.size() != 2 ||
            cohort.samples_with_likelihoods != 2 || cohort.iterations <= 0 ||
            !cohort.converged || cohort.qual <= 0.0 || cohort.log10_p_no_variant >= 0.0 ||
            cohort_count_sum < 3.9 || cohort_count_sum > 4.1 ||
            cohort.execution_space.empty() || cohort.seconds <= 0.0)
            throw std::runtime_error("Cohort allele-frequency Kokkos API smoke failed");
        const auto gq_only_cohort = fastgatk::kernels::calculate_allele_frequency_kokkos(
            {}, 1, 2, 2, {10.0, 0.01}, {30}, {0, 0});
        if (gq_only_cohort.approximate_gq_samples != 1 ||
            gq_only_cohort.samples_with_likelihoods != 1 ||
            gq_only_cohort.effective_allele_counts.size() != 2 ||
            gq_only_cohort.effective_allele_counts[0] <= 1.9 ||
            gq_only_cohort.effective_allele_counts[1] >= 0.1 ||
            gq_only_cohort.execution_space.empty() ||
            gq_only_cohort.seconds <= 0.0)
            throw std::runtime_error("GQ-only hom-ref AF approximation smoke failed");
        // GATK treats REF/REF and genotypes containing a spanning-deletion
        // allele (*) as the non-variant set when computing P(no-variant).
        // The explicit index keeps that rule in the shared Kokkos AF kernel
        // instead of collapsing '*' into an ordinary alternate allele.
        const auto without_spanning_deletion =
            fastgatk::kernels::calculate_allele_frequency_kokkos(
                {0, 0, 100, 100, 100, 100}, 1, 3, 2, {10.0, 0.01, 0.01}, {}, {}, -1);
        const auto with_spanning_deletion =
            fastgatk::kernels::calculate_allele_frequency_kokkos(
                {0, 0, 100, 100, 100, 100}, 1, 3, 2, {10.0, 0.01, 0.01}, {}, {}, 1);
        if (!(with_spanning_deletion.qual < without_spanning_deletion.qual) ||
            !(with_spanning_deletion.log10_p_no_variant >
              without_spanning_deletion.log10_p_no_variant) ||
            with_spanning_deletion.samples_with_likelihoods != 1)
            throw std::runtime_error("Spanning-deletion AF reduction smoke failed");
        const auto reference_confidence =
            fastgatk::kernels::calculate_reference_confidence_kokkos(
                {0, 0, 1}, {30, 30, 10}, {0, 1, 0}, 2);
        if (reference_confidence.depth != std::vector<std::uint32_t>({2, 1}) ||
            reference_confidence.reference_count != std::vector<std::uint32_t>({1, 1}) ||
            reference_confidence.non_ref_count != std::vector<std::uint32_t>({1, 0}) ||
            reference_confidence.hom_ref.size() != 2 ||
            reference_confidence.het.size() != 2 ||
            reference_confidence.hom_alt.size() != 2 ||
            !std::isfinite(reference_confidence.het[0]) ||
            reference_confidence.execution_space.empty() ||
            reference_confidence.seconds <= 0.0)
            throw std::runtime_error("Reference-confidence Kokkos API smoke failed");
        const auto reference_confidence_polyploid =
            fastgatk::kernels::calculate_reference_confidence_genotypes_kokkos(
                {0, 0, 1},
                {std::log10(0.9), std::log10(1.0e-3), std::log10(0.9)},
                {std::log10(0.1), std::log10(0.9), std::log10(0.1)},
                2, 3);
        if (reference_confidence_polyploid.ploidy != 3 ||
            reference_confidence_polyploid.locus_count != 2 ||
            reference_confidence_polyploid.genotype_pl.size() != 8 ||
            reference_confidence_polyploid.gq.size() != 2 ||
            reference_confidence_polyploid.genotype_pl[0] != 0 ||
            reference_confidence_polyploid.genotype_pl[4] != 0 ||
            reference_confidence_polyploid.genotype_pl[5] != 2 ||
            reference_confidence_polyploid.genotype_pl[6] != 4 ||
            reference_confidence_polyploid.genotype_pl[7] != 10 ||
            reference_confidence_polyploid.gq[1] != 2 ||
            reference_confidence_polyploid.execution_space.empty() ||
            reference_confidence_polyploid.prepare_seconds <= 0.0 ||
            reference_confidence_polyploid.seconds <= 0.0)
            throw std::runtime_error(
                "Polyploid reference-confidence genotype Kokkos API smoke failed");
        std::cout << "{\"status\":\"pass\",\"pairhmm_likelihood\":"
                  << result.likelihoods[0] << ",\"sw_score\":" << score
                  << ",\"sw_kokkos_score\":" << sw_batch.scores[0]
                  << ",\"sw_cigar\":\"" << alignment.cigar
                  << "\",\"pairhmm_matrix_pairs\":" << full_matrix.likelihoods.size() << "}\n";
        persistent_plan.clear();
        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
