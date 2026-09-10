#include "fastgatk/kernels/pairhmm_kokkos.hpp"
#include "fastgatk/kernels/smith_waterman.hpp"
#include "fastgatk/kernels/kmer_graph.hpp"
#include "fastgatk/kernels/activity_profile.hpp"
#include "fastgatk/kernels/bqsr.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "fastgatk/kernels/reference_confidence.hpp"
#include "fastgatk/kernels/read_error_correction.hpp"
#include "fastgatk/kernels/somatic.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

std::size_t positive_env(const char* name, std::size_t fallback, std::size_t maximum) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    try {
        std::size_t consumed = 0;
        const auto value = std::stoull(raw, &consumed);
        if (consumed != std::string(raw).size() || value == 0 || value > maximum)
            throw std::invalid_argument("range");
        return static_cast<std::size_t>(value);
    } catch (...) {
        throw std::invalid_argument(std::string(name) + " must be an integer in [1," +
                                    std::to_string(maximum) + "]");
    }
}

struct TimingSummary {
    std::vector<double> samples;
    double p50 = 0.0;
    double p95 = 0.0;
};

TimingSummary summarize(std::vector<double> samples) {
    if (samples.empty()) return {};
    std::sort(samples.begin(), samples.end());
    const auto quantile = [&](double q) {
        const double position = q * static_cast<double>(samples.size() - 1);
        const auto lower = static_cast<std::size_t>(position);
        const auto upper = std::min(samples.size() - 1, lower + 1);
        const double fraction = position - static_cast<double>(lower);
        return samples[lower] + fraction * (samples[upper] - samples[lower]);
    };
    const double p50 = quantile(0.50);
    const double p95 = quantile(0.95);
    return {std::move(samples), p50, p95};
}

void emit_timing(const char* name, const TimingSummary& summary) {
    std::cout << "\"" << name << "\":{\"p50_seconds\":" << summary.p50
              << ",\"p95_seconds\":" << summary.p95 << ",\"samples\":[";
    for (std::size_t index = 0; index < summary.samples.size(); ++index) {
        if (index != 0) std::cout << ',';
        std::cout << summary.samples[index];
    }
    std::cout << "]}";
}

}  // namespace

int main() {
    Kokkos::initialize();
    try {
        const auto repeats = positive_env("FASTGATK_BENCH_REPEATS", 5, 1000);
        const auto warmup = positive_env("FASTGATK_BENCH_WARMUP", 1, 100);
        constexpr std::size_t requests = 512;
        constexpr std::size_t somatic_candidates = 128;
        constexpr std::size_t somatic_reads = 64;
        constexpr std::size_t somatic_multiallelic_allele_count = 4;
        constexpr int joint_allele_count = 3;
        constexpr int joint_ploidy = 4;
        constexpr std::size_t bqsr_observations = 1U << 16;
        std::vector<fastgatk::kernels::SmithWatermanRequest> sw_requests;
        std::vector<fastgatk::kernels::SmithWatermanRequest> sw_uniform_requests;
        sw_requests.reserve(requests);
        sw_uniform_requests.reserve(requests);
        std::vector<fastgatk::pairhmm::PairHmmRead> reads;
        std::vector<fastgatk::pairhmm::PairHmmHaplotype> haplotypes;
        std::vector<fastgatk::pairhmm::PairHmmRequest> pairs;
        std::vector<double> normalization_values;
        std::vector<std::uint32_t> normalization_read_ids;
        std::vector<std::uint8_t> normalization_eligible;
        std::vector<fastgatk::pairhmm::AlleleMarginalizationRequest>
            marginalization_requests;
        std::vector<fastgatk::pairhmm::ReadAlleleUncertaintyRequest>
            uncertainty_requests;
        std::vector<fastgatk::pairhmm::ReadAlleleBestRequest>
            best_allele_requests;
        std::vector<fastgatk::pairhmm::FlowPairHmmRead> flow_reads;
        std::vector<fastgatk::pairhmm::FlowPairHmmHaplotype> flow_haplotypes;
        std::vector<fastgatk::pairhmm::FlowPairHmmRequest> flow_pairs;
        std::vector<std::int32_t> genotype_pl;
        std::vector<std::int32_t> genotype_alleles;
        std::vector<double> joint_likelihoods(
            static_cast<std::size_t>(joint_allele_count) * requests, -0.0);
        const std::vector<double> genotype_posterior_priors{
            0.0, -1.0, -2.0};
        const std::vector<double> genotype_prior_het{
            0.0, -3.9542425094393248, -4.0};
        const std::vector<double> genotype_prior_hom{
            0.0, -6.4771212547196626, -8.0};
        std::vector<std::uint32_t> rcm_locus;
        std::vector<std::uint8_t> rcm_quality;
        std::vector<std::uint8_t> rcm_alt;
        std::vector<double> rcm_poly_reference;
        std::vector<double> rcm_poly_non_ref;
        std::vector<double> somatic_reference_likelihoods(
            somatic_candidates * somatic_reads, -0.0);
        std::vector<double> somatic_alternate_likelihoods(
            somatic_candidates * somatic_reads, -0.0);
        std::vector<double> somatic_multiallelic_likelihoods(
            somatic_multiallelic_allele_count * somatic_reads, -0.0);
        std::vector<std::uint32_t> somatic_f1r2(somatic_candidates, 2);
        std::vector<std::uint32_t> somatic_r1f2(somatic_candidates, 2);
        fastgatk::kernels::BqsrObservationBatch bqsr_input;
        bqsr_input.qualities.resize(bqsr_observations);
        bqsr_input.mismatches.resize(bqsr_observations);
        fastgatk::kernels::BqsrCovariateObservationBatch bqsr_covariate_input;
        bqsr_covariate_input.covariate_ids.resize(bqsr_observations);
        bqsr_covariate_input.mismatches.resize(bqsr_observations);
        std::vector<std::int16_t> bqsr_deltas(bqsr_observations);
        for (std::size_t index = 0; index < bqsr_observations; ++index) {
            bqsr_input.qualities[index] = static_cast<std::uint8_t>(20 + (index % 60));
            bqsr_input.mismatches[index] = (index % 17 == 0) ? 1U : 0U;
            bqsr_covariate_input.covariate_ids[index] =
                static_cast<std::uint32_t>(index % 512);
            bqsr_covariate_input.mismatches[index] = bqsr_input.mismatches[index];
            bqsr_deltas[index] = static_cast<std::int16_t>(static_cast<int>(index % 7) - 3);
        }
        reads.reserve(requests);
        haplotypes.reserve(requests);
        pairs.reserve(requests);
        normalization_values.reserve(requests);
        normalization_read_ids.reserve(requests);
        normalization_eligible.reserve(requests);
        marginalization_requests.reserve(requests);
        uncertainty_requests.reserve(requests);
        best_allele_requests.reserve(requests);
        flow_reads.reserve(requests);
        flow_haplotypes.reserve(requests);
        flow_pairs.reserve(requests);
        genotype_pl.reserve(requests * 3);
        genotype_alleles.reserve(requests * 2);
        rcm_locus.reserve(requests);
        rcm_quality.reserve(requests);
        rcm_alt.reserve(requests);
        rcm_poly_reference.reserve(requests);
        rcm_poly_non_ref.reserve(requests);
        for (std::size_t i = 0; i < requests; ++i) {
            const auto read_length = 64 + (i % 4);
            const auto haplotype_length = 80 + (i % 5);
            fastgatk::kernels::SmithWatermanRequest sw;
            sw.read.resize(read_length);
            sw.reference.resize(haplotype_length);
            fastgatk::pairhmm::PairHmmRead pair_read;
            pair_read.bases.resize(read_length);
            pair_read.qualities.assign(read_length, 30);
            pair_read.insertion_gop.assign(read_length, 40);
            pair_read.deletion_gop.assign(read_length, 40);
            pair_read.gap_continuation.assign(read_length, 10);
            fastgatk::pairhmm::PairHmmHaplotype pair_haplotype;
            pair_haplotype.bases.resize(haplotype_length);
            for (std::size_t j = 0; j < read_length; ++j)
                sw.read[j] = pair_read.bases[j] = static_cast<std::uint8_t>((j + i) & 3U);
            for (std::size_t j = 0; j < haplotype_length; ++j)
                sw.reference[j] = pair_haplotype.bases[j] = static_cast<std::uint8_t>((j + i) & 3U);
            sw_requests.push_back(std::move(sw));
            fastgatk::kernels::SmithWatermanRequest uniform;
            uniform.read.resize(72);
            uniform.reference.resize(96);
            for (std::size_t j = 0; j < uniform.read.size(); ++j)
                uniform.read[j] = static_cast<std::uint8_t>((j + i * 3) & 3U);
            for (std::size_t j = 0; j < uniform.reference.size(); ++j)
                uniform.reference[j] = static_cast<std::uint8_t>((j + i * 5) & 3U);
            sw_uniform_requests.push_back(std::move(uniform));
            reads.push_back(std::move(pair_read));
            haplotypes.push_back(std::move(pair_haplotype));
            pairs.push_back({static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i)});
            normalization_values.push_back(-0.01 * static_cast<double>(i % 17));
            normalization_read_ids.push_back(static_cast<std::uint32_t>(i % 64));
            normalization_eligible.push_back((i % 7 == 0) ? 0U : 1U);
            marginalization_requests.push_back(
                fastgatk::pairhmm::AlleleMarginalizationRequest{
                    static_cast<std::uint32_t>(i / 2),
                    static_cast<std::uint8_t>(i & 1U),
                    normalization_values.back()});
            uncertainty_requests.push_back(
                fastgatk::pairhmm::ReadAlleleUncertaintyRequest{
                    static_cast<std::uint32_t>(i / 2), normalization_values.back()});
            best_allele_requests.push_back(
                fastgatk::pairhmm::ReadAlleleBestRequest{
                    static_cast<std::uint32_t>(i / 2),
                    static_cast<std::uint32_t>(i & 1U), normalization_values.back()});
            fastgatk::pairhmm::FlowPairHmmRead flow_read;
            const auto flow_read_length = 32 + (i % 2);
            flow_read.key.resize(flow_read_length);
            flow_read.flow_order.resize(flow_read_length);
            flow_read.insertion_gop.assign(flow_read_length, 40);
            flow_read.deletion_gop.assign(flow_read_length, 40);
            flow_read.gap_continuation.assign(flow_read_length, 10);
            flow_read.probabilities.assign(flow_read_length * 256, 0.001);
            for (std::size_t j = 0; j < flow_read_length; ++j) {
                flow_read.key[j] = 1 + static_cast<std::int32_t>((j + i) % 3);
                flow_read.flow_order[j] = static_cast<std::uint8_t>((j + i) & 3U);
                flow_read.probabilities[j * 256 + static_cast<std::size_t>(flow_read.key[j])] = 0.999;
            }
            fastgatk::pairhmm::FlowPairHmmHaplotype flow_haplotype;
            const auto flow_haplotype_length = 40 + (i % 3);
            flow_haplotype.key.resize(flow_haplotype_length);
            flow_haplotype.flow_order.resize(flow_haplotype_length);
            for (std::size_t j = 0; j < flow_haplotype_length; ++j) {
                flow_haplotype.key[j] = 1 + static_cast<std::int32_t>((j + i) % 3);
                flow_haplotype.flow_order[j] = static_cast<std::uint8_t>((j + i) & 3U);
            }
            flow_reads.push_back(std::move(flow_read));
            flow_haplotypes.push_back(std::move(flow_haplotype));
            flow_pairs.push_back({static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i)});
            if ((i & 1U) == 0)
                genotype_pl.insert(genotype_pl.end(), {50, 0, 80});
            else
                genotype_pl.insert(genotype_pl.end(), {100, 80, 0});
            if ((i & 1U) == 0)
                genotype_alleles.insert(genotype_alleles.end(), {0, 1});
            else
                genotype_alleles.insert(genotype_alleles.end(), {1, 1});
            for (int allele = 0; allele < joint_allele_count; ++allele) {
                joint_likelihoods[static_cast<std::size_t>(allele) * requests + i] =
                    -0.01 * static_cast<double>((i + static_cast<std::size_t>(allele) * 3) % 17);
            }
            rcm_locus.push_back(static_cast<std::uint32_t>(i % 128));
            const auto quality = static_cast<std::uint8_t>(30 + (i % 5));
            const auto alt = (i % 7 == 0) ? 1U : 0U;
            rcm_quality.push_back(quality);
            rcm_alt.push_back(alt);
            const auto error = std::pow(10.0, -0.1 * static_cast<double>(quality));
            const auto correct = std::max(1.0 - error, 1e-300);
            const auto incorrect = std::max(error / 3.0, 1e-300);
            rcm_poly_reference.push_back(std::log10(alt != 0 ? incorrect : correct));
            rcm_poly_non_ref.push_back(std::log10(alt != 0 ? correct : incorrect));
        }
        for (std::size_t candidate = 0; candidate < somatic_candidates; ++candidate) {
            for (std::size_t read = 0; read < somatic_reads; ++read) {
                const auto index = candidate * somatic_reads + read;
                somatic_reference_likelihoods[index] = -0.05 * static_cast<double>(read % 3);
                somatic_alternate_likelihoods[index] =
                    -0.01 * static_cast<double>((read + candidate) % 5);
                if ((read + candidate) % 17 == 0)
                    somatic_reference_likelihoods[index] = -1.0e300;
            }
        }
        for (std::size_t read = 0; read < somatic_reads; ++read) {
            somatic_multiallelic_likelihoods[read] =
                -0.05 * static_cast<double>(read % 3);
            somatic_multiallelic_likelihoods[somatic_reads + read] =
                -0.01 * static_cast<double>(read % 5);
            somatic_multiallelic_likelihoods[2 * somatic_reads + read] =
                -0.02 * static_cast<double>((read + 2) % 5);
            somatic_multiallelic_likelihoods[3 * somatic_reads + read] =
                -0.03 * static_cast<double>((read + 3) % 5);
        }
        // Warmup is deliberately outside reported samples.  It exercises
        // Kokkos initialization and first allocation without hiding that cost
        // in steady-state kernel timings.
        for (std::size_t iteration = 0; iteration < warmup; ++iteration) {
            (void)fastgatk::kernels::smith_waterman_score_kokkos(sw_requests);
            (void)fastgatk::kernels::smith_waterman_score_kokkos(sw_uniform_requests);
            (void)fastgatk::pairhmm::compute_kokkos_bucketed(reads, haplotypes, pairs, 1);
            (void)fastgatk::pairhmm::compute_kokkos_bucketed(
                reads, haplotypes, pairs, 1, fastgatk::pairhmm::PairHmmPrecision::Float32);
            (void)fastgatk::pairhmm::normalize_likelihoods_kokkos(
                normalization_values, normalization_read_ids, normalization_eligible, 64, 4.5);
            (void)fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos(
                marginalization_requests, requests / 2);
            (void)fastgatk::pairhmm::reduce_read_allele_uncertainty_kokkos(
                uncertainty_requests, requests / 2);
            (void)fastgatk::pairhmm::compute_kokkos_flow(flow_reads, flow_haplotypes, flow_pairs, 1);
            (void)fastgatk::kernels::derive_diploid_gt_gq_kokkos(genotype_pl, requests, 2);
            (void)fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                genotype_pl, requests, 2, 2, genotype_posterior_priors);
            (void)fastgatk::kernels::calculate_joint_genotype_pl_kokkos(
                joint_likelihoods, requests, joint_allele_count, joint_ploidy);
            (void)fastgatk::kernels::calculate_genotype_priors_kokkos(
                genotype_prior_het, genotype_prior_hom, joint_ploidy);
            (void)fastgatk::kernels::count_alleles_kokkos(genotype_alleles, requests, 2, 2);
            (void)fastgatk::kernels::calculate_site_posterior_kokkos(
                genotype_pl, requests, 2, 2,
                {0.0, -3.4771212547196626, -6.477121254719663});
            (void)fastgatk::kernels::calculate_cross_sample_reference_confidence_kokkos(
                genotype_pl, requests, 2, 2,
                {0.0, -3.4771212547196626, -6.477121254719663});
            (void)fastgatk::kernels::calculate_allele_frequency_kokkos(
                genotype_pl, requests, 2, 2, {10.0, 0.01});
            (void)fastgatk::kernels::calculate_reference_confidence_kokkos(
                rcm_locus, rcm_quality, rcm_alt, 128);
            (void)fastgatk::kernels::calculate_reference_confidence_genotypes_kokkos(
                rcm_locus, rcm_poly_reference, rcm_poly_non_ref, 128, joint_ploidy);
            (void)fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                somatic_reference_likelihoods, somatic_alternate_likelihoods,
                somatic_candidates, somatic_reads, 101);
            (void)fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                somatic_multiallelic_likelihoods,
                somatic_multiallelic_allele_count, somatic_reads);
            (void)fastgatk::kernels::calculate_somatic_posterior_kokkos(
                somatic_reference_likelihoods, somatic_alternate_likelihoods,
                {}, {}, somatic_f1r2, somatic_r1f2,
                somatic_candidates, somatic_reads, 0.05);
            (void)fastgatk::kernels::count_bqsr_quality_kokkos(bqsr_input);
            (void)fastgatk::kernels::count_bqsr_covariates_kokkos(bqsr_covariate_input);
            (void)fastgatk::kernels::apply_bqsr_quality_kokkos(
                bqsr_input.qualities, bqsr_deltas);
        }

        fastgatk::pairhmm::PersistentBucketPlan persistent_pairhmm(32);
        std::vector<double> sw_prepare_samples;
        std::vector<double> sw_kernel_samples;
        std::vector<double> sw_uniform_prepare_samples;
        std::vector<double> sw_uniform_kernel_samples;
        std::vector<double> pairhmm_prepare_samples;
        std::vector<double> pairhmm_kernel_samples;
        std::vector<double> float_pairhmm_prepare_samples;
        std::vector<double> float_pairhmm_kernel_samples;
        std::vector<double> normalization_prepare_samples;
        std::vector<double> normalization_kernel_samples;
        std::vector<double> marginalization_prepare_samples;
        std::vector<double> marginalization_kernel_samples;
        std::vector<double> uncertainty_prepare_samples;
        std::vector<double> uncertainty_kernel_samples;
        std::vector<double> best_allele_prepare_samples;
        std::vector<double> best_allele_kernel_samples;
        std::vector<double> flow_pairhmm_prepare_samples;
        std::vector<double> flow_pairhmm_kernel_samples;
        std::vector<double> persistent_prepare_samples;
        std::vector<double> persistent_kernel_samples;
        std::vector<double> genotype_prepare_samples;
        std::vector<double> genotype_kernel_samples;
        std::vector<double> joint_genotype_prepare_samples;
        std::vector<double> joint_genotype_kernel_samples;
        std::vector<double> genotype_prior_prepare_samples;
        std::vector<double> genotype_prior_kernel_samples;
        std::vector<double> posterior_assignment_prepare_samples;
        std::vector<double> posterior_assignment_kernel_samples;
        std::vector<double> allele_count_prepare_samples;
        std::vector<double> allele_count_kernel_samples;
        std::vector<double> posterior_prepare_samples;
        std::vector<double> posterior_kernel_samples;
        std::vector<double> cross_sample_reference_prepare_samples;
        std::vector<double> cross_sample_reference_kernel_samples;
        std::vector<double> cohort_prepare_samples;
        std::vector<double> cohort_kernel_samples;
        std::vector<double> rcm_prepare_samples;
        std::vector<double> rcm_kernel_samples;
        std::vector<double> rcm_poly_prepare_samples;
        std::vector<double> rcm_poly_kernel_samples;
        std::vector<double> somatic_prepare_samples;
        std::vector<double> somatic_kernel_samples;
        std::vector<double> somatic_posterior_prepare_samples;
        std::vector<double> somatic_posterior_kernel_samples;
        std::vector<double> somatic_multiallelic_prepare_samples;
        std::vector<double> somatic_multiallelic_kernel_samples;
        std::vector<double> bqsr_count_prepare_samples;
        std::vector<double> bqsr_count_kernel_samples;
        std::vector<double> bqsr_covariate_prepare_samples;
        std::vector<double> bqsr_covariate_kernel_samples;
        std::vector<double> bqsr_apply_prepare_samples;
        std::vector<double> bqsr_apply_kernel_samples;
        sw_prepare_samples.reserve(repeats);
        sw_kernel_samples.reserve(repeats);
        pairhmm_prepare_samples.reserve(repeats);
        pairhmm_kernel_samples.reserve(repeats);
        float_pairhmm_prepare_samples.reserve(repeats);
        float_pairhmm_kernel_samples.reserve(repeats);
        normalization_prepare_samples.reserve(repeats);
        normalization_kernel_samples.reserve(repeats);
        marginalization_prepare_samples.reserve(repeats);
        marginalization_kernel_samples.reserve(repeats);
        uncertainty_prepare_samples.reserve(repeats);
        uncertainty_kernel_samples.reserve(repeats);
        best_allele_prepare_samples.reserve(repeats);
        best_allele_kernel_samples.reserve(repeats);
        flow_pairhmm_prepare_samples.reserve(repeats);
        flow_pairhmm_kernel_samples.reserve(repeats);
        persistent_prepare_samples.reserve(repeats);
        persistent_kernel_samples.reserve(repeats);
        genotype_prepare_samples.reserve(repeats);
        genotype_kernel_samples.reserve(repeats);
        joint_genotype_prepare_samples.reserve(repeats);
        joint_genotype_kernel_samples.reserve(repeats);
        genotype_prior_prepare_samples.reserve(repeats);
        genotype_prior_kernel_samples.reserve(repeats);
        posterior_assignment_prepare_samples.reserve(repeats);
        posterior_assignment_kernel_samples.reserve(repeats);
        allele_count_prepare_samples.reserve(repeats);
        allele_count_kernel_samples.reserve(repeats);
        posterior_prepare_samples.reserve(repeats);
        posterior_kernel_samples.reserve(repeats);
        cross_sample_reference_prepare_samples.reserve(repeats);
        cross_sample_reference_kernel_samples.reserve(repeats);
        cohort_prepare_samples.reserve(repeats);
        cohort_kernel_samples.reserve(repeats);
        rcm_prepare_samples.reserve(repeats);
        rcm_kernel_samples.reserve(repeats);
        rcm_poly_prepare_samples.reserve(repeats);
        rcm_poly_kernel_samples.reserve(repeats);
        somatic_prepare_samples.reserve(repeats);
        somatic_kernel_samples.reserve(repeats);
        somatic_posterior_prepare_samples.reserve(repeats);
        somatic_posterior_kernel_samples.reserve(repeats);
        somatic_multiallelic_prepare_samples.reserve(repeats);
        somatic_multiallelic_kernel_samples.reserve(repeats);
        bqsr_count_prepare_samples.reserve(repeats);
        bqsr_count_kernel_samples.reserve(repeats);
        bqsr_covariate_prepare_samples.reserve(repeats);
        bqsr_covariate_kernel_samples.reserve(repeats);
        bqsr_apply_prepare_samples.reserve(repeats);
        bqsr_apply_kernel_samples.reserve(repeats);

        std::uint64_t sw_checksum = 0;
        std::uint64_t sw_uniform_checksum = 0;
        double pairhmm_checksum = 0.0;
        double float_pairhmm_checksum = 0.0;
        double float_pairhmm_max_abs_error = 0.0;
        std::uint64_t normalization_checksum = 0;
        std::uint64_t marginalization_checksum = 0;
        std::uint64_t uncertainty_checksum = 0;
        std::uint64_t best_allele_checksum = 0;
        double flow_pairhmm_checksum = 0.0;
        double persistent_checksum = 0.0;
        std::uint64_t genotype_checksum = 0;
        std::uint64_t joint_genotype_checksum = 0;
        std::uint64_t genotype_prior_checksum = 0;
        std::uint64_t posterior_assignment_checksum = 0;
        std::uint64_t allele_count_checksum = 0;
        std::uint64_t posterior_checksum = 0;
        std::uint64_t cross_sample_reference_checksum = 0;
        std::uint64_t cohort_checksum = 0;
        std::uint64_t rcm_checksum = 0;
        std::uint64_t rcm_poly_checksum = 0;
        std::uint64_t somatic_checksum = 0;
        std::uint64_t somatic_posterior_checksum = 0;
        std::uint64_t somatic_multiallelic_checksum = 0;
        std::uint64_t bqsr_count_checksum = 0;
        std::uint64_t bqsr_covariate_checksum = 0;
        std::uint64_t bqsr_apply_checksum = 0;
        std::size_t persistent_cached_shapes = 0;
        std::size_t persistent_cache_hits = 0;
        std::string execution_space;
        std::string normalization_execution_space;
        std::string marginalization_execution_space;
        std::string uncertainty_execution_space;
        std::string best_allele_execution_space;
        std::string flow_pairhmm_execution_policy;
        std::size_t simd_width = 1;
        std::size_t float_pairhmm_simd_width = 1;
        std::size_t sw_simd_width = 1;
        std::size_t sw_simd_groups = 0;
        std::string genotype_execution_space;
        std::string joint_genotype_execution_space;
        std::string genotype_prior_execution_space;
        std::string posterior_assignment_execution_space;
        std::string allele_count_execution_space;
        std::string posterior_execution_space;
        std::string cross_sample_reference_execution_space;
        std::string cohort_execution_space;
        std::string rcm_execution_space;
        std::string rcm_poly_execution_space;
        std::string somatic_execution_space;
        std::string somatic_posterior_execution_space;
        std::string somatic_multiallelic_execution_space;
        std::string bqsr_count_execution_space;
        std::string bqsr_count_execution_policy;
        std::string bqsr_covariate_execution_space;
        std::string bqsr_covariate_execution_policy;
        std::uint64_t bqsr_covariate_workspace_bytes = 0;
        bool bqsr_covariate_team_local_histogram = false;
        std::string bqsr_apply_execution_space;
        std::string bqsr_apply_execution_policy;
        for (std::size_t iteration = 0; iteration < repeats; ++iteration) {
            const auto sw = fastgatk::kernels::smith_waterman_score_kokkos(sw_requests);
            const auto sw_uniform = fastgatk::kernels::smith_waterman_score_kokkos(sw_uniform_requests);
            const auto pairhmm = fastgatk::pairhmm::compute_kokkos_bucketed(
                reads, haplotypes, pairs, 1);
            const auto float_pairhmm = fastgatk::pairhmm::compute_kokkos_bucketed(
                reads, haplotypes, pairs, 1, fastgatk::pairhmm::PairHmmPrecision::Float32);
            if (float_pairhmm.precision != "float32" ||
                float_pairhmm.likelihoods.size() != pairhmm.likelihoods.size())
                throw std::runtime_error("float PairHMM precision contract failed");
            auto normalized = fastgatk::pairhmm::normalize_likelihoods_kokkos(
                pairhmm.likelihoods, normalization_read_ids, normalization_eligible, 64, 4.5);
            for (std::size_t index = 0; index < marginalization_requests.size(); ++index)
                marginalization_requests[index].likelihood = pairhmm.likelihoods[index];
            for (std::size_t index = 0; index < uncertainty_requests.size(); ++index)
                uncertainty_requests[index].likelihood = pairhmm.likelihoods[index];
            for (std::size_t index = 0; index < best_allele_requests.size(); ++index)
                best_allele_requests[index].likelihood = pairhmm.likelihoods[index];
            const auto marginalized =
                fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos(
                    marginalization_requests, requests / 2);
            const auto uncertainty =
                fastgatk::pairhmm::reduce_read_allele_uncertainty_kokkos(
                    uncertainty_requests, requests / 2);
            const auto best_alleles =
                fastgatk::pairhmm::reduce_read_allele_best_kokkos(
                    best_allele_requests, requests / 2);
            const auto flow_pairhmm = fastgatk::pairhmm::compute_kokkos_flow(
                flow_reads, flow_haplotypes, flow_pairs, 1);
            const auto persistent = persistent_pairhmm.execute(reads, haplotypes, pairs, 1);
            const auto genotype = fastgatk::kernels::derive_diploid_gt_gq_kokkos(
                genotype_pl, requests, 2);
            const auto joint_genotype = fastgatk::kernels::calculate_joint_genotype_pl_kokkos(
                joint_likelihoods, requests, joint_allele_count, joint_ploidy);
            const auto genotype_prior = fastgatk::kernels::calculate_genotype_priors_kokkos(
                genotype_prior_het, genotype_prior_hom, joint_ploidy);
            const auto posterior_assignment =
                fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                    genotype_pl, requests, 2, 2, genotype_posterior_priors);
            const auto allele_counts = fastgatk::kernels::count_alleles_kokkos(
                genotype_alleles, requests, 2, 2);
            const auto posterior = fastgatk::kernels::calculate_site_posterior_kokkos(
                genotype_pl, requests, 2, 2,
                {0.0, -3.4771212547196626, -6.477121254719663});
            const auto cross_sample_reference =
                fastgatk::kernels::calculate_cross_sample_reference_confidence_kokkos(
                    genotype_pl, requests, 2, 2,
                    {0.0, -3.4771212547196626, -6.477121254719663});
            const auto cohort = fastgatk::kernels::calculate_allele_frequency_kokkos(
                genotype_pl, requests, 2, 2, {10.0, 0.01});
            const auto rcm = fastgatk::kernels::calculate_reference_confidence_kokkos(
                rcm_locus, rcm_quality, rcm_alt, 128);
            const auto rcm_poly =
                fastgatk::kernels::calculate_reference_confidence_genotypes_kokkos(
                    rcm_locus, rcm_poly_reference, rcm_poly_non_ref, 128, joint_ploidy);
            const auto somatic = fastgatk::kernels::calculate_somatic_likelihood_kokkos(
                somatic_reference_likelihoods, somatic_alternate_likelihoods,
                somatic_candidates, somatic_reads, 101);
            const auto somatic_posterior = fastgatk::kernels::calculate_somatic_posterior_kokkos(
                somatic_reference_likelihoods, somatic_alternate_likelihoods,
                {}, {}, somatic_f1r2, somatic_r1f2,
                somatic_candidates, somatic_reads, 0.05);
            const auto somatic_multiallelic =
                fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                    somatic_multiallelic_likelihoods,
                    somatic_multiallelic_allele_count, somatic_reads);
            const auto bqsr_count = fastgatk::kernels::count_bqsr_quality_kokkos(bqsr_input);
            const auto bqsr_covariate = fastgatk::kernels::count_bqsr_covariates_kokkos(
                bqsr_covariate_input);
            const auto bqsr_apply = fastgatk::kernels::apply_bqsr_quality_kokkos(
                bqsr_input.qualities, bqsr_deltas);
            sw_prepare_samples.push_back(sw.prepare_seconds);
            sw_kernel_samples.push_back(sw.seconds);
            sw_uniform_prepare_samples.push_back(sw_uniform.prepare_seconds);
            sw_uniform_kernel_samples.push_back(sw_uniform.seconds);
            sw_simd_width = sw_uniform.simd_width;
            sw_simd_groups = sw_uniform.simd_groups;
            pairhmm_prepare_samples.push_back(pairhmm.prepare_seconds);
            pairhmm_kernel_samples.push_back(pairhmm.seconds);
            float_pairhmm_prepare_samples.push_back(float_pairhmm.prepare_seconds);
            float_pairhmm_kernel_samples.push_back(float_pairhmm.seconds);
            float_pairhmm_simd_width = float_pairhmm.simd_width;
            normalization_prepare_samples.push_back(normalized.prepare_seconds);
            normalization_kernel_samples.push_back(normalized.seconds);
            marginalization_prepare_samples.push_back(marginalized.prepare_seconds);
            marginalization_kernel_samples.push_back(marginalized.seconds);
            uncertainty_prepare_samples.push_back(uncertainty.prepare_seconds);
            uncertainty_kernel_samples.push_back(uncertainty.seconds);
            best_allele_prepare_samples.push_back(best_alleles.prepare_seconds);
            best_allele_kernel_samples.push_back(best_alleles.seconds);
            flow_pairhmm_prepare_samples.push_back(flow_pairhmm.prepare_seconds);
            flow_pairhmm_kernel_samples.push_back(flow_pairhmm.seconds);
            flow_pairhmm_execution_policy = flow_pairhmm.execution_policy;
            persistent_prepare_samples.push_back(persistent.prepare_seconds);
            persistent_kernel_samples.push_back(persistent.seconds);
            genotype_prepare_samples.push_back(genotype.prepare_seconds);
            genotype_kernel_samples.push_back(genotype.seconds);
            joint_genotype_prepare_samples.push_back(joint_genotype.prepare_seconds);
            joint_genotype_kernel_samples.push_back(joint_genotype.seconds);
            genotype_prior_prepare_samples.push_back(genotype_prior.prepare_seconds);
            genotype_prior_kernel_samples.push_back(genotype_prior.seconds);
            posterior_assignment_prepare_samples.push_back(
                posterior_assignment.prepare_seconds);
            posterior_assignment_kernel_samples.push_back(
                posterior_assignment.seconds);
            allele_count_prepare_samples.push_back(allele_counts.prepare_seconds);
            allele_count_kernel_samples.push_back(allele_counts.seconds);
            posterior_prepare_samples.push_back(posterior.prepare_seconds);
            posterior_kernel_samples.push_back(posterior.seconds);
            cross_sample_reference_prepare_samples.push_back(
                cross_sample_reference.prepare_seconds);
            cross_sample_reference_kernel_samples.push_back(
                cross_sample_reference.seconds);
            cohort_prepare_samples.push_back(cohort.prepare_seconds);
            cohort_kernel_samples.push_back(cohort.seconds);
            rcm_prepare_samples.push_back(rcm.prepare_seconds);
            rcm_kernel_samples.push_back(rcm.seconds);
            rcm_poly_prepare_samples.push_back(rcm_poly.prepare_seconds);
            rcm_poly_kernel_samples.push_back(rcm_poly.seconds);
            somatic_prepare_samples.push_back(somatic.prepare_seconds);
            somatic_kernel_samples.push_back(somatic.seconds);
            somatic_posterior_prepare_samples.push_back(somatic_posterior.prepare_seconds);
            somatic_posterior_kernel_samples.push_back(somatic_posterior.seconds);
            somatic_multiallelic_prepare_samples.push_back(
                somatic_multiallelic.prepare_seconds);
            somatic_multiallelic_kernel_samples.push_back(somatic_multiallelic.seconds);
            bqsr_count_prepare_samples.push_back(bqsr_count.prepare_seconds);
            bqsr_count_kernel_samples.push_back(bqsr_count.execute_seconds);
            bqsr_covariate_prepare_samples.push_back(bqsr_covariate.prepare_seconds);
            bqsr_covariate_kernel_samples.push_back(bqsr_covariate.execute_seconds);
            bqsr_apply_prepare_samples.push_back(bqsr_apply.prepare_seconds);
            bqsr_apply_kernel_samples.push_back(bqsr_apply.execute_seconds);
            const auto current_sw_checksum =
                std::accumulate(sw.scores.begin(), sw.scores.end(), std::uint64_t{0});
            const auto current_sw_uniform_checksum =
                std::accumulate(sw_uniform.scores.begin(), sw_uniform.scores.end(), std::uint64_t{0});
            const auto current_pairhmm_checksum = pairhmm.checksum;
            const auto current_float_pairhmm_checksum = float_pairhmm.checksum;
            double current_float_pairhmm_max_abs_error = 0.0;
            for (std::size_t index = 0; index < pairhmm.likelihoods.size(); ++index)
                current_float_pairhmm_max_abs_error = std::max(
                    current_float_pairhmm_max_abs_error,
                    std::abs(pairhmm.likelihoods[index] - float_pairhmm.likelihoods[index]));
            if (iteration == 0) {
                sw_checksum = current_sw_checksum;
                sw_uniform_checksum = current_sw_uniform_checksum;
                pairhmm_checksum = current_pairhmm_checksum;
                float_pairhmm_checksum = current_float_pairhmm_checksum;
                float_pairhmm_max_abs_error = current_float_pairhmm_max_abs_error;
                normalization_checksum = 0;
                for (const auto value : normalized.likelihoods)
                    normalization_checksum += static_cast<std::uint64_t>(
                        std::llround((value + 20.0) * 1'000'000.0));
                marginalization_checksum = 0;
                for (const auto value : marginalized.best_by_row_allele)
                    if (std::isfinite(value))
                        marginalization_checksum += static_cast<std::uint64_t>(
                            std::llround((value + 20.0) * 1'000'000.0));
                uncertainty_checksum = 0;
                for (const auto value : uncertainty.best_second_by_row)
                    if (std::isfinite(value))
                        uncertainty_checksum += static_cast<std::uint64_t>(
                            std::llround((value + 20.0) * 1'000'000.0));
                best_allele_checksum = 0;
                for (const auto value : best_alleles.best_second_allele_by_row)
                    if (value != std::numeric_limits<std::uint32_t>::max())
                        best_allele_checksum += value + 1U;
                for (const auto value : best_alleles.best_second_likelihood_by_row)
                    if (std::isfinite(value))
                        best_allele_checksum += static_cast<std::uint64_t>(
                            std::llround((value + 20.0) * 1'000'000.0));
                flow_pairhmm_checksum = flow_pairhmm.checksum;
                persistent_checksum = persistent.checksum;
                genotype_checksum = 0;
                for (std::size_t sample = 0; sample < requests; ++sample)
                    genotype_checksum += static_cast<std::uint64_t>(
                        genotype.first_allele[sample] + genotype.second_allele[sample] +
                        genotype.gq[sample]);
                joint_genotype_checksum = 0;
                for (const auto value : joint_genotype.pl)
                    joint_genotype_checksum += static_cast<std::uint64_t>(value + 1);
                genotype_prior_checksum = 0;
                for (const auto value : genotype_prior.log10_priors)
                    genotype_prior_checksum += static_cast<std::uint64_t>(
                        std::llround((value + 20.0) * 1'000'000.0));
                posterior_assignment_checksum = 0;
                for (std::size_t sample = 0; sample < requests; ++sample)
                    posterior_assignment_checksum += static_cast<std::uint64_t>(
                        posterior_assignment.first_allele[sample] +
                        posterior_assignment.second_allele[sample] +
                        posterior_assignment.gq[sample]);
                allele_count_checksum = 0;
                for (const auto count : allele_counts.counts)
                    allele_count_checksum += static_cast<std::uint64_t>(count);
                posterior_checksum = static_cast<std::uint64_t>(
                    posterior.samples_with_likelihoods * 1000 + posterior.qual);
                cross_sample_reference_checksum = static_cast<std::uint64_t>(
                    cross_sample_reference.samples_with_likelihoods * 1000 +
                    cross_sample_reference.joint_qual);
                cohort_checksum = static_cast<std::uint64_t>(
                    cohort.samples_with_likelihoods * 1000 + cohort.integer_allele_counts[1] +
                    cohort.iterations);
                rcm_checksum = static_cast<std::uint64_t>(
                    rcm.depth[0] * 1000 + rcm.non_ref_count[0]);
                rcm_poly_checksum = 0;
                for (const auto value : rcm_poly.genotype_pl)
                    rcm_poly_checksum += static_cast<std::uint64_t>(value + 1);
                for (const auto value : rcm_poly.gq)
                    rcm_poly_checksum += value;
                somatic_checksum = static_cast<std::uint64_t>(
                    somatic.informative_reads[0] + somatic.tlod[0] * 1000.0);
                somatic_posterior_checksum = static_cast<std::uint64_t>(
                    somatic_posterior.informative_reads[0] +
                    somatic_posterior.somatic_probability[0] * 1000000.0 +
                    somatic_posterior.artifact_probability[0] * 1000000.0);
                somatic_multiallelic_checksum = 0;
                for (std::size_t alternate = 0;
                     alternate < somatic_multiallelic.tlod.size(); ++alternate) {
                    somatic_multiallelic_checksum += static_cast<std::uint64_t>(
                        somatic_multiallelic.informative_reads[alternate]);
                    somatic_multiallelic_checksum += static_cast<std::uint64_t>(
                        std::llround((somatic_multiallelic.tlod[alternate] + 100.0) * 1000.0));
                }
                bqsr_count_checksum = 0;
                for (const auto& bin : bqsr_count.bins)
                    bqsr_count_checksum += bin.count + bin.mismatches;
                bqsr_covariate_checksum = 0;
                for (const auto& bin : bqsr_covariate.covariates)
                    bqsr_covariate_checksum += bin.count + bin.mismatches;
                bqsr_apply_checksum = std::accumulate(
                    bqsr_apply.adjusted.begin(), bqsr_apply.adjusted.end(), std::uint64_t{0});
            } else if (sw_checksum != current_sw_checksum ||
                       sw_uniform_checksum != current_sw_uniform_checksum ||
                       pairhmm_checksum != current_pairhmm_checksum ||
                       float_pairhmm_checksum != current_float_pairhmm_checksum ||
                       float_pairhmm_max_abs_error != current_float_pairhmm_max_abs_error ||
                       persistent_checksum != persistent.checksum ||
                       flow_pairhmm_checksum != flow_pairhmm.checksum) {
                throw std::runtime_error("non-deterministic kernel benchmark checksum");
            }
            std::uint64_t current_normalization_checksum = 0;
            for (const auto value : normalized.likelihoods)
                current_normalization_checksum += static_cast<std::uint64_t>(
                    std::llround((value + 20.0) * 1'000'000.0));
            if (normalization_checksum != current_normalization_checksum)
                throw std::runtime_error("non-deterministic likelihood-normalization checksum");
            std::uint64_t current_marginalization_checksum = 0;
            for (const auto value : marginalized.best_by_row_allele)
                if (std::isfinite(value))
                    current_marginalization_checksum += static_cast<std::uint64_t>(
                        std::llround((value + 20.0) * 1'000'000.0));
            if (marginalization_checksum != current_marginalization_checksum)
                throw std::runtime_error(
                    "non-deterministic allele-marginalization checksum");
            std::uint64_t current_uncertainty_checksum = 0;
            for (const auto value : uncertainty.best_second_by_row)
                if (std::isfinite(value))
                    current_uncertainty_checksum += static_cast<std::uint64_t>(
                        std::llround((value + 20.0) * 1'000'000.0));
            if (uncertainty_checksum != current_uncertainty_checksum)
                throw std::runtime_error(
                    "non-deterministic read-allele-uncertainty checksum");
            std::uint64_t current_best_allele_checksum = 0;
            for (const auto value : best_alleles.best_second_allele_by_row)
                if (value != std::numeric_limits<std::uint32_t>::max())
                    current_best_allele_checksum += value + 1U;
            for (const auto value : best_alleles.best_second_likelihood_by_row)
                if (std::isfinite(value))
                    current_best_allele_checksum += static_cast<std::uint64_t>(
                        std::llround((value + 20.0) * 1'000'000.0));
            if (best_allele_checksum != current_best_allele_checksum)
                throw std::runtime_error(
                    "non-deterministic read-allele-best checksum");
            std::uint64_t current_genotype_checksum = 0;
            for (std::size_t sample = 0; sample < requests; ++sample)
                current_genotype_checksum += static_cast<std::uint64_t>(
                    genotype.first_allele[sample] + genotype.second_allele[sample] +
                    genotype.gq[sample]);
            if (genotype_checksum != current_genotype_checksum)
                throw std::runtime_error("non-deterministic genotype benchmark checksum");
            std::uint64_t current_joint_genotype_checksum = 0;
            for (const auto value : joint_genotype.pl)
                current_joint_genotype_checksum += static_cast<std::uint64_t>(value + 1);
            if (joint_genotype_checksum != current_joint_genotype_checksum)
                throw std::runtime_error("non-deterministic joint genotype benchmark checksum");
            std::uint64_t current_genotype_prior_checksum = 0;
            for (const auto value : genotype_prior.log10_priors)
                current_genotype_prior_checksum += static_cast<std::uint64_t>(
                    std::llround((value + 20.0) * 1'000'000.0));
            if (genotype_prior_checksum != current_genotype_prior_checksum)
                throw std::runtime_error("non-deterministic genotype-prior benchmark checksum");
            std::uint64_t current_posterior_assignment_checksum = 0;
            for (std::size_t sample = 0; sample < requests; ++sample)
                current_posterior_assignment_checksum += static_cast<std::uint64_t>(
                    posterior_assignment.first_allele[sample] +
                    posterior_assignment.second_allele[sample] +
                    posterior_assignment.gq[sample]);
            if (posterior_assignment_checksum != current_posterior_assignment_checksum)
                throw std::runtime_error(
                    "non-deterministic posterior assignment benchmark checksum");
            std::uint64_t current_allele_count_checksum = 0;
            for (const auto count : allele_counts.counts)
                current_allele_count_checksum += static_cast<std::uint64_t>(count);
            if (allele_count_checksum != current_allele_count_checksum)
                throw std::runtime_error("non-deterministic allele-count benchmark checksum");
            const auto current_posterior_checksum = static_cast<std::uint64_t>(
                posterior.samples_with_likelihoods * 1000 + posterior.qual);
            if (posterior_checksum != current_posterior_checksum)
                throw std::runtime_error("non-deterministic posterior benchmark checksum");
            const auto current_cross_sample_reference_checksum = static_cast<std::uint64_t>(
                cross_sample_reference.samples_with_likelihoods * 1000 +
                cross_sample_reference.joint_qual);
            if (cross_sample_reference_checksum != current_cross_sample_reference_checksum)
                throw std::runtime_error("non-deterministic cross-sample reference benchmark checksum");
            const auto current_cohort_checksum = static_cast<std::uint64_t>(
                cohort.samples_with_likelihoods * 1000 + cohort.integer_allele_counts[1] +
                cohort.iterations);
            if (cohort_checksum != current_cohort_checksum)
                throw std::runtime_error("non-deterministic cohort benchmark checksum");
            const auto current_rcm_checksum = static_cast<std::uint64_t>(
                rcm.depth[0] * 1000 + rcm.non_ref_count[0]);
            if (rcm_checksum != current_rcm_checksum)
                throw std::runtime_error("non-deterministic reference-confidence benchmark checksum");
            std::uint64_t current_rcm_poly_checksum = 0;
            for (const auto value : rcm_poly.genotype_pl)
                current_rcm_poly_checksum += static_cast<std::uint64_t>(value + 1);
            for (const auto value : rcm_poly.gq)
                current_rcm_poly_checksum += value;
            if (rcm_poly_checksum != current_rcm_poly_checksum)
                throw std::runtime_error(
                    "non-deterministic polyploid reference-confidence benchmark checksum");
            const auto current_somatic_checksum = static_cast<std::uint64_t>(
                somatic.informative_reads[0] + somatic.tlod[0] * 1000.0);
            if (somatic_checksum != current_somatic_checksum)
                throw std::runtime_error("non-deterministic somatic benchmark checksum");
            const auto current_somatic_posterior_checksum = static_cast<std::uint64_t>(
                somatic_posterior.informative_reads[0] +
                somatic_posterior.somatic_probability[0] * 1000000.0 +
                somatic_posterior.artifact_probability[0] * 1000000.0);
            if (somatic_posterior_checksum != current_somatic_posterior_checksum)
                throw std::runtime_error("non-deterministic somatic posterior benchmark checksum");
            std::uint64_t current_somatic_multiallelic_checksum = 0;
            for (std::size_t alternate = 0;
                 alternate < somatic_multiallelic.tlod.size(); ++alternate) {
                current_somatic_multiallelic_checksum += static_cast<std::uint64_t>(
                    somatic_multiallelic.informative_reads[alternate]);
                current_somatic_multiallelic_checksum += static_cast<std::uint64_t>(
                    std::llround((somatic_multiallelic.tlod[alternate] + 100.0) * 1000.0));
            }
            if (somatic_multiallelic_checksum != current_somatic_multiallelic_checksum)
                throw std::runtime_error(
                    "non-deterministic multiallelic somatic benchmark checksum");
            std::uint64_t current_bqsr_count_checksum = 0;
            for (const auto& bin : bqsr_count.bins)
                current_bqsr_count_checksum += bin.count + bin.mismatches;
            if (bqsr_count_checksum != current_bqsr_count_checksum)
                throw std::runtime_error("non-deterministic BQSR count benchmark checksum");
            std::uint64_t current_bqsr_covariate_checksum = 0;
            for (const auto& bin : bqsr_covariate.covariates)
                current_bqsr_covariate_checksum += bin.count + bin.mismatches;
            if (bqsr_covariate_checksum != current_bqsr_covariate_checksum)
                throw std::runtime_error("non-deterministic BQSR covariate benchmark checksum");
            const auto current_bqsr_apply_checksum = std::accumulate(
                bqsr_apply.adjusted.begin(), bqsr_apply.adjusted.end(), std::uint64_t{0});
            if (bqsr_apply_checksum != current_bqsr_apply_checksum)
                throw std::runtime_error("non-deterministic BQSR apply benchmark checksum");
            persistent_cached_shapes = persistent.cached_shapes;
            persistent_cache_hits = persistent.cache_hits;
            execution_space = pairhmm.execution_space;
            normalization_execution_space = normalized.execution_space;
            marginalization_execution_space = marginalized.execution_space;
            uncertainty_execution_space = uncertainty.execution_space;
            best_allele_execution_space = best_alleles.execution_space;
            simd_width = pairhmm.simd_width;
            genotype_execution_space = genotype.execution_space;
            joint_genotype_execution_space = joint_genotype.execution_space;
            genotype_prior_execution_space = genotype_prior.execution_space;
            posterior_assignment_execution_space = posterior_assignment.execution_space;
            allele_count_execution_space = allele_counts.execution_space;
            posterior_execution_space = posterior.execution_space;
            cross_sample_reference_execution_space = cross_sample_reference.execution_space;
            cohort_execution_space = cohort.execution_space;
            rcm_execution_space = rcm.execution_space;
            rcm_poly_execution_space = rcm_poly.execution_space;
            somatic_execution_space = somatic.execution_space;
            somatic_posterior_execution_space = somatic_posterior.execution_space;
            somatic_multiallelic_execution_space = somatic_multiallelic.execution_space;
            bqsr_count_execution_space = bqsr_count.execution_space;
            bqsr_count_execution_policy = bqsr_count.execution_policy;
            bqsr_covariate_execution_space = bqsr_covariate.execution_space;
            bqsr_covariate_execution_policy = bqsr_covariate.execution_policy;
            bqsr_covariate_workspace_bytes = bqsr_covariate.workspace_bytes;
            bqsr_covariate_team_local_histogram = bqsr_covariate.team_local_histogram;
            bqsr_apply_execution_space = bqsr_apply.execution_space;
            bqsr_apply_execution_policy = bqsr_apply.execution_policy;
        }

        const auto sw_prepare = summarize(std::move(sw_prepare_samples));
        const auto sw_kernel = summarize(std::move(sw_kernel_samples));
        const auto sw_uniform_prepare = summarize(std::move(sw_uniform_prepare_samples));
        const auto sw_uniform_kernel = summarize(std::move(sw_uniform_kernel_samples));
        const auto pairhmm_prepare = summarize(std::move(pairhmm_prepare_samples));
        const auto pairhmm_kernel = summarize(std::move(pairhmm_kernel_samples));
        const auto float_pairhmm_prepare =
            summarize(std::move(float_pairhmm_prepare_samples));
        const auto float_pairhmm_kernel =
            summarize(std::move(float_pairhmm_kernel_samples));
        const auto normalization_prepare =
            summarize(std::move(normalization_prepare_samples));
        const auto normalization_kernel =
            summarize(std::move(normalization_kernel_samples));
        const auto marginalization_prepare =
            summarize(std::move(marginalization_prepare_samples));
        const auto marginalization_kernel =
            summarize(std::move(marginalization_kernel_samples));
        const auto uncertainty_prepare = summarize(std::move(uncertainty_prepare_samples));
        const auto uncertainty_kernel = summarize(std::move(uncertainty_kernel_samples));
        const auto best_allele_prepare = summarize(std::move(best_allele_prepare_samples));
        const auto best_allele_kernel = summarize(std::move(best_allele_kernel_samples));
        const auto flow_pairhmm_prepare = summarize(std::move(flow_pairhmm_prepare_samples));
        const auto flow_pairhmm_kernel = summarize(std::move(flow_pairhmm_kernel_samples));
        const auto persistent_prepare = summarize(std::move(persistent_prepare_samples));
        const auto persistent_kernel = summarize(std::move(persistent_kernel_samples));
        const auto genotype_prepare = summarize(std::move(genotype_prepare_samples));
        const auto genotype_kernel = summarize(std::move(genotype_kernel_samples));
        const auto joint_genotype_prepare =
            summarize(std::move(joint_genotype_prepare_samples));
        const auto joint_genotype_kernel =
            summarize(std::move(joint_genotype_kernel_samples));
        const auto genotype_prior_prepare =
            summarize(std::move(genotype_prior_prepare_samples));
        const auto genotype_prior_kernel =
            summarize(std::move(genotype_prior_kernel_samples));
        const auto posterior_assignment_prepare =
            summarize(std::move(posterior_assignment_prepare_samples));
        const auto posterior_assignment_kernel =
            summarize(std::move(posterior_assignment_kernel_samples));
        const auto allele_count_prepare = summarize(std::move(allele_count_prepare_samples));
        const auto allele_count_kernel = summarize(std::move(allele_count_kernel_samples));
        const auto posterior_prepare = summarize(std::move(posterior_prepare_samples));
        const auto posterior_kernel = summarize(std::move(posterior_kernel_samples));
        const auto cross_sample_reference_prepare =
            summarize(std::move(cross_sample_reference_prepare_samples));
        const auto cross_sample_reference_kernel =
            summarize(std::move(cross_sample_reference_kernel_samples));
        const auto cohort_prepare = summarize(std::move(cohort_prepare_samples));
        const auto cohort_kernel = summarize(std::move(cohort_kernel_samples));
        const auto rcm_prepare = summarize(std::move(rcm_prepare_samples));
        const auto rcm_kernel = summarize(std::move(rcm_kernel_samples));
        const auto rcm_poly_prepare = summarize(std::move(rcm_poly_prepare_samples));
        const auto rcm_poly_kernel = summarize(std::move(rcm_poly_kernel_samples));
        const auto somatic_prepare = summarize(std::move(somatic_prepare_samples));
        const auto somatic_kernel = summarize(std::move(somatic_kernel_samples));
        const auto somatic_posterior_prepare = summarize(std::move(somatic_posterior_prepare_samples));
        const auto somatic_posterior_kernel = summarize(std::move(somatic_posterior_kernel_samples));
        const auto somatic_multiallelic_prepare =
            summarize(std::move(somatic_multiallelic_prepare_samples));
        const auto somatic_multiallelic_kernel =
            summarize(std::move(somatic_multiallelic_kernel_samples));
        const auto bqsr_count_prepare = summarize(std::move(bqsr_count_prepare_samples));
        const auto bqsr_count_kernel = summarize(std::move(bqsr_count_kernel_samples));
        const auto bqsr_covariate_prepare = summarize(std::move(bqsr_covariate_prepare_samples));
        const auto bqsr_covariate_kernel = summarize(std::move(bqsr_covariate_kernel_samples));
        const auto bqsr_apply_prepare = summarize(std::move(bqsr_apply_prepare_samples));
        const auto bqsr_apply_kernel = summarize(std::move(bqsr_apply_kernel_samples));
        fastgatk::kernels::KmerGraphInput graph_input;
        graph_input.offsets.push_back(0);
        for (const auto& read : reads) {
            graph_input.bases.insert(graph_input.bases.end(), read.bases.begin(), read.bases.end());
            graph_input.offsets.push_back(static_cast<std::uint32_t>(graph_input.bases.size()));
        }
        graph_input.reference_offsets = {0, 1024};
        graph_input.reference_bases.resize(1024);
        for (std::size_t i = 0; i < graph_input.reference_bases.size(); ++i)
            graph_input.reference_bases[i] = static_cast<std::uint8_t>(i & 3U);
        const auto graph = fastgatk::kernels::build_kmer_graph_kokkos(
            graph_input, fastgatk::kernels::KmerGraphOptions{5, 1});
        fastgatk::kernels::ReadErrorCorrectionInput correction_input;
        correction_input.offsets.push_back(0);
        for (const auto& read : reads) {
            correction_input.bases.insert(correction_input.bases.end(),
                                          read.bases.begin(), read.bases.end());
            correction_input.qualities.insert(correction_input.qualities.end(),
                                              read.qualities.begin(), read.qualities.end());
            correction_input.offsets.push_back(
                static_cast<std::uint32_t>(correction_input.bases.size()));
        }
        const auto correction = fastgatk::kernels::correct_read_errors_kokkos(
            correction_input, fastgatk::kernels::ReadErrorCorrectionOptions{5, 2, 2, 30, 1});
        fastgatk::kernels::ActivityProfileInput activity_input;
        activity_input.tids.assign(requests, 0);
        activity_input.positions.resize(requests);
        activity_input.counts.resize(requests * 4);
        activity_input.reference_bases.assign(requests, 0);
        for (std::size_t i = 0; i < requests; ++i) {
            activity_input.positions[i] = static_cast<std::int32_t>(i);
            activity_input.counts[i * 4] = 28;
            activity_input.counts[i * 4 + 1] = (i % 17 == 0) ? 4 : 2;
            activity_input.counts[i * 4 + 2] = 0;
            activity_input.counts[i * 4 + 3] = 0;
        }
        fastgatk::kernels::ActivityProfileOptions activity_options;
        activity_options.min_depth = 4;
        activity_options.min_activity = 0.002;
        activity_options.halo = 32;
        activity_options.max_region_size = 300;
        activity_options.max_prob_propagation_distance = 300;
        const auto activity = fastgatk::kernels::compute_activity_profile_kokkos(
            activity_input, activity_options);
        std::cout << "{\"schema_version\":7,\"status\":\"pass\",\"requests\":" << requests
                  << ",\"warmup\":" << warmup << ",\"repeats\":" << repeats
                  << ",\"backend\":\"" << execution_space << "\""
                  << ",\"simd_width\":" << simd_width
                  << ",\"compiled_simd_backend\":\""
                  << fastgatk::pairhmm::compiled_simd_backend_name() << "\""
                  << ",\"compiled_simd_width\":"
                  << fastgatk::pairhmm::compiled_simd_width()
                  << ",\"cpu_supports_avx2\":"
                  << (fastgatk::pairhmm::cpu_supports_avx2() ? "true" : "false")
                  << ",\"cpu_supports_avx512\":"
                  << (fastgatk::pairhmm::cpu_supports_avx512() ? "true" : "false")
                  << ",\"sw_simd_width\":" << sw_simd_width
                  << ",\"sw_simd_groups\":" << sw_simd_groups
                  << ",\"concurrency\":" << Kokkos::DefaultExecutionSpace{}.concurrency()
                  << ",\"compiler\":\"" << __VERSION__ << "\""
                  << ",\"flow_pairhmm_error_model\":\"flow\""
                  << ",\"flow_pairhmm_execution_policy\":\""
                  << flow_pairhmm_execution_policy << "\""
                  << ",\"hardware_concurrency\":" << std::thread::hardware_concurrency()
                  << ",\"sw_checksum\":" << sw_checksum
                  << ",\"sw_uniform_checksum\":" << sw_uniform_checksum
                  << ",\"pairhmm_checksum\":" << pairhmm_checksum
                  << ",\"float_pairhmm_checksum\":" << float_pairhmm_checksum
                  << ",\"float_pairhmm_max_abs_error\":" << float_pairhmm_max_abs_error
                  << ",\"float_pairhmm_simd_width\":" << float_pairhmm_simd_width
                  << ",\"pairhmm_likelihood_normalization_checksum\":"
                  << normalization_checksum
                  << ",\"pairhmm_allele_marginalization_checksum\":"
                  << marginalization_checksum
                  << ",\"pairhmm_read_allele_uncertainty_checksum\":"
                  << uncertainty_checksum
                  << ",\"pairhmm_read_allele_best_checksum\":"
                  << best_allele_checksum
                  << ",\"flow_pairhmm_checksum\":" << flow_pairhmm_checksum
                  << ",\"persistent_pairhmm_checksum\":" << persistent_checksum
                  << ",\"genotype_checksum\":" << genotype_checksum
                  << ",\"joint_genotype_checksum\":" << joint_genotype_checksum
                  << ",\"genotype_prior_checksum\":" << genotype_prior_checksum
                  << ",\"posterior_assignment_checksum\":" << posterior_assignment_checksum
                  << ",\"allele_count_checksum\":" << allele_count_checksum
                  << ",\"posterior_checksum\":" << posterior_checksum
                  << ",\"cross_sample_reference_checksum\":" << cross_sample_reference_checksum
                  << ",\"cohort_checksum\":" << cohort_checksum
                  << ",\"reference_confidence_checksum\":" << rcm_checksum
                  << ",\"reference_confidence_polyploid_checksum\":" << rcm_poly_checksum
                  << ",\"somatic_checksum\":" << somatic_checksum
                  << ",\"somatic_posterior_checksum\":" << somatic_posterior_checksum
                  << ",\"somatic_multiallelic_checksum\":"
                  << somatic_multiallelic_checksum
                  << ",\"bqsr_count_checksum\":" << bqsr_count_checksum
                  << ",\"bqsr_covariate_checksum\":" << bqsr_covariate_checksum
                  << ",\"bqsr_apply_checksum\":" << bqsr_apply_checksum
                  << ",\"persistent_pairhmm_cached_shapes\":" << persistent_cached_shapes
                  << ",\"persistent_pairhmm_cache_hits\":" << persistent_cache_hits
                  << ",\"timings\":{\"prepare_and_execute\":{";
        emit_timing("sw_prepare", sw_prepare);
        std::cout << ',';
        emit_timing("sw_kernel", sw_kernel);
        std::cout << ',';
        emit_timing("sw_uniform_prepare", sw_uniform_prepare);
        std::cout << ',';
        emit_timing("sw_uniform_kernel", sw_uniform_kernel);
        std::cout << ',';
        emit_timing("pairhmm_prepare", pairhmm_prepare);
        std::cout << ',';
        emit_timing("pairhmm_kernel", pairhmm_kernel);
        std::cout << ',';
        emit_timing("float_pairhmm_prepare", float_pairhmm_prepare);
        std::cout << ',';
        emit_timing("float_pairhmm_kernel", float_pairhmm_kernel);
        std::cout << ',';
        emit_timing("pairhmm_likelihood_normalization_prepare", normalization_prepare);
        std::cout << ',';
        emit_timing("pairhmm_likelihood_normalization_kernel", normalization_kernel);
        std::cout << ',';
        emit_timing("pairhmm_allele_marginalization_prepare", marginalization_prepare);
        std::cout << ',';
        emit_timing("pairhmm_allele_marginalization_kernel", marginalization_kernel);
        std::cout << ',';
        emit_timing("pairhmm_read_allele_uncertainty_prepare", uncertainty_prepare);
        std::cout << ',';
        emit_timing("pairhmm_read_allele_uncertainty_kernel", uncertainty_kernel);
        std::cout << ',';
        emit_timing("pairhmm_read_allele_best_prepare", best_allele_prepare);
        std::cout << ',';
        emit_timing("pairhmm_read_allele_best_kernel", best_allele_kernel);
        std::cout << ',';
        emit_timing("flow_pairhmm_prepare", flow_pairhmm_prepare);
        std::cout << ',';
        emit_timing("flow_pairhmm_kernel", flow_pairhmm_kernel);
        std::cout << ',';
        emit_timing("persistent_pairhmm_prepare", persistent_prepare);
        std::cout << ',';
        emit_timing("persistent_pairhmm_kernel", persistent_kernel);
        std::cout << ',';
        emit_timing("genotype_prepare", genotype_prepare);
        std::cout << ',';
        emit_timing("genotype_kernel", genotype_kernel);
        std::cout << ',';
        emit_timing("joint_genotype_prepare", joint_genotype_prepare);
        std::cout << ',';
        emit_timing("joint_genotype_kernel", joint_genotype_kernel);
        std::cout << ',';
        emit_timing("genotype_prior_prepare", genotype_prior_prepare);
        std::cout << ',';
        emit_timing("genotype_prior_kernel", genotype_prior_kernel);
        std::cout << ',';
        emit_timing("genotype_posterior_assignment_prepare", posterior_assignment_prepare);
        std::cout << ',';
        emit_timing("genotype_posterior_assignment_kernel", posterior_assignment_kernel);
        std::cout << ',';
        emit_timing("genotype_allele_count_prepare", allele_count_prepare);
        std::cout << ',';
        emit_timing("genotype_allele_count_kernel", allele_count_kernel);
        std::cout << ',';
        emit_timing("genotype_posterior_prepare", posterior_prepare);
        std::cout << ',';
        emit_timing("genotype_posterior_kernel", posterior_kernel);
        std::cout << ',';
        emit_timing("cross_sample_reference_prepare", cross_sample_reference_prepare);
        std::cout << ',';
        emit_timing("cross_sample_reference_kernel", cross_sample_reference_kernel);
        std::cout << ',';
        emit_timing("genotype_cohort_prepare", cohort_prepare);
        std::cout << ',';
        emit_timing("genotype_cohort_kernel", cohort_kernel);
        std::cout << ',';
        emit_timing("reference_confidence_prepare", rcm_prepare);
        std::cout << ',';
        emit_timing("reference_confidence_kernel", rcm_kernel);
        std::cout << ',';
        emit_timing("reference_confidence_polyploid_prepare", rcm_poly_prepare);
        std::cout << ',';
        emit_timing("reference_confidence_polyploid_kernel", rcm_poly_kernel);
        std::cout << ',';
        emit_timing("somatic_prepare", somatic_prepare);
        std::cout << ',';
        emit_timing("somatic_kernel", somatic_kernel);
        std::cout << ',';
        emit_timing("somatic_posterior_prepare", somatic_posterior_prepare);
        std::cout << ',';
        emit_timing("somatic_posterior_kernel", somatic_posterior_kernel);
        std::cout << ',';
        emit_timing("somatic_multiallelic_prepare", somatic_multiallelic_prepare);
        std::cout << ',';
        emit_timing("somatic_multiallelic_kernel", somatic_multiallelic_kernel);
        std::cout << ',';
        emit_timing("bqsr_count_prepare", bqsr_count_prepare);
        std::cout << ',';
        emit_timing("bqsr_count_kernel", bqsr_count_kernel);
        std::cout << ',';
        emit_timing("bqsr_covariate_prepare", bqsr_covariate_prepare);
        std::cout << ',';
        emit_timing("bqsr_covariate_kernel", bqsr_covariate_kernel);
        std::cout << ',';
        emit_timing("bqsr_apply_prepare", bqsr_apply_prepare);
        std::cout << ',';
        emit_timing("bqsr_apply_kernel", bqsr_apply_kernel);
        std::cout << "},\"pairhmm_likelihood_normalization_execution_space\":\""
                  << normalization_execution_space
                  << "\",\"pairhmm_allele_marginalization_execution_space\":\""
                  << marginalization_execution_space
                  << "\",\"pairhmm_read_allele_uncertainty_execution_space\":\""
                  << uncertainty_execution_space
                  << "\",\"pairhmm_read_allele_best_execution_space\":\""
                  << best_allele_execution_space
                  << "\",\"genotype_execution_space\":\"" << genotype_execution_space
                  << "\",\"joint_genotype_execution_space\":\""
                  << joint_genotype_execution_space
                  << "\",\"genotype_prior_execution_space\":\""
                  << genotype_prior_execution_space
                  << "\",\"allele_count_execution_space\":\"" << allele_count_execution_space
                  << "\",\"posterior_assignment_execution_space\":\""
                  << posterior_assignment_execution_space
                  << "\",\"posterior_execution_space\":\"" << posterior_execution_space
                  << "\",\"cross_sample_reference_execution_space\":\""
                  << cross_sample_reference_execution_space
                  << "\",\"cohort_execution_space\":\"" << cohort_execution_space
                  << "\",\"reference_confidence_execution_space\":\""
                  << rcm_execution_space
                  << "\",\"reference_confidence_polyploid_execution_space\":\""
                  << rcm_poly_execution_space << "\",\"somatic_execution_space\":\""
                  << somatic_execution_space << "\",\"somatic_posterior_execution_space\":\""
                  << somatic_posterior_execution_space
                  << "\",\"somatic_multiallelic_execution_space\":\""
                  << somatic_multiallelic_execution_space
                  << "\",\"bqsr_count_execution_space\":\""
                  << bqsr_count_execution_space
                  << "\",\"bqsr_count_execution_policy\":\""
                  << bqsr_count_execution_policy
                  << "\",\"bqsr_covariate_execution_space\":\""
                  << bqsr_covariate_execution_space
                  << "\",\"bqsr_covariate_execution_policy\":\""
                  << bqsr_covariate_execution_policy
                  << "\",\"bqsr_covariate_workspace_bytes\":"
                  << bqsr_covariate_workspace_bytes
                  << ",\"bqsr_covariate_team_local_histogram\":"
                  << (bqsr_covariate_team_local_histogram ? "true" : "false")
                  << ",\"bqsr_apply_execution_space\":\""
                  << bqsr_apply_execution_space
                  << "\",\"bqsr_apply_execution_policy\":\""
                  << bqsr_apply_execution_policy << "\"},\"graph_nodes\":" << graph.nodes
                  << ",\"graph_edges\":" << graph.edges
                  << ",\"graph_reference_nodes\":" << graph.reference_nodes
                  << ",\"graph_reference_edges\":" << graph.reference_edges
                  << ",\"graph_reference_connected_nodes\":" << graph.reference_connected_nodes
                  << ",\"graph_dangling_nodes\":" << graph.dangling_nodes
                  << ",\"graph_dangling_branch_paths\":" << graph.dangling_branch_paths
                  << ",\"graph_dangling_branch_bases\":" << graph.dangling_branch_bases
                  << ",\"graph_dangling_recovered_paths\":" << graph.dangling_recovered_paths
                  << ",\"graph_dangling_recovered_bases\":" << graph.dangling_recovered_bases
                  << ",\"graph_artificial_haplotype_recovery_paths\":"
                  << graph.artificial_haplotype_recovery_paths
                  << ",\"graph_artificial_haplotype_recovery_bases\":"
                  << graph.artificial_haplotype_recovery_bases
                  << ",\"graph_seqgraph_nodes\":" << graph.seqgraph_nodes
                  << ",\"graph_seqgraph_edges\":" << graph.seqgraph_edges
                  << ",\"graph_seqgraph_linear_chain_merges\":" << graph.seqgraph_linear_chain_merges
                  << ",\"graph_seqgraph_diamond_merges\":" << graph.seqgraph_diamond_merges
                  << ",\"graph_seqgraph_tail_merges\":" << graph.seqgraph_tail_merges
                  << ",\"graph_seqgraph_suffix_splits\":" << graph.seqgraph_suffix_splits
                  << ",\"graph_seqgraph_suffix_merges\":" << graph.seqgraph_suffix_merges
                  << ",\"graph_adaptive_pruned_nodes\":" << graph.adaptive_pruned_nodes
                  << ",\"graph_non_unique_kmers\":" << graph.non_unique_kmers.size()
                  << ",\"graph_reference_non_unique_kmers\":" << graph.reference_non_unique_kmers
                  << ",\"graph_reference_kmer_rejected\":"
                  << (graph.reference_kmer_rejected ? "true" : "false")
                  << ",\"graph_reference_paths\":" << graph.reference_path_count
                  << ",\"graph_haplotype_paths\":" << graph.haplotype_path_count
                  << ",\"graph_haplotype_sequences\":" << graph.haplotype_path_sequences.size()
                  << ",\"graph_kmer_size_selected\":" << graph.kmer_size
                  << ",\"graph_kmer_iterations\":" << graph.kmer_iterations
                  << ",\"graph_has_non_reference_cycles\":"
                  << (graph.has_non_reference_cycles ? "true" : "false")
                  << ",\"graph_prepare_seconds\":" << graph.prepare_seconds
                  << ",\"graph_seconds\":" << graph.seconds
                  << ",\"graph_execution_space\":\"" << graph.execution_space << "\""
                  << ",\"error_correction_solid_kmers\":" << correction.solid_kmers
                  << ",\"error_correction_corrected_kmers\":" << correction.corrected_kmers
                  << ",\"error_correction_uncorrectable_kmers\":" << correction.uncorrectable_kmers
                  << ",\"error_correction_corrected_bases\":" << correction.corrected_bases
                  << ",\"error_correction_prepare_seconds\":" << correction.prepare_seconds
                  << ",\"error_correction_seconds\":" << correction.seconds
                  << ",\"error_correction_execution_space\":\""
                  << correction.execution_space << "\""
                  << ",\"activity_loci\":" << activity.loci
                  << ",\"activity_reference_aware\":"
                  << (activity.reference_aware ? "true" : "false")
                  << ",\"active_loci\":" << activity.active_loci
                  << ",\"assembly_regions\":" << activity.regions.size()
                  << ",\"activity_max_region_size\":" << activity_options.max_region_size
                  << ",\"activity_max_prob_propagation_distance\":"
                  << activity_options.max_prob_propagation_distance
                  << ",\"activity_filter_size\":" << activity.filter_size
                  << ",\"activity_effective_max_prob_propagation_distance\":"
                  << activity.effective_max_prob_propagation_distance
                  << ",\"activity_prepare_seconds\":" << activity.prepare_seconds
                  << ",\"activity_seconds\":" << activity.seconds
                  << ",\"activity_execution_space\":\"" << activity.execution_space << "\"}\n";
        persistent_pairhmm.clear();
        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
