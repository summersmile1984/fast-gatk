#include "fastgatk/kernels/genotype.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;
using InputView = Kokkos::View<std::int32_t*, MemorySpace>;
using LikelihoodView = Kokkos::View<double*, MemorySpace>;
using OutputView = Kokkos::View<std::int32_t*, MemorySpace>;
using PriorView = Kokkos::View<double*, MemorySpace>;
using PosteriorView = Kokkos::View<double*, MemorySpace>;
using CountView = Kokkos::View<std::int32_t*, MemorySpace>;

// GATK's MathUtils.approximateLog10SumLog10 uses this quantized Jacobian
// table for the two-allele path in GenotypeLikelihoodCalculator.  Keep the
// table dimensions and rounding rule explicit here: it is part of the PL
// numerical contract, not merely a performance optimization.
constexpr double kJacobianMaximumDifference = 8.0;
constexpr double kJacobianStep = 0.0001;
constexpr int kJacobianEntryCount =
    static_cast<int>(kJacobianMaximumDifference / kJacobianStep) + 1;

KOKKOS_INLINE_FUNCTION std::size_t genotype_count(int allele_count, int ploidy) {
    // C(A+P-1,P), capped to avoid integer overflow in malformed inputs.
    std::size_t result = 1;
    for (int step = 1; step <= ploidy; ++step) {
        const auto numerator = static_cast<std::size_t>(allele_count + step - 1);
        if (result > std::numeric_limits<std::size_t>::max() / numerator) return 0;
        result *= numerator;
        result /= static_cast<std::size_t>(step);
    }
    return result;
}

KOKKOS_INLINE_FUNCTION void unrank_genotype(std::size_t index, int allele_count,
                                             int ploidy, int* alleles) {
    int available_alleles = allele_count;
    for (int position = ploidy - 1; position >= 0; --position) {
        int selected = available_alleles - 1;
        for (int candidate = 0; candidate < available_alleles; ++candidate) {
            const auto group = genotype_count(candidate + 1, position);
            if (index < group) {
                selected = candidate;
                break;
            }
            index -= group;
        }
        alleles[position] = selected;
        available_alleles = selected + 1;
    }
}

}  // namespace

JointGenotypePlResult calculate_joint_genotype_pl_kokkos(
    const std::vector<double>& log10_likelihoods,
    const std::size_t read_count,
    const int allele_count,
    const int ploidy,
    const int max_phred) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (read_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid joint genotype likelihood dimensions");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || width > 1'000'000 ||
        log10_likelihoods.size() != read_count * static_cast<std::size_t>(allele_count))
        throw std::invalid_argument("joint genotype likelihood dimensions do not match");

    const auto prepare_begin = std::chrono::steady_clock::now();
    LikelihoodView device_likelihoods("joint_genotype_log10_likelihoods",
                                      log10_likelihoods.size());
    LikelihoodView device_rescaled_likelihoods("joint_genotype_rescaled_likelihoods",
                                                log10_likelihoods.size());
    LikelihoodView device_jacobian("joint_genotype_jacobian", kJacobianEntryCount);
    InputView device_complete_reads("joint_genotype_complete_reads", read_count);
    OutputView device_pl("joint_genotype_pl", width);
    Kokkos::View<double*, MemorySpace> device_scores("joint_genotype_scores", width);
    auto host_likelihoods = Kokkos::create_mirror_view(device_likelihoods);
    auto host_rescaled_likelihoods = Kokkos::create_mirror_view(device_rescaled_likelihoods);
    auto host_jacobian = Kokkos::create_mirror_view(device_jacobian);
    auto host_complete_reads = Kokkos::create_mirror_view(device_complete_reads);
    for (std::size_t index = 0; index < log10_likelihoods.size(); ++index)
        host_likelihoods(index) = log10_likelihoods[index];
    // This is the GATK multi-allelic rescaledNonLogLikelihoods recurrence:
    // normalize every read column by its best allele, then enter probability
    // space and retain the summed scale factor.  It is intentionally prepared
    // in allele-major/read-minor order to match the source calculator.
    std::vector<double> per_read_maxima(read_count,
        -std::numeric_limits<double>::infinity());
    std::size_t complete_read_count = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        bool complete = true;
        for (int allele = 0; allele < allele_count; ++allele) {
            const auto value = log10_likelihoods[
                static_cast<std::size_t>(allele) * read_count + read];
            if (!std::isfinite(value)) {
                complete = false;
                break;
            }
            per_read_maxima[read] = std::max(per_read_maxima[read], value);
        }
        host_complete_reads(read) = complete ? 1 : 0;
        if (complete) ++complete_read_count;
    }
    double log10_rescaling = 0.0;
    for (std::size_t read = 0; read < read_count; ++read)
        if (host_complete_reads(read) != 0)
            log10_rescaling += per_read_maxima[read];
    for (int allele = 0; allele < allele_count; ++allele)
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto index = static_cast<std::size_t>(allele) * read_count + read;
            host_rescaled_likelihoods(index) = host_complete_reads(read) != 0
                ? std::pow(10.0, log10_likelihoods[index] - per_read_maxima[read])
                : 0.0;
        }
    for (int index = 0; index < kJacobianEntryCount; ++index)
        host_jacobian(index) = std::log10(
            1.0 + std::pow(10.0, -static_cast<double>(index) * kJacobianStep));
    Kokkos::deep_copy(device_likelihoods, host_likelihoods);
    Kokkos::deep_copy(device_rescaled_likelihoods, host_rescaled_likelihoods);
    Kokkos::deep_copy(device_jacobian, host_jacobian);
    Kokkos::deep_copy(device_complete_reads, host_complete_reads);
    Kokkos::deep_copy(device_scores, -1.0e300);
    Kokkos::deep_copy(device_pl, 999);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_joint_genotype_scores",
        Kokkos::RangePolicy<ExecSpace>(0, width),
        KOKKOS_LAMBDA(const std::size_t genotype) {
            int alleles[32];
            unrank_genotype(genotype, allele_count, ploidy, alleles);
            int distinct_alleles[32];
            int distinct_counts[32];
            int component_count = 0;
            for (int copy = 0; copy < ploidy; ++copy) {
                if (copy == 0 || alleles[copy] != alleles[copy - 1]) {
                    distinct_alleles[component_count] = alleles[copy];
                    distinct_counts[component_count] = 0;
                    ++component_count;
                }
                ++distinct_counts[component_count - 1];
            }
            double score = 0.0;
            if (component_count == 1) {
                // GATK: MathUtils.sum(log10LikelihoodsByAlleleAndRead[allele]).
                const auto allele = distinct_alleles[0];
                for (std::size_t read = 0; read < read_count; ++read)
                    if (device_complete_reads(read) != 0)
                    score += device_likelihoods[
                        static_cast<std::size_t>(allele) * read_count + read];
            } else if (component_count == 2) {
                // GATK's exact branch selection includes a quantized
                // Jacobian log-sum for the biallelic mixture.
                const auto first_allele = distinct_alleles[0];
                const auto second_allele = distinct_alleles[1];
                const auto log10_first_count = Kokkos::log10(
                    static_cast<double>(distinct_counts[0]));
                const auto log10_second_count = Kokkos::log10(
                    static_cast<double>(distinct_counts[1]));
                for (std::size_t read = 0; read < read_count; ++read) {
                    if (device_complete_reads(read) == 0) continue;
                    double left = device_likelihoods[
                        static_cast<std::size_t>(first_allele) * read_count + read] +
                        log10_first_count;
                    double right = device_likelihoods[
                        static_cast<std::size_t>(second_allele) * read_count + read] +
                        log10_second_count;
                    if (left > right) {
                        const auto swap = left;
                        left = right;
                        right = swap;
                    }
                    if (left <= -1.0e299) {
                        score += right;
                        continue;
                    }
                    const auto difference = right - left;
                    if (difference < kJacobianMaximumDifference) {
                        const auto table_index = static_cast<int>(
                            difference / kJacobianStep + 0.5);
                        score += right + device_jacobian(table_index);
                    } else {
                        score += right;
                    }
                }
                score -= static_cast<double>(complete_read_count) * Kokkos::log10(
                    static_cast<double>(ploidy));
            } else {
                // GATK's 3+-allele path combines the rescaled probability
                // rows in allele-major order, then restores the column scale.
                for (std::size_t read = 0; read < read_count; ++read) {
                    if (device_complete_reads(read) == 0) continue;
                    double probability = 0.0;
                    for (int component = 0; component < component_count; ++component)
                        probability += static_cast<double>(distinct_counts[component]) *
                            device_rescaled_likelihoods[
                                static_cast<std::size_t>(distinct_alleles[component]) *
                                read_count + read];
                    score += Kokkos::log10(probability);
                }
                score -= static_cast<double>(complete_read_count) * Kokkos::log10(
                    static_cast<double>(ploidy));
                score += log10_rescaling;
            }
            device_scores(genotype) = score;
        });
    ExecSpace{}.fence();
    double maximum = -1.0e300;
    Kokkos::parallel_reduce(
        "fastgatk_joint_genotype_maximum",
        Kokkos::RangePolicy<ExecSpace>(0, width),
        KOKKOS_LAMBDA(const std::size_t genotype, double& local_max) {
            local_max = Kokkos::fmax(local_max, device_scores(genotype));
        }, Kokkos::Max<double>(maximum));
    ExecSpace{}.fence();
    Kokkos::parallel_for(
        "fastgatk_joint_genotype_pl",
        Kokkos::RangePolicy<ExecSpace>(0, width),
        KOKKOS_LAMBDA(const std::size_t genotype) {
            const auto score = device_scores(genotype);
            if (!(score > -1.0e299) || !(maximum > -1.0e299)) return;
            const auto phred = Kokkos::fmax(0.0, Kokkos::round(-10.0 * (score - maximum)));
            const auto bounded = max_phred > 0
                ? Kokkos::fmin(static_cast<double>(max_phred), phred) : phred;
            device_pl(genotype) = static_cast<std::int32_t>(Kokkos::fmin(
                static_cast<double>(std::numeric_limits<std::int32_t>::max()), bounded));
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_pl = Kokkos::create_mirror_view(device_pl);
    Kokkos::deep_copy(host_pl, device_pl);
    JointGenotypePlResult result;
    result.pl.resize(width, 999);
    for (std::size_t index = 0; index < width; ++index)
        result.pl[index] = host_pl(index);
    result.ploidy = ploidy;
    result.allele_count = allele_count;
    result.read_count = read_count;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

GenotypePriorResult calculate_genotype_priors_kokkos(
    const std::vector<double>& log10_heterozygous_priors,
    const std::vector<double>& log10_homozygous_priors,
    const int ploidy) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (log10_heterozygous_priors.size() != log10_homozygous_priors.size() ||
        log10_heterozygous_priors.size() < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid genotype-prior dimensions");
    for (std::size_t allele = 0; allele < log10_heterozygous_priors.size(); ++allele) {
        if (!std::isfinite(log10_heterozygous_priors[allele]) ||
            !std::isfinite(log10_homozygous_priors[allele]))
            throw std::invalid_argument("genotype priors must be finite");
    }
    const auto allele_count = static_cast<int>(log10_heterozygous_priors.size());
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || width > 1'000'000)
        throw std::invalid_argument("genotype-prior dimensions exceed the bounded contract");

    const auto prepare_begin = std::chrono::steady_clock::now();
    PriorView device_het("genotype_prior_heterozygous", log10_heterozygous_priors.size());
    PriorView device_hom("genotype_prior_homozygous", log10_homozygous_priors.size());
    PriorView device_output("genotype_prior_output", width);
    auto host_het = Kokkos::create_mirror_view(device_het);
    auto host_hom = Kokkos::create_mirror_view(device_hom);
    for (std::size_t allele = 0; allele < log10_heterozygous_priors.size(); ++allele) {
        host_het(allele) = log10_heterozygous_priors[allele];
        host_hom(allele) = log10_homozygous_priors[allele];
    }
    Kokkos::deep_copy(device_het, host_het);
    Kokkos::deep_copy(device_hom, host_hom);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_genotype_priors",
        Kokkos::RangePolicy<ExecSpace>(0, width),
        KOKKOS_LAMBDA(const std::size_t genotype) {
            int alleles[32];
            int counts[32];
            for (int allele = 0; allele < allele_count; ++allele) counts[allele] = 0;
            unrank_genotype(genotype, allele_count, ploidy, alleles);
            for (int copy = 0; copy < ploidy; ++copy) ++counts[alleles[copy]];
            double prior = 0.0;
            for (int allele = 1; allele < allele_count; ++allele) {
                const int count = counts[allele];
                if (count == 0) continue;
                const double het = device_het(allele);
                const double hom = device_hom(allele);
                prior += count == 2
                    ? hom
                    : het + (hom - het) * static_cast<double>(count - 1);
            }
            device_output(genotype) = prior;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_output = Kokkos::create_mirror_view(device_output);
    Kokkos::deep_copy(host_output, device_output);
    GenotypePriorResult result;
    result.log10_priors.resize(width, 0.0);
    for (std::size_t genotype = 0; genotype < width; ++genotype)
        result.log10_priors[genotype] = host_output(genotype);
    result.ploidy = ploidy;
    result.allele_count = allele_count;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

BiallelicCallConfidenceResult calculate_biallelic_call_confidence_kokkos(
    const std::vector<std::int32_t>& pl,
    const std::vector<double>& prior_pseudocounts) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (pl.empty() || pl.size() % 3 != 0)
        throw std::invalid_argument("biallelic call-confidence PL rows must have width three");
    const auto candidate_count = pl.size() / 3;
    if (prior_pseudocounts.size() != candidate_count * 2)
        throw std::invalid_argument("biallelic call-confidence prior dimensions do not match");
    for (const auto prior : prior_pseudocounts)
        if (!(prior > 0.0) || !std::isfinite(prior))
            throw std::invalid_argument("biallelic call-confidence pseudocounts must be positive and finite");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("biallelic_call_confidence_pl", pl.size());
    PriorView device_priors("biallelic_call_confidence_priors", prior_pseudocounts.size());
    PosteriorView device_qual("biallelic_call_confidence_qual", candidate_count);
    PosteriorView device_log10_p_alt_absent(
        "biallelic_call_confidence_log10_p_alt_absent", candidate_count);
    PosteriorView device_log10_p_variant_present(
        "biallelic_call_confidence_log10_p_variant_present", candidate_count);
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    for (std::size_t index = 0; index < pl.size(); ++index) host_pl(index) = pl[index];
    for (std::size_t index = 0; index < prior_pseudocounts.size(); ++index)
        host_priors(index) = prior_pseudocounts[index];
    Kokkos::deep_copy(device_pl, host_pl);
    Kokkos::deep_copy(device_priors, host_priors);
    Kokkos::deep_copy(device_qual, 0.0);
    Kokkos::deep_copy(device_log10_p_alt_absent, 0.0);
    Kokkos::deep_copy(device_log10_p_variant_present, 0.0);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_biallelic_call_confidence",
        Kokkos::RangePolicy<ExecSpace>(0, candidate_count),
        KOKKOS_LAMBDA(const std::size_t candidate) {
            const auto base = candidate * 3;
            if (device_pl(base) < 0 || device_pl(base + 1) < 0 || device_pl(base + 2) < 0)
                return;
            // The implementation follows AlleleFrequencyCalculator's
            // one-sample diploid EM.  Genotype 0/1 has a log10 combination
            // count of log10(2); hom-ref and hom-alt have coefficient one.
            double ref_frequency = 0.5;
            double alt_frequency = 0.5;
            double previous_ref_count = 0.0;
            double previous_alt_count = 0.0;
            for (int iteration = 0; iteration < 100; ++iteration) {
                const double log10_ref = Kokkos::log10(ref_frequency);
                const double log10_alt = Kokkos::log10(alt_frequency);
                const double terms[3] = {
                    -0.1 * static_cast<double>(device_pl(base)) + 2.0 * log10_ref,
                    -0.1 * static_cast<double>(device_pl(base + 1)) +
                        0.30102999566398119521 + log10_ref + log10_alt,
                    -0.1 * static_cast<double>(device_pl(base + 2)) + 2.0 * log10_alt};
                const double maximum = Kokkos::fmax(terms[0], Kokkos::fmax(terms[1], terms[2]));
                const double p0 = Kokkos::pow(10.0, terms[0] - maximum);
                const double p1 = Kokkos::pow(10.0, terms[1] - maximum);
                const double p2 = Kokkos::pow(10.0, terms[2] - maximum);
                const double denominator = p0 + p1 + p2;
                if (!(denominator > 0.0)) return;
                const double new_ref_count = (2.0 * p0 + p1) / denominator;
                const double new_alt_count = (p1 + 2.0 * p2) / denominator;
                const double difference = Kokkos::fmax(
                    Kokkos::fabs(new_ref_count - previous_ref_count),
                    Kokkos::fabs(new_alt_count - previous_alt_count));
                previous_ref_count = new_ref_count;
                previous_alt_count = new_alt_count;
                const double ref_pseudocount = device_priors[candidate * 2];
                const double alt_pseudocount = device_priors[candidate * 2 + 1];
                const double total = ref_pseudocount + alt_pseudocount +
                    new_ref_count + new_alt_count;
                ref_frequency = (ref_pseudocount + new_ref_count) / total;
                alt_frequency = (alt_pseudocount + new_alt_count) / total;
                if (difference <= 0.1) break;
            }

            const double log10_ref = Kokkos::log10(ref_frequency);
            const double log10_alt = Kokkos::log10(alt_frequency);
            const double terms[3] = {
                -0.1 * static_cast<double>(device_pl(base)) + 2.0 * log10_ref,
                -0.1 * static_cast<double>(device_pl(base + 1)) +
                    0.30102999566398119521 + log10_ref + log10_alt,
                -0.1 * static_cast<double>(device_pl(base + 2)) + 2.0 * log10_alt};
            const double maximum = Kokkos::fmax(terms[0], Kokkos::fmax(terms[1], terms[2]));
            const double p0 = Kokkos::pow(10.0, terms[0] - maximum);
            const double denominator = p0 +
                Kokkos::pow(10.0, terms[1] - maximum) +
                Kokkos::pow(10.0, terms[2] - maximum);
            if (p0 > 0.0 && denominator > 0.0) {
                const double p_alt_absent = p0 / denominator;
                device_log10_p_alt_absent(candidate) = Kokkos::fmin(
                    0.0, Kokkos::log10(p_alt_absent));
                device_qual(candidate) = Kokkos::fmax(
                    0.0, -10.0 * Kokkos::log10(p_alt_absent));
                const double p_variant_present = 1.0 - p_alt_absent;
                if (p_variant_present > 0.0)
                    device_log10_p_variant_present(candidate) = Kokkos::fmin(
                        0.0, Kokkos::log10(p_variant_present));
            }
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_qual = Kokkos::create_mirror_view(device_qual);
    auto host_log10_p_alt_absent = Kokkos::create_mirror_view(device_log10_p_alt_absent);
    auto host_log10_p_variant_present =
        Kokkos::create_mirror_view(device_log10_p_variant_present);
    Kokkos::deep_copy(host_qual, device_qual);
    Kokkos::deep_copy(host_log10_p_alt_absent, device_log10_p_alt_absent);
    Kokkos::deep_copy(host_log10_p_variant_present, device_log10_p_variant_present);
    BiallelicCallConfidenceResult result;
    result.qual.resize(candidate_count, 0.0);
    result.log10_p_alt_absent.resize(candidate_count, 0.0);
    result.log10_p_variant_present.resize(candidate_count, 0.0);
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        result.qual[candidate] = host_qual(candidate);
        result.log10_p_alt_absent[candidate] = host_log10_p_alt_absent(candidate);
        result.log10_p_variant_present[candidate] = host_log10_p_variant_present(candidate);
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

GenotypePlResult derive_genotype_gt_gq_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid genotype PL dimensions/ploidy");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || pl.size() != sample_count * width)
        throw std::invalid_argument("PL vector does not match sample/allele dimensions: samples=" +
                                    std::to_string(sample_count) + ", alleles=" +
                                    std::to_string(allele_count) + ", ploidy=" +
                                    std::to_string(ploidy) + ", width=" +
                                    std::to_string(width) + ", values=" +
                                    std::to_string(pl.size()));

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("genotype_pl", pl.size());
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    for (std::size_t index = 0; index < pl.size(); ++index) host_pl(index) = pl[index];
    Kokkos::deep_copy(device_pl, host_pl);
    OutputView alleles("genotype_alleles", sample_count * static_cast<std::size_t>(ploidy));
    OutputView gq("genotype_gq", sample_count);
    Kokkos::deep_copy(alleles, -1);
    Kokkos::deep_copy(gq, -1);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    constexpr std::int32_t kBestSentinel = std::numeric_limits<std::int32_t>::max();
    Kokkos::parallel_for(
        "fastgatk_genotype_pl",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            std::int32_t best = kBestSentinel;
            std::int32_t runner_up = kBestSentinel;
            int best_genotype[32];
            for (int genotype = 0; genotype < static_cast<int>(width); ++genotype) {
                const auto value = device_pl(sample * width + static_cast<std::size_t>(genotype));
                // HTSlib uses INT32_MIN/vector-end sentinels; accepting
                // negative values as missing keeps this kernel independent
                // of HTSlib while preserving all legal non-negative PLs.
                if (value < 0) continue;
                if (value < best) {
                    runner_up = best;
                    best = value;
                    unrank_genotype(static_cast<std::size_t>(genotype), allele_count,
                                    ploidy, best_genotype);
                } else if (value < runner_up) {
                    runner_up = value;
                }
            }
            if (best == kBestSentinel) return;
            for (int copy = 0; copy < ploidy; ++copy)
                alleles[sample * static_cast<std::size_t>(ploidy) + static_cast<std::size_t>(copy)] =
                    best_genotype[copy];
            if (runner_up != kBestSentinel) {
                const auto difference = runner_up - best;
                gq(sample) = difference > 99 ? 99 : difference;
            } else {
                gq(sample) = 0;
            }
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_alleles = Kokkos::create_mirror_view(alleles);
    auto host_gq = Kokkos::create_mirror_view(gq);
    Kokkos::deep_copy(host_alleles, alleles);
    Kokkos::deep_copy(host_gq, gq);

    GenotypePlResult result;
    result.ploidy = ploidy;
    result.alleles.resize(sample_count * static_cast<std::size_t>(ploidy));
    result.gq.resize(sample_count);
    for (std::size_t index = 0; index < result.alleles.size(); ++index)
        result.alleles[index] = host_alleles(index);
    if (ploidy == 2) {
        result.first_allele.resize(sample_count);
        result.second_allele.resize(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            result.first_allele[sample] = result.alleles[sample * 2];
            result.second_allele[sample] = result.alleles[sample * 2 + 1];
        }
    }
    for (std::size_t sample = 0; sample < sample_count; ++sample)
        result.gq[sample] = host_gq(sample);
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

AlleleLikelihoodScoreResult calculate_allele_likelihood_scores_kokkos(
    const std::vector<std::int32_t>& pl,
    const std::size_t sample_count,
    const int allele_count,
    const int ploidy) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid allele-likelihood score dimensions");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || pl.size() != sample_count * width)
        throw std::invalid_argument("PL vector does not match allele-likelihood score dimensions");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("allele_likelihood_score_pl", pl.size());
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    for (std::size_t index = 0; index < pl.size(); ++index)
        host_pl(index) = pl[index];
    Kokkos::deep_copy(device_pl, host_pl);
    // A per-sample matrix avoids atomics and makes the final reduction exactly
    // sample ordered on every execution space.  The matrix is bounded by the
    // same Number=G contract as the caller's PL payload.
    Kokkos::View<double**, MemorySpace> device_sample_scores(
        "allele_likelihood_sample_scores", sample_count, allele_count);
    Kokkos::deep_copy(device_sample_scores, 0.0);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_allele_likelihood_scores",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            std::int32_t best = std::numeric_limits<std::int32_t>::max();
            int best_index = -1;
            const auto offset = sample * width;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(offset + genotype);
                // Missing/vector-end PLs are not likelihoods.  Strictly less
                // preserves MathUtils.maxElementIndex's first-index tie rule
                // after converting GL (-0.1 * PL) to the lowest PL.
                if (value < 0 || value >= best) continue;
                best = value;
                best_index = static_cast<int>(genotype);
            }
            if (best_index < 0) return;
            const auto hom_ref = device_pl(offset);
            if (hom_ref < 0) return;
            const auto difference = (static_cast<double>(hom_ref) -
                                     static_cast<double>(best)) / 10.0;
            if (!(difference >= 0.0) || !Kokkos::isfinite(difference)) return;
            int best_alleles[32];
            unrank_genotype(static_cast<std::size_t>(best_index), allele_count,
                            ploidy, best_alleles);
            for (int copy = 0; copy < ploidy; ++copy) {
                const int allele = best_alleles[copy];
                if (allele > 0) device_sample_scores(sample, allele) = difference;
            }
            // A homozygous or polyploid genotype can contain the same allele
            // more than once.  AlleleSubsettingUtils.containsAllele() counts
            // it once, so retain one score per allele rather than ploidy copies.
            for (int allele = 1; allele < allele_count; ++allele) {
                bool present = false;
                for (int copy = 0; copy < ploidy; ++copy)
                    if (best_alleles[copy] == allele) { present = true; break; }
                if (!present) device_sample_scores(sample, allele) = 0.0;
            }
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_sample_scores = Kokkos::create_mirror_view(device_sample_scores);
    Kokkos::deep_copy(host_sample_scores, device_sample_scores);
    AlleleLikelihoodScoreResult result;
    result.scores.assign(static_cast<std::size_t>(allele_count), 0.0);
    for (std::size_t sample = 0; sample < sample_count; ++sample)
        for (int allele = 1; allele < allele_count; ++allele)
            result.scores[static_cast<std::size_t>(allele)] +=
                host_sample_scores(sample, allele);
    result.sample_count = sample_count;
    result.allele_count = allele_count;
    result.ploidy = ploidy;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

GenotypePlResult derive_genotype_gt_gq_from_log10_priors_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& log10_genotype_priors) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid posterior genotype dimensions/ploidy");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || pl.size() != sample_count * width ||
        log10_genotype_priors.size() != width)
        throw std::invalid_argument("posterior PL/prior dimensions do not match");
    for (const auto value : log10_genotype_priors)
        if (!std::isfinite(value))
            throw std::invalid_argument("posterior genotype priors must be finite");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("posterior_genotype_pl", pl.size());
    PriorView device_priors("posterior_genotype_log10_priors", width);
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    for (std::size_t index = 0; index < pl.size(); ++index) host_pl(index) = pl[index];
    for (std::size_t index = 0; index < width; ++index) host_priors(index) = log10_genotype_priors[index];
    Kokkos::deep_copy(device_pl, host_pl);
    Kokkos::deep_copy(device_priors, host_priors);
    OutputView alleles("posterior_genotype_alleles",
                       sample_count * static_cast<std::size_t>(ploidy));
    OutputView gq("posterior_genotype_gq", sample_count);
    PosteriorView posterior_phred("posterior_genotype_phred",
                                  sample_count * width);
    Kokkos::deep_copy(alleles, -1);
    Kokkos::deep_copy(gq, -1);
    Kokkos::deep_copy(posterior_phred, -1.0);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_genotype_posterior_assignment",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            double maximum = -1.0e300;
            int best_index = -1;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double posterior = device_priors(genotype) -
                    0.1 * static_cast<double>(value);
                if (posterior > maximum) {
                    maximum = posterior;
                    best_index = static_cast<int>(genotype);
                }
            }
            if (best_index < 0) return;
            double denominator = 0.0;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double posterior = device_priors(genotype) -
                    0.1 * static_cast<double>(value);
                denominator += Kokkos::pow(10.0, posterior - maximum);
            }
            if (!(denominator > 0.0)) return;
            const double best_probability = 1.0 / denominator;
            const double error_probability = best_probability >= 1.0
                ? 0.0 : 1.0 - best_probability;
            int genotype_quality = 99;
            if (error_probability > 0.0) {
                genotype_quality = static_cast<int>(
                    -10.0 * Kokkos::log10(error_probability) + 0.5);
                if (genotype_quality > 99) genotype_quality = 99;
                if (genotype_quality < 0) genotype_quality = 0;
            }
            const double log10_denominator = Kokkos::log10(denominator);
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double posterior = device_priors(genotype) -
                    0.1 * static_cast<double>(value);
                const double normalized = posterior - maximum - log10_denominator;
                posterior_phred(sample * width + genotype) =
                    normalized == 0.0 ? 0.0 : -10.0 * normalized;
            }
            int best_genotype[32];
            unrank_genotype(static_cast<std::size_t>(best_index), allele_count,
                            ploidy, best_genotype);
            for (int copy = 0; copy < ploidy; ++copy)
                alleles[sample * static_cast<std::size_t>(ploidy) +
                        static_cast<std::size_t>(copy)] = best_genotype[copy];
            gq(sample) = genotype_quality;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_alleles = Kokkos::create_mirror_view(alleles);
    auto host_gq = Kokkos::create_mirror_view(gq);
    auto host_posterior_phred = Kokkos::create_mirror_view(posterior_phred);
    Kokkos::deep_copy(host_alleles, alleles);
    Kokkos::deep_copy(host_gq, gq);
    Kokkos::deep_copy(host_posterior_phred, posterior_phred);
    GenotypePlResult result;
    result.ploidy = ploidy;
    result.alleles.resize(sample_count * static_cast<std::size_t>(ploidy));
    result.gq.resize(sample_count);
    result.posterior_phred.resize(sample_count * width);
    result.prior_phred.resize(width);
    for (std::size_t index = 0; index < result.alleles.size(); ++index)
        result.alleles[index] = host_alleles(index);
    if (ploidy == 2) {
        result.first_allele.resize(sample_count);
        result.second_allele.resize(sample_count);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            result.first_allele[sample] = result.alleles[sample * 2];
            result.second_allele[sample] = result.alleles[sample * 2 + 1];
        }
    }
    for (std::size_t sample = 0; sample < sample_count; ++sample)
        result.gq[sample] = host_gq(sample);
    for (std::size_t index = 0; index < result.posterior_phred.size(); ++index)
        result.posterior_phred[index] = host_posterior_phred(index);
    for (std::size_t index = 0; index < width; ++index)
        result.prior_phred[index] = log10_genotype_priors[index] == 0.0
            ? 0.0 : -10.0 * log10_genotype_priors[index];
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

GenotypePlResult derive_diploid_gt_gq_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count) {
    return derive_genotype_gt_gq_kokkos(pl, sample_count, allele_count, 2);
}

RemapGenotypePlResult remap_genotype_pl_kokkos(
    const std::vector<std::int32_t>& source_pl,
    std::size_t sample_count,
    int source_allele_count,
    int target_allele_count,
    int ploidy,
    const std::vector<std::int32_t>& kept_alleles) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || source_allele_count < 2 || target_allele_count < 2 ||
        ploidy <= 0 || ploidy > 32 || kept_alleles.size() !=
        static_cast<std::size_t>(target_allele_count))
        throw std::invalid_argument("invalid genotype PL remap dimensions");
    for (const auto allele : kept_alleles)
        if (allele < -1 || allele >= source_allele_count)
            throw std::invalid_argument("genotype PL remap contains an invalid source allele");
    const auto source_width = genotype_count(source_allele_count, ploidy);
    const auto target_width = genotype_count(target_allele_count, ploidy);
    if (source_width == 0 || target_width == 0 || source_pl.size() !=
        sample_count * source_width)
        throw std::invalid_argument("source PL vector does not match remap dimensions");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_source("genotype_pl_remap_source", source_pl.size());
    InputView device_kept("genotype_pl_remap_alleles", kept_alleles.size());
    auto host_source = Kokkos::create_mirror_view(device_source);
    auto host_kept = Kokkos::create_mirror_view(device_kept);
    for (std::size_t index = 0; index < source_pl.size(); ++index)
        host_source(index) = source_pl[index];
    for (std::size_t index = 0; index < kept_alleles.size(); ++index)
        host_kept(index) = kept_alleles[index];
    Kokkos::deep_copy(device_source, host_source);
    Kokkos::deep_copy(device_kept, host_kept);
    OutputView device_target("genotype_pl_remap_target",
                             sample_count * target_width);
    Kokkos::deep_copy(device_target, std::numeric_limits<std::int32_t>::min());
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_genotype_pl_remap",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count * target_width),
        KOKKOS_LAMBDA(const std::size_t cell) {
            const auto sample = cell / target_width;
            const auto target_index = cell % target_width;
            int target_tuple[32];
            int source_tuple[32];
            unrank_genotype(target_index, target_allele_count, ploidy, target_tuple);
            bool invalid = false;
            for (int copy = 0; copy < ploidy; ++copy) {
                source_tuple[copy] = device_kept(target_tuple[copy]);
                if (source_tuple[copy] < 0) invalid = true;
            }
            if (invalid) {
                device_target(cell) = std::numeric_limits<std::int32_t>::min();
                return;
            }
            // A joint ALT union can reorder the source shard's alleles.  The
            // VCF rank is defined on a nondecreasing tuple, so normalize the
            // mapped tuple before computing its source row.
            for (int left = 1; left < ploidy; ++left) {
                const auto value = source_tuple[left];
                int right = left - 1;
                while (right >= 0 && source_tuple[right] > value) {
                    source_tuple[right + 1] = source_tuple[right];
                    --right;
                }
                source_tuple[right + 1] = value;
            }
            // The VCF combination order is the same nondecreasing tuple
            // ordering used by rank_genotype; compute it inline so this
            // remains backend-independent and does not allocate in a kernel.
            std::size_t source_index = 0;
            for (int position = ploidy - 1; position >= 0; --position) {
                const int selected = source_tuple[position];
                for (int candidate = 0; candidate < selected; ++candidate)
                    source_index += genotype_count(candidate + 1, position);
            }
            const auto value = device_source(sample * source_width + source_index);
            device_target(cell) = value;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    auto host_target = Kokkos::create_mirror_view(device_target);
    Kokkos::deep_copy(host_target, device_target);

    RemapGenotypePlResult result;
    result.ploidy = ploidy;
    result.pl.resize(sample_count * target_width);
    for (std::size_t index = 0; index < result.pl.size(); ++index)
        result.pl[index] = host_target(index);
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

RemapAlleleFieldResult remap_allele_field_kokkos(
    const std::vector<std::int32_t>& source_values,
    std::size_t sample_count,
    int source_allele_count,
    int target_allele_count,
    const std::vector<std::int32_t>& target_to_source) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || source_allele_count < 1 || target_allele_count < 1 ||
        target_to_source.size() != static_cast<std::size_t>(target_allele_count))
        throw std::invalid_argument("invalid allele-field remap dimensions");
    for (const auto source : target_to_source)
        if (source < -1 || source >= source_allele_count)
            throw std::invalid_argument("allele-field remap contains an invalid source allele");
    const auto source_size = sample_count * static_cast<std::size_t>(source_allele_count);
    const auto target_size = sample_count * static_cast<std::size_t>(target_allele_count);
    if (source_values.size() != source_size)
        throw std::invalid_argument("source allele field does not match remap dimensions");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_source("allele_field_remap_source", source_size);
    InputView device_map("allele_field_remap_map", target_to_source.size());
    auto host_source = Kokkos::create_mirror_view(device_source);
    auto host_map = Kokkos::create_mirror_view(device_map);
    for (std::size_t index = 0; index < source_size; ++index) host_source(index) = source_values[index];
    for (std::size_t index = 0; index < target_to_source.size(); ++index) host_map(index) = target_to_source[index];
    Kokkos::deep_copy(device_source, host_source);
    Kokkos::deep_copy(device_map, host_map);
    OutputView device_target("allele_field_remap_target", target_size);
    Kokkos::deep_copy(device_target, std::numeric_limits<std::int32_t>::min());
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_allele_field_remap",
        Kokkos::RangePolicy<ExecSpace>(0, target_size),
        KOKKOS_LAMBDA(const std::size_t cell) {
            const auto target_allele = static_cast<std::size_t>(
                cell % static_cast<std::size_t>(target_allele_count));
            const auto source_allele = device_map(target_allele);
            if (source_allele >= 0)
                device_target(cell) = device_source(
                    (cell / static_cast<std::size_t>(target_allele_count)) *
                    static_cast<std::size_t>(source_allele_count) +
                    static_cast<std::size_t>(source_allele));
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    auto host_target = Kokkos::create_mirror_view(device_target);
    Kokkos::deep_copy(host_target, device_target);

    RemapAlleleFieldResult result;
    result.values.resize(target_size);
    for (std::size_t index = 0; index < target_size; ++index) result.values[index] = host_target(index);
    result.prepare_seconds = std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds = std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

AlleleCountResult count_alleles_kokkos(
    const std::vector<std::int32_t>& allele_indices,
    std::size_t sample_count,
    int ploidy,
    int allele_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || ploidy <= 0 || allele_count < 2)
        throw std::invalid_argument("invalid allele-count dimensions");
    const auto expected = sample_count * static_cast<std::size_t>(ploidy);
    if (allele_indices.size() != expected)
        throw std::invalid_argument("allele vector does not match sample/ploidy dimensions");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_alleles("genotype_alleles", allele_indices.size());
    auto host_alleles = Kokkos::create_mirror_view(device_alleles);
    for (std::size_t index = 0; index < allele_indices.size(); ++index)
        host_alleles(index) = allele_indices[index];
    Kokkos::deep_copy(device_alleles, host_alleles);
    OutputView device_counts("genotype_allele_counts",
                            sample_count * static_cast<std::size_t>(allele_count));
    Kokkos::deep_copy(device_counts, 0);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    // One work item owns one (sample, allele) cell and scans only that
    // sample's fixed-ploidy row.  This avoids global atomics and preserves
    // deterministic integer semantics on CPU and accelerator backends.
    const auto cells = sample_count * static_cast<std::size_t>(allele_count);
    Kokkos::parallel_for(
        "fastgatk_genotype_allele_counts",
        Kokkos::RangePolicy<ExecSpace>(0, cells),
        KOKKOS_LAMBDA(const std::size_t cell) {
            const auto sample = cell / static_cast<std::size_t>(allele_count);
            const auto allele = static_cast<int>(cell % static_cast<std::size_t>(allele_count));
            std::int32_t count = 0;
            for (int copy = 0; copy < ploidy; ++copy) {
                const auto value = device_alleles[sample * static_cast<std::size_t>(ploidy) +
                                                  static_cast<std::size_t>(copy)];
                if (value == allele) ++count;
            }
            device_counts(cell) = count;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_counts = Kokkos::create_mirror_view(device_counts);
    Kokkos::deep_copy(host_counts, device_counts);
    AlleleCountResult result;
    result.counts.assign(static_cast<std::size_t>(allele_count), 0);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        for (int allele = 0; allele < allele_count; ++allele)
            result.counts[static_cast<std::size_t>(allele)] +=
                host_counts(sample * static_cast<std::size_t>(allele_count) +
                            static_cast<std::size_t>(allele));
    }
    for (const auto count : result.counts) result.an += count;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

GenotypePosteriorResult calculate_site_posterior_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& log10_priors,
    const int spanning_deletion_index) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid posterior genotype dimensions");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || pl.size() != sample_count * width ||
        log10_priors.size() != width)
        throw std::invalid_argument("posterior PL/prior dimensions do not match");
    if (spanning_deletion_index < -1 || spanning_deletion_index >= allele_count)
        throw std::invalid_argument("invalid posterior spanning-deletion allele index");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("genotype_posterior_pl", pl.size());
    PriorView device_priors("genotype_log10_priors", log10_priors.size());
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    for (std::size_t index = 0; index < pl.size(); ++index) host_pl(index) = pl[index];
    for (std::size_t index = 0; index < log10_priors.size(); ++index)
        host_priors(index) = log10_priors[index];
    Kokkos::deep_copy(device_pl, host_pl);
    Kokkos::deep_copy(device_priors, host_priors);
    PosteriorView device_log10_p0("genotype_log10_p0", sample_count);
    Kokkos::deep_copy(device_log10_p0, -1.0e300);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_genotype_posterior",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            // PL is a relative likelihood, so subtracting the largest finite
            // log10 term makes the normalization stable without changing the
            // posterior. Negative PL values are HTSlib missing/vector-end
            // sentinels and are ignored just as in the GT/GQ kernel.
            double maximum = -1.0e300;
            bool has_value = false;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double term = -0.1 * static_cast<double>(value) +
                                    device_priors(genotype);
                if (term > maximum) maximum = term;
                has_value = true;
            }
            if (!has_value) return;
            double sum = 0.0;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double term = -0.1 * static_cast<double>(value) +
                                    device_priors(genotype);
                sum += Kokkos::pow(10.0, term - maximum);
            }
            if (sum <= 0.0) return;
            double non_variant_scaled = 0.0;
            int genotype_alleles[32];
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                bool non_variant_genotype = genotype == 0;
                if (spanning_deletion_index >= 0) {
                    non_variant_genotype = true;
                    unrank_genotype(genotype, allele_count, ploidy, genotype_alleles);
                    for (int copy = 0; copy < ploidy; ++copy) {
                        const auto allele = genotype_alleles[copy];
                        if (allele != 0 && allele != spanning_deletion_index) {
                            non_variant_genotype = false;
                            break;
                        }
                    }
                }
                if (!non_variant_genotype) continue;
                const double term = -0.1 * static_cast<double>(value) +
                                    device_priors(genotype);
                non_variant_scaled += Kokkos::pow(10.0, term - maximum);
            }
            if (!(non_variant_scaled > 0.0)) return;
            double log10_p0 = Kokkos::log10(non_variant_scaled) - Kokkos::log10(sum);
            // Roundoff can make a probability infinitesimally greater than
            // one.  Keep the quality non-negative and deterministic.
            if (log10_p0 > 0.0) log10_p0 = 0.0;
            device_log10_p0(sample) = log10_p0;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_log10_p0 = Kokkos::create_mirror_view(device_log10_p0);
    Kokkos::deep_copy(host_log10_p0, device_log10_p0);
    GenotypePosteriorResult result;
    result.log10_no_variant = 0.0;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto value = host_log10_p0(sample);
        if (value <= -1.0e299) continue;
        result.log10_no_variant += value;
        ++result.samples_with_likelihoods;
    }
    result.qual = std::max(0.0, -10.0 * result.log10_no_variant);
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

CrossSampleReferenceConfidenceResult
calculate_cross_sample_reference_confidence_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& shared_log10_priors) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid cross-sample posterior dimensions");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || pl.size() != sample_count * width ||
        shared_log10_priors.size() != width)
        throw std::invalid_argument("cross-sample PL/prior dimensions do not match");
    for (const auto prior : shared_log10_priors)
        if (!std::isfinite(prior))
            throw std::invalid_argument("cross-sample genotype priors must be finite");

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("cross_sample_reference_pl", pl.size());
    PriorView device_priors("cross_sample_reference_priors", shared_log10_priors.size());
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    auto host_priors = Kokkos::create_mirror_view(device_priors);
    for (std::size_t index = 0; index < pl.size(); ++index) host_pl(index) = pl[index];
    for (std::size_t index = 0; index < shared_log10_priors.size(); ++index)
        host_priors(index) = shared_log10_priors[index];
    Kokkos::deep_copy(device_pl, host_pl);
    Kokkos::deep_copy(device_priors, host_priors);
    PosteriorView device_sample_log10_p_reference(
        "cross_sample_log10_p_reference", sample_count);
    Kokkos::deep_copy(device_sample_log10_p_reference, -1.0e300);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_cross_sample_reference_confidence",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            double maximum = -1.0e300;
            bool has_value = false;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double term = -0.1 * static_cast<double>(value) +
                                    device_priors(genotype);
                if (term > maximum) maximum = term;
                has_value = true;
            }
            if (!has_value) return;
            double denominator = 0.0;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                const double term = -0.1 * static_cast<double>(value) +
                                    device_priors(genotype);
                denominator += Kokkos::pow(10.0, term - maximum);
            }
            const auto ref_value = device_pl(sample * width);
            if (ref_value < 0 || denominator <= 0.0) return;
            const double ref_term = -0.1 * static_cast<double>(ref_value) +
                                    device_priors(0);
            double log10_probability = ref_term - maximum - Kokkos::log10(denominator);
            if (log10_probability > 0.0) log10_probability = 0.0;
            device_sample_log10_p_reference(sample) = log10_probability;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_sample = Kokkos::create_mirror_view(device_sample_log10_p_reference);
    Kokkos::deep_copy(host_sample, device_sample_log10_p_reference);
    CrossSampleReferenceConfidenceResult result;
    result.sample_log10_p_reference.resize(sample_count, -1.0e300);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto value = host_sample(sample);
        result.sample_log10_p_reference[sample] = value;
        if (value <= -1.0e299) continue;
        result.joint_log10_p_reference += std::min(0.0, value);
        ++result.samples_with_likelihoods;
    }
    if (result.samples_with_likelihoods == 0)
        result.joint_log10_p_reference = 0.0;
    result.joint_qual = std::max(0.0, -10.0 * result.joint_log10_p_reference);
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

AlleleFrequencyResult calculate_allele_frequency_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& prior_pseudocounts,
    const std::vector<std::int32_t>& sample_gq,
    const std::vector<std::int32_t>& sample_alleles,
    const int spanning_deletion_index) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (sample_count == 0 || allele_count < 2 || ploidy <= 0 || ploidy > 32)
        throw std::invalid_argument("invalid cohort genotype dimensions");
    if (spanning_deletion_index < -1 || spanning_deletion_index >= allele_count)
        throw std::invalid_argument("invalid spanning-deletion allele index");
    const auto width = genotype_count(allele_count, ploidy);
    if (width == 0 || (pl.size() != 0 && pl.size() != sample_count * width) ||
        prior_pseudocounts.size() != static_cast<std::size_t>(allele_count))
        throw std::invalid_argument("cohort PL/prior dimensions do not match");
    if ((!sample_gq.empty() && sample_gq.size() != sample_count) ||
        (!sample_alleles.empty() && sample_alleles.size() !=
            sample_count * static_cast<std::size_t>(ploidy)))
        throw std::invalid_argument("cohort fallback GT/GQ dimensions do not match");
    for (const auto prior : prior_pseudocounts)
        if (!(prior > 0.0) || !std::isfinite(prior))
            throw std::invalid_argument("cohort allele pseudocounts must be finite and positive");
    // Very large PL arrays make both GATK and this exact-combination model
    // impractical.  Fail closed instead of allowing a device allocation to
    // become an unbounded resource request.
    if (width > 1'000'000)
        throw std::invalid_argument("cohort genotype PL width exceeds safe limit");

    // GATK accepts a diploid hom-ref genotype with GQ but no PL as an
    // informative AF-calculation input.  Synthesize its PL row before the
    // device transfer so the rest of the EM/posterior path remains a single
    // backend-independent Kokkos implementation.  An all-missing PL row is
    // deliberately required; partially populated rows retain their supplied
    // likelihoods and are not silently rewritten.
    std::vector<std::int32_t> effective_pl = pl;
    if (effective_pl.empty())
        effective_pl.assign(sample_count * width, std::numeric_limits<std::int32_t>::min());
    std::size_t approximate_gq_samples = 0;
    if (ploidy == 2 && sample_gq.size() == sample_count &&
        sample_alleles.size() == sample_count * static_cast<std::size_t>(ploidy)) {
        constexpr int hom_var_scale_factor = 10;
        std::vector<int> genotype_alleles(static_cast<std::size_t>(ploidy), 0);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            bool missing_pl = true;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                if (effective_pl[sample * width + genotype] >= 0) {
                    missing_pl = false;
                    break;
                }
            }
            if (!missing_pl || sample_gq[sample] < 0 ||
                sample_alleles[sample * 2] != 0 ||
                sample_alleles[sample * 2 + 1] != 0)
                continue;
            const auto gq = sample_gq[sample];
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                unrank_genotype(genotype, allele_count, ploidy, genotype_alleles.data());
                const bool hom_ref = genotype == 0;
                const bool contains_ref = genotype_alleles[0] == 0 || genotype_alleles[1] == 0;
                effective_pl[sample * width + genotype] = hom_ref ? 0 :
                    (contains_ref ? gq : hom_var_scale_factor * gq);
            }
            ++approximate_gq_samples;
        }
    }

    const auto prepare_begin = std::chrono::steady_clock::now();
    InputView device_pl("cohort_genotype_pl", effective_pl.size());
    PriorView device_log10_af("cohort_log10_af", static_cast<std::size_t>(allele_count));
    PriorView device_log10_combo("cohort_log10_combo", width);
    CountView device_genotype_counts("cohort_genotype_counts",
                                    width * static_cast<std::size_t>(allele_count));
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    for (std::size_t index = 0; index < effective_pl.size(); ++index)
        host_pl(index) = effective_pl[index];
    Kokkos::deep_copy(device_pl, host_pl);

    std::vector<double> log10_combo(width, 0.0);
    std::vector<std::int32_t> genotype_counts(
        width * static_cast<std::size_t>(allele_count), 0);
    const double log10_e = std::log(10.0);
    std::vector<int> genotype_alleles;
    genotype_alleles.resize(static_cast<std::size_t>(ploidy));
    for (std::size_t genotype = 0; genotype < width; ++genotype) {
        unrank_genotype(genotype, allele_count, ploidy, genotype_alleles.data());
        std::vector<int> counts(static_cast<std::size_t>(allele_count), 0);
        for (const auto allele : genotype_alleles) ++counts[static_cast<std::size_t>(allele)];
        double log10_multinomial = std::lgamma(static_cast<double>(ploidy) + 1.0) / log10_e;
        for (int allele = 0; allele < allele_count; ++allele) {
            const auto count = counts[static_cast<std::size_t>(allele)];
            genotype_counts[genotype * static_cast<std::size_t>(allele_count) +
                            static_cast<std::size_t>(allele)] = count;
            log10_multinomial -= std::lgamma(static_cast<double>(count) + 1.0) / log10_e;
        }
        log10_combo[genotype] = log10_multinomial;
    }
    auto host_log10_combo = Kokkos::create_mirror_view(device_log10_combo);
    auto host_genotype_counts = Kokkos::create_mirror_view(device_genotype_counts);
    for (std::size_t index = 0; index < width; ++index) host_log10_combo(index) = log10_combo[index];
    for (std::size_t index = 0; index < genotype_counts.size(); ++index)
        host_genotype_counts(index) = genotype_counts[index];
    Kokkos::deep_copy(device_log10_combo, host_log10_combo);
    Kokkos::deep_copy(device_genotype_counts, host_genotype_counts);

    const auto prepare_end = std::chrono::steady_clock::now();
    const auto execute_begin = std::chrono::steady_clock::now();

    // Evaluate one EM effective-count step.  One work item owns a sample and
    // writes all allele counts for that sample, so the Host reduction remains
    // deterministic and no device atomics are needed.
    const auto evaluate_counts = [&](const std::vector<double>& log10_af) {
        auto host_af = Kokkos::create_mirror_view(device_log10_af);
        for (int allele = 0; allele < allele_count; ++allele)
            host_af(allele) = log10_af[static_cast<std::size_t>(allele)];
        Kokkos::deep_copy(device_log10_af, host_af);
        PriorView device_sample_counts(
            "cohort_sample_effective_counts",
            sample_count * static_cast<std::size_t>(allele_count));
        Kokkos::deep_copy(device_sample_counts, 0.0);
        Kokkos::parallel_for(
            "fastgatk_cohort_effective_allele_counts",
            Kokkos::RangePolicy<ExecSpace>(0, sample_count),
            KOKKOS_LAMBDA(const std::size_t sample) {
                double maximum = -1.0e300;
                bool has_value = false;
                for (std::size_t genotype = 0; genotype < width; ++genotype) {
                    const auto value = device_pl(sample * width + genotype);
                    if (value < 0) continue;
                    double term = device_log10_combo(genotype) -
                                  0.1 * static_cast<double>(value);
                    for (int allele = 0; allele < allele_count; ++allele)
                        term += static_cast<double>(device_genotype_counts[
                            genotype * static_cast<std::size_t>(allele_count) +
                            static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                    if (term > maximum) maximum = term;
                    has_value = true;
                }
                if (!has_value) return;
                double denominator = 0.0;
                for (std::size_t genotype = 0; genotype < width; ++genotype) {
                    const auto value = device_pl(sample * width + genotype);
                    if (value < 0) continue;
                    double term = device_log10_combo(genotype) -
                                  0.1 * static_cast<double>(value);
                    for (int allele = 0; allele < allele_count; ++allele)
                        term += static_cast<double>(device_genotype_counts[
                            genotype * static_cast<std::size_t>(allele_count) +
                            static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                    denominator += Kokkos::pow(10.0, term - maximum);
                }
                if (denominator <= 0.0) return;
                for (std::size_t genotype = 0; genotype < width; ++genotype) {
                    const auto value = device_pl(sample * width + genotype);
                    if (value < 0) continue;
                    double term = device_log10_combo(genotype) -
                                  0.1 * static_cast<double>(value);
                    for (int allele = 0; allele < allele_count; ++allele)
                        term += static_cast<double>(device_genotype_counts[
                            genotype * static_cast<std::size_t>(allele_count) +
                            static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                    const double probability = Kokkos::pow(10.0, term - maximum) / denominator;
                    for (int allele = 0; allele < allele_count; ++allele)
                        device_sample_counts(sample * static_cast<std::size_t>(allele_count) +
                                             static_cast<std::size_t>(allele)) +=
                            probability * static_cast<double>(device_genotype_counts[
                                genotype * static_cast<std::size_t>(allele_count) +
                                static_cast<std::size_t>(allele)]);
                }
            });
        ExecSpace{}.fence();
        auto host_sample_counts = Kokkos::create_mirror_view(device_sample_counts);
        Kokkos::deep_copy(host_sample_counts, device_sample_counts);
        std::vector<double> counts(static_cast<std::size_t>(allele_count), 0.0);
        for (std::size_t sample = 0; sample < sample_count; ++sample)
            for (int allele = 0; allele < allele_count; ++allele)
                counts[static_cast<std::size_t>(allele)] += host_sample_counts(
                    sample * static_cast<std::size_t>(allele_count) +
                    static_cast<std::size_t>(allele));
        return counts;
    };

    std::vector<double> allele_counts(static_cast<std::size_t>(allele_count), 0.0);
    std::vector<double> log10_af(static_cast<std::size_t>(allele_count),
                                 -std::log10(static_cast<double>(allele_count)));
    constexpr double convergence_threshold = 0.1;
    constexpr int max_iterations = 100;
    int iterations = 0;
    bool converged = false;
    for (; iterations < max_iterations; ++iterations) {
        const auto new_counts = evaluate_counts(log10_af);
        double maximum_difference = 0.0;
        for (int allele = 0; allele < allele_count; ++allele)
            maximum_difference = std::max(maximum_difference,
                std::abs(allele_counts[static_cast<std::size_t>(allele)] -
                         new_counts[static_cast<std::size_t>(allele)]));
        allele_counts = new_counts;
        double sum = 0.0;
        for (int allele = 0; allele < allele_count; ++allele)
            sum += prior_pseudocounts[static_cast<std::size_t>(allele)] +
                   allele_counts[static_cast<std::size_t>(allele)];
        if (!(sum > 0.0)) throw std::runtime_error("cohort Dirichlet sum is non-positive");
        for (int allele = 0; allele < allele_count; ++allele) {
            const auto posterior_pseudocount =
                prior_pseudocounts[static_cast<std::size_t>(allele)] +
                allele_counts[static_cast<std::size_t>(allele)];
            log10_af[static_cast<std::size_t>(allele)] =
                std::log10(posterior_pseudocount) - std::log10(sum);
        }
        if (maximum_difference <= convergence_threshold) {
            converged = true;
            ++iterations;
            break;
        }
    }
    if (!converged) iterations = max_iterations;

    // Evaluate the final normalized genotype posteriors once more.  This
    // supplies site P(no-variant), per-allele P(allele absent), and the
    // deterministic effective AC used for GATK-style annotations.
    auto host_af = Kokkos::create_mirror_view(device_log10_af);
    for (int allele = 0; allele < allele_count; ++allele)
        host_af(allele) = log10_af[static_cast<std::size_t>(allele)];
    Kokkos::deep_copy(device_log10_af, host_af);
    PriorView device_sample_p0("cohort_sample_p0", sample_count);
    PriorView device_sample_absent("cohort_sample_absent",
                                  sample_count * static_cast<std::size_t>(allele_count));
    // The site posteriors are accumulated in LOG space (log-sum-exp) rather than
    // as linear probabilities.  A confident homozygous-alternate sample has
    // P(hom-ref) far below the smallest positive double (PL differences of a few
    // thousand are common), so the linear form underflowed to exactly 0 and its
    // log10 became -inf, which the host then mistook for "this sample has no
    // likelihoods" -- leaving the cohort without any allele-absent information
    // and pruning a perfectly good ALT.  Measured on GATK's own chr20 corpus: 18
    // of 252 emitted loci (all confident hom-alt calls) disappeared that way.
    PriorView device_absent_log_max("cohort_absent_log_max",
                                   sample_count * static_cast<std::size_t>(allele_count));
    PriorView device_absent_log_sum("cohort_absent_log_sum",
                                   sample_count * static_cast<std::size_t>(allele_count));
    Kokkos::deep_copy(device_sample_p0, -1.0e300);
    Kokkos::deep_copy(device_sample_absent, 0.0);
    Kokkos::deep_copy(device_absent_log_max, -1.0e300);
    Kokkos::deep_copy(device_absent_log_sum, 0.0);
    Kokkos::parallel_for(
        "fastgatk_cohort_final_posteriors",
        Kokkos::RangePolicy<ExecSpace>(0, sample_count),
        KOKKOS_LAMBDA(const std::size_t sample) {
            double maximum = -1.0e300;
            bool has_value = false;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                double term = device_log10_combo(genotype) -
                              0.1 * static_cast<double>(value);
                for (int allele = 0; allele < allele_count; ++allele)
                    term += static_cast<double>(device_genotype_counts[
                        genotype * static_cast<std::size_t>(allele_count) +
                        static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                if (term > maximum) maximum = term;
                has_value = true;
            }
            if (!has_value) return;
            double denominator = 0.0;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                double term = device_log10_combo(genotype) -
                              0.1 * static_cast<double>(value);
                for (int allele = 0; allele < allele_count; ++allele)
                    term += static_cast<double>(device_genotype_counts[
                        genotype * static_cast<std::size_t>(allele_count) +
                        static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                denominator += Kokkos::pow(10.0, term - maximum);
            }
            if (denominator <= 0.0) return;
            const double log10_denominator = Kokkos::log10(denominator);
            double p0_log_max = -1.0e300;
            double p0_log_sum = 0.0;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = device_pl(sample * width + genotype);
                if (value < 0) continue;
                double term = device_log10_combo(genotype) -
                              0.1 * static_cast<double>(value);
                for (int allele = 0; allele < allele_count; ++allele)
                    term += static_cast<double>(device_genotype_counts[
                        genotype * static_cast<std::size_t>(allele_count) +
                        static_cast<std::size_t>(allele)]) * device_log10_af(allele);
                const double log10_probability = term - maximum - log10_denominator;
                bool non_variant_genotype = genotype == 0;
                if (spanning_deletion_index >= 0) {
                    non_variant_genotype = true;
                    for (int allele = 1; allele < allele_count; ++allele) {
                        if (allele == spanning_deletion_index) continue;
                        if (device_genotype_counts[
                                genotype * static_cast<std::size_t>(allele_count) +
                                static_cast<std::size_t>(allele)] > 0) {
                            non_variant_genotype = false;
                            break;
                        }
                    }
                }
                if (non_variant_genotype) {
                    if (log10_probability > p0_log_max) {
                        p0_log_sum = p0_log_sum *
                            Kokkos::pow(10.0, p0_log_max - log10_probability) + 1.0;
                        p0_log_max = log10_probability;
                    } else {
                        p0_log_sum += Kokkos::pow(10.0, log10_probability - p0_log_max);
                    }
                }
                for (int allele = 0; allele < allele_count; ++allele) {
                    if (device_genotype_counts[genotype * static_cast<std::size_t>(allele_count) +
                                               static_cast<std::size_t>(allele)] != 0)
                        continue;
                    const auto index = sample * static_cast<std::size_t>(allele_count) +
                                       static_cast<std::size_t>(allele);
                    const auto running_max = device_absent_log_max(index);
                    if (log10_probability > running_max) {
                        device_absent_log_sum(index) = device_absent_log_sum(index) *
                            Kokkos::pow(10.0, running_max - log10_probability) + 1.0;
                        device_absent_log_max(index) = log10_probability;
                    } else {
                        device_absent_log_sum(index) +=
                            Kokkos::pow(10.0, log10_probability - running_max);
                    }
                }
            }
            for (int allele = 0; allele < allele_count; ++allele) {
                const auto index = sample * static_cast<std::size_t>(allele_count) +
                                   static_cast<std::size_t>(allele);
                const auto total = device_absent_log_sum(index);
                device_sample_absent(index) = total > 0.0
                    ? device_absent_log_max(index) + Kokkos::log10(total)
                    : -1.0e300;
            }
            device_sample_p0(sample) = p0_log_sum > 0.0
                ? p0_log_max + Kokkos::log10(p0_log_sum)
                : -1.0e300;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_sample_p0 = Kokkos::create_mirror_view(device_sample_p0);
    auto host_sample_absent = Kokkos::create_mirror_view(device_sample_absent);
    Kokkos::deep_copy(host_sample_p0, device_sample_p0);
    Kokkos::deep_copy(host_sample_absent, device_sample_absent);
    AlleleFrequencyResult result;
    result.effective_allele_counts = allele_counts;
    result.approximate_gq_samples = approximate_gq_samples;
    result.log10_allele_frequencies = log10_af;
    result.integer_allele_counts.resize(static_cast<std::size_t>(allele_count), 0);
    result.log10_p_allele_absent.assign(static_cast<std::size_t>(allele_count), 0.0);
    result.log10_p_no_variant = 0.0;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto p0 = host_sample_p0(sample);
        if (p0 <= -1.0e299) continue;
        ++result.samples_with_likelihoods;
        result.log10_p_no_variant += std::min(0.0, p0);
        for (int allele = 0; allele < allele_count; ++allele) {
            const auto absent = host_sample_absent(
                sample * static_cast<std::size_t>(allele_count) +
                static_cast<std::size_t>(allele));
            result.log10_p_allele_absent[static_cast<std::size_t>(allele)] +=
                std::min(0.0, absent);
        }
    }
    for (int allele = 0; allele < allele_count; ++allele) {
        const auto rounded = std::llround(allele_counts[static_cast<std::size_t>(allele)]);
        result.integer_allele_counts[static_cast<std::size_t>(allele)] =
            static_cast<std::int32_t>(std::clamp<long long>(rounded, 0, std::numeric_limits<std::int32_t>::max()));
    }
    result.qual = std::max(0.0, -10.0 * result.log10_p_no_variant);
    result.iterations = iterations;
    result.converged = converged;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

AlleleFrequencyResult calculate_allele_frequency_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& prior_pseudocounts) {
    return calculate_allele_frequency_kokkos(
        pl, sample_count, allele_count, ploidy, prior_pseudocounts, {}, {});
}

}  // namespace fastgatk::kernels
