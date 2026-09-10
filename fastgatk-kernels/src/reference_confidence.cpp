#include "fastgatk/kernels/reference_confidence.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;
using IndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
using QualityView = Kokkos::View<std::uint8_t*, MemorySpace>;
using FlagView = Kokkos::View<std::uint8_t*, MemorySpace>;
using LikelihoodView = Kokkos::View<double*[3], MemorySpace>;
using ScalarLikelihoodView = Kokkos::View<double*, MemorySpace>;

KOKKOS_INLINE_FUNCTION double log10_sum10(const double first, const double second) {
    const auto high = first > second ? first : second;
    if (!Kokkos::isfinite(high)) return high;
    const auto scaled = Kokkos::pow(10.0, first - high) +
                        Kokkos::pow(10.0, second - high);
    return high + Kokkos::log10(Kokkos::fmax(scaled,
                                             std::numeric_limits<double>::min()));
}

KOKKOS_INLINE_FUNCTION double reference_genotype_mix_score(
    const double reference_likelihood,
    const double non_ref_likelihood,
    const int ploidy,
    const int non_ref_copies) {
    if (non_ref_copies == 0) return reference_likelihood;
    if (non_ref_copies == ploidy) return non_ref_likelihood;
    const auto reference_term = reference_likelihood + Kokkos::log10(
        static_cast<double>(ploidy - non_ref_copies));
    const auto non_ref_term = non_ref_likelihood + Kokkos::log10(
        static_cast<double>(non_ref_copies));
    const auto high = reference_term > non_ref_term ? reference_term : non_ref_term;
    if (!Kokkos::isfinite(high)) return high;
    // ReferenceConfidenceModel deliberately calls
    // MathUtils.approximateLog10SumLog10() here, rather than the exact
    // log-sum.  Its Jacobian table quantizes the difference to 1e-4 up to
    // eight log10 units; reproducing that in the Kokkos kernel is required
    // for high-ploidy Number=G reference blocks (p8 exposes the difference).
    const auto low = reference_term > non_ref_term ? non_ref_term : reference_term;
    const auto difference = high - low;
    double correction = 0.0;
    if (difference < 8.0) {
        const auto table_index = static_cast<int>(difference * 10000.0 + 0.5);
        const auto quantized_difference = static_cast<double>(table_index) * 1.0e-4;
        correction = Kokkos::log10(1.0 + Kokkos::pow(10.0, -quantized_difference));
    }
    return high + correction - Kokkos::log10(static_cast<double>(ploidy));
}

}  // namespace

ReferenceConfidenceResult calculate_reference_confidence_kokkos(
    const std::vector<std::uint32_t>& locus_indices,
    const std::vector<std::uint8_t>& qualities,
    const std::vector<std::uint8_t>& is_alt,
    const std::size_t locus_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (locus_indices.size() != qualities.size() ||
        locus_indices.size() != is_alt.size())
        throw std::invalid_argument("reference-confidence observation arrays have different lengths");
    if (locus_count == 0 && !locus_indices.empty())
        throw std::invalid_argument("reference-confidence locus count is zero");
    for (const auto locus : locus_indices)
        if (locus >= locus_count)
            throw std::invalid_argument("reference-confidence locus index is out of range");
    for (const auto flag : is_alt)
        if (flag > 1)
            throw std::invalid_argument("reference-confidence is_alt must be 0 or 1");

    const auto prepare_begin = std::chrono::steady_clock::now();
    IndexView device_locus("rcm_locus", locus_indices.size());
    QualityView device_quality("rcm_quality", qualities.size());
    FlagView device_alt("rcm_alt", is_alt.size());
    auto host_locus = Kokkos::create_mirror_view(device_locus);
    auto host_quality = Kokkos::create_mirror_view(device_quality);
    auto host_alt = Kokkos::create_mirror_view(device_alt);
    for (std::size_t index = 0; index < locus_indices.size(); ++index) {
        host_locus(index) = locus_indices[index];
        host_quality(index) = qualities[index];
        host_alt(index) = is_alt[index];
    }
    Kokkos::deep_copy(device_locus, host_locus);
    Kokkos::deep_copy(device_quality, host_quality);
    Kokkos::deep_copy(device_alt, host_alt);
    LikelihoodView device_likelihoods("rcm_likelihoods", locus_indices.size());
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_reference_confidence",
        Kokkos::RangePolicy<ExecSpace>(0, locus_indices.size()),
        KOKKOS_LAMBDA(const std::size_t observation) {
            const auto quality = device_quality(observation);
            const auto error = Kokkos::pow(10.0, -0.1 * static_cast<double>(quality));
            const auto correct = Kokkos::fmax(1.0 - error, 1e-300);
            const auto incorrect = Kokkos::fmax(error / 3.0, 1e-300);
            const bool alt = device_alt(observation) != 0;
            const auto reference_likelihood = Kokkos::log10(alt ? incorrect : correct);
            const auto non_ref_likelihood = Kokkos::log10(alt ? correct : incorrect);
            device_likelihoods(observation, 0) = reference_likelihood;
            device_likelihoods(observation, 1) =
                log10_sum10(reference_likelihood, non_ref_likelihood) - Kokkos::log10(2.0);
            device_likelihoods(observation, 2) = non_ref_likelihood;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_likelihoods = Kokkos::create_mirror_view(device_likelihoods);
    Kokkos::deep_copy(host_likelihoods, device_likelihoods);
    ReferenceConfidenceResult result;
    result.reference_likelihoods.resize(locus_indices.size());
    result.non_ref_likelihoods.resize(locus_indices.size());
    result.hom_ref.assign(locus_count, 0.0);
    result.het.assign(locus_count, 0.0);
    result.hom_alt.assign(locus_count, 0.0);
    result.depth.assign(locus_count, 0);
    result.reference_count.assign(locus_count, 0);
    result.non_ref_count.assign(locus_count, 0);
    // Host reduction in observation order is intentional: it makes the
    // floating-point result reproducible when the same input is dispatched to
    // OpenMP, CUDA, HIP or another Kokkos execution space.
    for (std::size_t observation = 0; observation < locus_indices.size(); ++observation) {
        const auto locus = static_cast<std::size_t>(host_locus(observation));
        result.reference_likelihoods[observation] = host_likelihoods(observation, 0);
        result.non_ref_likelihoods[observation] = host_likelihoods(observation, 2);
        result.hom_ref[locus] += host_likelihoods(observation, 0);
        result.het[locus] += host_likelihoods(observation, 1);
        result.hom_alt[locus] += host_likelihoods(observation, 2);
        ++result.depth[locus];
        if (host_alt(observation) != 0) ++result.non_ref_count[locus];
        else ++result.reference_count[locus];
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

ReferenceConfidenceGenotypeResult calculate_reference_confidence_genotypes_kokkos(
    const std::vector<std::uint32_t>& locus_indices,
    const std::vector<double>& reference_likelihoods,
    const std::vector<double>& non_ref_likelihoods,
    const std::size_t locus_count,
    const int ploidy) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (locus_indices.size() != reference_likelihoods.size() ||
        locus_indices.size() != non_ref_likelihoods.size())
        throw std::invalid_argument("reference-confidence genotype arrays have different lengths");
    if (ploidy <= 0 || ploidy > 32 || (locus_count == 0 && !locus_indices.empty()))
        throw std::invalid_argument("invalid reference-confidence genotype dimensions");
    for (const auto locus : locus_indices)
        if (locus >= locus_count)
            throw std::invalid_argument("reference-confidence genotype locus is out of range");
    for (std::size_t index = 0; index < reference_likelihoods.size(); ++index)
        if (!std::isfinite(reference_likelihoods[index]) ||
            !std::isfinite(non_ref_likelihoods[index]))
            throw std::invalid_argument("reference-confidence genotype likelihoods must be finite");

    const auto width = static_cast<std::size_t>(ploidy + 1);
    const auto total = locus_count * width;
    const auto prepare_begin = std::chrono::steady_clock::now();

    // The previous implementation launched one output item per (locus,
    // genotype), then scanned the complete observation array and discarded
    // all observations belonging to other loci.  That is O(loci * G * N)
    // work for a sparse AssemblyRegion.  Build stable per-locus segments once
    // on Host instead.  Stable sorting preserves the original observation
    // order inside each locus, so the strict floating-point reduction has the
    // same summation order and therefore the same PL/GQ bits on every backend.
    std::vector<std::size_t> order(locus_indices.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](const auto left, const auto right) {
        return locus_indices[left] < locus_indices[right];
    });
    std::vector<double> sorted_reference(reference_likelihoods.size());
    std::vector<double> sorted_non_ref(non_ref_likelihoods.size());
    std::vector<std::uint64_t> locus_offsets(locus_count + 1, 0);
    for (std::size_t sorted = 0; sorted < order.size(); ++sorted) {
        const auto original = order[sorted];
        const auto locus = static_cast<std::size_t>(locus_indices[original]);
        sorted_reference[sorted] = reference_likelihoods[original];
        sorted_non_ref[sorted] = non_ref_likelihoods[original];
        ++locus_offsets[locus + 1];
    }
    for (std::size_t locus = 1; locus < locus_offsets.size(); ++locus)
        locus_offsets[locus] += locus_offsets[locus - 1];

    using OffsetView = Kokkos::View<std::uint64_t*, MemorySpace>;
    ScalarLikelihoodView device_reference("rcm_genotype_reference", sorted_reference.size());
    ScalarLikelihoodView device_non_ref("rcm_genotype_non_ref", sorted_non_ref.size());
    OffsetView device_offsets("rcm_genotype_locus_offsets", locus_offsets.size());
    Kokkos::View<double*, MemorySpace> device_scores("rcm_genotype_scores", total);
    Kokkos::View<std::int32_t*, MemorySpace> device_pl("rcm_genotype_pl", total);
    Kokkos::View<std::uint8_t*, MemorySpace> device_gq("rcm_genotype_gq", locus_count);
    Kokkos::View<double*, MemorySpace> device_hom_ref_gq(
        "rcm_genotype_hom_ref_gq", locus_count);
    auto host_reference = Kokkos::create_mirror_view(device_reference);
    auto host_non_ref = Kokkos::create_mirror_view(device_non_ref);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    for (std::size_t index = 0; index < order.size(); ++index) {
        host_reference(index) = sorted_reference[index];
        host_non_ref(index) = sorted_non_ref[index];
    }
    for (std::size_t locus = 0; locus < locus_offsets.size(); ++locus)
        host_offsets(locus) = locus_offsets[locus];
    Kokkos::deep_copy(device_reference, host_reference);
    Kokkos::deep_copy(device_non_ref, host_non_ref);
    Kokkos::deep_copy(device_offsets, host_offsets);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_reference_confidence_genotype_scores",
        Kokkos::RangePolicy<ExecSpace>(0, total),
        KOKKOS_LAMBDA(const std::size_t flat) {
            const auto locus = flat / width;
            const auto non_ref_copies = static_cast<int>(flat % width);
            double score = 0.0;
            const auto begin = device_offsets(locus);
            const auto end = device_offsets(locus + 1);
            for (std::uint64_t observation = begin; observation < end; ++observation) {
                score += reference_genotype_mix_score(
                    device_reference(observation), device_non_ref(observation),
                    ploidy, non_ref_copies);
            }
            device_scores(flat) = score;
        });
    ExecSpace{}.fence();
    Kokkos::parallel_for(
        "fastgatk_reference_confidence_genotype_pl",
        Kokkos::RangePolicy<ExecSpace>(0, locus_count),
        KOKKOS_LAMBDA(const std::size_t locus) {
            const auto ref_score = device_scores[locus * width];
            int gq = 99;
            double raw_hom_ref_gq = 1.0e300;
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                auto score = device_scores[locus * width + genotype];
                if (genotype != 0) score = Kokkos::fmin(score, ref_score);
                const auto raw = genotype == 0 ? 0.0 : Kokkos::fmax(
                    0.0, -10.0 * (score - ref_score));
                // HTSJDK's GenotypeLikelihoods does not apply the 999 cap
                // used by the PairHMM candidate-envelope writer.  GVCF
                // ReferenceConfidenceModel therefore preserves large PLs
                // (for example 1080 on a high-depth hom-ref block); only the
                // signed integer representation bounds the materialized
                // value.  Keeping this distinction is required for direct
                // BP_RESOLUTION/GVCF replacement.
                const auto rounded = Kokkos::fmin(
                    2147483647.0, Kokkos::fmax(0.0, Kokkos::round(raw)));
                device_pl[locus * width + genotype] = static_cast<std::int32_t>(rounded);
                if (genotype != 0) {
                    gq = gq < static_cast<int>(rounded) ? gq : static_cast<int>(rounded);
                    raw_hom_ref_gq = Kokkos::fmin(raw_hom_ref_gq, raw);
                }
            }
            device_gq[locus] = static_cast<std::uint8_t>(gq > 99 ? 99 : gq);
            device_hom_ref_gq(locus) = raw_hom_ref_gq < 1.0e299
                ? raw_hom_ref_gq : 0.0;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();
    auto host_pl = Kokkos::create_mirror_view(device_pl);
    auto host_gq = Kokkos::create_mirror_view(device_gq);
    auto host_hom_ref_gq = Kokkos::create_mirror_view(device_hom_ref_gq);
    Kokkos::deep_copy(host_pl, device_pl);
    Kokkos::deep_copy(host_gq, device_gq);
    Kokkos::deep_copy(host_hom_ref_gq, device_hom_ref_gq);
    ReferenceConfidenceGenotypeResult result;
    result.genotype_pl.resize(total, 0);
    result.gq.resize(locus_count, 0);
    result.hom_ref_gq_phred.resize(locus_count, 0.0);
    for (std::size_t index = 0; index < total; ++index) result.genotype_pl[index] = host_pl(index);
    for (std::size_t locus = 0; locus < locus_count; ++locus) {
        result.gq[locus] = host_gq(locus);
        result.hom_ref_gq_phred[locus] = host_hom_ref_gq(locus);
    }
    result.ploidy = ploidy;
    result.locus_count = locus_count;
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

}  // namespace fastgatk::kernels
