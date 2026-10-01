#include "fastgatk/kernels/gpu_safety.hpp"
#include "fastgatk/kernels/activity_profile.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

using U32View = Kokkos::View<std::uint32_t*, MemorySpace>;
using U8View = Kokkos::View<std::uint8_t*, MemorySpace>;
using DView = Kokkos::View<double*, MemorySpace>;

constexpr double kMinProbabilityToKeepInFilter = 1.0e-5;
constexpr std::int32_t kMaxFilterSize = 50;
constexpr double kDefaultSigma = 17.0;

KOKKOS_INLINE_FUNCTION double activity_digamma_positive(const double input) {
    if (!(input > 0.0)) return -std::numeric_limits<double>::max();
    if (input <= 1.0e-5)
        return -0.5772156649015329 - 1.0 / input;
    double x = input;
    while (x < 49.0) x += 1.0;
    double result = Kokkos::log(x) - 0.5 / x;
    const double inverse = 1.0 / x;
    const double inverse2 = inverse * inverse;
    // Apache Commons Math Gamma.digamma (the overload used by
    // Mutect2Engine.logLikelihoodRatio) retains the 1/120 Bernoulli term.
    result -= inverse2 * (1.0 / 12.0 + inverse2 *
                          (1.0 / 120.0 - inverse2 / 252.0));
    while (x - 1.0 >= input) {
        x -= 1.0;
        result -= 1.0 / x;
    }
    return result;
}

KOKKOS_INLINE_FUNCTION double activity_fast_bernoulli_entropy(const double p) {
    const auto product = p * (1.0 - p);
    return product * (11.0 + 33.0 * product) / (2.0 + 20.0 * product);
}

KOKKOS_INLINE_FUNCTION std::uint64_t activity_gcd(
    std::uint64_t left, std::uint64_t right) {
    while (right != 0U) {
        const auto remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

// Commons Math's CombinatoricsUtils.binomialCoefficientLog uses exact
// integer coefficients for the small counts encountered in a pileup, a
// rounded multiplicative coefficient for the medium range, and a log-product
// path after double overflow.  Keep the branch structure in the Kokkos
// activity kernel instead of replacing it with a backend-specific lgamma.
KOKKOS_INLINE_FUNCTION double activity_binomial_log(
    const std::uint64_t n, std::uint64_t k) {
    if (k == 0U || k == n) return 0.0;
    if (k == 1U || k == n - 1U) return Kokkos::log(static_cast<double>(n));
    if (k > n / 2U) k = n - k;
    if (n < 67U) {
        std::uint64_t coefficient = 1U;
        std::uint64_t factor = n - k + 1U;
        for (std::uint64_t divisor = 1U; divisor <= k; ++divisor, ++factor) {
            const auto divisor_gcd = activity_gcd(factor, divisor);
            coefficient = (coefficient / (divisor / divisor_gcd)) *
                (factor / divisor_gcd);
        }
        return Kokkos::log(static_cast<double>(coefficient));
    }
    if (n < 1030U) {
        double coefficient = 1.0;
        for (std::uint64_t i = 1U; i <= k; ++i)
            coefficient *= static_cast<double>(n - k + i) /
                           static_cast<double>(i);
        coefficient = Kokkos::floor(coefficient + 0.5);
        return Kokkos::log(coefficient);
    }
    double log_sum = 0.0;
    for (std::uint64_t i = n - k + 1U; i <= n; ++i)
        log_sum += Kokkos::log(static_cast<double>(i));
    for (std::uint64_t i = 2U; i <= k; ++i)
        log_sum -= Kokkos::log(static_cast<double>(i));
    return log_sum;
}

KOKKOS_INLINE_FUNCTION double activity_log_likelihood_ratio(
    const U32View& quality_offsets,
    const U8View& quality_bases,
    const U8View& quality_values,
    const U8View& quality_alt_flags,
    const U32View& indel_quality_offsets,
    const U8View& indel_quality_values,
    const DView& digamma_values,
    const DView& binomial_log_values,
    const DView& quality_epsilon_values,
    const DView& quality_log_epsilon_values,
    const DView& quality_log_one_minus_epsilon_values,
    const std::size_t binomial_stride,
    const std::size_t locus,
    const std::uint8_t reference_base,
    const std::uint8_t pcr_snv_quality,
    const std::uint8_t quality_correction,
    std::uint32_t& total_observations) {
    // Overlapping mates are capped by the shared Host quality-correction
    // stage before these observations reach the activity kernel.  Keep the
    // parameter in the API so a future direct activity caller can supply the
    // same GATK pcr-snv-qual policy without changing the ABI.
    (void)pcr_snv_quality;
    const auto begin = static_cast<std::size_t>(quality_offsets(locus));
    const auto end = static_cast<std::size_t>(quality_offsets(locus + 1));
    const auto indel_begin = static_cast<std::size_t>(indel_quality_offsets(locus));
    const auto indel_end = static_cast<std::size_t>(indel_quality_offsets(locus + 1));
    total_observations = static_cast<std::uint32_t>(end - begin);
    if (reference_base >= 4 || end <= begin) return 0.0;

    std::uint64_t best_alt_count = 0;
    std::uint64_t best_alt_qual_sum = 0;
    std::uint8_t best_alt = reference_base;
    bool best_alt_is_indel = false;
    for (std::uint8_t base = 0; base < 4; ++base) {
        if (base == reference_base) continue;
        std::uint64_t count = 0;
        std::uint64_t qual_sum = 0;
        for (std::size_t index = begin; index < end; ++index) {
            // PileupQualBuffer places an aligned base next to a useful
            // soft clip in its indel bucket and does not subsequently add it
            // as a substitution.  The Host marks that mutually exclusive
            // classification in quality_alt_flags for the somatic path.
            if (quality_alt_flags(index) != 0U) continue;
            if (quality_bases(index) != base || quality_values(index) <= 6) continue;
            ++count;
            const auto corrected = Kokkos::min(
                static_cast<int>(quality_values(index)) +
                    static_cast<int>(quality_correction), 93);
            qual_sum += static_cast<std::uint64_t>(corrected);
        }
        // PileupQualBuffer.likeliestIndexAndQuals breaks ties by the first
        // base index (strictly greater quality sum).
        if (qual_sum > best_alt_qual_sum) {
            best_alt_qual_sum = qual_sum;
            best_alt_count = count;
            best_alt = base;
        }
    }
    std::uint64_t indel_qual_sum = 0;
    for (std::size_t index = indel_begin; index < indel_end; ++index)
        indel_qual_sum += static_cast<std::uint64_t>(indel_quality_values(index));
    if (indel_qual_sum > best_alt_qual_sum) {
        best_alt_qual_sum = indel_qual_sum;
        best_alt_count = indel_end - indel_begin;
        best_alt = 4U;
        best_alt_is_indel = true;
    }
    (void)best_alt;
    // `PileupQualBuffer.likeliestIndexAndQuals()` returns an empty ALT list
    // when every pileup element supports the reference (or its alternate
    // base has quality <= 6).  Mutect2 still calls logLikelihoodRatio with
    // nRef = pileup.size() and nAlt = 0 in that case.  For the default flat
    // Beta prior its beta-entropy term is -log(nRef + 1), which is strictly
    // negative for a non-empty pileup.  Returning zero here subtly differs
    // at --initial-tumor-lod 0: every reference-only locus becomes active,
    // joins into a large AssemblyRegion, and can manufacture graph calls.
    // Preserve the source calculation before the ALT-only early exit.
    if (best_alt_count == 0U)
        return -Kokkos::log(static_cast<double>(total_observations) + 1.0);

    const auto n_ref = static_cast<std::uint64_t>(total_observations) -
                       best_alt_count;
    const auto n_alt = best_alt_count;
    const auto f_tilde_ratio = Kokkos::exp(
        digamma_values(static_cast<std::size_t>(n_ref + 1U)) -
        digamma_values(static_cast<std::size_t>(n_alt + 1U)));
    double read_sum = 0.0;
    const auto accumulate_qual = [&](const int quality) {
        const auto quality_index = static_cast<std::size_t>(quality);
        const auto epsilon = quality_epsilon_values(quality_index);
        const auto z_bar_alt = (1.0 - epsilon) /
            (1.0 - epsilon + epsilon * f_tilde_ratio);
        const auto log_epsilon = quality_log_epsilon_values(quality_index);
        const auto log_one_minus_epsilon =
            quality_log_one_minus_epsilon_values(quality_index);
        read_sum += z_bar_alt * (log_one_minus_epsilon - log_epsilon) +
                    activity_fast_bernoulli_entropy(z_bar_alt);
    };
    if (best_alt_is_indel) {
        for (std::size_t index = indel_begin; index < indel_end; ++index)
            accumulate_qual(static_cast<int>(indel_quality_values(index)));
    } else {
        for (std::size_t index = begin; index < end; ++index) {
            if (quality_alt_flags(index) != 0U) continue;
            if (quality_bases(index) != best_alt || quality_values(index) <= 6) continue;
            const auto corrected = Kokkos::min(
                static_cast<int>(quality_values(index)) +
                    static_cast<int>(quality_correction), 93);
            accumulate_qual(corrected);
        }
    }
    const auto n = n_ref + n_alt;
    const auto beta_entropy = -Kokkos::log(static_cast<double>(n) + 1.0) -
                              binomial_log_values(static_cast<std::size_t>(n) * binomial_stride +
                                                  static_cast<std::size_t>(n_alt));
    return beta_entropy + read_sum;
}

KOKKOS_INLINE_FUNCTION double activity_log10_sum(const double left,
                                                  const double right) {
    const auto maximum = left > right ? left : right;
    const auto minimum = left > right ? right : left;
    if (!Kokkos::isfinite(maximum)) return maximum;
    const auto difference = maximum - minimum;
    if (difference >= 8.0) return maximum;
    // GATK's MathUtils.approximateLog10SumLog10 looks up the Jacobian term
    // after MathUtils.fastRound(difference * 10000).  Preserve the 1e-4
    // cache quantization instead of replacing it with an exact log1p: the
    // latter can alter a band-pass boundary by one reference base near the
    // activity threshold.
    const auto rounded_difference = Kokkos::floor(difference * 10000.0 + 0.5) / 10000.0;
    return maximum + Kokkos::log(1.0 + Kokkos::pow(10.0, -rounded_difference)) /
        Kokkos::log(10.0);
}

// Default-ploidy numerical port of HaplotypeCaller's isActive() fast path:
// ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny followed by
// AlleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior.
// `quality_alt_flags` carries the CIGAR-adjacent half of
// ReferenceConfidenceModel.isAltBeforeAssembly from the Host representation.
// Active-region discovery uses at least diploid ploidy in GATK; the current
// Host API supplies that default two-copy case.  The ordinary caller's
// arbitrary-ploidy genotyper remains a later stage and is deliberately not
// conflated with this activity decision.
KOKKOS_INLINE_FUNCTION double activity_hc_ref_vs_any_probability(
    const U32View& quality_offsets,
    const U8View& quality_bases,
    const U8View& quality_values,
    const U8View& quality_alt_flags,
    const std::size_t locus,
    const std::uint8_t reference_base,
    const std::uint8_t min_base_quality,
    const double snp_heterozygosity,
    const double heterozygosity_stdev) {
    constexpr double log10_three = 0.47712125471966243730;
    constexpr double log10_two = 0.30102999566398119521;
    if (reference_base >= 4 || !(snp_heterozygosity > 0.0) ||
        !(heterozygosity_stdev > 0.0))
        return 0.0;
    const auto begin = static_cast<std::size_t>(quality_offsets(locus));
    const auto end = static_cast<std::size_t>(quality_offsets(locus + 1));
    if (end <= begin) return 0.0;

    double hom_ref = 0.0;
    double het = 0.0;
    double hom_alt = 0.0;
    std::size_t usable = 0;
    for (std::size_t index = begin; index < end; ++index) {
        const auto quality = quality_values(index);
        // GATK's ReferenceConfidenceModel skips qualities <= its HC
        // minBaseQualityScore.  This differs deliberately from PairHMM's
        // input floor and from the legacy pileup candidate threshold.
        if (quality <= min_base_quality) continue;
        const double log10_error = -static_cast<double>(quality) / 10.0;
        const double error = Kokkos::pow(10.0, log10_error);
        const double log10_correct = Kokkos::log1p(-error) / Kokkos::log(10.0);
        const double log10_wrong = log10_error - log10_three;
        const bool alternate = quality_bases(index) != reference_base ||
            quality_alt_flags(index) != 0;
        const double reference_likelihood = alternate ? log10_wrong : log10_correct;
        const double non_ref_likelihood = alternate ? log10_correct : log10_wrong;
        hom_ref += reference_likelihood;
        hom_alt += non_ref_likelihood;
        // RefVsAny adds log10(1) for both diploid allele counts then divides
        // every genotype by ploidy after the pileup.  Express that directly
        // to retain the Java operation's normalized likelihood.
        het += activity_log10_sum(reference_likelihood, non_ref_likelihood) - log10_two;
        ++usable;
    }
    if (usable == 0) return 0.0;

    // `returnZeroIfRefIsMax` first examines raw genotype likelihoods.  GATK
    // resolves ties toward index zero, so use strict comparisons here.
    if (!(het > hom_ref) && !(hom_alt > hom_ref)) return 0.0;

    // AlleleFrequencyCalculator's diploid Dirichlet prior can be written as
    // stable ratios to the hom-ref term, avoiding a device lgamma dependency:
    // ref pseudocount = h / s^2 and SNP pseudocount = h * ref.
    const double reference_pseudocount =
        snp_heterozygosity / (heterozygosity_stdev * heterozygosity_stdev);
    const double snp_pseudocount = snp_heterozygosity * reference_pseudocount;
    if (!(reference_pseudocount > 0.0) || !(snp_pseudocount > 0.0)) return 0.0;
    const double het_prior_delta = log10_two + Kokkos::log(snp_pseudocount) /
        Kokkos::log(10.0) - Kokkos::log(reference_pseudocount + 1.0) / Kokkos::log(10.0);
    const double hom_alt_prior_delta =
        (Kokkos::log(snp_pseudocount) + Kokkos::log(snp_pseudocount + 1.0) -
         Kokkos::log(reference_pseudocount) - Kokkos::log(reference_pseudocount + 1.0)) /
        Kokkos::log(10.0);
    const double posterior_ref = hom_ref;
    const double posterior_het = het + het_prior_delta;
    const double posterior_alt = hom_alt + hom_alt_prior_delta;
    if (!(posterior_het > posterior_ref) && !(posterior_alt > posterior_ref)) return 0.0;
    const double maximum = posterior_ref > posterior_het
        ? (posterior_ref > posterior_alt ? posterior_ref : posterior_alt)
        : (posterior_het > posterior_alt ? posterior_het : posterior_alt);
    const double ref_weight = Kokkos::pow(10.0, posterior_ref - maximum);
    const double het_weight = Kokkos::pow(10.0, posterior_het - maximum);
    const double alt_weight = Kokkos::pow(10.0, posterior_alt - maximum);
    const double total = ref_weight + het_weight + alt_weight;
    return total > 0.0 ? 1.0 - ref_weight / total : 0.0;
}

double normal_distribution(const double mean, const double sigma, const double x) {
    if (sigma < 0.0) throw std::invalid_argument("bandpass sigma must be non-negative");
    if (sigma == 0.0) return x == mean ? 1.0 : 0.0;
    // Keep the operation order of MathUtils.normalDistribution: the host
    // libm exp is intentional because GATK calls Math.exp rather than a
    // strict-math implementation.
    constexpr double two_pi = 2.0 * 3.141592653589793238462643383279502884;
    return std::exp(-(x - mean) * (x - mean) / (2.0 * sigma * sigma)) /
           (sigma * std::sqrt(two_pi));
}

std::vector<double> make_bandpass_kernel(const std::int32_t filter_size,
                                         const double sigma) {
    if (filter_size < 0) throw std::invalid_argument("bandpass filter size must be non-negative");
    const auto width = static_cast<std::size_t>(2 * filter_size + 1);
    std::vector<double> kernel(width, 0.0);
    for (std::size_t i = 0; i < width; ++i)
        kernel[i] = normal_distribution(static_cast<double>(filter_size), sigma,
                                         static_cast<double>(i));
    double sum = 0.0;
    for (const auto value : kernel) sum += value;
    if (!(sum > 0.0) || !std::isfinite(sum))
        throw std::invalid_argument("bandpass kernel normalization failed");
    for (auto& value : kernel) value /= sum;
    return kernel;
}

std::int32_t determine_filter_size(const std::vector<double>& kernel,
                                   const double minimum_probability) {
    if (kernel.empty() || (kernel.size() % 2) == 0)
        throw std::invalid_argument("bandpass kernel must have odd non-zero width");
    const auto middle = static_cast<std::int32_t>((kernel.size() - 1) / 2);
    auto filter_end = middle;
    while (filter_end > 0 && kernel[static_cast<std::size_t>(filter_end - 1)] >=
                                  minimum_probability)
        --filter_end;
    return middle - filter_end;
}

struct ProfileState {
    std::int32_t start = 0;
    double probability = 0.0;
};

struct PoppedProfileRegion {
    std::int32_t start = 0;
    std::int32_t end = 0;
    bool active = false;
};

// Host-side port of ActivityProfile + BandPassActivityProfile.  The raw
// activity calculation above is device portable; this object owns the
// variable-length, add-and-accumulate state list whose exact boundaries feed
// graph and PairHMM windows.
class BandPassProfile {
public:
    BandPassProfile(const std::int32_t max_propagation,
                    const double active_threshold,
                    const std::int32_t max_filter_size,
                    const double sigma,
                    const bool adaptive)
        : active_threshold_(active_threshold),
          kernel_(make_bandpass_kernel(max_filter_size, sigma)) {
        filter_size_ = adaptive
            ? determine_filter_size(kernel_, kMinProbabilityToKeepInFilter)
            : max_filter_size;
        if (filter_size_ != max_filter_size)
            kernel_ = make_bandpass_kernel(filter_size_, sigma);
        effective_max_propagation_ = max_propagation + filter_size_;
    }

    std::int32_t filter_size() const { return filter_size_; }
    std::int32_t effective_max_propagation() const { return effective_max_propagation_; }
    bool empty() const { return states_.empty(); }
    std::int32_t end() const { return raw_stop_.value_or(0); }

    void add(const std::int32_t start, const double probability) {
        if (!raw_stop_.has_value()) {
            raw_start_ = start;
            raw_stop_ = start;
        } else {
            if (*raw_stop_ != start - 1)
                throw std::logic_error("activity profile loci must be contiguous");
            raw_stop_ = start;
        }
        if (probability <= 0.0 || filter_size_ == 0) {
            incorporate(start, probability);
            return;
        }
        for (std::int32_t offset = -filter_size_; offset <= filter_size_; ++offset) {
            const auto position = start + offset;
            if (position < 0) continue;
            const auto kernel_index = static_cast<std::size_t>(offset + filter_size_);
            incorporate(position, probability * kernel_[kernel_index]);
        }
    }

    std::vector<PoppedProfileRegion> pop_ready_regions(const std::int32_t min_region_size,
                                                        const std::int32_t max_region_size,
                                                        const bool force) {
        std::vector<PoppedProfileRegion> result;
        while (!states_.empty()) {
            auto next = pop_next(min_region_size, max_region_size, force);
            if (!next.has_value()) break;
            result.push_back(*next);
        }
        return result;
    }

private:
    std::optional<PoppedProfileRegion> pop_next(const std::int32_t min_region_size,
                                                const std::int32_t max_region_size,
                                                const bool force) {
        if (states_.empty()) return std::nullopt;
        if (force && raw_start_.has_value() && raw_stop_.has_value()) {
            const auto span = static_cast<std::int64_t>(*raw_stop_) - *raw_start_ + 1;
            if (span >= 0 && static_cast<std::size_t>(span) < states_.size())
                states_.resize(static_cast<std::size_t>(span));
        }
        const auto is_active = states_.front().probability > active_threshold_;
        if (!force && states_.size() <
            static_cast<std::size_t>(max_region_size + effective_max_propagation_))
            return std::nullopt;
        auto end = find_first_activity_boundary(is_active, max_region_size);
        if (is_active && end == static_cast<std::size_t>(max_region_size))
            end = find_best_cut_site(end, min_region_size);
        if (end == 0) return std::nullopt;
        const auto first_start = states_.front().start;
        const auto last_index = end - 1;
        const auto last_start = states_[last_index].start;
        // ActivityProfile drains one bounded AssemblyRegion at a time.  A
        // vector erase at the front shifts every unconsumed state; for an
        // unbounded 1 Mb traversal of mostly inactive loci that turns the
        // source walk into O(n^2) Host copies (one <=300 bp chunk at a
        // time).  The profile state is logically a FIFO, so pop exactly the
        // same prefix from a deque and retain GATK's boundary decisions.
        for (std::size_t index = 0; index < end; ++index)
            states_.pop_front();
        if (states_.empty()) {
            raw_start_.reset();
            raw_stop_.reset();
        } else {
            // GATK keeps the raw stop for the contiguous-add/end check, but
            // advances the profile start to the first retained filtered state.
            raw_start_ = states_.front().start;
        }
        return PoppedProfileRegion{first_start, last_start, is_active};
    }

    std::size_t find_first_activity_boundary(const bool is_active,
                                             const std::int32_t max_region_size) const {
        std::size_t end = 0;
        while (end < states_.size() && end < static_cast<std::size_t>(max_region_size)) {
            const bool state_active = states_[end].probability > active_threshold_;
            if (state_active != is_active) break;
            ++end;
        }
        return end;
    }

    bool is_minimum(const std::size_t index) const {
        if (index == 0 || index + 1 >= states_.size()) return false;
        const auto here = states_[index].probability;
        return here <= states_[index + 1].probability &&
               here < states_[index - 1].probability;
    }

    std::size_t find_best_cut_site(const std::size_t end_of_active_region,
                                   const std::int32_t min_region_size) const {
        auto minimum_index = end_of_active_region - 1;
        double minimum_probability = std::numeric_limits<double>::max();
        for (std::int64_t index = static_cast<std::int64_t>(end_of_active_region) - 1;
             index >= static_cast<std::int64_t>(min_region_size) - 1; --index) {
            const auto current = states_[static_cast<std::size_t>(index)].probability;
            if (current < minimum_probability && is_minimum(static_cast<std::size_t>(index))) {
                minimum_probability = current;
                minimum_index = static_cast<std::size_t>(index);
            }
        }
        return minimum_index + 1;
    }

    void incorporate(const std::int32_t start, const double probability) {
        if (!raw_start_.has_value()) throw std::logic_error("activity profile start is unset");
        const auto relative = static_cast<std::int64_t>(start) - *raw_start_;
        if (relative < 0) return;
        const auto position = static_cast<std::size_t>(relative);
        if (position < states_.size()) {
            states_[position].probability += probability;
            return;
        }
        if (position != states_.size())
            throw std::logic_error("activity profile filtered states are not contiguous");
        states_.push_back(ProfileState{start, probability});
    }

    double active_threshold_ = 0.0;
    std::vector<double> kernel_;
    std::int32_t filter_size_ = 0;
    std::int32_t effective_max_propagation_ = 0;
    std::optional<std::int32_t> raw_start_;
    std::optional<std::int32_t> raw_stop_;
    std::deque<ProfileState> states_;
};

}  // namespace

ActivityProfileResult compute_activity_profile_kokkos(
    const ActivityProfileInput& input, ActivityProfileOptions options) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (input.tids.size() != input.positions.size() ||
        input.counts.size() != input.tids.size() * 4)
        throw std::invalid_argument("activity profile arrays do not match");
    if (!input.reference_bases.empty() &&
        input.reference_bases.size() != input.tids.size())
        throw std::invalid_argument("activity profile reference array does not match loci");
    if (!input.indel_counts.empty() &&
        input.indel_counts.size() != input.tids.size())
        throw std::invalid_argument("activity profile indel array does not match loci");
    const bool high_quality_softclip_arrays_empty = input.high_quality_softclip_bases.empty() &&
        input.high_quality_softclip_events.empty();
    const bool high_quality_softclip_arrays_complete =
        input.high_quality_softclip_bases.size() == input.tids.size() &&
        input.high_quality_softclip_events.size() == input.tids.size();
    if (!high_quality_softclip_arrays_empty && !high_quality_softclip_arrays_complete)
        throw std::invalid_argument("activity profile high-quality soft-clip arrays do not match loci");
    if (options.min_depth == 0 || options.halo == 0 || options.max_region_size == 0 ||
        options.max_prob_propagation_distance == 0 || options.min_region_size == 0 ||
        options.min_region_size > options.max_region_size ||
        options.bandpass_max_filter_size == 0 ||
        options.bandpass_max_filter_size > static_cast<std::uint32_t>(kMaxFilterSize))
        throw std::invalid_argument("activity profile thresholds must be positive");
    if (options.min_activity < 0.0 || options.min_activity > 1.0)
        throw std::invalid_argument("activity profile min_activity must be in [0,1]");
    if (!std::isfinite(options.bandpass_sigma) || options.bandpass_sigma < 0.0)
        throw std::invalid_argument("activity profile bandpass sigma must be finite and non-negative");
    for (const auto& span : options.traversal_spans) {
        if (span.tid < 0 || span.start < 0 || span.end <= span.start)
            throw std::invalid_argument("activity profile traversal spans must be non-empty and non-negative");
    }
    const bool quality_arrays_empty = input.quality_offsets.empty() &&
        input.quality_bases.empty() && input.quality_values.empty() &&
        input.quality_alt_flags.empty();
    const bool quality_arrays_complete = input.quality_offsets.size() == input.tids.size() + 1 &&
        input.quality_bases.size() == input.quality_values.size() &&
        (input.quality_alt_flags.empty() ||
         input.quality_bases.size() == input.quality_alt_flags.size()) &&
        !input.quality_offsets.empty() && input.quality_offsets.front() == 0 &&
        input.quality_offsets.back() == input.quality_bases.size();
    if (!quality_arrays_empty && !quality_arrays_complete)
        throw std::invalid_argument("activity profile quality arrays do not match loci");
    if (quality_arrays_complete) {
        for (std::size_t i = 1; i < input.quality_offsets.size(); ++i)
            if (input.quality_offsets[i] < input.quality_offsets[i - 1])
                throw std::invalid_argument("activity profile quality offsets are not monotonic");
    }
    const bool indel_quality_arrays_empty = input.indel_quality_offsets.empty() &&
        input.indel_quality_values.empty();
    const bool indel_quality_arrays_complete = indel_quality_arrays_empty ||
        (input.indel_quality_offsets.size() == input.tids.size() + 1 &&
         !input.indel_quality_offsets.empty() && input.indel_quality_offsets.front() == 0 &&
         input.indel_quality_offsets.back() == input.indel_quality_values.size());
    if (!indel_quality_arrays_empty && !indel_quality_arrays_complete)
        throw std::invalid_argument("activity profile indel quality arrays do not match loci");
    if (indel_quality_arrays_complete) {
        for (std::size_t i = 1; i < input.indel_quality_offsets.size(); ++i)
            if (input.indel_quality_offsets[i] < input.indel_quality_offsets[i - 1])
                throw std::invalid_argument("activity profile indel quality offsets are not monotonic");
    }
    if (!input.somatic_normal_suppressed.empty() &&
        input.somatic_normal_suppressed.size() != input.tids.size())
        throw std::invalid_argument("activity profile normal suppression mask does not match loci");
    if (!input.somatic_feature_suppressed.empty() &&
        input.somatic_feature_suppressed.size() != input.tids.size())
        throw std::invalid_argument("activity profile feature suppression mask does not match loci");
    if (!input.forced_allele_active.empty() &&
        input.forced_allele_active.size() != input.tids.size())
        throw std::invalid_argument("activity profile forced-allele mask does not match loci");

    // The quality-aware somatic likelihood repeatedly evaluates the same
    // pure functions for the small set of pileup depths and Phred scores in
    // a batch.  Build those values in a Kokkos pre-kernel below so the main
    // numerical kernel retains GATK's formula and operation order per locus
    // without spending most of a high-depth tile in transcendental calls.
    std::size_t maximum_quality_depth = 0;
    if (quality_arrays_complete) {
        for (std::size_t i = 0; i < input.tids.size(); ++i) {
            maximum_quality_depth = std::max(
                maximum_quality_depth,
                static_cast<std::size_t>(input.quality_offsets[i + 1] -
                                         input.quality_offsets[i]));
        }
    }
    std::size_t maximum_quality_score = 93;
    for (const auto quality : input.indel_quality_values)
        maximum_quality_score = std::max(
            maximum_quality_score, static_cast<std::size_t>(quality));
    constexpr std::size_t kActivityQualityTableMinimumSize = 94;
    const auto activity_quality_table_size = std::max(
        kActivityQualityTableMinimumSize, maximum_quality_score + 1U);
    const auto digamma_table_size = maximum_quality_depth + 2U;
    const auto binomial_stride = maximum_quality_depth + 1U;

    ActivityProfileResult result;
    result.used = true;
    result.reference_aware = !input.reference_bases.empty();
    result.loci = input.tids.size();
    result.execution_space = ExecSpace::name();
    const auto configured_filter_size = static_cast<std::int32_t>(
        options.bandpass_max_filter_size);
    BandPassProfile profile(
        static_cast<std::int32_t>(options.max_prob_propagation_distance),
        options.min_activity, configured_filter_size, options.bandpass_sigma,
        options.bandpass_adaptive_filter);
    result.filter_size = static_cast<std::uint32_t>(profile.filter_size());
    result.effective_max_prob_propagation_distance = static_cast<std::uint32_t>(
        profile.effective_max_propagation());
    result.active.assign(result.loci, 0);
    result.activity.assign(result.loci, 0.0);
    // A selected traversal interval can legitimately have no reads at all.
    // AssemblyRegionIterator still supplies every such empty pileup to the
    // source ActivityProfile (and --force-active then processes the emitted
    // zero-activity regions).  Do not return early when traversal metadata is
    // present: the Host profile below can represent that case without staging
    // synthetic rows on Kokkos.
    if (result.loci == 0 && options.traversal_spans.empty()) return result;

    fastgatk::core::HostBatch host("activity-profile-v1");
    host.records = result.loci;
    host.bytes = input.tids.size() * sizeof(std::int32_t) * 2 +
                 input.counts.size() * sizeof(std::uint32_t) +
                 input.reference_bases.size() * sizeof(std::uint8_t) +
                 input.indel_counts.size() * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("activity-profile");
    plan.begin_prepare(host);
    Kokkos::View<std::uint32_t*, MemorySpace> depth("activity_depth", result.loci);
    Kokkos::View<std::uint32_t*, MemorySpace> max_count("activity_max_count", result.loci);
    Kokkos::View<std::uint32_t*, MemorySpace> counts("activity_counts", result.loci * 4);
    Kokkos::View<std::uint8_t*, MemorySpace> reference("activity_reference", result.loci);
    U32View indel_counts("activity_indel_counts", result.loci);
    Kokkos::View<double*, MemorySpace> activity("activity_signal", result.loci);
    Kokkos::View<std::uint8_t*, MemorySpace> active("activity_mask", result.loci);
    U32View quality_offsets("activity_quality_offsets", input.quality_offsets.size());
    U8View quality_bases("activity_quality_bases", input.quality_bases.size());
    U8View quality_values("activity_quality_values", input.quality_values.size());
    U8View quality_alt_flags("activity_quality_alt_flags", input.quality_bases.size());
    U32View indel_quality_offsets("activity_indel_quality_offsets",
                                  indel_quality_arrays_empty ? result.loci + 1 :
                                      input.indel_quality_offsets.size());
    U8View indel_quality_values("activity_indel_quality_values",
                                input.indel_quality_values.size());
    U8View normal_suppressed("activity_normal_suppressed", result.loci);
    U8View feature_suppressed("activity_feature_suppressed", result.loci);
    U8View forced_allele_active("activity_forced_allele_active", result.loci);
    DView digamma_values("activity_digamma_values", digamma_table_size);
    DView binomial_log_values("activity_binomial_log_values",
                              binomial_stride * binomial_stride);
    DView quality_epsilon_values("activity_quality_epsilon_values",
                                 activity_quality_table_size);
    DView quality_log_epsilon_values("activity_quality_log_epsilon_values",
                                     activity_quality_table_size);
    DView quality_log_one_minus_epsilon_values(
        "activity_quality_log_one_minus_epsilon_values", activity_quality_table_size);
    auto host_depth = Kokkos::create_mirror_view(depth);
    auto host_max = Kokkos::create_mirror_view(max_count);
    auto host_counts = Kokkos::create_mirror_view(counts);
    auto host_reference = Kokkos::create_mirror_view(reference);
    auto host_indel_counts = Kokkos::create_mirror_view(indel_counts);
    auto host_quality_offsets = Kokkos::create_mirror_view(quality_offsets);
    auto host_quality_bases = Kokkos::create_mirror_view(quality_bases);
    auto host_quality_values = Kokkos::create_mirror_view(quality_values);
    auto host_quality_alt_flags = Kokkos::create_mirror_view(quality_alt_flags);
    auto host_indel_quality_offsets = Kokkos::create_mirror_view(indel_quality_offsets);
    auto host_indel_quality_values = Kokkos::create_mirror_view(indel_quality_values);
    auto host_normal_suppressed = Kokkos::create_mirror_view(normal_suppressed);
    auto host_feature_suppressed = Kokkos::create_mirror_view(feature_suppressed);
    auto host_forced_allele_active = Kokkos::create_mirror_view(forced_allele_active);
    for (std::size_t i = 0; i < result.loci; ++i) {
        std::uint32_t total = 0;
        std::uint32_t maximum = 0;
        for (std::size_t base = 0; base < 4; ++base) {
            const auto count = input.counts[i * 4 + base];
            total += count;
            maximum = std::max(maximum, count);
            host_counts(i * 4 + base) = count;
        }
        host_depth(i) = total;
        host_max(i) = maximum;
        host_reference(i) = input.reference_bases.empty() ? 4 : input.reference_bases[i];
        host_indel_counts(i) = input.indel_counts.empty() ? 0U : input.indel_counts[i];
        host_normal_suppressed(i) = input.somatic_normal_suppressed.empty() ? 0U :
            input.somatic_normal_suppressed[i];
        host_feature_suppressed(i) = input.somatic_feature_suppressed.empty() ? 0U :
            input.somatic_feature_suppressed[i];
        host_forced_allele_active(i) = input.forced_allele_active.empty() ? 0U :
            input.forced_allele_active[i];
    }
    Kokkos::deep_copy(depth, host_depth);
    Kokkos::deep_copy(max_count, host_max);
    Kokkos::deep_copy(counts, host_counts);
    Kokkos::deep_copy(reference, host_reference);
    Kokkos::deep_copy(indel_counts, host_indel_counts);
    for (std::size_t i = 0; i < input.quality_offsets.size(); ++i)
        host_quality_offsets(i) = input.quality_offsets[i];
    for (std::size_t i = 0; i < input.quality_bases.size(); ++i) {
        host_quality_bases(i) = input.quality_bases[i];
        host_quality_values(i) = input.quality_values[i];
        host_quality_alt_flags(i) = input.quality_alt_flags.empty() ? 0U :
            input.quality_alt_flags[i];
    }
    for (std::size_t i = 0; i < indel_quality_offsets.extent(0); ++i)
        host_indel_quality_offsets(i) = indel_quality_arrays_empty ? 0U :
            input.indel_quality_offsets[i];
    for (std::size_t i = 0; i < input.indel_quality_values.size(); ++i)
        host_indel_quality_values(i) = input.indel_quality_values[i];
    Kokkos::deep_copy(quality_offsets, host_quality_offsets);
    Kokkos::deep_copy(quality_bases, host_quality_bases);
    Kokkos::deep_copy(quality_values, host_quality_values);
    Kokkos::deep_copy(quality_alt_flags, host_quality_alt_flags);
    Kokkos::deep_copy(indel_quality_offsets, host_indel_quality_offsets);
    Kokkos::deep_copy(indel_quality_values, host_indel_quality_values);
    Kokkos::deep_copy(normal_suppressed, host_normal_suppressed);
    Kokkos::deep_copy(feature_suppressed, host_feature_suppressed);
    Kokkos::deep_copy(forced_allele_active, host_forced_allele_active);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(result.loci);
    device.bind("depth", depth);
    device.bind("max_count", max_count);
    device.bind("counts", counts);
    device.bind("reference", reference);
    device.bind("indel_counts", indel_counts);
    device.bind("quality_offsets", quality_offsets);
    device.bind("quality_bases", quality_bases);
    device.bind("quality_values", quality_values);
    device.bind("quality_alt_flags", quality_alt_flags);
    device.bind("indel_quality_offsets", indel_quality_offsets);
    device.bind("indel_quality_values", indel_quality_values);
    device.bind("normal_suppressed", normal_suppressed);
    device.bind("feature_suppressed", feature_suppressed);
    device.bind("forced_allele_active", forced_allele_active);
    device.bind("digamma_values", digamma_values);
    device.bind("binomial_log_values", binomial_log_values);
    device.bind("quality_epsilon_values", quality_epsilon_values);
    device.bind("quality_log_epsilon_values", quality_log_epsilon_values);
    device.bind("quality_log_one_minus_epsilon_values",
                quality_log_one_minus_epsilon_values);
    device.bind("activity", activity);
    device.bind("active", active);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    const auto min_depth = options.min_depth;
    const auto min_activity = options.min_activity;
    const auto quality_aware_somatic = options.quality_aware_somatic && quality_arrays_complete &&
        indel_quality_arrays_complete;
    const auto initial_log10_odds = options.initial_tumor_log10_odds;
    const auto pcr_snv_quality = options.pcr_snv_quality;
    const auto quality_correction = options.multiple_substitution_quality_correction;
    const auto genotype_aware_hc = options.genotype_aware_hc && quality_arrays_complete;
    const auto hc_min_base_quality = options.hc_min_base_quality;
    const auto hc_snp_heterozygosity = options.hc_snp_heterozygosity;
    const auto hc_heterozygosity_stdev = options.hc_heterozygosity_stdev;
    const auto lookup_work = std::max({
        digamma_table_size, binomial_log_values.extent(0),
        activity_quality_table_size});
    Kokkos::parallel_for("activity_profile_lookup",
        Kokkos::RangePolicy<ExecSpace>(0, lookup_work),
        KOKKOS_LAMBDA(const std::size_t index) {
            if (index < digamma_table_size) {
                digamma_values(index) = activity_digamma_positive(
                    static_cast<double>(index));
            }
            if (index < binomial_log_values.extent(0)) {
                const auto n = index / binomial_stride;
                const auto k = index % binomial_stride;
                binomial_log_values(index) = k <= n
                    ? activity_binomial_log(static_cast<std::uint64_t>(n),
                                            static_cast<std::uint64_t>(k))
                    : 0.0;
            }
            if (index < activity_quality_table_size) {
                const auto quality = static_cast<double>(index);
                const auto log_epsilon = quality * (-Kokkos::log(10.0) / 10.0);
                quality_epsilon_values(index) = Kokkos::pow(10.0, -quality / 10.0);
                quality_log_epsilon_values(index) = log_epsilon;
                quality_log_one_minus_epsilon_values(index) =
                    log_epsilon < -0.69314718055994530942
                        ? Kokkos::log(1.0 - Kokkos::exp(log_epsilon))
                        : Kokkos::log(-Kokkos::expm1(log_epsilon));
            }
        });
    // Mutect2's quality evidence is highly sparse: a long low-coverage
    // prefix can sit beside a short high-depth interval.  The per-locus LOD
    // work is therefore not uniform (the high-depth locus scans its complete
    // PileupQualBuffer several times).  Static RangePolicy partitioning can
    // leave one OpenMP worker with every expensive locus while the other
    // Kokkos workers wait.  Dynamic scheduling changes no per-locus numeric
    // operation or output ordering, but keeps the Kokkos numerical kernel
    // balanced on real high-coverage assembly tiles.
    using ActivityRangePolicy = Kokkos::RangePolicy<
        ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>>;
    Kokkos::parallel_for("activity_profile", ActivityRangePolicy(0, result.loci),
        KOKKOS_LAMBDA(const std::size_t i) {
            const auto total = depth(i);
            const auto ref = reference(i);
            // With a known reference, count every non-reference observation.
            // For an unknown/N reference retain the old diversity signal so
            // no-reference smoke inputs do not silently change behavior.
            double signal = total == 0 ? 0.0 :
                (ref < 4
                    ? 1.0 - static_cast<double>(counts(i * 4 + ref)) /
                        static_cast<double>(total)
                    : 1.0 - static_cast<double>(max_count(i)) /
                        static_cast<double>(total));
            bool is_active = total >= min_depth && signal >= min_activity;
            // Retain the compact indel seed only for the legacy count-only
            // activity API.  HaplotypeCaller instead supplies its complete
            // ReferenceConfidenceModel pileup (including Q30 synthetic
            // deletion elements) through the quality arrays below; forcing
            // every I/D anchor to one would diverge from GATK's posterior in
            // deep reference-supporting pileups.
            if (!genotype_aware_hc && !quality_aware_somatic && indel_counts(i) > 0) {
                signal = 1.0;
                is_active = true;
            }
            if (quality_aware_somatic && ref < 4) {
                std::uint32_t quality_total = 0;
                const auto lod = activity_log_likelihood_ratio(
                    quality_offsets, quality_bases, quality_values, quality_alt_flags,
                    indel_quality_offsets, indel_quality_values,
                    digamma_values, binomial_log_values, quality_epsilon_values,
                    quality_log_epsilon_values, quality_log_one_minus_epsilon_values,
                    binomial_stride, i, ref,
                    pcr_snv_quality, quality_correction, quality_total);
                // Mutect2Engine.isActive() does not apply the caller's
                // callable/min-depth threshold here.  It returns an active
                // state solely from the pileup log-likelihood ratio (an
                // empty pileup naturally has no ALT qualities and therefore
                // stays inactive).  `min_depth` belongs to the count-only HC
                // activity contract; reusing it for this branch silently
                // drops low-depth, high-quality somatic seeds and changes
                // AssemblyRegion ownership.
                (void)quality_total;
                is_active = lod >= initial_log10_odds * Kokkos::log(10.0) &&
                    normal_suppressed(i) == 0U && feature_suppressed(i) == 0U;
                signal = is_active ? 1.0 : 0.0;
            } else if (genotype_aware_hc && ref < 4) {
                signal = activity_hc_ref_vs_any_probability(
                    quality_offsets, quality_bases, quality_values, quality_alt_flags, i, ref,
                    hc_min_base_quality, hc_snp_heterozygosity,
                    hc_heterozygosity_stdev);
                is_active = signal > min_activity;
            }
            // HaplotypeCallerEngine.isActive() checks a supplied --alleles
            // feature before the normal Ref-vs-Any calculation.  Preserve
            // that exact priority: a force-call-filtered eligible event is
            // active even with zero ALT pileup evidence, while it does not
            // alter the downstream PairHMM likelihood observations.
            if (forced_allele_active(i) != 0U) {
                signal = 1.0;
                is_active = true;
            }
            activity(i) = signal;
            active(i) = is_active ? 1 : 0;
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_activity = Kokkos::create_mirror_view(activity);
    auto host_active = Kokkos::create_mirror_view(active);
    Kokkos::deep_copy(host_activity, activity);
    Kokkos::deep_copy(host_active, active);
    for (std::size_t i = 0; i < result.loci; ++i) {
        result.activity[i] = host_activity(i);
        result.active[i] = host_active(i);
        if (result.active[i]) ++result.active_loci;
    }
    std::vector<std::size_t> order(result.loci);
    for (std::size_t i = 0; i < result.loci; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](const auto left, const auto right) {
        if (input.tids[left] != input.tids[right]) return input.tids[left] < input.tids[right];
        return input.positions[left] < input.positions[right];
    });
    const auto halo = static_cast<std::int32_t>(options.halo);
    const auto min_region_size = static_cast<std::int32_t>(options.min_region_size);
    const auto max_region_size = static_cast<std::int32_t>(options.max_region_size);
    const auto gap_limit = static_cast<std::int64_t>(max_region_size) +
                           profile.effective_max_propagation();

    // GATK's ActivityProfile.processState replaces a high-quality-soft-clip
    // state with copies over +/- the truncated running average of its
    // high-quality soft-clipped bases.  BandPassActivityProfile then applies
    // every copy's Gaussian at the original `justAddedState` coordinate (not
    // the copied coordinate), so the net effect is an exact central-weight
    // multiplier of 2*span+1.  Preserve that Java control-flow detail here;
    // moving the Gaussian centres would falsely extend the profile tails.
    std::map<std::pair<std::int32_t, std::int32_t>, double> profile_probabilities;
    for (const auto index : order) {
        const auto tid = input.tids[index];
        const auto position = input.positions[index];
        if (tid < 0 || position < 0) continue;
        std::int32_t softclip_span = 0;
        if (high_quality_softclip_arrays_complete &&
            input.high_quality_softclip_events[index] != 0U) {
            const double average = static_cast<double>(input.high_quality_softclip_bases[index]) /
                static_cast<double>(input.high_quality_softclip_events[index]);
            // HaplotypeCallerEngine's threshold is strictly greater than six
            // bases and ActivityProfile truncates the resulting Double via
            // intValue before applying its propagation cap.
            if (average > 6.0) {
                const auto truncated = static_cast<std::int32_t>(average);
                softclip_span = std::min(truncated, profile.effective_max_propagation());
            }
        }
        if (softclip_span == 0) {
            profile_probabilities[{tid, position}] += result.activity[index];
            continue;
        }
        profile_probabilities[{tid, position}] += result.activity[index] *
            static_cast<double>(2 * softclip_span + 1);
    }
    const auto profile_probability = [&](const std::int32_t tid, const std::int32_t position) {
        const auto found = profile_probabilities.find({tid, position});
        return found == profile_probabilities.end() ? 0.0 : found->second;
    };

    // Active raw loci are retained only for assigning the exact active span
    // to a popped region.  The profile's region span itself may include
    // Gaussian tails and is the bounded graph/PairHMM window.
    std::vector<std::pair<std::int32_t, std::int32_t>> raw_active;
    raw_active.reserve(result.active_loci);
    for (const auto index : order)
        if (result.active[index]) raw_active.emplace_back(input.tids[index], input.positions[index]);

    auto append_popped = [&](const std::int32_t tid,
                             const PoppedProfileRegion& popped) {
        // AssemblyRegion's active span is the complete band-pass segment,
        // not merely the sparse raw loci whose unfiltered PileupQualBuffer
        // score crossed the threshold.  The latter remain counted below for
        // telemetry, while EventMap ownership and the padded graph window
        // must use the former exactly as AssemblyRegion(activeSpan, padding).
        const std::int32_t active_start = popped.start;
        const std::int32_t active_end = popped.end;
        std::size_t active_count = 0;
        for (const auto& raw : raw_active) {
            if (raw.first != tid) continue;
            if (raw.second < popped.start) continue;
            if (raw.second > popped.end) break;
            ++active_count;
        }
        // Mutect2's --force-active does not remove ActivityProfile segment
        // boundaries; it only forces each segment through the active-region
        // walker.  Preserve raw coordinates/counts for the debug track while
        // exposing the effective active state to both the writer and caller.
        const auto effective_active = popped.active || options.force_active;
        result.profile_regions.push_back(ActivityRegion{
            tid, popped.start, popped.end, active_start, active_end,
            active_count, effective_active});
        if (!effective_active) return;  // calling graph materializes active entries only.
        const auto region_start = std::max<std::int32_t>(0, popped.start - halo);
        result.regions.push_back(ActivityRegion{
            tid, region_start, popped.end + halo, active_start, active_end,
            active_count, true});
    };
    auto flush_profile = [&](const std::int32_t tid, const bool force) {
        for (const auto& popped : profile.pop_ready_regions(min_region_size, max_region_size, force))
            append_popped(tid, popped);
    };
    const auto append_zero_profile_regions = [&](const std::int32_t tid,
                                                 std::int64_t begin,
                                                 const std::int64_t end) {
        // A constant zero ActivityProfile has no activity boundary, so
        // ActivityProfile.findFirstActivityBoundary() cuts it only at the
        // hard max-region boundary.  Materialize those final Host regions
        // directly rather than iterating each uncovered base.  This is
        // required for source-equivalent --force-active traversal, and keeps
        // an empty chromosome interval out of the Kokkos evidence payload.
        while (begin < end) {
            const auto inclusive_end = std::min<std::int64_t>(
                end - 1, begin + static_cast<std::int64_t>(max_region_size) - 1);
            append_popped(tid, PoppedProfileRegion{
                static_cast<std::int32_t>(begin),
                static_cast<std::int32_t>(inclusive_end), false});
            begin = inclusive_end + 1;
        }
    };

    if (options.traversal_spans.empty()) {
        // The compact API historically has no traversal-domain metadata. Keep
        // that sparse behavior for callers that do not provide it.
        std::int32_t current_tid = -1;
        std::int32_t last_position = -1;
        bool have_locus = false;
        for (const auto index : order) {
            const auto tid = input.tids[index];
            const auto position = input.positions[index];
            if (position < 0) continue;
            if (!have_locus || tid != current_tid) {
                if (have_locus) flush_profile(current_tid, true);
                current_tid = tid;
                last_position = position - 1;
                have_locus = true;
            }
            const auto gap = static_cast<std::int64_t>(position) - last_position - 1;
            if (gap > gap_limit) {
                flush_profile(current_tid, true);
                last_position = position - 1;
            } else {
                // Missing contexts inside a local interval are explicit zero
                // activity states in GATK.  Filling only bounded gaps avoids
                // allocating a chromosome-sized profile for sparse API inputs.
                for (std::int64_t missing = 1; missing <= gap; ++missing) {
                    const auto locus = last_position + static_cast<std::int32_t>(missing);
                    flush_profile(current_tid, false);
                    profile.add(locus, profile_probability(current_tid, locus));
                }
            }
            flush_profile(current_tid, false);
            profile.add(position, profile_probability(tid, position));
            last_position = position;
        }
        if (have_locus) flush_profile(current_tid, true);
        return result;
    }

    // AssemblyRegionIterator wraps the source pileup iterator in an
    // IntervalAlignmentContextIterator, so every coordinate inside -L is
    // added to ActivityProfile, including uncovered loci.  Materializing a
    // chromosome-worth of zero states here would defeat the compact Host /
    // Kokkos interface.  Instead, retain the exact finite Gaussian shoulders
    // on both sides of observed evidence and force-convert only across a gap
    // that is already too long to affect a source max-region decision.
    std::vector<ActivityProfileSpan> traversal_spans = options.traversal_spans;
    std::sort(traversal_spans.begin(), traversal_spans.end(),
              [](const auto& left, const auto& right) {
                  if (left.tid != right.tid) return left.tid < right.tid;
                  if (left.start != right.start) return left.start < right.start;
                  return left.end < right.end;
              });
    std::vector<ActivityProfileSpan> merged_spans;
    merged_spans.reserve(traversal_spans.size());
    for (const auto& span : traversal_spans) {
        if (!merged_spans.empty() && merged_spans.back().tid == span.tid &&
            span.start <= merged_spans.back().end) {
            merged_spans.back().end = std::max(merged_spans.back().end, span.end);
        } else {
            merged_spans.push_back(span);
        }
    }
    const auto add_profile_state = [&](const std::int32_t tid,
                                       const std::int32_t position) {
        flush_profile(tid, false);
        profile.add(position, profile_probability(tid, position));
    };
    const auto finite_shoulder = static_cast<std::int64_t>(profile.filter_size());
    std::size_t next_order = 0;
    for (const auto& span : merged_spans) {
        while (next_order < order.size()) {
            const auto index = order[next_order];
            const auto tid = input.tids[index];
            const auto position = input.positions[index];
            if (tid < span.tid || (tid == span.tid && position < span.start))
                throw std::invalid_argument(
                    "activity profile traversal spans do not cover an input locus");
            break;
        }
        const auto span_begin = next_order;
        while (next_order < order.size()) {
            const auto index = order[next_order];
            if (input.tids[index] != span.tid || input.positions[index] >= span.end)
                break;
            ++next_order;
        }
        if (span_begin == next_order) {
            // There is no compact Kokkos evidence row in this selected
            // interval, but GATK still emits its inactive ActivityProfile
            // chunks.  In force-active mode those become real assembly work.
            append_zero_profile_regions(span.tid, span.start, span.end);
            continue;
        }

        // Mutect2's quality-aware kernel emits a binary activity signal.  If
        // every compact pileup in a selected interval is inactive, the source
        // ActivityProfile is a constant zero sequence even when reads cover
        // part of that interval.  Do not split that sequence at the compact
        // evidence shoulder: AssemblyRegionIterator only cuts an all-zero
        // interval at the hard max-region boundary.  Besides avoiding work,
        // this preserves the IGV AssemblyRegion contract for a zero initial
        // LOD with reference-only pileups.
        bool span_has_nonzero_profile = false;
        for (std::size_t ordered = span_begin; ordered < next_order; ++ordered) {
            if (result.activity[order[ordered]] != 0.0) {
                span_has_nonzero_profile = true;
                break;
            }
        }
        if (!span_has_nonzero_profile) {
            append_zero_profile_regions(span.tid, span.start, span.end);
            continue;
        }

        const auto first_position = input.positions[order[span_begin]];
        const auto prefix_start = std::max<std::int64_t>(
            span.start, static_cast<std::int64_t>(first_position) - finite_shoulder);
        append_zero_profile_regions(span.tid, span.start, prefix_start);
        std::int32_t last_position = static_cast<std::int32_t>(prefix_start) - 1;
        for (std::int64_t position = prefix_start; position < first_position; ++position) {
            add_profile_state(span.tid, static_cast<std::int32_t>(position));
            last_position = static_cast<std::int32_t>(position);
        }
        for (std::size_t ordered = span_begin; ordered < next_order; ++ordered) {
            const auto index = order[ordered];
            const auto position = input.positions[index];
            const auto gap = static_cast<std::int64_t>(position) - last_position - 1;
            if (gap > gap_limit) {
                // Preserve the outgoing source Gaussian shoulder before an
                // inactive desert, then start the incoming shoulder exactly
                // where it can first influence a later active segment.
                const auto tail_end = std::min<std::int64_t>(
                    static_cast<std::int64_t>(position) - 1,
                    static_cast<std::int64_t>(last_position) + finite_shoulder);
                for (std::int64_t zero = static_cast<std::int64_t>(last_position) + 1;
                     zero <= tail_end; ++zero)
                    add_profile_state(span.tid, static_cast<std::int32_t>(zero));
                flush_profile(span.tid, true);
                const auto incoming_start = std::max<std::int64_t>(
                    span.start, static_cast<std::int64_t>(position) - finite_shoulder);
                append_zero_profile_regions(
                    span.tid, static_cast<std::int64_t>(tail_end) + 1, incoming_start);
                last_position = static_cast<std::int32_t>(incoming_start) - 1;
                for (std::int64_t zero = incoming_start; zero < position; ++zero) {
                    add_profile_state(span.tid, static_cast<std::int32_t>(zero));
                    last_position = static_cast<std::int32_t>(zero);
                }
            } else {
                for (std::int64_t zero = static_cast<std::int64_t>(last_position) + 1;
                     zero < position; ++zero) {
                    add_profile_state(span.tid, static_cast<std::int32_t>(zero));
                    last_position = static_cast<std::int32_t>(zero);
                }
            }
            add_profile_state(span.tid, position);
            last_position = position;
        }
        const auto shoulder_end = std::min<std::int64_t>(
            static_cast<std::int64_t>(span.end) - 1,
            static_cast<std::int64_t>(last_position) + finite_shoulder);
        for (std::int64_t zero = static_cast<std::int64_t>(last_position) + 1;
             zero <= shoulder_end; ++zero)
            add_profile_state(span.tid, static_cast<std::int32_t>(zero));
        flush_profile(span.tid, true);
        append_zero_profile_regions(
            span.tid, static_cast<std::int64_t>(shoulder_end) + 1, span.end);
    }
    if (next_order != order.size())
        throw std::invalid_argument("activity profile traversal spans do not cover an input locus");
    return result;
}

SomaticReferenceConfidenceResult calculate_somatic_reference_confidence_kokkos(
    const ActivityProfileInput& input, const std::uint8_t minimum_base_quality) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    const auto loci = input.tids.size();
    if (input.positions.size() != loci || input.reference_bases.size() != loci)
        throw std::invalid_argument("somatic reference-confidence loci are malformed");
    if (input.quality_offsets.size() != loci + 1U ||
        input.quality_bases.size() != input.quality_values.size() ||
        input.quality_bases.size() != input.quality_alt_flags.size() ||
        (!input.quality_offsets.empty() &&
         input.quality_offsets.back() != input.quality_bases.size()))
        throw std::invalid_argument("somatic reference-confidence observations are malformed");

    SomaticReferenceConfidenceResult result;
    result.log10_lods.assign(loci, 0.0);
    result.depths.assign(loci, 0U);
    result.reference_depths.assign(loci, 0U);
    result.non_reference_depths.assign(loci, 0U);
    result.execution_space = ExecSpace::name();
    if (loci == 0U) return result;

    std::size_t maximum_depth = 0U;
    for (std::size_t locus = 0; locus < loci; ++locus) {
        const auto begin = input.quality_offsets[locus];
        const auto end = input.quality_offsets[locus + 1U];
        if (end < begin) throw std::invalid_argument("somatic reference-confidence offsets are unsorted");
        maximum_depth = std::max(maximum_depth, static_cast<std::size_t>(end - begin));
    }
    std::size_t maximum_quality = 93U;
    for (const auto quality : input.quality_values)
        maximum_quality = std::max(maximum_quality, static_cast<std::size_t>(quality));
    const auto quality_table_size = maximum_quality + 1U;
    const auto binomial_stride = maximum_depth + 1U;

    const auto prepare_begin = std::chrono::steady_clock::now();
    U32View quality_offsets("somatic_rcm_quality_offsets", input.quality_offsets.size());
    U8View quality_bases("somatic_rcm_quality_bases", input.quality_bases.size());
    U8View quality_values("somatic_rcm_quality_values", input.quality_values.size());
    U8View quality_alt_flags("somatic_rcm_quality_alt_flags", input.quality_alt_flags.size());
    U8View reference_bases("somatic_rcm_reference_bases", input.reference_bases.size());
    DView digamma_values("somatic_rcm_digamma", maximum_depth + 2U);
    DView binomial_log_values("somatic_rcm_binomial", binomial_stride * binomial_stride);
    DView quality_epsilon_values("somatic_rcm_quality_epsilon", quality_table_size);
    DView quality_log_epsilon_values("somatic_rcm_quality_log_epsilon", quality_table_size);
    DView quality_log_one_minus_epsilon_values(
        "somatic_rcm_quality_log_one_minus_epsilon", quality_table_size);
    DView lods("somatic_rcm_lods", loci);
    U32View depths("somatic_rcm_depths", loci);
    U32View reference_depths("somatic_rcm_reference_depths", loci);
    U32View non_reference_depths("somatic_rcm_non_reference_depths", loci);

    auto host_offsets = Kokkos::create_mirror_view(quality_offsets);
    auto host_bases = Kokkos::create_mirror_view(quality_bases);
    auto host_values = Kokkos::create_mirror_view(quality_values);
    auto host_alt_flags = Kokkos::create_mirror_view(quality_alt_flags);
    auto host_reference = Kokkos::create_mirror_view(reference_bases);
    for (std::size_t index = 0; index < input.quality_offsets.size(); ++index)
        host_offsets(index) = input.quality_offsets[index];
    for (std::size_t index = 0; index < input.quality_bases.size(); ++index) {
        host_bases(index) = input.quality_bases[index];
        host_values(index) = input.quality_values[index];
        host_alt_flags(index) = input.quality_alt_flags[index];
    }
    for (std::size_t locus = 0; locus < loci; ++locus)
        host_reference(locus) = input.reference_bases[locus];
    Kokkos::deep_copy(quality_offsets, host_offsets);
    Kokkos::deep_copy(quality_bases, host_bases);
    Kokkos::deep_copy(quality_values, host_values);
    Kokkos::deep_copy(quality_alt_flags, host_alt_flags);
    Kokkos::deep_copy(reference_bases, host_reference);
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    const auto lookup_work = std::max({maximum_depth + 2U,
        binomial_log_values.extent(0), quality_table_size});
    Kokkos::parallel_for("somatic_reference_confidence_lookup",
        Kokkos::RangePolicy<ExecSpace>(0, lookup_work),
        KOKKOS_LAMBDA(const std::size_t index) {
            if (index < maximum_depth + 2U)
                digamma_values(index) = activity_digamma_positive(static_cast<double>(index));
            if (index < binomial_log_values.extent(0)) {
                const auto n = index / binomial_stride;
                const auto k = index % binomial_stride;
                binomial_log_values(index) = k <= n
                    ? activity_binomial_log(static_cast<std::uint64_t>(n),
                                            static_cast<std::uint64_t>(k))
                    : 0.0;
            }
            if (index < quality_table_size) {
                const auto quality = static_cast<double>(index);
                const auto log_epsilon = quality * (-Kokkos::log(10.0) / 10.0);
                quality_epsilon_values(index) = Kokkos::pow(10.0, -quality / 10.0);
                quality_log_epsilon_values(index) = log_epsilon;
                quality_log_one_minus_epsilon_values(index) =
                    log_epsilon < -0.69314718055994530942
                        ? Kokkos::log(1.0 - Kokkos::exp(log_epsilon))
                        : Kokkos::log(-Kokkos::expm1(log_epsilon));
            }
        });
    Kokkos::parallel_for("somatic_reference_confidence",
        Kokkos::RangePolicy<ExecSpace>(0, loci),
        KOKKOS_LAMBDA(const std::size_t locus) {
            const auto reference = reference_bases(locus);
            const auto begin = static_cast<std::size_t>(quality_offsets(locus));
            const auto end = static_cast<std::size_t>(quality_offsets(locus + 1U));
            std::uint32_t n_ref = 0U;
            std::uint32_t n_alt = 0U;
            for (std::size_t index = begin; index < end; ++index) {
                const auto quality = quality_values(index);
                if (quality <= minimum_base_quality) continue;
                const bool alternate = quality_alt_flags(index) != 0U ||
                    quality_bases(index) != reference;
                if (alternate) ++n_alt;
                else ++n_ref;
            }
            const auto n = static_cast<std::size_t>(n_ref) + static_cast<std::size_t>(n_alt);
            depths(locus) = static_cast<std::uint32_t>(n);
            reference_depths(locus) = n_ref;
            non_reference_depths(locus) = n_alt;
            if (n == 0U) {
                lods(locus) = 0.0;
                return;
            }
            const auto f_tilde_ratio = Kokkos::exp(
                digamma_values(static_cast<std::size_t>(n_ref) + 1U) -
                digamma_values(static_cast<std::size_t>(n_alt) + 1U));
            double read_sum = 0.0;
            for (std::size_t index = begin; index < end; ++index) {
                const auto quality = quality_values(index);
                if (quality <= minimum_base_quality) continue;
                const bool alternate = quality_alt_flags(index) != 0U ||
                    quality_bases(index) != reference;
                if (!alternate) continue;
                const auto quality_index = static_cast<std::size_t>(quality);
                const auto epsilon = quality_epsilon_values(quality_index);
                const auto z_bar_alt = (1.0 - epsilon) /
                    (1.0 - epsilon + epsilon * f_tilde_ratio);
                read_sum += z_bar_alt *
                    (quality_log_one_minus_epsilon_values(quality_index) -
                     quality_log_epsilon_values(quality_index)) +
                    activity_fast_bernoulli_entropy(z_bar_alt);
            }
            const auto beta_entropy = -Kokkos::log(static_cast<double>(n) + 1.0) -
                binomial_log_values(n * binomial_stride + static_cast<std::size_t>(n_alt));
            lods(locus) = (beta_entropy + read_sum) / Kokkos::log(10.0);
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_lods = Kokkos::create_mirror_view(lods);
    auto host_depths = Kokkos::create_mirror_view(depths);
    auto host_reference_depths = Kokkos::create_mirror_view(reference_depths);
    auto host_non_reference_depths = Kokkos::create_mirror_view(non_reference_depths);
    Kokkos::deep_copy(host_lods, lods);
    Kokkos::deep_copy(host_depths, depths);
    Kokkos::deep_copy(host_reference_depths, reference_depths);
    Kokkos::deep_copy(host_non_reference_depths, non_reference_depths);
    for (std::size_t locus = 0; locus < loci; ++locus) {
        result.log10_lods[locus] = host_lods(locus);
        result.depths[locus] = host_depths(locus);
        result.reference_depths[locus] = host_reference_depths(locus);
        result.non_reference_depths[locus] = host_non_reference_depths(locus);
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds = std::chrono::duration<double>(execute_end - execute_begin).count();
    return result;
}

}  // namespace fastgatk::kernels
