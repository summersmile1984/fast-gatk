#include "fastgatk/kernels/somatic.hpp"

#include <Kokkos_Core.hpp>

#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;
using ValueView = Kokkos::View<double*, MemorySpace>;
using OutputView = Kokkos::View<double*, MemorySpace>;
using CountView = Kokkos::View<std::size_t*, MemorySpace>;
using OffsetView = Kokkos::View<std::uint32_t*, MemorySpace>;

constexpr double kMissing = -1.0e299;
constexpr double kLog10ToNatural = 2.3025850929940456840179914546843642;
// These are the exact cutoffs used by SomaticLikelihoodsEngine.  A likelihood
// contribution is dropped only when its responsibility is strictly below
// 1e-10; x*log(x) is dropped only below 1e-8.  Keeping the inequalities named
// avoids accidentally changing the boundary behavior while refactoring the
// biallelic and multiallelic paths.
constexpr double kNegligibleResponsibility = 1.0e-10;
constexpr double kEntropyResponsibility = 1.0e-8;

KOKKOS_INLINE_FUNCTION bool valid_log_likelihood(const double value) {
    // GATK's NaturalLogUtils accepts -Infinity as a zero-probability allele
    // (and skips it in logSumExp), but a row containing one finite allele must
    // still contribute.  The native sentinel is also treated as missing.  A
    // NaN/+Infinity is not a valid likelihood and is ignored at this boundary
    // rather than allowing it to poison a whole candidate.
    return value > kMissing && Kokkos::isfinite(value);
}

KOKKOS_INLINE_FUNCTION double mix_log10(const double reference,
                                        const double alternate,
                                        const double fraction) {
    if (fraction <= 0.0) return reference;
    if (fraction >= 1.0) return alternate;
    const auto high = reference > alternate ? reference : alternate;
    if (high <= kMissing) return kMissing;
    const auto scaled = (1.0 - fraction) * Kokkos::pow(10.0, reference - high) +
                        fraction * Kokkos::pow(10.0, alternate - high);
    return high + Kokkos::log10(scaled > 0.0 ? scaled : Kokkos::pow(10.0, -300.0));
}

KOKKOS_INLINE_FUNCTION double orientation_log10(const std::uint32_t forward,
                                                const std::uint32_t reverse) {
    const auto total = static_cast<double>(forward + reverse);
    if (total <= 0.0) return 0.0;
    constexpr double error = 0.05;
    const auto balanced = total * Kokkos::log10(0.5);
    const auto forward_biased = static_cast<double>(forward) * Kokkos::log10(1.0 - error) +
                                static_cast<double>(reverse) * Kokkos::log10(error);
    const auto reverse_biased = static_cast<double>(forward) * Kokkos::log10(error) +
                                static_cast<double>(reverse) * Kokkos::log10(1.0 - error);
    const auto biased = forward_biased > reverse_biased ? forward_biased : reverse_biased;
    return biased > balanced ? biased - balanced : 0.0;
}

// Commons-Math 3.x Gamma.digamma for the positive alpha values used by the
// Mutect2 Dirichlet update.  Keep the same branch points and operation order
// as Gamma.digamma: the Java implementation recurs to x >= 49 and then uses
// the 1/12, 1/120 and 1/252 Bernoulli terms.  The 1/120 term is small at the
// branch point (~1.4e-9), but omitting it moves every variational responsibility
// and eventually the TLOD on sufficiently informative loci.
KOKKOS_INLINE_FUNCTION double digamma_positive(double value) {
    double x = value;
    double result = 0.0;
    if (!(x > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (x <= 1.0e-5)
        return -0.5772156649015329 - 1.0 / x;
    while (x < 49.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    // This is intentionally written in the same parenthesization as Commons
    // Math 3.6 Gamma.digamma (rather than as a generic asymptotic expansion):
    //   log(x) - 0.5/x - inv * (1/12 + inv * (1/120 - inv/252)),
    // where inv = 1/(x*x).
    const double inverse = 1.0 / (x * x);
    return result + Kokkos::log(x) - 0.5 / x - inverse * (
        1.0 / 12.0 + inverse * (1.0 / 120.0 - inverse / 252.0));
}

// Commons-Math's logGamma implementation.  SomaticLikelihoodsEngine evaluates
// Gamma.logGamma for every Dirichlet normalization; relying on the host libm
// lgamma gives backend- and libc-dependent last bits.  For x <= 8 Commons
// Math uses its NSWC logGamma1p approximation, rather than the Lanczos branch.
// Keeping that branch here matters for the small pseudocounts at the start of
// the variational update and avoids a systematic ~1e-15 drift in every
// evidence evaluation.  For x > 8 we use the pinned Lanczos coefficient table
// and operation order.  The alpha values produced by the variational update
// are positive, so reflection and pole handling are intentionally not needed.
KOKKOS_INLINE_FUNCTION double inv_gamma1pm1(const double x) {
    // These coefficients are Commons Math 3.x's port of the NSWC DGAM1
    // routine.  The Horner order below deliberately matches Gamma.java.
    constexpr double a0 = .611609510448141581788e-08;
    constexpr double a1 = .624730830116465516210e-08;
    constexpr double b1 = .203610414066806987300e+00;
    constexpr double b2 = .266205348428949217746e-01;
    constexpr double b3 = .493944979382446875238e-03;
    constexpr double b4 = -.851419432440314906588e-05;
    constexpr double b5 = -.643045481779353022248e-05;
    constexpr double b6 = .992641840672773722196e-06;
    constexpr double b7 = -.607761895722825260739e-07;
    constexpr double b8 = .195755836614639731882e-09;
    constexpr double p0 = .6116095104481415817861e-08;
    constexpr double p1 = .6871674113067198736152e-08;
    constexpr double p2 = .6820161668496170657918e-09;
    constexpr double p3 = .4686843322948848031080e-10;
    constexpr double p4 = .1572833027710446286995e-11;
    constexpr double p5 = -.1249441572276366213222e-12;
    constexpr double p6 = .4343529937408594255178e-14;
    constexpr double q1 = .3056961078365221025009e+00;
    constexpr double q2 = .5464213086042296536016e-01;
    constexpr double q3 = .4956830093825887312020e-02;
    constexpr double q4 = .2692369466186361192876e-03;
    constexpr double c = -.422784335098467139393487909917598e+00;
    constexpr double c0 = .577215664901532860606512090082402e+00;
    constexpr double c1 = -.655878071520253881077019515145390e+00;
    constexpr double c2 = -.420026350340952355290039348754298e-01;
    constexpr double c3 = .166538611382291489501700795102105e+00;
    constexpr double c4 = -.421977345555443367482083012891874e-01;
    constexpr double c5 = -.962197152787697356211492167234820e-02;
    constexpr double c6 = .721894324666309954239501034044657e-02;
    constexpr double c7 = -.116516759185906511211397108401839e-02;
    constexpr double c8 = -.215241674114950972815729963053648e-03;
    constexpr double c9 = .128050282388116186153198626328164e-03;
    constexpr double c10 = -.201348547807882386556893914210218e-04;
    constexpr double c11 = -.125049348214267065734535947383309e-05;
    constexpr double c12 = .113302723198169588237412962033074e-05;
    constexpr double c13 = -.205633841697760710345015413002057e-06;

    const double t = x <= 0.5 ? x : (x - 0.5) - 0.5;
    double ret;
    if (t < 0.0) {
        const double a = a0 + t * a1;
        double b = b8;
        b = b7 + t * b;
        b = b6 + t * b;
        b = b5 + t * b;
        b = b4 + t * b;
        b = b3 + t * b;
        b = b2 + t * b;
        b = b1 + t * b;
        b = 1.0 + t * b;

        double c_value = c13 + t * (a / b);
        c_value = c12 + t * c_value;
        c_value = c11 + t * c_value;
        c_value = c10 + t * c_value;
        c_value = c9 + t * c_value;
        c_value = c8 + t * c_value;
        c_value = c7 + t * c_value;
        c_value = c6 + t * c_value;
        c_value = c5 + t * c_value;
        c_value = c4 + t * c_value;
        c_value = c3 + t * c_value;
        c_value = c2 + t * c_value;
        c_value = c1 + t * c_value;
        c_value = c + t * c_value;
        if (x > 0.5)
            ret = t * c_value / x;
        else
            ret = x * ((c_value + 0.5) + 0.5);
    } else {
        double p = p6;
        p = p5 + t * p;
        p = p4 + t * p;
        p = p3 + t * p;
        p = p2 + t * p;
        p = p1 + t * p;
        p = p0 + t * p;

        double q = q4;
        q = q3 + t * q;
        q = q2 + t * q;
        q = q1 + t * q;
        q = 1.0 + t * q;

        double c_value = c13 + (p / q) * t;
        c_value = c12 + t * c_value;
        c_value = c11 + t * c_value;
        c_value = c10 + t * c_value;
        c_value = c9 + t * c_value;
        c_value = c8 + t * c_value;
        c_value = c7 + t * c_value;
        c_value = c6 + t * c_value;
        c_value = c5 + t * c_value;
        c_value = c4 + t * c_value;
        c_value = c3 + t * c_value;
        c_value = c2 + t * c_value;
        c_value = c1 + t * c_value;
        c_value = c0 + t * c_value;
        if (x > 0.5)
            ret = (t / x) * ((c_value - 0.5) - 0.5);
        else
            ret = x * c_value;
    }
    return ret;
}

KOKKOS_INLINE_FUNCTION double log_gamma1p(const double x) {
    return -Kokkos::log1p(inv_gamma1pm1(x));
}

KOKKOS_INLINE_FUNCTION double log_gamma_positive(double value) {
    if (!(value > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (value <= 8.0) {
        if (value < 0.5)
            return log_gamma1p(value) - Kokkos::log(value);
        if (value <= 2.5)
            return log_gamma1p((value - 0.5) - 0.5);
        const auto count = static_cast<int>(Kokkos::floor(value - 1.5));
        double product = 1.0;
        for (int index = 1; index <= count; ++index)
            product *= value - static_cast<double>(index);
        return log_gamma1p(value - static_cast<double>(count + 1)) +
               Kokkos::log(product);
    }
    constexpr double lanczos[] = {
        0.9999999999999971, 57.15623566586292, -59.59796035547549,
        14.136097974741746, -0.4919138160976202, 3.399464998481189e-5,
        4.652362892704858e-5, -9.837447530487956e-5,
        1.580887032249125e-4, -2.1026444172410488e-4,
        2.1743961811521265e-4, -1.643181065367639e-4,
        8.441822398385275e-5, -2.6190838401581408e-5,
        3.6899182659531625e-6};
    double x = value;
    double correction = 0.0;
    while (x < 8.0) {
        correction -= Kokkos::log(x);
        x += 1.0;
    }
    double lanczos_sum = 0.0;
    for (int index = 14; index > 0; --index)
        lanczos_sum += lanczos[index] / (x + static_cast<double>(index));
    lanczos_sum += lanczos[0];
    constexpr double g = 4.7421875;
    const double t = x + g + 0.5;
    constexpr double half_log_two_pi = 0.9189385332046727;
    return correction + (x + 0.5) * Kokkos::log(t) - t +
           half_log_two_pi + Kokkos::log(lanczos_sum / x);
}

KOKKOS_INLINE_FUNCTION double log_dirichlet_normalization_two(
    const double reference_alpha, const double alternate_alpha) {
    return log_gamma_positive(reference_alpha + alternate_alpha) -
           log_gamma_positive(reference_alpha) - log_gamma_positive(alternate_alpha);
}

KOKKOS_INLINE_FUNCTION double log_dirichlet_normalization(
    const ValueView& alpha, const std::size_t alpha_offset,
    const std::size_t allele_count, const std::size_t excluded_allele) {
    double sum = 0.0;
    double denominator = 0.0;
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        if (allele == excluded_allele) continue;
        const auto value = alpha(alpha_offset + allele);
        sum += value;
        denominator += log_gamma_positive(value);
    }
    return log_gamma_positive(sum) - denominator;
}

KOKKOS_INLINE_FUNCTION double log_dirichlet_prior(
    const std::size_t allele_count,
    const std::size_t excluded_allele,
    const double reference_pseudocount,
    const double alternate_pseudocount) {
    double sum = 0.0;
    double denominator = 0.0;
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        if (allele == excluded_allele) continue;
        const double value = allele == 0 ? reference_pseudocount : alternate_pseudocount;
        sum += value;
        denominator += log_gamma_positive(value);
    }
    return log_gamma_positive(sum) - denominator;
}

// General all-alleles-vs-without-alt evidence used by GATK's
// SomaticGenotypingEngine.  Alpha and effective-count vectors are held in
// caller-provided workspace so the device path never allocates dynamically.
// `excluded_allele == allele_count` denotes the all-alleles model.
KOKKOS_INLINE_FUNCTION void variational_multiallelic_evidence(
    const ValueView& likelihoods,
    const std::size_t likelihood_offset,
    const std::size_t allele_count,
    const std::size_t read_count,
    const std::size_t excluded_allele,
    const std::size_t tracked_allele,
    const ValueView& alpha_workspace,
    const std::size_t alpha_offset,
    const std::size_t counts_offset,
    const double reference_pseudocount,
    const double alternate_pseudocount,
    double& reference_log10,
    double& evidence_log10,
    double& tracked_fraction,
    std::size_t& covered) {
    reference_log10 = 0.0;
    evidence_log10 = 0.0;
    tracked_fraction = 0.0;
    covered = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        const auto value = likelihoods(likelihood_offset + read);
        if (valid_log_likelihood(value)) reference_log10 += value;
    }

    std::size_t active_alleles = 0;
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        if (allele == excluded_allele) {
            alpha_workspace(alpha_offset + allele) = 0.0;
            alpha_workspace(counts_offset + allele) = 0.0;
        } else {
            alpha_workspace(alpha_offset + allele) = 1.0;
            alpha_workspace(counts_offset + allele) = 0.0;
            ++active_alleles;
        }
    }
    if (active_alleles == 0) {
        evidence_log10 = reference_log10;
        return;
    }

    for (std::size_t read = 0; read < read_count; ++read) {
        bool any = false;
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (allele == excluded_allele) continue;
            const auto value = likelihoods(
                likelihood_offset + allele * read_count + read);
            if (valid_log_likelihood(value)) {
                any = true;
                break;
            }
        }
        if (any) ++covered;
    }
    if (covered == 0) return;

    for (int iteration = 0; iteration < 1000; ++iteration) {
        double sum_alpha = 0.0;
        for (std::size_t allele = 0; allele < allele_count; ++allele)
            if (allele != excluded_allele) sum_alpha += alpha_workspace(alpha_offset + allele);
        for (std::size_t allele = 0; allele < allele_count; ++allele)
            if (allele != excluded_allele) alpha_workspace(counts_offset + allele) = 0.0;

        for (std::size_t read = 0; read < read_count; ++read) {
            double high = -std::numeric_limits<double>::infinity();
            bool any = false;
            std::size_t max_allele = allele_count;
            for (std::size_t allele = 0; allele < allele_count; ++allele) {
                if (allele == excluded_allele) continue;
                const auto value = likelihoods(
                    likelihood_offset + allele * read_count + read);
                if (!valid_log_likelihood(value)) continue;
                const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                    digamma_positive(sum_alpha) + value * kLog10ToNatural;
                if (!any || log_value > high) {
                    high = log_value;
                    max_allele = allele;
                    any = true;
                }
            }
            if (!any) continue;
            // Commons Math/GATK NaturalLogUtils.logSumExp starts at the
            // selected maximum's unit contribution and then adds the other
            // terms.  Keep that exact operation order (including the
            // first-max tie break) instead of summing exp(0) in an arbitrary
            // allele position; this is observable in Strict raw-bit oracles.
            double denominator = 1.0;
            for (std::size_t allele = 0; allele < allele_count; ++allele) {
                if (allele == excluded_allele) continue;
                const auto value = likelihoods(
                    likelihood_offset + allele * read_count + read);
                if (!valid_log_likelihood(value)) continue;
                const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                    digamma_positive(sum_alpha) + value * kLog10ToNatural;
                if (allele != max_allele)
                    denominator += Kokkos::exp(log_value - high);
            }
            if (!(denominator > 0.0)) continue;
            for (std::size_t allele = 0; allele < allele_count; ++allele) {
                if (allele == excluded_allele) continue;
                const auto value = likelihoods(
                    likelihood_offset + allele * read_count + read);
                if (!valid_log_likelihood(value)) continue;
                const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                    digamma_positive(sum_alpha) + value * kLog10ToNatural;
                const auto log_sum = high + (denominator != 1.0 ? Kokkos::log(denominator) : 0.0);
                alpha_workspace(counts_offset + allele) +=
                    Kokkos::exp(log_value - log_sum);
            }
        }

        double distance = 0.0;
        double next_sum = 0.0;
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (allele == excluded_allele) continue;
            const auto prior = allele == 0 ? reference_pseudocount : alternate_pseudocount;
            const auto next = prior + alpha_workspace(counts_offset + allele);
            distance += Kokkos::abs(next - alpha_workspace(alpha_offset + allele));
            next_sum += next;
            alpha_workspace(alpha_offset + allele) = next;
        }
        if (distance / next_sum < 0.001) break;
    }

    double posterior_sum = 0.0;
    for (std::size_t allele = 0; allele < allele_count; ++allele)
        if (allele != excluded_allele) posterior_sum += alpha_workspace(alpha_offset + allele);
    if (tracked_allele < allele_count && tracked_allele != excluded_allele)
        tracked_fraction = alpha_workspace(alpha_offset + tracked_allele) / posterior_sum;
    if (covered == 0) {
        evidence_log10 = reference_log10;
        return;
    }

    double likelihood_entropy = 0.0;
    for (std::size_t read = 0; read < read_count; ++read) {
        double high = -std::numeric_limits<double>::infinity();
        bool any = false;
        std::size_t max_allele = allele_count;
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (allele == excluded_allele) continue;
            const auto value = likelihoods(
                likelihood_offset + allele * read_count + read);
            if (!valid_log_likelihood(value)) continue;
            const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                digamma_positive(posterior_sum) + value * kLog10ToNatural;
            if (!any || log_value > high) {
                high = log_value;
                max_allele = allele;
                any = true;
            }
        }
        if (!any) continue;
        double denominator = 1.0;
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (allele == excluded_allele) continue;
            const auto value = likelihoods(
                likelihood_offset + allele * read_count + read);
            if (!valid_log_likelihood(value)) continue;
            const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                digamma_positive(posterior_sum) + value * kLog10ToNatural;
            if (allele != max_allele)
                denominator += Kokkos::exp(log_value - high);
        }
        if (!(denominator > 0.0)) continue;
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (allele == excluded_allele) continue;
            const auto value = likelihoods(
                likelihood_offset + allele * read_count + read);
            if (!valid_log_likelihood(value)) continue;
            const auto log_value = digamma_positive(alpha_workspace(alpha_offset + allele)) -
                digamma_positive(posterior_sum) + value * kLog10ToNatural;
            const auto log_sum = high + (denominator != 1.0 ? Kokkos::log(denominator) : 0.0);
            const auto responsibility = Kokkos::exp(log_value - log_sum);
            if (responsibility >= kNegligibleResponsibility)
                likelihood_entropy += responsibility * value * kLog10ToNatural;
            if (responsibility >= kEntropyResponsibility)
                likelihood_entropy -= responsibility * Kokkos::log(responsibility);
        }
    }
    const auto prior_contribution = log_dirichlet_prior(
        allele_count, excluded_allele, reference_pseudocount, alternate_pseudocount);
    const auto posterior_contribution = -log_dirichlet_normalization(
        alpha_workspace, alpha_offset, allele_count, excluded_allele);
    evidence_log10 = (prior_contribution + posterior_contribution + likelihood_entropy) /
                     kLog10ToNatural;
}

// GATK's SomaticLikelihoodsEngine uses a two-state Dirichlet variational
// posterior for each biallelic candidate.  The input PairHMM values are
// log10 likelihoods; this helper performs the variational update and the
// evidence/entropy calculation in natural-log space, then returns log10
// evidence through the output references.  It is deliberately fixed-order:
// candidate and read loops are serial inside one Kokkos work item so Strict
// mode is independent of the execution backend.
KOKKOS_INLINE_FUNCTION void variational_biallelic_evidence(
    const ValueView& reference,
    const ValueView& alternate,
    const std::size_t offset,
    const std::size_t read_count,
    const double reference_pseudocount,
    const double alternate_pseudocount,
    double& reference_log10,
    double& evidence_log10,
    double& alternate_fraction,
    std::size_t& covered) {
    double reference_sum = 0.0;
    covered = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        const auto index = offset + read;
        const auto ref = reference(index);
        const auto alt = alternate(index);
        const bool valid_ref = valid_log_likelihood(ref);
        const bool valid_alt = valid_log_likelihood(alt);
        if (!valid_ref && !valid_alt) continue;
        if (valid_ref) reference_sum += ref * kLog10ToNatural;
        ++covered;
    }
    reference_log10 = reference_sum / kLog10ToNatural;
    evidence_log10 = reference_log10;
    alternate_fraction = 0.0;
    if (covered == 0) return;

    // minAF=0 in GATK gives one pseudocount for each state.
    double reference_alpha = 1.0;
    double alternate_alpha = 1.0;
    for (int iteration = 0; iteration < 1000; ++iteration) {
        const double sum_alpha = reference_alpha + alternate_alpha;
        const double log_weight_ref = digamma_positive(reference_alpha) -
                                       digamma_positive(sum_alpha);
        const double log_weight_alt = digamma_positive(alternate_alpha) -
                                       digamma_positive(sum_alpha);
        double reference_count = 0.0;
        double alternate_count = 0.0;
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto index = offset + read;
            const auto ref = reference(index);
            const auto alt = alternate(index);
            const bool valid_ref = valid_log_likelihood(ref);
            const bool valid_alt = valid_log_likelihood(alt);
            if (!valid_ref && !valid_alt) continue;
            const double ref_log = valid_ref
                ? log_weight_ref + ref * kLog10ToNatural
                : -std::numeric_limits<double>::infinity();
            const double alt_log = valid_alt
                ? log_weight_alt + alt * kLog10ToNatural
                : -std::numeric_limits<double>::infinity();
            const double high = ref_log > alt_log ? ref_log : alt_log;
            // Match NaturalLogUtils.logSumExp: the first (max) term is the
            // unit contribution and only the non-max term is exponentiated.
            const double denominator = 1.0 + Kokkos::exp(
                (ref_log > alt_log ? alt_log : ref_log) - high);
            if (!(denominator > 0.0)) continue;
            const double log_sum = high +
                (denominator != 1.0 ? Kokkos::log(denominator) : 0.0);
            const double alternate_responsibility = Kokkos::exp(
                alt_log - log_sum);
            alternate_count += alternate_responsibility;
            reference_count += 1.0 - alternate_responsibility;
        }
        const double next_reference_alpha = reference_pseudocount + reference_count;
        const double next_alternate_alpha = alternate_pseudocount + alternate_count;
        const double relative_l1 =
            (Kokkos::abs(next_reference_alpha - reference_alpha) +
             Kokkos::abs(next_alternate_alpha - alternate_alpha)) /
            (next_reference_alpha + next_alternate_alpha);
        reference_alpha = next_reference_alpha;
        alternate_alpha = next_alternate_alpha;
        if (relative_l1 < 0.001) break;
    }

    alternate_fraction = alternate_alpha / (reference_alpha + alternate_alpha);
    const double log_posterior_normalization =
        log_dirichlet_normalization_two(reference_alpha, alternate_alpha);
    double likelihood_entropy = 0.0;
    const double sum_alpha = reference_alpha + alternate_alpha;
    const double log_weight_ref = digamma_positive(reference_alpha) -
                                  digamma_positive(sum_alpha);
    const double log_weight_alt = digamma_positive(alternate_alpha) -
                                  digamma_positive(sum_alpha);
    for (std::size_t read = 0; read < read_count; ++read) {
        const auto index = offset + read;
        const auto ref = reference(index);
        const auto alt = alternate(index);
        const bool valid_ref = valid_log_likelihood(ref);
        const bool valid_alt = valid_log_likelihood(alt);
        if (!valid_ref && !valid_alt) continue;
        const double ref_log = valid_ref
            ? log_weight_ref + ref * kLog10ToNatural
            : -std::numeric_limits<double>::infinity();
        const double alt_log = valid_alt
            ? log_weight_alt + alt * kLog10ToNatural
            : -std::numeric_limits<double>::infinity();
        const double high = ref_log > alt_log ? ref_log : alt_log;
        const double denominator = 1.0 + Kokkos::exp(
            (ref_log > alt_log ? alt_log : ref_log) - high);
        if (!(denominator > 0.0)) continue;
        const double log_sum = high +
            (denominator != 1.0 ? Kokkos::log(denominator) : 0.0);
        const double ref_responsibility = Kokkos::exp(ref_log - log_sum);
        const double alt_responsibility = Kokkos::exp(alt_log - log_sum);
        const double entropy =
            (ref_responsibility >= kEntropyResponsibility ?
                ref_responsibility * Kokkos::log(ref_responsibility) : 0.0) +
            (alt_responsibility >= kEntropyResponsibility ?
                alt_responsibility * Kokkos::log(alt_responsibility) : 0.0);
        // This mirrors SomaticLikelihoodsEngine.likelihoodsContribution: in
        // particular, never evaluate 0 * -Infinity for a missing allele.
        const double reference_contribution =
            (!valid_ref || ref_responsibility < kNegligibleResponsibility) ? 0.0 :
                ref_responsibility * ref * kLog10ToNatural;
        const double alternate_contribution =
            (!valid_alt || alt_responsibility < kNegligibleResponsibility) ? 0.0 :
                alt_responsibility * alt * kLog10ToNatural;
        likelihood_entropy += reference_contribution + alternate_contribution - entropy;
    }
    const double prior_normalization = log_dirichlet_normalization_two(
        reference_pseudocount, alternate_pseudocount);
    const double evidence_natural = prior_normalization - log_posterior_normalization +
                                    likelihood_entropy;
    evidence_log10 = evidence_natural / kLog10ToNatural;
}

}  // namespace

SomaticLikelihoodResult calculate_somatic_likelihood_kokkos(
    const std::vector<double>& reference_read_likelihoods,
    const std::vector<double>& alternate_read_likelihoods,
    std::size_t candidate_count,
    std::size_t read_count,
    std::size_t allele_fraction_grid,
    const double minimum_allele_fraction) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (candidate_count == 0 || read_count == 0 || allele_fraction_grid < 2 ||
        allele_fraction_grid > 10001 || !std::isfinite(minimum_allele_fraction) ||
        minimum_allele_fraction < 0.0 || minimum_allele_fraction >= 1.0)
        throw std::invalid_argument("invalid somatic likelihood dimensions/grid");
    const auto expected = candidate_count * read_count;
    if (reference_read_likelihoods.size() != expected ||
        alternate_read_likelihoods.size() != expected)
        throw std::invalid_argument("somatic likelihood matrices do not match dimensions");
    // GATK's minimum-allele-fraction prior is Beta(1+epsilon, 1), where
    // epsilon is chosen so that minAF^epsilon = 1/2.  The Dirichlet
    // pseudocounts therefore remain [1,1] at the default minAF=0 and use the
    // release-compatible ALT pseudocount for an explicit positive value.
    const double alternate_pseudocount = minimum_allele_fraction == 0.0
        ? 1.0 : 1.0 - Kokkos::log(2.0) / Kokkos::log(minimum_allele_fraction);
    const double reference_pseudocount = 1.0;

    const auto prepare_begin = std::chrono::steady_clock::now();
    ValueView reference("somatic_reference_likelihoods", expected);
    ValueView alternate("somatic_alternate_likelihoods", expected);
    auto host_reference = Kokkos::create_mirror_view(reference);
    auto host_alternate = Kokkos::create_mirror_view(alternate);
    for (std::size_t index = 0; index < expected; ++index) {
        host_reference(index) = reference_read_likelihoods[index];
        host_alternate(index) = alternate_read_likelihoods[index];
    }
    Kokkos::deep_copy(reference, host_reference);
    Kokkos::deep_copy(alternate, host_alternate);
    OutputView tlod("somatic_tlod", candidate_count);
    OutputView normal_log10_odds("somatic_normal_log10_odds", candidate_count);
    OutputView best_fraction("somatic_best_af", candidate_count);
    OutputView reference_sum("somatic_reference_sum", candidate_count);
    OutputView best_sum("somatic_best_sum", candidate_count);
    CountView informative("somatic_informative_reads", candidate_count);
    Kokkos::deep_copy(tlod, 0.0);
    Kokkos::deep_copy(normal_log10_odds, 0.0);
    Kokkos::deep_copy(best_fraction, 0.0);
    Kokkos::deep_copy(reference_sum, 0.0);
    Kokkos::deep_copy(best_sum, 0.0);
    Kokkos::deep_copy(informative, static_cast<std::size_t>(0));
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_somatic_likelihood",
        Kokkos::RangePolicy<ExecSpace>(0, candidate_count),
        KOKKOS_LAMBDA(const std::size_t candidate) {
            double reference_likelihood = 0.0;
            double heterozygous_likelihood = 0.0;
            double evidence_likelihood = 0.0;
            double alternate_fraction = 0.0;
            std::size_t covered = 0;
            variational_biallelic_evidence(
                reference, alternate, candidate * read_count, read_count,
                reference_pseudocount, alternate_pseudocount,
                reference_likelihood, evidence_likelihood, alternate_fraction, covered);
            // NLOD is not derived from the variational somatic evidence.
            // GATK evaluates the fixed diploid ref/alt heterozygous model for
            // every retained normal fragment, then subtracts it from hom-ref.
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto index = candidate * read_count + read;
                const auto ref = reference(index);
                const auto alt = alternate(index);
                if (!valid_log_likelihood(ref)) continue;
                if (valid_log_likelihood(alt)) {
                    heterozygous_likelihood += mix_log10(ref, alt, 0.5);
                } else if (!Kokkos::isfinite(alt) && alt < 0.0) {
                    // NaturalLogUtils.logSumExp(ref, -Infinity) + log(1/2).
                    heterozygous_likelihood += ref + Kokkos::log10(0.5);
                } else {
                    continue;
                }
            }
            reference_sum(candidate) = reference_likelihood;
            best_sum(candidate) = evidence_likelihood;
            best_fraction(candidate) = alternate_fraction;
            informative(candidate) = covered;
            tlod(candidate) = evidence_likelihood - reference_likelihood;
            normal_log10_odds(candidate) = reference_likelihood - heterozygous_likelihood;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_tlod = Kokkos::create_mirror_view(tlod);
    auto host_normal_log10_odds = Kokkos::create_mirror_view(normal_log10_odds);
    auto host_fraction = Kokkos::create_mirror_view(best_fraction);
    auto host_reference_sum = Kokkos::create_mirror_view(reference_sum);
    auto host_best_sum = Kokkos::create_mirror_view(best_sum);
    auto host_informative = Kokkos::create_mirror_view(informative);
    Kokkos::deep_copy(host_tlod, tlod);
    Kokkos::deep_copy(host_normal_log10_odds, normal_log10_odds);
    Kokkos::deep_copy(host_fraction, best_fraction);
    Kokkos::deep_copy(host_reference_sum, reference_sum);
    Kokkos::deep_copy(host_best_sum, best_sum);
    Kokkos::deep_copy(host_informative, informative);

    SomaticLikelihoodResult result;
    result.tlod.resize(candidate_count);
    result.normal_log10_odds.resize(candidate_count);
    result.best_allele_fraction.resize(candidate_count);
    result.reference_log10_likelihood.resize(candidate_count);
    result.best_log10_likelihood.resize(candidate_count);
    result.informative_reads.resize(candidate_count);
    result.evidence_groups = read_count;
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        result.tlod[candidate] = host_tlod(candidate);
        result.normal_log10_odds[candidate] = host_normal_log10_odds(candidate);
        result.best_allele_fraction[candidate] = host_fraction(candidate);
        result.reference_log10_likelihood[candidate] = host_reference_sum(candidate);
        result.best_log10_likelihood[candidate] = host_best_sum(candidate);
        result.informative_reads[candidate] = host_informative(candidate);
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

SomaticLikelihoodResult calculate_somatic_likelihood_sparse_kokkos(
    const std::vector<double>& reference_group_likelihoods,
    const std::vector<double>& alternate_group_likelihoods,
    const std::vector<std::uint32_t>& candidate_offsets,
    const std::size_t candidate_count,
    const std::size_t evidence_group_count,
    const std::size_t allele_fraction_grid,
    const double minimum_allele_fraction) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (candidate_count == 0 || evidence_group_count == 0 || allele_fraction_grid < 2 ||
        allele_fraction_grid > 10001 || !std::isfinite(minimum_allele_fraction) ||
        minimum_allele_fraction < 0.0 || minimum_allele_fraction >= 1.0)
        throw std::invalid_argument("invalid sparse somatic likelihood dimensions/grid");
    if (candidate_offsets.size() != candidate_count + 1 ||
        reference_group_likelihoods.size() != alternate_group_likelihoods.size() ||
        candidate_offsets.front() != 0 ||
        candidate_offsets.back() != reference_group_likelihoods.size())
        throw std::invalid_argument("sparse somatic likelihood rows do not match offsets");
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        if (candidate_offsets[candidate] > candidate_offsets[candidate + 1])
            throw std::invalid_argument("sparse somatic likelihood offsets are not monotonic");
    }
    const double alternate_pseudocount = minimum_allele_fraction == 0.0
        ? 1.0 : 1.0 - Kokkos::log(2.0) / Kokkos::log(minimum_allele_fraction);
    const double reference_pseudocount = 1.0;

    const auto prepare_begin = std::chrono::steady_clock::now();
    ValueView reference("sparse_somatic_reference_likelihoods",
                        reference_group_likelihoods.size());
    ValueView alternate("sparse_somatic_alternate_likelihoods",
                        alternate_group_likelihoods.size());
    OffsetView offsets("sparse_somatic_candidate_offsets", candidate_offsets.size());
    auto host_reference = Kokkos::create_mirror_view(reference);
    auto host_alternate = Kokkos::create_mirror_view(alternate);
    auto host_offsets = Kokkos::create_mirror_view(offsets);
    for (std::size_t index = 0; index < reference_group_likelihoods.size(); ++index) {
        host_reference(index) = reference_group_likelihoods[index];
        host_alternate(index) = alternate_group_likelihoods[index];
    }
    for (std::size_t index = 0; index < candidate_offsets.size(); ++index)
        host_offsets(index) = candidate_offsets[index];
    Kokkos::deep_copy(reference, host_reference);
    Kokkos::deep_copy(alternate, host_alternate);
    Kokkos::deep_copy(offsets, host_offsets);
    OutputView tlod("sparse_somatic_tlod", candidate_count);
    OutputView normal_log10_odds("sparse_somatic_normal_log10_odds", candidate_count);
    OutputView best_fraction("sparse_somatic_best_af", candidate_count);
    OutputView reference_sum("sparse_somatic_reference_sum", candidate_count);
    OutputView best_sum("sparse_somatic_best_sum", candidate_count);
    CountView informative("sparse_somatic_informative_reads", candidate_count);
    Kokkos::deep_copy(tlod, 0.0);
    Kokkos::deep_copy(normal_log10_odds, 0.0);
    Kokkos::deep_copy(best_fraction, 0.0);
    Kokkos::deep_copy(reference_sum, 0.0);
    Kokkos::deep_copy(best_sum, 0.0);
    Kokkos::deep_copy(informative, static_cast<std::size_t>(0));
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_sparse_somatic_likelihood",
        Kokkos::RangePolicy<ExecSpace>(0, candidate_count),
        KOKKOS_LAMBDA(const std::size_t candidate) {
            const auto offset = static_cast<std::size_t>(offsets(candidate));
            const auto read_count = static_cast<std::size_t>(
                offsets(candidate + 1) - offsets(candidate));
            double reference_likelihood = 0.0;
            double heterozygous_likelihood = 0.0;
            double evidence_likelihood = 0.0;
            double alternate_fraction = 0.0;
            std::size_t covered = 0;
            variational_biallelic_evidence(
                reference, alternate, offset, read_count,
                reference_pseudocount, alternate_pseudocount,
                reference_likelihood, evidence_likelihood, alternate_fraction, covered);
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto index = offset + read;
                const auto ref = reference(index);
                const auto alt = alternate(index);
                if (!valid_log_likelihood(ref)) continue;
                if (valid_log_likelihood(alt)) {
                    heterozygous_likelihood += mix_log10(ref, alt, 0.5);
                } else if (!Kokkos::isfinite(alt) && alt < 0.0) {
                    heterozygous_likelihood += ref + Kokkos::log10(0.5);
                }
            }
            reference_sum(candidate) = reference_likelihood;
            best_sum(candidate) = evidence_likelihood;
            best_fraction(candidate) = alternate_fraction;
            informative(candidate) = covered;
            tlod(candidate) = evidence_likelihood - reference_likelihood;
            normal_log10_odds(candidate) = reference_likelihood - heterozygous_likelihood;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_tlod = Kokkos::create_mirror_view(tlod);
    auto host_normal_log10_odds = Kokkos::create_mirror_view(normal_log10_odds);
    auto host_fraction = Kokkos::create_mirror_view(best_fraction);
    auto host_reference_sum = Kokkos::create_mirror_view(reference_sum);
    auto host_best_sum = Kokkos::create_mirror_view(best_sum);
    auto host_informative = Kokkos::create_mirror_view(informative);
    Kokkos::deep_copy(host_tlod, tlod);
    Kokkos::deep_copy(host_normal_log10_odds, normal_log10_odds);
    Kokkos::deep_copy(host_fraction, best_fraction);
    Kokkos::deep_copy(host_reference_sum, reference_sum);
    Kokkos::deep_copy(host_best_sum, best_sum);
    Kokkos::deep_copy(host_informative, informative);

    SomaticLikelihoodResult result;
    result.tlod.resize(candidate_count);
    result.normal_log10_odds.resize(candidate_count);
    result.best_allele_fraction.resize(candidate_count);
    result.reference_log10_likelihood.resize(candidate_count);
    result.best_log10_likelihood.resize(candidate_count);
    result.informative_reads.resize(candidate_count);
    result.evidence_groups = evidence_group_count;
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        result.tlod[candidate] = host_tlod(candidate);
        result.normal_log10_odds[candidate] = host_normal_log10_odds(candidate);
        result.best_allele_fraction[candidate] = host_fraction(candidate);
        result.reference_log10_likelihood[candidate] = host_reference_sum(candidate);
        result.best_log10_likelihood[candidate] = host_best_sum(candidate);
        result.informative_reads[candidate] = host_informative(candidate);
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

SomaticLikelihoodResult calculate_somatic_multiallelic_likelihood_kokkos(
    const std::vector<double>& log10_likelihoods,
    const std::size_t allele_count,
    const std::size_t read_count,
    const double minimum_allele_fraction) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (allele_count < 2 || read_count == 0 || !std::isfinite(minimum_allele_fraction) ||
        minimum_allele_fraction < 0.0 || minimum_allele_fraction >= 1.0)
        throw std::invalid_argument("invalid multiallelic somatic dimensions");
    if (allele_count > std::numeric_limits<std::size_t>::max() / read_count ||
        log10_likelihoods.size() != allele_count * read_count)
        throw std::invalid_argument("multiallelic somatic matrix does not match dimensions");
    const double alternate_pseudocount = minimum_allele_fraction == 0.0
        ? 1.0 : 1.0 - Kokkos::log(2.0) / Kokkos::log(minimum_allele_fraction);
    const double reference_pseudocount = 1.0;
    const auto alternate_count = allele_count - 1;
    if (allele_count > std::numeric_limits<std::size_t>::max() / 4 ||
        alternate_count > std::numeric_limits<std::size_t>::max() /
            (4 * allele_count))
        throw std::invalid_argument("multiallelic somatic workspace is too large");

    const auto prepare_begin = std::chrono::steady_clock::now();
    ValueView likelihoods("somatic_multiallelic_likelihoods", log10_likelihoods.size());
    auto host_likelihoods = Kokkos::create_mirror_view(likelihoods);
    for (std::size_t index = 0; index < log10_likelihoods.size(); ++index)
        host_likelihoods(index) = log10_likelihoods[index];
    Kokkos::deep_copy(likelihoods, host_likelihoods);

    // Each ALT gets one all-alleles and one without-that-ALT evaluation. Each
    // evaluation owns an alpha and an effective-count row.
    const auto workspace_size = alternate_count * 4 * allele_count;
    ValueView workspace("somatic_multiallelic_workspace", workspace_size);
    OutputView tlod("somatic_multiallelic_tlod", alternate_count);
    OutputView normal_log10_odds("somatic_multiallelic_normal_log10_odds", alternate_count);
    OutputView best_fraction("somatic_multiallelic_af", alternate_count);
    OutputView reference_sum("somatic_multiallelic_reference", alternate_count);
    OutputView best_sum("somatic_multiallelic_best", alternate_count);
    CountView informative("somatic_multiallelic_informative", alternate_count);
    Kokkos::deep_copy(tlod, 0.0);
    Kokkos::deep_copy(normal_log10_odds, 0.0);
    Kokkos::deep_copy(best_fraction, 0.0);
    Kokkos::deep_copy(reference_sum, 0.0);
    Kokkos::deep_copy(best_sum, 0.0);
    Kokkos::deep_copy(informative, static_cast<std::size_t>(0));
    const auto prepare_end = std::chrono::steady_clock::now();

    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_somatic_multiallelic_likelihood",
        Kokkos::RangePolicy<ExecSpace>(0, alternate_count),
        KOKKOS_LAMBDA(const std::size_t alternate_index) {
            const auto tracked = alternate_index + 1;
            const auto base = alternate_index * 4 * allele_count;
            double all_reference = 0.0;
            double heterozygous_likelihood = 0.0;
            double all_evidence = 0.0;
            double all_fraction = 0.0;
            std::size_t all_covered = 0;
            variational_multiallelic_evidence(
                likelihoods, 0, allele_count, read_count, allele_count, tracked,
                workspace, base, base + allele_count,
                reference_pseudocount, alternate_pseudocount,
                all_reference, all_evidence, all_fraction, all_covered);
            double without_reference = 0.0;
            double without_evidence = 0.0;
            double unused_fraction = 0.0;
            std::size_t without_covered = 0;
            variational_multiallelic_evidence(
                likelihoods, 0, allele_count, read_count, tracked, allele_count,
                workspace, base + 2 * allele_count, base + 3 * allele_count,
                reference_pseudocount, alternate_pseudocount,
                without_reference, without_evidence, unused_fraction, without_covered);
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto ref = likelihoods(read);
                const auto alt = likelihoods(tracked * read_count + read);
                if (!valid_log_likelihood(ref)) continue;
                if (valid_log_likelihood(alt)) {
                    heterozygous_likelihood += mix_log10(ref, alt, 0.5);
                } else if (!Kokkos::isfinite(alt) && alt < 0.0) {
                    heterozygous_likelihood += ref + Kokkos::log10(0.5);
                }
            }
            reference_sum(alternate_index) = all_reference;
            best_sum(alternate_index) = all_evidence;
            best_fraction(alternate_index) = all_fraction;
            informative(alternate_index) = all_covered;
            tlod(alternate_index) = all_evidence - without_evidence;
            normal_log10_odds(alternate_index) = all_reference - heterozygous_likelihood;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_tlod = Kokkos::create_mirror_view(tlod);
    auto host_normal_log10_odds = Kokkos::create_mirror_view(normal_log10_odds);
    auto host_fraction = Kokkos::create_mirror_view(best_fraction);
    auto host_reference_sum = Kokkos::create_mirror_view(reference_sum);
    auto host_best_sum = Kokkos::create_mirror_view(best_sum);
    auto host_informative = Kokkos::create_mirror_view(informative);
    Kokkos::deep_copy(host_tlod, tlod);
    Kokkos::deep_copy(host_normal_log10_odds, normal_log10_odds);
    Kokkos::deep_copy(host_fraction, best_fraction);
    Kokkos::deep_copy(host_reference_sum, reference_sum);
    Kokkos::deep_copy(host_best_sum, best_sum);
    Kokkos::deep_copy(host_informative, informative);

    SomaticLikelihoodResult result;
    result.tlod.resize(alternate_count);
    result.normal_log10_odds.resize(alternate_count);
    result.best_allele_fraction.resize(alternate_count);
    result.reference_log10_likelihood.resize(alternate_count);
    result.best_log10_likelihood.resize(alternate_count);
    result.informative_reads.resize(alternate_count);
    result.evidence_groups = read_count;
    for (std::size_t alternate = 0; alternate < alternate_count; ++alternate) {
        result.tlod[alternate] = host_tlod(alternate);
        result.normal_log10_odds[alternate] = host_normal_log10_odds(alternate);
        result.best_allele_fraction[alternate] = host_fraction(alternate);
        result.reference_log10_likelihood[alternate] = host_reference_sum(alternate);
        result.best_log10_likelihood[alternate] = host_best_sum(alternate);
        result.informative_reads[alternate] = host_informative(alternate);
    }
    result.prepare_seconds =
        std::chrono::duration<double>(prepare_end - prepare_begin).count();
    result.seconds =
        std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

SomaticPosteriorResult calculate_somatic_posterior_kokkos(
    const std::vector<double>& tumor_reference_read_likelihoods,
    const std::vector<double>& tumor_alternate_read_likelihoods,
    const std::vector<double>& normal_reference_read_likelihoods,
    const std::vector<double>& normal_alternate_read_likelihoods,
    const std::vector<std::uint32_t>& f1r2,
    const std::vector<std::uint32_t>& r1f2,
    const std::size_t candidate_count,
    const std::size_t read_count,
    const double contamination,
    const double somatic_prior,
    const double germline_prior,
    const double artifact_prior,
    const std::size_t allele_fraction_grid) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (candidate_count == 0 || read_count == 0 || allele_fraction_grid < 2 ||
        allele_fraction_grid > 10001)
        throw std::invalid_argument("invalid somatic posterior dimensions/grid");
    const auto expected = candidate_count * read_count;
    if (tumor_reference_read_likelihoods.size() != expected ||
        tumor_alternate_read_likelihoods.size() != expected)
        throw std::invalid_argument("tumor somatic likelihood matrices do not match dimensions");
    const bool has_normal = !normal_reference_read_likelihoods.empty() ||
                            !normal_alternate_read_likelihoods.empty();
    // Tumor and normal are independent AlleleLikelihoods matrices in GATK;
    // their evidence counts need not match (and usually do not).  Derive the
    // normal row width from its candidate-major payload instead of silently
    // treating a different-width normal as all missing.  This keeps one
    // candidate's normal evidence aligned while allowing separate coverage.
    std::size_t normal_read_count = 0;
    if (has_normal) {
        if (normal_reference_read_likelihoods.empty() ||
            normal_alternate_read_likelihoods.empty() ||
            normal_reference_read_likelihoods.size() != normal_alternate_read_likelihoods.size() ||
            normal_reference_read_likelihoods.size() % candidate_count != 0)
            throw std::invalid_argument("normal somatic likelihood matrices do not match dimensions");
        normal_read_count = normal_reference_read_likelihoods.size() / candidate_count;
        if (normal_read_count == 0)
            throw std::invalid_argument("normal somatic likelihood matrices have no evidence rows");
    }
    if (f1r2.size() != candidate_count || r1f2.size() != candidate_count)
        throw std::invalid_argument("orientation counts do not match candidate dimensions");
    if (!std::isfinite(contamination) || contamination < 0.0 || contamination >= 1.0 ||
        !std::isfinite(somatic_prior) || !std::isfinite(germline_prior) ||
        !std::isfinite(artifact_prior) || somatic_prior <= 0.0 ||
        germline_prior <= 0.0 || artifact_prior <= 0.0)
        throw std::invalid_argument("somatic posterior contamination/priors are invalid");

    const auto prepare_begin = std::chrono::steady_clock::now();
    ValueView tumor_reference("somatic_posterior_tumor_reference", expected);
    ValueView tumor_alternate("somatic_posterior_tumor_alternate", expected);
    ValueView normal_reference("somatic_posterior_normal_reference",
                               has_normal ? candidate_count * normal_read_count : 1);
    ValueView normal_alternate("somatic_posterior_normal_alternate",
                               has_normal ? candidate_count * normal_read_count : 1);
    Kokkos::View<std::uint32_t*, MemorySpace> forward_counts(
        "somatic_posterior_f1r2", candidate_count);
    Kokkos::View<std::uint32_t*, MemorySpace> reverse_counts(
        "somatic_posterior_r1f2", candidate_count);
    auto host_tumor_reference = Kokkos::create_mirror_view(tumor_reference);
    auto host_tumor_alternate = Kokkos::create_mirror_view(tumor_alternate);
    auto host_normal_reference = Kokkos::create_mirror_view(normal_reference);
    auto host_normal_alternate = Kokkos::create_mirror_view(normal_alternate);
    auto host_forward = Kokkos::create_mirror_view(forward_counts);
    auto host_reverse = Kokkos::create_mirror_view(reverse_counts);
    for (std::size_t index = 0; index < expected; ++index) {
        host_tumor_reference(index) = tumor_reference_read_likelihoods[index];
        host_tumor_alternate(index) = tumor_alternate_read_likelihoods[index];
    }
    if (has_normal) {
        for (std::size_t index = 0; index < candidate_count * normal_read_count; ++index) {
            host_normal_reference(index) = normal_reference_read_likelihoods[index];
            host_normal_alternate(index) = normal_alternate_read_likelihoods[index];
        }
    } else {
        host_normal_reference(0) = kMissing;
        host_normal_alternate(0) = kMissing;
    }
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        host_forward(candidate) = f1r2[candidate];
        host_reverse(candidate) = r1f2[candidate];
    }
    Kokkos::deep_copy(tumor_reference, host_tumor_reference);
    Kokkos::deep_copy(tumor_alternate, host_tumor_alternate);
    Kokkos::deep_copy(normal_reference, host_normal_reference);
    Kokkos::deep_copy(normal_alternate, host_normal_alternate);
    Kokkos::deep_copy(forward_counts, host_forward);
    Kokkos::deep_copy(reverse_counts, host_reverse);

    OutputView somatic_probability("somatic_posterior_somatic", candidate_count);
    OutputView germline_probability("somatic_posterior_germline", candidate_count);
    OutputView artifact_probability("somatic_posterior_artifact", candidate_count);
    OutputView adjusted_fraction("somatic_posterior_adjusted_af", candidate_count);
    OutputView orientation_probability("somatic_posterior_orientation", candidate_count);
    OutputView somatic_evidence("somatic_posterior_somatic_log10", candidate_count);
    OutputView germline_evidence("somatic_posterior_germline_log10", candidate_count);
    OutputView artifact_evidence("somatic_posterior_artifact_log10", candidate_count);
    CountView informative("somatic_posterior_informative", candidate_count);
    Kokkos::deep_copy(somatic_probability, 0.0);
    Kokkos::deep_copy(germline_probability, 0.0);
    Kokkos::deep_copy(artifact_probability, 0.0);
    Kokkos::deep_copy(adjusted_fraction, 0.0);
    Kokkos::deep_copy(orientation_probability, 0.0);
    Kokkos::deep_copy(somatic_evidence, kMissing);
    Kokkos::deep_copy(germline_evidence, kMissing);
    Kokkos::deep_copy(artifact_evidence, kMissing);
    Kokkos::deep_copy(informative, static_cast<std::size_t>(0));
    const auto log_somatic_prior = Kokkos::log10(somatic_prior);
    const auto log_germline_prior = Kokkos::log10(germline_prior);
    const auto log_artifact_prior = Kokkos::log10(artifact_prior);
    const auto log_one_minus_contamination = Kokkos::log10(1.0 - contamination);
    const auto normal_available = has_normal;
    const auto execute_begin = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "fastgatk_somatic_posterior",
        Kokkos::RangePolicy<ExecSpace>(0, candidate_count),
        KOKKOS_LAMBDA(const std::size_t candidate) {
            double tumor_reference_sum = 0.0;
            double best_tumor = kMissing;
            double best_fraction = 0.0;
            double tumor_germline = kMissing;
            double normal_reference_sum = 0.0;
            double normal_germline = 0.0;
            double normal_hom_alt = 0.0;
            std::size_t normal_covered = 0;
            std::size_t covered = 0;
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto index = candidate * read_count + read;
                const auto tumor_ref = tumor_reference(index);
                const auto tumor_alt = tumor_alternate(index);
                if (!valid_log_likelihood(tumor_ref) ||
                    !valid_log_likelihood(tumor_alt)) continue;
                tumor_reference_sum += tumor_ref;
                ++covered;
            }
            informative(candidate) = covered;
            if (covered == 0) return;
            const auto effective_heterozygous_fraction = 0.5 * (1.0 - contamination);
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto index = candidate * read_count + read;
                const auto tumor_ref = tumor_reference(index);
                const auto tumor_alt = tumor_alternate(index);
                if (!valid_log_likelihood(tumor_ref) ||
                    !valid_log_likelihood(tumor_alt)) continue;
                tumor_germline += mix_log10(tumor_ref, tumor_alt,
                                            effective_heterozygous_fraction);
            }
            if (normal_available) {
                for (std::size_t normal_read = 0;
                     normal_read < normal_read_count; ++normal_read) {
                    const auto normal_index = candidate * normal_read_count + normal_read;
                    const auto normal_ref = normal_reference(normal_index);
                    const auto normal_alt = normal_alternate(normal_index);
                    if (!valid_log_likelihood(normal_ref) ||
                        !valid_log_likelihood(normal_alt)) continue;
                    normal_reference_sum += normal_ref;
                    normal_germline += mix_log10(normal_ref, normal_alt, 0.5);
                    normal_hom_alt += mix_log10(normal_ref, normal_alt, 1.0);
                    ++normal_covered;
                }
            }
            for (std::size_t grid = 0; grid < allele_fraction_grid; ++grid) {
                const auto fraction = static_cast<double>(grid) /
                                       static_cast<double>(allele_fraction_grid - 1);
                const auto effective_fraction = fraction * (1.0 - contamination);
                double mixture = 0.0;
                for (std::size_t read = 0; read < read_count; ++read) {
                    const auto index = candidate * read_count + read;
                    const auto tumor_ref = tumor_reference(index);
                    const auto tumor_alt = tumor_alternate(index);
                    if (!valid_log_likelihood(tumor_ref) ||
                        !valid_log_likelihood(tumor_alt)) continue;
                    mixture += mix_log10(tumor_ref, tumor_alt, effective_fraction);
                }
                if (mixture > best_tumor) {
                    best_tumor = mixture;
                    best_fraction = fraction;
                }
            }
            const auto normal_alt_evidence = normal_available && normal_covered > 0
                ? (normal_hom_alt > normal_germline ? normal_hom_alt : normal_germline) -
                  normal_reference_sum : 0.0;
            const auto normal_penalty = normal_alt_evidence > 0.0 ? normal_alt_evidence : 0.0;
            const auto somatic_log = best_tumor - normal_penalty + log_somatic_prior;
            const auto germline_log = tumor_germline +
                (normal_available && normal_covered > 0
                    ? normal_germline + log_germline_prior : log_germline_prior);
            const auto orientation_odds = orientation_log10(
                forward_counts(candidate), reverse_counts(candidate));
            // A balanced alternate signal is unlikely under the artifact
            // component; the orientation term can overcome this penalty only
            // when F1R2/R1F2 is strongly one-sided.
            const auto artifact_log = best_tumor + orientation_odds - 3.0 + log_artifact_prior;
            const auto high = somatic_log > germline_log
                ? (somatic_log > artifact_log ? somatic_log : artifact_log)
                : (germline_log > artifact_log ? germline_log : artifact_log);
            const auto denominator = Kokkos::pow(10.0, somatic_log - high) +
                                    Kokkos::pow(10.0, germline_log - high) +
                                    Kokkos::pow(10.0, artifact_log - high);
            const auto log_denominator = high +
                Kokkos::log10(denominator > 0.0 ? denominator : Kokkos::pow(10.0, -300.0));
            somatic_probability(candidate) = Kokkos::pow(10.0, somatic_log - log_denominator);
            germline_probability(candidate) = Kokkos::pow(10.0, germline_log - log_denominator);
            artifact_probability(candidate) = Kokkos::pow(10.0, artifact_log - log_denominator);
            const auto orientation_denominator = 1.0 + Kokkos::pow(10.0, -orientation_odds);
            orientation_probability(candidate) = 1.0 / orientation_denominator;
            const auto correction = Kokkos::pow(10.0, log_one_minus_contamination);
            const auto corrected = best_fraction / (correction > 1.0e-12 ? correction : 1.0e-12);
            adjusted_fraction(candidate) = corrected > 1.0 ? 1.0 : corrected;
            somatic_evidence(candidate) = best_tumor - normal_penalty;
            germline_evidence(candidate) = tumor_germline +
                (normal_available && normal_covered > 0 ? normal_germline : 0.0);
            artifact_evidence(candidate) = orientation_odds;
        });
    ExecSpace{}.fence();
    const auto execute_end = std::chrono::steady_clock::now();

    auto host_somatic = Kokkos::create_mirror_view(somatic_probability);
    auto host_germline = Kokkos::create_mirror_view(germline_probability);
    auto host_artifact = Kokkos::create_mirror_view(artifact_probability);
    auto host_adjusted = Kokkos::create_mirror_view(adjusted_fraction);
    auto host_orientation = Kokkos::create_mirror_view(orientation_probability);
    auto host_somatic_evidence = Kokkos::create_mirror_view(somatic_evidence);
    auto host_germline_evidence = Kokkos::create_mirror_view(germline_evidence);
    auto host_artifact_evidence = Kokkos::create_mirror_view(artifact_evidence);
    auto host_informative = Kokkos::create_mirror_view(informative);
    Kokkos::deep_copy(host_somatic, somatic_probability);
    Kokkos::deep_copy(host_germline, germline_probability);
    Kokkos::deep_copy(host_artifact, artifact_probability);
    Kokkos::deep_copy(host_adjusted, adjusted_fraction);
    Kokkos::deep_copy(host_orientation, orientation_probability);
    Kokkos::deep_copy(host_somatic_evidence, somatic_evidence);
    Kokkos::deep_copy(host_germline_evidence, germline_evidence);
    Kokkos::deep_copy(host_artifact_evidence, artifact_evidence);
    Kokkos::deep_copy(host_informative, informative);
    SomaticPosteriorResult result;
    result.somatic_probability.resize(candidate_count);
    result.germline_probability.resize(candidate_count);
    result.artifact_probability.resize(candidate_count);
    result.contamination_adjusted_allele_fraction.resize(candidate_count);
    result.orientation_bias_probability.resize(candidate_count);
    result.somatic_log10_evidence.resize(candidate_count);
    result.germline_log10_evidence.resize(candidate_count);
    result.artifact_log10_evidence.resize(candidate_count);
    result.informative_reads.resize(candidate_count);
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        result.somatic_probability[candidate] = host_somatic(candidate);
        result.germline_probability[candidate] = host_germline(candidate);
        result.artifact_probability[candidate] = host_artifact(candidate);
        result.contamination_adjusted_allele_fraction[candidate] = host_adjusted(candidate);
        result.orientation_bias_probability[candidate] = host_orientation(candidate);
        result.somatic_log10_evidence[candidate] = host_somatic_evidence(candidate);
        result.germline_log10_evidence[candidate] = host_germline_evidence(candidate);
        result.artifact_log10_evidence[candidate] = host_artifact_evidence(candidate);
        result.informative_reads[candidate] = host_informative(candidate);
    }
    result.prepare_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - prepare_begin).count();
    result.seconds = std::chrono::duration<double>(execute_end - execute_begin).count();
    result.execution_space = ExecSpace::name();
    return result;
}

}  // namespace fastgatk::kernels
