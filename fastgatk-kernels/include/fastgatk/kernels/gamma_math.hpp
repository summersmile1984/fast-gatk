// Shared Commons-Math 3.x gamma/digamma kernels for FilterMutectCalls and the
// somatic likelihood engines.  These device functions mirror Gamma.java
// exactly (NSWC logGamma1p branch for x <= 8 and the pinned Lanczos table
// above 8; Gamma.digamma recursion to x >= 49) so every numeric caller shares
// one implementation and no host libm lgamma/digamma drift can enter.
#pragma once

#include <Kokkos_Core.hpp>

#include <limits>

namespace fastgatk::kernels::detail {

// Commons-Math 3.x Gamma.digamma for positive values.
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
    const double inverse = 1.0 / (x * x);
    return result + Kokkos::log(x) - 0.5 / x - inverse * (
        1.0 / 12.0 + inverse * (1.0 / 120.0 - inverse / 252.0));
}

// NSWC DGAM1 coefficients ported from Commons Math 3.x Gamma.java.
KOKKOS_INLINE_FUNCTION double inv_gamma1pm1(const double x) {
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
           log_gamma_positive(reference_alpha) -
           log_gamma_positive(alternate_alpha);
}

}  // namespace fastgatk::kernels::detail
