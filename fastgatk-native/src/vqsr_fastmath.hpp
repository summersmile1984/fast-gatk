// Bit-exact ports of the Apache Commons Math 3.5 routines that GATK 4.6.2.0's
// VQSR VBEM path uses.  GATK's MultivariateGaussian calls Gamma.digamma, whose
// 3.5 implementation uses FastMath.log (NOT Math.log); both are ported here
// with the exact Java operation order and the exact table bits dumped from the
// bundled jar, so results are raw-bit identical to the Java runtime.
//
// Only the hiPrec == null paths are ported (the digamma/Normal-Wishart use
// case); the hiPrec variants are intentionally omitted.
#pragma once

#include <bit>
#include <cstdint>
#include <limits>

#include "vqsr_fastmath_tables.hpp"

namespace fastgatk::vqsr_fastmath {

static constexpr double HEX_40000000 = 1073741824.0;
static constexpr double TWO_POWER_52 = 4503599627370496.0;

// Exact port of commons-math3 3.5 FastMath.log(double x) == log(x, null).
inline double fastmath_log(double x) {
    if (x == 0.0)  // Handle special case of +0/-0
        return -std::numeric_limits<double>::infinity();
    std::int64_t bits = static_cast<std::int64_t>(std::bit_cast<std::uint64_t>(x));

    // Handle special cases of negative input, and NaN
    if (((bits & static_cast<std::int64_t>(0x8000000000000000ULL)) != 0 || x != x) && x != 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    // Handle special cases of positive infinity.
    if (x == std::numeric_limits<double>::infinity())
        return std::numeric_limits<double>::infinity();

    // Extract the exponent
    int exp = static_cast<int>(bits >> 52) - 1023;

    if ((bits & 0x7ff0000000000000LL) == 0) {
        // Subnormal! (x == 0 was handled above)
        bits <<= 1;
        while ((bits & 0x0010000000000000LL) == 0) {
            --exp;
            bits <<= 1;
        }
    }

    if ((exp == -1 || exp == 0) && x < 1.01 && x > 0.99) {
        // Straight polynomial expansion in higher precision (hiPrec == null).
        double xa = x - 1.0;
        double xb = xa - x + 1.0;
        double tmp = xa * HEX_40000000;
        double aa = xa + tmp - tmp;
        double ab = xa - aa;
        xa = aa;
        xb = ab;

        const double* lnCoef_last = LN_QUICK_COEF[LN_QUICK_COEF_LEN - 1];
        double ya = lnCoef_last[0];
        double yb = lnCoef_last[1];

        for (int i = LN_QUICK_COEF_LEN - 2; i >= 0; --i) {
            // Multiply a = y * x
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            // split, so now y = a
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;

            // Add  a = y + lnQuickCoef
            const double* lnCoef_i = LN_QUICK_COEF[i];
            aa = ya + lnCoef_i[0];
            ab = yb + lnCoef_i[1];
            // Split y = a
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }

        // Multiply a = y * x
        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;
        // split, so now y = a
        tmp = aa * HEX_40000000;
        ya = aa + tmp - tmp;
        yb = aa - ya + ab;

        return ya + yb;
    }

    // lnm is a log of a number in the range of 1.0 - 2.0, so 0 <= lnm < ln(2)
    const double* lnm = LN_MANT[(bits & 0x000ffc0000000000LL) >> 42];

    // y is the most significant 10 bits of the mantissa; epsilon = (x - y) / y
    const double epsilon = static_cast<double>(bits & 0x3ffffffffffLL) /
        static_cast<double>(static_cast<std::uint64_t>(TWO_POWER_52) +
                            static_cast<std::uint64_t>(bits & 0x000ffc0000000000LL));

    double lnza = 0.0;
    double lnzb = 0.0;

    // hiPrec == null: Remez polynomial in standard double precision.
    lnza = -0.16624882440418567;
    lnza = lnza * epsilon + 0.19999954120254515;
    lnza = lnza * epsilon + -0.2499999997677497;
    lnza = lnza * epsilon + 0.3333333333332802;
    lnza = lnza * epsilon + -0.5;
    lnza = lnza * epsilon + 1.0;
    lnza *= epsilon;

    // Compute the following sum:
    //  lnzb + lnm[1] + ln2B*exp + lnza + lnm[0] + ln2A*exp;
    double a = LN_2_A * exp;
    double b = 0.0;
    double c = a + lnm[0];
    double d = -(c - a - lnm[0]);
    a = c;
    b += d;

    c = a + lnza;
    d = -(c - a - lnza);
    a = c;
    b += d;

    c = a + LN_2_B * exp;
    d = -(c - a - LN_2_B * exp);
    a = c;
    b += d;

    c = a + lnm[1];
    d = -(c - a - lnm[1]);
    a = c;
    b += d;

    c = a + lnzb;
    d = -(c - a - lnzb);
    a = c;
    b += d;

    return a + b;
}

// Exact port of commons-math3 3.5 Gamma.digamma (the version bundled with
// GATK 4.6.2.0).  The recurrence is recursive in Java (digamma(x+1) - 1/x),
// so the subtraction association differs from an iterative accumulation; the
// recursion is kept for raw-bit parity.  S_LIMIT = 1e-5, C_LIMIT = 49 and the
// asymptotic series divides by 252.0 (rather than multiplying by 1.0/252.0).
inline double gatk_digamma_exact(double x) {
    if (x > 0.0 && x <= 1.0e-5)
        return -0.5772156649015329 - 1.0 / x;
    if (x >= 49.0) {
        const double inv = 1.0 / (x * x);
        // FastMath.log(x) - 0.5/x - inv * (1/12 + inv * (1/120 - inv/252.0))
        return fastmath_log(x) - 0.5 / x -
            inv * (0.08333333333333333 + inv * (0.008333333333333333 - inv / 252.0));
    }
    if (x <= 0.0 || x != x)
        // commons-math3 3.5: x <= 0 falls into the recurrence, producing NaN
        // through digamma(x+1) chains only for negative integers; for other
        // non-positive x the recurrence still applies.  Keep the same control
        // flow as the bytecode: any x < 49 that is not in (0, 1e-5] recurses.
        return gatk_digamma_exact(x + 1.0) - 1.0 / x;
    return gatk_digamma_exact(x + 1.0) - 1.0 / x;
}

}  // namespace fastgatk::vqsr_fastmath
