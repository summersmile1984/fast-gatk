// Bit-exact port of JAMA 1.0.3 LUDecomposition as bundled in GATK 4.6.2.0
// (Jama/LUDecomposition.class and Jama/Matrix.class, 2012-11-14).  GATK's
// MultivariateGaussian calls sigma.inverse() == new LUDecomposition(sigma)
// .solve(identity) and sigma.det() == new LUDecomposition(sigma).det(), so
// VQSR raw-bit parity requires the same left-looking Crout/Doolittle loop,
// the same Math.abs column pivoting, and the same solve/identity pass.  This
// header is deliberately standalone (host-only) so both the refresh_precision
// path and scoring can share one definition.
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace fastgatk::vqsr_jama {

// Returns true when the factorization succeeded (all diagonal factors
// non-zero).  On success `inverse` holds sigma^-1 computed with the exact
// JAMA row-pivoted forward/backward substitution against an identity RHS and
// `det` holds pivsign * product(LU[j][j]).
inline bool jama_lu_inverse_det(const std::vector<double>& matrix,
                                std::size_t n, std::vector<double>& inverse,
                                double& det) {
    // LU = row-major copy of the matrix.
    std::vector<double> lu(matrix);
    std::vector<int> piv(n);
    for (std::size_t i = 0; i < n; ++i) piv[i] = static_cast<int>(i);
    int pivsign = 1;

    std::vector<double> lucolj(n);
    for (std::size_t j = 0; j < n; ++j) {
        // Make a copy of the j-th column to localize references.
        for (std::size_t i = 0; i < n; ++i) lucolj[i] = lu[i * n + j];

        // Apply previous transformations.
        for (std::size_t i = 0; i < n; ++i) {
            // Most of the time is spent in the following dot product.
            const std::size_t kmax = i < j ? i : j;
            double s = 0.0;
            for (std::size_t k = 0; k < kmax; ++k)
                s += lu[i * n + k] * lucolj[k];
            lucolj[i] -= s;
            lu[i * n + j] = lucolj[i];
        }

        // Find pivot and exchange if necessary.
        std::size_t p = j;
        for (std::size_t i = j + 1; i < n; ++i) {
            if (std::abs(lucolj[i]) > std::abs(lucolj[p])) p = i;
        }
        if (p != j) {
            for (std::size_t k = 0; k < n; ++k) {
                const double t = lu[p * n + k];
                lu[p * n + k] = lu[j * n + k];
                lu[j * n + k] = t;
            }
            const int t = piv[p];
            piv[p] = piv[j];
            piv[j] = t;
            pivsign = -pivsign;
        }

        // Compute multipliers.
        if (j < n && lu[j * n + j] != 0.0) {
            for (std::size_t i = j + 1; i < n; ++i) lu[i * n + j] /= lu[j * n + j];
        }
    }

    // Singularity check mirrors JAMA isNonsingular().
    for (std::size_t j = 0; j < n; ++j)
        if (lu[j * n + j] == 0.0) return false;

    det = static_cast<double>(pivsign);
    for (std::size_t j = 0; j < n; ++j) det *= lu[j * n + j];

    // solve(identity): X = identity rows permuted by piv.
    inverse.assign(n * n, 0.0);
    for (std::size_t row = 0; row < n; ++row) {
        const std::size_t target = static_cast<std::size_t>(piv[row]);
        inverse[row * n + target] = 1.0;
    }
    // Solve L*Y = P*B (row permutation is already applied to X).
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = k + 1; i < n; ++i) {
            for (std::size_t column = 0; column < n; ++column)
                inverse[i * n + column] -= inverse[k * n + column] * lu[i * n + k];
        }
    }
    // Solve U*X = Y.
    for (std::size_t k = n; k-- > 0;) {
        for (std::size_t column = 0; column < n; ++column)
            inverse[k * n + column] /= lu[k * n + k];
        for (std::size_t i = 0; i < k; ++i) {
            for (std::size_t column = 0; column < n; ++column)
                inverse[i * n + column] -= inverse[k * n + column] * lu[i * n + k];
        }
    }
    return true;
}

}  // namespace fastgatk::vqsr_jama
