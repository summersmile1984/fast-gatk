#include "fastgatk/contamination_kernel_segmenter.hpp"
#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>

namespace fastgatk::contamination {
namespace {

using Matrix = std::vector<std::vector<double>>;
constexpr double EPSILON = 1.0e-10;
constexpr double KERNEL_VARIANCE = 0.025;
constexpr std::int64_t RANDOM_SEED = 1216;

void record_segmenter_plan(SegmenterTelemetry* telemetry,
                           const fastgatk::core::PlanTelemetry& plan,
                           std::size_t observations) {
    if (telemetry == nullptr) return;
    ++telemetry->batches;
    telemetry->observations += observations;
    telemetry->prepare_seconds += plan.prepare_seconds;
    telemetry->execute_seconds += plan.execute_seconds;
    telemetry->execution_space = Kokkos::DefaultExecutionSpace::name();
}

int fast_exponent(double value) {
    return static_cast<int>((std::bit_cast<std::uint64_t>(value) >> 52) & 0x7ffULL) - 1023;
}

double fast_scalb(double value, int exponent) {
    if (exponent > -1023 && exponent < 1024)
        return value * std::bit_cast<double>(static_cast<std::uint64_t>(exponent + 1023) << 52);
    if (std::isnan(value) || std::isinf(value) || value == 0.0) return value;
    if (exponent < -2098) return std::signbit(value) ? -0.0 : 0.0;
    if (exponent > 2097) return std::signbit(value) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto sign = bits & 0x8000000000000000ULL;
    const int raw_exponent = static_cast<int>((bits >> 52) & 0x7ffULL);
    std::uint64_t mantissa = bits & 0x000fffffffffffffULL;
    int scaled_exponent = raw_exponent + exponent;
    if (exponent < 0) {
        if (scaled_exponent > 0) return std::bit_cast<double>(sign | (static_cast<std::uint64_t>(scaled_exponent) << 52) | mantissa);
        if (scaled_exponent > -53) {
            mantissa |= 1ULL << 52;
            const auto lost = mantissa & (1ULL << (-scaled_exponent));
            mantissa >>= 1 - scaled_exponent;
            if (lost != 0) ++mantissa;
            return std::bit_cast<double>(sign | mantissa);
        }
        return sign == 0 ? 0.0 : -0.0;
    }
    if (raw_exponent == 0) {
        while ((mantissa >> 52) != 1) { mantissa <<= 1; --scaled_exponent; }
        ++scaled_exponent;
        mantissa &= 0x000fffffffffffffULL;
    }
    if (scaled_exponent < 2047) return std::bit_cast<double>(sign | (static_cast<std::uint64_t>(scaled_exponent) << 52) | mantissa);
    return sign == 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
}

double fast_hypot(double x, double y) {
    if (std::isinf(x) || std::isinf(y)) return std::numeric_limits<double>::infinity();
    if (std::isnan(x) || std::isnan(y)) return std::numeric_limits<double>::quiet_NaN();
    const int exponent_x = fast_exponent(x), exponent_y = fast_exponent(y);
    if (exponent_x > exponent_y + 27) return std::abs(x);
    if (exponent_y > exponent_x + 27) return std::abs(y);
    const int middle = (exponent_x + exponent_y) / 2;
    const double scaled_x = fast_scalb(x, -middle), scaled_y = fast_scalb(y, -middle);
    return fast_scalb(std::sqrt(scaled_x * scaled_x + scaled_y * scaled_y), middle);
}

// java.util.Random.nextInt(bound), including its rejection step.  The
// ContaminationSegmenter uses this exact seeded generator for sampling with
// replacement, so using std::uniform_int_distribution would change the
// reduced kernel matrix.
class JavaRandom {
public:
    explicit JavaRandom(std::int64_t seed)
        : state((static_cast<std::uint64_t>(seed) ^ 0x5DEECE66DULL) & ((1ULL << 48) - 1)) {}

    std::int32_t next_int(std::int32_t bound) {
        if (bound <= 0) return 0;
        const auto m = bound - 1;
        if ((bound & m) == 0)
            return static_cast<std::int32_t>((static_cast<std::int64_t>(bound) * next(31)) >> 31);
        std::int32_t bits;
        std::int32_t value;
        do {
            bits = next(31);
            value = bits % bound;
        } while (bits - value + m < 0);
        return value;
    }

private:
    std::uint64_t state;

    std::int32_t next(int bits) {
        state = (state * 0x5DEECE66DULL + 0xBULL) & ((1ULL << 48) - 1);
        return static_cast<std::int32_t>(state >> (48 - bits));
    }
};

KOKKOS_INLINE_FUNCTION double kernel(double first, double second) {
    const double first_complement = 1.0 - first;
    const double second_complement = 1.0 - second;
    const double maf_first = first < first_complement ? first : first_complement;
    const double maf_second = second < second_complement ? second : second_complement;
    const double difference = maf_first - maf_second;
    return Kokkos::exp(-(difference * difference) / (2.0 * KERNEL_VARIANCE));
}

Matrix symmetric_kernel_matrix(const std::vector<double>& subsample,
                               SegmenterTelemetry* telemetry) {
    const std::size_t n = subsample.size();
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("contamination-segmenter-symmetric-v1");
    host_batch.records = n * n;
    host_batch.bytes = n * sizeof(double) + n * n * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("contamination-segmenter-symmetric");
    kernel_plan.begin_prepare(host_batch);
    Kokkos::View<double*> subsample_view("contamination_subsample", n);
    auto host_subsample = Kokkos::create_mirror_view(subsample_view);
    for (std::size_t i = 0; i < n; ++i) host_subsample(i) = subsample[i];
    Kokkos::deep_copy(subsample_view, host_subsample);
    Kokkos::View<double**> matrix_view("contamination_sub_kernel", n, n);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(n * n);
    device_batch.bind("subsample", subsample_view);
    device_batch.bind("matrix", matrix_view);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "contamination_sub_kernel_fill",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0, 0}, {n, n}),
        KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
            matrix_view(i, j) = kernel(subsample_view(i), subsample_view(j));
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    record_segmenter_plan(telemetry, kernel_plan.telemetry(), n * n);
    Matrix result(n, std::vector<double>(n));
    auto copied = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrix_view);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            result[i][j] = copied(i, j);
    return result;
}

struct EigenDecomposition {
    std::vector<double> values;
    Matrix vectors;
};

// Deterministic Jacobi decomposition of the symmetric Gaussian kernel.  The
// kernel is positive semidefinite, so its eigen-decomposition is also the SVD
// needed by KernelSegmenter.  Keeping the matrix host-resident here preserves
// Java's deterministic accumulation order while the matrix construction and
// projection remain Kokkos operations.
[[maybe_unused]] EigenDecomposition symmetric_eigen(Matrix matrix) {
    const std::size_t n = matrix.size();
    Matrix vectors(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) vectors[i][i] = 1.0;
    if (n == 0) return {{}, {}};
    const std::size_t max_sweeps = std::max<std::size_t>(32, 12 * n * n);
    for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
        std::size_t p = 0, q = 0;
        double largest = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const double candidate = std::abs(matrix[i][j]);
                if (candidate > largest) { largest = candidate; p = i; q = j; }
            }
        }
        if (!(largest > 1.0e-14)) break;
        const double angle = 0.5 * std::atan2(2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        const double app = matrix[p][p];
        const double aqq = matrix[q][q];
        const double apq = matrix[p][q];
        matrix[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = matrix[q][p] = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            if (k == p || k == q) continue;
            const double mkp = matrix[k][p];
            const double mkq = matrix[k][q];
            matrix[k][p] = matrix[p][k] = cosine * mkp - sine * mkq;
            matrix[k][q] = matrix[q][k] = sine * mkp + cosine * mkq;
        }
        for (std::size_t k = 0; k < n; ++k) {
            const double vkp = vectors[k][p];
            const double vkq = vectors[k][q];
            vectors[k][p] = cosine * vkp - sine * vkq;
            vectors[k][q] = sine * vkp + cosine * vkq;
        }
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t first, std::size_t second) {
        return matrix[first][first] > matrix[second][second];
    });
    Matrix sorted_vectors(n, std::vector<double>(n));
    std::vector<double> values(n);
    for (std::size_t column = 0; column < n; ++column) {
        values[column] = std::max(0.0, matrix[order[column]][order[column]]);
        for (std::size_t row = 0; row < n; ++row)
            sorted_vectors[row][column] = vectors[row][order[column]];
    }
    return {std::move(values), std::move(sorted_vectors)};
}

// Commons-Math/JAMA Golub-Reinsch SVD for the square kernel matrix.  The
// kernel is symmetric positive semidefinite, but using the bidiagonal QR path
// rather than a generic eigensolver preserves the singular-value ordering and
// the U-column scaling expected by KernelSegmenter.
EigenDecomposition commons_math_svd(const Matrix& input) {
    const int m = static_cast<int>(input.size());
    const int n = m == 0 ? 0 : static_cast<int>(input.front().size());
    Matrix a = input;
    std::vector<double> singular(static_cast<std::size_t>(n), 0.0);
    std::vector<double> e(static_cast<std::size_t>(n), 0.0);
    std::vector<double> work(static_cast<std::size_t>(m), 0.0);
    Matrix u(static_cast<std::size_t>(m), std::vector<double>(static_cast<std::size_t>(n), 0.0));
    Matrix v(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n), 0.0));
    const int nct = std::min(m - 1, n);
    const int nrt = std::max(0, std::min(n - 2, m));
    const int limit = std::max(nct, nrt);
    for (int k = 0; k < limit; ++k) {
        if (k < nct) {
            singular[static_cast<std::size_t>(k)] = 0.0;
            for (int i = k; i < m; ++i)
                    singular[static_cast<std::size_t>(k)] = fast_hypot(singular[static_cast<std::size_t>(k)], a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]);
            if (singular[static_cast<std::size_t>(k)] != 0.0) {
                if (a[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)] < 0.0)
                    singular[static_cast<std::size_t>(k)] = -singular[static_cast<std::size_t>(k)];
                for (int i = k; i < m; ++i) a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] /= singular[static_cast<std::size_t>(k)];
                a[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)] += 1.0;
            }
            singular[static_cast<std::size_t>(k)] = -singular[static_cast<std::size_t>(k)];
        }
        for (int j = k + 1; j < n; ++j) {
            if (k < nct && singular[static_cast<std::size_t>(k)] != 0.0) {
                double t = 0.0;
                for (int i = k; i < m; ++i) t += a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                t = -t / a[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
                for (int i = k; i < m; ++i) a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += t * a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            }
            e[static_cast<std::size_t>(j)] = a[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
        }
        if (k < nct)
            for (int i = k; i < m; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
        if (k < nrt) {
            e[static_cast<std::size_t>(k)] = 0.0;
            for (int i = k + 1; i < n; ++i) e[static_cast<std::size_t>(k)] = fast_hypot(e[static_cast<std::size_t>(k)], e[static_cast<std::size_t>(i)]);
            if (e[static_cast<std::size_t>(k)] != 0.0) {
                if (e[static_cast<std::size_t>(k + 1)] < 0.0) e[static_cast<std::size_t>(k)] = -e[static_cast<std::size_t>(k)];
                for (int i = k + 1; i < n; ++i) e[static_cast<std::size_t>(i)] /= e[static_cast<std::size_t>(k)];
                e[static_cast<std::size_t>(k + 1)] += 1.0;
            }
            e[static_cast<std::size_t>(k)] = -e[static_cast<std::size_t>(k)];
            if (k + 1 < m && e[static_cast<std::size_t>(k)] != 0.0) {
                std::fill(work.begin() + k + 1, work.end(), 0.0);
                for (int j = k + 1; j < n; ++j)
                    for (int i = k + 1; i < m; ++i) work[static_cast<std::size_t>(i)] += e[static_cast<std::size_t>(j)] * a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                for (int j = k + 1; j < n; ++j) {
                    const double t = -e[static_cast<std::size_t>(j)] / e[static_cast<std::size_t>(k + 1)];
                    for (int i = k + 1; i < m; ++i) a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += t * work[static_cast<std::size_t>(i)];
                }
            }
            for (int i = k + 1; i < n; ++i) v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = e[static_cast<std::size_t>(i)];
        }
    }
    if (nct < n) singular[static_cast<std::size_t>(nct)] = a[static_cast<std::size_t>(nct)][static_cast<std::size_t>(nct)];
    if (m < n) singular[static_cast<std::size_t>(n - 1)] = 0.0;
    if (nrt + 1 < n) e[static_cast<std::size_t>(nrt)] = a[static_cast<std::size_t>(nrt)][static_cast<std::size_t>(n - 1)];
    if (n > 0) e[static_cast<std::size_t>(n - 1)] = 0.0;
    for (int j = nct; j < n; ++j) {
        for (int i = 0; i < m; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0.0;
        if (j < m) u[static_cast<std::size_t>(j)][static_cast<std::size_t>(j)] = 1.0;
    }
    for (int k = nct - 1; k >= 0; --k) {
        if (singular[static_cast<std::size_t>(k)] != 0.0) {
            for (int j = k + 1; j < n; ++j) {
                double t = 0.0;
                for (int i = k; i < m; ++i) t += u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                t = -t / u[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
                for (int i = k; i < m; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += t * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            }
            for (int i = k; i < m; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = -u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            u[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)] = 1.0 + u[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
            for (int i = 0; i < k; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = 0.0;
        } else {
            for (int i = 0; i < m; ++i) u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = 0.0;
            u[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)] = 1.0;
        }
    }
    for (int k = n - 1; k >= 0; --k) {
        if (k < nrt && e[static_cast<std::size_t>(k)] != 0.0) {
            for (int j = k + 1; j < n; ++j) {
                double t = 0.0;
                for (int i = k + 1; i < n; ++i) t += v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                t = -t / v[static_cast<std::size_t>(k + 1)][static_cast<std::size_t>(k)];
                for (int i = k + 1; i < n; ++i) v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += t * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            }
        }
        for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = 0.0;
        v[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)] = 1.0;
    }
    const double tiny = 1.6033346880071782e-291;
    const double eps = 2.220446049250313e-16;
    int p = n;
    while (p > 0) {
        int k = p - 2;
        while (k >= 0) {
            const double threshold = tiny + eps * (std::abs(singular[static_cast<std::size_t>(k)]) + std::abs(singular[static_cast<std::size_t>(k + 1)]));
            if (!(std::abs(e[static_cast<std::size_t>(k)]) > threshold)) { e[static_cast<std::size_t>(k)] = 0.0; break; }
            --k;
        }
        int kase = 0;
        if (k == p - 2) kase = 4;
        else {
            int ks = p - 1;
            while (ks >= k) {
                if (ks == k) break;
                const double t = (ks != p - 1 ? std::abs(e[static_cast<std::size_t>(ks)]) : 0.0) +
                                 (ks != k + 1 ? std::abs(e[static_cast<std::size_t>(ks - 1)]) : 0.0);
                if (std::abs(singular[static_cast<std::size_t>(ks)]) <= tiny + eps * t) { singular[static_cast<std::size_t>(ks)] = 0.0; break; }
                --ks;
            }
            if (ks == k) kase = 3;
            else if (ks == p - 1) kase = 1;
            else { kase = 2; k = ks; }
        }
        ++k;
        if (kase == 1) {
            double f = e[static_cast<std::size_t>(p - 2)]; e[static_cast<std::size_t>(p - 2)] = 0.0;
            for (int j = p - 2; j >= k; --j) {
                double t = fast_hypot(singular[static_cast<std::size_t>(j)], f);
                const double cs = singular[static_cast<std::size_t>(j)] / t, sn = f / t;
                singular[static_cast<std::size_t>(j)] = t;
                if (j != k) { f = -sn * e[static_cast<std::size_t>(j - 1)]; e[static_cast<std::size_t>(j - 1)] = cs * e[static_cast<std::size_t>(j - 1)]; }
                for (int i = 0; i < n; ++i) { t = cs * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + sn * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(p - 1)]; v[static_cast<std::size_t>(i)][static_cast<std::size_t>(p - 1)] = -sn * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + cs * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(p - 1)]; v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = t; }
            }
        } else if (kase == 2) {
            double f = e[static_cast<std::size_t>(k - 1)]; e[static_cast<std::size_t>(k - 1)] = 0.0;
            for (int j = k; j < p; ++j) {
                double t = fast_hypot(singular[static_cast<std::size_t>(j)], f);
                const double cs = singular[static_cast<std::size_t>(j)] / t, sn = f / t;
                singular[static_cast<std::size_t>(j)] = t; f = -sn * e[static_cast<std::size_t>(j)]; e[static_cast<std::size_t>(j)] = cs * e[static_cast<std::size_t>(j)];
                for (int i = 0; i < m; ++i) { t = cs * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + sn * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k - 1)]; u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k - 1)] = -sn * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + cs * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k - 1)]; u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = t; }
            }
        } else if (kase == 3) {
            const double scale = std::max({std::abs(singular[static_cast<std::size_t>(p - 1)]), std::abs(singular[static_cast<std::size_t>(p - 2)]), std::abs(e[static_cast<std::size_t>(p - 2)]), std::abs(singular[static_cast<std::size_t>(k)]), std::abs(e[static_cast<std::size_t>(k)])});
            const double sp = singular[static_cast<std::size_t>(p - 1)] / scale, spm1 = singular[static_cast<std::size_t>(p - 2)] / scale, epm1 = e[static_cast<std::size_t>(p - 2)] / scale, sk = singular[static_cast<std::size_t>(k)] / scale, ek = e[static_cast<std::size_t>(k)] / scale;
            const double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
            const double c = (sp * epm1) * (sp * epm1);
            double shift = 0.0;
            if (b != 0.0 || c != 0.0) { double root = std::sqrt(b * b + c); if (b < 0.0) root = -root; shift = c / (b + root); }
            double f = (sk + sp) * (sk - sp) + shift, g = sk * ek;
            for (int j = k; j < p - 1; ++j) {
                double t = fast_hypot(f, g); const double cs = f / t, sn = g / t;
                if (j != k) e[static_cast<std::size_t>(j - 1)] = t;
                f = cs * singular[static_cast<std::size_t>(j)] + sn * e[static_cast<std::size_t>(j)]; e[static_cast<std::size_t>(j)] = cs * e[static_cast<std::size_t>(j)] - sn * singular[static_cast<std::size_t>(j)]; g = sn * singular[static_cast<std::size_t>(j + 1)]; singular[static_cast<std::size_t>(j + 1)] = cs * singular[static_cast<std::size_t>(j + 1)];
                for (int i = 0; i < n; ++i) { t = cs * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + sn * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)]; v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)] = -sn * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + cs * v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)]; v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = t; }
                t = fast_hypot(f, g); const double cs2 = f / t, sn2 = g / t;
                singular[static_cast<std::size_t>(j)] = t; f = cs2 * e[static_cast<std::size_t>(j)] + sn2 * singular[static_cast<std::size_t>(j + 1)]; singular[static_cast<std::size_t>(j + 1)] = -sn2 * e[static_cast<std::size_t>(j)] + cs2 * singular[static_cast<std::size_t>(j + 1)]; g = sn2 * e[static_cast<std::size_t>(j + 1)]; e[static_cast<std::size_t>(j + 1)] = cs2 * e[static_cast<std::size_t>(j + 1)];
                if (j < m - 1) for (int i = 0; i < m; ++i) { t = cs2 * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + sn2 * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)]; u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)] = -sn2 * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] + cs2 * u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)]; u[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = t; }
            }
            e[static_cast<std::size_t>(p - 2)] = f;
        } else {
            if (singular[static_cast<std::size_t>(k)] <= 0.0) {
                if (singular[static_cast<std::size_t>(k)] < 0.0) singular[static_cast<std::size_t>(k)] = -singular[static_cast<std::size_t>(k)]; else singular[static_cast<std::size_t>(k)] = 0.0;
                for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = -v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            }
            while (k < p - 1 && singular[static_cast<std::size_t>(k)] < singular[static_cast<std::size_t>(k + 1)]) {
                std::swap(singular[static_cast<std::size_t>(k)], singular[static_cast<std::size_t>(k + 1)]);
                for (int i = 0; i < n; ++i) std::swap(v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)], v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k + 1)]);
                for (int i = 0; i < m; ++i) std::swap(u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)], u[static_cast<std::size_t>(i)][static_cast<std::size_t>(k + 1)]);
                ++k;
            }
            --p;
        }
    }
    return {std::move(singular), std::move(u)};
}

Matrix reduced_observation_matrix(const std::vector<double>& data,
                                  std::size_t dimension,
                                  SegmenterTelemetry* telemetry) {
    const std::size_t n = data.size();
    const std::size_t p = std::min(dimension, n);
    if (p == 0) return {};
    JavaRandom random(RANDOM_SEED);
    std::vector<double> subsample;
    subsample.reserve(p);
    if (p == n) subsample = data;
    else for (std::size_t i = 0; i < p; ++i) subsample.push_back(data[random.next_int(static_cast<std::int32_t>(n))]);

    const Matrix sub_kernel = symmetric_kernel_matrix(subsample, telemetry);
    const auto eigen = commons_math_svd(sub_kernel);
    Kokkos::View<double*> subsample_view("contamination_subsample_projection", p);
    auto host_subsample_view = Kokkos::create_mirror_view(subsample_view);
    for (std::size_t i = 0; i < p; ++i) host_subsample_view(i) = subsample[i];
    Kokkos::deep_copy(subsample_view, host_subsample_view);
    Kokkos::View<double*> data_view("contamination_data", n);
    Kokkos::View<double*> singular_view("contamination_singular_values", p);
    Kokkos::View<double**> u_view("contamination_svd_u", p, p);
    auto host_data = Kokkos::create_mirror_view(data_view);
    auto host_singular = Kokkos::create_mirror_view(singular_view);
    auto host_u = Kokkos::create_mirror_view(u_view);
    for (std::size_t i = 0; i < n; ++i) host_data(i) = data[i];
    for (std::size_t i = 0; i < p; ++i) {
        host_singular(i) = eigen.values[i];
        for (std::size_t j = 0; j < p; ++j) host_u(i, j) = eigen.vectors[i][j];
    }
    Kokkos::deep_copy(data_view, host_data);
    Kokkos::deep_copy(singular_view, host_singular);
    Kokkos::deep_copy(u_view, host_u);
    Kokkos::View<double**> reduced_view("contamination_reduced_kernel", n, p);
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("contamination-segmenter-reduced-v1");
    host_batch.records = n * p;
    host_batch.bytes = n * sizeof(double) + p * sizeof(double) +
                       p * sizeof(double) + p * p * sizeof(double) + n * p * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("contamination-segmenter-reduced");
    kernel_plan.begin_prepare(host_batch);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(n * p);
    device_batch.bind("data", data_view);
    device_batch.bind("subsample", subsample_view);
    device_batch.bind("singular", singular_view);
    device_batch.bind("u", u_view);
    device_batch.bind("reduced", reduced_view);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    Kokkos::parallel_for(
        "contamination_reduced_kernel_fill",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0, 0}, {n, p}),
        KOKKOS_LAMBDA(const std::size_t row, const std::size_t column) {
            double value = 0.0;
            const double scale = 1.0 / (Kokkos::sqrt(singular_view(column)) + EPSILON);
            for (std::size_t i = 0; i < p; ++i)
                value += kernel(data_view(row), subsample_view(i)) * u_view(i, column) * scale;
            reduced_view(row, column) = value;
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    record_segmenter_plan(telemetry, kernel_plan.telemetry(), n * p);
    auto copied = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, reduced_view);
    Matrix result(n, std::vector<double>(p));
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < p; ++column)
            result[row][column] = copied(row, column);
    return result;
}

std::vector<double> approximation_diagonal(const Matrix& reduced) {
    std::vector<double> result(reduced.size(), 0.0);
    for (std::size_t row = 0; row < reduced.size(); ++row) {
        double sum = 0.0;
        for (const double value : reduced[row]) sum += value * value;
        result[row] = std::sqrt(sum) * std::sqrt(sum);
    }
    return result;
}

struct Cost { double d; std::vector<double> w; double v; double c; };

Cost segment_cost(std::size_t start, std::size_t end,
                  const Matrix& reduced, const std::vector<double>& diagonal) {
    const std::size_t n = reduced.size();
    const std::size_t p = reduced.empty() ? 0 : reduced.front().size();
    double d = diagonal[start];
    std::vector<double> w = reduced[start];
    double v = 0.0;
    for (const double value : w) v += value * value;
    auto add = [&](std::size_t index) {
        d += diagonal[index];
        double dot = 0.0;
        for (std::size_t j = 0; j < p; ++j) {
            dot += reduced[index][j] * w[j];
            w[j] += reduced[index][j];
        }
        v += 2.0 * dot + diagonal[index];
    };
    if (start <= end) for (std::size_t index = start + 1; index <= end; ++index) add(index);
    else {
        for (std::size_t index = start + 1; index < n; ++index) add(index);
        for (std::size_t index = 0; index <= end; ++index) add(index);
    }
    return {d, std::move(w), v, d - v / (start <= end ? end - start + 1.0 : n - start + end + 1.0)};
}

std::vector<double> window_costs(const Matrix& reduced,
                                 const std::vector<double>& diagonal,
                                 std::size_t window_size) {
    const std::size_t n = reduced.size();
    const std::size_t p = reduced.empty() ? 0 : reduced.front().size();
    std::size_t center = 0;
    std::size_t start = (center + n - window_size + 1) % n;
    std::size_t end = (center + window_size) % n;
    auto left = segment_cost(start, center, reduced, diagonal);
    auto right = segment_cost(center + 1, end, reduced, diagonal);
    auto total = segment_cost(start, end, reduced, diagonal);
    std::vector<double> costs(n);
    costs[center] = left.c + right.c - total.c;
    const double reciprocal = 1.0 / static_cast<double>(window_size);
    for (std::size_t iteration = 0; iteration < n; ++iteration) {
        center = iteration;
        const std::size_t center_next = (center + 1) % n;
        const std::size_t end_next = (end + 1) % n;
        auto update = [&](double& d, std::vector<double>& w, double& v,
                          std::size_t remove, std::size_t add_index) {
            d -= diagonal[remove];
            double dot = 0.0;
            for (std::size_t j = 0; j < p; ++j) {
                dot += reduced[remove][j] * w[j];
                w[j] -= reduced[remove][j];
            }
            v += -2.0 * dot + diagonal[remove];
            d += diagonal[add_index];
            dot = 0.0;
            for (std::size_t j = 0; j < p; ++j) {
                dot += reduced[add_index][j] * w[j];
                w[j] += reduced[add_index][j];
            }
            v += 2.0 * dot + diagonal[add_index];
        };
        update(left.d, left.w, left.v, start, center_next);
        left.c = left.d - left.v * reciprocal;
        update(right.d, right.w, right.v, center_next, end_next);
        right.c = right.d - right.v * reciprocal;
        update(total.d, total.w, total.v, start, end_next);
        total.c = total.d - 0.5 * total.v * reciprocal;
        costs[center_next] = left.c + right.c - total.c;
        start = (start + 1) % n;
        end = end_next;
    }
    return costs;
}

// Java's PersistenceOptimizer is a watershed over a one-dimensional cost
// curve.  This implementation keeps its stable value/index ordering and
// component merge rule, including plateau minima.
struct Component { std::size_t left, right, minimum; double value; };
struct Pair { std::size_t minimum; double persistence; };

bool java_less(double first, double second) {
    if (first < second) return true;
    if (first > second) return false;
    if (std::isnan(first) || std::isnan(second)) return std::isnan(second) && !std::isnan(first);
    return std::signbit(first) != std::signbit(second) ? std::signbit(first) : false;
}

std::vector<std::size_t> persistence_minima(const std::vector<double>& data) {
    if (data.empty()) return {};
    std::vector<std::size_t> order(data.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return java_less(data[a], data[b]); });
    if (data.size() == 1) return {0};
    std::vector<Component> components;
    std::vector<std::ptrdiff_t> colors(data.size(), -1);
    std::vector<Pair> pairs;
    auto create = [&](std::size_t index) {
        colors[index] = static_cast<std::ptrdiff_t>(components.size());
        components.push_back({index, index, index, data[index]});
    };
    auto extend = [&](std::size_t color, std::size_t index) {
        if (index + 1 == components[color].left) components[color].left = index;
        else if (components[color].right + 1 == index) components[color].right = index;
        colors[index] = static_cast<std::ptrdiff_t>(color);
    };
    auto merge = [&](std::size_t left_color, std::size_t right_color, std::size_t index) {
        std::size_t keep = left_color, remove = right_color;
        if (java_less(components[right_color].value, components[left_color].value)) {
            keep = right_color; remove = left_color;
        }
        const auto old_left = components[remove].left;
        const auto old_right = components[remove].right;
        const auto old_min = components[remove].minimum;
        colors[old_left] = static_cast<std::ptrdiff_t>(keep);
        colors[old_right] = static_cast<std::ptrdiff_t>(keep);
        if (components[keep].minimum > old_min) components[keep].left = old_left;
        else components[keep].right = old_right;
        colors[index] = colors[index - 1];
    };
    for (const auto index : order) {
        if (index == 0) {
            if (colors[1] < 0) create(index); else extend(static_cast<std::size_t>(colors[1]), index);
        } else if (index + 1 == data.size()) {
            if (colors[index - 1] < 0) create(index); else extend(static_cast<std::size_t>(colors[index - 1]), index);
        } else {
            const auto left = colors[index - 1];
            const auto right = colors[index + 1];
            if (left < 0 && right < 0) create(index);
            else if (left >= 0 && right < 0) extend(static_cast<std::size_t>(left), index);
            else if (left < 0 && right >= 0) extend(static_cast<std::size_t>(right), index);
            else {
                const auto left_color = static_cast<std::size_t>(left);
                const auto right_color = static_cast<std::size_t>(right);
                // PersistenceOptimizer's ExtremaPair breaks equal minima by
                // the lower data index (and mergeComponents keeps the lower
                // component index).  Do not let a floating-point tie move the
                // watershed minimum to the right-hand component.
                const auto minimum = java_less(components[right_color].value,
                                               components[left_color].value)
                    ? components[left_color].minimum
                    : java_less(components[left_color].value,
                                components[right_color].value)
                        ? components[right_color].minimum
                        : std::min(components[left_color].minimum,
                                   components[right_color].minimum);
                const auto maximum = index;
                pairs.push_back({minimum, data[maximum] - data[minimum]});
                merge(left_color, right_color, index);
            }
        }
    }
    std::stable_sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
        return java_less(b.persistence, a.persistence);
    });
    std::vector<std::size_t> minima{order.front()};
    for (const auto& pair : pairs) minima.push_back(pair.minimum);
    return minima;
}

std::vector<std::size_t> select_changepoints(const std::vector<std::size_t>& candidates,
                                             std::size_t max_num,
                                             double linear_factor, double log_factor,
                                             const Matrix& reduced,
                                             const std::vector<double>& diagonal) {
    const std::size_t n = reduced.size();
    std::vector<std::size_t> sorted = candidates;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::vector<std::size_t> starts{0};
    for (const auto point : sorted) starts.push_back(std::min(point + 1, n - 1));
    std::vector<std::size_t> ends = sorted;
    ends.push_back(n - 1);
    struct Segment { std::size_t start, end; double cost; };
    std::vector<Segment> segments;
    for (std::size_t i = 0; i < starts.size(); ++i)
        segments.push_back({starts[i], ends[i], segment_cost(starts[i], ends[i], reduced, diagonal).c});
    const std::size_t num_segments = segments.size();
    std::vector<double> total_costs{0.0};
    for (const auto& segment : segments) total_costs.front() += segment.cost;
    std::vector<double> pair_costs, merged_costs, merging_costs;
    for (std::size_t i = 0; i + 1 < segments.size(); ++i) {
        pair_costs.push_back(segments[i].cost + segments[i + 1].cost);
        merged_costs.push_back(segment_cost(segments[i].start, segments[i + 1].end, reduced, diagonal).c);
        merging_costs.push_back(pair_costs.back() - merged_costs.back());
    }
    std::vector<std::size_t> changepoints;
    for (std::size_t iteration = 0; iteration + 1 < num_segments; ++iteration) {
        const auto index = static_cast<std::size_t>(std::distance(merging_costs.begin(),
            std::max_element(merging_costs.begin(), merging_costs.end())));
        const auto new_start = segments[index].start;
        const auto mergepoint = segments[index].end;
        const auto new_end = segments[index + 1].end;
        const auto new_cost = merged_costs[index];
        segments.erase(segments.begin() + static_cast<std::ptrdiff_t>(index), segments.begin() + static_cast<std::ptrdiff_t>(index + 2));
        segments.insert(segments.begin() + static_cast<std::ptrdiff_t>(index), {new_start, new_end, new_cost});
        pair_costs.erase(pair_costs.begin() + static_cast<std::ptrdiff_t>(index));
        merged_costs.erase(merged_costs.begin() + static_cast<std::ptrdiff_t>(index));
        merging_costs.erase(merging_costs.begin() + static_cast<std::ptrdiff_t>(index));
        if (index > 0) {
            pair_costs[index - 1] = segments[index - 1].cost + segments[index].cost;
            merged_costs[index - 1] = segment_cost(segments[index - 1].start, new_end, reduced, diagonal).c;
            merging_costs[index - 1] = pair_costs[index - 1] - merged_costs[index - 1];
        }
        if (index + 1 < segments.size()) {
            pair_costs[index] = segments[index].cost + segments[index + 1].cost;
            merged_costs[index] = segment_cost(new_start, segments[index + 1].end, reduced, diagonal).c;
            merging_costs[index] = pair_costs[index] - merged_costs[index];
        }
        double total = 0.0;
        for (const auto& segment : segments) total += segment.cost;
        total_costs.insert(total_costs.begin(), total);
        changepoints.insert(changepoints.begin(), mergepoint);
    }
    const std::size_t effective = std::min(max_num, changepoints.size());
    std::size_t optimal = 0;
    double minimum = std::numeric_limits<double>::infinity();
    for (std::size_t count = 0; count <= effective; ++count) {
        const double penalty = linear_factor * static_cast<double>(count) +
            log_factor * static_cast<double>(count) *
            std::log(static_cast<double>(n) / (static_cast<double>(count) + EPSILON));
        const double value = total_costs[count] + penalty;
        if (value < minimum) { minimum = value; optimal = count; }
    }
    changepoints.resize(optimal);
    return changepoints;
}

}  // namespace

std::vector<std::size_t> find_changepoints_with_windows(
    const std::vector<double>& data,
    std::size_t max_num_changepoints,
    std::size_t kernel_approximation_dimension,
    const std::vector<std::size_t>& window_sizes,
    double linear_penalty,
    double log_linear_penalty,
    SegmenterTelemetry* telemetry) {
    if (data.empty() || max_num_changepoints == 0) return {};
    const Matrix reduced = reduced_observation_matrix(data, kernel_approximation_dimension, telemetry);
    const auto diagonal = approximation_diagonal(reduced);
    std::vector<std::size_t> candidates;
    for (const auto window_size : window_sizes) {
        if (window_size == 0 || 2 * window_size > data.size()) continue;
        const auto costs = window_costs(reduced, diagonal, window_size);
        const auto minima = persistence_minima(costs);
        for (const auto index : minima)
            if (index != 0 && index + 1 != data.size() &&
                candidates.size() < max_num_changepoints)
                candidates.push_back(index);
    }
    auto selected = select_changepoints(candidates, max_num_changepoints,
                                        linear_penalty, log_linear_penalty,
                                        reduced, diagonal);
    // ContaminationSegmenter requests INDEX order; BACKWARD_SELECTION is an
    // internal ranking only and must not become an interval-order bug.
    std::sort(selected.begin(), selected.end());
    return selected;
}

std::vector<std::size_t> find_changepoints(const std::vector<double>& data,
                                           std::size_t max_num_changepoints,
                                           std::size_t kernel_approximation_dimension,
                                           std::size_t window_size,
                                           double linear_penalty,
                                           double log_linear_penalty,
                                           SegmenterTelemetry* telemetry) {
    return find_changepoints_with_windows(data, max_num_changepoints,
                                          kernel_approximation_dimension,
                                          {window_size}, linear_penalty,
                                          log_linear_penalty, telemetry);
}

}  // namespace fastgatk::contamination
