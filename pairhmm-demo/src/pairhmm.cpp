#include "pairhmm.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <thread>

#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#include <immintrin.h>
#define FASTGATK_X86_SIMD 1
#else
#define FASTGATK_X86_SIMD 0
#endif

namespace fastgatk::pairhmm {
namespace {

constexpr int kAvx2Lanes = 4;  // AVX2: four independent double lanes.
constexpr int kAvx512Lanes = 8;  // AVX-512: eight independent double lanes.

#if defined(__GNUC__) || defined(__clang__)
#define FASTGATK_TARGET_AVX2 __attribute__((target("avx2")))
#define FASTGATK_TARGET_AVX512 __attribute__((target("avx512f")))
#else
#define FASTGATK_TARGET_AVX2
#define FASTGATK_TARGET_AVX512
#endif
// Mirror GATK's `Math.log10(INITIAL_CONDITION)` expression, rather than algebraically
// rewriting it; the last few bits are part of the compatibility contract.
const double kInitialConditionLog10 = std::log10(std::pow(2.0, 1020.0));
constexpr double kTristateCorrection = 3.0;

// The strict contract follows the ordinary Java floating-point environment:
// round-to-nearest, with subnormals preserved. GKL's fast path intentionally
// uses FTZ, so this is only applied by the scalar compatibility backend.
inline void configure_strict_fp_environment() {
#if defined(__x86_64__) || defined(__i386__)
    unsigned csr = _mm_getcsr();
    csr &= ~(_MM_FLUSH_ZERO_MASK | _MM_DENORMALS_ZERO_MASK | _MM_ROUND_MASK);
    csr |= _MM_ROUND_NEAREST;
    _mm_setcsr(csr);
#endif
}

inline double approximate_log10_sum(double a, double b);

struct GatkTables {
    std::array<double, 256> quality{};
    std::array<double, 256 * 256> match{};
};

const GatkTables& gatk_tables() {
    static const GatkTables table = [] {
        GatkTables values{};
        for (std::size_t q = 0; q < values.quality.size(); ++q) {
            values.quality[q] = std::pow(10.0, -0.1 * static_cast<double>(q));
        }
        for (std::size_t i = 0; i < 256; ++i) {
            for (std::size_t d = 0; d < 256; ++d) {
                values.match[i * 256 + d] = 1.0 - std::pow(10.0,
                    approximate_log10_sum(-0.1 * static_cast<double>(i),
                                          -0.1 * static_cast<double>(d)));
            }
        }
        // A Java-generated table can be supplied for strict GATK compatibility. The
        // text format is: 256 lines "q <raw-bits>" then 65536 lines
        // "m <ins> <del> <raw-bits>". Production builds can embed these values;
        // the file hook keeps the demo auditable without copying Java into the kernel.
        if (const char* path = std::getenv("FAST_GATK_PAIRHMM_TABLES")) {
            std::ifstream in(path);
            std::string kind;
            std::uint64_t bits = 0;
            bool ok = static_cast<bool>(in);
            for (std::size_t q = 0; ok && q < 256; ++q) {
                in >> kind >> std::hex >> bits >> std::dec;
                ok = in && kind == "q";
                if (ok) std::memcpy(&values.quality[q], &bits, sizeof(bits));
            }
            for (std::size_t i = 0; ok && i < 256; ++i) {
                for (std::size_t d = 0; ok && d < 256; ++d) {
                    std::size_t file_i = 0, file_d = 0;
                    in >> kind >> file_i >> file_d >> std::hex >> bits >> std::dec;
                    ok = in && kind == "m" && file_i == i && file_d == d;
                    if (ok) std::memcpy(&values.match[i * 256 + d], &bits, sizeof(bits));
                }
            }
            if (!ok) throw std::runtime_error("invalid FAST_GATK_PAIRHMM_TABLES file");
        }
        return values;
    }();
    return table;
}

const std::array<double, 256>& quality_error_table() { return gatk_tables().quality; }
const std::array<double, 256 * 256>& match_to_match_table() { return gatk_tables().match; }

inline void validate(const PairInput& p) {
    const std::size_t n = p.read.size();
    if (p.haplotype.empty() || n == 0 ||
        p.read_qual.size() != n || p.insertion_gop.size() != n ||
        p.deletion_gop.size() != n || p.gap_continuation.size() != n) {
        throw std::invalid_argument("PairHMM input lengths are inconsistent or empty");
    }
}

inline double quality_error(std::uint8_t q) {
    return quality_error_table()[q];
}

// GATK's PairHMMModel uses a 0.0001-quantized Jacobian log table when building
// match-to-match transitions. Keeping this approximation (rather than replacing it
// with an exact log-sum) is important for matching LoglessPairHMM.
inline double approximate_log10_sum(double a, double b) {
    if (a > b) std::swap(a, b);
    if (std::isinf(a) && a < 0.0) return b;
    const double diff = b - a;
    if (diff >= 8.0) return b;
    const int index = static_cast<int>(diff * 10000.0 + 0.5);
    const double quantized = static_cast<double>(index) * 0.0001;
    return b + std::log10(1.0 + std::pow(10.0, -quantized));
}

inline double match_to_match(std::uint8_t insertion, std::uint8_t deletion) {
    return match_to_match_table()[static_cast<std::size_t>(insertion) * 256 + deletion];
}

inline double prior(std::uint8_t read_base, std::uint8_t hap_base, std::uint8_t qual,
                    bool tristate) {
    const double error = quality_error(qual);
    if (read_base == hap_base || read_base == 'N' || hap_base == 'N') return 1.0 - error;
    return error / (tristate ? kTristateCorrection : 1.0);
}

double scalar_once(const PairInput& p, bool tristate) {
    validate(p);
    const std::size_t r = p.read.size();
    const std::size_t h = p.haplotype.size();
    const std::size_t stride = h + 1;
    const std::size_t cells = (r + 1) * stride;
    std::vector<double> m(cells, 0.0), ins(cells, 0.0), del(cells, 0.0);

    const double initial = std::pow(2.0, 1020.0) / static_cast<double>(h);
    for (std::size_t j = 0; j <= h; ++j) del[j] = initial;

    for (std::size_t i = 1; i <= r; ++i) {
        const std::size_t row = i * stride;
        const std::size_t previous = (i - 1) * stride;
        const double mm = match_to_match(p.insertion_gop[i - 1], p.deletion_gop[i - 1]);
        const double im = 1.0 - quality_error(p.gap_continuation[i - 1]);
        const double mi = quality_error(p.insertion_gop[i - 1]);
        const double ii = quality_error(p.gap_continuation[i - 1]);
        const double md = quality_error(p.deletion_gop[i - 1]);
        const double dd = ii;
        for (std::size_t j = 1; j <= h; ++j) {
            const double p0 = prior(p.read[i - 1], p.haplotype[j - 1], p.read_qual[i - 1], tristate);
            m[row + j] = p0 * (m[previous + j - 1] * mm + ins[previous + j - 1] * im +
                               del[previous + j - 1] * im);
            ins[row + j] = m[previous + j] * mi + ins[previous + j] * ii;
            del[row + j] = m[row + j - 1] * md + del[row + j - 1] * dd;
        }
    }

    double sum = 0.0;
    const std::size_t row = r * stride;
    for (std::size_t j = 1; j <= h; ++j) sum += m[row + j] + ins[row + j];
    return std::log10(sum) - kInitialConditionLog10;
}

#if FASTGATK_X86_SIMD

struct alignas(32) Lane4 {
    double v[kAvx2Lanes];
};

struct ScalarWorkspace {
    std::vector<double> m, ins, del;
    std::size_t stride = 0;
    void ensure(std::size_t read_length, std::size_t hap_length) {
        stride = hap_length + 1;
        const std::size_t cells = (read_length + 1) * stride;
        if (m.size() != cells) {
            m.resize(cells); ins.resize(cells); del.resize(cells);
        }
        std::fill(m.begin(), m.end(), 0.0);
        std::fill(ins.begin(), ins.end(), 0.0);
        std::fill(del.begin(), del.end(), 0.0);
    }
};

double scalar_once(const PairInput& p, bool tristate, ScalarWorkspace& ws) {
    validate(p);
    const std::size_t r = p.read.size();
    const std::size_t h = p.haplotype.size();
    ws.ensure(r, h);
    const std::size_t stride = ws.stride;
    const double initial = std::pow(2.0, 1020.0) / static_cast<double>(h);
    for (std::size_t j = 0; j <= h; ++j) ws.del[j] = initial;

    for (std::size_t i = 1; i <= r; ++i) {
        const std::size_t row = i * stride, previous = (i - 1) * stride;
        const double mm = match_to_match(p.insertion_gop[i - 1], p.deletion_gop[i - 1]);
        const double im = 1.0 - quality_error(p.gap_continuation[i - 1]);
        const double mi = quality_error(p.insertion_gop[i - 1]);
        const double ii = quality_error(p.gap_continuation[i - 1]);
        const double md = quality_error(p.deletion_gop[i - 1]);
        const double dd = ii;
        for (std::size_t j = 1; j <= h; ++j) {
            const double p0 = prior(p.read[i - 1], p.haplotype[j - 1], p.read_qual[i - 1], tristate);
            ws.m[row + j] = p0 * (ws.m[previous + j - 1] * mm + ws.ins[previous + j - 1] * im + ws.del[previous + j - 1] * im);
            ws.ins[row + j] = ws.m[previous + j] * mi + ws.ins[previous + j] * ii;
            ws.del[row + j] = ws.m[row + j - 1] * md + ws.del[row + j - 1] * dd;
        }
    }
    double sum = 0.0;
    const std::size_t row = r * stride;
    for (std::size_t j = 1; j <= h; ++j) sum += ws.m[row + j] + ws.ins[row + j];
    return std::log10(sum) - kInitialConditionLog10;
}

FASTGATK_TARGET_AVX2 inline __m256d load(const Lane4& x) { return _mm256_load_pd(x.v); }
FASTGATK_TARGET_AVX2 inline void store(Lane4& x, __m256d v) { _mm256_store_pd(x.v, v); }

struct SimdWorkspace {
    std::vector<Lane4> m, ins, del;
    std::size_t stride = 0;
    void ensure(std::size_t read_length, std::size_t hap_length) {
        (void)read_length;
        stride = hap_length + 1;
        // PairHMM only consumes the previous and current rows. Keeping the full
        // matrix made an AVX-512 worker use ~4.7 MiB for 150x160 input and caused
        // severe shared-cache pressure at high thread counts.
        const std::size_t cells = 2 * stride;
        if (m.size() != cells) { m.resize(cells); ins.resize(cells); del.resize(cells); }
        std::fill(m.begin(), m.end(), Lane4{});
        std::fill(ins.begin(), ins.end(), Lane4{});
        std::fill(del.begin(), del.end(), Lane4{});
    }
};

FASTGATK_TARGET_AVX2
void simd_batch_once_avx2(const PairInput* p, std::size_t count, double* output, bool tristate, SimdWorkspace& ws) {
    const std::size_t max_r = [&] { std::size_t n = 0; for (std::size_t i = 0; i < count; ++i) n = std::max(n, p[i].read.size()); return n; }();
    const std::size_t max_h = [&] { std::size_t n = 0; for (std::size_t i = 0; i < count; ++i) n = std::max(n, p[i].haplotype.size()); return n; }();
    const std::size_t stride = max_h + 1;
    ws.ensure(max_r, max_h);
    auto& m = ws.m; auto& ins = ws.ins; auto& del = ws.del;
    Lane4 initial{};
    for (int lane = 0; lane < kAvx2Lanes; ++lane) {
        if (lane < static_cast<int>(count))
            initial.v[lane] = std::pow(2.0, 1020.0) / static_cast<double>(p[lane].haplotype.size());
    }
    for (std::size_t j = 0; j <= max_h; ++j) {
        Lane4 x{};
        for (int lane = 0; lane < kAvx2Lanes; ++lane) {
            x.v[lane] = (lane < static_cast<int>(count) && j <= p[lane].haplotype.size())
                ? initial.v[lane] : 0.0;
        }
        del[j] = x;
    }

    for (std::size_t i = 1; i <= max_r; ++i) {
        alignas(32) double mm_a[kAvx2Lanes]{}, im_a[kAvx2Lanes]{}, mi_a[kAvx2Lanes]{}, ii_a[kAvx2Lanes]{}, md_a[kAvx2Lanes]{}, dd_a[kAvx2Lanes]{};
        alignas(32) double prior_match_a[kAvx2Lanes]{}, prior_mismatch_a[kAvx2Lanes]{};
        for (int lane = 0; lane < kAvx2Lanes; ++lane) {
            if (lane >= static_cast<int>(count) || i > p[lane].read.size()) continue;
            const auto& q = p[lane];
            mm_a[lane] = match_to_match(q.insertion_gop[i - 1], q.deletion_gop[i - 1]);
            im_a[lane] = 1.0 - quality_error(q.gap_continuation[i - 1]);
            mi_a[lane] = quality_error(q.insertion_gop[i - 1]);
            ii_a[lane] = quality_error(q.gap_continuation[i - 1]);
            md_a[lane] = quality_error(q.deletion_gop[i - 1]);
            dd_a[lane] = ii_a[lane];
            const double error = quality_error(q.read_qual[i - 1]);
            prior_match_a[lane] = 1.0 - error;
            prior_mismatch_a[lane] = tristate ? error / kTristateCorrection : error;
        }
        const __m256d mm = _mm256_load_pd(mm_a), im = _mm256_load_pd(im_a), mi = _mm256_load_pd(mi_a);
        const __m256d ii = _mm256_load_pd(ii_a), md = _mm256_load_pd(md_a), dd = _mm256_load_pd(dd_a);
        const std::size_t row = (i & 1U) * stride, previous = ((i - 1) & 1U) * stride;
        m[row] = Lane4{}; ins[row] = Lane4{}; del[row] = Lane4{};
        for (std::size_t j = 1; j <= max_h; ++j) {
            alignas(32) double prior_a[kAvx2Lanes]{};
            for (int lane = 0; lane < kAvx2Lanes; ++lane) {
                if (lane < static_cast<int>(count) && i <= p[lane].read.size() && j <= p[lane].haplotype.size()) {
                    const bool matches = p[lane].read[i - 1] == p[lane].haplotype[j - 1] ||
                                         p[lane].read[i - 1] == 'N' || p[lane].haplotype[j - 1] == 'N';
                    prior_a[lane] = matches ? prior_match_a[lane] : prior_mismatch_a[lane];
                }
            }
            const __m256d pv = _mm256_load_pd(prior_a);
            // Match the scalar Java/C++ left-associative order: (A + B) + C.
            // Writing B + C first changes one-ulp rounding and breaks strict SIMD conformance.
            const __m256d mv_ab = _mm256_add_pd(
                _mm256_mul_pd(load(m[previous + j - 1]), mm),
                _mm256_mul_pd(load(ins[previous + j - 1]), im));
            const __m256d mv = _mm256_add_pd(
                mv_ab, _mm256_mul_pd(load(del[previous + j - 1]), im));
            const __m256d mcur = _mm256_mul_pd(pv, mv);
            const __m256d icur = _mm256_add_pd(_mm256_mul_pd(load(m[previous + j]), mi),
                                               _mm256_mul_pd(load(ins[previous + j]), ii));
            const __m256d dcur = _mm256_add_pd(_mm256_mul_pd(load(m[row + j - 1]), md),
                                               _mm256_mul_pd(load(del[row + j - 1]), dd));
            store(m[row + j], mcur); store(ins[row + j], icur); store(del[row + j], dcur);
        }
    }

    const std::size_t row = (max_r & 1U) * stride;
    for (std::size_t lane = 0; lane < count; ++lane) {
        // For mixed-length batches the last row includes zero-padded cells. Re-run only
        // those lanes through the reference; equal-length benchmark batches stay SIMD-only.
        if (p[lane].read.size() != max_r || p[lane].haplotype.size() != max_h) {
            output[lane] = scalar_once(p[lane], tristate);
        } else {
            // Keep the Java final accumulation order explicitly scalar. The matrix
            // recurrence is SIMD; this tiny tail avoids a compiler/vector reduction
            // changing the last bit of a GATK likelihood.
            double sum = 0.0;
            for (std::size_t j = 1; j <= max_h; ++j) sum += m[row + j].v[lane] + ins[row + j].v[lane];
            output[lane] = std::log10(sum) - kInitialConditionLog10;
        }
    }
}

struct alignas(64) Lane8 {
    double v[kAvx512Lanes];
};

struct Simd512Workspace {
    std::vector<Lane8> m, ins, del;
    std::size_t stride = 0;
    void ensure(std::size_t read_length, std::size_t hap_length) {
        (void)read_length;
        stride = hap_length + 1;
        const std::size_t cells = 2 * stride;
        if (m.size() != cells) { m.resize(cells); ins.resize(cells); del.resize(cells); }
        std::fill(m.begin(), m.end(), Lane8{});
        std::fill(ins.begin(), ins.end(), Lane8{});
        std::fill(del.begin(), del.end(), Lane8{});
    }
};

FASTGATK_TARGET_AVX512
void simd_batch_once_avx512(const PairInput* p, std::size_t count, double* output,
                            bool tristate, Simd512Workspace& ws) {
    const std::size_t max_r = [&] {
        std::size_t n = 0;
        for (std::size_t i = 0; i < count; ++i) n = std::max(n, p[i].read.size());
        return n;
    }();
    const std::size_t max_h = [&] {
        std::size_t n = 0;
        for (std::size_t i = 0; i < count; ++i) n = std::max(n, p[i].haplotype.size());
        return n;
    }();
    const std::size_t stride = max_h + 1;
    ws.ensure(max_r, max_h);
    auto& m = ws.m;
    auto& ins = ws.ins;
    auto& del = ws.del;

    Lane8 initial{};
    for (int lane = 0; lane < kAvx512Lanes; ++lane) {
        if (lane < static_cast<int>(count))
            initial.v[lane] = std::pow(2.0, 1020.0) / static_cast<double>(p[lane].haplotype.size());
    }
    for (std::size_t j = 0; j <= max_h; ++j) {
        Lane8 x{};
        for (int lane = 0; lane < kAvx512Lanes; ++lane) {
            x.v[lane] = (lane < static_cast<int>(count) && j <= p[lane].haplotype.size())
                ? initial.v[lane] : 0.0;
        }
        del[j] = x;
    }

    for (std::size_t i = 1; i <= max_r; ++i) {
        alignas(64) double mm_a[kAvx512Lanes]{}, im_a[kAvx512Lanes]{}, mi_a[kAvx512Lanes]{},
            ii_a[kAvx512Lanes]{}, md_a[kAvx512Lanes]{}, dd_a[kAvx512Lanes]{};
        alignas(64) double prior_match_a[kAvx512Lanes]{}, prior_mismatch_a[kAvx512Lanes]{};
        for (int lane = 0; lane < kAvx512Lanes; ++lane) {
            if (lane >= static_cast<int>(count) || i > p[lane].read.size()) continue;
            const auto& q = p[lane];
            mm_a[lane] = match_to_match(q.insertion_gop[i - 1], q.deletion_gop[i - 1]);
            im_a[lane] = 1.0 - quality_error(q.gap_continuation[i - 1]);
            mi_a[lane] = quality_error(q.insertion_gop[i - 1]);
            ii_a[lane] = quality_error(q.gap_continuation[i - 1]);
            md_a[lane] = quality_error(q.deletion_gop[i - 1]);
            dd_a[lane] = ii_a[lane];
            const double error = quality_error(q.read_qual[i - 1]);
            prior_match_a[lane] = 1.0 - error;
            prior_mismatch_a[lane] = tristate ? error / kTristateCorrection : error;
        }
        const __m512d mm = _mm512_load_pd(mm_a), im = _mm512_load_pd(im_a), mi = _mm512_load_pd(mi_a);
        const __m512d ii = _mm512_load_pd(ii_a), md = _mm512_load_pd(md_a), dd = _mm512_load_pd(dd_a);
        const std::size_t row = (i & 1U) * stride, previous = ((i - 1) & 1U) * stride;
        m[row] = Lane8{}; ins[row] = Lane8{}; del[row] = Lane8{};
        for (std::size_t j = 1; j <= max_h; ++j) {
            alignas(64) double prior_a[kAvx512Lanes]{};
            for (int lane = 0; lane < kAvx512Lanes; ++lane) {
                if (lane < static_cast<int>(count) && i <= p[lane].read.size() && j <= p[lane].haplotype.size()) {
                    const bool matches = p[lane].read[i - 1] == p[lane].haplotype[j - 1] ||
                                         p[lane].read[i - 1] == 'N' || p[lane].haplotype[j - 1] == 'N';
                    prior_a[lane] = matches ? prior_match_a[lane] : prior_mismatch_a[lane];
                }
            }
            const __m512d pv = _mm512_load_pd(prior_a);
            // Match the scalar Java/C++ left-associative order: (A + B) + C.
            const __m512d mv_ab = _mm512_add_pd(
                _mm512_mul_pd(_mm512_load_pd(m[previous + j - 1].v), mm),
                _mm512_mul_pd(_mm512_load_pd(ins[previous + j - 1].v), im));
            const __m512d mv = _mm512_add_pd(
                mv_ab, _mm512_mul_pd(_mm512_load_pd(del[previous + j - 1].v), im));
            const __m512d mcur = _mm512_mul_pd(pv, mv);
            const __m512d icur = _mm512_add_pd(
                _mm512_mul_pd(_mm512_load_pd(m[previous + j].v), mi),
                _mm512_mul_pd(_mm512_load_pd(ins[previous + j].v), ii));
            const __m512d dcur = _mm512_add_pd(
                _mm512_mul_pd(_mm512_load_pd(m[row + j - 1].v), md),
                _mm512_mul_pd(_mm512_load_pd(del[row + j - 1].v), dd));
            _mm512_store_pd(m[row + j].v, mcur);
            _mm512_store_pd(ins[row + j].v, icur);
            _mm512_store_pd(del[row + j].v, dcur);
        }
    }

    const std::size_t row = (max_r & 1U) * stride;
    for (std::size_t lane = 0; lane < count; ++lane) {
        if (p[lane].read.size() != max_r || p[lane].haplotype.size() != max_h) {
            output[lane] = scalar_once(p[lane], tristate);
        } else {
            double sum = 0.0;
            for (std::size_t j = 1; j <= max_h; ++j) sum += m[row + j].v[lane] + ins[row + j].v[lane];
            output[lane] = std::log10(sum) - kInitialConditionLog10;
        }
    }
}

#endif  // FASTGATK_X86_SIMD

template <typename Fn>
BatchResult timed_parallel(const std::vector<PairInput>& inputs, int threads, int iterations, Fn fn) {
    if (threads < 1) throw std::invalid_argument("threads must be positive");
    BatchResult result; result.likelihoods.resize(inputs.size());
    const auto start = std::chrono::steady_clock::now();
    for (int it = 0; it < iterations; ++it) {
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        const int actual = std::min<int>(threads, std::max<std::size_t>(1, inputs.size()));
        workers.reserve(actual);
        for (int t = 0; t < actual; ++t) {
            workers.emplace_back([&] {
                while (true) {
                    const std::size_t begin = next.fetch_add(1);
                    if (begin >= inputs.size()) break;
                    result.likelihoods[begin] = fn(inputs[begin]);
                }
            });
        }
        for (auto& worker : workers) worker.join();
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    result.seconds = seconds;
    result.pairs_per_second = static_cast<double>(inputs.size()) * iterations / seconds;
    for (double x : result.likelihoods) result.checksum += x;
    return result;
}

}  // namespace

double log10_likelihood_scalar(const PairInput& input, bool tristate_correction) {
    configure_strict_fp_environment();
    return scalar_once(input, tristate_correction);
}

BatchResult compute_scalar(const std::vector<PairInput>& inputs, int threads, int iterations) {
    if (threads < 1) throw std::invalid_argument("threads must be positive");
    configure_strict_fp_environment();
    BatchResult result; result.likelihoods.resize(inputs.size());
    if (threads == 1) {
        ScalarWorkspace ws;
        for (std::size_t i = 0; i < inputs.size(); ++i)
            result.likelihoods[i] = scalar_once(inputs[i], true, ws);
        const auto start = std::chrono::steady_clock::now();
        for (int it = 0; it < iterations; ++it)
            for (std::size_t i = 0; i < inputs.size(); ++i)
                result.likelihoods[i] = scalar_once(inputs[i], true, ws);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        result.seconds = seconds;
        result.pairs_per_second = static_cast<double>(inputs.size()) * iterations / seconds;
        for (double x : result.likelihoods) result.checksum += x;
        return result;
    }
    const int actual = std::min<int>(threads, std::max<std::size_t>(1, inputs.size()));
    std::atomic<std::size_t> next{0};
    std::barrier warmup_barrier(actual + 1), start_barrier(actual + 1), done_barrier(actual + 1);
    std::vector<std::thread> workers;
    for (int t = 0; t < actual; ++t) {
        workers.emplace_back([&] {
            configure_strict_fp_environment();
            ScalarWorkspace ws;
            while (true) {
                const std::size_t index = next.fetch_add(1);
                if (index >= inputs.size()) break;
                result.likelihoods[index] = scalar_once(inputs[index], true, ws);
            }
            warmup_barrier.arrive_and_wait();
            for (int it = 0; it < iterations; ++it) {
                start_barrier.arrive_and_wait();
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= inputs.size()) break;
                    result.likelihoods[index] = scalar_once(inputs[index], true, ws);
                }
                done_barrier.arrive_and_wait();
            }
        });
    }
    warmup_barrier.arrive_and_wait();
    const auto start = std::chrono::steady_clock::now();
    for (int it = 0; it < iterations; ++it) {
        next.store(0);
        start_barrier.arrive_and_wait();
        done_barrier.arrive_and_wait();
    }
    for (auto& worker : workers) worker.join();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    result.seconds = seconds;
    result.pairs_per_second = static_cast<double>(inputs.size()) * iterations / seconds;
    for (double x : result.likelihoods) result.checksum += x;
    return result;
}

template <typename Workspace, typename Fn>
BatchResult compute_simd_impl(const std::vector<PairInput>& inputs, int threads, int iterations,
                              std::size_t lanes, Fn fn) {
    if (threads < 1) throw std::invalid_argument("threads must be positive");
    BatchResult result; result.likelihoods.resize(inputs.size());
    if (threads == 1) {
        Workspace ws;
        for (std::size_t group = 0; group < inputs.size(); group += lanes) {
            const std::size_t count = std::min<std::size_t>(lanes, inputs.size() - group);
            fn(inputs.data() + group, count, result.likelihoods.data() + group, true, ws);
        }
        const auto start = std::chrono::steady_clock::now();
        for (int it = 0; it < iterations; ++it)
            for (std::size_t group = 0; group < inputs.size(); group += lanes) {
                const std::size_t count = std::min<std::size_t>(lanes, inputs.size() - group);
                fn(inputs.data() + group, count, result.likelihoods.data() + group, true, ws);
            }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        result.seconds = seconds;
        result.pairs_per_second = static_cast<double>(inputs.size()) * iterations / seconds;
        for (double x : result.likelihoods) result.checksum += x;
        return result;
    }
    const int actual = std::min<int>(threads, std::max<std::size_t>(1, inputs.size()));
    std::atomic<std::size_t> next{0};
    std::barrier warmup_barrier(actual + 1), start_barrier(actual + 1), done_barrier(actual + 1);
    std::vector<std::thread> workers;
    for (int t = 0; t < actual; ++t) {
        workers.emplace_back([&] {
            Workspace ws;
            while (true) {
                const std::size_t group = next.fetch_add(lanes);
                if (group >= inputs.size()) break;
                const std::size_t count = std::min<std::size_t>(lanes, inputs.size() - group);
                fn(inputs.data() + group, count, result.likelihoods.data() + group, true, ws);
            }
            warmup_barrier.arrive_and_wait();
            for (int it = 0; it < iterations; ++it) {
                start_barrier.arrive_and_wait();
                while (true) {
                    const std::size_t group = next.fetch_add(lanes);
                    if (group >= inputs.size()) break;
                    const std::size_t count = std::min<std::size_t>(lanes, inputs.size() - group);
                    fn(inputs.data() + group, count, result.likelihoods.data() + group, true, ws);
                }
                done_barrier.arrive_and_wait();
            }
        });
    }
    warmup_barrier.arrive_and_wait();
    const auto start = std::chrono::steady_clock::now();
    for (int it = 0; it < iterations; ++it) {
        next.store(0);
        start_barrier.arrive_and_wait();
        done_barrier.arrive_and_wait();
    }
    for (auto& worker : workers) worker.join();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    result.seconds = seconds;
    result.pairs_per_second = static_cast<double>(inputs.size()) * iterations / seconds;
    for (double x : result.likelihoods) result.checksum += x;
    return result;
}

bool cpu_supports_avx2_impl() {
#if defined(__x86_64__) || defined(__i386__)
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
#else
    return false;
#endif
}

bool cpu_supports_avx512_impl() {
#if defined(__x86_64__) || defined(__i386__)
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f");
#else
    return false;
#endif
#else
    return false;
#endif
}

bool cpu_is_amd_family_19h() {
#if FASTGATK_X86_SIMD
    unsigned max_leaf = 0, ebx = 0, ecx = 0, edx = 0;
    max_leaf = __get_cpuid_max(0, nullptr);
    if (max_leaf < 1) return false;
    __cpuid(0, max_leaf, ebx, ecx, edx);
    char vendor[13]{};
    std::memcpy(vendor, &ebx, 4);
    std::memcpy(vendor + 4, &edx, 4);
    std::memcpy(vendor + 8, &ecx, 4);
    if (std::strcmp(vendor, "AuthenticAMD") != 0) return false;
    unsigned eax = 0;
    __cpuid(1, eax, ebx, ecx, edx);
    unsigned family = (eax >> 8U) & 0x0fU;
    if (family == 0x0fU) family += (eax >> 20U) & 0xffU;
    return family == 0x19U;
#else
    return false;
#endif
}

SimdBackend preferred_simd_backend_impl() {
    const bool avx512 = cpu_supports_avx512_impl();
    const bool avx2 = cpu_supports_avx2_impl();
    // Zen 4 executes 512-bit operations through narrower internal datapaths. This
    // kernel is ~39% faster with AVX2 on the measured family-19h host, so Auto
    // selects by architecture policy while explicit --backend=avx512 remains available.
    if (avx2 && avx512 && cpu_is_amd_family_19h()) return SimdBackend::Avx2;
    if (avx512) return SimdBackend::Avx512;
    if (avx2) return SimdBackend::Avx2;
    return SimdBackend::Auto;
}

BatchResult compute_simd(const std::vector<PairInput>& inputs, int threads, int iterations,
                         SimdBackend backend) {
#if FASTGATK_X86_SIMD
    const bool avx512 = cpu_supports_avx512_impl();
    const bool avx2 = cpu_supports_avx2_impl();
    if (backend == SimdBackend::Avx512 && !avx512)
        throw std::runtime_error("AVX-512 backend requested but avx512f is unavailable");
    if (backend == SimdBackend::Avx2 && !avx2)
        throw std::runtime_error("AVX2 backend requested but avx2 is unavailable");
    const SimdBackend selected = backend == SimdBackend::Auto ? preferred_simd_backend_impl() : backend;
    if (selected == SimdBackend::Avx512) {
        return compute_simd_impl<Simd512Workspace>(inputs, threads, iterations, kAvx512Lanes,
                                                    simd_batch_once_avx512);
    }
    if (selected == SimdBackend::Avx2) {
        return compute_simd_impl<SimdWorkspace>(inputs, threads, iterations, kAvx2Lanes,
                                                simd_batch_once_avx2);
    }
    // Auto mode remains usable on non-x86/non-SIMD hosts; strict scalar is the
    // compatibility fallback rather than an illegal instruction.
    return compute_scalar(inputs, threads, iterations);
#else
    if (backend != SimdBackend::Auto)
        throw std::runtime_error("requested x86 SIMD backend on a non-x86 build");
    return compute_scalar(inputs, threads, iterations);
#endif
}

std::string cpu_simd_description() {
    const bool avx512 = cpu_supports_avx512_impl();
    const bool avx2 = cpu_supports_avx2_impl();
    const SimdBackend preferred = preferred_simd_backend_impl();
    if (avx512 && preferred == SimdBackend::Avx2)
        return "AVX2 auto-selected; AVX-512 available (AMD family-19h policy)";
    if (avx512) return "AVX-512 auto-selected, 8x f64 lanes (runtime)";
    if (avx2) return "AVX2 auto-selected, 4x f64 lanes (runtime)";
    return "scalar fallback (runtime)";
}

bool cpu_supports_avx2() { return cpu_supports_avx2_impl(); }
bool cpu_supports_avx512() { return cpu_supports_avx512_impl(); }

}  // namespace fastgatk::pairhmm
