#include "fastgatk/kernels/pairhmm_kokkos.hpp"
#include "fastgatk/kernels/gpu_safety.hpp"
#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_BitManipulation.hpp>
#include <Kokkos_SIMD.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <tuple>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <string>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

#if defined(__linux__)
#include <system_error>
#endif

namespace fastgatk::pairhmm {
namespace {

// `fragment_aggregation_*`, `allele_marginalization_*` and the
// `likelihood_*` reduce kernels allocate transient device Views plus
// their host mirrors.  Kokkos OpenMP's default memory pool keeps slab
// growth after the View's C++ scope exits, which inflates cumulative
// peak RSS across the dozens of region/fragment/per-read reductions
// invoked by a single Mutect2 invocation.  This helper explicitly
// releases each View back to the slab before the function returns;
// subsequent callers that need a smaller View see the slab reused
// instead of accumulated.
//
// `Kokkos::resize(view, 0)` shrinks the live extent and returns the
// underlying allocation to the pool on the OpenMP backend (the only
// backend currently used).  On other backends the call is a no-op.
//
// `Kokkos::view_alloc(WithoutInitializing)` removes the implicit
// `memset(0)` cost for Views that are fully overwritten by the host
// fill loop or by `Kokkos::deep_copy` before any read; this is purely
// a wall-clock reduction (no semantic effect).
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

template <typename View>
inline void shrink_release(View& view) noexcept {
    if (view.span() != 0) Kokkos::resize(view, 0);
}

template <typename... Views>
inline void release_views(Views&... views) noexcept {
    (shrink_release(views), ...);
}

using Simd = Kokkos::Experimental::simd<double>;
using View2D = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;
using IndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
using ResultView = Kokkos::View<double*, MemorySpace>;
using HostView2D = Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>;
using HostIndexView = Kokkos::View<std::uint32_t*, Kokkos::HostSpace>;
using WorkspaceView = Kokkos::View<Simd**, Kokkos::LayoutRight, MemorySpace>;

constexpr double kTristateCorrection = 3.0;

// Java 17 Math.log10 delegates to StrictMath.log10, whose native reference is
// OpenJDK fdlibm.  System libm implementations are allowed to choose the other
// adjacent double, which showed up in 3/4096 matrix results on this host.  This
// device-callable formulation preserves the fdlibm operation order and exact
// constants so compatibility does not depend on the host C library.
// Algorithm/provenance: OpenJDK 17u libfdlibm e_log.c and e_log10.c.
KOKKOS_INLINE_FUNCTION double gatk_strict_log(double x) {
    constexpr double ln2_hi = 6.93147180369123816490e-01;
    constexpr double ln2_lo = 1.90821492927058770002e-10;
    constexpr double two54 = 1.80143985094819840000e+16;
    constexpr double lg1 = 6.666666666666735130e-01;
    constexpr double lg2 = 3.999999999940941908e-01;
    constexpr double lg3 = 2.857142874366239149e-01;
    constexpr double lg4 = 2.222219843214978396e-01;
    constexpr double lg5 = 1.818357216161805012e-01;
    constexpr double lg6 = 1.531383769920937332e-01;
    constexpr double lg7 = 1.479819860511658591e-01;
    std::uint64_t bits = Kokkos::bit_cast<std::uint64_t>(x);
    std::int32_t hx = static_cast<std::int32_t>(bits >> 32);
    const std::uint32_t lx = static_cast<std::uint32_t>(bits);
    int k = 0;
    if (hx < 0x00100000) {
        if (((hx & 0x7fffffff) | lx) == 0) return -two54 / 0.0;
        if (hx < 0) return (x - x) / 0.0;
        k -= 54;
        x *= two54;
        bits = Kokkos::bit_cast<std::uint64_t>(x);
        hx = static_cast<std::int32_t>(bits >> 32);
    }
    if (hx >= 0x7ff00000) return x + x;
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    const int i_normalize = (hx + 0x95f64) & 0x100000;
    const std::uint32_t normalized_hi = static_cast<std::uint32_t>(
        hx | (i_normalize ^ 0x3ff00000));
    bits = (bits & 0x00000000ffffffffULL) |
           (static_cast<std::uint64_t>(normalized_hi) << 32);
    x = Kokkos::bit_cast<double>(bits);
    k += i_normalize >> 20;
    const double f = x - 1.0;
    if ((0x000fffff & (2 + hx)) < 3) {
        if (f == 0.0) {
            if (k == 0) return 0.0;
            const double dk = static_cast<double>(k);
            return dk * ln2_hi + dk * ln2_lo;
        }
        const double r = f * f * (0.5 - 0.33333333333333333 * f);
        if (k == 0) return f - r;
        const double dk = static_cast<double>(k);
        return dk * ln2_hi - ((r - dk * ln2_lo) - f);
    }
    const double s = f / (2.0 + f);
    const double dk = static_cast<double>(k);
    const double z = s * s;
    const int i_range = hx - 0x6147a;
    const double w = z * z;
    const int j_range = 0x6b851 - hx;
    const double t1 = w * (lg2 + w * (lg4 + w * lg6));
    const double t2 = z * (lg1 + w * (lg3 + w * (lg5 + w * lg7)));
    const double r = t2 + t1;
    if ((i_range | j_range) > 0) {
        const double hfsq = 0.5 * f * f;
        if (k == 0) return f - (hfsq - s * (hfsq + r));
        return dk * ln2_hi - ((hfsq - (s * (hfsq + r) + dk * ln2_lo)) - f);
    }
    if (k == 0) return f - s * (f - r);
    return dk * ln2_hi - ((s * (f - r) - dk * ln2_lo) - f);
}

KOKKOS_INLINE_FUNCTION double gatk_strict_log10(double x) {
    constexpr double two54 = 1.80143985094819840000e+16;
    constexpr double ivln10 = 4.34294481903251816668e-01;
    constexpr double log10_2hi = 3.01029995663611771306e-01;
    constexpr double log10_2lo = 3.69423907715893078616e-13;
    std::uint64_t bits = Kokkos::bit_cast<std::uint64_t>(x);
    std::int32_t hx = static_cast<std::int32_t>(bits >> 32);
    const std::uint32_t lx = static_cast<std::uint32_t>(bits);
    int k = 0;
    if (hx < 0x00100000) {
        if (((hx & 0x7fffffff) | lx) == 0) return -two54 / 0.0;
        if (hx < 0) return (x - x) / 0.0;
        k -= 54;
        x *= two54;
        bits = Kokkos::bit_cast<std::uint64_t>(x);
        hx = static_cast<std::int32_t>(bits >> 32);
    }
    if (hx >= 0x7ff00000) return x + x;
    k += (hx >> 20) - 1023;
    const int i = static_cast<int>((static_cast<std::uint32_t>(k) & 0x80000000U) >> 31);
    const std::uint32_t normalized_hi = static_cast<std::uint32_t>(
        (hx & 0x000fffff) | ((0x3ff - i) << 20));
    bits = (bits & 0x00000000ffffffffULL) |
           (static_cast<std::uint64_t>(normalized_hi) << 32);
    x = Kokkos::bit_cast<double>(bits);
    const double y = static_cast<double>(k + i);
    const double z = y * log10_2lo + ivln10 * gatk_strict_log(x);
    return z + y * log10_2hi;
}

struct GatkTables {
    std::array<double, 256> quality{};
    std::array<double, 256 * 256> match{};
};

// GATK's default native PairHMM implementation is Intel GKL, rather than the
// Java LoglessPairHMM.  Its Float32 Context builds a *different* transition
// cache from the Java Double cache: qualities come from powf(), the Jacobian
// correction table is stored as float, and the historical GKL constants are
// deliberately shortened (0.434294 instead of 1 / log(10)).  Feeding a
// Float32 recurrence with Java's double transition table and casting it late
// gives a numerically plausible answer, but it is not the answer that Mutect2
// gets through VectorLoglessPairHMM by default.  Keep this source-compatible
// cache separate from the pinned Java table above: Float64 remains the strict
// Java-compatible mode, while Float32 models GKL 0.8.11.
struct GklFloatTables {
    std::array<float, 128> quality{};
    std::array<float, ((255 * 256) >> 1)> match{};
};

float gkl_float_approximate_log10_sum(float first, float second,
                                      const std::array<float, 80001>& jacobian) {
    if (first > second) std::swap(first, second);
    if (std::isinf(first) || std::isinf(second)) return second;
    const float difference = second - first;
    if (difference >= 8.0F) return second;
    // Context<float>::fastRound: positive values are rounded by adding 0.5f
    // before the integer cast.  Its input is necessarily non-negative here.
    const int index = static_cast<int>(difference * 10000.0F + 0.5F);
    return second + jacobian[static_cast<std::size_t>(index)];
}

const GklFloatTables& gkl_float_tables() {
    static const GklFloatTables tables = [] {
        GklFloatTables values{};
        std::array<float, 80001> jacobian{};
        for (std::size_t index = 0; index < jacobian.size(); ++index) {
            // ContextBase<float>::initializeJacobianLogTable() computes this
            // expression in double and stores the rounded Float32 value.
            jacobian[index] = static_cast<float>(std::log10(
                1.0 + std::pow(10.0, -static_cast<double>(index) * 0.0001)));
        }
        for (std::size_t quality = 0; quality < values.quality.size(); ++quality) {
            // Match Context<float>::ph2pr initialization, including powf().
            values.quality[quality] = std::pow(
                10.0F, -static_cast<float>(quality) / 10.0F);
        }
        constexpr double kGklInverseLn10 = 0.434294;
        for (std::size_t maximum = 0; maximum <= 254; ++maximum) {
            for (std::size_t minimum = 0; minimum <= maximum; ++minimum) {
                const float log10_sum = gkl_float_approximate_log10_sum(
                    -0.1F * static_cast<float>(minimum),
                    -0.1F * static_cast<float>(maximum), jacobian);
                const double match_to_match_log10 = std::log1p(-std::min(
                    1.0, std::pow(10.0, static_cast<double>(log10_sum)))) *
                    kGklInverseLn10;
                values.match[(maximum * (maximum + 1) >> 1) + minimum] =
                    static_cast<float>(std::pow(10.0, match_to_match_log10));
            }
        }
        return values;
    }();
    return tables;
}

float gkl_float_match_to_match(const GklFloatTables& tables,
                               const std::uint8_t insertion_quality,
                               const std::uint8_t deletion_quality) {
    // GKL's testcase bridge masks every native quality byte to the same
    // seven-bit domain before Context<float> looks it up.
    std::size_t insertion = static_cast<std::size_t>(insertion_quality & 127U);
    std::size_t deletion = static_cast<std::size_t>(deletion_quality & 127U);
    const std::size_t minimum = std::min(insertion, deletion);
    const std::size_t maximum = std::max(insertion, deletion);
    return tables.match[(maximum * (maximum + 1) >> 1) + minimum];
}

// Keep the Jacobian table construction in the same operation order as
// MathUtils.JacobianLogTable in GATK.  In particular, use the strict log10
// implementation for the outer logarithm instead of the host libm's log10;
// this prevents a platform-dependent adjacent-double choice from leaking into
// every PairHMM transition table entry.
double gatk_jacobian_correction(const int index) {
    constexpr double table_step = 0.0001;
    return gatk_strict_log10(
        1.0 + std::pow(10.0, -static_cast<double>(index) * table_step));
}

double approximate_log10_sum(double a, double b) {
    if (a > b) std::swap(a, b);
    if (std::isinf(a) && a < 0.0) return b;
    const double diff = b - a;
    if (diff >= 8.0) return b;
    const int index = static_cast<int>(diff * 10000.0 + 0.5);
    return b + gatk_jacobian_correction(index);
}

bool load_gatk_tables_file(const char* path, GatkTables& values, const bool required) {
    if (path == nullptr || *path == '\0') return false;
    std::ifstream input(path);
    if (!input) {
        if (required)
            throw std::runtime_error(std::string("cannot open FAST_GATK_PAIRHMM_TABLES file: ") + path);
        return false;
    }
    std::string kind;
    std::uint64_t bits = 0;
    bool ok = true;
    for (std::size_t q = 0; ok && q < 256; ++q) {
        input >> kind >> std::hex >> bits >> std::dec;
        ok = input && kind == "q";
        if (ok) std::memcpy(&values.quality[q], &bits, sizeof(bits));
    }
    for (std::size_t ins = 0; ok && ins < 256; ++ins) {
        for (std::size_t del = 0; ok && del < 256; ++del) {
            std::size_t file_ins = 0, file_del = 0;
            input >> kind >> file_ins >> file_del >> std::hex >> bits >> std::dec;
            ok = input && kind == "m" && file_ins == ins && file_del == del;
            if (ok) std::memcpy(&values.match[ins * 256 + del], &bits, sizeof(bits));
        }
    }
    if (!ok)
        throw std::runtime_error(std::string("invalid FAST_GATK_PAIRHMM_TABLES file: ") + path);
    return true;
}

#if defined(FASTGATK_INSTALLED_PAIRHMM_TABLES_RELATIVE_PATH)
bool load_relative_installed_gatk_tables(GatkTables& values) {
#if defined(__linux__)
    std::error_code error;
    const auto executable = std::filesystem::read_symlink("/proc/self/exe", error);
    if (!error && !executable.empty()) {
        const auto candidate = executable.parent_path() /
            FASTGATK_INSTALLED_PAIRHMM_TABLES_RELATIVE_PATH;
        if (load_gatk_tables_file(candidate.c_str(), values, false)) return true;
    }
#else
    (void)values;
#endif
    return false;
}
#endif

#if defined(FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH) || \
    defined(FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH)
bool load_compiled_default_gatk_tables(GatkTables& values) {
#ifdef FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH
    if (load_gatk_tables_file(FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH, values, false))
        return true;
#endif
#ifdef FASTGATK_INSTALLED_PAIRHMM_TABLES_RELATIVE_PATH
    if (load_relative_installed_gatk_tables(values)) return true;
#endif
#ifdef FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH
    // A package build may not retain the source tree.  Probe the configured
    // install location after the auditable source-tree path, while preserving
    // the old math fallback for intentionally minimal kernel-only builds.
#if defined(FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH)
    if (std::string(FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH) !=
        std::string(FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH) &&
        load_gatk_tables_file(FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH, values, false))
        return true;
#else
    if (load_gatk_tables_file(FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH, values, false))
        return true;
#endif
#endif
    return false;
}
#endif

const GatkTables& gatk_tables() {
    static const GatkTables tables = [] {
        GatkTables values{};
        for (std::size_t q = 0; q < 256; ++q)
            values.quality[q] = std::pow(10.0, -0.1 * static_cast<double>(q));
        for (std::size_t ins = 0; ins < 256; ++ins) {
            for (std::size_t del = 0; del < 256; ++del) {
                // PairHMMModel initializes matchToMatchProb by first forming
                // matchToMatchLog10 with log1p(-min(1, 10^sum)) / ln(10), and
                // then converting that log10 value back to linear scale.  The
                // seemingly equivalent `1 - 10^sum` is not equivalent at the
                // low-quality boundary (for Q=0 it is even negative).
                const double log10_sum = approximate_log10_sum(
                    -0.1 * static_cast<double>(ins),
                    -0.1 * static_cast<double>(del));
                const double event_probability = std::min(
                    1.0, std::pow(10.0, log10_sum));
                if (event_probability >= 1.0) {
                    values.match[ins * 256 + del] = 0.0;
                } else {
                    const double no_event_log10 =
                        std::log1p(-event_probability) /
                        2.3025850929940456840179914546843642;
                    values.match[ins * 256 + del] =
                        std::pow(10.0, no_event_log10);
                }
            }
        }
        if (const char* path = std::getenv("FAST_GATK_PAIRHMM_TABLES")) {
            // An explicit override is fail-closed: a typo must not silently
            // switch a strict deployment back to platform-dependent libm.
            load_gatk_tables_file(path, values, true);
        } else {
#ifdef FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH
            // Source-tree/package builds may provide the pinned Java table as
            // the default.  If the data file is not staged (for example a
            // deliberately minimal kernel-only install), retain the portable
            // math fallback rather than making library initialization fail.
            load_compiled_default_gatk_tables(values);
#elif defined(FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH)
            load_compiled_default_gatk_tables(values);
#endif
        }
        return values;
    }();
    return tables;
}

void validate(const std::vector<PairInput>& records, const PairIndexBatch& pairs,
              std::size_t& read_length, std::size_t& haplotype_length) {
    if (records.empty() || pairs.read_ids.empty() ||
        pairs.read_ids.size() != pairs.haplotype_ids.size())
        throw std::invalid_argument("empty or inconsistent Kokkos PairHMM batch");
    read_length = records.front().read.size();
    haplotype_length = records.front().haplotype.size();
    if (read_length == 0 || haplotype_length == 0)
        throw std::invalid_argument("empty read or haplotype");
    for (const PairInput& record : records) {
        if (record.read.size() != read_length || record.haplotype.size() != haplotype_length ||
            record.read_qual.size() != read_length || record.insertion_gop.size() != read_length ||
            record.deletion_gop.size() != read_length || record.gap_continuation.size() != read_length)
            throw std::invalid_argument(
                "rectangular Kokkos PairHMM kernel requires uniform batch lengths; "
                "use compute_kokkos_bucketed or PersistentBucketPlan for ragged input");
    }
    for (std::size_t i = 0; i < pairs.read_ids.size(); ++i) {
        if (pairs.read_ids[i] >= records.size() || pairs.haplotype_ids[i] >= records.size())
            throw std::invalid_argument("pair index exceeds record count");
    }
}

template <class View>
auto mirror(const View& view) {
    // Host mirrors are filled in the loop below before the deep_copy,
    // so the implicit memset(0) inside `create_mirror_view` is pure
    // overhead.  Skip it via `WithoutInitializing`.  GPU-safe: a
    // Kokkos View construction runtime property; works on every backend.
    return Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), view);
}

}  // namespace

LikelihoodNormalizationResult normalize_likelihoods_kokkos(
    const std::vector<double>& likelihoods,
    const std::vector<std::uint32_t>& read_ids,
    const std::vector<std::uint8_t>& eligible_for_best,
    const std::size_t read_count,
    const double maximum_likelihood_difference_cap) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (likelihoods.size() != read_ids.size() ||
        likelihoods.size() != eligible_for_best.size())
        throw std::invalid_argument("likelihood normalization dimensions do not match");
    const bool cap_disabled = std::isinf(maximum_likelihood_difference_cap) &&
        maximum_likelihood_difference_cap < 0.0;
    if (std::isnan(maximum_likelihood_difference_cap) ||
        (maximum_likelihood_difference_cap < 0.0 && !cap_disabled) ||
        (std::isinf(maximum_likelihood_difference_cap) && !cap_disabled))
        throw std::invalid_argument(
            "likelihood normalization cap must be finite and non-negative, or negative infinity to disable");
    for (const auto read : read_ids)
        if (read >= read_count)
            throw std::invalid_argument("likelihood normalization read id exceeds read count");

    LikelihoodNormalizationResult result;
    result.likelihoods = likelihoods;
    result.best_by_read.assign(read_count, -std::numeric_limits<double>::infinity());
    result.execution_space = ExecSpace::name();
    using OffsetView = Kokkos::View<std::size_t*, MemorySpace>;
    using SortedIndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using ReadIdView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using EligibleView = Kokkos::View<std::uint8_t*, MemorySpace>;
    using ValueView = Kokkos::View<double*, MemorySpace>;

    const auto request_count = likelihoods.size();
    std::vector<std::uint32_t> sorted(request_count);
    std::iota(sorted.begin(), sorted.end(), 0U);
    std::stable_sort(sorted.begin(), sorted.end(), [&](const auto left, const auto right) {
        if (read_ids[left] != read_ids[right]) return read_ids[left] < read_ids[right];
        return left < right;
    });
    std::vector<std::size_t> offsets(read_count + 1, 0);
    for (const auto index : sorted) ++offsets[read_ids[index] + 1];
    for (std::size_t read = 0; read < read_count; ++read)
        offsets[read + 1] += offsets[read];

    fastgatk::core::HostBatch host_batch("pairhmm-likelihood-normalization-v1");
    host_batch.records = request_count;
    host_batch.bytes = likelihoods.size() * sizeof(double) +
        read_ids.size() * sizeof(std::uint32_t) + eligible_for_best.size() * sizeof(std::uint8_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm_likelihood_normalization");
    plan.begin_prepare(host_batch);
    // All six device Views are fully overwritten before any kernel read.
    OffsetView device_offsets(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_read_offsets"),
        offsets.size());
    SortedIndexView device_sorted(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_sorted_indices"),
        sorted.size());
    ReadIdView device_read_ids(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_read_ids"),
        read_ids.size());
    EligibleView device_eligible(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_eligible"),
        eligible_for_best.size());
    ValueView device_values(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_values"),
        likelihoods.size());
    ValueView device_best(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "likelihood_best"),
        read_count);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    auto host_sorted = Kokkos::create_mirror_view(device_sorted);
    auto host_read_ids = Kokkos::create_mirror_view(device_read_ids);
    auto host_eligible = Kokkos::create_mirror_view(device_eligible);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t index = 0; index < offsets.size(); ++index) host_offsets(index) = offsets[index];
    for (std::size_t index = 0; index < sorted.size(); ++index) host_sorted(index) = sorted[index];
    for (std::size_t index = 0; index < read_ids.size(); ++index) {
        host_read_ids(index) = read_ids[index];
        host_eligible(index) = eligible_for_best[index];
        host_values(index) = likelihoods[index];
    }
    Kokkos::deep_copy(device_offsets, host_offsets);
    Kokkos::deep_copy(device_sorted, host_sorted);
    Kokkos::deep_copy(device_read_ids, host_read_ids);
    Kokkos::deep_copy(device_eligible, host_eligible);
    Kokkos::deep_copy(device_values, host_values);
    Kokkos::deep_copy(device_best, -std::numeric_limits<double>::infinity());
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(request_count);
    device_batch.bind("read_offsets", device_offsets);
    device_batch.bind("sorted_indices", device_sorted);
    device_batch.bind("read_ids", device_read_ids);
    device_batch.bind("eligible", device_eligible);
    device_batch.bind("values", device_values);
    device_batch.bind("best", device_best);
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const auto execute_start = std::chrono::steady_clock::now();
    Kokkos::parallel_for("fastgatk_likelihood_best", RangePolicy(0, read_count),
        KOKKOS_LAMBDA(const std::size_t read) {
            double best = -std::numeric_limits<double>::infinity();
            for (std::size_t cursor = device_offsets(read);
                 cursor < device_offsets(read + 1); ++cursor) {
                const auto request = device_sorted(cursor);
                if (device_eligible(request) != 0 && device_values(request) > best)
                    best = device_values(request);
            }
            device_best(read) = best;
        });
    // GATK uses Double.NEGATIVE_INFINITY as an explicit no-cap mode and
    // returns the PairHMM values unchanged. Do not replace it with a very
    // large finite width: that converts valid -Infinity zero-probability
    // cells into finite likelihoods and changes subsequent allele
    // marginalization. The best reduction above is retained for telemetry.
    if (!cap_disabled) {
        Kokkos::parallel_for("fastgatk_likelihood_cap", RangePolicy(0, request_count),
            KOKKOS_LAMBDA(const std::size_t request) {
                const auto read = device_read_ids(request);
                const auto best = device_best(read);
                if (best == -std::numeric_limits<double>::infinity()) return;
                const auto floor = best - maximum_likelihood_difference_cap;
                if (device_values(request) < floor) device_values(request) = floor;
            });
    }
    ExecSpace().fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(execute_end - execute_start).count();
    auto host_best = Kokkos::create_mirror_view(device_best);
    auto host_result = Kokkos::create_mirror_view(device_values);
    Kokkos::deep_copy(host_best, device_best);
    Kokkos::deep_copy(host_result, device_values);
    for (std::size_t read = 0; read < read_count; ++read) result.best_by_read[read] = host_best(read);
    for (std::size_t request = 0; request < request_count; ++request) result.likelihoods[request] = host_result(request);
    // Return device allocations to the OpenMP memory pool before this
    // function's scope exits.  See fragment_aggregation comment for the
    // slab-growth rationale.
    release_views(device_offsets, device_sorted, device_read_ids,
                  device_eligible, device_values, device_best);
    return result;
}

AlleleMarginalizationResult marginalize_read_allele_likelihoods_kokkos(
    const std::vector<AlleleMarginalizationRequest>& requests,
    const std::size_t row_count,
    const std::size_t allele_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (allele_count == 0)
        throw std::invalid_argument("allele marginalization requires at least one allele");
    if (row_count > std::numeric_limits<std::size_t>::max() / allele_count)
        throw std::invalid_argument("allele marginalization output dimensions overflow");
    for (const auto& request : requests) {
        if (request.row_id >= row_count)
            throw std::invalid_argument("allele marginalization row id exceeds row count");
        if (static_cast<std::size_t>(request.allele) >= allele_count)
            throw std::invalid_argument("allele marginalization allele index exceeds allele count");
    }

    AlleleMarginalizationResult result;
    result.allele_count = allele_count;
    result.best_by_row_allele.assign(
        row_count * allele_count, -std::numeric_limits<double>::infinity());
    result.execution_space = ExecSpace::name();
    using OffsetView = Kokkos::View<std::size_t*, MemorySpace>;
    using SortedIndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using ValueView = Kokkos::View<double*, MemorySpace>;

    const auto request_count = requests.size();
    std::vector<std::uint32_t> sorted(request_count);
    std::iota(sorted.begin(), sorted.end(), 0U);
    std::stable_sort(sorted.begin(), sorted.end(), [&](const auto left, const auto right) {
        if (requests[left].row_id != requests[right].row_id)
            return requests[left].row_id < requests[right].row_id;
        if (requests[left].allele != requests[right].allele)
            return requests[left].allele < requests[right].allele;
        return left < right;
    });
    const auto cell_count = row_count * allele_count;
    std::vector<std::size_t> offsets(cell_count + 1, 0);
    for (const auto index : sorted)
        ++offsets[static_cast<std::size_t>(requests[index].row_id) * allele_count +
                  requests[index].allele + 1];
    for (std::size_t cell = 0; cell < cell_count; ++cell)
        offsets[cell + 1] += offsets[cell];

    fastgatk::core::HostBatch host_batch("pairhmm-allele-marginalization-v1");
    host_batch.records = request_count;
    host_batch.bytes = request_count * (sizeof(double) + sizeof(std::uint32_t) +
                                        sizeof(std::uint8_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm_allele_marginalization");
    plan.begin_prepare(host_batch);
    // All four device Views are fully overwritten by `Kokkos::deep_copy`
    // before any kernel read; skip the implicit `memset(0)` cost.
    OffsetView device_offsets(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "allele_marginalization_offsets"),
        offsets.size());
    SortedIndexView device_sorted(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "allele_marginalization_sorted"),
        sorted.size());
    ValueView device_values(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "allele_marginalization_values"),
        request_count);
    ValueView device_best(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "allele_marginalization_best"),
        cell_count);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    auto host_sorted = Kokkos::create_mirror_view(device_sorted);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t cell = 0; cell < offsets.size(); ++cell)
        host_offsets(cell) = offsets[cell];
    for (std::size_t index = 0; index < sorted.size(); ++index)
        host_sorted(index) = sorted[index];
    for (std::size_t index = 0; index < requests.size(); ++index)
        host_values(index) = requests[index].likelihood;
    Kokkos::deep_copy(device_offsets, host_offsets);
    Kokkos::deep_copy(device_sorted, host_sorted);
    Kokkos::deep_copy(device_values, host_values);
    Kokkos::deep_copy(device_best, -std::numeric_limits<double>::infinity());
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(request_count);
    device_batch.bind("offsets", device_offsets);
    device_batch.bind("sorted", device_sorted);
    device_batch.bind("values", device_values);
    device_batch.bind("best", device_best);
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const auto execute_start = std::chrono::steady_clock::now();
    Kokkos::parallel_for("fastgatk_allele_marginalization", RangePolicy(0, cell_count),
        KOKKOS_LAMBDA(const std::size_t cell) {
            double best = -std::numeric_limits<double>::infinity();
            for (std::size_t cursor = device_offsets(cell);
                 cursor < device_offsets(cell + 1); ++cursor) {
                const auto request = device_sorted(cursor);
                if (device_values(request) > best) best = device_values(request);
            }
            device_best(cell) = best;
        });
    ExecSpace().fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(execute_end - execute_start).count();
    auto host_best = Kokkos::create_mirror_view(device_best);
    Kokkos::deep_copy(host_best, device_best);
    for (std::size_t cell = 0; cell < cell_count; ++cell)
        result.best_by_row_allele[cell] = host_best(cell);
    // Return device allocations to the OpenMP memory pool before this
    // function's scope exits.  See fragment_aggregation comment for the
    // slab-growth rationale.
    release_views(device_offsets, device_sorted, device_values, device_best);
    return result;
}

AlleleMarginalizationResult marginalize_read_allele_likelihoods_kokkos(
    const std::vector<AlleleMarginalizationRequest>& requests,
    const std::size_t row_count) {
    return marginalize_read_allele_likelihoods_kokkos(requests, row_count, 2);
}

FragmentHaplotypeAggregationResult aggregate_fragment_haplotype_likelihoods_kokkos(
    const std::vector<FragmentHaplotypeAggregationRequest>& requests,
    const std::size_t cell_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    for (const auto& request : requests)
        if (request.cell_id >= cell_count)
            throw std::invalid_argument("fragment aggregation cell id exceeds cell count");

    FragmentHaplotypeAggregationResult result;
    result.sums_by_cell.assign(cell_count, 0.0);
    result.execution_space = ExecSpace::name();
    using OffsetView = Kokkos::View<std::size_t*, MemorySpace>;
    using SortedIndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using ValueView = Kokkos::View<double*, MemorySpace>;

    const auto request_count = requests.size();
    std::vector<std::uint32_t> sorted(request_count);
    std::iota(sorted.begin(), sorted.end(), 0U);
    std::stable_sort(sorted.begin(), sorted.end(), [&](const auto left, const auto right) {
        if (requests[left].cell_id != requests[right].cell_id)
            return requests[left].cell_id < requests[right].cell_id;
        return left < right;
    });
    std::vector<std::size_t> offsets(cell_count + 1, 0);
    for (const auto index : sorted) ++offsets[requests[index].cell_id + 1];
    for (std::size_t cell = 0; cell < cell_count; ++cell)
        offsets[cell + 1] += offsets[cell];

    fastgatk::core::HostBatch host_batch("pairhmm-fragment-aggregation-v1");
    host_batch.records = request_count;
    host_batch.bytes = request_count * (sizeof(double) + sizeof(std::uint32_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm_fragment_aggregation");
    plan.begin_prepare(host_batch);
    // All four device Views are fully overwritten by the host fill loop or
    // by `Kokkos::deep_copy` before any kernel read.  Skip the implicit
    // `memset(0)` cost via `WithoutInitializing`.
    OffsetView device_offsets(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "fragment_aggregation_offsets"),
        offsets.size());
    SortedIndexView device_sorted(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "fragment_aggregation_sorted"),
        sorted.size());
    ValueView device_values(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "fragment_aggregation_values"),
        request_count);
    ValueView device_sums(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "fragment_aggregation_sums"),
        cell_count);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    auto host_sorted = Kokkos::create_mirror_view(device_sorted);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t cell = 0; cell < offsets.size(); ++cell)
        host_offsets(cell) = offsets[cell];
    for (std::size_t index = 0; index < sorted.size(); ++index)
        host_sorted(index) = sorted[index];
    for (std::size_t index = 0; index < requests.size(); ++index)
        host_values(index) = requests[index].likelihood;
    Kokkos::deep_copy(device_offsets, host_offsets);
    Kokkos::deep_copy(device_sorted, host_sorted);
    Kokkos::deep_copy(device_values, host_values);
    Kokkos::deep_copy(device_sums, 0.0);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(request_count);
    device_batch.bind("offsets", device_offsets);
    device_batch.bind("sorted", device_sorted);
    device_batch.bind("values", device_values);
    device_batch.bind("sums", device_sums);
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const auto execute_start = std::chrono::steady_clock::now();
    Kokkos::parallel_for("fastgatk_fragment_aggregation", RangePolicy(0, cell_count),
        KOKKOS_LAMBDA(const std::size_t cell) {
            double sum = 0.0;
            // This loop is intentionally scalar and ordered. Reassociation
            // would change Java double addition and can also hide -Infinity.
            for (std::size_t cursor = device_offsets(cell);
                 cursor < device_offsets(cell + 1); ++cursor)
                sum += device_values(device_sorted(cursor));
            device_sums(cell) = sum;
        });
    ExecSpace().fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(execute_end - execute_start).count();
    auto host_sums = Kokkos::create_mirror_view(device_sums);
    Kokkos::deep_copy(host_sums, device_sums);
    for (std::size_t cell = 0; cell < cell_count; ++cell)
        result.sums_by_cell[cell] = host_sums(cell);
    // Return device allocations to the OpenMP memory pool before this
    // function's scope exits.  Without these explicit shrinks the
    // underlying slabs are retained until program exit, inflating
    // cumulative peak RSS across many region/fragment reductions.
    release_views(device_offsets, device_sorted, device_values, device_sums);
    return result;
}

ReadAlleleUncertaintyResult reduce_read_allele_uncertainty_kokkos(
    const std::vector<ReadAlleleUncertaintyRequest>& requests,
    const std::size_t row_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (row_count > std::numeric_limits<std::size_t>::max() / 2)
        throw std::invalid_argument("read uncertainty row count overflows output dimensions");
    for (const auto& request : requests)
        if (request.row_id >= row_count)
            throw std::invalid_argument("read uncertainty row id exceeds row count");

    ReadAlleleUncertaintyResult result;
    result.best_second_by_row.assign(
        row_count * 2, -std::numeric_limits<double>::infinity());
    result.execution_space = ExecSpace::name();
    using OffsetView = Kokkos::View<std::size_t*, MemorySpace>;
    using SortedIndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using ValueView = Kokkos::View<double*, MemorySpace>;
    using BestSecondView = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

    const auto request_count = requests.size();
    std::vector<std::uint32_t> sorted(request_count);
    std::iota(sorted.begin(), sorted.end(), 0U);
    std::stable_sort(sorted.begin(), sorted.end(), [&](const auto left, const auto right) {
        if (requests[left].row_id != requests[right].row_id)
            return requests[left].row_id < requests[right].row_id;
        return left < right;
    });
    std::vector<std::size_t> offsets(row_count + 1, 0);
    for (const auto index : sorted)
        ++offsets[requests[index].row_id + 1];
    for (std::size_t row = 0; row < row_count; ++row)
        offsets[row + 1] += offsets[row];

    fastgatk::core::HostBatch host_batch("pairhmm-read-allele-uncertainty-v1");
    host_batch.records = request_count;
    host_batch.bytes = request_count * (sizeof(double) + sizeof(std::uint32_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm_read_allele_uncertainty");
    plan.begin_prepare(host_batch);
    OffsetView device_offsets("read_uncertainty_offsets", offsets.size());
    SortedIndexView device_sorted("read_uncertainty_sorted", sorted.size());
    ValueView device_values("read_uncertainty_values", request_count);
    BestSecondView device_best_second("read_uncertainty_best_second", row_count, 2);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    auto host_sorted = Kokkos::create_mirror_view(device_sorted);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t row = 0; row < offsets.size(); ++row)
        host_offsets(row) = offsets[row];
    for (std::size_t index = 0; index < sorted.size(); ++index)
        host_sorted(index) = sorted[index];
    for (std::size_t index = 0; index < requests.size(); ++index)
        host_values(index) = requests[index].likelihood;
    Kokkos::deep_copy(device_offsets, host_offsets);
    Kokkos::deep_copy(device_sorted, host_sorted);
    Kokkos::deep_copy(device_values, host_values);
    Kokkos::deep_copy(device_best_second,
                      -std::numeric_limits<double>::infinity());
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(request_count);
    device_batch.bind("offsets", device_offsets);
    device_batch.bind("sorted", device_sorted);
    device_batch.bind("values", device_values);
    device_batch.bind("best_second", device_best_second);
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const auto execute_start = std::chrono::steady_clock::now();
    Kokkos::parallel_for("fastgatk_read_allele_uncertainty", RangePolicy(0, row_count),
        KOKKOS_LAMBDA(const std::size_t row) {
            double best = -std::numeric_limits<double>::infinity();
            double second = -std::numeric_limits<double>::infinity();
            for (std::size_t cursor = device_offsets(row);
                 cursor < device_offsets(row + 1); ++cursor) {
                const auto request = device_sorted(cursor);
                const auto value = device_values(request);
                // Keep the same strict comparison order as Host
                // AlleleLikelihoods: equal values preserve the first allele
                // as best and become the second-best tie.
                if (value > best) {
                    second = best;
                    best = value;
                } else if (value > second) {
                    second = value;
                }
            }
            device_best_second(row, 0) = best;
            device_best_second(row, 1) = second;
        });
    ExecSpace().fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(execute_end - execute_start).count();
    auto host_best_second = Kokkos::create_mirror_view(device_best_second);
    Kokkos::deep_copy(host_best_second, device_best_second);
    for (std::size_t row = 0; row < row_count; ++row) {
        result.best_second_by_row[row * 2] = host_best_second(row, 0);
        result.best_second_by_row[row * 2 + 1] = host_best_second(row, 1);
    }
    return result;
}

ReadAlleleBestResult reduce_read_allele_best_kokkos(
    const std::vector<ReadAlleleBestRequest>& requests,
    const std::size_t row_count) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (row_count > std::numeric_limits<std::size_t>::max() / 2)
        throw std::invalid_argument("read allele best row count overflows output dimensions");
    if (row_count > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("read allele best row count exceeds uint32 row ids");
    for (const auto& request : requests) {
        if (request.row_id >= row_count)
            throw std::invalid_argument("read allele best row id exceeds row count");
        if (!std::isfinite(request.likelihood))
            throw std::invalid_argument("read allele best likelihood must be finite");
    }

    constexpr auto kMissing = std::numeric_limits<std::uint32_t>::max();
    ReadAlleleBestResult result;
    result.best_second_allele_by_row.assign(row_count * 2, kMissing);
    result.best_second_likelihood_by_row.assign(
        row_count * 2, -std::numeric_limits<double>::infinity());
    result.execution_space = ExecSpace::name();

    using OffsetView = Kokkos::View<std::size_t*, MemorySpace>;
    using SortedIndexView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using AlleleView = Kokkos::View<std::uint32_t*, MemorySpace>;
    using ValueView = Kokkos::View<double*, MemorySpace>;
    using BestIndexView = Kokkos::View<std::uint32_t**, Kokkos::LayoutRight, MemorySpace>;
    using BestValueView = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

    const auto request_count = requests.size();
    std::vector<std::uint32_t> sorted(request_count);
    std::iota(sorted.begin(), sorted.end(), 0U);
    // Stable allele order is part of the GATK contract: REF is emitted first,
    // followed by concrete ALT order.  The original index is the final tie
    // breaker so duplicate requests remain deterministic as well.
    std::stable_sort(sorted.begin(), sorted.end(), [&](const auto left, const auto right) {
        if (requests[left].row_id != requests[right].row_id)
            return requests[left].row_id < requests[right].row_id;
        if (requests[left].allele_index != requests[right].allele_index)
            return requests[left].allele_index < requests[right].allele_index;
        return left < right;
    });
    std::vector<std::size_t> offsets(row_count + 1, 0);
    for (const auto index : sorted)
        ++offsets[requests[index].row_id + 1];
    for (std::size_t row = 0; row < row_count; ++row)
        offsets[row + 1] += offsets[row];

    fastgatk::core::HostBatch host_batch("pairhmm-read-allele-best-v1");
    host_batch.records = request_count;
    host_batch.bytes = request_count *
        (sizeof(double) + 2 * sizeof(std::uint32_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm_read_allele_best");
    plan.begin_prepare(host_batch);
    OffsetView device_offsets("read_allele_best_offsets", offsets.size());
    SortedIndexView device_sorted("read_allele_best_sorted", sorted.size());
    AlleleView device_alleles("read_allele_best_alleles", request_count);
    ValueView device_values("read_allele_best_values", request_count);
    BestIndexView device_best("read_allele_best_indices", row_count, 2);
    BestValueView device_best_values("read_allele_best_values_out", row_count, 2);
    auto host_offsets = Kokkos::create_mirror_view(device_offsets);
    auto host_sorted = Kokkos::create_mirror_view(device_sorted);
    auto host_alleles = Kokkos::create_mirror_view(device_alleles);
    auto host_values = Kokkos::create_mirror_view(device_values);
    for (std::size_t row = 0; row < offsets.size(); ++row)
        host_offsets(row) = offsets[row];
    for (std::size_t index = 0; index < sorted.size(); ++index)
        host_sorted(index) = sorted[index];
    for (std::size_t index = 0; index < requests.size(); ++index) {
        host_alleles(index) = requests[index].allele_index;
        host_values(index) = requests[index].likelihood;
    }
    Kokkos::deep_copy(device_offsets, host_offsets);
    Kokkos::deep_copy(device_sorted, host_sorted);
    Kokkos::deep_copy(device_alleles, host_alleles);
    Kokkos::deep_copy(device_values, host_values);
    Kokkos::deep_copy(device_best, kMissing);
    Kokkos::deep_copy(device_best_values,
                      -std::numeric_limits<double>::infinity());
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(request_count);
    device_batch.bind("offsets", device_offsets);
    device_batch.bind("sorted", device_sorted);
    device_batch.bind("alleles", device_alleles);
    device_batch.bind("values", device_values);
    device_batch.bind("best", device_best);
    device_batch.bind("best_values", device_best_values);
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const auto execute_start = std::chrono::steady_clock::now();
    Kokkos::parallel_for("fastgatk_read_allele_best", RangePolicy(0, row_count),
        KOKKOS_LAMBDA(const std::size_t row) {
            double best = -std::numeric_limits<double>::infinity();
            double second = -std::numeric_limits<double>::infinity();
            std::uint32_t best_allele = kMissing;
            std::uint32_t second_allele = kMissing;
            for (std::size_t cursor = device_offsets(row);
                 cursor < device_offsets(row + 1); ++cursor) {
                const auto request = device_sorted(cursor);
                const auto value = device_values(request);
                const auto allele = device_alleles(request);
                // Strict comparisons plus sorted allele order reproduce the
                // Java searchBestAllele() tie behavior: equal likelihoods
                // retain the earlier REF/ALT entry as best and become the
                // second-best value.
                if (value > best) {
                    second = best;
                    second_allele = best_allele;
                    best = value;
                    best_allele = allele;
                } else if (value > second) {
                    second = value;
                    second_allele = allele;
                }
            }
            device_best(row, 0) = best_allele;
            device_best(row, 1) = second_allele;
            device_best_values(row, 0) = best;
            device_best_values(row, 1) = second;
        });
    ExecSpace().fence();
    const auto execute_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(execute_end - execute_start).count();
    auto host_best = Kokkos::create_mirror_view(device_best);
    auto host_best_values = Kokkos::create_mirror_view(device_best_values);
    Kokkos::deep_copy(host_best, device_best);
    Kokkos::deep_copy(host_best_values, device_best_values);
    for (std::size_t row = 0; row < row_count; ++row) {
        result.best_second_allele_by_row[row * 2] = host_best(row, 0);
        result.best_second_allele_by_row[row * 2 + 1] = host_best(row, 1);
        result.best_second_likelihood_by_row[row * 2] = host_best_values(row, 0);
        result.best_second_likelihood_by_row[row * 2 + 1] = host_best_values(row, 1);
    }
    return result;
}

KokkosBatchResult compute_kokkos(const std::vector<PairInput>& records,
                                 const PairIndexBatch& pairs, int iterations,
                                 bool tristate_correction) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");

    std::size_t read_length = 0, haplotype_length = 0;
    validate(records, pairs, read_length, haplotype_length);
    const std::size_t record_count = records.size();
    const std::size_t pair_count = pairs.read_ids.size();
    const auto& tables = gatk_tables();
    fastgatk::core::HostBatch host_batch("pairhmm-v1");
    host_batch.records = pair_count;
    for (const auto& input : records) {
        host_batch.bytes += input.haplotype.size() * sizeof(std::uint8_t);
        host_batch.bytes += input.read.size() * sizeof(std::uint8_t);
        host_batch.bytes += (input.read_qual.size() + input.insertion_gop.size() +
                             input.deletion_gop.size() + input.gap_continuation.size()) * sizeof(std::uint8_t);
    }
    host_batch.bytes += (pairs.read_ids.size() + pairs.haplotype_ids.size()) * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm");
    plan.begin_prepare(host_batch);

    View2D read_bases("read_bases", read_length, record_count);
    View2D hap_bases("hap_bases", haplotype_length, record_count);
    View2D mm_values("mm", read_length, record_count);
    View2D im_values("im", read_length, record_count);
    View2D mi_values("mi", read_length, record_count);
    View2D ii_values("ii", read_length, record_count);
    View2D md_values("md", read_length, record_count);
    View2D prior_match("prior_match", read_length, record_count);
    View2D prior_mismatch("prior_mismatch", read_length, record_count);
    IndexView read_ids("pair_read_ids", pair_count);
    IndexView hap_ids("pair_haplotype_ids", pair_count);
    ResultView likelihoods("likelihoods", pair_count);
    ResultView scaled_sums("scaled_sums", pair_count);

    auto h_read_bases = mirror(read_bases);
    auto h_hap_bases = mirror(hap_bases);
    auto h_mm = mirror(mm_values); auto h_im = mirror(im_values);
    auto h_mi = mirror(mi_values); auto h_ii = mirror(ii_values);
    auto h_md = mirror(md_values); auto h_pm = mirror(prior_match);
    auto h_px = mirror(prior_mismatch);
    auto h_read_ids = mirror(read_ids); auto h_hap_ids = mirror(hap_ids);

    for (std::size_t record = 0; record < record_count; ++record) {
        const PairInput& input = records[record];
        for (std::size_t j = 0; j < haplotype_length; ++j)
            h_hap_bases(j, record) = static_cast<double>(input.haplotype[j]);
        for (std::size_t i = 0; i < read_length; ++i) {
            h_read_bases(i, record) = static_cast<double>(input.read[i]);
            const double gcp = tables.quality[input.gap_continuation[i]];
            const double error = tables.quality[input.read_qual[i]];
            h_mm(i, record) = tables.match[static_cast<std::size_t>(input.insertion_gop[i]) * 256 + input.deletion_gop[i]];
            h_im(i, record) = 1.0 - gcp;
            h_mi(i, record) = tables.quality[input.insertion_gop[i]];
            h_ii(i, record) = gcp;
            h_md(i, record) = tables.quality[input.deletion_gop[i]];
            h_pm(i, record) = 1.0 - error;
            h_px(i, record) = error / (tristate_correction ? kTristateCorrection : 1.0);
        }
    }
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        h_read_ids(pair) = pairs.read_ids[pair];
        h_hap_ids(pair) = pairs.haplotype_ids[pair];
    }

    Kokkos::deep_copy(read_bases, h_read_bases); Kokkos::deep_copy(hap_bases, h_hap_bases);
    Kokkos::deep_copy(mm_values, h_mm); Kokkos::deep_copy(im_values, h_im);
    Kokkos::deep_copy(mi_values, h_mi); Kokkos::deep_copy(ii_values, h_ii);
    Kokkos::deep_copy(md_values, h_md); Kokkos::deep_copy(prior_match, h_pm);
    Kokkos::deep_copy(prior_mismatch, h_px); Kokkos::deep_copy(read_ids, h_read_ids);
    Kokkos::deep_copy(hap_ids, h_hap_ids);
    ExecSpace().fence();

    fastgatk::core::DeviceBatch<ExecSpace> device_batch(pair_count);
    device_batch.bind("read_bases", read_bases); device_batch.bind("hap_bases", hap_bases);
    device_batch.bind("mm", mm_values); device_batch.bind("im", im_values);
    device_batch.bind("mi", mi_values); device_batch.bind("ii", ii_values);
    device_batch.bind("md", md_values); device_batch.bind("prior_match", prior_match);
    device_batch.bind("prior_mismatch", prior_mismatch);
    device_batch.bind("pair_read_ids", read_ids); device_batch.bind("pair_haplotype_ids", hap_ids);

    KokkosBatchResult result;
    result.simd_width = Simd::size();
    result.execution_space = ExecSpace::name();
    result.likelihoods.resize(pair_count);
    result.scaled_sums.resize(pair_count);

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    using WorkspaceView = Kokkos::View<Simd**, Kokkos::LayoutRight, MemorySpace>;
    const std::size_t stride = haplotype_length + 1;
    const std::size_t groups = (pair_count + Simd::size() - 1) / Simd::size();
    const double initial = std::pow(2.0, 1020.0) / static_cast<double>(haplotype_length);
    const double initial_log10 = gatk_strict_log10(std::pow(2.0, 1020.0));
    // M/I/D are updated in place from left to right. The overwritten previous-row
    // cell is kept in a SIMD register for the next diagonal dependency, cutting
    // the workspace from six rows to three without changing arithmetic order.
    WorkspaceView workspace("pairhmm_workspace", groups, 3 * stride);
    device_batch.bind("workspace", workspace);
    device_batch.bind("likelihoods", likelihoods);
    device_batch.bind("scaled_sums", scaled_sums);
    RangePolicy policy(0, groups);
    ExecSpace().fence();
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    auto launch = [&] {
        Kokkos::parallel_for("fastgatk_pairhmm", policy, KOKKOS_LAMBDA(const std::size_t group) {
            const std::size_t base = group * Simd::size();
            const bool full_group = base + Simd::size() <= pair_count;
            const std::uint32_t first_read = read_ids(base);
            const std::uint32_t first_hap = hap_ids(base);
            bool uniform_read = full_group;
            bool contiguous_read = full_group;
            bool uniform_hap = full_group;
            bool contiguous_hap = full_group;
            for (std::size_t lane = 1; lane < Simd::size() && full_group; ++lane) {
                uniform_read = uniform_read && read_ids(base + lane) == first_read;
                contiguous_read = contiguous_read && read_ids(base + lane) == first_read + lane;
                uniform_hap = uniform_hap && hap_ids(base + lane) == first_hap;
                contiguous_hap = contiguous_hap && hap_ids(base + lane) == first_hap + lane;
            }
            // These loaders recognize the two dominant layouts without exposing
            // an ISA: independent pairs are contiguous, while a matrix8 SIMD
            // group broadcasts one read across consecutive haplotypes.
            const auto load_read = [&](const View2D& source, std::size_t row) -> Simd {
                if (uniform_read) return Simd(source(row, first_read));
                if (contiguous_read)
                    return Simd(&source(row, first_read), Kokkos::Experimental::simd_flag_default);
                return Simd([=](auto lane) -> double {
                    const std::size_t pair = base + static_cast<std::size_t>(lane);
                    return pair < pair_count ? source(row, read_ids(pair)) : 0.0;
                });
            };
            const auto load_hap = [&](const View2D& source, std::size_t row) -> Simd {
                if (uniform_hap) return Simd(source(row, first_hap));
                if (contiguous_hap)
                    return Simd(&source(row, first_hap), Kokkos::Experimental::simd_flag_default);
                return Simd([=](auto lane) -> double {
                    const std::size_t pair = base + static_cast<std::size_t>(lane);
                    return pair < pair_count ? source(row, hap_ids(pair)) : 0.0;
                });
            };
            Simd* const dp = &workspace(group, 0);
            const Simd zero(0.0);
            const Simd initial_v(initial);
            for (std::size_t j = 0; j <= haplotype_length; ++j) {
                dp[j] = zero;
                dp[stride + j] = zero;
                dp[2 * stride + j] = initial_v;
            }

            for (std::size_t i = 1; i <= read_length; ++i) {
                Simd previous_m = dp[0];
                Simd previous_i = dp[stride];
                Simd previous_d = dp[2 * stride];
                dp[0] = zero;
                dp[stride] = zero;
                dp[2 * stride] = zero;

                const Simd mm = load_read(mm_values, i - 1);
                const Simd im = load_read(im_values, i - 1);
                const Simd mi = load_read(mi_values, i - 1);
                const Simd ii = load_read(ii_values, i - 1);
                const Simd md = load_read(md_values, i - 1);
                const Simd p_match = load_read(prior_match, i - 1);
                const Simd p_mismatch = load_read(prior_mismatch, i - 1);
                const Simd read_base_v = load_read(read_bases, i - 1);

                for (std::size_t j = 1; j <= haplotype_length; ++j) {
                    const Simd hap_base_v = load_hap(hap_bases, j - 1);
                    const auto matches = (read_base_v == hap_base_v) |
                                         (read_base_v == Simd(static_cast<double>('N'))) |
                                         (hap_base_v == Simd(static_cast<double>('N')));
                    const Simd prior = Kokkos::Experimental::condition(matches, p_match, p_mismatch);
                    const Simd old_m = dp[j];
                    const Simd old_i = dp[stride + j];
                    const Simd old_d = dp[2 * stride + j];
                    const Simd mv_ab = previous_m * mm + previous_i * im;
                    dp[j] = prior * (mv_ab + previous_d * im);
                    dp[stride + j] = old_m * mi + old_i * ii;
                    dp[2 * stride + j] = dp[j - 1] * md + dp[2 * stride + j - 1] * ii;
                    previous_m = old_m;
                    previous_i = old_i;
                    previous_d = old_d;
                }
            }

            Simd sum(0.0);
            for (std::size_t j = 1; j <= haplotype_length; ++j)
                sum = sum + (dp[j] + dp[stride + j]);
            const Simd values([&](auto lane) -> double {
                return gatk_strict_log10(sum[static_cast<std::size_t>(lane)]) - initial_log10;
            });
            for (std::size_t lane = 0; lane < Simd::size(); ++lane) {
                const std::size_t pair = base + lane;
                if (pair < pair_count) {
                    scaled_sums(pair) = sum[lane];
                    likelihoods(pair) = values[lane];
                }
            }
        });
    };

    plan.begin_execute();
    for (int iteration = 0; iteration < iterations; ++iteration) launch();
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    result.pairs_per_second = static_cast<double>(pair_count) * iterations / result.seconds;

    auto h_likelihoods = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), likelihoods);
    auto h_scaled_sums = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scaled_sums);
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        result.scaled_sums[pair] = h_scaled_sums(pair);
        result.likelihoods[pair] = h_likelihoods(pair);
        result.checksum += h_likelihoods(pair);
    }
    // See fragment_aggregation comment for the slab-growth rationale.
    release_views(read_bases, hap_bases, mm_values, im_values, mi_values, ii_values,
                  md_values, prior_match, prior_mismatch, read_ids, hap_ids,
                  likelihoods, scaled_sums, workspace);
    return result;
}

// GKL-compatible Float32 path used by the regular base-space Mutect2
// pipeline.  It stores DP state and transition tables as float and uses
// Kokkos::Experimental::simd<float>, matching GKL's default Float32 context.
// The initial DP scale is reduced to 2^120 (rather than Java's 2^1020) and
// removed in the same Float32 operation at the end; this keeps the recurrence
// finite in IEEE-754 float while preserving its scale invariance.
KokkosBatchResult compute_kokkos_float_records(
    const std::vector<PairInput>& records,
    const PairIndexBatch& pairs,
    int iterations,
    bool tristate_correction = true) {
    using FloatSimd = Kokkos::Experimental::simd<float>;
    using FloatView2D = Kokkos::View<float**, Kokkos::LayoutRight, MemorySpace>;
    using FloatWorkspaceView = Kokkos::View<FloatSimd**, Kokkos::LayoutRight, MemorySpace>;
    using FloatResultView = Kokkos::View<float*, MemorySpace>;
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");

    std::size_t read_length = 0, haplotype_length = 0;
    validate(records, pairs, read_length, haplotype_length);
    const std::size_t record_count = records.size();
    const std::size_t pair_count = pairs.read_ids.size();
    const auto& tables = gkl_float_tables();
    fastgatk::core::HostBatch host_batch("pairhmm-float32-v1");
    host_batch.records = pair_count;
    for (const auto& input : records) {
        host_batch.bytes += input.haplotype.size() * sizeof(std::uint8_t);
        host_batch.bytes += input.read.size() * sizeof(std::uint8_t);
        host_batch.bytes += (input.read_qual.size() + input.insertion_gop.size() +
                             input.deletion_gop.size() + input.gap_continuation.size()) * sizeof(std::uint8_t);
    }
    host_batch.bytes += (pairs.read_ids.size() + pairs.haplotype_ids.size()) * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm-float32");
    plan.begin_prepare(host_batch);

    FloatView2D read_bases("float_read_bases", read_length, record_count);
    FloatView2D hap_bases("float_hap_bases", haplotype_length, record_count);
    FloatView2D mm_values("float_mm", read_length, record_count);
    FloatView2D im_values("float_im", read_length, record_count);
    FloatView2D mi_values("float_mi", read_length, record_count);
    FloatView2D ii_values("float_ii", read_length, record_count);
    FloatView2D md_values("float_md", read_length, record_count);
    FloatView2D prior_match("float_prior_match", read_length, record_count);
    FloatView2D prior_mismatch("float_prior_mismatch", read_length, record_count);
    IndexView read_ids("float_pair_read_ids", pair_count);
    IndexView hap_ids("float_pair_hap_ids", pair_count);
    FloatResultView likelihoods("float_likelihoods", pair_count);
    FloatResultView scaled_sums("float_scaled_sums", pair_count);

    // Host mirrors are filled in the loop below before the deep_copy,
    // so the implicit memset(0) inside `create_mirror_view` is pure
    // overhead.  Use `WithoutInitializing` to skip it.  GPU-safe: a
    // Kokkos View construction with a runtime property; works on every
    // backend.  Saves `read_length × record_count` writes per mirror
    // × 9 mirrors per bucket × ~30 buckets per region.
    auto h_read_bases = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), read_bases);
    auto h_hap_bases = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), hap_bases);
    auto h_mm = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), mm_values);
    auto h_im = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), im_values);
    auto h_mi = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), mi_values);
    auto h_ii = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), ii_values);
    auto h_md = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), md_values);
    auto h_pm = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), prior_match);
    auto h_px = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing), prior_mismatch);
    auto h_read_ids = Kokkos::create_mirror_view(read_ids);
    auto h_hap_ids = Kokkos::create_mirror_view(hap_ids);
    for (std::size_t record = 0; record < record_count; ++record) {
        const PairInput& input = records[record];
        for (std::size_t j = 0; j < haplotype_length; ++j)
            h_hap_bases(j, record) = static_cast<float>(input.haplotype[j]);
        for (std::size_t i = 0; i < read_length; ++i) {
            h_read_bases(i, record) = static_cast<float>(input.read[i]);
            const auto gcp = tables.quality[input.gap_continuation[i] & 127U];
            const auto error = tables.quality[input.read_qual[i] & 127U];
            h_mm(i, record) = gkl_float_match_to_match(
                tables, input.insertion_gop[i], input.deletion_gop[i]);
            // Do the transitions in Float32, as Context<float> and GKL's
            // AVX arithmetic do.  Computing them in double then casting is
            // observably different at the final Mutect2 NALOD boundary.
            h_im(i, record) = 1.0F - gcp;
            h_mi(i, record) = tables.quality[input.insertion_gop[i] & 127U];
            h_ii(i, record) = gcp;
            h_md(i, record) = tables.quality[input.deletion_gop[i] & 127U];
            h_pm(i, record) = 1.0F - error;
            h_px(i, record) = error / (tristate_correction ? 3.0F : 1.0F);
        }
    }
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        h_read_ids(pair) = pairs.read_ids[pair];
        h_hap_ids(pair) = pairs.haplotype_ids[pair];
    }
    Kokkos::deep_copy(read_bases, h_read_bases);
    Kokkos::deep_copy(hap_bases, h_hap_bases);
    Kokkos::deep_copy(mm_values, h_mm);
    Kokkos::deep_copy(im_values, h_im);
    Kokkos::deep_copy(mi_values, h_mi);
    Kokkos::deep_copy(ii_values, h_ii);
    Kokkos::deep_copy(md_values, h_md);
    Kokkos::deep_copy(prior_match, h_pm);
    Kokkos::deep_copy(prior_mismatch, h_px);
    Kokkos::deep_copy(read_ids, h_read_ids);
    Kokkos::deep_copy(hap_ids, h_hap_ids);
    ExecSpace().fence();

    fastgatk::core::DeviceBatch<ExecSpace> device_batch(pair_count);
    device_batch.bind("float_read_bases", read_bases);
    device_batch.bind("float_hap_bases", hap_bases);
    device_batch.bind("float_mm", mm_values);
    device_batch.bind("float_im", im_values);
    device_batch.bind("float_mi", mi_values);
    device_batch.bind("float_ii", ii_values);
    device_batch.bind("float_md", md_values);
    device_batch.bind("float_prior_match", prior_match);
    device_batch.bind("float_prior_mismatch", prior_mismatch);
    device_batch.bind("float_pair_read_ids", read_ids);
    device_batch.bind("float_pair_hap_ids", hap_ids);

    KokkosBatchResult result;
    result.precision = "float32";
    result.simd_width = FloatSimd::size();
    result.execution_space = ExecSpace::name();
    result.likelihoods.resize(pair_count);
    result.scaled_sums.resize(pair_count);

    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const std::size_t stride = haplotype_length + 1;
    const std::size_t groups = (pair_count + FloatSimd::size() - 1) / FloatSimd::size();
    constexpr float initial_scale = 1.329227995784916e36F;  // 2^120
    const float initial = initial_scale / static_cast<float>(haplotype_length);
    // GKL returns `(double)(log10f(sum) - log10f(2^120))` for its normal
    // Float32 path.  Keep both logarithms and their subtraction in Float32.
    const float initial_log10 = Kokkos::log10(initial_scale);
    FloatWorkspaceView workspace("float_pairhmm_workspace", groups, 3 * stride);
    device_batch.bind("float_workspace", workspace);
    device_batch.bind("float_likelihoods", likelihoods);
    device_batch.bind("float_scaled_sums", scaled_sums);
    RangePolicy policy(0, groups);
    ExecSpace().fence();
    plan.end_prepare(device_batch);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    auto launch = [&] {
        Kokkos::parallel_for("fastgatk_pairhmm_float32", policy,
            KOKKOS_LAMBDA(const std::size_t group) {
                const std::size_t base = group * FloatSimd::size();
                const bool full_group = base + FloatSimd::size() <= pair_count;
                const std::uint32_t first_read = read_ids(base);
                const std::uint32_t first_hap = hap_ids(base);
                bool uniform_read = full_group;
                bool contiguous_read = full_group;
                bool uniform_hap = full_group;
                bool contiguous_hap = full_group;
                for (std::size_t lane = 1; lane < FloatSimd::size() && full_group; ++lane) {
                    uniform_read = uniform_read && read_ids(base + lane) == first_read;
                    contiguous_read = contiguous_read && read_ids(base + lane) == first_read + lane;
                    uniform_hap = uniform_hap && hap_ids(base + lane) == first_hap;
                    contiguous_hap = contiguous_hap && hap_ids(base + lane) == first_hap + lane;
                }
                const auto load_read = [&](const FloatView2D& source, std::size_t row) -> FloatSimd {
                    if (uniform_read) return FloatSimd(source(row, first_read));
                    if (contiguous_read)
                        return FloatSimd(&source(row, first_read), Kokkos::Experimental::simd_flag_default);
                    return FloatSimd([=](auto lane) -> float {
                        const std::size_t pair = base + static_cast<std::size_t>(lane);
                        return pair < pair_count ? source(row, read_ids(pair)) : 0.0F;
                    });
                };
                const auto load_hap = [&](const FloatView2D& source, std::size_t row) -> FloatSimd {
                    if (uniform_hap) return FloatSimd(source(row, first_hap));
                    if (contiguous_hap)
                        return FloatSimd(&source(row, first_hap), Kokkos::Experimental::simd_flag_default);
                    return FloatSimd([=](auto lane) -> float {
                        const std::size_t pair = base + static_cast<std::size_t>(lane);
                        return pair < pair_count ? source(row, hap_ids(pair)) : 0.0F;
                    });
                };
                FloatSimd* const dp = &workspace(group, 0);
                const FloatSimd zero(0.0F);
                const FloatSimd initial_v(initial);
                for (std::size_t j = 0; j <= haplotype_length; ++j) {
                    dp[j] = zero;
                    dp[stride + j] = zero;
                    dp[2 * stride + j] = initial_v;
                }
                for (std::size_t i = 1; i <= read_length; ++i) {
                    FloatSimd previous_m = dp[0];
                    FloatSimd previous_i = dp[stride];
                    FloatSimd previous_d = dp[2 * stride];
                    dp[0] = zero;
                    dp[stride] = zero;
                    dp[2 * stride] = zero;
                    const FloatSimd mm = load_read(mm_values, i - 1);
                    const FloatSimd im = load_read(im_values, i - 1);
                    const FloatSimd mi = load_read(mi_values, i - 1);
                    const FloatSimd ii = load_read(ii_values, i - 1);
                    const FloatSimd md = load_read(md_values, i - 1);
                    const FloatSimd p_match = load_read(prior_match, i - 1);
                    const FloatSimd p_mismatch = load_read(prior_mismatch, i - 1);
                    const FloatSimd read_base_v = load_read(read_bases, i - 1);
                    for (std::size_t j = 1; j <= haplotype_length; ++j) {
                        const FloatSimd hap_base_v = load_hap(hap_bases, j - 1);
                        const auto matches = (read_base_v == hap_base_v) |
                            (read_base_v == FloatSimd(static_cast<float>('N'))) |
                            (hap_base_v == FloatSimd(static_cast<float>('N')));
                        const FloatSimd prior = Kokkos::Experimental::condition(matches, p_match, p_mismatch);
                        const FloatSimd old_m = dp[j];
                        const FloatSimd old_i = dp[stride + j];
                        const FloatSimd old_d = dp[2 * stride + j];
                        const FloatSimd mv_ab = previous_m * mm + previous_i * im;
                        dp[j] = prior * (mv_ab + previous_d * im);
                        dp[stride + j] = old_m * mi + old_i * ii;
                        dp[2 * stride + j] = dp[j - 1] * md + dp[2 * stride + j - 1] * ii;
                        previous_m = old_m;
                        previous_i = old_i;
                        previous_d = old_d;
                    }
                }
                FloatSimd sum(0.0F);
                for (std::size_t j = 1; j <= haplotype_length; ++j)
                    sum = sum + (dp[j] + dp[stride + j]);
                const FloatSimd values([&](auto lane) -> float {
                    return Kokkos::log10(sum[static_cast<std::size_t>(lane)]) - initial_log10;
                });
                for (std::size_t lane = 0; lane < FloatSimd::size(); ++lane) {
                    const std::size_t pair = base + lane;
                    if (pair < pair_count) {
                        scaled_sums(pair) = sum[lane];
                        likelihoods(pair) = values[lane];
                    }
                }
            });
    };

    plan.begin_execute();
    for (int iteration = 0; iteration < iterations; ++iteration) launch();
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    result.pairs_per_second = static_cast<double>(pair_count) * iterations / result.seconds;

    auto h_likelihoods = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), likelihoods);
    auto h_scaled_sums = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scaled_sums);
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        result.scaled_sums[pair] = static_cast<double>(h_scaled_sums(pair));
        result.likelihoods[pair] = static_cast<double>(h_likelihoods(pair));
        result.checksum += result.likelihoods[pair];
    }
    // See fragment_aggregation comment for the slab-growth rationale.
    release_views(read_bases, hap_bases, mm_values, im_values, mi_values, ii_values,
                  md_values, prior_match, prior_mismatch, read_ids, hap_ids,
                  likelihoods, scaled_sums, workspace);
    return result;
}

KokkosBatchResult compute_kokkos_bucketed(
    const std::vector<PairHmmRead>& reads,
    const std::vector<PairHmmHaplotype>& haplotypes,
    const std::vector<PairHmmRequest>& requests,
    int iterations,
    PairHmmPrecision precision) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");

    KokkosBatchResult aggregate;
    aggregate.precision = precision == PairHmmPrecision::Float32 ? "float32" : "double";
    aggregate.simd_width = precision == PairHmmPrecision::Float32
        ? Kokkos::Experimental::simd<float>::size() : Simd::size();
    aggregate.execution_space = ExecSpace::name();
    aggregate.likelihoods.assign(requests.size(), 0.0);
    aggregate.scaled_sums.assign(requests.size(), 0.0);
    if (requests.empty()) return aggregate;
    using BucketKey = std::pair<std::size_t, std::size_t>;
    std::map<BucketKey, std::vector<std::size_t>> buckets;
    for (std::size_t request = 0; request < requests.size(); ++request) {
        const auto& item = requests[request];
        if (item.read_id >= reads.size() || item.haplotype_id >= haplotypes.size())
            throw std::invalid_argument("bucketed PairHMM request index exceeds input");
        const auto& read = reads[item.read_id];
        const auto& haplotype = haplotypes[item.haplotype_id];
        if (read.bases.empty() || haplotype.bases.empty())
            throw std::invalid_argument("bucketed PairHMM input contains an empty sequence");
        if (read.qualities.size() != read.bases.size() ||
            read.insertion_gop.size() != read.bases.size() ||
            read.deletion_gop.size() != read.bases.size() ||
            read.gap_continuation.size() != read.bases.size())
            throw std::invalid_argument("bucketed PairHMM read quality arrays differ");
        buckets[{read.bases.size(), haplotype.bases.size()}].push_back(request);
    }

    for (const auto& [key, request_indices] : buckets) {
        const auto read_length = key.first;
        const auto haplotype_length = key.second;
        std::map<std::uint32_t, std::uint32_t> read_local;
        std::map<std::uint32_t, std::uint32_t> haplotype_local;
        std::vector<std::uint32_t> bucket_reads;
        std::vector<std::uint32_t> bucket_haplotypes;
        for (const auto request_index : request_indices) {
            const auto& request = requests[request_index];
            if (!read_local.contains(request.read_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_reads.size());
                read_local.emplace(request.read_id, local);
                bucket_reads.push_back(request.read_id);
            }
            if (!haplotype_local.contains(request.haplotype_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_haplotypes.size());
                haplotype_local.emplace(request.haplotype_id, local);
                bucket_haplotypes.push_back(request.haplotype_id);
            }
        }

        // Reorder `request_indices` so that consecutive SIMD-width lanes in
        // the Kokkos kernel share the same read_id and occupy consecutive
        // haplotype_ids.  This maximizes the `uniform_read` and
        // `contiguous_hap` fast paths in `load_read`/`load_hap`, which
        // let the kernel fall back to a single 16-lane broadcast per row
        // instead of a per-lane scalar gather.  Sort is stable on the
        // request order so the bucket's likelihood / checksum ordering
        // remains deterministic for callers (GATK parity oracles in
        // particular match by row index within a bucket).  Done before
        // building `read_local` / `haplotype_local` so the per-bucket
        // record indices stay consistent with the sorted request order.
        std::vector<std::size_t> sorted_indices(request_indices.begin(), request_indices.end());
        std::stable_sort(sorted_indices.begin(), sorted_indices.end(),
            [&](const std::size_t left, const std::size_t right) {
                const auto& l = requests[left];
                const auto& r = requests[right];
                if (l.read_id != r.read_id) return l.read_id < r.read_id;
                return l.haplotype_id < r.haplotype_id;
            });

        // One PairInput record can act as a read source or a haplotype source;
        // the unused side is filled with deterministic N/quality defaults.
        const auto record_count = bucket_reads.size() + bucket_haplotypes.size();
        std::vector<PairInput> records;
        records.reserve(record_count);
        for (const auto read_id : bucket_reads) {
            const auto& source = reads[read_id];
            PairInput input;
            input.haplotype.assign(haplotype_length, static_cast<std::uint8_t>('N'));
            input.read = source.bases;
            input.read_qual = source.qualities;
            input.insertion_gop = source.insertion_gop;
            input.deletion_gop = source.deletion_gop;
            input.gap_continuation = source.gap_continuation;
            records.push_back(std::move(input));
        }
        const auto haplotype_base = bucket_reads.size();
        for (const auto haplotype_id : bucket_haplotypes) {
            const auto& source = haplotypes[haplotype_id];
            PairInput input;
            input.haplotype = source.bases;
            input.read.assign(read_length, static_cast<std::uint8_t>('N'));
            input.read_qual.assign(read_length, 30);
            input.insertion_gop.assign(read_length, 40);
            input.deletion_gop.assign(read_length, 40);
            input.gap_continuation.assign(read_length, 10);
            records.push_back(std::move(input));
        }

        PairIndexBatch pairs;
        pairs.read_ids.reserve(sorted_indices.size());
        pairs.haplotype_ids.reserve(sorted_indices.size());
        for (const auto request_index : sorted_indices) {
            const auto& request = requests[request_index];
            pairs.read_ids.push_back(read_local.at(request.read_id));
            pairs.haplotype_ids.push_back(static_cast<std::uint32_t>(haplotype_base +
                                                                       haplotype_local.at(request.haplotype_id)));
        }
        auto bucket_result = precision == PairHmmPrecision::Float32
            ? compute_kokkos_float_records(records, pairs, iterations)
            : compute_kokkos(records, pairs, iterations);
        if (precision == PairHmmPrecision::Float32) {
            // GKL's JNI wrapper first evaluates every pair in Context<float>,
            // but it never returns a Float32 result whose *scaled* sum is
            // below MIN_ACCEPTED (1e-28f).  Those pairs are recomputed in
            // Context<double> before `log10` is taken, which avoids an
            // otherwise silent Float32 underflow on very poor alignments.
            // Preserve that per-pair selection in the Kokkos path: the Host
            // only compacts request IDs, while both recurrences remain
            // Kokkos kernels and no read/assembly ownership changes.
            constexpr double kGklMinimumAcceptedScaledSum = 1.0e-28;
            PairIndexBatch fallback_pairs;
            std::vector<std::size_t> fallback_indices;
            fallback_pairs.read_ids.reserve(pairs.read_ids.size());
            fallback_pairs.haplotype_ids.reserve(pairs.haplotype_ids.size());
            fallback_indices.reserve(pairs.read_ids.size());
            for (std::size_t local = 0; local < pairs.read_ids.size(); ++local) {
                if (bucket_result.scaled_sums[local] >= kGklMinimumAcceptedScaledSum) continue;
                fallback_pairs.read_ids.push_back(pairs.read_ids[local]);
                fallback_pairs.haplotype_ids.push_back(pairs.haplotype_ids[local]);
                fallback_indices.push_back(local);
            }
            if (!fallback_indices.empty()) {
                const auto fallback = compute_kokkos(records, fallback_pairs, iterations);
                bucket_result.prepare_seconds += fallback.prepare_seconds;
                bucket_result.seconds += fallback.seconds;
                for (std::size_t fallback_index = 0;
                     fallback_index < fallback_indices.size(); ++fallback_index) {
                    const auto local = fallback_indices[fallback_index];
                    bucket_result.likelihoods[local] = fallback.likelihoods[fallback_index];
                    bucket_result.scaled_sums[local] = fallback.scaled_sums[fallback_index];
                }
                bucket_result.checksum = 0.0;
                for (const auto value : bucket_result.likelihoods)
                    bucket_result.checksum += value;
            }
        }
        aggregate.simd_width = bucket_result.simd_width;
        aggregate.prepare_seconds += bucket_result.prepare_seconds;
        aggregate.seconds += bucket_result.seconds;
        // `bucket_result.likelihoods` is indexed by the sorted order
        // passed to `compute_kokkos(_float_records)` via `pairs`.  Map
        // each `local` position back to its original request index, then
        // write into the aggregate's request-aligned output arrays.
        for (std::size_t local = 0; local < sorted_indices.size(); ++local) {
            const auto original = sorted_indices[local];
            aggregate.likelihoods[original] = bucket_result.likelihoods[local];
            aggregate.scaled_sums[original] = bucket_result.scaled_sums[local];
            aggregate.checksum += bucket_result.likelihoods[local];
        }
    }
    aggregate.pairs_per_second = aggregate.seconds > 0.0
        ? static_cast<double>(requests.size()) * iterations / aggregate.seconds : 0.0;
    // The per-bucket Kokkos Views were released back to the OpenMP memory
    // pool by `release_views`.  glibc keeps the pool's top chunks in the
    // process RSS even after that release; `malloc_trim(0)` returns them
    // to the OS so the running peak accurately reflects the working set
    // after each region.  A no-op on non-glibc platforms.
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
    return aggregate;
}

KokkosBatchResult compute_kokkos_full_matrix(
    const std::vector<PairHmmRead>& reads,
    const std::vector<PairHmmHaplotype>& haplotypes,
    int iterations) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");

    KokkosBatchResult result;
    result.simd_width = Simd::size();
    result.execution_space = ExecSpace::name();
    result.matrix_rows = reads.size();
    result.matrix_columns = haplotypes.size();
    if (reads.empty() || haplotypes.empty()) return result;
    if (reads.size() > std::numeric_limits<std::uint32_t>::max() ||
        haplotypes.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("PairHMM full matrix dimensions exceed uint32 request IDs");
    if (reads.size() > std::numeric_limits<std::size_t>::max() / haplotypes.size())
        throw std::invalid_argument("PairHMM full matrix dimensions overflow");

    const auto pair_count = reads.size() * haplotypes.size();
    std::vector<PairHmmRequest> requests;
    requests.reserve(pair_count);
    for (std::size_t read = 0; read < reads.size(); ++read)
        for (std::size_t haplotype = 0; haplotype < haplotypes.size(); ++haplotype)
            requests.push_back(PairHmmRequest{
                static_cast<std::uint32_t>(read), static_cast<std::uint32_t>(haplotype)});

    // The full matrix is often invoked repeatedly for one assembly region.  A
    // local persistent plan keeps the implementation on the same lifecycle as
    // callers that explicitly own PersistentBucketPlan, while still allowing
    // this convenience API to remain value-oriented and thread-safe.
    PersistentBucketPlan plan(64);
    result = plan.execute(reads, haplotypes, requests, iterations);
    result.matrix_rows = reads.size();
    result.matrix_columns = haplotypes.size();
    return result;
}

namespace {

struct PersistentState {
    const std::size_t read_length;
    const std::size_t haplotype_length;
    const std::size_t record_count;
    const std::size_t pair_count;
    const std::size_t stride;
    const std::size_t groups;

    View2D read_bases;
    View2D hap_bases;
    View2D mm_values;
    View2D im_values;
    View2D mi_values;
    View2D ii_values;
    View2D md_values;
    View2D prior_match;
    View2D prior_mismatch;
    IndexView read_ids;
    IndexView hap_ids;
    ResultView likelihoods;
    ResultView scaled_sums;
    WorkspaceView workspace;
    bool warmed = false;

    PersistentState(std::size_t read_length_value, std::size_t haplotype_length_value,
                    std::size_t record_count_value, std::size_t pair_count_value)
        : read_length(read_length_value), haplotype_length(haplotype_length_value),
          record_count(record_count_value), pair_count(pair_count_value),
          stride(haplotype_length_value + 1),
          groups((pair_count_value + Simd::size() - 1) / Simd::size()),
          read_bases("persistent_read_bases", read_length_value, record_count_value),
          hap_bases("persistent_hap_bases", haplotype_length_value, record_count_value),
          mm_values("persistent_mm", read_length_value, record_count_value),
          im_values("persistent_im", read_length_value, record_count_value),
          mi_values("persistent_mi", read_length_value, record_count_value),
          ii_values("persistent_ii", read_length_value, record_count_value),
          md_values("persistent_md", read_length_value, record_count_value),
          prior_match("persistent_prior_match", read_length_value, record_count_value),
          prior_mismatch("persistent_prior_mismatch", read_length_value, record_count_value),
          read_ids("persistent_pair_read_ids", pair_count_value),
          hap_ids("persistent_pair_hap_ids", pair_count_value),
          likelihoods("persistent_likelihoods", pair_count_value),
          scaled_sums("persistent_scaled_sums", pair_count_value),
          workspace("persistent_pairhmm_workspace", groups, 3 * stride) {}
};

void prepare_persistent_state(PersistentState& state,
                              const std::vector<PairInput>& records,
                              const PairIndexBatch& pairs,
                              KokkosBatchResult& result) {
    const auto& tables = gatk_tables();
    fastgatk::core::HostBatch host_batch("pairhmm-v1-persistent");
    host_batch.records = pairs.read_ids.size();
    for (const auto& input : records) {
        host_batch.bytes += input.haplotype.size() * sizeof(std::uint8_t);
        host_batch.bytes += input.read.size() * sizeof(std::uint8_t);
        host_batch.bytes += (input.read_qual.size() + input.insertion_gop.size() +
                             input.deletion_gop.size() + input.gap_continuation.size()) * sizeof(std::uint8_t);
    }
    host_batch.bytes += (pairs.read_ids.size() + pairs.haplotype_ids.size()) * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm-persistent");
    plan.begin_prepare(host_batch);

    auto h_read_bases = Kokkos::create_mirror_view(state.read_bases);
    auto h_hap_bases = Kokkos::create_mirror_view(state.hap_bases);
    auto h_mm = Kokkos::create_mirror_view(state.mm_values);
    auto h_im = Kokkos::create_mirror_view(state.im_values);
    auto h_mi = Kokkos::create_mirror_view(state.mi_values);
    auto h_ii = Kokkos::create_mirror_view(state.ii_values);
    auto h_md = Kokkos::create_mirror_view(state.md_values);
    auto h_pm = Kokkos::create_mirror_view(state.prior_match);
    auto h_px = Kokkos::create_mirror_view(state.prior_mismatch);
    auto h_read_ids = Kokkos::create_mirror_view(state.read_ids);
    auto h_hap_ids = Kokkos::create_mirror_view(state.hap_ids);
    for (std::size_t record = 0; record < records.size(); ++record) {
        const PairInput& input = records[record];
        for (std::size_t j = 0; j < state.haplotype_length; ++j)
            h_hap_bases(j, record) = static_cast<double>(input.haplotype[j]);
        for (std::size_t i = 0; i < state.read_length; ++i) {
            h_read_bases(i, record) = static_cast<double>(input.read[i]);
            const double gcp = tables.quality[input.gap_continuation[i]];
            const double error = tables.quality[input.read_qual[i]];
            h_mm(i, record) = tables.match[static_cast<std::size_t>(input.insertion_gop[i]) * 256 + input.deletion_gop[i]];
            h_im(i, record) = 1.0 - gcp;
            h_mi(i, record) = tables.quality[input.insertion_gop[i]];
            h_ii(i, record) = gcp;
            h_md(i, record) = tables.quality[input.deletion_gop[i]];
            h_pm(i, record) = 1.0 - error;
            h_px(i, record) = error / kTristateCorrection;
        }
    }
    for (std::size_t pair = 0; pair < pairs.read_ids.size(); ++pair) {
        h_read_ids(pair) = pairs.read_ids[pair];
        h_hap_ids(pair) = pairs.haplotype_ids[pair];
    }
    Kokkos::deep_copy(state.read_bases, h_read_bases);
    Kokkos::deep_copy(state.hap_bases, h_hap_bases);
    Kokkos::deep_copy(state.mm_values, h_mm);
    Kokkos::deep_copy(state.im_values, h_im);
    Kokkos::deep_copy(state.mi_values, h_mi);
    Kokkos::deep_copy(state.ii_values, h_ii);
    Kokkos::deep_copy(state.md_values, h_md);
    Kokkos::deep_copy(state.prior_match, h_pm);
    Kokkos::deep_copy(state.prior_mismatch, h_px);
    Kokkos::deep_copy(state.read_ids, h_read_ids);
    Kokkos::deep_copy(state.hap_ids, h_hap_ids);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(state.pair_count);
    device.bind("read_bases", state.read_bases);
    device.bind("hap_bases", state.hap_bases);
    device.bind("mm", state.mm_values);
    device.bind("im", state.im_values);
    device.bind("mi", state.mi_values);
    device.bind("ii", state.ii_values);
    device.bind("md", state.md_values);
    device.bind("prior_match", state.prior_match);
    device.bind("prior_mismatch", state.prior_mismatch);
    device.bind("pair_read_ids", state.read_ids);
    device.bind("pair_haplotype_ids", state.hap_ids);
    device.bind("workspace", state.workspace);
    device.bind("likelihoods", state.likelihoods);
    device.bind("scaled_sums", state.scaled_sums);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
}

KokkosBatchResult execute_persistent_state(PersistentState& state, int iterations) {
    KokkosBatchResult result;
    result.simd_width = Simd::size();
    result.execution_space = ExecSpace::name();
    result.likelihoods.resize(state.pair_count);
    result.scaled_sums.resize(state.pair_count);
    const double initial = std::pow(2.0, 1020.0) / static_cast<double>(state.haplotype_length);
    const double initial_log10 = gatk_strict_log10(std::pow(2.0, 1020.0));
    const auto read_bases = state.read_bases;
    const auto hap_bases = state.hap_bases;
    const auto mm_values = state.mm_values;
    const auto im_values = state.im_values;
    const auto mi_values = state.mi_values;
    const auto ii_values = state.ii_values;
    const auto md_values = state.md_values;
    const auto prior_match = state.prior_match;
    const auto prior_mismatch = state.prior_mismatch;
    const auto read_ids = state.read_ids;
    const auto hap_ids = state.hap_ids;
    const auto likelihoods = state.likelihoods;
    const auto scaled_sums = state.scaled_sums;
    const auto workspace = state.workspace;
    const auto read_length = state.read_length;
    const auto haplotype_length = state.haplotype_length;
    const auto pair_count = state.pair_count;
    const auto stride = state.stride;
    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
    const RangePolicy policy(0, state.groups);
    auto launch = [&] {
        Kokkos::parallel_for("fastgatk_pairhmm_persistent", policy,
            KOKKOS_LAMBDA(const std::size_t group) {
                const std::size_t base = group * Simd::size();
                const bool full_group = base + Simd::size() <= pair_count;
                const std::uint32_t first_read = read_ids(base);
                const std::uint32_t first_hap = hap_ids(base);
                bool uniform_read = full_group;
                bool contiguous_read = full_group;
                bool uniform_hap = full_group;
                bool contiguous_hap = full_group;
                for (std::size_t lane = 1; lane < Simd::size() && full_group; ++lane) {
                    uniform_read = uniform_read && read_ids(base + lane) == first_read;
                    contiguous_read = contiguous_read && read_ids(base + lane) == first_read + lane;
                    uniform_hap = uniform_hap && hap_ids(base + lane) == first_hap;
                    contiguous_hap = contiguous_hap && hap_ids(base + lane) == first_hap + lane;
                }
                const auto load_read = [&](const View2D& source, std::size_t row) -> Simd {
                    if (uniform_read) return Simd(source(row, first_read));
                    if (contiguous_read)
                        return Simd(&source(row, first_read), Kokkos::Experimental::simd_flag_default);
                    return Simd([=](auto lane) -> double {
                        const std::size_t pair = base + static_cast<std::size_t>(lane);
                        return pair < pair_count ? source(row, read_ids(pair)) : 0.0;
                    });
                };
                const auto load_hap = [&](const View2D& source, std::size_t row) -> Simd {
                    if (uniform_hap) return Simd(source(row, first_hap));
                    if (contiguous_hap)
                        return Simd(&source(row, first_hap), Kokkos::Experimental::simd_flag_default);
                    return Simd([=](auto lane) -> double {
                        const std::size_t pair = base + static_cast<std::size_t>(lane);
                        return pair < pair_count ? source(row, hap_ids(pair)) : 0.0;
                    });
                };
                Simd* const dp = &workspace(group, 0);
                const Simd zero(0.0);
                const Simd initial_v(initial);
                for (std::size_t j = 0; j <= haplotype_length; ++j) {
                    dp[j] = zero;
                    dp[stride + j] = zero;
                    dp[2 * stride + j] = initial_v;
                }
                for (std::size_t i = 1; i <= read_length; ++i) {
                    Simd previous_m = dp[0];
                    Simd previous_i = dp[stride];
                    Simd previous_d = dp[2 * stride];
                    dp[0] = zero;
                    dp[stride] = zero;
                    dp[2 * stride] = zero;
                    const Simd mm = load_read(mm_values, i - 1);
                    const Simd im = load_read(im_values, i - 1);
                    const Simd mi = load_read(mi_values, i - 1);
                    const Simd ii = load_read(ii_values, i - 1);
                    const Simd md = load_read(md_values, i - 1);
                    const Simd p_match = load_read(prior_match, i - 1);
                    const Simd p_mismatch = load_read(prior_mismatch, i - 1);
                    const Simd read_base_v = load_read(read_bases, i - 1);
                    for (std::size_t j = 1; j <= haplotype_length; ++j) {
                        const Simd hap_base_v = load_hap(hap_bases, j - 1);
                        const auto matches = (read_base_v == hap_base_v) |
                            (read_base_v == Simd(static_cast<double>('N'))) |
                            (hap_base_v == Simd(static_cast<double>('N')));
                        const Simd prior = Kokkos::Experimental::condition(matches, p_match, p_mismatch);
                        const Simd old_m = dp[j];
                        const Simd old_i = dp[stride + j];
                        const Simd old_d = dp[2 * stride + j];
                        const Simd mv_ab = previous_m * mm + previous_i * im;
                        dp[j] = prior * (mv_ab + previous_d * im);
                        dp[stride + j] = old_m * mi + old_i * ii;
                        dp[2 * stride + j] = dp[j - 1] * md + dp[2 * stride + j - 1] * ii;
                        previous_m = old_m;
                        previous_i = old_i;
                        previous_d = old_d;
                    }
                }
                Simd sum(0.0);
                for (std::size_t j = 1; j <= haplotype_length; ++j)
                    sum = sum + (dp[j] + dp[stride + j]);
                const Simd values([&](auto lane) -> double {
                    return gatk_strict_log10(sum[static_cast<std::size_t>(lane)]) - initial_log10;
                });
                for (std::size_t lane = 0; lane < Simd::size(); ++lane) {
                    const std::size_t pair = base + lane;
                    if (pair < pair_count) {
                        scaled_sums(pair) = sum[lane];
                        likelihoods(pair) = values[lane];
                    }
                }
            });
    };
    if (!state.warmed) {
        state.warmed = true;
    }
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm-persistent-execute");
    plan.begin_execute();
    for (int iteration = 0; iteration < iterations; ++iteration) launch();
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    result.pairs_per_second = result.seconds > 0.0
        ? static_cast<double>(state.pair_count) * iterations / result.seconds : 0.0;
    auto host_likelihoods = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), likelihoods);
    auto host_scaled_sums = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scaled_sums);
    for (std::size_t pair = 0; pair < state.pair_count; ++pair) {
        result.likelihoods[pair] = host_likelihoods(pair);
        result.scaled_sums[pair] = host_scaled_sums(pair);
        result.checksum += result.likelihoods[pair];
    }
    return result;
}

}  // namespace

struct PersistentBucketPlan::Impl {
    using Shape = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;
    std::size_t max_cached_shapes = 8;
    std::size_t hits = 0;
    std::map<Shape, std::unique_ptr<PersistentState>> states;

    explicit Impl(std::size_t max_shapes) : max_cached_shapes(std::max<std::size_t>(1, max_shapes)) {}

    void clear() {
        states.clear();
        hits = 0;
    }
};

PersistentBucketPlan::PersistentBucketPlan(std::size_t max_cached_shapes)
    : impl_(std::make_unique<Impl>(max_cached_shapes)) {}

PersistentBucketPlan::~PersistentBucketPlan() = default;
PersistentBucketPlan::PersistentBucketPlan(PersistentBucketPlan&&) noexcept = default;
PersistentBucketPlan& PersistentBucketPlan::operator=(PersistentBucketPlan&&) noexcept = default;

KokkosBatchResult PersistentBucketPlan::execute(
    const std::vector<PairHmmRead>& reads,
    const std::vector<PairHmmHaplotype>& haplotypes,
    const std::vector<PairHmmRequest>& requests,
    int iterations) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");
    KokkosBatchResult aggregate;
    aggregate.simd_width = Simd::size();
    aggregate.execution_space = ExecSpace::name();
    aggregate.likelihoods.assign(requests.size(), 0.0);
    aggregate.scaled_sums.assign(requests.size(), 0.0);
    if (requests.empty()) return aggregate;
    using LengthKey = std::pair<std::size_t, std::size_t>;
    std::map<LengthKey, std::vector<std::size_t>> buckets;
    for (std::size_t request = 0; request < requests.size(); ++request) {
        const auto& item = requests[request];
        if (item.read_id >= reads.size() || item.haplotype_id >= haplotypes.size())
            throw std::invalid_argument("persistent PairHMM request index exceeds input");
        const auto& read = reads[item.read_id];
        const auto& haplotype = haplotypes[item.haplotype_id];
        if (read.bases.empty() || haplotype.bases.empty())
            throw std::invalid_argument("persistent PairHMM input contains an empty sequence");
        if (read.qualities.size() != read.bases.size() ||
            read.insertion_gop.size() != read.bases.size() ||
            read.deletion_gop.size() != read.bases.size() ||
            read.gap_continuation.size() != read.bases.size())
            throw std::invalid_argument("persistent PairHMM read quality arrays differ");
        buckets[{read.bases.size(), haplotype.bases.size()}].push_back(request);
    }

    for (const auto& [lengths, request_indices] : buckets) {
        const auto read_length = lengths.first;
        const auto haplotype_length = lengths.second;
        std::map<std::uint32_t, std::uint32_t> read_local;
        std::map<std::uint32_t, std::uint32_t> haplotype_local;
        std::vector<std::uint32_t> bucket_reads;
        std::vector<std::uint32_t> bucket_haplotypes;
        for (const auto request_index : request_indices) {
            const auto& request = requests[request_index];
            if (!read_local.contains(request.read_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_reads.size());
                read_local.emplace(request.read_id, local);
                bucket_reads.push_back(request.read_id);
            }
            if (!haplotype_local.contains(request.haplotype_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_haplotypes.size());
                haplotype_local.emplace(request.haplotype_id, local);
                bucket_haplotypes.push_back(request.haplotype_id);
            }
        }
        const auto record_count = bucket_reads.size() + bucket_haplotypes.size();
        const auto haplotype_base = bucket_reads.size();
        std::vector<PairInput> records;
        records.reserve(record_count);
        for (const auto read_id : bucket_reads) {
            const auto& source = reads[read_id];
            PairInput input;
            input.haplotype.assign(haplotype_length, static_cast<std::uint8_t>('N'));
            input.read = source.bases;
            input.read_qual = source.qualities;
            input.insertion_gop = source.insertion_gop;
            input.deletion_gop = source.deletion_gop;
            input.gap_continuation = source.gap_continuation;
            records.push_back(std::move(input));
        }
        for (const auto haplotype_id : bucket_haplotypes) {
            const auto& source = haplotypes[haplotype_id];
            PairInput input;
            input.haplotype = source.bases;
            input.read.assign(read_length, static_cast<std::uint8_t>('N'));
            input.read_qual.assign(read_length, 30);
            input.insertion_gop.assign(read_length, 40);
            input.deletion_gop.assign(read_length, 40);
            input.gap_continuation.assign(read_length, 10);
            records.push_back(std::move(input));
        }
        PairIndexBatch pairs;
        pairs.read_ids.reserve(request_indices.size());
        pairs.haplotype_ids.reserve(request_indices.size());
        for (const auto request_index : request_indices) {
            const auto& request = requests[request_index];
            pairs.read_ids.push_back(read_local.at(request.read_id));
            pairs.haplotype_ids.push_back(static_cast<std::uint32_t>(
                haplotype_base + haplotype_local.at(request.haplotype_id)));
        }
        const Impl::Shape shape{read_length, haplotype_length, record_count, request_indices.size()};
        auto found = impl_->states.find(shape);
        if (found == impl_->states.end()) {
            if (impl_->states.size() >= impl_->max_cached_shapes) impl_->states.erase(impl_->states.begin());
            found = impl_->states.emplace(shape, std::make_unique<PersistentState>(
                read_length, haplotype_length, record_count, request_indices.size())).first;
        } else {
            ++impl_->hits;
        }
        KokkosBatchResult bucket_result;
        prepare_persistent_state(*found->second, records, pairs, bucket_result);
        const auto execution = execute_persistent_state(*found->second, iterations);
        bucket_result.likelihoods = execution.likelihoods;
        bucket_result.scaled_sums = execution.scaled_sums;
        bucket_result.seconds = execution.seconds;
        bucket_result.pairs_per_second = execution.pairs_per_second;
        bucket_result.checksum = execution.checksum;
        bucket_result.simd_width = execution.simd_width;
        bucket_result.execution_space = execution.execution_space;
        aggregate.prepare_seconds += bucket_result.prepare_seconds;
        aggregate.seconds += bucket_result.seconds;
        for (std::size_t local = 0; local < request_indices.size(); ++local) {
            const auto original = request_indices[local];
            aggregate.likelihoods[original] = bucket_result.likelihoods[local];
            aggregate.scaled_sums[original] = bucket_result.scaled_sums[local];
            aggregate.checksum += bucket_result.likelihoods[local];
        }
    }
    aggregate.cache_hits = impl_->hits;
    aggregate.cached_shapes = impl_->states.size();
    aggregate.pairs_per_second = aggregate.seconds > 0.0
        ? static_cast<double>(requests.size()) * iterations / aggregate.seconds : 0.0;
    return aggregate;
}

void PersistentBucketPlan::clear() { impl_->clear(); }
std::size_t PersistentBucketPlan::cached_shapes() const { return impl_->states.size(); }
std::size_t PersistentBucketPlan::cache_hits() const { return impl_->hits; }

namespace {

constexpr std::size_t kFlowSize = 4;
constexpr std::size_t kFlowProbabilityStride = 256;

KOKKOS_INLINE_FUNCTION std::size_t flow_cell_index(
    std::size_t pair, std::size_t state, std::size_t row,
    std::size_t column, std::size_t stride) {
    return pair * (3 * 5 * stride) + state * (5 * stride) +
           (row % 5) * stride + column;
}

KokkosBatchResult compute_kokkos_flow_bucket(
    const std::vector<FlowPairHmmRead>& reads,
    const std::vector<FlowPairHmmHaplotype>& haplotypes,
    const std::vector<FlowPairHmmRequest>& requests,
    int iterations) {
    const auto read_length = reads.front().key.size();
    const auto haplotype_length = haplotypes.front().key.size();
    const auto pair_count = requests.size();
    const auto read_padded = read_length + 1 + kFlowSize;
    const auto haplotype_padded = haplotype_length + 1 + kFlowSize;
    const auto stride = haplotype_padded;
    if (pair_count > std::numeric_limits<std::size_t>::max() / (3 * 5 * stride))
        throw std::invalid_argument("flow PairHMM workspace dimensions overflow");
    const auto workspace_size = pair_count * 3 * 5 * stride;
    const auto& tables = gatk_tables();

    fastgatk::core::HostBatch host("pairhmm-flow-v1");
    host.records = pair_count;
    host.bytes = pair_count * (read_length * sizeof(std::int32_t) +
                               haplotype_length * sizeof(std::int32_t) +
                               read_length * kFlowProbabilityStride * sizeof(double));
    fastgatk::core::KernelPlan<ExecSpace> plan("pairhmm-flow");
    plan.begin_prepare(host);

    Kokkos::View<std::int32_t*> read_key("flow_read_key", pair_count * read_length);
    Kokkos::View<std::int32_t*> haplotype_key("flow_haplotype_key", pair_count * haplotype_length);
    Kokkos::View<std::uint8_t*> read_order("flow_read_order", pair_count * read_length);
    Kokkos::View<std::uint8_t*> haplotype_order("flow_haplotype_order", pair_count * haplotype_length);
    Kokkos::View<double*> probabilities("flow_probabilities",
                                        pair_count * read_length * kFlowProbabilityStride);
    Kokkos::View<double*> match_to_match("flow_match_to_match", pair_count * read_length);
    Kokkos::View<double*> match_to_insertion("flow_match_to_insertion", pair_count * read_length);
    Kokkos::View<double*> match_to_deletion("flow_match_to_deletion", pair_count * read_length);
    Kokkos::View<double*> indel_to_match("flow_indel_to_match", pair_count * read_length);
    Kokkos::View<double*> insertion_to_insertion("flow_insertion_to_insertion", pair_count * read_length);
    Kokkos::View<double*> deletion_to_deletion("flow_deletion_to_deletion", pair_count * read_length);
    Kokkos::View<std::uint32_t*> hap_start("flow_hap_start", pair_count);
    Kokkos::View<double*> workspace("flow_workspace", workspace_size);
    Kokkos::View<double*> likelihoods("flow_likelihoods", pair_count);
    Kokkos::View<double*> scaled_sums("flow_scaled_sums", pair_count);

    auto h_read_key = Kokkos::create_mirror_view(read_key);
    auto h_haplotype_key = Kokkos::create_mirror_view(haplotype_key);
    auto h_read_order = Kokkos::create_mirror_view(read_order);
    auto h_haplotype_order = Kokkos::create_mirror_view(haplotype_order);
    auto h_probabilities = Kokkos::create_mirror_view(probabilities);
    auto h_match_to_match = Kokkos::create_mirror_view(match_to_match);
    auto h_match_to_insertion = Kokkos::create_mirror_view(match_to_insertion);
    auto h_match_to_deletion = Kokkos::create_mirror_view(match_to_deletion);
    auto h_indel_to_match = Kokkos::create_mirror_view(indel_to_match);
    auto h_insertion_to_insertion = Kokkos::create_mirror_view(insertion_to_insertion);
    auto h_deletion_to_deletion = Kokkos::create_mirror_view(deletion_to_deletion);
    auto h_hap_start = Kokkos::create_mirror_view(hap_start);

    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        const auto& read = reads[requests[pair].read_id];
        const auto& haplotype = haplotypes[requests[pair].haplotype_id];
        for (std::size_t i = 0; i < read_length; ++i) {
            h_read_key(pair * read_length + i) = read.key[i];
            h_read_order(pair * read_length + i) = read.flow_order[i];
            const auto source_probability = i * kFlowProbabilityStride;
            const auto destination_probability = pair * read_length * kFlowProbabilityStride +
                                                  source_probability;
            for (std::size_t value = 0; value < kFlowProbabilityStride; ++value)
                h_probabilities(destination_probability + value) = read.probabilities[source_probability + value];
            const auto ins = read.insertion_gop[i];
            const auto del = read.deletion_gop[i];
            const auto gcp = read.gap_continuation[i];
            h_match_to_match(pair * read_length + i) =
                tables.match[static_cast<std::size_t>(ins) * 256 + del];
            h_match_to_insertion(pair * read_length + i) = tables.quality[ins];
            h_match_to_deletion(pair * read_length + i) = tables.quality[del];
            h_indel_to_match(pair * read_length + i) = 1.0 - tables.quality[gcp];
            h_insertion_to_insertion(pair * read_length + i) = tables.quality[gcp];
            h_deletion_to_deletion(pair * read_length + i) = tables.quality[gcp];
        }
        for (std::size_t j = 0; j < haplotype_length; ++j) {
            h_haplotype_key(pair * haplotype_length + j) = haplotype.key[j];
            h_haplotype_order(pair * haplotype_length + j) = haplotype.flow_order[j];
        }
        std::uint32_t start = 0;
        if (!read.flow_order.empty()) {
            for (std::size_t j = 0; j < haplotype_length; ++j) {
                if (haplotype.flow_order[j] == read.flow_order.front()) {
                    start = static_cast<std::uint32_t>(j);
                    break;
                }
            }
        }
        h_hap_start(pair) = start;
    }
    Kokkos::deep_copy(read_key, h_read_key);
    Kokkos::deep_copy(haplotype_key, h_haplotype_key);
    Kokkos::deep_copy(read_order, h_read_order);
    Kokkos::deep_copy(haplotype_order, h_haplotype_order);
    Kokkos::deep_copy(probabilities, h_probabilities);
    Kokkos::deep_copy(match_to_match, h_match_to_match);
    Kokkos::deep_copy(match_to_insertion, h_match_to_insertion);
    Kokkos::deep_copy(match_to_deletion, h_match_to_deletion);
    Kokkos::deep_copy(indel_to_match, h_indel_to_match);
    Kokkos::deep_copy(insertion_to_insertion, h_insertion_to_insertion);
    Kokkos::deep_copy(deletion_to_deletion, h_deletion_to_deletion);
    Kokkos::deep_copy(hap_start, h_hap_start);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(pair_count);
    device.bind("flow_read_key", read_key);
    device.bind("flow_haplotype_key", haplotype_key);
    device.bind("flow_read_order", read_order);
    device.bind("flow_haplotype_order", haplotype_order);
    device.bind("flow_probabilities", probabilities);
    device.bind("flow_match_to_match", match_to_match);
    device.bind("flow_match_to_insertion", match_to_insertion);
    device.bind("flow_match_to_deletion", match_to_deletion);
    device.bind("flow_indel_to_match", indel_to_match);
    device.bind("flow_insertion_to_insertion", insertion_to_insertion);
    device.bind("flow_deletion_to_deletion", deletion_to_deletion);
    device.bind("flow_hap_start", hap_start);
    device.bind("flow_workspace", workspace);
    device.bind("flow_likelihoods", likelihoods);
    device.bind("flow_scaled_sums", scaled_sums);
    plan.end_prepare(device);

    const double initial = std::pow(2.0, 1020.0) / static_cast<double>(haplotype_length);
    const double initial_log10 = gatk_strict_log10(std::pow(2.0, 1020.0));
    // One league entry owns one independent flow-pair workspace.  TeamPolicy
    // is deliberately used here even while the recurrence itself remains
    // serial in (read-flow, haplotype-flow) order: it gives CUDA/HIP/SYCL a
    // backend-native launch and a stable place to add anti-diagonal/team
    // tiling later, without introducing a second device-only algorithm.  A
    // single team member preserves the Strict recurrence order and therefore
    // the CPU oracle's floating-point bits.
    const auto launch = [&] {
        using FlowTeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
        using FlowTeamMember = FlowTeamPolicy::member_type;
        Kokkos::parallel_for("fastgatk_pairhmm_flow_team",
            FlowTeamPolicy(pair_count, 1),
            KOKKOS_LAMBDA(const FlowTeamMember& team) {
                if (team.team_rank() != 0) return;
                const std::size_t pair = team.league_rank();
                const auto stride_local = haplotype_padded;
                for (std::size_t state = 0; state < 3; ++state)
                    for (std::size_t row = 0; row < 5; ++row)
                        for (std::size_t column = 0; column < stride_local; ++column)
                            workspace(flow_cell_index(pair, state, row, column, stride_local)) = 0.0;
                for (std::size_t row = 0; row < 5; ++row)
                    for (std::size_t column = 0; column < stride_local; ++column)
                        workspace(flow_cell_index(pair, 2, row, column, stride_local)) = initial;

                for (std::size_t i = kFlowSize + 1; i < read_padded; ++i) {
                    const auto row = i % 5;
                    for (std::size_t state = 0; state < 3; ++state)
                        for (std::size_t column = 0; column < stride_local; ++column)
                            workspace(flow_cell_index(pair, state, row, column, stride_local)) = 0.0;
                    const auto read_flow = i - (kFlowSize + 1);
                    const auto transition_offset = pair * read_length + read_flow;
                    const auto phase = (hap_start(pair) + i) % kFlowSize;
                    for (std::size_t j = phase + kFlowSize; j < haplotype_padded; j += kFlowSize) {
                        double prior = 0.0;
                        if (j >= kFlowSize + 1) {
                            const auto hap_flow = j - (kFlowSize + 1);
                            if (hap_flow < haplotype_length &&
                                read_order(pair * read_length + read_flow) ==
                                    haplotype_order(pair * haplotype_length + hap_flow)) {
                                const auto probability_offset =
                                    pair * read_length * kFlowProbabilityStride +
                                    read_flow * kFlowProbabilityStride;
                                const auto hap_value = haplotype_key(pair * haplotype_length + hap_flow);
                                const auto probability_index = hap_value < 0 ? 0U :
                                    (hap_value >= static_cast<std::int32_t>(kFlowProbabilityStride)
                                         ? static_cast<std::uint32_t>(kFlowProbabilityStride - 1)
                                         : static_cast<std::uint32_t>(hap_value));
                                prior = probabilities(probability_offset + probability_index);
                            }
                        }
                        const auto previous_row = (i - 1) % 5;
                        const auto insertion_row = (i - kFlowSize) % 5;
                        const auto current_match = flow_cell_index(pair, 0, row, j, stride_local);
                        const auto current_insertion = flow_cell_index(pair, 1, row, j, stride_local);
                        const auto current_deletion = flow_cell_index(pair, 2, row, j, stride_local);
                        const auto previous_match = flow_cell_index(pair, 0, previous_row, j - 1, stride_local);
                        const auto previous_insertion = flow_cell_index(pair, 1, previous_row, j - 1, stride_local);
                        const auto previous_deletion = flow_cell_index(pair, 2, previous_row, j - 1, stride_local);
                        const auto insertion_match = flow_cell_index(pair, 0, insertion_row, j, stride_local);
                        const auto insertion_insertion = flow_cell_index(pair, 1, insertion_row, j, stride_local);
                        const auto deletion_match = flow_cell_index(pair, 0, row, j - kFlowSize, stride_local);
                        const auto deletion_deletion = flow_cell_index(pair, 2, row, j - kFlowSize, stride_local);
                        workspace(current_match) = prior *
                            (workspace(previous_match) * match_to_match(transition_offset) +
                             workspace(previous_insertion) * indel_to_match(transition_offset) +
                             workspace(previous_deletion) * indel_to_match(transition_offset));
                        workspace(current_insertion) =
                            workspace(insertion_match) * match_to_insertion(transition_offset) +
                            workspace(insertion_insertion) * insertion_to_insertion(transition_offset);
                        workspace(current_deletion) =
                            workspace(deletion_match) * match_to_deletion(transition_offset) +
                            workspace(deletion_deletion) * deletion_to_deletion(transition_offset);
                    }
                }
                const auto final_row = (read_padded - 1) % 5;
                double sum = 0.0;
                for (std::size_t j = 1; j < haplotype_padded; ++j)
                    sum += workspace(flow_cell_index(pair, 0, final_row, j, stride_local)) +
                           workspace(flow_cell_index(pair, 1, final_row, j, stride_local));
                scaled_sums(pair) = sum;
                likelihoods(pair) = gatk_strict_log10(sum) - initial_log10;
            });
    };
    plan.begin_execute();
    for (int iteration = 0; iteration < iterations; ++iteration) launch();
    ExecSpace().fence();
    plan.end_execute();

    KokkosBatchResult result;
    result.error_model = "flow";
    result.simd_width = 1;
    result.execution_space = ExecSpace::name();
    result.execution_policy = "TeamPolicy";
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    result.seconds = plan.telemetry().execute_seconds;
    result.likelihoods.resize(pair_count);
    result.scaled_sums.resize(pair_count);
    auto host_likelihoods = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), likelihoods);
    auto host_scaled_sums = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scaled_sums);
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        result.likelihoods[pair] = host_likelihoods(pair);
        result.scaled_sums[pair] = host_scaled_sums(pair);
        result.checksum += result.likelihoods[pair];
    }
    result.pairs_per_second = result.seconds > 0.0
        ? static_cast<double>(pair_count) * iterations / result.seconds : 0.0;
    return result;
}

}  // namespace

KokkosBatchResult compute_kokkos_flow(
    const std::vector<FlowPairHmmRead>& reads,
    const std::vector<FlowPairHmmHaplotype>& haplotypes,
    const std::vector<FlowPairHmmRequest>& requests,
    int iterations) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");
    KokkosBatchResult aggregate;
    aggregate.error_model = "flow";
    aggregate.simd_width = 1;
    aggregate.execution_space = ExecSpace::name();
    aggregate.execution_policy = "TeamPolicy";
    aggregate.likelihoods.assign(requests.size(), 0.0);
    aggregate.scaled_sums.assign(requests.size(), 0.0);
    if (requests.empty()) return aggregate;
    using LengthKey = std::pair<std::size_t, std::size_t>;
    std::map<LengthKey, std::vector<std::size_t>> buckets;
    for (std::size_t index = 0; index < requests.size(); ++index) {
        const auto& request = requests[index];
        if (request.read_id >= reads.size() || request.haplotype_id >= haplotypes.size())
            throw std::invalid_argument("flow PairHMM request index exceeds input");
        const auto& read = reads[request.read_id];
        const auto& haplotype = haplotypes[request.haplotype_id];
        if (read.key.empty() || haplotype.key.empty())
            throw std::invalid_argument("flow PairHMM input contains an empty key");
        if (read.flow_order.size() != read.key.size() ||
            read.insertion_gop.size() != read.key.size() ||
            read.deletion_gop.size() != read.key.size() ||
            read.gap_continuation.size() != read.key.size() ||
            read.key.size() > std::numeric_limits<std::size_t>::max() /
                kFlowProbabilityStride ||
            read.probabilities.size() != read.key.size() * kFlowProbabilityStride)
            throw std::invalid_argument("flow PairHMM read arrays/probability table differ");
        if (haplotype.flow_order.size() != haplotype.key.size())
            throw std::invalid_argument("flow PairHMM haplotype key/order arrays differ");
        for (const auto probability : read.probabilities)
            if (!std::isfinite(probability) || probability < 0.0 || probability > 1.0)
                throw std::invalid_argument("flow PairHMM probability is outside [0,1]");
        buckets[{read.key.size(), haplotype.key.size()}].push_back(index);
    }
    for (const auto& [lengths, indices] : buckets) {
        std::vector<FlowPairHmmRead> bucket_reads;
        std::vector<FlowPairHmmHaplotype> bucket_haplotypes;
        std::vector<FlowPairHmmRequest> bucket_requests;
        std::map<std::uint32_t, std::uint32_t> read_local;
        std::map<std::uint32_t, std::uint32_t> haplotype_local;
        for (const auto index : indices) {
            const auto& request = requests[index];
            if (!read_local.contains(request.read_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_reads.size());
                read_local.emplace(request.read_id, local);
                bucket_reads.push_back(reads[request.read_id]);
            }
            if (!haplotype_local.contains(request.haplotype_id)) {
                const auto local = static_cast<std::uint32_t>(bucket_haplotypes.size());
                haplotype_local.emplace(request.haplotype_id, local);
                bucket_haplotypes.push_back(haplotypes[request.haplotype_id]);
            }
            bucket_requests.push_back(FlowPairHmmRequest{
                read_local.at(request.read_id), haplotype_local.at(request.haplotype_id)});
        }
        const auto bucket_result = compute_kokkos_flow_bucket(
            bucket_reads, bucket_haplotypes, bucket_requests, iterations);
        aggregate.prepare_seconds += bucket_result.prepare_seconds;
        aggregate.seconds += bucket_result.seconds;
        for (std::size_t local = 0; local < indices.size(); ++local) {
            aggregate.likelihoods[indices[local]] = bucket_result.likelihoods[local];
            aggregate.scaled_sums[indices[local]] = bucket_result.scaled_sums[local];
            aggregate.checksum += bucket_result.likelihoods[local];
        }
    }
    aggregate.pairs_per_second = aggregate.seconds > 0.0
        ? static_cast<double>(requests.size()) * iterations / aggregate.seconds : 0.0;
    return aggregate;
}

KokkosBatchResult compute_kokkos_flow_alignment(
    const std::vector<FlowPairHmmRead>& reads,
    const std::vector<FlowPairHmmHaplotype>& haplotypes,
    const std::vector<FlowPairHmmRequest>& requests,
    int iterations) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (iterations < 1) throw std::invalid_argument("iterations must be positive");

    KokkosBatchResult result;
    result.error_model = "flow-alignment";
    result.simd_width = 1;
    result.execution_space = ExecSpace::name();
    result.execution_policy = "RangePolicy";
    result.likelihoods.assign(requests.size(), -std::numeric_limits<double>::infinity());
    result.scaled_sums.assign(requests.size(), -std::numeric_limits<double>::infinity());
    if (requests.empty()) return result;

    std::size_t read_flow_count = 0;
    std::size_t haplotype_flow_count = 0;
    for (const auto& read : reads) {
        if (read.key.empty() || read.flow_order.size() != read.key.size() ||
            read.key.size() > std::numeric_limits<std::size_t>::max() / kFlowProbabilityStride ||
            read.probabilities.size() != read.key.size() * kFlowProbabilityStride)
            throw std::invalid_argument("flow alignment read arrays/probability table differ");
        for (const auto probability : read.probabilities)
            if (!std::isfinite(probability) || probability < 0.0)
                throw std::invalid_argument("flow alignment probability is invalid");
        if (read.key.size() > std::numeric_limits<std::size_t>::max() - read_flow_count)
            throw std::invalid_argument("flow alignment read offsets overflow");
        read_flow_count += read.key.size();
    }
    for (const auto& haplotype : haplotypes) {
        if (haplotype.key.empty() || haplotype.flow_order.size() != haplotype.key.size())
            throw std::invalid_argument("flow alignment haplotype key/order arrays differ");
        if (haplotype.key.size() > std::numeric_limits<std::size_t>::max() - haplotype_flow_count)
            throw std::invalid_argument("flow alignment haplotype offsets overflow");
        haplotype_flow_count += haplotype.key.size();
    }
    for (const auto& request : requests)
        if (request.read_id >= reads.size() || request.haplotype_id >= haplotypes.size())
            throw std::invalid_argument("flow alignment request index exceeds input");

    fastgatk::core::HostBatch host("flow-alignment-v1");
    host.records = requests.size();
    host.bytes = read_flow_count * (sizeof(std::int32_t) + sizeof(std::uint8_t) +
                                     kFlowProbabilityStride * sizeof(double)) +
                 haplotype_flow_count * (sizeof(std::int32_t) + sizeof(std::uint8_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("flow-alignment");
    plan.begin_prepare(host);

    Kokkos::View<std::int32_t*> haplotype_keys("flow_alignment_haplotype_keys",
                                                haplotype_flow_count);
    Kokkos::View<std::uint8_t*> read_orders("flow_alignment_read_orders", read_flow_count);
    Kokkos::View<std::uint8_t*> haplotype_orders("flow_alignment_haplotype_orders",
                                                  haplotype_flow_count);
    Kokkos::View<double*> probabilities("flow_alignment_probabilities",
                                        read_flow_count * kFlowProbabilityStride);
    Kokkos::View<std::size_t*> read_offsets("flow_alignment_read_offsets", reads.size() + 1);
    Kokkos::View<std::size_t*> haplotype_offsets("flow_alignment_haplotype_offsets",
                                                  haplotypes.size() + 1);
    Kokkos::View<std::uint32_t*> request_reads("flow_alignment_request_reads", requests.size());
    Kokkos::View<std::uint32_t*> request_haplotypes("flow_alignment_request_haplotypes",
                                                    requests.size());
    Kokkos::View<double*> likelihoods("flow_alignment_likelihoods", requests.size());

    auto h_haplotype_keys = Kokkos::create_mirror_view(haplotype_keys);
    auto h_read_orders = Kokkos::create_mirror_view(read_orders);
    auto h_haplotype_orders = Kokkos::create_mirror_view(haplotype_orders);
    auto h_probabilities = Kokkos::create_mirror_view(probabilities);
    auto h_read_offsets = Kokkos::create_mirror_view(read_offsets);
    auto h_haplotype_offsets = Kokkos::create_mirror_view(haplotype_offsets);
    auto h_request_reads = Kokkos::create_mirror_view(request_reads);
    auto h_request_haplotypes = Kokkos::create_mirror_view(request_haplotypes);

    std::size_t read_offset = 0;
    h_read_offsets(0) = 0;
    for (std::size_t read_id = 0; read_id < reads.size(); ++read_id) {
        const auto& read = reads[read_id];
        for (std::size_t flow = 0; flow < read.key.size(); ++flow) {
            h_read_orders(read_offset + flow) = read.flow_order[flow];
            const auto source = flow * kFlowProbabilityStride;
            const auto destination = (read_offset + flow) * kFlowProbabilityStride;
            for (std::size_t hmer = 0; hmer < kFlowProbabilityStride; ++hmer)
                h_probabilities(destination + hmer) = read.probabilities[source + hmer];
        }
        read_offset += read.key.size();
        h_read_offsets(read_id + 1) = read_offset;
    }
    std::size_t haplotype_offset = 0;
    h_haplotype_offsets(0) = 0;
    for (std::size_t haplotype_id = 0; haplotype_id < haplotypes.size(); ++haplotype_id) {
        const auto& haplotype = haplotypes[haplotype_id];
        for (std::size_t flow = 0; flow < haplotype.key.size(); ++flow) {
            h_haplotype_keys(haplotype_offset + flow) = haplotype.key[flow];
            h_haplotype_orders(haplotype_offset + flow) = haplotype.flow_order[flow];
        }
        haplotype_offset += haplotype.key.size();
        h_haplotype_offsets(haplotype_id + 1) = haplotype_offset;
    }
    for (std::size_t request_id = 0; request_id < requests.size(); ++request_id) {
        h_request_reads(request_id) = requests[request_id].read_id;
        h_request_haplotypes(request_id) = requests[request_id].haplotype_id;
    }
    Kokkos::deep_copy(haplotype_keys, h_haplotype_keys);
    Kokkos::deep_copy(read_orders, h_read_orders);
    Kokkos::deep_copy(haplotype_orders, h_haplotype_orders);
    Kokkos::deep_copy(probabilities, h_probabilities);
    Kokkos::deep_copy(read_offsets, h_read_offsets);
    Kokkos::deep_copy(haplotype_offsets, h_haplotype_offsets);
    Kokkos::deep_copy(request_reads, h_request_reads);
    Kokkos::deep_copy(request_haplotypes, h_request_haplotypes);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(requests.size());
    device.bind("flow_alignment_haplotype_keys", haplotype_keys);
    device.bind("flow_alignment_read_orders", read_orders);
    device.bind("flow_alignment_haplotype_orders", haplotype_orders);
    device.bind("flow_alignment_probabilities", probabilities);
    device.bind("flow_alignment_read_offsets", read_offsets);
    device.bind("flow_alignment_haplotype_offsets", haplotype_offsets);
    device.bind("flow_alignment_request_reads", request_reads);
    device.bind("flow_alignment_request_haplotypes", request_haplotypes);
    device.bind("flow_alignment_likelihoods", likelihoods);
    plan.end_prepare(device);

    const auto launch = [&] {
        Kokkos::parallel_for("fastgatk_flow_alignment",
            Kokkos::RangePolicy<ExecSpace>(0, requests.size()),
            KOKKOS_LAMBDA(const std::size_t request_id) {
                const auto read_id = request_reads(request_id);
                const auto haplotype_id = request_haplotypes(request_id);
                const auto read_begin = read_offsets(read_id);
                const auto read_end = read_offsets(read_id + 1);
                const auto haplotype_begin = haplotype_offsets(haplotype_id);
                const auto haplotype_end = haplotype_offsets(haplotype_id + 1);
                const auto read_length = read_end - read_begin;
                const auto haplotype_length = haplotype_end - haplotype_begin;
                double best = -std::numeric_limits<double>::infinity();
                std::size_t starting_point = 0;
                const auto first_read_flow = read_orders(read_begin);
                for (std::size_t flow = 0; flow < haplotype_length; ++flow) {
                    if (haplotype_orders(haplotype_begin + flow) == first_read_flow) {
                        starting_point = flow;
                        break;
                    }
                }
                // FlowBasedAlignmentLikelihoodEngine.flowLikelihoodScore():
                // align only cycle-compatible starts, collect the log10
                // product for every full read placement, and retain its max.
                for (std::size_t start = starting_point;
                     start + read_length <= haplotype_length;
                     start += kFlowSize) {
                    double score = 0.0;
                    for (std::size_t flow = 0; flow < read_length; ++flow) {
                        const auto hmer = haplotype_keys(haplotype_begin + start + flow);
                        const auto probability_index = hmer < 0 ? 0U :
                            (hmer >= static_cast<std::int32_t>(kFlowProbabilityStride)
                                ? static_cast<std::uint32_t>(kFlowProbabilityStride - 1)
                                : static_cast<std::uint32_t>(hmer));
                        auto probability = probabilities(
                            (read_begin + flow) * kFlowProbabilityStride + probability_index);
                        // FlowBasedRead.getProb() limits an accumulated cell
                        // to one before Math.log10 observes it.
                        if (probability > 1.0) probability = 1.0;
                        score += gatk_strict_log10(probability);
                    }
                    if (score > best) best = score;
                }
                likelihoods(request_id) = best;
            });
    };
    plan.begin_execute();
    for (int iteration = 0; iteration < iterations; ++iteration) launch();
    ExecSpace().fence();
    plan.end_execute();

    result.prepare_seconds = plan.telemetry().prepare_seconds;
    result.seconds = plan.telemetry().execute_seconds;
    result.pairs_per_second = result.seconds > 0.0
        ? static_cast<double>(requests.size()) * iterations / result.seconds : 0.0;
    auto h_likelihoods = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), likelihoods);
    for (std::size_t request_id = 0; request_id < requests.size(); ++request_id) {
        result.likelihoods[request_id] = h_likelihoods(request_id);
        // This engine has no HMM scaled-sum state.  Retain its final score as
        // the diagnostic companion value so callers never mistake it for a
        // base/flow PairHMM underflow signal.
        result.scaled_sums[request_id] = h_likelihoods(request_id);
        result.checksum += h_likelihoods(request_id);
    }
    return result;
}

std::string kokkos_backend_description() {
    return std::string(ExecSpace::name()) + ", simd_width=" + std::to_string(Simd::size());
}

}  // namespace fastgatk::pairhmm
