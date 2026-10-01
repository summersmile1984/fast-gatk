#include "fastgatk/kernels/gpu_safety.hpp"
#include "fastgatk/kernels/smith_waterman.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;
constexpr int kNegativeInfinity = -1000000000;
using SimdInt = Kokkos::Experimental::simd<int>;

// Score a homogeneous batch with Kokkos' portable SIMD API.  SIMD lanes are
// independent alignments; the affine DP recurrence remains scalar in the
// read/reference dimensions, so this does not change GATK's tie/overhang
// policy.  Heterogeneous batches use the regular Kokkos kernel below.
SmithWatermanBatchResult smith_waterman_score_kokkos_simd_uniform(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy overhang) {
    constexpr std::size_t width = SimdInt::size();
    SmithWatermanBatchResult result;
    result.execution_space = ExecSpace::name();
    result.simd_width = width;
    result.simd_groups = (requests.size() + width - 1) / width;
    result.scores.assign(requests.size(), 0);
    if (requests.empty()) return result;

    const std::size_t read_length = requests.front().read.size();
    const std::size_t reference_length = requests.front().reference.size();
    if (read_length == 0 || reference_length == 0)
        throw std::invalid_argument("Smith-Waterman request must be non-empty");
    for (const auto& request : requests) {
        if (request.read.size() != read_length ||
            request.reference.size() != reference_length)
            throw std::invalid_argument("SIMD SW helper requires uniform lengths");
    }

    fastgatk::core::HostBatch host("smith-waterman-simd-v1");
    host.records = requests.size();
    host.bytes = requests.size() * (read_length + reference_length) * sizeof(int);
    fastgatk::core::KernelPlan<ExecSpace> plan("smith-waterman-simd");
    plan.begin_prepare(host);

    // The request dimension is the contiguous LayoutRight dimension, allowing
    // the SIMD constructor to load each lane with one portable Kokkos load.
    using BaseView = Kokkos::View<int**, Kokkos::LayoutRight, MemorySpace>;
    using ScoreView = Kokkos::View<int*, MemorySpace>;
    using SimdMatrix = Kokkos::View<SimdInt***, Kokkos::LayoutRight, MemorySpace>;
    const std::size_t groups = result.simd_groups;
    // SIMD loads are contiguous across the request (lane) dimension.  The
    // final group may contain fewer than ``width`` real requests; allocate a
    // padded lane region so the portable simd constructor never reads past
    // the View extent.  Padded lanes are inactive and their scores are never
    // copied back to the public result.
    const std::size_t padded_requests = groups * width;
    BaseView reads("sw_simd_reads", read_length, padded_requests);
    BaseView references("sw_simd_references", reference_length, padded_requests);
    ScoreView scores("sw_simd_scores", requests.size());
    const std::size_t stride = read_length + 1;
    SimdMatrix h_rows("sw_simd_h_rows", groups, 2, stride);
    SimdMatrix vertical_rows("sw_simd_vertical_rows", groups, 2, stride);
    SimdMatrix best_vertical_rows("sw_simd_best_vertical_rows", groups, 1, stride);
    auto host_reads = Kokkos::create_mirror_view(reads);
    auto host_references = Kokkos::create_mirror_view(references);
    for (std::size_t request = 0; request < requests.size(); ++request) {
        for (std::size_t index = 0; index < read_length; ++index)
            host_reads(index, request) = requests[request].read[index];
        for (std::size_t index = 0; index < reference_length; ++index)
            host_references(index, request) = requests[request].reference[index];
    }
    for (std::size_t request = requests.size(); request < padded_requests; ++request) {
        for (std::size_t index = 0; index < read_length; ++index)
            host_reads(index, request) = 0;
        for (std::size_t index = 0; index < reference_length; ++index)
            host_references(index, request) = 0;
    }
    Kokkos::deep_copy(reads, host_reads);
    Kokkos::deep_copy(references, host_references);
    ExecSpace().fence();

    fastgatk::core::DeviceBatch<ExecSpace> device(requests.size());
    device.bind("reads", reads);
    device.bind("references", references);
    device.bind("h_rows", h_rows);
    device.bind("vertical_rows", vertical_rows);
    device.bind("best_vertical_rows", best_vertical_rows);
    device.bind("scores", scores);
    plan.end_prepare(device);

    plan.begin_execute();
    const int match = parameters.match;
    const int mismatch = parameters.mismatch;
    const int gap_open = parameters.gap_open;
    const int gap_extend = parameters.gap_extend;
    const int strategy = static_cast<int>(overhang);
    Kokkos::parallel_for(
        "fastgatk_smith_waterman_simd",
        Kokkos::RangePolicy<ExecSpace>(0, groups),
        KOKKOS_LAMBDA(const std::size_t group) {
            const std::size_t base = group * width;
            const std::size_t active =
                base + width <= requests.size() ? width : requests.size() - base;
            const int negative = kNegativeInfinity;
            const bool gap_edges =
                strategy == static_cast<int>(SmithWatermanOverhangStrategy::Indel) ||
                strategy == static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel);
            const SimdInt zero(0);
            const SimdInt negative_v(negative);
            for (std::size_t j = 0; j <= read_length; ++j) {
                h_rows(group, 0, j) = zero;
                h_rows(group, 1, j) = zero;
                vertical_rows(group, 0, j) = negative_v;
                vertical_rows(group, 1, j) = negative_v;
                best_vertical_rows(group, 0, j) = negative_v;
            }
            if (gap_edges) {
                if (read_length >= 1) h_rows(group, 0, 1) = SimdInt(gap_open);
                for (std::size_t j = 2; j <= read_length; ++j)
                    h_rows(group, 0, j) = h_rows(group, 0, j - 1) + gap_extend;
            }
            int right_score[width];
            std::size_t right_row[width];
            int bottom_score[width];
            std::size_t bottom_column[width];
            int bottom_right_score[width];
            for (std::size_t lane = 0; lane < width; ++lane) {
                right_score[lane] = negative;
                right_row[lane] = 0;
                bottom_score[lane] = negative;
                bottom_column[lane] = 0;
                bottom_right_score[lane] = 0;
            }
            for (std::size_t i = 1; i <= reference_length; ++i) {
                const int current = static_cast<int>(i & 1U);
                const int previous = current ^ 1;
                h_rows(group, current, 0) = zero;
                vertical_rows(group, current, 0) = negative_v;
                if (gap_edges)
                    h_rows(group, current, 0) = i == 1
                        ? SimdInt(gap_open)
                        : h_rows(group, previous, 0) + gap_extend;
                SimdInt best_horizontal = negative_v;
                for (std::size_t j = 1; j <= read_length; ++j) {
                    const SimdInt read_v(&reads(j - 1, base),
                                         Kokkos::Experimental::simd_flag_default);
                    const SimdInt reference_v(&references(i - 1, base),
                                              Kokkos::Experimental::simd_flag_default);
                    const auto substitution = Kokkos::Experimental::condition(
                        read_v == reference_v, SimdInt(match), SimdInt(mismatch));
                    const SimdInt diagonal =
                        h_rows(group, previous, j - 1) + substitution;
                    const SimdInt previous_vertical =
                        h_rows(group, previous, j) + gap_open;
                    best_vertical_rows(group, 0, j) += gap_extend;
                    best_vertical_rows(group, 0, j) = Kokkos::max(
                        previous_vertical, best_vertical_rows(group, 0, j));
                    const SimdInt vertical = best_vertical_rows(group, 0, j);
                    const SimdInt previous_horizontal =
                        h_rows(group, current, j - 1) + gap_open;
                    best_horizontal += gap_extend;
                    best_horizontal = Kokkos::max(previous_horizontal, best_horizontal);
                    const SimdInt right = best_horizontal;
                    const SimdInt cell = Kokkos::max(
                        negative_v, Kokkos::max(diagonal, Kokkos::max(vertical, right)));
                    h_rows(group, current, j) = cell;
                    vertical_rows(group, current, j) = vertical;
                    for (std::size_t lane = 0; lane < active; ++lane) {
                        const int value = cell[lane];
                        if (j == read_length) {
                            if (value >= right_score[lane]) {
                                right_score[lane] = value;
                                right_row[lane] = i;
                            }
                            if (i == reference_length) {
                                bottom_right_score[lane] = value;
                                if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                                    const auto distance = reference_length > j
                                        ? reference_length - j : j - reference_length;
                                    const auto old_distance = reference_length > bottom_column[lane]
                                        ? reference_length - bottom_column[lane]
                                        : bottom_column[lane] - reference_length;
                                    if (value > bottom_score[lane] ||
                                        (value == bottom_score[lane] && distance < old_distance)) {
                                        bottom_score[lane] = value;
                                        bottom_column[lane] = j;
                                    }
                                }
                            }
                        } else if (i == reference_length &&
                                   strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                            const auto distance = reference_length > j
                                ? reference_length - j : j - reference_length;
                            const auto old_distance = reference_length > bottom_column[lane]
                                ? reference_length - bottom_column[lane]
                                : bottom_column[lane] - reference_length;
                            if (value > bottom_score[lane] ||
                                (value == bottom_score[lane] && distance < old_distance)) {
                                bottom_score[lane] = value;
                                bottom_column[lane] = j;
                            }
                        }
                    }
                }
            }
            for (std::size_t lane = 0; lane < active; ++lane) {
                int selected = bottom_right_score[lane];
                if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::Indel)) {
                    selected = right_score[lane];
                    const auto bottom_distance = reference_length > bottom_column[lane]
                        ? reference_length - bottom_column[lane]
                        : bottom_column[lane] - reference_length;
                    const auto right_distance = right_row[lane] > read_length
                        ? right_row[lane] - read_length : read_length - right_row[lane];
                    if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel) &&
                        (bottom_score[lane] > selected ||
                         (bottom_score[lane] == selected && bottom_distance < right_distance)))
                        selected = bottom_score[lane];
                }
                scores(base + lane) = selected;
            }
        });
    ExecSpace().fence();
    plan.end_execute();
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    result.seconds = plan.telemetry().execute_seconds;
    auto host_scores = Kokkos::create_mirror_view(scores);
    Kokkos::deep_copy(host_scores, scores);
    for (std::size_t index = 0; index < requests.size(); ++index)
        result.scores[index] = host_scores(index);
    return result;
}

}  // namespace

SmithWatermanBatchResult smith_waterman_score_kokkos(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy overhang) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    SmithWatermanBatchResult result;
    result.execution_space = ExecSpace::name();
    result.simd_width = SimdInt::size();
    result.scores.assign(requests.size(), 0);
    if (requests.empty()) return result;
    // Bound the ragged DP workspace.  A single AssemblyRegion can produce
    // thousands of read/haplotype requests; allocating all rows at once makes
    // Host-space deallocation dominate wall time on OpenMP and creates a large
    // transient footprint on accelerators.  Chunking preserves request order
    // and the exact recurrence while keeping one Kokkos workspace below the
    // allocator's large-object threshold.
    constexpr std::size_t request_chunk = 512;
    if (requests.size() > request_chunk) {
        result.scores.resize(requests.size());
        for (std::size_t begin = 0; begin < requests.size(); begin += request_chunk) {
            const auto end = std::min(requests.size(), begin + request_chunk);
            std::vector<SmithWatermanRequest> chunk(
                requests.begin() + static_cast<std::ptrdiff_t>(begin),
                requests.begin() + static_cast<std::ptrdiff_t>(end));
            const auto partial = smith_waterman_score_kokkos(chunk, parameters, overhang);
            std::copy(partial.scores.begin(), partial.scores.end(),
                      result.scores.begin() + static_cast<std::ptrdiff_t>(begin));
            result.prepare_seconds += partial.prepare_seconds;
            result.seconds += partial.seconds;
            result.simd_width = std::max(result.simd_width, partial.simd_width);
            result.simd_groups += partial.simd_groups;
        }
        return result;
    }

    // Use the SIMD API when the batch is rectangular and large enough to fill
    // at least one vector.  The normal path below remains the deterministic
    // fallback for ragged request sets and for scalar-only Kokkos backends.
    if constexpr (SimdInt::size() > 1) {
        if (requests.size() >= SimdInt::size()) {
            const auto read_length = requests.front().read.size();
            const auto reference_length = requests.front().reference.size();
            const bool uniform = std::all_of(requests.begin(), requests.end(),
                [&](const auto& request) {
                    return !request.read.empty() && !request.reference.empty() &&
                           request.read.size() == read_length &&
                           request.reference.size() == reference_length;
                });
            if (uniform)
                return smith_waterman_score_kokkos_simd_uniform(requests, parameters, overhang);
        }
    }

    std::size_t max_read = 0, max_reference = 0;
    std::size_t host_bytes = 0;
    std::vector<std::uint32_t> read_offsets(requests.size() + 1, 0);
    std::vector<std::uint32_t> reference_offsets(requests.size() + 1, 0);
    std::vector<std::uint8_t> reads, references;
    for (std::size_t i = 0; i < requests.size(); ++i) {
        const auto& request = requests[i];
        if (request.read.empty() || request.reference.empty())
            throw std::invalid_argument("Smith-Waterman Kokkos requests must be non-empty");
        if (request.read.size() > UINT32_MAX || request.reference.size() > UINT32_MAX)
            throw std::invalid_argument("Smith-Waterman request exceeds uint32 offsets");
        max_read = std::max(max_read, request.read.size());
        max_reference = std::max(max_reference, request.reference.size());
        read_offsets[i + 1] = read_offsets[i] + static_cast<std::uint32_t>(request.read.size());
        reference_offsets[i + 1] = reference_offsets[i] + static_cast<std::uint32_t>(request.reference.size());
        reads.insert(reads.end(), request.read.begin(), request.read.end());
        references.insert(references.end(), request.reference.begin(), request.reference.end());
    }
    host_bytes = reads.size() + references.size() +
                 (read_offsets.size() + reference_offsets.size()) * sizeof(std::uint32_t);
    fastgatk::core::HostBatch host("smith-waterman-v1");
    host.records = requests.size();
    host.bytes = host_bytes;
    fastgatk::core::KernelPlan<ExecSpace> plan("smith-waterman");
    plan.begin_prepare(host);

    Kokkos::View<std::uint8_t*, MemorySpace> read_bases("sw_read_bases", reads.size());
    Kokkos::View<std::uint8_t*, MemorySpace> reference_bases("sw_reference_bases", references.size());
    Kokkos::View<std::uint32_t*, MemorySpace> read_index("sw_read_offsets", read_offsets.size());
    Kokkos::View<std::uint32_t*, MemorySpace> reference_index("sw_reference_offsets", reference_offsets.size());
    Kokkos::View<int*, MemorySpace> scores("sw_scores", requests.size());
    using Matrix = Kokkos::View<int***, Kokkos::LayoutRight, MemorySpace>;
    // Two H rows and two affine vertical-gap rows per request.  Each request
    // owns its workspace, so requests can execute independently on CPU/GPU.
    // The Java/GATK matrix uses reference rows and alternate(read) columns.
    // Keep the same orientation here so the endpoint/tie policy is shared.
    const std::size_t stride = max_read + 1;
    Matrix h_rows("sw_h_rows", requests.size(), 2, stride);
    Matrix vertical_rows("sw_vertical_rows", requests.size(), 2, stride);
    Matrix best_vertical_rows("sw_best_vertical_rows", requests.size(), 1, stride);
    auto h_reads = Kokkos::create_mirror_view(read_bases);
    auto h_references = Kokkos::create_mirror_view(reference_bases);
    auto h_read_index = Kokkos::create_mirror_view(read_index);
    auto h_reference_index = Kokkos::create_mirror_view(reference_index);
    for (std::size_t i = 0; i < reads.size(); ++i) h_reads(i) = reads[i];
    for (std::size_t i = 0; i < references.size(); ++i) h_references(i) = references[i];
    for (std::size_t i = 0; i < read_offsets.size(); ++i) h_read_index(i) = read_offsets[i];
    for (std::size_t i = 0; i < reference_offsets.size(); ++i) h_reference_index(i) = reference_offsets[i];
    Kokkos::deep_copy(read_bases, h_reads);
    Kokkos::deep_copy(reference_bases, h_references);
    Kokkos::deep_copy(read_index, h_read_index);
    Kokkos::deep_copy(reference_index, h_reference_index);
    ExecSpace().fence();

    fastgatk::core::DeviceBatch<ExecSpace> device(requests.size());
    device.bind("read_bases", read_bases);
    device.bind("reference_bases", reference_bases);
    device.bind("read_offsets", read_index);
    device.bind("reference_offsets", reference_index);
    device.bind("h_rows", h_rows);
    device.bind("vertical_rows", vertical_rows);
    device.bind("best_vertical_rows", best_vertical_rows);
    device.bind("scores", scores);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    plan.begin_execute();
    const int match = parameters.match;
    const int mismatch = parameters.mismatch;
    const int gap_open = parameters.gap_open;
    const int gap_extend = parameters.gap_extend;
    const int strategy = static_cast<int>(overhang);
    Kokkos::parallel_for("fastgatk_smith_waterman", Kokkos::RangePolicy<ExecSpace>(0, requests.size()),
        KOKKOS_LAMBDA(const std::size_t request) {
            const auto read_begin = read_index(request);
            const auto read_end = read_index(request + 1);
            const auto reference_begin = reference_index(request);
            const auto reference_end = reference_index(request + 1);
            const auto read_length = static_cast<std::size_t>(read_end - read_begin);
            const auto reference_length = static_cast<std::size_t>(reference_end - reference_begin);
            const int negative = kNegativeInfinity;
            const bool gap_edges = strategy == static_cast<int>(SmithWatermanOverhangStrategy::Indel) ||
                strategy == static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel);
            for (std::size_t j = 0; j <= read_length; ++j) {
                h_rows(request, 0, j) = 0;
                h_rows(request, 1, j) = 0;
                vertical_rows(request, 0, j) = negative;
                vertical_rows(request, 1, j) = negative;
                best_vertical_rows(request, 0, j) = negative;
            }
            if (gap_edges) {
                if (read_length >= 1) h_rows(request, 0, 1) = gap_open;
                for (std::size_t j = 2; j <= read_length; ++j)
                    h_rows(request, 0, j) = h_rows(request, 0, j - 1) + gap_extend;
            }
            int right_score = std::numeric_limits<int>::min();
            std::size_t right_row = 0;
            int bottom_score = std::numeric_limits<int>::min();
            std::size_t bottom_column = 0;
            int bottom_right_score = 0;
            for (std::size_t i = 1; i <= reference_length; ++i) {
                const int current = static_cast<int>(i & 1U);
                const int previous = current ^ 1;
                h_rows(request, current, 0) = 0;
                vertical_rows(request, current, 0) = negative;
                if (gap_edges) {
                    h_rows(request, current, 0) = i == 1 ? gap_open :
                        h_rows(request, previous, 0) + gap_extend;
                }
                int best_horizontal = negative;
                for (std::size_t j = 1; j <= read_length; ++j) {
                    const int diagonal = h_rows(request, previous, j - 1) +
                        (reference_bases(reference_begin + i - 1) == read_bases(read_begin + j - 1)
                             ? match : mismatch);
                    const int previous_vertical = h_rows(request, previous, j) + gap_open;
                    best_vertical_rows(request, 0, j) += gap_extend;
                    if (previous_vertical > best_vertical_rows(request, 0, j))
                        best_vertical_rows(request, 0, j) = previous_vertical;
                    const int vertical = best_vertical_rows(request, 0, j);
                    const int previous_horizontal = h_rows(request, current, j - 1) + gap_open;
                    best_horizontal += gap_extend;
                    if (previous_horizontal > best_horizontal)
                        best_horizontal = previous_horizontal;
                    const int right = best_horizontal;
                    int cell = diagonal;
                    if (!(diagonal >= vertical && diagonal >= right))
                        cell = right >= vertical ? right : vertical;
                    cell = Kokkos::max(kNegativeInfinity, cell);
                    h_rows(request, current, j) = cell;
                    vertical_rows(request, current, j) = vertical;
                    if (j == read_length) {
                        if (cell >= right_score) {
                            right_score = cell;
                            right_row = i;
                        }
                        if (i == reference_length) {
                            bottom_right_score = cell;
                            if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                                const auto candidate_distance = reference_length > j
                                    ? reference_length - j : j - reference_length;
                                if (cell > bottom_score ||
                                    (cell == bottom_score && candidate_distance <
                                     (reference_length > bottom_column ? reference_length - bottom_column :
                                      bottom_column - reference_length))) {
                                    bottom_score = cell;
                                    bottom_column = j;
                                }
                            }
                        }
                    } else if (i == reference_length &&
                               strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                        if (cell > bottom_score ||
                            (cell == bottom_score &&
                             (reference_length > j ? reference_length - j : j - reference_length) <
                             (reference_length > bottom_column ? reference_length - bottom_column :
                              bottom_column - reference_length))) {
                            bottom_score = cell;
                            bottom_column = j;
                        }
                    }
                }
            }
            int selected = bottom_right_score;
            if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::Indel)) {
                selected = right_score;
                if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel) &&
                    (bottom_score > selected ||
                     (bottom_score == selected &&
                      (reference_length > bottom_column ? reference_length - bottom_column :
                       bottom_column - reference_length) <
                      (right_row > read_length ? right_row - read_length : read_length - right_row))))
                    selected = bottom_score;
            }
            scores(request) = selected;
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_scores = Kokkos::create_mirror_view(scores);
    Kokkos::deep_copy(host_scores, scores);
    for (std::size_t i = 0; i < requests.size(); ++i) result.scores[i] = host_scores(i);
    return result;
}

// Kokkos batch alignment: single parallel DP fill that produces both
// score and backtrack matrices for every request in the chunk, then
// scalar traceback walks on the deep-copied data to build CIGARs.
//
// The DP recurrence (and tie policy: diagonal >= vertical, diagonal >=
// right, then right >= vertical) mirrors the scalar
// `smith_waterman::calculate_matrix` so the resulting scores and
// backtrack values are bit-identical to the scalar path.  The scalar
// `smith_waterman::traceback` walk then consumes the deep-copied backtrack
// matrix.  Host-side CIGAR string building remains the same as the
// scalar path; the only savings are the elimination of the per-request
// `std::vector::assign(rows*cols, 0)` for score + backtrack and the
// per-request DP-fill re-execution that the existing
// "score then traceback" two-pass does.
//
// Local helper: same semantics as the `release_views` defined in
// `pairhmm_kokkos.cpp`, inlined here so this TU doesn't have to reach
// into another translation unit's anonymous namespace.  Hands the
// View back to the OpenMP memory pool before its C++ scope exits; a
// no-op on non-OpenMP backends where `Kokkos::resize` has no slab.
namespace align_kokkos_detail {
template <typename View>
inline void shrink_release(View& view) noexcept {
    if (view.span() != 0) Kokkos::resize(view, 0);
}
template <typename... Views>
inline void release_views(Views&... views) noexcept {
    (shrink_release(views), ...);
}
}  // namespace align_kokkos_detail

// Forward declaration of the chunk helper, which is defined inside
// the chunked-implementation anonymous namespace below.  The public
// `smith_waterman_align_kokkos` calls into it through this declaration
// so the chunking loop stays out-of-line.
std::vector<SmithWatermanAlignment> smith_waterman_align_kokkos_chunk_impl(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy overhang);

std::vector<SmithWatermanAlignment> smith_waterman_align_kokkos(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy overhang) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (requests.empty()) return {};
    for (const auto& request : requests) {
        if (request.read.empty() || request.reference.empty())
            throw std::invalid_argument("Smith-Waterman request must be non-empty");
        if (request.read.size() > UINT32_MAX || request.reference.size() > UINT32_MAX)
            throw std::invalid_argument("Smith-Waterman request exceeds uint32 offsets");
    }

    // Chunk the batch the same way `smith_waterman_score_kokkos` does, so
    // workspace peak matches the score path.  Each chunk is processed
    // independently and results are concatenated in input order.
    constexpr std::size_t request_chunk = 512;
    std::vector<SmithWatermanAlignment> alignments(requests.size());
    for (std::size_t begin = 0; begin < requests.size(); begin += request_chunk) {
        const auto end = std::min(requests.size(), begin + request_chunk);
        std::vector<SmithWatermanRequest> chunk(
            requests.begin() + static_cast<std::ptrdiff_t>(begin),
            requests.begin() + static_cast<std::ptrdiff_t>(end));
        const auto chunk_alignments =
            smith_waterman_align_kokkos_chunk_impl(chunk, parameters, overhang);
        std::copy(chunk_alignments.begin(), chunk_alignments.end(),
                  alignments.begin() + static_cast<std::ptrdiff_t>(begin));
    }
    return alignments;
}

// Matches GATK's `SmithWatermanJavaAligner.traceback`.  Walks the backtrack
// matrix from the endpoint produced by the Kokkos DP fill, accumulates runs
// of M / I / D / S, and produces the CIGAR string and the
// `SmithWatermanAlignment` result for one request.  Mirrors the scalar
// `smith_waterman.cpp::traceback` so the byte stream matches for the same
// backtrack values, which is what guarantees bit-identical CIGARs vs the
// scalar `smith_waterman_align_reference` path.
SmithWatermanAlignment traceback_one(
    const std::uint8_t* read, std::size_t read_length,
    const std::uint8_t* /*reference*/, std::size_t reference_length,
    std::size_t endpoint_row, std::size_t endpoint_col,
    int endpoint_score,
    const Kokkos::View<const int***, Kokkos::LayoutRight, Kokkos::HostSpace>& trace_matrix,
    std::size_t request_index,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy strategy) {
    enum class TraceState : std::uint8_t { Match, Insertion, Deletion, Clip };
    SmithWatermanAlignment result;
    result.score = endpoint_score;
    std::size_t row = endpoint_row;
    std::size_t column = endpoint_col;
    std::vector<std::pair<TraceState, std::size_t>> reverse_segments;
    const auto trailing_overhang = (reference_length) - endpoint_col;
    std::size_t segment_length = trailing_overhang;
    if (strategy == SmithWatermanOverhangStrategy::Softclip && trailing_overhang != 0) {
        reverse_segments.emplace_back(TraceState::Clip, trailing_overhang);
        segment_length = 0;
    }
    TraceState state = TraceState::Match;
    while (row > 0 && column > 0) {
        const int backtrack = trace_matrix(request_index, row, column);
        TraceState next_state = TraceState::Match;
        std::size_t step_length = 1;
        if (backtrack > 0) {
            next_state = TraceState::Deletion;
            step_length = static_cast<std::size_t>(backtrack);
        } else if (backtrack < 0) {
            next_state = TraceState::Insertion;
            step_length = static_cast<std::size_t>(-backtrack);
        }
        if (next_state == state) {
            segment_length += step_length;
        } else {
            if (segment_length != 0) reverse_segments.emplace_back(state, segment_length);
            segment_length = step_length;
            state = next_state;
        }
        switch (next_state) {
            case TraceState::Match: --row; --column; break;
            case TraceState::Insertion: column -= step_length; break;
            case TraceState::Deletion: row -= step_length; break;
            case TraceState::Clip: break;
        }
    }
    if (strategy == SmithWatermanOverhangStrategy::Softclip) {
        reverse_segments.emplace_back(state, segment_length);
        if (column > 0) reverse_segments.emplace_back(TraceState::Clip, column);
        result.read_start = row;
        result.reference_start = row;
    } else if (strategy == SmithWatermanOverhangStrategy::Ignore) {
        reverse_segments.emplace_back(state, segment_length + column);
        result.read_start = row;
        result.reference_start = row;
    } else {
        reverse_segments.emplace_back(state, segment_length);
        if (row > 0) {
            reverse_segments.emplace_back(TraceState::Deletion, row);
            // Indel/LeadingIndel strategies consume the remaining alternate
            // (read) prefix as a leading insertion, matching the scalar
            // `traceback` semantics that push `Insertion(column)` when
            // `row == 0 && column > 0`.  When both row and column are zero
            // here we land at (0, 0) and the for-loop above already
            // captured any residue in `segment_length`.
        } else if (column > 0) {
            reverse_segments.emplace_back(TraceState::Insertion, column);
        }
        result.read_start = 0;
        result.reference_start = 0;
    }
    std::reverse(reverse_segments.begin(), reverse_segments.end());

    // Compact CIGAR.  Mirrors the scalar `smith_waterman.cpp::compact_cigar`.
    std::string cigar;
    TraceState previous = TraceState::Clip;
    std::size_t run = 0;
    auto flush = [&]() {
        if (run == 0) return;
        cigar += std::to_string(run);
        switch (previous) {
            case TraceState::Match: cigar.push_back('M'); break;
            case TraceState::Insertion: cigar.push_back('I'); break;
            case TraceState::Deletion: cigar.push_back('D'); break;
            case TraceState::Clip: cigar.push_back('S'); break;
        }
        run = 0;
    };
    for (const auto& [seg_state, length] : reverse_segments) {
        if (length == 0) continue;
        if (run != 0 && seg_state != previous) flush();
        previous = seg_state;
        run += length;
    }
    flush();
    result.cigar = std::move(cigar);
    result.read_end = endpoint_col;
    result.reference_end = endpoint_row;
    if (strategy == SmithWatermanOverhangStrategy::Softclip)
        result.alignment_offset = static_cast<int>(result.reference_start);
    else if (strategy == SmithWatermanOverhangStrategy::Ignore)
        result.alignment_offset = static_cast<int>(result.reference_start) -
            static_cast<int>(result.read_start);
    (void)parameters; (void)read; (void)read_length;
    return result;
}

// Chunked implementation.  Allocates one set of Kokkos Views per chunk,
// runs the parallel DP fill, deep-copies the backtrack matrix, runs
// scalar traceback walks on the deep-copied data, and returns the
// resulting `SmithWatermanAlignment` vector.  External linkage so the
// public `smith_waterman_align_kokkos` and the HC caller
// (calling_pipeline.cpp) can link against it.  Implementation mirrors
// the scalar `smith_waterman::calculate_matrix` + `traceback` so the
// resulting scores and CIGARs are bit-identical to the scalar path.
std::vector<SmithWatermanAlignment> smith_waterman_align_kokkos_chunk_impl(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters,
    SmithWatermanOverhangStrategy overhang) {
    if (requests.empty()) return {};
    std::size_t max_read = 0, max_reference = 0;
    std::size_t host_bytes = 0;
    std::vector<std::uint32_t> read_offsets(requests.size() + 1, 0);
    std::vector<std::uint32_t> reference_offsets(requests.size() + 1, 0);
    std::vector<std::uint8_t> reads, references;
    for (std::size_t i = 0; i < requests.size(); ++i) {
        const auto& request = requests[i];
        max_read = std::max(max_read, request.read.size());
        max_reference = std::max(max_reference, request.reference.size());
        read_offsets[i + 1] = read_offsets[i] +
            static_cast<std::uint32_t>(request.read.size());
        reference_offsets[i + 1] = reference_offsets[i] +
            static_cast<std::uint32_t>(request.reference.size());
        reads.insert(reads.end(), request.read.begin(), request.read.end());
        references.insert(references.end(),
            request.reference.begin(), request.reference.end());
    }
    host_bytes = reads.size() + references.size() +
                 (read_offsets.size() + reference_offsets.size()) * sizeof(std::uint32_t);
    fastgatk::core::HostBatch host("smith-waterman-align-v1");
    host.records = requests.size();
    host.bytes = host_bytes;
    fastgatk::core::KernelPlan<ExecSpace> plan("smith-waterman-align");
    plan.begin_prepare(host);

    Kokkos::View<std::uint8_t*, MemorySpace> read_bases(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "sw_align_read_bases"),
        reads.size());
    Kokkos::View<std::uint8_t*, MemorySpace> reference_bases(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "sw_align_reference_bases"),
        references.size());
    Kokkos::View<std::uint32_t*, MemorySpace> read_index("sw_align_read_offsets",
        read_offsets.size());
    Kokkos::View<std::uint32_t*, MemorySpace> reference_index("sw_align_reference_offsets",
        reference_offsets.size());
    Kokkos::View<int*, MemorySpace> scores("sw_align_scores", requests.size());
    using ScoreMatrix = Kokkos::View<int***, Kokkos::LayoutRight, MemorySpace>;
    const std::size_t stride = max_read + 1;
    ScoreMatrix h_rows("sw_align_h_rows", requests.size(), 2, stride);
    ScoreMatrix vertical_rows("sw_align_vertical_rows", requests.size(), 2, stride);
    ScoreMatrix best_vertical_rows("sw_align_best_vertical_rows", requests.size(), 1, stride);
    // One backtrack entry per DP cell, row-major by (request, ref_row, read_col).
    // Sourced from the same `vertical_gap_size`/`horizontal_gap_size`
    // bookkeeping as `smith_waterman::calculate_matrix`, so the
    // traceback walk yields the same CIGAR as the scalar path.
    using TraceMatrix = Kokkos::View<int***, Kokkos::LayoutRight, MemorySpace>;
    const std::size_t trace_stride = max_reference + 1;
    TraceMatrix trace_matrix("sw_align_trace_matrix", requests.size(), trace_stride, stride);
    auto h_reads = Kokkos::create_mirror_view(read_bases);
    auto h_references = Kokkos::create_mirror_view(reference_bases);
    auto h_read_index = Kokkos::create_mirror_view(read_index);
    auto h_reference_index = Kokkos::create_mirror_view(reference_index);
    for (std::size_t i = 0; i < reads.size(); ++i) h_reads(i) = reads[i];
    for (std::size_t i = 0; i < references.size(); ++i) h_references(i) = references[i];
    for (std::size_t i = 0; i < read_offsets.size(); ++i) h_read_index(i) = read_offsets[i];
    for (std::size_t i = 0; i < reference_offsets.size(); ++i) h_reference_index(i) = reference_offsets[i];
    Kokkos::deep_copy(read_bases, h_reads);
    Kokkos::deep_copy(reference_bases, h_references);
    Kokkos::deep_copy(read_index, h_read_index);
    Kokkos::deep_copy(reference_index, h_reference_index);
    ExecSpace().fence();

    fastgatk::core::DeviceBatch<ExecSpace> device(requests.size());
    device.bind("read_bases", read_bases);
    device.bind("reference_bases", reference_bases);
    device.bind("read_offsets", read_index);
    device.bind("reference_offsets", reference_index);
    device.bind("h_rows", h_rows);
    device.bind("vertical_rows", vertical_rows);
    device.bind("best_vertical_rows", best_vertical_rows);
    device.bind("scores", scores);
    device.bind("trace_matrix", trace_matrix);
    plan.end_prepare(device);

    // Endpoint selection on the device; mirrors `choose_endpoint` in
    // smith_waterman.cpp.  Endpoint selection feeds back to the
    // Host-side traceback walks.
    Kokkos::View<int*, MemorySpace> endpoint_score("sw_align_endpoint_score", requests.size());
    Kokkos::View<std::uint32_t*, MemorySpace> endpoint_row("sw_align_endpoint_row", requests.size());
    Kokkos::View<std::uint32_t*, MemorySpace> endpoint_col("sw_align_endpoint_col", requests.size());
    device.bind("endpoint_score", endpoint_score);
    device.bind("endpoint_row", endpoint_row);
    device.bind("endpoint_col", endpoint_col);
    Kokkos::deep_copy(endpoint_score, std::numeric_limits<int>::min());

    plan.begin_execute();
    const int match = parameters.match;
    const int mismatch = parameters.mismatch;
    const int gap_open = parameters.gap_open;
    const int gap_extend = parameters.gap_extend;
    const int strategy = static_cast<int>(overhang);
    Kokkos::parallel_for("fastgatk_smith_waterman_align",
        Kokkos::RangePolicy<ExecSpace>(0, requests.size()),
        KOKKOS_LAMBDA(const std::size_t request) {
            const auto read_begin = read_index(request);
            const auto read_end = read_index(request + 1);
            const auto reference_begin = reference_index(request);
            const auto reference_end = reference_index(request + 1);
            const auto read_length = static_cast<std::size_t>(read_end - read_begin);
            const auto reference_length = static_cast<std::size_t>(reference_end - reference_begin);
            const int negative = kNegativeInfinity;
            const bool gap_edges =
                strategy == static_cast<int>(SmithWatermanOverhangStrategy::Indel) ||
                strategy == static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel);
            for (std::size_t j = 0; j <= read_length; ++j) {
                h_rows(request, 0, j) = 0;
                h_rows(request, 1, j) = 0;
                vertical_rows(request, 0, j) = negative;
                vertical_rows(request, 1, j) = negative;
                best_vertical_rows(request, 0, j) = negative;
            }
            if (gap_edges) {
                if (read_length >= 1) h_rows(request, 0, 1) = gap_open;
                for (std::size_t j = 2; j <= read_length; ++j)
                    h_rows(request, 0, j) = h_rows(request, 0, j - 1) + gap_extend;
            }
            int right_score = std::numeric_limits<int>::min();
            std::size_t right_row = 0;
            int bottom_score = std::numeric_limits<int>::min();
            std::size_t bottom_column = 0;
            int bottom_right_score = 0;
            std::size_t bottom_right_row = reference_length;
            std::size_t bottom_right_col = read_length;
            for (std::size_t i = 1; i <= reference_length; ++i) {
                const int current = static_cast<int>(i & 1U);
                const int previous = current ^ 1;
                h_rows(request, current, 0) = 0;
                vertical_rows(request, current, 0) = negative;
                if (gap_edges) {
                    h_rows(request, current, 0) = i == 1 ? gap_open :
                        h_rows(request, previous, 0) + gap_extend;
                }
                int best_horizontal = negative;
                int best_horizontal_size = 0;
                int previous_horizontal_size = 0;
                for (std::size_t j = 1; j <= read_length; ++j) {
                    const int diagonal = h_rows(request, previous, j - 1) +
                        (reference_bases(reference_begin + i - 1) == read_bases(read_begin + j - 1)
                             ? match : mismatch);
                    const int previous_vertical = h_rows(request, previous, j) + gap_open;
                    int best_vertical_local = best_vertical_rows(request, 0, j) + gap_extend;
                    int best_vertical_local_size = previous_horizontal_size + 1;
                    if (previous_vertical > best_vertical_local) {
                        best_vertical_local = previous_vertical;
                        best_vertical_local_size = 1;
                    }
                    best_vertical_rows(request, 0, j) = best_vertical_local;
                    vertical_rows(request, current, j) = best_vertical_local;
                    const int vertical = best_vertical_local;
                    const int previous_horizontal = h_rows(request, current, j - 1) + gap_open;
                    best_horizontal += gap_extend;
                    int local_best_horizontal_size = previous_horizontal_size + 1;
                    if (previous_horizontal > best_horizontal) {
                        best_horizontal = previous_horizontal;
                        local_best_horizontal_size = 1;
                    }
                    previous_horizontal_size = local_best_horizontal_size;
                    const int right = best_horizontal;
                    int cell = diagonal;
                    int trace_value = 0;
                    if (diagonal >= vertical && diagonal >= right) {
                        cell = diagonal;
                        trace_value = 0;
                    } else if (right >= vertical) {
                        cell = right;
                        trace_value = -best_horizontal_size;
                    } else {
                        cell = vertical;
                        trace_value = best_vertical_local_size;
                    }
                    cell = Kokkos::max(kNegativeInfinity, cell);
                    h_rows(request, current, j) = cell;
                    trace_matrix(request, i, j) = trace_value;
                    if (j == read_length) {
                        if (cell >= right_score) {
                            right_score = cell;
                            right_row = i;
                        }
                        if (i == reference_length) {
                            bottom_right_score = cell;
                            bottom_right_row = i;
                            bottom_right_col = j;
                            if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                                const auto candidate_distance = reference_length > j
                                    ? reference_length - j : j - reference_length;
                                if (cell > bottom_score ||
                                    (cell == bottom_score && candidate_distance <
                                     (reference_length > bottom_column ? reference_length - bottom_column :
                                      bottom_column - reference_length))) {
                                    bottom_score = cell;
                                    bottom_column = j;
                                }
                            }
                        }
                    } else if (i == reference_length &&
                               strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel)) {
                        if (cell > bottom_score ||
                            (cell == bottom_score &&
                             (reference_length > j ? reference_length - j : j - reference_length) <
                             (reference_length > bottom_column ? reference_length - bottom_column :
                              bottom_column - reference_length))) {
                            bottom_score = cell;
                            bottom_column = j;
                        }
                    }
                }
                best_horizontal_size = previous_horizontal_size;
            }
            int selected = bottom_right_score;
            std::size_t selected_row = bottom_right_row;
            std::size_t selected_col = bottom_right_col;
            if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::Indel)) {
                selected = right_score;
                selected_row = right_row;
                selected_col = read_length;
                if (strategy != static_cast<int>(SmithWatermanOverhangStrategy::LeadingIndel) &&
                    (bottom_score > selected ||
                     (bottom_score == selected &&
                      (reference_length > bottom_column ? reference_length - bottom_column :
                       bottom_column - reference_length) <
                      (selected_row > read_length ? selected_row - read_length : read_length - selected_row))))
                    selected = bottom_score;
            }
            scores(request) = selected;
            endpoint_score(request) = selected;
            endpoint_row(request) = static_cast<std::uint32_t>(selected_row);
            endpoint_col(request) = static_cast<std::uint32_t>(selected_col);
        });
    ExecSpace().fence();
    plan.end_execute();

    auto host_scores = Kokkos::create_mirror_view(scores);
    auto host_endpoint_score = Kokkos::create_mirror_view(endpoint_score);
    auto host_endpoint_row = Kokkos::create_mirror_view(endpoint_row);
    auto host_endpoint_col = Kokkos::create_mirror_view(endpoint_col);
    auto host_trace = Kokkos::create_mirror_view(trace_matrix);
    Kokkos::deep_copy(host_scores, scores);
    Kokkos::deep_copy(host_endpoint_score, endpoint_score);
    Kokkos::deep_copy(host_endpoint_row, endpoint_row);
    Kokkos::deep_copy(host_endpoint_col, endpoint_col);
    Kokkos::deep_copy(host_trace, trace_matrix);
    // Bind lifetime: every Kokkos View allocated in this function goes
    // out of scope when we return; release them eagerly back to the slab
    // to avoid cumulative growth across many HC regions.
    align_kokkos_detail::release_views(
        read_bases, reference_bases, read_index, reference_index,
        scores, h_rows, vertical_rows, best_vertical_rows,
        trace_matrix, endpoint_score, endpoint_row, endpoint_col);

    // Scalar traceback walks on the deep-copied backtrack matrix.  Bit
    // match vs `smith_waterman::traceback` is the contract; any drift
    // will surface as a numerical-contract failure at the HC caller.
    std::vector<SmithWatermanAlignment> alignments(requests.size());
    for (std::size_t request = 0; request < requests.size(); ++request) {
        const auto& req = requests[request];
        alignments[request] = traceback_one(
            req.read.data(), req.read.size(),
            req.reference.data(), req.reference.size(),
            host_endpoint_row(request),
            host_endpoint_col(request),
            host_endpoint_score(request),
            host_trace, request,
            parameters, overhang);
    }
    // Numerical-contract assertion.  HC's existing score-then-traceback
    // path asserts the same invariant; we keep it for parity.
    for (std::size_t request = 0; request < requests.size(); ++request) {
        if (alignments[request].score != host_scores(request))
            throw std::runtime_error(
                "NUMERICAL_CONTRACT_FAILURE: SW batch score/traceback mismatch");
    }
    return alignments;
}

}  // namespace fastgatk::kernels
