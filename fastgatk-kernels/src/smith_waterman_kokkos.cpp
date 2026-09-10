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

}  // namespace fastgatk::kernels
