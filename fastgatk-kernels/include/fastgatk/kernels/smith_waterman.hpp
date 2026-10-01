#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

struct SmithWatermanParameters {
    int match = 10;
    int mismatch = -15;
    int gap_open = -30;
    int gap_extend = -5;
};

// Matches GATK's SWOverhangStrategy.  The first argument to the public API is
// the alternate/read sequence and the second is the reference/haplotype;
// consequently I/D operators use the usual read-vs-reference CIGAR meaning.
enum class SmithWatermanOverhangStrategy : std::uint8_t {
    Softclip,
    Indel,
    LeadingIndel,
    Ignore,
};

// Host-side traceback result.  The score recurrence is reusable by device
// kernels, while variable-length CIGAR construction intentionally stays on
// Host.  Coordinates use half-open [start,end) intervals.
struct SmithWatermanAlignment {
    int score = 0;
    std::size_t read_start = 0;
    std::size_t read_end = 0;
    std::size_t reference_start = 0;
    std::size_t reference_end = 0;
    // GATK's alignment offset: SOFTCLIP uses the reference start, IGNORE
    // uses reference_start-read_start, and indel modes use zero.
    int alignment_offset = 0;
    std::string cigar;
};

int smith_waterman_score_reference(const std::uint8_t* read, std::size_t read_length,
                                   const std::uint8_t* reference, std::size_t reference_length,
                                   SmithWatermanParameters parameters = {},
                                   SmithWatermanOverhangStrategy overhang =
                                       SmithWatermanOverhangStrategy::Softclip);

SmithWatermanAlignment smith_waterman_align_reference(
    const std::uint8_t* read, std::size_t read_length,
    const std::uint8_t* reference, std::size_t reference_length,
    SmithWatermanParameters parameters = {},
    SmithWatermanOverhangStrategy overhang = SmithWatermanOverhangStrategy::Softclip);

// Kokkos batch score API for variable-length requests.  The affine DP score
// runs through the selected Kokkos execution space; variable-length CIGAR
// traceback remains on Host in smith_waterman_align_reference().
struct SmithWatermanRequest {
    std::vector<std::uint8_t> read;
    std::vector<std::uint8_t> reference;
};

struct SmithWatermanBatchResult {
    std::vector<int> scores;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    // Width selected by Kokkos::Experimental::simd for the score kernel.
    // A value of one is the portable scalar ABI (including accelerator
    // execution spaces where Kokkos intentionally maps SIMD to scalar).
    std::size_t simd_width = 1;
    std::size_t simd_groups = 0;
};

SmithWatermanBatchResult smith_waterman_score_kokkos(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters = {},
    SmithWatermanOverhangStrategy overhang = SmithWatermanOverhangStrategy::Softclip);

// Kokkos batch score + traceback + CIGAR.  Backed by a single parallel DP
// fill over the request set: the affine DP matrix and its backtrack matrix
// live in Kokkos Views, and only the variable-length CIGAR traceback walk
// runs on Host after a deep_copy.  Equivalent to calling
// `smith_waterman_align_reference` on every request but avoids the
// per-request `std::vector::assign(rows*cols, 0)` allocation that the
// scalar traceback uses, so the per-request overhead drops from O(mn) of
// `memset(0)` and reallocation to a single Kokkos View deep_copy.
//
// Score `result.scores[i]` is bit-identical to what
// `smith_waterman_align_reference(...)` returns for request `i`.  CIGAR and
// alignment offsets follow GATK's `SmithWatermanJavaAligner` tie policy:
// diagonal >= right/down, then right >= down, then tie-broken by the
// closest-to-diagonal endpoint.
//
// Bounded by `requests.size() <= std::numeric_limits<std::uint32_t>::max()`.
std::vector<SmithWatermanAlignment> smith_waterman_align_kokkos(
    const std::vector<SmithWatermanRequest>& requests,
    SmithWatermanParameters parameters = {},
    SmithWatermanOverhangStrategy overhang = SmithWatermanOverhangStrategy::Softclip);

}  // namespace fastgatk::kernels
