#pragma once

#include "fastgatk/kernels/pairhmm.hpp"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace fastgatk::pairhmm {

struct PairIndexBatch {
    std::vector<std::uint32_t> read_ids;
    std::vector<std::uint32_t> haplotype_ids;
};

// Variable-length Host input for real callers.  The low-level SIMD kernel
// uses one rectangular length bucket; the bucketed entry point below groups
// requests by (read_length, haplotype_length), preserving original request
// order while avoiding padding-induced likelihood changes.
struct PairHmmRead {
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    std::vector<std::uint8_t> insertion_gop;
    std::vector<std::uint8_t> deletion_gop;
    std::vector<std::uint8_t> gap_continuation;
};

struct PairHmmHaplotype {
    std::vector<std::uint8_t> bases;
};

struct PairHmmRequest {
    std::uint32_t read_id = 0;
    std::uint32_t haplotype_id = 0;
};

// Arithmetic mode for the portable Kokkos PairHMM kernel.  Float32 mirrors
// GKL's default execution: it runs the Float32 recurrence first and reruns
// only scaled sums below GKL's 1e-28 threshold in Float64.  Float64 requests
// the direct Java-compatible double recurrence for callers that explicitly
// select --native-pair-hmm-use-double-precision.
enum class PairHmmPrecision { Float32, Float64 };

struct KokkosBatchResult : BatchResult {
    std::vector<double> scaled_sums;
    double prepare_seconds = 0.0;
    std::size_t simd_width = 1;
    std::size_t cached_shapes = 0;
    std::size_t cache_hits = 0;
    // Non-zero only for compute_kokkos_full_matrix().  Likelihoods and
    // scaled_sums are then row-major [read_id * matrix_columns + haplotype_id].
    std::size_t matrix_rows = 0;
    std::size_t matrix_columns = 0;
    std::string execution_space;
    // Kokkos execution policy used by the numerical kernel.  Keeping this in
    // the result makes backend/benchmark reports auditable instead of
    // inferring a GPU launch policy from the execution-space name.
    std::string execution_policy = "RangePolicy";
    // Explicit model tag.  Regular and flow-space PairHMM are different
    // recurrences; callers and manifests must never infer one from a generic
    // boolean or silently run regular bases through a flow request.
    std::string error_model = "regular";
    // Arithmetic provenance for audit/benchmark output.  "float32" denotes
    // GKL-compatible Float32-first execution, including its per-pair double
    // underflow fallback; "double" denotes a direct Float64 recurrence.
    std::string precision = "double";
};

// Host-facing GATK AlleleLikelihoods normalization boundary.  The numerical
// work (per-read best reduction and global-mismapping cap) is executed through
// Kokkos so CPU and device callers share one implementation; CIGAR/allele
// semantics remain on Host.  `eligible_for_best` mirrors the caller's retained
// haplotype set (for example, graph paths can be scored but excluded from the
// poorly-modelled-read baseline until their pruning contract is complete).
struct LikelihoodNormalizationResult {
    std::vector<double> likelihoods;
    std::vector<double> best_by_read;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "RangePolicy";
};

LikelihoodNormalizationResult normalize_likelihoods_kokkos(
    const std::vector<double>& likelihoods,
    const std::vector<std::uint32_t>& read_ids,
    const std::vector<std::uint8_t>& eligible_for_best,
    std::size_t read_count,
    // Positive values are the log10-width of GATK's lower likelihood floor.
    // Negative infinity is the release-compatible disabled-cap sentinel;
    // finite negative values and NaN remain invalid.
    double maximum_likelihood_difference_cap = 4.5);

// Reduce the read-by-haplotype matrix to concrete alleles for a candidate-read
// row. Host supplies the CIGAR/graph-filtered rows and a zero-based concrete
// allele index (0=REF, 1..N=ALT); Kokkos performs the deterministic
// max-marginalization that mirrors AlleleLikelihoods.marginalize().
struct AlleleMarginalizationRequest {
    std::uint32_t row_id = 0;
    std::uint8_t allele = 0;
    double likelihood = 0.0;
};

struct AlleleMarginalizationResult {
    // Row-major [row_id * allele_count + allele], initialized to -infinity
    // when no retained haplotype supports that allele.
    std::vector<double> best_by_row_allele;
    std::size_t allele_count = 2;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "RangePolicy";
};

AlleleMarginalizationResult marginalize_read_allele_likelihoods_kokkos(
    const std::vector<AlleleMarginalizationRequest>& requests,
    std::size_t row_count,
    std::size_t allele_count);

// Backward-compatible biallelic convenience overload.  New callers that
// retain more than one concrete ALT must pass the concrete allele count so a
// read-by-allele matrix is not silently collapsed to REF/ALT.
AlleleMarginalizationResult marginalize_read_allele_likelihoods_kokkos(
    const std::vector<AlleleMarginalizationRequest>& requests,
    std::size_t row_count);

// Sum read log-likelihoods into fragment x haplotype cells before allele
// marginalization.  Requests for one cell are added in their original Host
// order, matching AlleleLikelihoods.groupEvidence().  In particular,
// -infinity is a real zero-probability value and is deliberately propagated;
// it must not be treated as a missing read.
struct FragmentHaplotypeAggregationRequest {
    std::uint32_t cell_id = 0;
    double likelihood = 0.0;
};

struct FragmentHaplotypeAggregationResult {
    // One sum per fragment x haplotype cell. Empty cells retain 0.0, matching
    // Java's freshly allocated likelihood matrix before evidence is added.
    std::vector<double> sums_by_cell;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "RangePolicy";
};

FragmentHaplotypeAggregationResult aggregate_fragment_haplotype_likelihoods_kokkos(
    const std::vector<FragmentHaplotypeAggregationRequest>& requests,
    std::size_t cell_count);

// Reduce the concrete allele likelihoods for one (locus, read) row to the
// best and second-best values used by GATK's BestAllele confidence margin.
// Host owns locus/read row construction and stable tie ordering; Kokkos owns
// the deterministic top-two reduction so regular and flow PairHMM callers
// share the same uncertainty kernel.
struct ReadAlleleUncertaintyRequest {
    std::uint32_t row_id = 0;
    double likelihood = 0.0;
};

struct ReadAlleleUncertaintyResult {
    // One pair per row: [best, second-best]. Missing entries remain -infinity.
    std::vector<double> best_second_by_row;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "RangePolicy";
};

ReadAlleleUncertaintyResult reduce_read_allele_uncertainty_kokkos(
    const std::vector<ReadAlleleUncertaintyRequest>& requests,
    std::size_t row_count);

// Select the best concrete allele for each (locus, read) row after allele
// max-marginalization.  Requests are sorted by row and allele index before
// they enter the kernel, so strict `>` comparisons preserve GATK's stable
// REF-first/ALT-order tie behavior.  The returned indices are local to the
// caller's row (REF=0, concrete ALT=1..N); missing rows use UINT32_MAX.
struct ReadAlleleBestRequest {
    std::uint32_t row_id = 0;
    std::uint32_t allele_index = 0;
    double likelihood = 0.0;
};

struct ReadAlleleBestResult {
    // Row-major [row_id * 2 + {best, second}] allele indices.
    std::vector<std::uint32_t> best_second_allele_by_row;
    // Row-major [row_id * 2 + {best, second}] likelihoods.
    std::vector<double> best_second_likelihood_by_row;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "RangePolicy";
};

ReadAlleleBestResult reduce_read_allele_best_kokkos(
    const std::vector<ReadAlleleBestRequest>& requests,
    std::size_t row_count);

KokkosBatchResult compute_kokkos(const std::vector<PairInput>& records,
                                 const PairIndexBatch& pairs,
                                 int iterations = 1,
                                 bool tristate_correction = true);

// Execute variable-length read/haplotype requests through deterministic
// rectangular buckets.  Returned likelihoods/scaled sums use request order.
KokkosBatchResult compute_kokkos_bucketed(
    const std::vector<PairHmmRead>& reads,
    const std::vector<PairHmmHaplotype>& haplotypes,
    const std::vector<PairHmmRequest>& requests,
    int iterations = 1,
    PairHmmPrecision precision = PairHmmPrecision::Float64);

// Evaluate the complete ragged read-by-haplotype matrix.  Requests are
// internally length-bucketed and executed through the same persistent Kokkos
// workspace plan; no sequence padding is introduced.  The returned vectors
// are row-major and therefore retain the natural GATK matrix ordering.
KokkosBatchResult compute_kokkos_full_matrix(
    const std::vector<PairHmmRead>& reads,
    const std::vector<PairHmmHaplotype>& haplotypes,
    int iterations = 1);

// Flow-space PairHMM input.  `key` contains homopolymer flow lengths and
// `flow_order` contains the corresponding flow base/order symbols.  The
// caller supplies the calibrated per-read-flow probability table flattened
// as [read_flow][256], matching FlowBasedRead.getProb(flow, haplotype_flow).
// This keeps calibration/HTS tags on Host while the fixed-stride DP remains a
// portable Kokkos kernel.
struct FlowPairHmmRead {
    std::vector<std::int32_t> key;
    std::vector<std::uint8_t> flow_order;
    std::vector<std::uint8_t> insertion_gop;
    std::vector<std::uint8_t> deletion_gop;
    std::vector<std::uint8_t> gap_continuation;
    std::vector<double> probabilities;
};

struct FlowPairHmmHaplotype {
    std::vector<std::int32_t> key;
    std::vector<std::uint8_t> flow_order;
};

struct FlowPairHmmRequest {
    std::uint32_t read_id = 0;
    std::uint32_t haplotype_id = 0;
};

// Execute the independent flow-space recurrence.  Requests are length
// bucketed, preserve input order, and use a separate Kokkos kernel tag from
// regular PairHMM.  No regular-base approximation is used.
KokkosBatchResult compute_kokkos_flow(
    const std::vector<FlowPairHmmRead>& reads,
    const std::vector<FlowPairHmmHaplotype>& haplotypes,
    const std::vector<FlowPairHmmRequest>& requests,
    int iterations = 1);

// Execute GATK's FlowBasedAlignmentLikelihoodEngine (the `FlowBased`
// likelihood-calculation-engine).  Unlike FlowBasedHMM above, this is not a
// PairHMM recurrence: Host has already trimmed each haplotype to its
// alignment-aware flow window and the Kokkos kernel evaluates every legal
// phase-aligned placement, returning the greatest product of calibrated
// per-flow probabilities.  Keeping it separate is deliberate: FlowBased and
// FlowBasedHMM are distinct GATK engines and must not silently share the HMM
// scoring path.
KokkosBatchResult compute_kokkos_flow_alignment(
    const std::vector<FlowPairHmmRead>& reads,
    const std::vector<FlowPairHmmHaplotype>& haplotypes,
    const std::vector<FlowPairHmmRequest>& requests,
    int iterations = 1);

// Long-lived bucket dispatcher.  A plan retains Kokkos Views and the three-row
// DP workspace for each observed (read length, haplotype length, record count,
// pair count) shape.  Repeated execute calls still repack input values, but do
// not reallocate the device workspace; this is the production lifecycle used
// by region-level callers to amortize prepare/allocation costs.
class PersistentBucketPlan {
public:
    explicit PersistentBucketPlan(std::size_t max_cached_shapes = 8);
    ~PersistentBucketPlan();
    PersistentBucketPlan(PersistentBucketPlan&&) noexcept;
    PersistentBucketPlan& operator=(PersistentBucketPlan&&) noexcept;
    PersistentBucketPlan(const PersistentBucketPlan&) = delete;
    PersistentBucketPlan& operator=(const PersistentBucketPlan&) = delete;

    KokkosBatchResult execute(
        const std::vector<PairHmmRead>& reads,
        const std::vector<PairHmmHaplotype>& haplotypes,
        const std::vector<PairHmmRequest>& requests,
        int iterations = 1);
    void clear();
    std::size_t cached_shapes() const;
    std::size_t cache_hits() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

std::string kokkos_backend_description();

}  // namespace fastgatk::pairhmm
