#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// Host-side observations after CIGAR/reference projection and known-sites
// masking.  The kernel deliberately receives only flat primitive arrays: all
// read-group/cycle/context semantics remain in the Host tool.
struct BqsrObservationBatch {
    std::vector<std::uint8_t> qualities;
    std::vector<std::uint8_t> mismatches;
};

// Integer identifiers for Host-resolved BQSR covariate keys.  The Host owns
// string/context/cycle construction and assigns a dense id for the current
// decoded batch; the device only sees primitive arrays and therefore remains
// independent of HTSlib and allocation-heavy std::map state.
struct BqsrCovariateObservationBatch {
    std::vector<std::uint32_t> covariate_ids;
    std::vector<std::uint8_t> mismatches;
};

struct BqsrQualityCounts {
    std::uint64_t count = 0;
    std::uint64_t mismatches = 0;
};

struct BqsrCountResult {
    std::vector<BqsrQualityCounts> bins;
    std::size_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space;
    std::string execution_policy;
};

struct BqsrCovariateCountResult {
    std::vector<BqsrQualityCounts> covariates;
    std::size_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::uint64_t workspace_bytes = 0;
    bool team_local_histogram = false;
    std::string execution_space;
    std::string execution_policy;
};

struct BqsrQualityTransformResult {
    std::vector<std::uint8_t> adjusted;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space;
    std::string execution_policy;
};

// Count raw-quality and mismatch observations with Kokkos.  Each quality bin
// uses TeamPolicy-local integer tables and a Host league-order merge, so
// backend parallel order cannot alter the result. The API is batch-oriented
// and can be called once per decoded read batch; it does not own HTSlib
// objects or reference strings.
BqsrCountResult count_bqsr_quality_kokkos(const BqsrObservationBatch& input);

// Aggregate Host-assigned covariate ids with deterministic integer counters.
// The implementation uses a bounded TeamPolicy-local histogram when its
// workspace fits the safety cap, then merges team rows in league order.  Very
// wide/sparse batches use a bounded RangePolicy atomic fallback instead of
// allocating an unbounded per-team table.  Addition is exact, so the result is
// independent of execution-space scheduling; the Host merges dense rows into
// its stable CovariateKey map.
BqsrCovariateCountResult count_bqsr_covariates_kokkos(
    const BqsrCovariateObservationBatch& input);

// Apply already-resolved Host covariate deltas to a flat quality batch.  The
// Host decides which covariate row wins; this kernel only performs the exact
// integer clamp used by the output tool.
BqsrQualityTransformResult apply_bqsr_quality_kokkos(
    const std::vector<std::uint8_t>& qualities,
    const std::vector<std::int16_t>& deltas);

}  // namespace fastgatk::kernels
