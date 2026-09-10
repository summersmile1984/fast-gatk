#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace fastgatk::contamination {

// Telemetry for the two Kokkos stages used by KernelSegmenter. The SVD and
// change-point selection stay host-side for deterministic Java-compatible
// ordering, while matrix construction/projection use the shared lifecycle.
struct SegmenterTelemetry {
    std::size_t batches = 0;
    std::size_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space;
    std::string execution_policy = "MDRangePolicy";
};

// KernelSegmenter-compatible changepoints.  A returned index belongs to the
// segment on its left, matching GATK's ContaminationSegmenter convention.
std::vector<std::size_t> find_changepoints(const std::vector<double>& minor_allele_fractions,
                                           std::size_t max_num_changepoints = 10,
                                           std::size_t kernel_approximation_dimension = 100,
                                           std::size_t window_size = 50,
                                           double linear_penalty = 1.0,
                                           double log_linear_penalty = 1.0,
                                           SegmenterTelemetry* telemetry = nullptr);

std::vector<std::size_t> find_changepoints_with_windows(
    const std::vector<double>& minor_allele_fractions,
    std::size_t max_num_changepoints,
    std::size_t kernel_approximation_dimension,
    const std::vector<std::size_t>& window_sizes,
    double linear_penalty,
    double log_linear_penalty,
    SegmenterTelemetry* telemetry = nullptr);

}  // namespace fastgatk::contamination
