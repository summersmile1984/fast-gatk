#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

struct ReadErrorCorrectionInput {
    std::vector<std::uint8_t> bases;
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint8_t> qualities;
};

struct ReadErrorCorrectionOptions {
    std::uint32_t kmer_length = 25;
    std::uint32_t min_observations_for_kmer_to_be_solid = 20;
    // NearbyKmerErrorCorrector's production defaults search solid neighbors
    // up to two Hamming mismatches.  The field remains explicit so callers can
    // tighten the policy for low-coverage or highly repetitive regions.
    std::uint32_t max_mismatches = 2;
    std::uint8_t corrected_base_quality = 30;
    // Only very sparse kmers are eligible for correction.  This matches the
    // Java corrector's default maxObservationsForKmerToBeCorrectable=1 and
    // prevents a real low-frequency allele from being rewritten wholesale.
    std::uint32_t max_observations_for_kmer_to_be_correctable = 1;
};

struct ReadErrorCorrectionResult {
    bool used = false;
    std::size_t input_reads = 0;
    std::size_t solid_kmers = 0;
    std::size_t corrected_kmers = 0;
    std::size_t uncorrectable_kmers = 0;
    std::size_t corrected_reads = 0;
    std::size_t corrected_bases = 0;
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Build a deterministic exact-k-mer table on Host, map sparse kmers to the
// nearest solid Hamming neighbor (up to max_mismatches), and collect the
// resulting per-base edits with strict consensus across overlapping windows.
// Host owns the irregular correction map; Kokkos applies the final immutable
// edit mask and quality updates so all execution spaces share one contract.
ReadErrorCorrectionResult correct_read_errors_kokkos(
    const ReadErrorCorrectionInput& input,
    ReadErrorCorrectionOptions options = {});

}  // namespace fastgatk::kernels
