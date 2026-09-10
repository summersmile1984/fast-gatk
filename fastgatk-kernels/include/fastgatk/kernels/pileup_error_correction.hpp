#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// A flat, HTSlib-free representation of the observations needed by GATK's
// PileupReadErrorCorrector.  cigar_ops uses BAM's packed (length << 4)|opcode
// encoding; opcode values are the SAM CIGAR operation numbers.
struct PileupReadErrorCorrectionInput {
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    std::vector<std::int32_t> tids;
    std::vector<std::int32_t> positions;
    std::vector<std::uint32_t> cigar_offsets;
    std::vector<std::uint32_t> cigar_ops;
};

struct PileupReadErrorCorrectionOptions {
    // GATK exposes this as a hidden argument and defaults it to -infinity,
    // which disables pileup correction.  A finite value enables the path.
    double log_odds_threshold = -std::numeric_limits<double>::infinity();
    std::uint8_t corrected_base_quality = 30;
    std::uint32_t indel_span = 15;
    std::uint32_t indel_mismatches = 3;
};

struct PileupReadErrorCorrectionResult {
    bool used = false;
    std::size_t input_reads = 0;
    std::size_t loci = 0;
    std::size_t corrected_loci = 0;
    std::size_t skipped_indel_adjacent_bases = 0;
    std::size_t corrected_reads = 0;
    std::size_t corrected_bases = 0;
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Reproduce the GATK pileup correction policy.  Irregular locus/CIGAR work is
// deterministic Host preprocessing; the final base/quality edit is performed
// by a Kokkos kernel shared by all execution spaces.
PileupReadErrorCorrectionResult correct_reads_by_pileup_kokkos(
    const PileupReadErrorCorrectionInput& input,
    PileupReadErrorCorrectionOptions options = {});

}  // namespace fastgatk::kernels
