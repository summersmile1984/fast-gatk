#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::somatic {

// Host-side representation of one standard GATK CollectF1R2Counts sample.
// The arrays are intentionally plain containers: Kokkos kernels produce the
// counts, while this boundary owns the HTSJDK metrics/tar serialization.
struct F1R2ArchiveAltRow {
    std::string context;
    int ref_count = 0;
    int alt_count = 0;
    int ref_f1r2 = 0;
    int alt_f1r2 = 0;
    int depth = 0;
    char alt = 'N';
};

struct F1R2ArchiveSample {
    std::string sample;
    std::array<std::vector<std::uint64_t>, 64> ref_hist{};
    std::array<std::vector<std::uint64_t>, 384> alt_hist{};
    std::vector<F1R2ArchiveAltRow> alt_rows;
    int max_depth = 200;
};

// Write the standard three-member CollectF1R2Counts archive consumed by
// LearnReadOrientationModel.  This function does not perform pileup or
// filtering; callers provide deterministic, already-filtered counts.
void write_f1r2_archive(const std::string& output,
                        const std::vector<F1R2ArchiveSample>& samples);

}  // namespace fastgatk::somatic
