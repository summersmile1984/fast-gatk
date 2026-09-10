#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace fastgatk::kernels {

struct ReadFilterInput {
    std::vector<std::uint32_t> offsets;
    std::vector<std::int32_t> tids;
    std::vector<std::int32_t> positions;
    std::vector<std::uint8_t> mapq;
    std::vector<std::uint16_t> flags;
    // Optional packed BAM CIGAR metadata. Empty vectors retain the legacy
    // filter contract; CIGAR-aware filters fail closed when requested but the
    // metadata is absent or malformed.
    std::vector<std::uint32_t> cigar_offsets;
    std::vector<std::uint32_t> cigar_ops;
    // Optional per-record RG and OA/XM auxiliary-tag payloads. Wellformed
    // filtering requires a non-empty RG span; NonChimericOriginalAlignment
    // uses the OA/XM spans. Each offset array is records()+1 when present.
    std::vector<std::uint32_t> read_group_offsets;
    std::vector<std::uint8_t> read_groups;
    std::vector<std::uint32_t> original_alignment_offsets;
    std::vector<std::uint8_t> original_alignments;
    std::vector<std::uint32_t> mate_contig_offsets;
    std::vector<std::uint8_t> mate_contigs;
    // Appended after the legacy payload fields to preserve aggregate-
    // initializer/source compatibility for callers that provide OA/XM
    // offsets and payloads positionally.
    std::vector<std::uint8_t> original_alignment_present;
    std::vector<std::uint8_t> mate_contig_present;
};

struct ReadFilterOptions {
    std::uint8_t min_mapq = 0;
    bool exclude_unmapped = true;
    bool exclude_secondary = true;
    bool exclude_supplementary = true;
    bool exclude_duplicates = false;
    bool exclude_qcfail = true;
    std::uint32_t max_reads_per_locus = 0;
    std::uint64_t seed = 0x9e3779b97f4a7c15ULL;
    // Mutect2 replaces AssemblyRegionWalker's ordinary positional sampler
    // with MutectDownsampler.  Once a stride exceeds its coverage budget it
    // reservoirs only reads whose MAPQ is strictly greater than 50.  This is
    // Host-side traversal policy; the Kokkos mask remains the common fixed
    // read-filter implementation for all callers.
    bool mutect2_downsampling = false;
    // M2ArgumentCollection controls the MutectDownsampler pool width and its
    // optional all-stride rejection threshold.  A non-positive suspicious
    // threshold disables that rejection, as in GATK.
    std::uint32_t mutect2_downsampling_stride = 1;
    std::int32_t mutect2_max_suspicious_reads_per_alignment_start = 0;
    bool require_good_cigar = false;
    bool require_nonzero_reference_span = false;
    bool require_non_chimeric_original_alignment = false;
    // Parameterized GATK ReadLengthReadFilter.  The class is enabled only
    // when require_read_length is true; keeping the bounds in the typed
    // options lets the same Kokkos mask run on Serial/OpenMP/CUDA/HIP/SYCL.
    bool require_read_length = false;
    std::uint32_t min_read_length = 1;
    std::uint32_t max_read_length = std::numeric_limits<std::uint32_t>::max();
    // GATK MappingQualityAvailableReadFilter rejects the sentinel MAPQ=255;
    // MappingQualityNotZeroReadFilter rejects MAPQ=0. Keep these trailing
    // for aggregate-initializer source compatibility.
    // Keep this trailing for aggregate-initializer source compatibility.
    bool exclude_mapping_quality_unavailable = false;
    bool exclude_mapping_quality_zero = false;
    // GATK WellformedReadFilter rejects reference-skip (N) CIGAR operators
    // and records without a read-group tag. Keep these trailing for source
    // compatibility with older aggregate initializers.
    bool require_no_n_cigar = false;
    bool require_read_group = false;
};

struct ReadFilterResult {
    bool used = false;
    std::size_t input_reads = 0;
    std::size_t passed_reads = 0;
    std::size_t filtered_reads = 0;
    std::size_t downsampled_reads = 0;
    std::vector<std::uint8_t> keep;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Apply the fixed-order read filters on Kokkos and perform optional stable
// per-locus downsampling on Host.  The returned mask is ordinal-stable and can
// be used by HC, BQSR and somatic tools without carrying HTSlib objects onto a
// device.
ReadFilterResult filter_reads_kokkos(
    const ReadFilterInput& input, ReadFilterOptions options = {});

}  // namespace fastgatk::kernels
