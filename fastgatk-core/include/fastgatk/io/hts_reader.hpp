#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <map>
#include <string>
#include <vector>

namespace fastgatk::io {

// BAM/CRAM stores a CIGAR element as (length << 4) | operation.  Keep that
// packed representation in ReadBatch so it can be copied to a Kokkos View
// without introducing a vector of host-only objects.  The public helpers do
// not expose HTSlib types and are therefore usable by all backends.
enum class CigarOpCode : std::uint8_t {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
    ReferenceSkip = 3,
    SoftClip = 4,
    HardClip = 5,
    Padding = 6,
    SequenceMatch = 7,
    SequenceMismatch = 8,
};

struct CigarOp {
    std::uint32_t length = 0;
    CigarOpCode code = CigarOpCode::Match;

    static constexpr CigarOp unpack(std::uint32_t packed) noexcept {
        return CigarOp{packed >> 4, static_cast<CigarOpCode>(packed & 0x0fU)};
    }

    constexpr std::uint32_t pack() const noexcept {
        return (length << 4) | static_cast<std::uint8_t>(code);
    }

    constexpr bool valid() const noexcept {
        return length != 0 && static_cast<std::uint8_t>(code) <=
            static_cast<std::uint8_t>(CigarOpCode::SequenceMismatch);
    }

    constexpr bool consumes_read() const noexcept {
        return code == CigarOpCode::Match || code == CigarOpCode::Insertion ||
               code == CigarOpCode::SoftClip || code == CigarOpCode::SequenceMatch ||
               code == CigarOpCode::SequenceMismatch;
    }

    constexpr bool consumes_reference() const noexcept {
        return code == CigarOpCode::Match || code == CigarOpCode::Deletion ||
               code == CigarOpCode::ReferenceSkip || code == CigarOpCode::SequenceMatch ||
               code == CigarOpCode::SequenceMismatch;
    }

    // Only these operations map an observed read base to one reference base.
    constexpr bool projects_base() const noexcept {
        return code == CigarOpCode::Match || code == CigarOpCode::SequenceMatch ||
               code == CigarOpCode::SequenceMismatch;
    }
};

struct ReadProjection {
    std::int64_t reference_position = -1;
    std::uint32_t read_offset = 0;
    CigarOpCode operation = CigarOpCode::Match;
};

struct ReadBatch {
    // Legacy flat read payload.  These fields remain unchanged for existing
    // smoke callers and are indexed by offsets[record:record+1].
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    std::vector<std::int32_t> positions;
    std::vector<std::int32_t> tids;
    std::vector<std::uint8_t> mapq;

    // CIGAR is a second flat payload.  cigar_offsets has records()+1 entries
    // when metadata is present; cigar_ops contains packed CigarOp values.
    std::vector<std::uint32_t> cigar_offsets;
    std::vector<std::uint32_t> cigar_ops;

    // Fixed-width BAM core metadata and flat optional string tags.  Names and
    // RG values use [offsets[i], offsets[i+1]) and have no trailing NUL.
    std::vector<std::uint16_t> flags;
    std::vector<std::int32_t> mate_tids;
    std::vector<std::int32_t> mate_positions;
    std::vector<std::int32_t> template_lengths;
    std::vector<std::uint32_t> name_offsets;
    std::vector<std::uint8_t> names;
    std::vector<std::uint32_t> read_group_offsets;
    std::vector<std::uint8_t> read_groups;

    // Optional chimeric-alignment provenance used by GATK's
    // NonChimericOriginalAlignmentReadFilter. OA (original alignment) and
    // XM (mate contig) are kept as independent flat Z-tag payloads.  The
    // presence vectors are distinct from the payload offsets so a present but
    // empty tag remains distinguishable from a missing tag, matching
    // GATKRead.hasAttribute().
    std::vector<std::uint32_t> original_alignment_offsets;
    std::vector<std::uint8_t> original_alignments;
    std::vector<std::uint32_t> mate_contig_offsets;
    std::vector<std::uint8_t> mate_contigs;

    // Optional GATK/BQSR per-base indel qualities.  BI/BD are SAM Z tags
    // encoded as FASTQ characters; the Host reader stores decoded Phred bytes
    // with one offset span per record so missing tags remain aligned.
    std::vector<std::uint32_t> insertion_quality_offsets;
    std::vector<std::uint8_t> insertion_qualities;
    std::vector<std::uint32_t> deletion_quality_offsets;
    std::vector<std::uint8_t> deletion_qualities;

    // Optional flow-based read payload decoded from BAM/CRAM aux tags.  Each
    // offsets array has records()+1 entries when the corresponding metadata is
    // present; empty spans are retained so batches can be concatenated without
    // losing record alignment.  `flow_orders` stores the repeated per-flow
    // symbols from the read group's FO tag, while tp/t0 remain Host raw tag
    // values for `decode_flow_read`.
    std::vector<std::uint32_t> flow_tp_offsets;
    std::vector<std::int8_t> flow_tp;
    std::vector<std::uint32_t> flow_t0_offsets;
    std::vector<std::uint8_t> flow_t0_phred;
    std::vector<std::uint32_t> flow_order_offsets;
    std::vector<std::uint8_t> flow_orders;
    std::vector<std::uint16_t> flow_max_hmer;

    // Appended after all legacy fields to preserve aggregate-initializer and
    // source compatibility for callers that construct ReadBatch positionally.
    // These vectors distinguish a present-but-empty OA/XM tag from a missing
    // tag (GATKRead.hasAttribute()).
    std::vector<std::uint8_t> original_alignment_present;
    std::vector<std::uint8_t> mate_contig_present;
    // Optional resolved-sample ordinal for each record.  HTS decoding leaves
    // this empty for ordinary single-sample batches; a multi-sample Host
    // caller populates it before shared AssemblyRegion processing so fragment
    // cleanup and graph pruning cannot cross sample boundaries merely because
    // two inputs reuse a read name or RG identifier.
    std::vector<std::int32_t> sample_ids;

    // Optional Host-only ordinal in the outer decoded ReadBatch.  A compact
    // slice (for example an AssemblyRegion) retains this identity so later
    // annotation code can reattach a local PairHMM row to its exact original
    // BAM record, even when two records share a name and alignment span.  It
    // is deliberately not part of any Kokkos input contract.
    std::vector<std::uint32_t> source_records;

    // Optional Host-only ordinal of the outer input that supplied each
    // record.  Read names are fragment-local within one alignment input, but
    // separately supplied inputs may legitimately reuse a read name (and may
    // even be the same BAM listed twice).  Device kernels never consume this.
    std::vector<std::uint32_t> source_input_ids;

    void clear();
    std::size_t records() const;
    std::size_t bases_count() const;
    // Approximate resident bytes owned by this flat batch.  This includes
    // vector payloads and offset arrays, but intentionally excludes allocator
    // capacity and HTSlib's internal decode buffers.  It is used for runtime
    // backpressure, never as a serialization size.
    std::uint64_t bytes() const noexcept;

    // An empty cigar_offsets means a legacy batch without CIGAR metadata.
    // Such a batch deliberately retains the old contiguous-coordinate
    // projection behavior.  A non-empty, malformed layout is never accepted
    // as valid metadata.
    bool has_cigar() const noexcept;
    bool cigar_layout_valid() const noexcept;
    // Validate only one record's CIGAR span and base-consumption contract.
    // This deliberately does not inspect sibling records, so a malformed
    // record can be rejected by a per-read filter without dropping an entire
    // decoded batch.
    bool cigar_record_layout_valid(std::size_t record) const noexcept;
    bool indel_quality_layout_valid() const noexcept;
};

struct HeaderSummary {
    std::vector<std::string> contigs;
    std::vector<std::int64_t> contig_lengths;
    // Canonical SAM text for the sequence dictionary as emitted by
    // HTSJDK's SAMTextHeaderCodec when SimpleCountCollection serializes its
    // locatable metadata.  This intentionally contains only @HD (VN:1.6)
    // and @SQ records; read groups and other BAM header lines are not part of
    // a GATK HDF5SimpleCountCollection sequence dictionary.
    std::string sequence_dictionary;
    // Optional @SQ AS (assembly) tags, kept parallel to contigs.  Most BAMs
    // omit this tag; writers must therefore only render it when present.
    // Keeping the value in the Host header contract lets VCF/GVCF writers
    // preserve reference assembly metadata without exposing HTSlib objects to
    // Kokkos kernels.
    std::vector<std::string> contig_assemblies;
    // Unique sample names discovered from @RG SM tags, in header order.
    // Native writers use this to preserve the GATK/HTSJDK sample column.
    std::vector<std::string> samples;
    // Read-group to sample mapping from @RG ID/SM.  Keeping this mapping in
    // the Host header contract lets walkers (for example DepthOfCoverage)
    // partition a multi-sample BAM without pulling HTSlib objects into their
    // kernels.  A missing RG/SM is deliberately absent and is handled by the
    // caller's explicit unknown-sample policy.
    std::map<std::string, std::string> read_group_samples;
};

// Normalized Host interval used by BAM/CRAM readers. Coordinates are zero-based
// half-open, matching HTSlib and the rest of the native traversal code.
struct HtsInterval {
    std::int32_t tid = -1;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// GATK interval arguments support a union (the default) or an intersection
// across repeatable -L/--intervals selectors. Keep the rule in the shared Host
// reader so indexed and sequential traversal expose identical coordinates to
// every caller.
enum class HtsIntervalSetRule : std::uint8_t {
    Union = 0,
    Intersection = 1,
};

class HtsReader {
public:
    HtsReader(const std::string& path,
              const std::string& reference = {},
              const std::string& region = {},
              std::size_t batch_records = 1024,
              const std::string& sample = {},
              bool include_unmapped = false);
    HtsReader(const std::string& path,
              const std::string& reference,
              const std::vector<std::string>& regions,
              std::size_t batch_records = 1024,
              const std::string& sample = {},
              const std::vector<std::string>& exclusions = {},
              std::int64_t interval_padding = 0,
              std::int64_t exclusion_padding = 0,
              HtsIntervalSetRule interval_set_rule = HtsIntervalSetRule::Union,
              bool include_unmapped = false);
    ~HtsReader();

    HtsReader(HtsReader&&) noexcept;
    HtsReader& operator=(HtsReader&&) noexcept;
    HtsReader(const HtsReader&) = delete;
    HtsReader& operator=(const HtsReader&) = delete;

    bool next(ReadBatch& batch);
    // Change the decode batch at a safe record boundary.  The next call to
    // next() observes this value; an in-flight batch is never split.  Native
    // tools use this hook with AdaptiveController so byte pressure can reduce
    // staging without changing input order.
    void set_batch_records(std::size_t batch_records);
    std::size_t batch_records() const noexcept;
    bool compiled() const;
    // True when HTSlib found an index and interval traversal uses direct
    // indexed iteration.  Without an index the reader remains correct but
    // falls back to sequential decode-and-filter semantics.
    bool indexed() const noexcept;
    // True when HTSlib can load an index for this input, even when this reader
    // was opened without intervals and therefore remains a sequential reader.
    // Region-streaming callers use this capability check before constructing
    // their per-tile interval readers.
    bool has_index() const noexcept;
    // Number of records decoded from the underlying HTSlib stream, including
    // records rejected by interval/sample selection. This is Host telemetry
    // and does not change the filtered ReadBatch contract.
    std::size_t records_read() const noexcept;
    const std::string& backend_description() const;
    const HeaderSummary& header() const;
    // Raw selectors after parsing, before the reader's traversal union is
    // normalized.  Most walkers use intervals() for efficient indexed
    // traversal; tools that must reproduce GATK's requested interval-merging
    // rule use this view to preserve overlap/adjacency provenance.
    const std::vector<HtsInterval>& requested_intervals() const;
    const std::vector<HtsInterval>& intervals() const;
    // Normalized exclusion intervals.  Read traversal remains include-based;
    // callers apply this mask at the projected-locus boundary so a read that
    // spans an excluded and an included locus is not discarded wholesale.
    const std::vector<HtsInterval>& exclusion_intervals() const;
    std::size_t interval_file_inputs() const;
    std::size_t interval_file_records() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Decode and safely project a read base to a reference coordinate.  Returns
// false for insertions/soft clips, unmapped reads, out-of-range offsets, or an
// invalid CIGAR.  For legacy batches without CIGAR metadata, the historical
// contiguous mapping position + read_offset is used.
bool project_read_offset(const ReadBatch& batch, std::size_t record,
                         std::size_t read_offset, ReadProjection& projection) noexcept;

// Return the exclusive reference end of a read, or -1 for an invalid/unmapped
// record.  CIGAR-aware callers should use this for interval overlap instead of
// position + read length.
std::int64_t reference_end(const ReadBatch& batch, std::size_t record) noexcept;

}  // namespace fastgatk::io
