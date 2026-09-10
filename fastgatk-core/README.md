# fastgatk-core

`fastgatk-core` is the common C++ library boundary used by the Kokkos kernels and
native tools. It does not own GATK tool semantics or HTSlib; it standardizes the
batch/plan lifetime that the runtime uses to amortize preparation.

The public pieces are:

- `fastgatk::core::HostBatch`: schema, logical record count and host byte accounting;
- `fastgatk::core::DeviceBatch<ExecSpace>`: typed Kokkos View binding and device-byte accounting;
  `bind(label, view, logical_span)` accounts a logical batch prefix when a persistent
  capacity view is reused;
- `fastgatk::core::KernelPlan<ExecSpace>`: prepare/execute lifetime and telemetry.
- `fastgatk::io::ReadBatch`: flat read arrays plus packed CIGAR offsets/operations and BAM
  core metadata; `project_read_offset` and `reference_end` provide bounds-checked reference
  projection without exposing HTSlib types to kernels.
- `fastgatk::io::decode_flow_read`: Host-only decoder for GATK FlowBasedRead `tp`/`t0`
  tags and read-group flow order. It emits the flat `[flow][256]` calibrated table consumed
  by the independent flow-space PairHMM; tag parsing and boundary-flow policy stay out of
  Kokkos kernels.

`ReadBatch` keeps the original `offsets/bases/qualities/positions/tids/mapq` fields for
existing smoke callers. A batch with an empty `cigar_offsets` uses the legacy contiguous
projection; HTSlib-backed batches always populate CIGAR and metadata offsets.
`ReadBatch::bytes()` is the single saturating byte-accounting boundary for all flat
payloads (including CIGAR, names/RG, mate metadata, BI/BD and flow tags).  An
`HtsReader::set_batch_records()` update is observed only by the next `next()` call, so
adaptive controllers can change decode pressure without reordering records or mutating
an in-flight batch.

When a normalized interval is supplied and a BAM/CRAM index is available,
`HtsReader::indexed()` reports true and HTSlib iterates the index directly.  Without
an index it deliberately falls back to sequential decode-and-filter semantics; callers
that require fixed-memory region streaming must fail closed rather than silently
rescanning the whole input for every tile.

Module-specific plans still own strongly typed Views. The common plan deliberately
does not erase those types or expose raw device pointers.
