# fastgatk-runtime

Read-only resource boundary shared by native tools. `ResourceSnapshot::probe()`
reads the process CPU affinity, cgroup v2 memory limit/current usage, SLURM CPU
and memory allocation, scratch directory, job id, and launcher-visible
CUDA/HIP/SYCL device selectors.
When available it also reports scratch free inodes and the process open-file
limit from the Linux allocation boundary.
The launcher may additionally provide `FASTGATK_DEVICE_MEMORY_BYTES`,
`FASTGATK_DEVICE_FREE_BYTES`, `FASTGATK_LOCAL_SSD`, and
`FASTGATK_REMOTE_INPUT`; these values are surfaced as explicit hard/free
budget fields rather than inferred from a path. The JSON separates
`reported_device_*` (launcher provenance), `device_backend_compiled`, and
`device_assignment_available` from the controller-usable `device_*` budget.
Thus a CUDA-visible allocation cannot make an OpenMP/Serial binary claim a GPU
budget: it is surfaced as `device_telemetry_rejected:true` and usable device
fields remain zero. Unknown device limits remain zero and never become an
unbounded allocation.
`effective_threads(requested)` caps a requested Kokkos thread count to the
visible allocation, while `effective_memory_limit_bytes()` exposes the
conservative cgroup/SLURM limit. `safe_memory_budget_bytes()` reserves 20%
headroom for runtime/HTSlib overhead; native readers use it to reject
unbounded staging with `RESOURCE_EXHAUSTED` before entering a kernel.

`AdaptiveController` derives an initial batch limit at 60% of known host,
device and in-flight budgets, halves limits under pressure, and grows an idle
queue in 20% steps. `BoundedByteQueue<T>` provides blocking byte-capacity
backpressure for decoded reads, compute batches and encoded results. The native
HC reader applies the controller's initial safe batch before HTSlib decoding;
actual staging remains checked on every batch.

`ThreeStagePipeline<Decoded, Computed, Encoded>` composes those queues into a
single decode -> compute -> encode -> sink lifecycle. Decode, compute and encode
callbacks run on separate Host threads, while the sink is drained by the caller
thread to preserve deterministic output order. Each queue is byte-bounded,
stage exceptions close all queues and are rethrown after a clean join, and
metrics expose item/byte counts plus per-stage peak queue occupancy. The
`CountReads`, `FlagStat`, `CollectReadCounts`, `GetPileupSummaries`,
`CollectAllelicCounts`, `DepthOfCoverage`, and `CollectF1R2Counts` native
tools use this executor;
their manifests record the pipeline counters alongside the Kokkos `KernelPlan`
telemetry.

The native CTest `fastgatk-resource-limits-contract` replays this boundary in
both Kokkos host backends. It injects the active `fastgatk-hc-call` target via
`FASTGATK_RESOURCE_LIMITS_BINARY`, caps an eight-thread/requested-4096 batch to
the simulated `SLURM_CPUS_PER_TASK=2`/`SLURM_MEM_PER_NODE=1M` allocation, and
requires a tiny `1K` allocation to fail closed with `RESOURCE_EXHAUSTED`.
The same script remains runnable directly, with the OpenMP binary as its
default or an explicit target path in that environment variable.

## Output publication contract

`fastgatk/runtime/output.hpp` provides the reusable `OutputBundle` boundary for
tools that produce a primary file, an optional/required index, and an
`OutputManifest`. `validate_output_bundle()` is read-only: it requires regular,
non-empty artifacts, a bounded JSON-object manifest with `schema_version`,
`outputs`, and `complete:true` entries, and (by default) rejects conventional
`.tmp`, `.partial`, and `.incomplete` sibling residue. `require_complete_output()`
converts any incomplete state into the stable `OUTPUT_CONTRACT_FAILURE` class.

`publish_output_bundle(staging, destination)` validates the staging set before
renaming anything, rejects pre-existing destinations, renames the manifest last
as the completion marker, and rolls already-renamed files back to staging if a
rename or final validation fails. Staging and destination paths must be on the
same filesystem. Existing native tools are unchanged; callers can adopt this
API at their output boundary incrementally. `fastgatk-runtime-output-smoke`
covers complete success, incomplete `outputs[]`, missing-index failure,
temporary-residue failure, and failed-commit rollback in both OpenMP and Serial
builds.

Every new algorithm follows the same contract rather than adding a backend-
specific implementation:

1. Host decodes BAM/CRAM/VCF/FASTA and owns strings, ordering, CIGAR parsing,
   interval semantics, file formats, spill/checkpoint and fallback decisions.
2. A flat, bounded `HostBatch` is validated and passed to
   `KernelPlan::prepare`; device data is represented only by Kokkos `View`s.
3. Numeric work is launched through a Kokkos execution-space policy
   (`RangePolicy`, `TeamPolicy` or `MDRangePolicy`). CPU vectorization uses
   `Kokkos::Experimental::simd`; production sources must not include
   `immintrin.h` or call `_mm*` intrinsics.
4. `execute` fences at the contract boundary, then Host `collect` restores
   deterministic order and performs the format/annotation write.  Strict mode
   fixes traversal, reduction order, seeds and floating-point policy; Fast mode
   reports its tolerance and is never silently presented as bit-identical.

The `fastgatk-kokkos-api-boundary` CTest enforces the no-raw-intrinsics rule
for production `fastgatk-kernels` and native tool sources.  A tool is promoted
from `prototype` to `native-compatible` only after this lifecycle, its own
GATK oracle, format/index/sidecar checks, resource-pressure test and
Nextflow/SLURM path all pass.

`ReadBatch::bytes()` accounts the flat payload (including CIGAR, tags and
offsets), and `HtsReader::set_batch_records()` changes the decode size only at a
safe batch boundary. CollectReadCounts and GetPileupSummaries use these hooks
with `AdaptiveController::next()` to reduce the next batch under observed host
byte pressure while preserving read order. Their OutputManifests report the
initial/final batch sizes and reduction count.

The probe does not submit jobs, create a device context, or change scheduler
state. It only consumes launcher-provided device telemetry; a CUDA/HIP/SYCL-
specific discovery plugin can populate the same environment fields without
changing the C++ runtime ABI. A real device-context/throughput check remains a
GPU-runner responsibility and is never inferred from these fields.
Spill/retry remains an explicit caller policy, and an unknown limit is never
converted into a host-wide allocation.
