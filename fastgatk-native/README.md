# Native HTSlib / HaplotypeCaller smoke adapter

`AnnotateIntervals` 现接受 GATK 4.6.2.0 的 `--feature-query-lookahead`（含 0/负值），并
将其写入 OutputManifest/telemetry；pinned Java/native interval oracle 在 OpenMP/Serial
均通过，annotation 数值保持不变。云端 FeatureManager track/query 语义仍按 registry
显式 fallback。

`VariantRecalibrator` 现按 GATK `VariantDataManager` 在归一化后的 training/non-training 均值偏移
稳定重排 annotation dimensions，再进入 Kokkos VBEM；该顺序影响多高斯随机协方差初始化。
Pinned GATK 4.6.2.0 annotation-order oracle 在 OpenMP/Serial 均通过（VQSLOD 最大差约 4e-5）。
双高斯 VBEM oracle 另验证 Normal-Wishart mean outer-product 只乘一次 shrinkage factor：
PMix/mean/covariance 最大差 `5.23e-9`、VQSLOD 最大差 `5e-5`；修复前为 `1.22744`。
完整 raw-bit model convergence 与 resource calibration 仍保持显式 fallback。

`BaseRecalibrator` 现支持长读 CycleCovariate 的 Java 边界：默认 absolute cycle 上限为
500，超限读段与 GATK 一样 fail-closed；`--maximum-cycle-value`/`--max-cycle` 可显式
选择更大域，并写入 checkpoint、GATKReport Arguments 和 manifest。CIGAR insertion 的
read bases 也按 GATK `isSNP=0` 纳入 substitution covariates。pinned GATK 4.6.2.0
1,500M + insertion synthetic 在 OpenMP/Serial 比较 Quantized、RecalTable0/1/2 均
bit-identical；完整云端分布式 BQSR 仍按 registry 保持显式 fallback。

`ModelSegments` 的多样本 KernelSegmenter 现在按完整 sample vector 计算（线性 kernel
拼接 feature、Gaussian kernel 求和各样本 anchor），与 GATK 4.6.2.0 对齐并避免反相关
事件被均值归约抵消；双 OpenMP/Serial pinned oracle 已通过。完整多样本 posterior/SVD
raw-bit 一致性仍保持显式 fallback。

HC/Mutect2 的 PairHMM Host preparation 也已修正 GATK 无 `BI/BD` 标签时的质量：raw
默认 Q45，Conservative PCR 以原始 read 坐标仅降低非终位；clipping 后末位若原本是
内部重复位点可为 Q20/Q40，`NONE` 模式保持 Q45。706-row pinned GATK 4.6.2.0 输入
oracle 在 OpenMP/Serial 双后端验证，显式 BI/BD 标签不受影响。

`LearnReadOrientationModel` 现接受 GATK 4.6.2.0 的 `--QUIET/--quiet`、`--tmp-dir`、
`--verbosity/--VERBOSITY` 与 JDK deflater/inflater Boolean aliases，并在 manifest telemetry
记录选择。两样本 `.orientation_priors` tar writer pinned Java/native oracle 在 OpenMP/
Serial 均通过；member/gzip bytes 仍不宣称 bit-identical。

`CollectReadCounts` 现支持 GATK 4.6.2.0 的 `--QUIET`、`--read-validation-stringency/-VS`、
`--disable-bam-index-caching/-DBIC`、`--tmp-dir` 与 JDK deflater/inflater Boolean/value
controls，并在 manifest telemetry 记录选择。TSV interval/exclusion pinned Java/native
contract+oracle 在 OpenMP/Serial 均通过；HDF5/cloud 仍按 registry 显式 fallback。

`GatherBQSRReports` 现接受 GATK 4.6.2.0 的 `--QUIET/--quiet`、`--tmp-dir`、
`--verbosity` 与 JDK deflater/inflater Boolean/value controls，并在 manifest telemetry
记录选择。pinned 五表 Java/native oracle 在这些 CLI 参数下保持 bit-identical；空报告和
不兼容 report dimensions 仍 fail-closed。

`SortSam` 现支持 Picard/GATK 4.6.2.0 的大写 CLI aliases（`--INPUT/--OUTPUT`、
`--REFERENCE_SEQUENCE`、`-SO/--SORT_ORDER`、`--MAX_RECORDS_IN_RAM/--TMP_DIR`）、
`--CREATE_INDEX`、`--QUIET`、JDK codec switches、`--COMPRESSION_LEVEL`、
`--VALIDATION_STRINGENCY` 与 `--VERBOSITY`。pinned Java/native oracle 在 OpenMP/Serial
验证坐标排序记录/header、bounded spill 以及 manifest 的 compression/validation telemetry；
完整 Picard sort orders、cloud/MD5 与 byte-identical codec 输出仍按 registry 显式 fallback。

`GatherVcfs` 现接受 GATK/Picard 的大写 `--INPUT/--OUTPUT/--REFERENCE_SEQUENCE`、
`--QUIET`、`--VERBOSITY` 别名，并支持 comment、compression-level 与 `CREATE_INDEX`
writer controls；pinned GATK 4.6.2.0 CLI oracle 已在 OpenMP/Serial 通过。完整 cloud、
TMP 和 MD5 语义仍由 registry 显式 fallback。

`LeftAlignAndTrimVariants` 现支持 GATK 4.6.2.0 的 `-XL/--exclude-intervals`、
`-ip/--interval-padding`、`-ixp/--interval-exclusion-padding`，并按记录 END span
执行区间筛选；`-isr/-imr/-OVI` aliases 及 trim/split/keep-original 的 separated
Boolean 也可直接替换 GATK 调用。pinned Java/native CLI oracle 已在 OpenMP/Serial
通过，manifest 会记录 selector、padding、skip 与 index policy。复杂 symbolic/cloud/
bit-identical normalization 仍是 explicit fallback。

GenotypeGVCFs now accepts GATK's repeatable `-XL/--exclude-intervals`. The Host
traversal first constructs the `-L` include set and then removes records whose
gVCF span intersects the exclusion union; aggregate and `--stream-by-locus`
share this predicate. A pinned GATK 4.6.2.0 oracle compares concrete-variant
rows and OutputManifest telemetry on both OpenMP and Serial builds.

GenotypeGVCFs also implements GATK 4.6.2.0's deprecated
`--only-output-calls-starting-in-intervals` Boolean as a final writer boundary.
Input traversal remains span-overlap based, so a deletion beginning before
`-L` still participates in joint genotyping; only the finalized call's POS is
tested against the normalized include intervals. A self-contained pinned Java
oracle proves default overlap output versus STARTS_IN output and compares exact
Java/native rows for aggregate and `--stream-by-locus` on OpenMP and Serial.
Manifest telemetry exposes the selected mode and skipped/output record counts.

`LeftAlignAndTrimVariants` 支持 GATK 4.6.2.0 的 `--sites-only-vcf-output` optional
Boolean：native 在完成 normalization、multiallelic split 及 genotype/annotation
投影后，于独立 writer header 边界用 `bcf_subset_format` 移除 FORMAT/sample，输出严格
8 列 site-only VCF；关闭时保留完整样本字段。pinned Java/native oracle 覆盖 bare、
separated true/false，OpenMP/Serial 均通过。

`ValidateVariants` 的 ALT usage check 也与 GATK 4.6.2.0 对齐：只有 concrete ALT
必须出现在 called genotype 中，`<DEL>`、`<CNV>` 等 symbolic ALT 即使未被调用也不
会误报。pinned symbolic/mixed/concrete-failure oracle 及新增 CTest 在 OpenMP/Serial
均通过；完整 htsjdk/Tribble validation 仍按 registry 保持显式 fallback。

`VariantFiltration` 还支持 `--set-filtered-genotype-to-no-call`：带非-PASS FORMAT/FT
的 called genotype 会变为 unphased `./.`，只有实际发生转换时才重算 AC/AN/AF；GATK
4.6.2.0 的 bare/true/false、existing FT 和 embedded `=` fail-closed 边界由 pinned
oracle 与 OpenMP/Serial CTest 覆盖。

`fastgatk-hc-smoke` is the first real file-boundary native adapter. It reads SAM/BAM/CRAM
with HTSlib, batches decoded reads into flat host arrays, copies them into Kokkos Views,
and runs a deterministic base/quality/mapQ reduction kernel. It accepts the GATK
HaplotypeCaller input/output aliases (`-I/-R/-L/-O`) and common launcher options, emits
a valid VCF container plus `OutputManifest` and telemetry JSON, and is designed to be
called from an unchanged Nextflow process or SLURM allocation. Its VCF is intentionally
empty: this is still an adapter/smoke boundary, not a biological variant caller.

The native build accepts an optional FASTGATK_KOKKOS_ARCH environment variable
and forwards it to Kokkos (for example ZEN3/ZEN4 on AMD Zen hosts or SKX/ICL on
Intel AVX-512 hosts). This changes only Kokkos' compile-time ISA/execution
policy; all callers continue to use the same Kokkos APIs. Use separate build
directories or containers for portable scalar and ISA-tuned builds, and never
enable an ISA unavailable on the target CPU.
`FASTGATK_KOKKOS_BACKEND=cpu|serial|cuda|hip|sycl` (plus explicit
`*_openmp` variants) selects the Kokkos execution-space artifact in
`scripts/build_native.sh`. Device builds require the corresponding compiler
and runtime on the worker; the full correctness gate runs on CPU/OpenMP and a
serial-only artifact is compile-verified. CUDA/HIP/SYCL remain compiler/runtime
dependent and are not claimed as executed here.
For direct CMake use, `-DFASTGATK_KOKKOS_BACKEND=SERIAL|OPENMP|CUDA|HIP|SYCL`
selects one backend before Kokkos is configured; `CUDA_OPENMP`, `HIP_OPENMP`,
and `SYCL_OPENMP` select a device plus an OpenMP host execution space
(`AUTO` preserves toolchain cache settings). Device-only aliases retain
Kokkos Serial as the portable Host decode/collect space, while the device
remains the default execution space for numeric kernels. The optional
`-DFASTGATK_KOKKOS_ARCH=ZEN4|SKX|AMPERE80|...` alias forwards architecture
selection to Kokkos, so the same source does not grow ISA-specific branches.
This is a build-time selection only: production kernels keep the same Kokkos
API and no backend-specific algorithm branch is introduced.

Every configured build also exposes the `fastgatk-kokkos-backend-config` CTest.
It validates the CMake cache, generated `KokkosConfigCommon.cmake`, selected
execution space, architecture flags, and compiler provenance without creating
a GPU context. Therefore CPU-only CI can catch a stale/mixed backend cache;
actual CUDA/HIP/SYCL device availability remains a separate runtime gate and
must fail closed to the documented fallback when the worker has no device.

The PairHMM compatibility API exposes `compiled_simd_backend()` and
`compiled_simd_width()` so a launcher can audit the ABI embedded in the
binary. `compute_simd(..., Avx2/Avx512)` rejects a request that does not match
that compile-time Kokkos ABI (and still checks host CPU flags); it never
silently executes an AVX-512-labelled request through an AVX2 or scalar build.
The kernel benchmark and API smoke persist this ABI/CPU contract in JSON.

## Architecture baseline (locked)

The implementation baseline is now fixed for all subsequent algorithms:

* C++ is the only native Host language; HTSlib, file formats, string/object
  semantics, ordering, traceback, and fallback state machines stay on Host.
* Kokkos is the only production kernel API. Kernels consume flat, typed Views
  through `KernelPlan.prepare/execute/collect`; CPU SIMD uses Kokkos SIMD and
  the same source selects Serial/OpenMP/CUDA/HIP/SYCL execution spaces.
* Every tool must preserve the GATK command/output/index/sidecar contract,
  expose deterministic Strict versus explicitly-tolerant Fast behavior, and
  record backend, resource, timing, and fallback state in `OutputManifest`.
* Unsupported semantics fail closed to the registered Java/Picard/Spark
  fallback. A faster result cannot promote a tool past an unmet oracle or
  output/resource contract.

The same rule is observable inside a native tool: HC telemetry and its
`OutputManifest.compatibility` record `pairhmm_fallback_reason` (for example
`executed`, `no-reference`, `no-candidates`, or
`no-read-haplotype-overlap`).  A pileup result therefore cannot be mistaken for
a PairHMM result merely because both satisfy the VCF container contract.

HC's file-boundary paths also use the shared runtime `OutputBundle` validator
after writing the manifest. A primary output, requested index, and complete
manifest must all be present before success is reported; stdout remains the
intentional stream-only exception.

This architecture is validated as an implementation gate, not as a claim that
all GATK algorithms are already complete. The active CPU/OpenMP and Serial CTest
inventory contains 209 tests in each configured build, including CountReads, FlagStat, SplitIntervals,
CountBasesInReference, CompareReferences,
CheckReferenceCompatibility and GatherTranches/GatherVcfs
merge contracts, including the bundled GATK 4.6.2.0 Oracle tests).
On 2026-09-01, the 194-test matrix completed continuously with 194/194 passes
on both OpenMP (1198.16 s) and Serial (1184.40 s) builds; the current inventory
is 209 after subsequent bounded additions, including the GenotypeGVCFs
STARTS_IN writer, BQSR context-size, SortSam duplicate-comparator, Mutect2 NLOD,
PairHMM fragment-aggregation, and GenomicsDB native-storage-boundary oracles; these gates passed independently on both backends. This is
a regression result, not a claim of complete GATK algorithm parity.
The current bounded oracles also cover ValidateVariants Boolean/fail-closed behavior,
GenomicsDBImport TileDB workspace publication/reopen plus the legacy bridge,
LeftAlignAndTrimVariants genotype/annotation projection, VariantEval CountVariants
reports, and Java-readable HDF5 SimpleCountCollection metadata; these remain scoped
slices rather than claims of complete assembly, evaluator, or CNV model parity.
The PairHMM fragment-first path now preserves Java's `-Infinity` zero-probability
mate evidence through ordered Kokkos fragment×haplotype aggregation; the pinned
six-cell oracle is raw-bit exact on OpenMP and Serial. ModelSegments manifests
also distinguish GATK total/burn-in/retained iteration semantics from the native
bounded conditional-point estimate, and explicitly mark full Gibbs/MCMC and Java
RNG identity as unavailable. Native GenomicsDB workspaces are deliberately
declared `fastgatk-portable-sparse-v1` (not TileDB and not query-compatible);
real TileDB publication/query remains the external adapter/bridge path.
`fastgatk-progress-score-contract` validates the weighted snapshot and registry
arithmetic, `fastgatk-fixture-digest-contract` binds the oracle inputs to exact
SHA-256 bytes, and `fastgatk-resource-limits-contract` exercises SLURM-aware
adaptive batching plus fail-closed `RESOURCE_EXHAUSTED` behavior against the
active `fastgatk-hc-call` target. The historical pre-gate full OpenMP matrix passed 84/84;
the new resource gate passes 1/1 on both OpenMP and Serial, and the default
Serial focused matrix now includes it (54 tests), including GatherBQSRReports,
FilterMutectCalls, CollectReadCounts, scheduler retry, VariantFiltration Java, and
MarkDuplicates/Picard oracles. The ReblockGVCF real GATK oracle
passes 2/2 on both OpenMP and Serial; its supported non-overlapping same-sample
shard oracle passes 1/1 on both backends. The targeted progress-score/
CountBases pair is 2/2 on each backend and fixture digest is 1/1. Set
`FASTGATK_VERIFY_SERIAL_FULL=1` to add the complete current 198-test Serial release
gate. The latest Mutect2/SomaticLikelihoods changes were rechecked
with the targeted HC/Mutect2/API gate on both backends. CTest injects the active
`FASTGATK_NATIVE_BUILD` and expected execution-space label, so a portability
run cannot silently execute a sibling OpenMP binary. This is a reusable
architecture gate, while
full HC assembly semantics, release-specific Mutect2 posterior, native
TileDB/GenomicsDB import/storage and high-throughput query (the current bridge
covers local workspace export), complete JEXL/Picard/cloud behavior, and GPU
runtime CI remain explicit follow-up work.

The current portability/compatibility increment also has independent evidence:
the BQSR Java oracle compares all four report tables and 493 Java-`PrintReads`
records (including the `OQ`/`--use-original-qualities` path), and the HC
reverse-strand indel oracle checks reverse-complement, anchored CIGAR and
left-normalization behavior. The dispatcher contract covers Java launcher
options, properties-file fallback, recursive `@args`, `--` separation and
fail-closed unknown options. Runtime gates cover cgroup/SLURM hard and free
host budgets, visible-device limits, scratch/local-SSD and remote-input
telemetry, byte-bounded backpressure and `RESOURCE_EXHAUSTED`; remote inputs
are staged atomically before native execution. A shared atomic Tribble writer
publishes ordinary VCF `.idx` sidecars only after the primary output is valid.
The SelectVariants Java 4.6.2.0 semantic oracle checks type/sample/JEXL/
concordance/allele-subset records and stable GT/AD/PL/GQ/AC/AN/AF payloads.
`fastgatk-runtime-output-smoke` checks complete primary/index/manifest
publication, missing-index and temporary-residue fail-closed behavior, and
atomic rollback. The focused 52-test slice for these boundaries passed on
both OpenMP and Serial; the latest 125-test release matrix passes 125/125
(OpenMP 766.81 s, Serial 1024.16 s), with newly added bounded gates
run separately on both backends. The current bounded
oracles additionally cover CRAM/CRAI BQSR round-trip, BI/BD indel report parity,
ApplyVQSR multi-ALT
Number=A/culprit projection, and Picard/GATK SortSam ordering with the
`CREATE_INDEX` default/explicit boundary.

The latest bounded gates also cover GenomicsDBImport native sparse workspace
input-order preservation (including reopen/query missing-genotype behavior),
HaplotypeCaller `--dont-use-soft-clipped-bases` semantics, the Kokkos
GermlineFilter posterior/GERMQ boundary in FilterMutectCalls, and CombineGVCFs
`--sites-only-vcf-output`.  CombineGVCFs keeps full sample FORMAT state through
allele/DP merging and only then uses an independent writer header plus
`bcf_subset_format`; the pinned GATK 4.6.2.0 oracle confirms the 8-column
site-only shape and reference-block INFO behavior on both host backends.

FilterMutectCalls also has a pinned orientation writer oracle against GATK 4.6.2.0.
Its canonical/reverse-complement prior and `F1R2/F2R1 Number=R` MNV fixture produce
`ROQ=[93,93,1]` and the same orientation-filter pattern on both paths.  Native now
uses the GATK `QualityUtils.errorProbToQual()`-compatible `[1,93]` bounded phred
encoding, so zero/one artifact posteriors cannot serialize as the former `1000/0`.
This is a quality-serialization boundary; the complete FilterMutectCalls joint
posterior and Mutect2 release-specific model remain explicit fallback scope.
The contamination-table path now matches GATK's sample-aware lookup: a row for
the tumor sample overrides the fallback, while an unmatched sample uses the
configured contamination estimate (0.0 by default) instead of inheriting the
maximum value from an unrelated row.  The pinned contamination oracle includes
this unmatched-sample case on both OpenMP and Serial.
The GATK `--create-output-variant-index` short alias `-OVI <Boolean>` is also
accepted by the native parser, dispatcher, and help text.  A pinned Java/native
oracle fixes `-OVI false` to no `.tbi` sidecar and `-OVI true` to a `.tbi` sidecar;
GATK 4.6.2.0 rejects inline `-OVI=true`, so that spelling remains outside the
direct-replacement syntax contract.

ReblockGVCF's low-quality deletion path also clears stale `FORMAT/AD` while
preserving the deletion `END` span and materializing GATK-compatible
`GT=0/0`, `PL=0,0,0`, `GQ=0`, and `DP/MIN_DP`; the OpenMP/Serial contract asserts
these fields. This is a bounded annotation guarantee, not full annotation or
posterior parity for every ReblockGVCF mode.

The ReblockGVCF writer now also handles the Java buffered overlap boundary:
concrete variant/deletion spans trim or split an input hom-ref block before GQ-band
merging, and split tails use the indexed reference to update REF/POS. The focused
`scripts/verify_reblock_overlap_gatk_oracle.py` compares the pinned GATK 4.6.2.0
coverage `69000-69004 / 69005 / 69006-69010` on both OpenMP and Serial. Full
annotation/posterior semantics remain explicit fallback.

The decode boundary now uses `ReadBatch::bytes()` for HC and Mutect2 as well as
CollectReadCounts/GetPileupSummaries.  `AdaptiveController` may lower the next
reader batch only at a safe boundary, and `adaptive_batch_reductions` is recorded
in summaries/manifests.  HC/Mutect2 still retain an aggregate HostBatch for their
current assembly contract, so this is byte-aware decode backpressure—not a claim
of fully streaming whole-genome assembly.  HC also exposes opt-in
`--stream-by-contig`: coordinate-sorted input is flushed at contig boundaries,
VCF/GVCF bodies are written incrementally (including BGZF/Tabix), and the
per-contig likelihood/reference-confidence matrices are released after output.
It is guarded against contig re-entry and coordinate reversal; it is a contig
bounded path, not yet a fixed-memory AssemblyRegion spill implementation.
For indexed BAM/CRAM, `--stream-by-region N` adds a fixed-size core/halo VCF path:
each tile is fetched with HTSlib indexed iteration; the halo feeds activity/graph/
PairHMM while the core interval remains the output domain, and is released before
the next tile.  `streamed_regions`,
`streamed_peak_host_bytes`, and `stream_indexed` are recorded in the manifest.  The
one-tile VCF/GVCF path is byte-exact to aggregate HC; multi-tile VCF locus identity is
covered by a contract.  Multi-tile GVCF now carries a cross-tile HomRefBlock stitcher:
adjacent blocks in the same GQ band merge with element-wise PL/minimum-MIN_DP and a
rounded median DP reconstructed from per-site depths; candidate records remain tile-core
records.  A bounded chr17 fixture compares normal and 50 bp tiled GVCF output byte-for-byte.
This proves the block/output contract, while full GATK assembly/graph context remains a
separate whole-genome validation item.

Mutect2 uses the same core/halo reader and `ThreeStagePipeline` boundary.  Its
`calling::Options::intervals` remains the serialization core while
`locus_intervals` carries the halo into activity/graph/PairHMM, preventing
cross-tile evidence loss and duplicate writes.  A 500 bp tile is byte-exact on
the bundled fixture; deliberately smaller tiles are checked for ordered,
core-only, duplicate-free output because release-specific somatic posterior
calibration is not associative across independently computed tiles.

For whole-reference GVCF emission, the reference-confidence writer builds one Host
coordinate index for observed loci and a merged per-contig interval index for candidate
REF spans. Membership is answered with binary search instead of rescanning all loci and
candidates for every reference base; this changes only Host scheduling complexity (to
O(log N) lookup) and leaves Kokkos likelihood/PL bytes unchanged. The telemetry field
`reference_block_lookup_indexed=true` makes the optimized path auditable in large-contig
and streaming runs.
On the bundled `17:69000-70000` broad fixture, `benchmark_hc_broad.py` measured native
OpenMP p50 ≈1.47 s versus Java GATK p50 ≈5.41 s; the Serial artifact measured ≈1.39 s
with the same output/telemetry contract. These are fixture-level wall times (decode,
Host assembly, Kokkos kernels, and writing included), not a whole-genome speedup claim.

The dispatcher now stages HTTP(S), presigned object-store, and `file://` input
URIs into a content-addressed local cache before invoking native tools. Multiple
independent input URIs are prefetched with a bounded thread pool (override with
`FASTGATK_REMOTE_DOWNLOAD_THREADS`, capped at 64); results are consumed in
argument order, so network overlap does not change deterministic manifests. The
download path uses argv-based `curl` (no shell interpolation), retry/atomic
commit, required FASTA `.fai` checking, and opportunistic canonical `.dict`/
HTSlib sidecar caching. Remote outputs default to `BACKEND_UNAVAILABLE`; setting
`FASTGATK_REMOTE_OUTPUT_MODE=upload` enables HTTP(S)/presigned upload with
sidecars committed before the primary object and manifest URI rewriting.
`--fallback` preserves the original URI for Java GATK. The contract is
`fastgatk-native/scripts/verify_remote_staging.py`.

`fastgatk-hc-call` is the next-stage executable. It exercises a real, composable
`assembly -> likelihood -> genotyping` data path over the same HTSlib/Kokkos input:

1. Kokkos counts high-quality bases per `(tid, position)` and emits deterministic
   multi-ALT REF/ALT local haplotypes with shared observation ranges.
2. The reusable Kokkos k-mer graph stage counts stable nodes/edges over bounded
   AssemblyRegion reference windows, preserves `(tid,start,end)` path
   coordinates, read-support spans and alternate-edge metadata, prunes
   disconnected branches, and materializes deterministic path sequences;
   anchored CIGAR insertions/deletions and terminal soft-clip insertions
   (including reverse-strand orientation) become explicit left-anchored
   REF/ALT candidates. Materialized graph paths are aligned back with a
   Kokkos-scored Host SW traceback, and each retained haplotype owns a Host
   EventMap (SNP/MNP/indel/compound) that is the first graph candidate source.
   The bounded graph is not yet a full ReadThreadingGraph, so a documented
   transitional CIGAR/pileup supplement remains until those EventMaps cover
   the entire source haplotype population.
   Graph input is now restricted to reads whose CIGAR-aware reference span
   overlaps the merged AssemblyRegion windows (while preserving input order and
   BAM orientation flags), so unrelated contigs/windows cannot inflate local
   k-mer support or graph workspace.
3. Kokkos SW score and persistent bucketed PairHMM evaluate bounded local
   haplotype combinations plus coordinate-compatible graph paths, and
   marginalize REF/ALT-compatible paths.
4. A Kokkos genotyping kernel applies explicit configurable SNP/indel priors,
   selects the maximum posterior genotype/GQ, and emits PL in VCF/GVCF.

The same pipeline now runs a fixed-order Kokkos read-filter mask (flag/MAPQ/
coordinate/empty-read) and uses the one ordinal-stable mask for assembly,
graph, activity, SW/PairHMM and reference-confidence depth. `--disable-read-filter`
requires a filter class and is applied selectively; `--read-filter` supports the
the built-in mask classes, including CIGAR-validity, MAPQ-sentinel, non-zero
reference-span and `NonChimericOriginalAlignmentReadFilter` checks, plus
parameterized `ReadLengthReadFilter` with
`--min-read-length`/`--max-read-length` bounds. Unsupported/custom classes fail closed and are
eligible for explicit Java fallback rather than being silently ignored.
The same mask implements GATK's `MappingQualityAvailableReadFilter` and
`MappingQualityNotZeroReadFilter`, rejecting MAPQ=255 and MAPQ=0 sentinels when enabled;
their effective state is recorded in telemetry.
The default policy now matches GATK's HaplotypeCaller/Mutect2 filter floor
(MAPQ >= 20, duplicate reads excluded); `--minimum-mapping-quality` and
`--include-duplicates` are explicit overrides, and the selected policy is
recorded in OutputManifest telemetry. `--disable-tool-default-read-filters`
clears only implicit defaults while preserving explicit `--read-filter` names,
independent of argument order, and records the choice in HC/Mutect2 manifests.
Mutect2's implicit mask additionally enables `NonChimericOriginalAlignmentReadFilter`:
the HTSlib Host reader carries OA/XM auxiliary tags through batch merges and the
Kokkos predicate compares OA's comma-delimited contig prefix with XM; reads missing
either tag pass, matching GATK. Presence is retained independently of payload length,
so present-but-empty OA/XM tags follow GATK's `hasAttribute()` semantics. The effective
OA/XM policy is reported in telemetry.
It also applies GATK's default `ReadLengthReadFilter(30, Integer.MAX_VALUE)`; the
effective bounds are emitted in the manifest and telemetry.
Deterministic per-locus downsampling is
available through `--max-reads-per-locus` (or the GATK-shaped
`--max-reads-per-alignment-start` alias) plus `--downsampling-seed`; the cap is
applied before every downstream kernel and its retained/dropped counts are
recorded. Anchored compound/graph candidates that cross disjoint AssemblyRegion
boundaries are suppressed and counted, while candidates outside active regions
remain an explicit fallback.

Before assembly, the native Host applies GATK's `cleanOverlappingReadPairs`
contract to the filtered batch. Paired reads with the same name within the
resolved sample batch and reciprocal mate coordinates are CIGAR-projected to
their shared reference coordinates; concordant overlap qualities are capped at
Phred 20. HaplotypeCaller
uses GATK's conflict-to-zero mode, while Mutect2 preserves disagreeing mate
qualities (its Java call passes `setConflictingToZero=false`). The corrected flat
qualities feed assembly, graph, SW/PairHMM and reference confidence (the
GATK-compatible RCM evidence map retains its filtered-region qualities). GATK's
`--do-not-correct-overlapping-quality` (plus the descriptive
`--do-not-correct-overlapping-base-qualities` alias) disables the transform for
compatibility experiments, while pair/base/conflict/cap counts and metadata
availability are recorded in telemetry and OutputManifest.

When a reference is available, the activity-profile kernel also receives the
reference-projected base for every locus and computes the non-reference signal
from `1 - ref_count/depth`; `activity_reference_projection` in the manifest
distinguishes this path from the no-reference diversity fallback.
The AssemblyRegion state machine exposes GATK-shaped
`--active-probability-threshold`, `--assembly-region-padding`,
`--max-assembly-region-size`, and `--max-prob-propagation-distance` controls.
The active-span cap and propagation distance are independent, and padded
regions are never re-merged into an unbounded graph/PairHMM window; all four
values and the resulting region count are included in telemetry.

HC, Mutect2 and BaseRecalibrator accept repeatable `-L/--intervals`/`--region`
selectors through the same HTSlib Host reader. Selectors may be literal
`contig:start-end`, Picard interval-list/interval files, or zero-based
half-open BED files (plain or `.gz`; UCSC `track`/`browser` headers are
ignored); overlapping ranges are normalized before read selection,
and the original file and record counts are written to the manifest. The same
interval set is passed into HC reference-confidence block generation, so
disjoint gVCF intervals do not silently expand to the whole contig.

It emits candidate calls as JSON telemetry, VCF/VCF.GZ, or a reference-backed
reference-block plus candidate-site gVCF/gVCF.GZ (intervals include zero-depth
blocks; whole-contig mode emits merged uncovered blocks) (all with an `OutputManifest` sidecar) and an integer
`signature`. The input boundary now
preserves packed CIGAR operations, BAM flags, mate metadata, read names and RG
tags; the local SNP pileup uses the safe CIGAR reference projection and excludes
insertions/soft clips from reference-coordinate counts. This remains a bounded
local caller, not a claim of full GATK HaplotypeCaller equivalence: reference-backed calls now run the reusable Kokkos
Smith-Waterman score and bounded local haplotype-set PairHMM over each
halo-expanded AssemblyRegion (never a whole-contig candidate set), including
graph-materialized haplotypes only when their coordinate span overlaps that
region; the manifest records the region-group count and any explicit
out-of-region fallback candidates. Before a read is attached to an
AssemblyRegion PairHMM request, Host walks its CIGAR projection and applies
the region hard-clip span to bases, qualities, and BI/BD indel-quality arrays.
`pairhmm_reads_clipped` and `pairhmm_reads_dropped_after_clipping` make this
boundary auditable; Kokkos receives only the resulting flat read/haplotype
arrays; before genotype reduction, Host applies GATK's
`filterPoorlyModeledEvidence` threshold
`min(2,ceil(read_length*expectedErrorRatePerBase))*(-4.0)` (default rate 0.02).
Disqualified reads remain available to assembly/RCM but are excluded from
PairHMM genotype evidence and AD/DP; `pairhmm_reads_disqualified` and the
threshold are recorded in the manifest. The GATK-shaped
`--expected-mismatch-rate-for-read-disqualification` override is accepted by
HC/Mutect2. Genotyping windows use SNP/INDEL AssemblyRegion trimming, and
reference-backed VCF PLs remain prior-free while priors affect GT selection.
Graph SW substitution runs can be emitted as one phased MNP with
`--max-mnp-distance N` (or GATK's short `--mnp-dist N`; default `0`, preserving independent SNP records);
the run is accepted only when all mismatch columns remain in one ungapped
alignment segment, and the graph/read support floor is reused. `graph_mnp_candidates`
and the configured distance are written to telemetry/OutputManifest.
graph haplotypes are now scored and deterministically pruned by a configurable
`--haplotype-pruning-log10` delta relative to the best local haplotype, with the
best graph path retained and considered/kept/pruned counts in the manifest. SW
`M` segments are also compared base-by-base; graph-derived SNP mismatches are
promoted when the bounded graph/SW path has sufficient read support. If an
independent high-quality pileup count exists, it caps the graph support; an
empty/low-quality pileup therefore does not discard a graph-only SNP. The
promoted candidate then uses the shared PairHMM/genotype path. The manifest records these as
`graph_snp_candidates`/`graph_snp_candidate_extraction`.
The assembly boundary also records two GATK-shaped evidence gates before
PairHMM: terminal soft-clip-only alternate observations are suppressed, and
an alternate supported only by one overlapping read fragment is not allowed to
seed a haplotype.  `candidate_softclip_suppressed` and
`candidate_fragment_suppressed` make those decisions auditable.  The
`candidate_low_support_suppressed` counter additionally records the default
`min-pruning=2` boundary: an isolated singleton pileup SNP remains in the
PairHMM matrix but is folded back into reference confidence when the graph is
available and its site QUAL is below 30.  Explicit `--min-pruning 1`, somatic
mode, and short-read no-graph fallback retain singleton candidates.  The pinned
1 kb broad oracle is runnable with
`python3 fastgatk-native/scripts/verify_hc_broad_gatk_oracle.py`; its companion
file-boundary benchmark is
`python3 fastgatk-native/scripts/benchmark_hc_broad.py`.  In the latest pinned
Ryzen 9 7945HX run (same fixture and Java 17 invocation), p50 wall time was
1.48 s native versus 5.24 s Java GATK (about 3.5x); this is a small-I/O
reference point rather than a general speedup claim.  The manifest records
the support floor and pruning policy.  The broad oracle now requires call-set,
The 2026-08-30 rerun after the RCM/assembly-boundary update measured native
OpenMP p50 1.406 s (p95 1.425 s) versus Java GATK p50 4.552 s on the same
fixture; explicit Float32 was 1.388 s and explicit double 1.401 s.  These
numbers remain fixture-level measurements, not a whole-genome throughput claim.

The pinned `verify_hc_min_pruning_gatk_oracle.py` adds the strict assembly
boundary: on the same fixture, explicit GATK 4.6.2.0
`--min-pruning 3 --num-pruning-samples 1` removes the low-support singleton and
the native normal VCF is byte-identical with provenance disabled.  Native now
applies that explicit stricter floor even when its bounded reference-k-mer
graph falls back without materialized alternate paths; the oracle reports the
fold-back telemetry.  This slice does not claim GVCF block segmentation or
`--min-pruning 1` parity.

The pinned `verify_hc_softclip_gatk_oracle.py` covers the complementary
read-threading switch `--dont-use-soft-clipped-bases`: GATK 4.6.2.0's chr17
fixture emits the aligned 17:69298 SNP while terminal soft-clip evidence is
excluded, and native matches the four-site call set plus QUAL/GT/AD/DP/GQ/PL
and stable annotations.  Native no longer suppresses an already aligned base
merely because it lies near a terminal clip.  Coverage INFO/DP remains an
explicit residual (17:69298 is 4 native versus 6 Java), so this oracle does
not claim full Coverage/annotation-population parity.

The opt-in `--use-haplotype-realignment-for-rcm` switch now exercises the
read-to-haplotype boundary: Host selects each read's best normalized PairHMM
haplotype, reuses its SW traceback, and projects matched bases through a
haplotype/reference CIGAR before RCM.  Unmappable bases fall back per read to
the original projected observation.  The pinned fixture reports 3,265
realigned observations and produces the same decompressed GVCF hash as the
compatibility path; telemetry exposes both realigned and fallback counts.
With `FASTGATK_BENCHMARK_RCM_REALIGNMENT=1`, the same fixture reports a 1.39 s
p50 for the opt-in path (two repeats), so Host traceback projection is not a
dominant cost at this scale.
QUAL, INFO/DP, FORMAT/GT/AD/DP/GQ and PL identity, plus exact `MLEAC/MLEAF`,
`MQ/QD/FS/SOR` and all three rank-sum annotations, for all three pinned calls.
INFO/DP uses a separate GATK Coverage-style Host evidence count (FORMAT/DP
remains the retained informative likelihood depth); MQ additionally consumes
the PairHMM poorly-modelled-read source-record mask. The focused result is
core-plus-annotation bit-identical; the complete VCF is still not byte-identical
because provenance/command-line/date metadata, contig assembly metadata and
broader AF/annotation-engine semantics remain on the compatibility roadmap. The
three normal VCF record lines in the pinned broad oracle are now text-identical
to GATK, and all schema-bearing `INFO`/`FORMAT`/`FILTER` definitions plus stable
header order are exact. With explicit `--add-output-vcf-command-line false`,
the ordinary VCF is now byte-identical to GATK; the default true mode retains
truthful native provenance and therefore differs in its command-line/date line.

The HaplotypeCaller writer now honors `--add-output-vcf-command-line` on all
aggregate, contig-streaming, region-streaming, VCF, and GVCF paths.  When true
(the default), it emits a deterministic, implementation-labelled
`##GATKCommandLine=<ID=HaplotypeCaller,Version=fastgatk-native,...>` header;
when false, no such header is emitted.  The native version is intentional: it
keeps the provenance truthful while remaining valid for HTSJDK readers.  The
release-specific Java Date/Epoch fields are not fabricated; the explicit false
mode is the reproducible bit-identical profile.

`--native-pair-hmm-use-double-precision [BOOL]` is also accepted with GATK's
optional-boolean syntax.  The request and the effective kernel precision are
recorded independently in summary/OutputManifest.  HaplotypeCaller's strict
default remains direct `Kokkos::Experimental::simd<double>` because the
current Float32 path, while within the pairwise GKL replay bound, still changes
one pinned assembled-HC annotation payload.  Explicit `false` selects the
GKL-compatible Float32-first path: every pair starts in
`Kokkos::Experimental::simd<float>`, and a scaled sum below GKL's `1e-28f`
acceptance threshold is recomputed in the Kokkos Float64 recurrence before it
is returned.  Explicit `true` selects direct Float64.  The focused GKL
fallback oracle and real GATK debug-row replay verify the fallback boundary
separately.
Flow-space inputs use the calibrated double recurrence and explicitly fall back
to regular base-space Float32 when this mode is requested.

`benchmark_hc_broad.py` now runs aggregate strict/default, explicit strict
double, and explicit Float32 HC cases on the same fixture.  It reports the
effective precision, Kokkos SIMD width, output hash/size, and end-to-end timing
separately; Float32 is never compared as a bit-identical oracle.
Full GATK pruning/traversal parameters and reference-block priors are not yet
implemented. The native graph now exposes GATK-shaped `--min-pruning` in
addition to `--min-kmer-count`; it counts read k-mer observations for
non-reference branch retention while keeping the reference scaffold. The
genotype-prior transform matches GATK's `assumingHW` convention (each
alternate allele receives the configured SNP/indel heterozygosity in log10
scale; there is no extra SNP `log10(3)` normalization),
distinguishes OTHER alleles using `max(snpHet,indelHet)`, and applies the
resulting Number=G joint prior to optional VCF/gVCF GT/GQ assignment through
the shared Kokkos posterior-assignment primitive. GATK-compatible default
assignment remains PL-based; pass
`--genotype-assignment-method USE_POSTERIOR_PROBABILITIES` to enable the
prior-aware path, while `--no-genotype-priors` disables the prior model. The default `1e-3`
SNP and `1/8000` indel inputs plus `genotype_priors`/`joint_genotype_priors_used`
are recorded in telemetry. The graph stage now also materializes a bounded set of deterministic
path sequences with coordinates and activity/AssemblyRegion telemetry; graph
haplotypes are counted in PairHMM telemetry; the shared SW kernel now
implements GATK's four overhang strategies and Java tie priority. Concrete
Variable-length SW score requests are processed in deterministic chunks of 512
to bound ragged DP workspace and allocator pressure while preserving request
order and exact scores across Kokkos execution spaces.
multi-ALT VCF and candidate-site gVCF records now use one shared PairHMM
per-read likelihood matrix and the reusable Kokkos joint-genotype-PL kernel
to evaluate the full triangular genotype set, including GATK's log10
half-mixture for heterozygous and cross-ALT genotypes, with deterministic
GT/AC/AN/GQ selection; VCF/gVCF GT/GQ additionally
consumes the assumingHW joint prior only in explicit posterior-assignment mode.
The arbitrary-ploidy `GenotypePriorCalculator.assumingHW` Number=G enumeration
now uses the shared `calculate_genotype_priors_kokkos` API in both HC and
GenotypeGVCFs; Host retains only allele-type classification and the configured
heterozygosity values, while prior prepare/execute time and execution space are
recorded in the result/manifest. The kernel benchmark includes a deterministic
prior checksum and timing so Serial/OpenMP/CUDA/HIP paths cannot silently
diverge in genotype ordering.
In gVCF candidate loci,
`<NON_REF>` is now computed per read using GATK's qualified-concrete median rule
and included in the full triangular genotype matrix; the fallback envelope is
used only when no concrete per-read matrix is available. The complete GATK
assembly haplotype set, unseen-haplotype priors and cross-sample posterior
remain pending. For reference-backed diploid calls, the candidate call list
now uses the same PairHMM-derived PL rows and shared Kokkos GT/GQ selection as
the VCF writer; candidates without complete PairHMM evidence retain the
documented aggregate/pileup fallback. This removes a previous call-list versus
VCF disagreement at low depth. Reference-only blocks now carry REF-vs-any
AD/PL/GQ derived from the same filtered/projected observations and are merged
using GVCF GQ bands with `MIN_DP`; terminal soft-clip candidate discovery is
covered by the calling-indel smoke and recorded as `softclip_candidates` in
telemetry/manifest. A graph-derived internal soft-clip regression now covers
the bounded SW/CIGAR projection path; CIGAR/soft-clip indels are merged by
deterministic reference-repeat left-normalization, mixed `M/=/X/I/D` evidence
is covered, and short-distance compound `I/D` events on one read now produce a
single anchored ref/alt candidate; bounded same-read `I/D` + internal
soft-clip combinations now also produce one anchored ref/alt candidate with
an independent support floor; graph deletion branches use the same
candidate contract. The Host Number=G ordering is the VCF/GATK colexicographic
order (`0/0,0/1,1/1,0/2,...`), and candidate-site GT/GQ assignment consumes
only the concrete PL prefix; `<NON_REF>` rows remain an emitted likelihood
envelope. FORMAT/AD keeps the trailing zero for that symbolic allele. The
chr17 oracle now checks candidate-site GT/GQ/AD/PL record text, QUAL/FILTER,
standard INFO order and `FORMAT/SB` exactly; candidate PL is intentionally
uncapped like GATK's Number=G path. Reference-only block FORMAT/DP is GATK's rounded
median depth and FORMAT/MIN_DP is the block minimum, not the first-site depth.
HC now also exposes GATK's call-confidence contract: ordinary VCF mode defaults
to `--standard-min-confidence-threshold-for-calling 30` (alias
`--stand-call-conf`), while `-ERC GVCF`/`BP_RESOLUTION` forces the threshold to
zero so GenotypeGVCFs owns cohort filtering. Candidate-site QUAL is computed by
the shared batched Kokkos biallelic AF/Dirichlet EM kernel rather than a Host
floating-point shortcut, and is written to VCF QUAL/manifest telemetry. The
threshold gate is semantically wired, but native PairHMM/assembly likelihoods
remain numerically different from Java GATK on broader corpora; the bounded
oracle therefore invokes the explicit zero-threshold mode when testing evidence
and output shape.
For normal multi-ALT VCF records, site QUAL now uses the same shared Kokkos
`calculate_allele_frequency_kokkos` AFCalculator/Dirichlet-EM posterior over the
complete Number=G PL vector (instead of the maximum independent ALT QUAL).
Coverage INFO/DP is restricted to PairHMM-retained evidence; the tetra-ploid
multi-ALT oracle now locks the complete candidate row, including QUAL, INFO,
RankSum and FORMAT text, to GATK for this corpus. Broader annotation-engine
provenance and assembly coverage remain explicit expansion work.
Arbitrary-ploidy VCF/GVCF candidate records omit GATK's diploid-only `ExcessHet`
field; diploid records retain the exact-test value.
Full assembly-region haplotype combinations,
unseen-haplotype posterior semantics,
filters outside the nine native classes and the remaining release-specific
graph pruning/traversal semantics remain pending; bounded graph k-mer size,
support floor, path count and path depth are now explicit HC parameters and
are recorded in telemetry. GATK-shaped `--min-dangling-branch-length` now
removes short non-rejoining graph paths at the bounded traversal endpoint;
`--recover-all-dangling-branches` only relaxes the Host walk through forks in
an already reference-connected graph, selecting its heaviest supported edge
as GATK's `findPath` does.  It never preserves or globally SW-promotes a
disconnected component: that component is removed by the matching
`removePathsNotConnectedToRef` stage, using the GATK directed
reference-source-to-reference-sink intersection after tail/head graph
recovery. `--min-dangling-matching-bases`
(default 3) remains the exact CIGAR suffix/prefix requirement for a connected
tail/head merge. The manifest records pruned and CIGAR-recovered branch
counts/base spans.

The graph wrapper also implements the GATK cyclic-graph policy: when a
read-supported directed cycle is observed, it retries with `k+10` up to `31`
(bounded by the local read/reference length). `--dont-increase-kmer-sizes-for-cycles`
disables that retry, and `--max-num-haplotypes-in-population` is accepted as the
GATK alias of `--max-haplotype-paths`. Selected k-mer size, retry count, and
whether any retry observed a non-reference cycle are emitted in telemetry.
The same graph stage reports SeqGraph-style non-branching-chain compression
(`graph_seqgraph_nodes/edges`) and now applies an actual sequence-level
SeqGraph rewrite: shared-prefix/suffix diamonds, GATK-shaped common-suffix
tails (default minimum 10 bases), and linear-chain zipping are materialized as
shared vertices. Rewrite counts (`graph_seqgraph_*_merges` and
`graph_seqgraph_suffix_*`) and the reconstructed path language are exposed to
Host callers; linear-chain, diamond/tail, SNP-rethreading and two-ALT
support-order fixtures are part of the native/kernel API smoke. Alternate
outgoing edges remain ordered by observed read multiplicity with stable node-id
tie-breaks.
Graph paths now also carry reference-connected provenance: the shared
Smith--Waterman stage uses GATK's `NEW_SW_PARAMETERS`
`(match=200,mismatch=-150,gap-open=-260,gap-extend=-11)` with the shared
Kokkos score recurrence, then emits one stable reference CIGAR through the
Host traceback,
alignment score, and absolute offset per materialized path. These values are
available in `graph_haplotype_path_details`, together with a deterministic
`graph_haplotype_cigar_signature`; the real assembly-graph contract rejects a
path that has sequence/coordinates but no valid traceback. This is a bounded
provenance layer; `graph_haplotype_sw_*` telemetry records the execution space,
SIMD batch width/groups, and score timing so the graph path cannot silently
fall back to a separate CPU-only scoring implementation. For the actual
read×haplotype requests staged for PairHMM,
the manifest also records valid/indel CIGAR-pair counts and the numbers of
uncertain/informative read-at-locus observations after haplotype
max-marginalization under GATK's 0.2 log10 BestAllele margin, together with a
deterministic uncertainty signature. The shared Kokkos `reduce_read_allele_best_kokkos`
also returns the concrete REF/ALT owner for each informative `(locus, read)` row;
HC and Mutect2 use that one-read-one-allele result for multi-ALT depth/counts, so
the same read is not counted independently against every biallelic ALT. The Host then projects each candidate
interval through that same read×haplotype CIGAR: associations whose local
alignment explicitly soft-clips away the event are excluded from allele
marginalization (unclipped mismatch/indel haplotypes retain their reference
likelihood), with soft-clip-filter counts exposed in telemetry. This makes the boundary
auditable, but
does not claim full GATK read-haplotype posterior, all
soft-clip/hmer semantics, or cross-AssemblyRegion assembly parity.
Before graph construction and PairHMM, overlapping/contiguous raw activity
windows are now coalesced into a bounded calling-region union. A compound
candidate that crosses only a raw halo boundary is therefore retained in one
local graph; genuinely separated active islands remain independent. The
manifest exposes `assembly_region_union_count`, `assembly_region_union_merges`,
and `assembly_cross_region_rescued` so this behavior is auditable rather than
silently changing the region partition.
The reusable kernel API additionally returns the sorted per-sequence
`non_unique_kmers` set used by ReadThreadingGraph-style dynamic-k decisions;
HC/Mutect2 surface combined/reference repeat counts and the GATK-compatible
skip-and-increase-k decision as `graph_non_unique_kmers`,
`graph_reference_non_unique_kmers` and `graph_reference_kmer_rejected` in
telemetry. Production callers reject a repeated reference k-mer by default;
the low-level kernel retains an explicit compatibility override.
The compatibility surface also accepts GATK's `--num-pruning-samples`,
adaptive-pruning controls, and linked-de-Bruijn toggles. Linked mode is marked
experimental and deliberately keeps the bounded raw topology (no SeqGraph
rewrite). When enabled, uncovered high-support pivotal edges are stapled onto
the best completed path; `graph_artificial_haplotype_recovery_paths/bases`
records this work. `--disable-artificial-haplotype-recovery` and the legacy
cycle toggle are recorded in the manifest so a workflow can pass through the
same command line without silently changing modes.
The graph input also carries BAM flags, so reverse-strand reads are oriented
before k-mer extraction instead of becoming disconnected branches.

The hidden GATK `--error-correct-reads` path is now available as an explicit
assembly preprocessing mode. A Kokkos RangePolicy performs conservative
nearest-solid-k-mer correction (up to two Hamming mismatches) using exact
Host-counted k-mers, followed by strict consensus across overlapping windows;
the controls are `--kmer-length-for-read-error-correction`/
`--min-observations-for-kmer-to-be-solid`. Corrected bases and Q30 qualities are used only for candidate assembly and
graph construction; original filtered reads continue into PairHMM and the
reference-confidence model. Correction counts and timings are recorded in
the manifest, and the default remains off to preserve ordinary GATK behavior.

The hidden `--error-correction-log-odds` selector now exposes GATK's
`PileupReadErrorCorrector` path. A finite threshold computes the same flat-
Beta Mutect2 pileup log-likelihood ratio, excludes bases adjacent to CIGAR
indels/soft clips, applies GATK's three-mismatch indel guard, and takes
precedence over nearby-k-mer correction. The corrected assembly copy is still
the only consumer; PairHMM/RCM retain the original filtered reads. `-infinity`
keeps this mode disabled and both correction modes report their mode and
pileup-specific telemetry in the manifest.

Reference-backed VCF output also accepts `--sample-ploidy 1..8`: the writer
uses the VCF Number=G non-decreasing genotype ordering to generate arbitrary-
ploidy PL from the per-read PairHMM matrix, then uses the shared Kokkos
genotype-combination API for GT/GQ and the corresponding AC/AN values. The
manifest records the selected ploidy. Reference-backed `-ERC/--gvcf` now uses
the same ploidy-specific REF/`<NON_REF>` Number=G block model and candidate
writer; opaque/no-reference inputs still fail closed rather than inventing a
polyploid posterior.

`-ERC BP_RESOLUTION` is also implemented: the Host reuses the same Kokkos
reference-confidence likelihoods but disables GQ-band coalescing and emits one
REF/`<NON_REF>` record per reference base (`GT:AD:DP:GQ:PL`, `INFO=.`), while
candidate sites remain separate records. The BAM/CRAM `@RG SM` sample name is
preserved in both VCF and gVCF headers, with `FASTGATK` only as the no-tag
fallback, so Nextflow/SLURM outputs can be consumed directly by HTSJDK/GATK.
When present, each `@SQ AS` assembly tag is also carried into the VCF/GVCF
`##contig` line, matching HTSJDK's reference-dictionary metadata contract.

Both native executables link the reusable top-level `fastgatk-kernels` CMake
target. It owns the Kokkos PairHMM API and the affine-gap Smith-Waterman score
plus Host traceback/CIGAR API, allowing callers to share one Kokkos
configuration and implementation. The regular PairHMM recurrence uses
`Kokkos::Experimental::simd<double>` for same-shape read/haplotype lanes (the
benchmark reports the selected width, currently 8 on the OpenMP build), while
the flow recurrence uses Kokkos `TeamPolicy` for device-portable parallelism.
The SW score kernel uses
`Kokkos::Experimental::simd<int>` across independent same-length alignments
and reports the selected width/groups; ragged batches use the identical Kokkos
scalar recurrence. The API smoke includes the Java-GATK overhang/tie corpus
(including original-default, indel, substring and long homopolymer vectors) and
checks SIMD scores against the Host reference; the remaining HC gap is full
graph/pruning and caller semantics.

Flow-based inputs now have a Host-only `FlowBasedRead` tag boundary: HTSlib
preserves per-read `tp`/`t0` and read-group `FO`/`mc`, and
`fastgatk::io::decode_flow_read` produces the calibrated `[flow][256]` table
for the separate flow-space PairHMM API. The decoder includes GATK's
`clipProbs`, boundary-flow and max-hmer rules; `encode_flow_key` also exposes
the `FlowBasedHaplotype` key-to-base/reverse mappings and clipping pairs. HC/
Mutect2 now automatically select the flow recurrence when every participating
read has complete, consistent `tp`/flow-order metadata; ordinary or mixed
batches stay on regular PairHMM with an explicit fallback reason. The checked-in
GATK `sample.bam` and `sample.t0.bam` key/matrix corpus is now an opt-in CTest
oracle (including the fixture's phase-normalized leading-zero representation),
and four fixed read/haplotype cases compare the Kokkos flow likelihoods
bit-for-bit with GATK `FlowBasedPairHMM`.
The flow recurrence now launches through a Kokkos `TeamPolicy` (one independent
team/workspace per read-haplotype pair) on Serial/OpenMP and device-capable
backends; the Strict scalar recurrence order is unchanged. Host
read/reference-span plus a bounded local haplotype-vs-reference SW-CIGAR now
materialize request-local flow haplotype windows and record clipping/fallback
telemetry. Flow reads also apply the same AssemblyRegion base hard-clip to the
decoded `FlowBasedRead` key and `[flow][256]` probability matrix on Host
(`pairhmm_flow_reads_clipped`), including boundary hmer shifts and optional
probability spreading before the Kokkos TeamPolicy recurrence. Full GATK
haplotype-CIGAR coordinate/hmer/uncertainty parity,
anti-diagonal GPU wavefront tiling and GKL
flow raw-bit corpus parity remain explicit follow-ups; the checked-in
GATK FlowBasedPairHMM corpus is already a Strict raw-bit gate.
显式 `--flow-assembly-collapse-hmer-size N` 现已接入 HC/Mutect2：Host 在 flow
haplotype key 编码前复现 GATK `collapseBases()` 的确定性策略（首个 homopolymer
保留，后续 run 截断到 `N`），并用 `pairhmm_flow_haplotypes_collapsed` 记录受影响
haplotype 数；`0` 为默认关闭，`-1` 从当前 batch 的最大 read-group `mc` 自动选择阈值。
shared calling library 现在还提供 `uncollapse_flow_hmers`：按 GATK 的双向
SW/CIGAR 规则恢复 reference deletion 中的 hmer，支持 partial/limit-to-threshold
模式；flow batch materialize post-likelihood restoration view，并以
`pairhmm_flow_haplotypes_uncollapsed` 记录恢复数。随后按恢复后的 bounded
`(tid,start,end,sequence)` 建立 canonical group，合并 candidate ownership/mask，
并让 Kokkos max-marginalization 消费 canonical entry，避免 identical collapsed
haplotype 重复计入；`pairhmm_flow_haplotype_remaps` 和
`pairhmm_flow_identical_haplotype_groups` 写入 HC/Mutect2 manifest。跨
AssemblyRegion posterior 仍是后续兼容性门禁。
启用前还会在当前 bounded reference window 复现 GATK `needsCollapsing` 判断，普通
reference 区域不会因 alternate haplotype 的偶然长 run 被误截断。
Regular PairHMM preparation now follows GATK's default Host semantics: base qualities
are capped by MAPQ (255 is treated as unavailable) and values below Q18 are mapped to
Q6; missing BI/BD tags start at raw insertion/deletion Q45 with gap continuation Q10.
The Conservative PCR curve lowers only non-terminal positions (effective Q40 at
repeat length zero); the terminal position remains Q45. Per-base
BI/BD indel-quality Z tags are decoded into `ReadBatch`, preserved through filtering
and batch aggregation, and consumed by the PairHMM; missing tags use the documented
the model default when tags are absent, while length-mismatched tags fail
closed. The pinned PairHMM default-quality oracle checks this terminal-Q45
boundary against GATK 4.6.2.0; explicit BI/BD tags retain their existing decode
and PCR-adjustment path.
The GATK-shaped `--pcr-indel-model NONE|HOSTILE|AGGRESSIVE|CONSERVATIVE` switch is accepted by
HaplotypeCaller and Mutect2; `NONE` restores the explicit Q45 no-tag behavior
for legacy runs. `HOSTILE`/`AGGRESSIVE`/`CONSERVATIVE` apply the GATK-shaped
tandem-repeat quality curve in Host preparation; selected model, rate factor,
and adjusted read positions are emitted in manifest/stats telemetry.
The PairHMM floor is configurable with the GATK-shaped
`--base-quality-score-threshold` (default 18, bounded to Q6..Q255);
`--disable-cap-base-qualities-to-map-quality` disables the MAPQ cap while retaining
the threshold floor. Both controls are shared by HaplotypeCaller and Mutect2.
The GATK-compatible `--dont-use-soft-clipped-bases` switch is also shared by both
callers: it removes terminal/internal soft-clip-derived assembly evidence before
graph construction while retaining explicit CIGAR I/D events. The selected policy
is recorded in stats and manifests.
Mutect2 also accepts GATK's `--min-base-quality-score`/`-mbq` (plus the
`--min-base-quality` compatibility alias); the default Q10 floor is applied in
the shared Host assembly/activity/PairHMM input path.
After PairHMM, GATK-shaped per-read likelihood normalization now runs through
the shared Kokkos `normalize_likelihoods_kokkos` API: Host supplies read IDs and
the graph-path eligibility mask, while Kokkos performs the stable best reduction
and applies the 4.5-log10 global-mismapping cap. A negative Q is passed as
GATK's `-Infinity` disabled-cap sentinel, preserving valid `-Infinity`
zero-probability cells instead of replacing them with a large finite floor; Q=0
remains rejected by the GATK boundary. The Host retains CIGAR/allele semantics
and the clipped-read poorly-modelled gate; normalization prepare, execute, and
execution-space telemetry are recorded as
`pairhmm_normalization_*`.
The `--phred-scaled-global-read-mismapping-rate` boundary is release-pinned as
well: positive Q enables the cap, a negative Q disables it, and zero is rejected
by GATK (it is not an alternative spelling for disabled). Native rejects zero
before assembly so an accidental zero cannot flatten every read's likelihoods
and change TLOD; `verify_mutect2_mismapping_rate_boundary.py` exercises the same
failure boundary against GATK 4.6.2.0 and directly checks the Java disabled-cap
probe `[0.0,-Infinity,-3.0]`.
The REF/ALT max-marginalization is centralized in Kokkos through
`marginalize_read_allele_likelihoods_kokkos`; Host retains only CIGAR/graph
eligibility and stable row construction, with independent
`pairhmm_marginalization_*` telemetry. HC now keys rows by stable
`(AssemblyRegion,locus,source-read)` and passes the concrete ALT index
`1..N` into the explicit allele-count overload, so a multi-ALT read is
marginalized once and its per-ALT likelihoods feed AD/DP/PL without a hidden
biallelic truncation. The rectangular Kokkos stride is the maximum ALT count
among loci in that batch; unused cells remain `-inf`. Mutect2's full somatic
cross-ALT posterior and cross-AssemblyRegion ownership remain explicit follow-up
gates.
Region streaming is now allocation-aware: when an indexed halo tile exceeds
the runtime's conservative Host-memory budget, the core interval is split and
requeried in stable left-to-right order until the tile fits. A small fixed
allowance covers the decoded-batch object/HTS metadata; at a single-base core
the tile is accepted when it still fits the hard cgroup/SLURM allocation and
otherwise fails with `RESOURCE_EXHAUSTED`. Halo reads are recomputed for each
child, while only core records are emitted; therefore cross-boundary evidence
is retained without duplicate output. The manifest and summary expose
`streamed_region_splits` and the peak staged bytes.
Read×haplotype uncertainty now follows the same boundary: Host supplies stable
(locus, read) rows and finite likelihoods, while
`reduce_read_allele_uncertainty_kokkos` performs the deterministic best/second
reduction through Kokkos `RangePolicy`. Regular and flow PairHMM callers share
this path, with `pairhmm_uncertainty_*` prepare/execute/execution-space telemetry.
Reference-connected graph SNP paths now enter this reduction only when their
bounded SW CIGAR covers the candidate without a soft-clip; ambiguous
repeat-kmer projections remain fail-closed and are counted by
`pairhmm_graph_snp_posterior_pairs`. Dangling graph paths whose metadata span
extends beyond the emitted sequence are clipped to the available suffix; the
missing tail remains `unknown` provenance and cannot contribute synthetic
REF/ALT posterior evidence.

## Build

```bash
bash fastgatk-native/scripts/build_native.sh
python3 fastgatk-native/scripts/verify_native.py
# The verifier also runs fastgatk-hts-reader-test, including synthetic
# insertion/deletion/soft-clip projection edge cases.
```

When the pinned GATK package and JDK are present, CMake also builds
`fastgatk-genomicsdb-export` and extracts only its matching native library:

```bash
third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake -S fastgatk-native \
  -B fastgatk-native/build -DFASTGATK_ENABLE_GENOMICSDB_BRIDGE=ON
third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake --build fastgatk-native/build \
  --target fastgatk-genomicsdb-export fastgatk-genotype-gvcf
```

For a separately installed GATK release, set
`FASTGATK_GENOMICSDB_LIBRARY` and `FASTGATK_GENOMICSDB_JVM_LIBRARY` to matching
files. The bridge is optional and version-coupled by design; it is never used
to guess at an incompatible workspace.

The build script uses HTSlib 1.22.1 and zlib 1.3.1 from `third_party/sources/`, with
optional BZ2/LZMA/libcurl features disabled because this environment does not provide
their development headers. BAM and CRAM decoding are both exercised against GATK test
fixtures; CRAM requires its matching reference FASTA.

Example:

```bash
fastgatk-native/build/fastgatk-hc-smoke \
	-I reads.bam -L 17:69000-69100 \
	--batch-records=1024 --threads=16 -O calls.vcf \
	--output-manifest calls.vcf.manifest.json --telemetry calls.telemetry.json

fastgatk-native/build/fastgatk-hc-call \
  -I reads.bam -R reference.fa -L 17:69000-69100 \
  --batch-records=1024 --native-pair-hmm-threads=16 -O calling.vcf \
  --output-manifest calling.vcf.manifest.json --telemetry calling.json

# Repeatable interval selectors support GATK UNION (default) or INTERSECTION
fastgatk-native/build/fastgatk-hc-call \
  -I reads.bam -R reference.fa -L 17:69000-70000 -L 17:69300-69400 \
  -isr INTERSECTION -O calling-intersection.vcf

# GATK-compatible repeatable inputs: shards are merged before one calling pass
fastgatk-native/build/fastgatk-hc-call \
  -I shard-1.bam -I shard-2.bam -R reference.fa -O calling.vcf

# Explicit approximate genotype-prior policy (recorded in telemetry/manifest)
fastgatk-native/build/fastgatk-hc-call \
  -I reads.bam -R reference.fa -O calling.vcf \
  --heterozygosity 0.001 --indel-heterozygosity 0.000125

# GenotypeGVCFs posterior QUAL is explicit opt-in; GT/GQ/PL remain PL-based
fastgatk-native/build/fastgatk-genotype-gvcf -V sample.g.vcf.gz \
  -O joint.vcf.gz --gp-qual --heterozygosity 0.001 \
  --indel-heterozygosity 0.000125
```

`-O *.vcf` (or `*.vcf.gz`) writes a valid VCF header and contig dictionary. If
`--output-manifest` is omitted for a VCF output, the default sidecar is
`<output>.manifest.json`; BGZF-compressed VCF/GVCF also receives a tabix
`<output>.tbi` index when HTSlib is enabled, while plain VCF receives a
HTSJDK/Tribble LinearIndex v3 `<output>.idx`; `--output=summary.json` retains
JSON telemetry output.
`--create-output-variant-index` follows GATK optional-boolean syntax (flag,
`true/false`, or `=true/=false`); passing `false` suppresses both `.tbi` and
`.idx` sidecars as appropriate.
Repeated `-I/--input` values are combined in order with sequence-dictionary and
(unless `--sample-name` is explicit) sample-set validation. `--stream-by-contig`
uses a bounded coordinate k-way merge across shards, while
`--stream-by-region` merges each indexed halo tile; both retain the same memory
budget enforcement. A single region tile is equivalent to aggregate output;
multiple bounded tiles independently rebuild activity/graph context, so local
depth/candidate evidence may differ from aggregate. The region contract
enforces coordinate ownership, sorted/unique records, index validation and
backpressure telemetry rather than claiming whole-VCF byte identity.
With `-R`, the Host validates FASTA/FAI contig order and lengths against the
BAM dictionary before calling when a readable `.fai` is present. Unindexed
plain FASTA remains a documented smoke mode without an authoritative
dictionary; a present but malformed `.fai` fails closed. The explicit optional
boolean `--disable-sequence-dictionary-validation` skips this check, and its
value is recorded in the summary/manifest.
The smoke manifest records `status=smoke`; the calling prototype records
`status=prototype`. Both explicitly report output completeness, HTSlib/Kokkos
telemetry, and that bit identity with GATK is not yet claimed.
Reference-backed HC candidate records now also emit standard site annotations
`MLEAC/MLEAF`, `MQ`, `QD`, `FS`, `SOR`, `MQRankSum`, `ReadPosRankSum`, and
`BaseQRankSum`, plus the diploid single-sample `ExcessHet=0.0000` exact-test
result, calculated from the same filtered
CIGAR-projected evidence used by the call. `MLEAC/MLEAF` use the bounded GATK
`AlleleFrequencyCalculator` EM/Dirichlet reduction rather than copying the
final GT for current diploid candidate groups. Unsupported indel rank
statistics and complex likelihood-best-allele cases fail closed rather than
fabricating values; small-sample exact permutation is bounded to 18
observations and implemented on Host. The contract is covered by
`scripts/verify_hc_variant_annotations.py`.

GenotypeGVCFs 还对二倍体双等位多样本记录计算 GATK `ExcessHet`：采用
Wigginton exact HWE right-sided tail，四样本全 het 的 contract fixture 输出
`6.4098`；polyploid/multiallelic 记录 fail-closed，单样本退化记录输出
`0.0000`。完整 pedigree/founder-aware 语义仍未宣称 bit-identical。
ActivityProfile 使用同一套 Kokkos kernel：HC 继续使用 reference-aware
non-reference signal；Mutect2 CLI 额外传入 per-locus base-quality 列表，按 GATK
`Mutect2Engine.logLikelihoodRatio` 的 flat-Beta/digamma/Phred 模型生成二值活动状态，
并支持 `--pcr-snv-qual` 与 `--base-qual-correction-factor`。Host 复现
`BandPassActivityProfile` 的 Gaussian 累积、adaptive filter size、延迟 pop、
local-minimum cut 和 forced flush；默认阈值为 GATK 的 `0.002`，默认
`minAssemblyRegionSize=50`、`MAX_FILTER_SIZE=50`、`sigma=17`。manifest 会记录
`activity_filter_size`、有效 propagation distance 和 quality-aware 参数；profile 的
变长状态机和 gap/halo 控制保持在 Host。Mutect2 保留独立 AssemblyRegion，即使 halo
重叠也不合并，以免改变 GATK 的局部 graph/PairHMM 边界；HC 仍使用共享 API 的有界
union 策略。这样 Serial/OpenMP/CUDA/HIP 可以复用相同 Kokkos kernel。
When the calling path reaches Smith--Waterman, its manifest/summary also records
`sw_simd_width` and `sw_simd_groups`: rectangular read/haplotype batches use the
portable `Kokkos::Experimental::simd<int>` score kernel, while ragged batches
use the same Kokkos scalar kernel and preserve the Host traceback contract.

The file-boundary benchmark is separate from kernel-only timing:

```bash
python3 fastgatk-native/scripts/benchmark_end_to_end.py > e2e-benchmark.json
```

It records HC/Mutect2 wall time, output bytes, Kokkos prepare/execute telemetry,
compression/Tabix indexing and `/usr/bin/time` RSS/I/O counters in one JSON
artifact (schema v4 also records HC `--stream-by-contig` and indexed region-stream cases,
contig/region/peak-byte
telemetry, PairHMM request-link deduplication, bounded combination blocks and SW
SIMD width/groups from each output manifest). When the
pinned Java 17/GATK 4.6.2.0 jar is present it also runs a
same-fixture Java HaplotypeCaller baseline; set
`FASTGATK_SKIP_GATK_BENCHMARK=1` to omit it. The baseline is a performance
reference only—call-set equivalence still requires the separate GATK oracles.
When a Nextflow executable or the pinned framework JAR is available,
`verify_pipeline_local.sh` also executes the actual DSL2 smoke, scatter/gather,
and generic dispatcher workflows; otherwise it runs the DSL2-equivalent process
and local SLURM allocation simulation. The current validation environment has
Nextflow 26.04.6 available and has exercised all three real workflows. This is
local workflow evidence, not a real SLURM allocation or multi-node production
run.

The kernel benchmark (schema v7) also reports the compiled Kokkos SIMD ABI
(`compiled_simd_backend`, `compiled_simd_width`) and host AVX2/AVX-512
capability flags, in addition to
`joint_genotype_prepare`/`joint_genotype_kernel` timings for the shared
arbitrary-ploidy Number=G PL materializer, alongside
`reference_confidence_polyploid_prepare`/`reference_confidence_polyploid_kernel`
for the shared biallelic arbitrary-ploidy (including diploid) reference-confidence
PL/GQ materializer,
`genotype_posterior_prepare` and `genotype_posterior_kernel` for the
GenotypeGVCFs posterior path, plus
selected graph k-mer size, cyclic-retry iterations, and cycle-observed status;
this keeps the new ReadThreadingAssembler retry policy visible in the same
OpenMP/Serial benchmark artifact:

```bash
python3 fastgatk-native/scripts/verify_kernel_benchmark.py
python3 fastgatk-native/scripts/verify_kokkos_backend_matrix.py
```

On the current AMD Ryzen 9 7945HX (Zen 4), the OpenMP artifact is configured
with `Kokkos_ARCH_ZEN4` and reports `compiled_simd_backend=avx512`, eight
double lanes (sixteen integer lanes for SW), and `cpu_supports_avx512=true`.
With `OMP_NUM_THREADS=1`, the mixed-kernel benchmark measured PairHMM p50
about 1.80 ms versus 8.03 ms for the scalar Serial artifact, and uniform SW
about 1.26 ms versus 6.59 ms; all PairHMM/SW/genotype checksums were identical.
These are kernel-only same-host measurements, not a whole-pipeline or
Java-GATK speedup claim.

For an existing `gatk HaplotypeCaller ...` process, use the compatibility launcher:

```bash
fastgatk-native/scripts/gatk_compat.sh --java-options '-Xmx8g' HaplotypeCaller \
  -I reads.bam -L 17:69000-69100 -O calls.vcf
```

For a literal executable-name replacement, `fastgatk-native/dispatcher/gatk`
is also provided. Prepend `fastgatk-native/dispatcher` to `PATH` (or install
that file as `gatk`) and existing `gatk <Tool> ...` commands use the same
registry validation, remote staging, and explicit Java fallback behavior.

For an install-tree deployment, configure and install with CMake:

```bash
cmake --install fastgatk-native/build --prefix /opt/fastgatk
export PATH=/opt/fastgatk/bin:$PATH
gatk --version
```

The install includes production native binaries, the dispatcher registry,
and the Java-generated PairHMM table under `share/fastgatk`; the generated
launchers set the installed binary/registry roots without requiring source
paths or `FAST_GATK_PAIRHMM_TABLES`.

The HC parser accepts both `-ERC/--gvcf` and GATK's
`--emit-ref-confidence` spelling; unsupported combinations are reported as an
explicit backend failure so the compatibility dispatcher can route them to
Java. `-XL/--exclude-intervals` is applied after `-L/--intervals` at projected
locus granularity, with `-ip/-ixp` padding shared by indexed and sequential
HtsReader traversal; a read spanning excluded and included loci retains its
included evidence.

For multi-sample BAM/CRAM inputs, `--sample-name NAME` now matches GATK's
sample selector: HTSlib resolves `@RG ID→SM`, drops reads from other samples
before Host batching, and emits the selected sample as the sole VCF/GVCF
column. The selected name and `sample_name_selection` contract are recorded in
the summary and OutputManifest; an unknown name is rejected rather than mixed
silently with another sample.

The local pipeline contract test exercises GATK aliases, a two-interval
scatter/gather with indexed output, and the SLURM resource environment without
requiring a cluster. The second scatter invocation reads the shard
`OutputManifest` sidecars and verifies that complete shards are resumed rather
than recomputed. It also runs the actual Nextflow DSL2 files when a `nextflow`
executable or the pinned framework JAR is available:

```bash
bash fastgatk-native/scripts/verify_pipeline_local.sh
```

The scheduler-neutral runner can be used directly from a CI job or a SLURM
allocation:

```bash
fastgatk-native/workflow/scatter_gather_smoke.sh \
  -I reads.bam -L '17:69000-69050,17:69051-69100' \
  --binary fastgatk-native/build/fastgatk-hc-smoke \
  --gather-binary fastgatk-native/build/fastgatk-gather-vcfs \
  --outdir scatter-results --threads 2
```

The equivalent Nextflow DSL2 workflow is
`fastgatk-native/workflow/nextflow_scatter_gather.nf`; invoke it with
`-resume`. Shards keep GATK `-I/-R/-L/-O` names, and the gather stage emits
`gathered.vcf.gz`, `.tbi`, and an `OutputManifest` sidecar. The adjacent
`nextflow.config` has a default `local` profile and a production `slurm`
profile:

```bash
nextflow run fastgatk-native/workflow/nextflow_scatter_gather.nf -profile slurm -resume \
  --input reads.bam --intervals '17:69000-69050,17:69051-69100' \
  --binary fastgatk-native/build/fastgatk-hc-smoke \
  --gather_binary fastgatk-native/build/fastgatk-gather-vcfs \
  --slurm_queue short --slurm_cluster_options '--account=genomics'
```

For a generic drop-in process that should keep the tool name and launcher
stable, use `fastgatk-native/workflow/nextflow_gatk_compat.nf`. It calls the
same `dispatcher/gatk` executable used by an existing `gatk <Tool> ...`
process, keeps `-I/-R/-L/-O` plus `--extra_args`, and publishes both the
primary output and its OutputManifest sidecar:

```bash
nextflow run fastgatk-native/workflow/nextflow_gatk_compat.nf -resume \
  --tool CountReads --input reads.bam --output result.txt \
  --launcher fastgatk-native/dispatcher/gatk \
  --extra_args '--batch-records 4096 --threads 8' \
  --outdir nextflow-results
```

The same process can select any registry tool or an explicit Java fallback;
unsupported native parameters remain visible to the dispatcher instead of
being silently discarded. The local pipeline verifier runs this generic
workflow when a Nextflow executable or the pinned framework JAR is available.

In this mode Nextflow owns `sbatch/squeue` submission and the in-allocation
`slurm_smoke.sh` wrapper does not submit a nested job. When the wrapper is
called outside Nextflow, `FASTGATK_USE_SBATCH=1` enables its direct
`sbatch --wait` path. `SLURM_CPUS_PER_TASK`, `SLURM_MEM_PER_NODE`/
`SLURM_MEM_PER_CPU` and optional `FASTGATK_SLURM_PARTITION`,
`FASTGATK_SLURM_ACCOUNT`, `FASTGATK_SLURM_TIME` are translated at this
boundary; without it the wrapper simulates the SLURM resource environment
locally. The native process always caps its batch/threads to the allocation and
returns `RESOURCE_EXHAUSTED` when the safe lower bound cannot be met.

The generic Nextflow compatibility process declares its published files from
the selected `--output` stem (not a fixed `result*` glob), so existing
processes may retain custom output names. `nextflow.config` additionally
provides optional `container`, `cpus`, `memory`, and `slurm_gpu` settings. The
GPU profile emits a scheduler `--gpus` request only; it does not claim a
working CUDA/Kokkos device. `slurm_smoke.sh` accepts
`FASTGATK_GPU_COUNT/TYPE` and an explicit `FASTGATK_CUDA_VISIBLE_DEVICES` for
site wrappers. The fake-sbatch verifier covers this mapping and rejects
invalid counts; real container/GPU runtime and multi-node/preemption remain
deployment gates.

The local GATK oracle check runs the same interval through GATK 4.6.2.0 and the
native prototype. It verifies a shared sentinel call and reports call-set overlap;
it also validates GATK/native `-ERC GVCF` interval coverage, non-overlap, GQ/PL
shape, GQ-band headers and the shared sentinel candidate. It deliberately asserts
that the current prototype is not bit-identical:

```bash
FASTGATK_REQUIRE_GATK_ORACLE=1 \
  python3 fastgatk-native/scripts/verify_gatk_oracle.py
```

`verify_gatk_genotype_gvcf.py` additionally reports the exact INFO/FORMAT key-set
and record-text differences. With `--gatk-compatible-annotations`, the bundled
chr17 `17:69000-70000` window has three raw-text-exact variant rows, 1,001
reference-only rows, and exact semantic header definitions (including rank-sum
formatting and uncovered-site `./.` behavior); the default native profile
deliberately retains `RAW_MQandDP/RCQ/RCP/SB` diagnostics. Remaining full-file
with explicit `--add-output-vcf-command-line false`, ordinary VCF headers and
records are byte-identical; the default native profile retains truthful
provenance and GVCF reference-confidence diagnostics, so full GVCF identity is
still not claimed. Remaining differences are RCM block semantics and
provenance, not an AVX/fp32-vs-fp64 mismatch.

The DSL2, profile, and SLURM examples are in
`fastgatk-native/workflow/nextflow_smoke.nf`,
`fastgatk-native/workflow/nextflow.config`, and
`fastgatk-native/workflow/slurm_smoke.sh`.

The BQSR regression runs the native report/apply loop and, when the bundled
Java 17/GATK 4.6.2.0 oracle is present, compares all four report tables and
consumes the Java report directly:

~~~text
FASTGATK_REQUIRE_GATK_ORACLE=1 python3 fastgatk-native/scripts/verify_bqsr.py
~~~

CTest also registers the focused `fastgatk-bqsr-gatk-oracle` slice.  It uses
the pinned GATK report as an ApplyBQSR input and compares all four report
tables plus 493 Java-`PrintReads` records, including the
`--use-original-qualities`/`OQ` path, so this regression can run independently
of the broader native BQSR option matrix on both OpenMP and Serial builds.

`fastgatk-bqsr-context-size-gatk-oracle` pins ownership of the mismatch
ContextCovariate dimension. BaseRecalibrator serializes
`mismatches_context_size` in the GATKReport `Arguments` table, and ApplyBQSR
reconstructs contexts with that width instead of assuming the default 2-mer.
The deterministic size-3 fixture compares `Arguments`, `Quantized`, and all
three recalibration tables, then requires 200/200 Java-decoded ApplyBQSR
records to match. The fixture deliberately produces a context-specific QUAL
change, so ignoring the report argument cannot pass by posterior rounding.

`fastgatk-bqsr-report-roundtrip-gatk-oracle` covers the producer-side direct
replacement boundary: Java GATK 4.6.2.0 `ApplyBQSR` consumes the report emitted
by native `BaseRecalibrator`. The native writer emits GATK v1.1 fixed-width
columns (including integer `Quantized` metadata), and Java/native outputs
compare as 493 exact records through `PrintReads`. This is report
serialization compatibility, not a claim of complete BQSR-model parity.

`fastgatk-bqsr-indel-gatk-oracle` is the focused BI/BD event-table gate. It
materializes the pinned reads with deterministic per-base `BI`/`BD` SAM Z tags,
then compares `Quantized`, `RecalTable0`, `RecalTable1`, and `RecalTable2` with
GATK 4.6.2.0. The gate fixes the I/D reported-quality rows (including their
quality-keyed empirical priors), checks the four insertion/deletion qualities,
and requires truncated non-empty tags to fail closed. Missing tags retain
GATK's Q45 (or configured default) fallback. OpenMP and Serial both pass this
bounded event-model slice; it does not claim full BQSR-model parity on large or
cloud inputs.

`fastgatk-bqsr-preserve-gatk-oracle` isolates ApplyBQSR's
`--preserve-qscores-less-than` boundary. It feeds a deterministic mixture of
Q0..Q40 qualities through the same GATK 4.6.2.0 report and compares all Java
and native records after pinned `PrintReads`; qualities below a non-default Q10
threshold are explicitly asserted unchanged. The C++ decode/encode stages use
one shared predicate, so the low-quality bases bypass both Bayesian recalibration
and quantization. OpenMP and Serial both pass this bounded slice.

`fastgatk-apply-bqsr-alias-gatk-oracle` covers the direct-replacement report
selector alias: GATK 4.6.2.0 accepts `-bqsr` as the short spelling of
`--bqsr-recal-file`, while the native parser previously rejected it. The gate
uses a pinned Java report and compares long/short spellings for both Java and
native ApplyBQSR, including 493 Java-`PrintReads` records. OpenMP and Serial
both pass; this is a CLI/report-selection boundary, not a claim of complete
BQSR-model, cloud, or MD5 parity.

`fastgatk-bqsr-read-filter-gatk-oracle` covers BaseRecalibrator's generic GATK
read-filter plugin controls. The native tool resolves the seven 4.6.2.0
defaults (`MappingQualityNotZeroReadFilter`, `MappingQualityAvailableReadFilter`,
`MappedReadFilter`, `NotSecondaryAlignmentReadFilter`, `NotDuplicateReadFilter`,
`PassesVendorQualityCheckReadFilter`, and `WellformedReadFilter`) and supports
repeatable `-RF/--read-filter`, `-DF/--disable-read-filter`, plus optional-boolean
`--disable-tool-default-read-filters`. Resolution is deterministic: default set
(or empty set), `-DF` removal, then `-RF` addition; unsupported and inverted
filters fail closed. A synthetic SAM with MAPQ/flag cases is compared against
pinned GATK 4.6.2.0 across four scenarios, including all four report tables and
the effective 2/3/6/1 passing-read cardinalities. OpenMP and Serial both pass;
this is a filter/plugin boundary, not a claim of complete BQSR-model parity.

`fastgatk-bqsr-cram-gatk-oracle` is a separate CRAM gate.  It consumes the
bundled `dictFix.cram`/CRAI with known sites, compares the four
BaseRecalibrator tables, writes CRAM from both implementations, and reopens
both outputs through pinned HTSJDK.  ApplyBQSR selects CRAM 3.0 at the HTSlib
writer boundary (GATK 4.6.2.0 rejects CRAM 3.1) and records `.crai`; Java's
generic index option may produce `.bai`.  The gate compares decoded SAM records,
not container bytes, and does not claim full BQSR-model parity.

The native BQSR prototype is a separate report/apply loop:

~~~text
fastgatk-native/build/fastgatk-bqsr -I reads.bam -R reference.fa \
  --known-sites sites.vcf --known-sites sites-2.vcf -O recal.tsv
fastgatk-native/build/fastgatk-bqsr -I reads.bam -R reference.fa \
  --batch-records 4096 --checkpoint recal.checkpoint -O recal.tsv
fastgatk-native/build/fastgatk-bqsr -I reads.bam -R reference.fa \
  --batch-records 4096 --resume-checkpoint recal.checkpoint -O recal.retry.tsv
fastgatk-native/build/fastgatk-apply-bqsr -I reads.bam -R reference.fa \
  --bqsr-recal-file recal.tsv -O recalibrated.bam
fastgatk-native/build/fastgatk-gather-bqsr-reports -I shard-1.tsv -I shard-2.tsv \
  -O merged-recal.tsv
fastgatk-native/build/fastgatk-analyze-covariates -before merged-recal.tsv \
  -csv recalibration.csv -plots recalibration.pdf \
  --output-manifest recalibration.manifest.json
fastgatk-native/build/fastgatk-genotype-gvcf -V sample.g.vcf.gz \
  -O cohort.vcf.gz --output-manifest cohort.vcf.gz.manifest.json
fastgatk-native/build/fastgatk-combine-gvcfs -V shard-1.g.vcf.gz -V shard-2.g.vcf.gz \
  -O combined.g.vcf.gz --output-manifest combined.g.vcf.gz.manifest.json
~~~

`BaseRecalibrator` emits a GATKReport v1.1-compatible report containing the
Arguments, Quantized, RecalTable0/1/2 tables, plus a
`<report>.covariates.tsv` sidecar keyed by PU (Platform Unit, falling back to
RG ID), GATK's 2-base read context, strand-aware cycle and raw quality.
Repeated `--known-sites` VCF/BCF inputs are decoded by HTSlib, mapped to the BAM
dictionary, and unioned deterministically as merged `VariantContext` REF/INFO-END
intervals (so multi-base indels are masked across their full reference span rather
than only at the record start); their ordered signature and count are recorded in
the checkpoint/manifest. BaseRecalibrator applies GATK's standard BQSR
MAPQ/mapped/secondary/duplicate/QC/Wellformed filter composition and hard-clips
CIGAR soft-clips before covariate projection. Host-side CIGAR/reference/known-sites projection feeds the reusable Kokkos
`count_bqsr_quality_kokkos` integer histogram kernel. Covariate keys are
resolved on Host and compressed to dense IDs per batch; their count/mismatch
reduction uses `count_bqsr_covariates_kokkos` with a bounded `TeamPolicy`
team-local histogram and deterministic Host league-order merge, falling back
to a `RangePolicy` atomic histogram for oversized key spaces. The device path
does not mutate a Host map. The quality count path uses the same bounded
team-local pattern, while ApplyBQSR uses a `RangePolicy` transform; selected
policies, histogram workspace bytes, prepare/execute timings and execution
spaces are recorded in the manifest. `GatherBQSRReports`
reads and writes the same report schema and merges sidecars deterministically.
It also parses Java GATK reports without a FASTGATK sidecar: `Arguments`,
`Quantized`, `RecalTable0`, `RecalTable1`, and the separate Context/Cycle rows in
`RecalTable2` are merged with GATK's expected-error reported-quality rule and
the Java `QualQuantizer` greedy merge using the first input report's
`quantizing_levels` value (including non-default levels such as 4).  Structural
recalibration-table dimensions are checked before merging, matching Java's
fail-closed behavior for incompatible covariate/indel-context reports.  The
bundled oracle compares all five tables against `gatk GatherBQSRReports` for
disjoint reports and a non-default quantization level; the manifest records
`direct_gatk_report_inputs`, `quantizing_levels`, strict Host execution, output
completeness, output sizes and wall time.
`ApplyBQSR` implements GATK's hierarchical Bayesian RecalDatum estimate
(read-group → reported-quality → context/cycle), preserves qualities below
`--preserve-qscores-less-than`, and defaults to GATK's `--quantize-quals 0`
(no dynamic quantization; `-1` consumes the report map). A positive
`--quantize-quals N` now rebuilds the N-level greedy `QualQuantizer` map from
the report's empirical histogram instead of reusing a map with a different
level count. The
`--global-qscore-prior` override is propagated through the same hierarchy. If
the optional
sidecar is absent, the same hierarchy is reconstructed from RecalTable1/2.
ApplyBQSR also accepts repeatable `--static-quantized-quals` bins and
`--round-down-quantized`, using GATK's probability-space nearest/lower-bin mapping
after recalibration; these options are mutually exclusive with dynamic
`--quantize-quals`. Missing read groups fail closed by default and can be handled
only with explicit `--allow-missing-read-group`, which bypasses model lookup and
uses the selected quantizer on the original/OQ-reset qualities.
`--use-original-qualities` follows GATK's OQ semantics: a valid `OQ:Z` FASTQ
string is decoded on Host for model lookup while the original OQ tag is
preserved in output; malformed tags fail closed. The pinned chr17 oracle
compares all 493 output records, including tags and per-base qualities, after
Java `PrintReads`. `--emit-original-quals` adds an OQ tag from the pre-transform
QUAL only when one is absent, matching the non-overwrite behavior of HTSJDK.
`--create-output-bam-index` follows GATK's optional-boolean contract: indexing
is enabled by default, a bare flag enables it, and `false`/`0` suppresses the
`.bai`/`.crai` sidecar. The selected artifact policy is recorded in the
checkpoint signature and OutputManifest, so resumed jobs cannot silently
change their published file set.
Both BQSR stages expose the same manifest contract as the other native tools:
Kokkos execution space, strict determinism, complete output/sidecar/index
checks, byte counts, and wall-clock telemetry are recorded and asserted by the
GATK oracle. Before the first decode, both stages cap the requested batch using
the conservative SLURM/cgroup Host budget; BaseRecalibrator may reduce the next
batch after observing `ReadBatch::bytes()` during a non-checkpointed stream, while
checkpoint/resume keeps stable record boundaries. `effective_batch_records` and
`adaptive_batch_reductions` make this backpressure auditable.
`AnalyzeCovariates` reads one or more Java/native GATKReport v1.1 files via
`-bqsr` (the GATK 4.6.2.0 short alias; the legacy `-BQSR` spelling remains
accepted), `-before`, or `-after`, emits GATK's intermediate CSV schema (including
the expected-error reported-quality merge and encoded Context/Cycle ordering),
and writes a dependency-free vector PDF when `-plots` is requested. The batch
accuracy transform uses the shared `HostBatch -> KernelPlan.prepare -> Kokkos
Views -> execute -> collect`/`RangePolicy` lifecycle. When
multiple reports are supplied, shared recalibration arguments are checked
before merging so incompatible before/after models fail closed. The CSV
numeric rows are bit-identical to GATK 4.6.2.0 on the pinned chr17 report; the
batch accuracy transform records execution space/policy/timing and row count in
the manifest. The manifest also records strict Host execution,
output completeness/byte counts and wall time; the reproducible file-boundary
benchmark is `fastgatk-native/scripts/benchmark_analyze_covariates.py`.
The focused `verify_analyze_covariates_bqsr_alias_gatk_oracle.py` gate pins
Java's lower-case `-bqsr`/`--bqsr-recal-file` equivalence, native lower-case /
legacy uppercase equivalence, and rejection of Barclay's embedded `-bqsr=`
form on both OpenMP and Serial builds.
For HaplotypeCaller `-ERC GVCF`, the reference-confidence path exposes
`--indel-size-to-eliminate-in-ref-model` (default `10`) and evaluates each
observation's read/reference suffix against insertion/deletion shifts before
adding it to indel-informative depth; the chosen bound is recorded in the
Manifest. The ref-vs-any pileup likelihood stage is now the reusable
`calculate_reference_confidence_kokkos` kernel (OpenMP/CUDA/HIP-capable through
Kokkos); Host code performs the deterministic locus reduction and retains the
GATK indel model. Kernel prepare/execute time and execution space are emitted
in HC telemetry/OutputManifest. Candidate-site `<NON_REF>` uses the
qualified-concrete median matrix; complete assembly-haplotype and cross-sample
posterior RCM semantics remain explicit fallback scope.
普通 `-ERC GVCF` 在最终 genotype selection 后重建 reference-confidence blocks：未成
call 的 SNP/MNP candidate 回收到 `<NON_REF>` block，必须保留 reference span 的
indel candidate 仍作为 candidate record 输出；零支持的 assembly-only indel 在
`--min-alt-support` 支持度门下回收到 block，避免产生 GATK 不会写出的空候选行。输出
字段采用 GATK-shaped `GT:DP:GQ:MIN_DP:PL`（reference block）和
`GT:AD:DP:GQ:PL:SB`（candidate），candidate 的 QUAL/FILTER、
`BaseQRankSum/DP/ExcessHet/MLEAC/MLEAF/MQRankSum/RAW_MQandDP/ReadPosRankSum`
顺序及未封顶的 Number=G PL 与 GATK 保持一致；普通 block 使用 `FILTER=.`、`INFO=END`；
BP_RESOLUTION 保持逐碱基 `GT:AD:DP:GQ:PL`/`INFO=.`。这保证 candidate-site 集合和
block 不重叠，并通过 manifest 的 `gvcf_candidate_rebuild_after_genotyping`、
`gvcf_standard_fields`、`gvcf_candidate_standard_fields` 显式记录；RCM 数值/边界仍
以 oracle 门禁报告，尚不宣称 bit-identical。
非二倍体 RCM 也已与 GATK 4.6.2.0 对齐：indel 模型的 least-confidence 比较使用
完整 `ploidy+1` likelihood vector，polyploid SNP likelihood 在转换 PL 前按
hom-ref cap；candidate-site `MLEAC/MLEAF` 对 concrete + `<NON_REF>` 的完整
Number=G 矩阵运行 bounded EM。固定 triploid chr17 oracle 对 69000 reference
block 和 69067 candidate 的 sample/ML 字段 exact；更大范围的 GATK block
合并边界与完整 reference-block prior corpus 仍是后续扩展项。
`verify_hc_bp_resolution_gatk_oracle.py` 另外对固定 chr17 区间的 21 个逐碱基
记录执行 GATK 4.6.2.0 exact GT/AD/DP/GQ/PL 对照，避免把 BP_RESOLUTION
误认为只完成了 writer 层。
自定义 GVCF 分箱也已覆盖：`-GQB/--GVCF-GQ-bands` 为可重复的上界参数，严格要求
正数、递增且不超过 100，省略终止 100 时自动补齐；例如 `-GQB 10 -GQB 20
-GQB 50` 生成 `0-10/10-20/20-50/50-100`，reference-block 合并和
`##GVCFBlock` header 使用同一分箱。`verify_hc_gq_bands_gatk_oracle.py` 在
OpenMP/Serial 上与 GATK 4.6.2.0 的 header 及记录均 exact。
同一 writer 还支持 GATK optional Boolean `--floor-blocks`：它不改变 RCM 或
candidate-site 选择，而是把普通 `REF/<NON_REF>` block 的 GQ 向下取到实际 GQ
分箱下界，并将该 block 序列化为 `GT:DP:GQ`（不写 `MIN_DP/PL`）。candidate rows 和
`BP_RESOLUTION` 保持既有字段。上述 pinned GQ-band oracle 同时比较此模式的 Java/native
header、records、reference-block FORMAT 和 OutputManifest，OpenMP/Serial 均通过。
`benchmark_bqsr.py` adds a reproducible file-boundary benchmark for both
BaseRecalibrator and ApplyBQSR; it reports p50/p95 wall time, report/BAM/index
bytes, and the Kokkos prepare/execute policy. `--compute-indel-bqsr-tables` adds GATK-compatible I/D event tables with the
independent `--indels-context-size`, configurable default qualities, CIGAR
event errors and the GATK four-base indel-cycle cushion. The bundled chr17
oracle compares all four tables for this mode as well; long-read/large-cohort
corpora and object-store/cloud staging remain explicit follow-up scope. For the
same pinned fixture, `verify_bqsr.py` additionally converts GATK/native ApplyBQSR
outputs through Java `PrintReads` and compares all stable SAM record fields
(flag, CIGAR, sequence, and every quality character): 493/493 records are
bit-identical, not merely equal in aggregate quality sum. The
optional BaseRecalibrator checkpoint is an atomic, batch-boundary text state
file: `--resume-checkpoint` validates input size/mtime, reference/known-sites,
interval and all relevant BQSR options before skipping completed batches. The
report and covariate sidecar are regenerated byte-for-byte. ApplyBQSR now
records the last fully flushed batch, validates the input/report and all
quality-model options. On resume it requires the durable BAM/CRAM output
prefix, copies the completed records, transforms only the remaining input,
then atomically publishes the replacement and rebuilds its index; this keeps
recovery portable without depending on codec-specific append behavior.
Checkpoint state, prefix reuse, and restart status are included in
OutputManifest telemetry.

## Dispatcher and tool registry

The P0 compatibility dispatcher is
`fastgatk-native/dispatcher/fastgatk`. It reads the versioned
`fastgatk-native/dispatcher/tool_registry.json`, expands GATK-style `@args`
files, supports `HaplotypeCaller`, `BaseRecalibrator`, `ApplyBQSR`,
`GatherBQSRReports`, `AnalyzeCovariates`, `Mutect2`,
`FilterMutectCalls`, `CollectF1R2Counts`, `LearnReadOrientationModel`, `CollectReadCounts`, `DenoiseReadCounts`, `CallCopyRatioSegments`, `CollectAllelicCounts`, `ModelSegments`, `GatherTranches`, `VariantRecalibrator`, `ApplyVQSR`, `GetPileupSummaries`, `CalculateContamination`, `GatherPileupSummaries`, `ReblockGVCF`, `SelectVariants`, `VariantsToTable`, `VariantEval`, `ValidateVariants`, `GatherVcfs`, `LeftAlignAndTrimVariants`, `VariantFiltration`, `SortSam`, `MarkDuplicates`, `CombineGVCFs`, and the scoped `GenotypeGVCFs` contract-compatible materialization path with
their GATK parameter aliases. `GenomicsDBImport` is registered as an explicit
`adapter` tool: a C++ Host layer adapts batch size, reader threads and
`SLURM_TMPDIR`, then delegates sparse workspace creation to the configured
GATK GenomicsDB backend and validates the resulting workspace/manifest. The
explicit `--fastgatk-native-workspace` (or
`FASTGATK_GENOMICSDB_BACKEND=native`) mode instead writes a portable
`fastgatk-inputs.tsv` + `fastgatk-workspace.json` sparse index for native
GenotypeGVCFs; it also materializes local VCF/GVCF inputs under
`native-inputs/` and writes `fastgatk-native-inputs.tsv`, so the workspace can
be consumed after the original files are moved. Its manifest is marked
non-TileDB/non-bit-identical.
The native adapter preserves the explicit `-V`/sample-map callset order with
stable deduplication; it does not sort paths by filename. A reopened native
`gendb://` query therefore keeps the same sample columns as the import order.
The GenomicsDB contract verifier covers a reversed/duplicate input list and a
chr2-only shard queried on chr1 (the non-overlapping sample remains in the
header with a missing genotype). This is a native sparse-index boundary only;
external TileDB/GenomicsDB workspaces continue to use the configured GATK
backend or the version-coupled legacy bridge.
Every workspace write (external backend, native materialization, and manifest)
is guarded by an atomic sibling-directory lock named
`<workspace>.fastgatk.lock`; concurrent retries fail closed with
`RESOURCE_EXHAUSTED`, and normal/error paths release the lock automatically.
The owner PID/thread and lock path are recorded in the manifest, so stale locks
can be inspected before manual cleanup.
Native materialization is restart-safe: each copied input is recorded in an
atomic `fastgatk-native-workspace.checkpoint` with size/mtime fingerprints.
If a SLURM/Nextflow task is preempted, rerun with
`--resume-native-workspace`; completed copies are validated and reused, changed
inputs fail closed, and the checkpoint is removed only after the record index
and workspace metadata publish successfully. Checkpoint reuse counts and the
resume flag are emitted in both the summary and OutputManifest telemetry.
Native sparse workspaces also accept `--genomicsdb-update-workspace-path` when
the existing workspace has a valid `fastgatk-inputs.tsv`: the adapter merges
the prior and new canonical inputs, atomically rebuilds the materialized VCF
and record index, and carries forward sample-name-map rewrites. This is a
deterministic no-Java incremental boundary for native GenotypeGVCFs, explicitly
reported as `native_incremental_update`; it is not a TileDB append and opaque
or partial workspaces fail closed.
Native sparse imports also apply `-L/--intervals` before materialization using
record-span overlap (including interval-list/BED selectors and UNION/
INTERSECTION). A reference-confidence block is retained in full whenever its
span overlaps the imported interval, and rewritten inputs do not inherit stale
source VCF/CSI indexes. The native checkpoint records the interval selector
set and incremental updates restore the initial selection, matching GATK's
update boundary. The pinned
`verify_genomicsdb_import_native_interval_gatk_oracle.py` compares two
GATK 4.6.2.0 bridge-exported boundaries on both OpenMP and Serial; complete
TileDB storage/query compatibility remains with the external backend.
adapter preserves release-specific arguments; `--fallback` remains available
when the adapter or backend is unavailable. The dispatcher provides
`--help`, `--version`, `--list`, `--dry-run`, and explicit `--fallback` routing.
Unknown parameters remain structured errors and are never silently discarded.
An unregistered long-tail tool is fail-closed by default, but with explicit
`--fallback` it is forwarded unchanged as `gatk <Tool> ...`; this keeps the
fallback-first inventory usable while a native contract is being developed.
The fallback launcher can be overridden with `FASTGATK_GATK_BINARY` (the tool
name is still inserted after the launcher).
Deployments with a nonstandard artifact directory can set
`FASTGATK_NATIVE_BUILD=/path/to/fastgatk-native/build`; the dispatcher resolves
the registry binary suffix there, and automatically probes the sibling
`build-serial` artifact when the conventional OpenMP binary is absent.

The launcher contract is explicit in the registry's `launcher_contract` object
and is checked by `verify_launcher_contract.py`. `--java-options` (including
the legacy `-java-options` spelling) is consumed before or after the tool and
is forwarded only when an explicit Java fallback is selected; native binaries
never receive JVM flags. `--gatk-config-file PATH` is likewise accepted in
either position, but routes to Java because a GATK properties file can alter
codec/cloud/tool defaults outside the native registry. Nested `@args` files
are expanded recursively up to depth 16, `@@` escapes a literal `@`, and
relative paths are resolved against the launcher working directory. A leading
`--` separates launcher options from the tool. A post-tool `--` group is
preserved for explicit Java fallback (including Spark arguments); native
execution fails closed instead of discarding that group. Unknown launcher or
tool options remain fail-closed unless `--fallback` is explicit.

`Mutect2` now has a native prototype (`fastgatk-mutect2`) that reuses the HC
kernel path for tumor/normal inputs and evaluates a Kokkos per-read REF/ALT
mixture grid for deterministic TLOD and best-AF evidence. A second Kokkos
posterior stage now combines tumor/normal evidence with somatic/germline/
artifact priors, F1R2/R1F2 orientation evidence and contamination correction;
VCF INFO includes `PSOMATIC`, `PGERMLINE`, `PARTIFACT`, `OBP` and corrected
`CONTAM` AF. It also emits deterministic `AF`/`AD` annotations, orientation
counts and F1R2 sidecar, and writes BGZF/Tabix VCF. Concrete sibling ALTs at
one `(contig,position,REF)` are materialized as one multi-ALT VariantBlock:
`TLOD/AF/F1R2/R1F2` and posterior INFO vectors use `Number=A`, while tumor
`AD/FAD` retain REF plus every concrete ALT in the same order. Tumor records expose the
GATK-shaped `GT:AD:AF:DP:F1R2:F2R1:FAD:SB` sample fields, and local shared-read
evidence adds GATK-shaped `0|1:...:PGT:PID:PS:SB` phasing fields when a phase
link is provable; normal records use explicit missing orientation values when the
native normal-read batch is not retained. Release-specific GATK
prior calibration, read-orientation model fitting and full FilterMutectCalls
joint inference remain explicit fallback scope.
The VCF now also carries the common Mutect2-shaped `AS_SB_TABLE`, `DP`,
`ECNT`, `ECNTH`, `MBQ`, `MFRL`, `MMQ`, `MPOS` and `POPAF` INFO fields.  The
strand/depth and per-allele quality summaries are computed deterministically
from projected Host reads; release-specific event-count and germline-resource
calibration remain explicitly approximate until their corresponding oracle
gates are added. `--af-of-alleles-not-in-resource` (and the
`--population-af` alias) controls the `POPAF` fallback and defaults to GATK's
`5e-8`.
`--minimum-allele-fraction` (or GATK's `-min-AF` alias) now reaches the
Kokkos SomaticLikelihoodsEngine prior. The default `0` keeps the flat `[1,1]`
Dirichlet prior; positive values use GATK's `1-log(2)/log(minAF)` ALT
pseudocount and are recorded in stats/manifest.
Tumor and normal evidence retain independent fragment widths, matching GATK's
separate `AlleleLikelihoods` collections: the Kokkos posterior infers the
normal row width from its own candidate-major matrix instead of discarding the
normal whenever tumor and normal coverage differ. The pinned
`fastgatk-somatic-posterior-normal-count-gatk-oracle` uses four tumor and two
normal evidence units and passes on OpenMP and Serial; the broader
release-specific somatic posterior remains explicit fallback scope.
`--mitochondria-mode` is accepted using GATK's optional-boolean syntax (bare,
`=true|false`, or a following boolean token). When enabled, the native path
applies the GATK mitochondrial defaults `initial-tumor-lod=0`,
`initial-tumor-lod=0`, `tumor-lod-to-emit=0`, `population-af=4e-3`,
`graph-pruning-log-odds-threshold=-4*ln(10)` and
`recover-all-dangling-branches`; explicit values still take precedence. The
effective mode/defaults are written to Mutect2 stats and OutputManifest. This
is a compatibility boundary for the current prototype: full GATK mitochondrial
`OriginalAlignment` annotation and release-specific somatic posterior
calibration remain explicit fallback work.
For loci with multiple concrete ALTs, the Host groups the post-filter fragment
likelihood rows and dispatches the shared Kokkos
`calculate_somatic_multiallelic_likelihood_kokkos` kernel. It follows GATK's
all-alleles-vs-without-each-ALT Dirichlet evidence calculation rather than
independently normalizing each ALT; singleton loci keep the biallelic path for
regression parity. The public kernel also accepts a final `<NON_REF>` row for
callers that explicitly run reference-confidence mode; Mutect2's default
emission path does not add that row.
Before the somatic likelihood reduction, the Host boundary preserves the
post-filter read-name map and assigns stable fragment/haplotype cells. Kokkos
then applies GATK `AlleleLikelihoods.groupEvidence` arithmetic: reads from one
fragment are summed per haplotype in encounter order (including propagation of
`-Infinity` zero-probability cells) and only then max-marginalized to REF/ALT;
missing names fail closed to singleton evidence. The pinned
`fastgatk-fragment-aggregation-gatk-oracle` checks the exact cell and allele
bits against GATK 4.6.2.0 and is registered in both Serial and OpenMP builds.
For a locus with two or more concrete ALTs, a haplotype carrying a
sibling ALT is marked unknown rather than reused as REF for the current ALT.
`somatic_evidence_groups` and
`somatic_evidence_grouping` are recorded in stats/manifest and included in the
Mutect2 benchmark output.
The variational kernel's positive-argument digamma now follows the pinned
Commons Math `Gamma.digamma` branch points and Bernoulli-term order.  Local
combination paths marked state=2 are requested for the complete per-read
normalization baseline of that local block but are not projected into the candidate allele row;
this preserves GATK's normalize-before-marginalize ordering without treating an
unrelated local path as REF evidence.  Fragment interval selection also keeps
the Java `Fragment.createAndAvoidFailure` limit for malformed names while still
summing all group likelihoods as GATK does.
The same boundary now retains compact post-filter tid/start/end metadata and
keeps a grouped fragment only when its span overlaps the merged candidate
expanded by GATK's default `informativeReadOverlapMargin=2`; both mates remain
summed when either mate overlaps. Mutect2's shared overlap-quality stage also
preserves conflicting mate qualities, matching its Java
`setConflictingToZero=false` call (HC keeps conflict-to-zero). On the pinned
chr17 broad oracle, the complete ordered VCF (including TLOD) is now exact to
GATK 4.6.2.0 after excluding only execution provenance. This is fixture-scoped
evidence, not a claim of parity for every release-specific input combination.
The Java/native integration check is registered as
`fastgatk-mutect2-gatk-oracle` in CTest and runs the pinned GATK 4.6.2.0 jar
when its corpus is available. The independent
`fastgatk-somatic-likelihood-gatk-oracle` matrix check isolates
`SomaticLikelihoodsEngine`: its Python reference mirrors the release's
Dirichlet update, convergence threshold, log-evidence entropy cutoffs,
all-ALT-vs-without-ALT calculation and the non-flat minAF prior, and passes on
both OpenMP and Serial.
This separates true TLOD-input differences (assembly/PairHMM) from a
somatic-evidence formula regression. GATK's
`--minimum-allele-fraction`/`--min-AF` is accepted by the native CLI and uses
the corresponding non-flat ALT pseudocount in the Kokkos kernel; the value is
recorded in stats/manifest. Other release-specific posterior calibration
remains explicit fallback scope.
The Kokkos SW score path pads the final SIMD request group when a batch is not
a multiple of the selected `Kokkos::Experimental::simd` width. Inactive lanes
are discarded after execution, so ragged production batches cannot read beyond
the lane `View` while preserving the Java score/tie contract.
When `--f1r2-tar-gz` ends in `.tar.gz`, the sidecar is a standard three-member
POSIX archive (`.ref_histogram`, `.alt_histogram`, `.alt_table`) consumable by
`LearnReadOrientationModel`; an explicit `.tsv` path remains the legacy fastgatk
sidecar. The output manifest records which format was selected.
The CLI also accepts GATK's `--tumor-sample` and `--normal-sample` aliases;
HTSlib resolves each selected name through `@RG ID→SM` and filters reads before
Host batching, assembly, PairHMM, fragment grouping and F1R2 counting.  When
tumor is omitted, the first sample in the tumor input header is selected and
enforced; a tumor header with no `@RG SM` fails closed before decoding, matching
GATK's `samples cannot be empty` contract, and unknown names fail closed.  The
chosen names and filtered read counts are recorded in stats/manifest.  A standard
single-input tumor/normal command (`-I multi-sample.bam --tumor-sample T
--normal-sample N`) reuses the same HTSlib stream for both filtered ReadBatch
objects; `--normal-input` remains available for separate tumor/normal files.
Aggregate mode also accepts repeatable `-I` shards, validates their sequence
dictionaries/assembly tags, skips shards without the selected sample, and merges
the remaining ReadBatch objects with a stable coordinate merge (argv order is
the tie-break). Multi-shard
`--stream-by-region` now probes every index, validates the merged dictionary,
and performs a stable coordinate merge for each bounded core/halo tile; an
unindexed shard still fails closed instead of being silently dropped. In
reference-confidence mode, streaming invokes the same somatic RCM writer as
aggregate traversal and assigns each complete GVCF block to the tile containing
its start, so a block crossing a core boundary is serialized exactly once. The
pinned chr17 EventMap oracle verifies aggregate, a one-transaction direct
interval, and 400 bp streamed execution in both GVCF and BP-resolution modes;
all are exact to GATK 4.6.2.0 after execution provenance is excluded.
Somatic depth aggregation updates both the candidate telemetry view and the
materialized VCF-call view before serialization, so FORMAT `DP/AD` cannot retain
the pre-grouped annotation depth. The pinned 4.6.2.0 fixture reports exact
ordered records, schema, site set, GT, DP, AD, and TLOD (6/6) on both OpenMP
and Serial.
Repeated `-L`/`--intervals` selectors also accept GATK's
`-isr/--interval-set-rule UNION|INTERSECTION`. The selected rule is applied
by the shared HTSlib Host interval engine in both aggregate and indexed
`--stream-by-region` modes and is recorded in stats/OutputManifest.
The common `--sites-only-vcf-output true|false` switch is also supported: it
keeps the same caller, posterior and F1R2 computation but makes the final VCF
site-only (8 columns, no FORMAT/sample fields).  BGZF/Tabix, stats and
manifest outputs remain enabled and record the selected serialization mode.
GATK's repeatable `-normal`/`--normal-sample` selectors are supported as well:
each `@RG ID→SM` sample is decoded into its own Host ReadBatch, all selected
normal evidence is combined for somatic inference, and the VCF retains one
independent normal genotype column per selector in argv order.  The aggregated
`normal_samples` list and read count are recorded in stats/manifest; the
single-normal and `--normal-input` forms remain backward compatible.
Mutect2 also accepts GATK's optional boolean
`--create-output-variant-index[=true|false]` (default `true`).  For compressed
VCF output, setting it to `false` suppresses the sibling `.tbi` without changing
caller/posterior, VCF records, or F1R2/stats; stats and OutputManifest record both
the requested switch and the actual `vcf_index` state.  OpenMP/Serial verifier,
targeted CTest, and dispatcher registry contracts cover this boundary.  For an
uncompressed VCF, the same enabled switch writes an HTSJDK/Tribble LinearIndex v3
sibling `.idx`, which is read by a GATK SelectVariants interval query; disabling
the switch suppresses both `.tbi` and `.idx` as appropriate.
`CombineGVCFs` has the same bounded host/Kokkos merge boundary: the pinned
two-sample oracle checks indexed inputs, concrete-ALT/sample union,
reference-block coordinates, output indexing, and full sample FORMAT payload.
The default path is GATK-compatible for ALT order, no-call GT, GQ/AD/PL
projection, and FORMAT order; `--fastgatk-materialize-genotypes` is an
explicit native GT/GQ experiment mode.
Run `python3 fastgatk-native/scripts/benchmark_mutect2.py` for a repeatable
file-boundary timing that includes HTSlib decode, Kokkos calling, BGZF/Tabix and
standard F1R2 tar serialization.
The shared graph policy also exposes bounded `--kmer-size`,
`--min-kmer-count`, `--min-pruning`, `--max-haplotype-paths` (or GATK's
`--max-num-haplotypes-in-population`) and `--max-haplotype-depth`
controls for cluster memory budgets, plus
`--max-haplotype-combination-alleles` (1..16) for the complete local
`2^N` haplotype combination cap; effective values are recorded in manifest
telemetry. Above the cap, the caller uses deterministic coordinate blocks,
retains each block's complete `2^N` combinations within a bounded
candidate-to-haplotype link budget, shares the reference haplotype across
blocks, and deduplicates identical read×haplotype PairHMM/SW requests.
Both HC and Mutect2 accept `--max-mnp-distance`; positive values merge
phased, same-length graph substitutions from one ungapped SW segment and
record `graph_mnp_candidates`, while the default zero preserves independent
SNP output.
`CombineGVCFs` now has a native site-level shard-merge prototype with
concrete ALT union (keeping `<NON_REF>` last), duplicate removal, strict
REF-only `<NON_REF>` block coalescing, and deterministic PL/AD remapping
for disjoint sample inputs; blocks merge on true coordinate adjacency even when
their depth metadata differs, while disjoint sample FORMAT calls are retained.
The default output preserves GATK's no-call GT and input GQ, while
`--call-genotypes` enables GATK-compatible arbitrary-ploidy GT/GQ derivation;
the historical `--fastgatk-materialize-genotypes` spelling remains available
as an explicitly native experimental mode.  Standard
reference-band controls are also accepted: `--convert-to-base-pair-resolution`
and `--break-bands-at-multiples-of N`; split boundaries are kept through the
merge stage and are never silently coalesced again.  If `-R` is supplied,
interior split REF bases are fetched from the FASTA index when available.
The manifest records the selected format mode, strict determinism, execution
space, output/index bytes, wall time, and Kokkos remap telemetry. The native
contract exercises both modes, while `verify_combine_gvcfs_gatk_oracle.py`
requires full sample-payload bit identity for the default path.
Repeated `-L`/`--intervals`/`--region` selectors use GATK's
`--interval-set-rule UNION` (default) or `INTERSECTION` semantics before site
grouping and block coalescing; the selected rule is recorded in the summary
and OutputManifest and is shared by aggregate and `--stream-merge` traversal.
The same reference-band controls are supported by `--stream-merge`: a cursor
retains only the original block and its next segment, and reports
`lazy_reference_blocks`/`lazy_reference_segments` instead of pre-expanding a
large base-pair-resolution block.
Full joint-genotyping/GenomicsDB semantics
remain explicit fallback. `FilterMutectCalls` now has a native TLOD/germline
filter prototype with stats and BGZF/Tabix output; contamination uses a
GATK-shaped AD/POPAF posterior (per tumour sample, ALT-depth weighted median)
through a Kokkos batch kernel and honours `--max-contamination-probability`;
when POPAF/AD is absent it retains the legacy approximate AF contamination
floor. Minimum allele fraction, minimum reads per strand, unique ALT-read
support, concrete-ALT-count, and F1R2/R1F2 orientation-balance filters are
exposed. GATK-shaped contamination tables are parsed per tumor sample (with a
global estimate retained for samples not present in the table). The native
output also exposes this posterior as `CONTQ` and
per-ALT `AS_FilterStatus`; the pinned GATK 4.6 minimal allele-filter path does
not materialize `CONTQ` consistently, so the native annotation is marked as an
explicit extension in the oracle.
Release-specific model fitting and full joint filtering still use explicit
fallback. The common `--sites-only-vcf-output true|false` writer option is
also supported: all AD/F1R2 and filter calculations run before HTSlib strips
FORMAT/sample fields, and the final output is an 8-column site-only VCF while
retaining INFO/FILTER, index, stats, and manifest sidecars. Clustered-event
annotation filtering is covered below.
The native path also implements a bounded GATK `GermlineFilter` posterior
surface: NLOD/POPAF/tumor AD are consumed by a Kokkos likelihood kernel and the
result is emitted as Number=A `GERMQ` using GATK-compatible phred bounds
`[1,93]`. Missing POPAF fails closed, and the optional `PGERMLINE` annotation is
not substituted for this calculation. A pinned GATK 4.6.2.0 oracle passes on
both OpenMP and Serial; release-calibrated full joint posterior/filter learning
remains explicit fallback.
The empirical `OPTIMAL_F_SCORE` pre-pass now reuses the same Kokkos AD/POPAF
contamination batch used by final filtering and merges its per-ALT posterior
with GermlineFilter as Java `ErrorProbabilities` does (`NON_SOMATIC=max(...)`).
The pinned `fastgatk-filter-mutect-calls-contamination-joint-oracle` fixes the
10-record AD-gradient contamination FILTER boundary and verifies that the
learned threshold is finite and no longer the pre-fix `1.0` sentinel on both
OpenMP and Serial; remaining release-specific somatic calibration is explicit
fallback.
When `--orientation-bias-artifact-priors` is supplied, the Kokkos
F1R2/F2R1 weighted-median `ReadOrientationFilter` posterior is also fed into
the empirical `OPTIMAL_F_SCORE` ARTIFACT reduction and the final ALT joint
posterior. The pinned
`fastgatk-filter-mutect-calls-orientation-joint-oracle` checks the orientation
FILTER pattern, GATK threshold-zero boundary and no-prior baseline on OpenMP
and Serial; complete release calibration remains explicit fallback.
`FilterMutectCalls` also accepts median-evidence thresholds
(`--min-median-base-quality`, `--min-median-mapping-quality`,
`--min-median-read-position`) from `MBQ`/`MMQ`/`MPOS` INFO annotations;
their FILTER names, stats counters and manifest compatibility flags are
deterministic. It can also consume Mutect2's `PSOMATIC`/`PGERMLINE`/`PARTIFACT`
annotations through `--min-somatic-probability`,
`--max-germline-probability` and `--max-artifact-probability`; enabled
posterior thresholds fail closed when annotations are missing. The
`--tumor-sample` selector makes AD-based hard filters independent of sample
column order (default is a `TUMOR` sample when present, otherwise the first
sample); the selected sample is recorded in stats and the OutputManifest. The
clustered-event hard filter now consumes Mutect2's
`ECNT`/`ECNTH` annotations with GATK-compatible defaults (`3`/`2`) and the
`--max-events-in-region`/`--max-events-in-haplotype` overrides; missing event
annotations fail closed (no filter) and are reported in the stats/manifest.
Optional `MFRL` fragment-length and `NCount`/tumor-AD ratio hard filters are
available through `--max-median-fragment-length-difference` and `--max-n-ratio`;
they emit GATK's `fragment` and `n_ratio` FILTER names and record thresholds and
counters without treating missing annotations as fabricated failures. The
`--max-n-ratio` default is GATK's `+Infinity` (an explicit negative value keeps
the native disable override); records with zero summed ALT AD are skipped just
like Java's `NRatioFilter`.
`--mitochondria-mode` and `--microbial-mode` are also accepted (as bare
switches or optional `=true|false`/next-token boolean values). They disable
GATK's genomic-only clustered/multiallelic/fragment/haplotype filters; microbial
mode retains polymerase-slippage filtering and uses implicit MMQ=20, while
mitochondrial mode suppresses slippage as well and applies GATK's default
log-SNV/log-indel priors of `-2.5*ln(10)`/`-3.75*ln(10)` unless explicitly
overridden. The selected mode, effective priors and disabled filter set are
recorded in the OutputManifest and stats telemetry.
When scalar orientation fields are absent, the strand-support filter also
consumes GATK's allele-specific `AS_SB_TABLE` ALT pair; and
`AS_UNIQ_ALT_READ_COUNT` takes precedence over total AD for
`--unique-alt-read-count`, matching the Rust/GATK allele-specific boundary.
The native boundary now emits `INFO/AS_FilterStatus` (`Number=A,Type=String`)
with one Rust/GATK hard-filter list per ALT.  `MBQ`/`MMQ` skip the reference
entry, `MPOS` is ALT-only, long insertions (>=3 bases) use the reference MMQ
while deletions do not, strict strand checks zero counts in `AS_SB_TABLE`, and
unique ALT support uses the GATK `count <= threshold` boundary.  Multiallelic
support counts `TLOD` after log10-to-natural-log conversion against the hard
5.0 threshold; `PON` is presence-based and `NCount` ratios aggregate ALT AD
across all samples.  Site FILTER names are retained for existing callers and
the per-ALT status is recorded in stats/manifest telemetry.  Full
ErrorProbabilities joint posterior calibration remains an explicit fallback.
Number=A `TLOD`/posterior vectors are reduced by concrete ALT rather than only
ALT1: site TLOD uses the maximum, somatic posterior the minimum, and
germline/artifact posteriors the maximum; contamination and minimum-AF use all
tumour ALT fractions. Scalar legacy annotations broadcast across ALTs and
missing vector elements retain their ALT position, so mixed multi-ALT records
cannot shift evidence between sibling alleles.
When `NALOD` and `##normal_sample` metadata are present, it also runs the
Kokkos-backed NormalArtifactFilter subset: aggregate matched-normal/tumour AD,
the 10% normal-AF ratio gate, normal pileup binomial tail and NALOD posterior;
the posterior is compared with the effective FilterMutectCalls error threshold
before emitting `normal_artifact` and being counted in the manifest. Empirical
multi-pass threshold learning and the remaining joint filter model remain
explicit fallback. The pinned GATK 4.6.2.0 normal-artifact oracle covers the
positive-but-below-threshold posterior boundary on OpenMP and Serial.
`weak_evidence` now has a deterministic initial-model subset: per-ALT TLOD/AD
is evaluated with the flat/high-AF SomaticClusteringModel and a Kokkos log-beta
kernel.  `slippage` similarly consumes integer `RPA` plus `RU`, checks the
one-repeat/minimum-length gate, evaluates the regularized-beta artifact model
through Kokkos, emits `STRQ`, and records `slippage` telemetry.  The explicit
`--filter-error-probability-threshold` (default `0.1`), `--min-slippage-length`
and `--slippage-rate` options make this boundary auditable; GATK's four-pass
empirical model/threshold learning and cross-filter ErrorProbabilities remain
fallback.
When `--threshold-strategy` is explicitly selected, a streaming pre-pass also
reconstructs the exposed `PSOMATIC`/`PGERMLINE`/`PARTIFACT`/`OBP` error types
(falling back to TLOD+AD for sequencing error), adds deterministic
filter-specific probabilities available in the same pass, reduces correlated
ARTIFACT posteriors with per-ALT max, then combines independent types with the
GATK product rule, learns CONSTANT/FDR/F-score thresholds, and applies the
result per ALT as `error_probability`.  This is the VCF-observable joint
posterior subset; filter-specific Java ErrorProbabilities maps and the complete
empirical SomaticClusteringModel remain explicit fallback.
The same pass fits a deterministic SomaticClusteringModel with background/high-AF
beta-binomial clusters and GATK-shaped probability-weighted AF peak splitting
(up to five additional fuzzy-binomial clusters, BIC-gated, with five EM updates
per accepted split).  Its Kokkos likelihood is used for TLOD+AD fallback and STR
slippage; cluster count/weights/means are emitted in filtering-stats JSON.  The
background/high-AF updates use the Java-shaped digamma gradient (10 epochs,
rate 0.01 and the same lower bounds); fuzzy-binomial peaks use weighted AF
updates.  Filter-specific ErrorProbabilities maps and release-specific
calibration still differ, so this is not yet claimed bit-identical.
The native boundary now accepts GATK's `--log-snv-prior`, `--log-indel-prior`,
`--log-artifact-prior` (and `--log-somatic-prior` compatibility alias), plus
`--normal-p-value-threshold`, `--min-pcr-slippage-size`, and
`--pcr-slippage-rate` aliases.  A readable GATK `statistic/value` stats table
is parsed for `callable`; learned per-indel priors and callable-site provenance
are emitted in filtering-stats JSON without overwriting the input table.
For pipeline compatibility, an existing `--stats` file that is recognized as a
GATK Mutect2 stats table is treated as read-only input.  Native filtering stats
are written to `<output>.stats.json`, or to the explicit `--filtering-stats`
path, and the input path is recorded in the JSON/manifest.  The imported
threshold/model is preserved as provenance; applying GATK's full empirical
threshold learning remains an explicit fallback.
When `--stats` is a recognized GATK Mutect2 table, the native path also enables
the GATK default `OPTIMAL_F_SCORE` empirical pass when
`--threshold-strategy` is omitted; legacy VCFs without that table retain the
explicit adapter mode.
FilteringOutputStats-compatible per-filter FP/FDR/FN/FNR reductions are now
included under `filtering_output_stats` in native JSON.  The stats pass reduces
available filter probabilities with Java's max-within-`ErrorType` rule before
the independent-type product, while preserving original ALT indexing.  If
`--filtering-stats` ends in `.table` or `.tsv`, the requested path is emitted as
the GATK `filter FP FDR FN FNR` table with metadata, while a `.json` telemetry
sidecar is retained. Per-ALT probabilities are keyed by the original VCF ALT
index, so a symbolic `<NON_REF>` allele interleaved before a concrete ALT does
not shift the FilteringOutputStats attribution.  When Mutect2 exposes both
`PARTIFACT` and `OBP`, the native joint error reduction treats them as
correlated `ARTIFACT` filters and takes their per-ALT maximum; it then combines
that type independently with `NON_SOMATIC` and `SEQUENCING` using
`1-prod(1-p)`, matching Java `ErrorProbabilities`.
`benchmark_filter_mutect_calls.py` also reports a 95-record empirical-model
workload alongside the ordinary 128-record file-boundary throughput, including
the AD/POPAF contamination posterior batch (`contamination_posterior_records`
and `contamination_posterior_alleles`), learned cluster count and model
records/second.  This directly exercises the
quantile peak split, BIC gate and EM path rather than timing VCF I/O alone.
When standard FORMAT `PGT` and `PID` are present, FilterMutectCalls also runs
the GATK-shaped two-pass `FilteredHaplotypeFilter` subset: it learns the
maximum `PARTIFACT` posterior per `PGT+PID` key and applies it to records within
the inclusive `--max-intra-haplotype-distance` window (default 100bp), writing
the `haplotype` FILTER and learning counters to stats/manifest. The VCF
boundary now exports an observable per-filter stats subset, but the VCF does not
carry all Java ErrorProbabilities inputs; normal-artifact exclusion and the
complete joint posterior remain explicit fallbacks.
It now also consumes standard `LearnReadOrientationModel` tar.gz artifacts via
`--orientation-bias-artifact-priors`/`--ob-priors`: FASTA 3-mer context and
standard FORMAT `F1R2`/`F2R1` evidence are evaluated by a Kokkos 12-state
beta-binomial posterior kernel, with GATK-shaped `ROQ` and `orientation`
output. With priors present, all samples not declared by a ##normal_sample
header line are evaluated and aggregated with GATK's alternate-AD-depth weighted
median; --tumor-sample explicitly limits this prior filter to one sample. The
deterministic adapter defaults to a 0.5 posterior threshold when priors are
supplied; GATK's learned multi-pass threshold and full joint filter model remain
explicit fallback. The compatible threshold options
--threshold-strategy CONSTANT/FALSE_DISCOVERY_RATE/OPTIMAL_F_SCORE,
--initial-threshold, --false-discovery-rate and --f-score-beta perform a
streaming two-pass learned threshold over the native orientation posterior
subset; complete cross-filter joint learning remains fallback. Equal-length MNPs
use the GATK per-base maximum posterior rule; indels remain fail-closed in this
orientation-prior subset. The beta-binomial kernel keeps AD allele depth
separate from F1R2/F2R1 orientation trials and records any depth disagreement
in the OutputManifest telemetry. Standard Number=R F1R2/F2R1 vectors provide
REF+ALT trial depth exactly as GATK; Number=A/scalar legacy native records use
an explicit AD-REF denominator fallback, while AD never leaks into the
standard posterior. Verify this hand-off with
`python3 fastgatk-native/scripts/verify_filter_mutect_calls.py`, and measure it
with `python3 fastgatk-native/scripts/benchmark_filter_mutect_calls.py`.
All FilterMutectCalls scalar posterior bridges now use the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> RangePolicy -> execute ->
collect` lifecycle; the manifest and benchmark expose execution space, policy,
batch/observation counts and prepare/execute time. Multi-ALT records use a
single tumor-evidence `RangePolicy` launch, while single-ALT records keep the
lower-overhead scalar path; cross-record shard fusion remains a follow-up.
`-L/--intervals`/`--region` also reuse the shared HTSlib interval selector and
support `--interval-set-rule UNION|INTERSECTION`; the same subset is applied to
orientation, haplotype, empirical-learning and final filtering passes, with rule,
input statistics and skipped-record counts recorded in stats/manifest. OpenMP and
Serial contracts cover this boundary.
`CollectF1R2Counts` now closes the upstream orientation hand-off: HTSlib reads
BAM/CRAM, applies the Mutect2 default read-filter surface and CIGAR-aware
pileup projection, while a Kokkos atomic kernel aggregates per-locus base and
F1R2/F2R1 counts through the shared `HostBatch -> KernelPlan.prepare -> Kokkos
Views -> execute -> collect`/`RangePolicy` lifecycle. It emits the standard
three-member
`.ref_histogram/.alt_histogram/.alt_table` tar.gz contract, including sample
metadata and configurable `--f1r2-median-mq`, `--f1r2-min-bq` and
`--f1r2-max-depth`. The chr17 fixture matches GATK 4.6.2.0 histogram counts
and all alt-table rows; OA/chimeric edge cases, cloud input and report byte
identity remain explicit fallback scope. Run
`python3 fastgatk-native/scripts/verify_collect_f1r2_counts.py` or
`python3 fastgatk-native/scripts/benchmark_collect_f1r2_counts.py` for the
reproducible contract and file-boundary timing. When `-L` intervals and a BAM/
CRAM index are available, HTSlib uses indexed region iterators; missing indexes
fall back to a sequential scan and are reported in telemetry rather than
silently changing the result contract. Multi-input aggregation keys use the
sample name plus canonical contig name and position, so different BAM/CRAM
reference-header orders merge correctly; the verifier includes an indexed plus
sequential reordered-header regression.
`LearnReadOrientationModel` now has a native Kokkos path for the common Mutect2
hand-off: it accepts repeatable fastgatk F1R2/R1F2 TSV sidecars (plain or
`.tsv.gz`, decoded by zlib on Host), reduces
strand/ALT support with a Kokkos atomic kernel, and emits a real `.tar.gz`
containing a 64-context GATK-shaped `.orientation_priors` table. For the standard
`CollectF1R2Counts` tar it parses `.ref_histogram`, `.alt_histogram` and
`.alt_table`, runs the Java 12-state Beta-Binomial EM with canonical/reverse-
complement context merging and prior pseudocounts, and matches the bundled GATK
4.6.2.0 oracle within 4.5e-14 on the chr17 fixture (OpenMP/Serial). The 64-context
row order follows Java `HashMap` iteration and no-data rows use Java
`Double.toString` formatting, so those rows are byte-identical to GATK; active EM
rows retain the documented special-function tolerance. Standard-tar
EM observations are evaluated through the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle
using `RangePolicy` batches; the final weighted sum stays in deterministic Host
input order. The legacy TSV aggregate uses the same lifecycle and atomic
`RangePolicy`. OutputManifest records execution space, policy, batch/observation
counts, prepare/execute seconds, EM counters, sample metadata and archive
completeness.
Legacy TSV remains a compatibility approximation; tar member/gzip byte identity and
release-specific downstream filtering are not claimed. When one standard
CollectF1R2Counts archive contains multiple samples, native keeps each sample's
histogram/table namespace separate and writes one `<sample>.orientation_priors`
member per sample, matching GATK's multi-sample writer. The pinned two-sample oracle is
`verify_learn_read_orientation_model_multisample_gatk_oracle.py`.
That oracle also covers `--max-depth 2` against source histograms containing
higher-depth bins: native preserves those bins for Java-compatible example
counts and limits only the EM observation loop (maximum probability delta about
1.9e-12 on the pinned fixture).
For a reproducible file-boundary timing on synthetic F1R2 records:
`python3 fastgatk-native/scripts/benchmark_learn_read_orientation_model.py`.
`CollectReadCounts` now has a native HTSlib/Kokkos TSV/HDF5 path for CNV workflows
(the native default is HDF5, matching GATK; pass `--format TSV` for text interchange):
it consumes repeatable literal/Picard/BED intervals, applies the GATK default
mapped/non-duplicate/MAPQ filter, applies GATK-compatible
`--interval-merging-rule OVERLAPPING_ONLY|ALL`, and counts read starts with a Kokkos overlap kernel,
and emits the SAM-header-plus-`CONTIG/START/END/COUNT` SimpleCountCollection
format. Manifest/JSON telemetry records both requested intervals and merge events,
so overlap and adjacency behavior is reproducible. Each read batch uses the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle;
the manifest records `RangePolicy`, execution space, batch/record counts and
prepare/execute seconds, while fixed interval metadata views are reused across
batches. HDF5 output follows GATK
`HDF5SimpleCountCollection` paths for sample,
sequence dictionary, intervals and counts, then performs a native read-back
validation of the sample, interval dimensions and count values.  The
sequence-dictionary dataset is the HTSJDK `SAMTextHeaderCodec` subset
(`@HD\tVN:1.6` plus every normalized `@SQ` record, including AS/M5/UR/SP
tags), not a reduced name/length projection and not the complete BAM header.
The focused `verify_hdf5_simple_count_collection.py` gate writes both native
and GATK 4.6.2.0 HDF5 files, compares dictionary metadata, and uses Java's
`DenoiseReadCounts` reader to round-trip the interval/count row.  The manifest
records strict execution, output bytes, wall time and the native round-trip
result.  This closes only the SimpleCountCollection metadata slice; HDF5 SVD
PoN/model state and release-specific byte identity remain separate prototype
scope. Cloud indexed metadata remains
explicit fallback. The pinned GATK 4.6.2.0 read-start/filter oracle is in
`fastgatk-native/scripts/verify_collect_read_counts_gatk_oracle.py`; the broad
native contract is `fastgatk-native/scripts/verify_collect_read_counts.py`, and
the reproducible file-boundary benchmark is
`fastgatk-native/scripts/benchmark_collect_read_counts.py`.
`DenoiseReadCounts` now consumes that TSV/HDF5 boundary and requires both
`--denoised-copy-ratios` (or `-O/--output`) and `--standardized-copy-ratios`,
matching GATK's two-output writer contract; omitting either fails before input
processing. It implements the
fractional-coverage path with Kokkos map kernels: positive sample median,
epsilon-floored log2 and centered log2 median. An interval-matched TSV
`--panel-of-normals` is supported by taking sample/PoN fractional-coverage
ratios before the same kernel. A wide TSV panel (`CONTIG/START/END` followed by
sample columns) uses Kokkos GEMM-style Gram construction and deterministic
low-rank SVD projection; `--number-of-eigensamples` (or `--svd-rank`) controls
the removed components. The standardized and denoised files are identical when
no panel is supplied. GATK-style `GC_CONTENT` annotations from
`AnnotateIntervals` are consumed through the same 101-bin correction path;
when a panel is supplied, `--annotated-intervals` is ignored as in GATK (and
the manifest reports `gc_content=false` unless the HDF5 PoN contains its own
GC model); cloud metadata remains explicit fallback. The standardization, HDF5 PoN
standardization, Gram, projection and reconstruction stages all use the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle
with `RangePolicy`/`MDRangePolicy`; its manifest records strict determinism,
execution space, policy, kernel batch/observation counts, prepare/execute time,
complete output entries, output byte counts and wall time. The contract and
file-boundary benchmark are in
`fastgatk-native/scripts/verify_denoise_read_counts.py` and
`fastgatk-native/scripts/benchmark_denoise_read_counts.py`.
For HDF5 SimpleCountCollection input, the complete
`/locatable_metadata/sequence_dictionary` is now carried into both output
headers, including all `@SQ` AS/M5/UR/SP tags, matching GATK's
CopyRatioCollection writer. The pinned 4.6.2.0 oracle
`fastgatk-native/scripts/verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py`
compares exact headers and data rows on OpenMP and Serial; other
release-specific HDF5 metadata and Spark/cloud behavior remain explicit
fallback.
`AnnotateIntervals` is the native GC annotation producer. It reads FASTA/FAIDX and
GATK interval selectors, computes `GC_CONTENT` with a Kokkos batch kernel, and
optionally adds htsjdk-compatible length-weighted `MAPPABILITY` and
`SEGMENTAL_DUPLICATION_CONTENT` values from BED tracks (including score and NaN
handling). Track rows are indexed by contig and binary-searched in the Kokkos
`RangePolicy` kernel, so large WGS tracks do not trigger an interval-by-track full
scan. It writes an annotated-interval TSV directly consumable by the CNV tools.
All annotation stages follow `HostBatch -> KernelPlan.prepare -> Kokkos Views ->
execute -> collect`; the manifest records execution space, record/track counts,
index strategy and prepare/execute seconds for cross-backend comparisons. The
reproducible file-boundary benchmark is
`fastgatk-native/scripts/benchmark_annotate_intervals.py --with-tracks`；加入
`--include-java` 可在 pinned GATK 4.6.2.0 上做同输入输出 byte-identity 和 wall-time
对照。
`CountBasesInReference` provides the corresponding reference-walker boundary: Host/FAIDX
loads the FASTA and literal, Picard interval-list, and BED intervals (including `.gz`
compressed text; overlaps are merged), while a Kokkos `RangePolicy` atomic histogram counts raw reference bytes. The stdout/`-O`
report follows GATK's ascending-byte format, and the OutputManifest records the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle, execution
space, records and timing. The contract/oracle and benchmark are
`fastgatk-native/scripts/verify_count_bases_in_reference.py` and
`fastgatk-native/scripts/benchmark_count_bases_in_reference.py`; compressed interval files,
cloud reference and release-specific logging remain explicit fallback.
`CompareReferences` and `CheckReferenceCompatibility` cover the reference-tool boundary:
the first reproduces the MD5-keyed table, pair status analysis and `FIND_SNPS_ONLY` Kokkos
mismatch mask using adjacent `.dict` M5 values or recalculated MD5; the second reads
BAM/CRAM/SAM or VCF/BCF dictionaries and emits GATK-shaped compatibility rows. Both have
GATK 4.6.2.0 fixture oracles, CTest, dispatcher, OutputManifest and benchmarks. FULL_ALIGNMENT/
MUMmer, cloud inputs and full htsjdk/Tribble warning semantics remain explicit fallback.
`FastaReferenceMaker` and `FastaAlternateReferenceMaker` now provide the common FASTA-walker
boundary: dictionary-ordered merged literal/interval-list/BED output, GATK sequence IDs and
descriptions, uppercase sequence bytes, configurable line width and `.fai`/Picard `.dict` (M5)
generation; the alternate path applies simple
SNP/anchored indel records, SNP masks, mask priority and diploid IUPAC. Fixed-byte copies use
the shared Kokkos `RangePolicy` lifecycle, with GATK 4.6.2.0 oracles, CTest, dispatcher,
manifests and file-boundary benchmarks. Complex symbolic/overlapping indels, cloud inputs and
release-specific feature validation remain explicit fallback. For
`--use-iupac-sample`, homozygous genotypes now select the genotype allele exactly as GATK
(`getIUPACBase`), including ALT[1] for a multi-ALT 2/2 call; the pinned
`verify_fasta_alternate_iupac_hom_gatk_oracle.py` covers hom-ref/hom-ALT projection on
both Kokkos backends.
`ShiftFasta` extends the same boundary for circular references: each contig is rotated by the
requested (or default half-length) offset through a Kokkos `RangePolicy`, with GATK-compatible
FASTA wrapping, `.fai`/`.dict`, UCSC chain, regular/shifted interval sidecars and an OutputManifest.
This kernel is intentionally a byte-permutation/memory-bandwidth operation, so it uses the
portable Kokkos execution-policy API rather than `Kokkos::Experimental::simd`; arithmetic-heavy
paths continue to use the Kokkos SIMD API where vector lanes materially help.
Its GATK 4.6.2.0 oracle, round-trip contract, CTest, dispatcher entry and file-boundary
benchmark are `fastgatk-native/scripts/verify_shift_fasta.py` and
`fastgatk-native/scripts/benchmark_shift_fasta.py`; remote/cloud storage and release-specific
validation remain explicit fallback.
`IndexFeatureFile` now covers the HTSlib-indexable feature-file boundary: BGZF VCF/GVCF/BED
inputs produce standard `.tbi` indexes and BCF inputs produce `.csi`, with configurable HTSlib
decoder threads and a Kokkos `RangePolicy` metadata lifecycle recorded in the manifest. Plain
VCF paths use HTSJDK's `DynamicIndexCreator(FOR_SEEK_TIME)` choice: sparse ordinary VCFs use
the adaptive 2000 bp LinearIndex candidate, dense VCFs use the 75-record interval-tree
candidate, and empty VCFs retain the empty interval-tree form. `.g.vcf` paths additionally
use GATK's 128000 bp starting bins and adaptive occupied-contig bin merging, and the pinned
GVCF corpus is byte-identical to GATK 4.6.2.0. Native ordinary indexes are HTSJDK query
compatible; the pinned sparse one-record VCF is also byte-identical after reproducing
HTSJDK's leading-empty-block optimization, and the dense 300-record fixture is byte-identical
after reproducing the HTSJDK red-black-tree insertion/rotation preorder. Uncompressed BED
also uses the same DynamicIndexCreator/LinearIndex v3 path and matches the pinned GATK BED
index. Other uncompressed Tribble dynamic indices
and unsupported codecs remain explicit fallback. Contract and file-boundary benchmark are
`fastgatk-native/scripts/verify_index_feature_file.py` and
`fastgatk-native/scripts/benchmark_index_feature_file.py` (pass
`--input path/to/sample.g.vcf` to benchmark the linear `.idx` path).
`CountReads` and `FlagStat` now cover the common read-QC boundary with bounded HTSlib
SAM/BAM/CRAM batches. Literal, Picard interval-list, and BED selectors may be plain or
`.interval_list.gz/.bed.gz`; HTSlib decodes compressed text before the shared filter path.
CountReads uses a Kokkos `RangePolicy` reduction; FlagStat uses
Kokkos atomic counters for all 12 GATK fields and treats missing read/mate coordinates as
unmapped like HTSJDK. Both accept repeatable `-I/--input` shard files and aggregate records/
counters in input order, while OutputManifest keeps the single-input `input` field and adds
`inputs`/`input_count`. Both accept `-R/--reference` for CRAM reference injection and repeatable
`-L/--intervals`/`-XL` selectors (indexed HTSlib iteration with sequential filtering fallback),
including GATK `-isr/--interval-set-rule UNION|INTERSECTION` across repeatable selectors
(an interval-list/BED file is treated as one selector set),
`-ip/--interval-padding` and `-ixp/--interval-exclusion-padding` with contig-boundary
clamping, and `-imr/--interval-merging-rule ALL|OVERLAPPING_ONLY` (multi-region indexed
traversal de-duplicates alignments just as GATK does),
including exclusion-only traversals where a read is retained when its alignment span still
touches non-excluded territory,
and the 30 native `-RF/--read-filter`, `-XRF/--inverted-read-filter` and `-DF/--disable-read-filter` predicates, including
parameterized `MappingQualityReadFilter` with GATK-compatible inclusive
`--minimum-mapping-quality`/`--maximum-mapping-quality` bounds, plus exact
fragment-length, mate-contig, CIGAR-N/GoodCigar, read-group and mate-strand predicates.
`ReadNameReadFilter` is also supported with repeatable `--read-name` values (and enforces
the same filter/argument dependency as GATK); `ReadGroupReadFilter` similarly accepts
repeatable `--keep-read-group` RG identifiers.
`ReadGroupBlackListReadFilter` accepts repeatable `--read-group-black-list ATTR:VALUE`
expressions and resolves RG-header attributes through HTSlib.
Numeric `ReadTagValueFilter` comparisons are available through
`--read-filter-tag`, `--read-filter-tag-comp`, and `--read-filter-tag-op` with the six
GATK comparison operators.
`-XRF` evaluates the same bounded native predicate under GATK `InvertedReadFilter`
semantics after default and ordinary `-RF` predicates; the pinned parameterized MAPQ
FlagStat boundary is exact on OpenMP and Serial.

`GetPileupSummaries` likewise accepts `-XRF/--inverted-read-filter` for its 14 explicit
LocusWalker predicates.  It applies the inverse after the ordinary default/`-RF` mask;
the pinned low-MAPQ table-row oracle is byte-exact against GATK 4.6.2.0 on both backends.
They provide
optional text/OutputManifest output, GATK 4.6.2.0 chr17 full/interval/filter oracles,
OpenMP/Serial contracts, dispatcher entries and `fastgatk-native/scripts/benchmark_read_metrics.py`.
As in GATK's `ReadWalker`, `WellformedReadFilter` is enabled by default (header/RG,
coordinate, sequence/CIGAR-length, GoodCigar, and N-op checks); `--disable-tool-default-read-filters`
removes only this implicit filter, while explicit predicates remain active. The malformed
SAM regression is included in both backend contracts.
Passing `--include-java` to that benchmark adds the pinned GATK 4.6.2.0 Java
baseline; tiny fixtures are startup-dominated and should not be read as a WGS
throughput claim.
The dispatcher keeps this larger read-QC surface scoped to CountReads/FlagStat via
their registry `native_read_filters` lists; HC/Mutect2 apply the GATK standard
GoodCigar/non-zero-reference mask at the CLI boundary, while unknown classes still
fail closed to explicit Java fallback.
`SplitIntervals` provides deterministic FASTA-dictionary scatter shards for SLURM/Nextflow
(`INTERVAL_SUBDIVISION`, `INTERVAL_COUNT`, distributed remainder/overflow modes,
exclusions, `--dont-mix-contigs`, and reference-derived `--min-contig-size`). Literal,
Picard interval-list, and BED selectors may be plain or `.interval_list.gz/.bed.gz`; the
compressed text path is decoded by HTSlib on Host before entering the common kernel. Its base
summary uses the shared `HostBatch -> KernelPlan.prepare -> Kokkos View -> RangePolicy -> collect`
lifecycle and records execution-space/policy/timing telemetry in the manifest; the
contract/oracle and benchmark are `verify_split_intervals.py` and
`benchmark_split_intervals.py`.
`FilterIntervals` now consumes that annotated TSV and applies the GATK inclusive
GC-content, mappability, and segmental-duplication bounds with a Kokkos mask
kernel. Interval-list/BED/literal selectors, exclusion intervals, deterministic
contig ordering, GATK's `--interval-set-rule UNION`, solitary-contig protection,
Picard interval-list output, and an OutputManifest are covered by the contract. TSV count collections now support
GATK low-count and Apache Commons Math percentile-based extreme-count filters;
native GATK SimpleCountCollection HDF5 inputs are also supported through the
shared HDF5 reader; the mask stage follows `HostBatch -> KernelPlan.prepare ->
Kokkos Views -> execute -> collect` with `RangePolicy`, and records execution
space, record count and prepare/execute seconds; cloud/object-store inputs remain explicit fallback. The GATK short aliases
`-imr`, `-isr`, `-ip`, and `-ixp` are accepted; default/explicit non-`OVERLAPPING_ONLY`
merging, non-zero `-ip`/`-ixp`, and `-isr INTERSECTION` are rejected exactly as GATK
4.6.2.0. Its reproducible benchmark is
`fastgatk-native/scripts/benchmark_filter_intervals.py`.
`PreprocessIntervals` now provides the CNV interval-preparation boundary. It reads
FASTA/FAIDX and optional interval-list/BED/literal selectors (omitted `-L` means
whole-reference), merges overlapping intervals, subtracts repeatable
`-XL/--exclude-intervals` selectors before padding, applies GATK-compatible
padding, generates fixed bins (or preserves padded intervals when `--bin-length 0`), and
uses a bounded Kokkos batch kernel to remove bins containing only N bases. The
filter stage follows `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute ->
collect` with `RangePolicy`; manifest/benchmark telemetry records execution space,
batch/record counts and prepare/execute seconds. `PreprocessIntervals` keeps GATK's
validation that `--interval-exclusion-padding` must be zero. The result is a Picard
interval-list directly consumable by `CollectReadCounts`; the
contract is `fastgatk-native/scripts/verify_preprocess_intervals.py`.
The file-boundary benchmark is `fastgatk-native/scripts/benchmark_preprocess_intervals.py`.
The native `CreateReadCountPanelOfNormals` is the upstream PoN writer. It accepts
repeatable TSV or SimpleCountCollection HDF5 inputs, applies the GATK-ordered
interval/sample zero filters, fractional medians, imputation and outlier
truncation, and uses Kokkos Gram/power-iteration SVD to write the
`HDF5SVDReadCountPanelOfNormals` v7 paths, with bounded row-chunked matrix
persistence controlled by `--maximum-chunk-size`. Gram construction and left-vector
projection use the shared `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`
(`MDRangePolicy`/`RangePolicy`) lifecycle, with kernel telemetry in the manifest. The resulting HDF5 file is read by
both native DenoiseReadCounts and GATK 4.6.2.0; GC annotated intervals are
consumed through the GATK 101-bin GC-bias correction when annotations are
supplied. Cloud and Spark/report extensions remain explicit fallback. See
`fastgatk-native/scripts/verify_create_read_count_panel_of_normals.py`.
The reproducible file-boundary benchmark is
`fastgatk-native/scripts/benchmark_create_read_count_panel_of_normals.py`.
`CallCopyRatioSegments` consumes the segment table emitted by a segmentation
backend and applies the GATK SimpleCopyRatioCaller two-pass, length-weighted
copy-ratio statistics. The `2^log2` transform runs through a portable
`Kokkos::Experimental::simd<double>` batch kernel (SIMD width/groups and kernel
time are recorded in the manifest), with scalar fallback provided by Kokkos on
non-SIMD backends. It emits the called segment table and optional IGV legacy
`.seg` file. The SIMD batch uses the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle
with `RangePolicy`; manifest telemetry includes policy, batch/observation and
prepare/execute times. The deterministic TSV ModelSegments boundary is implemented below,
while the probabilistic model and HDF5/PoN semantics stay explicit fallback. See
`fastgatk-native/scripts/verify_call_copy_ratio_segments.py` and
`fastgatk-native/scripts/benchmark_call_copy_ratio_segments.py`.
The native CLI accepts GATK's full 4.6.2.0 z-score spellings
`--outlier-neutral-segment-copy-ratio-z-score-threshold` and
`--calling-copy-ratio-z-score-threshold`; the shorter legacy spellings remain
accepted for existing workflows. The pinned Java/native oracle covers both.
`CollectAllelicCounts` now provides the complementary BAM/CRAM/SAM-to-TSV CNV
path: HTSlib decodes reads and reference, Kokkos atomically aggregates ACGT
observations, empty loci are emitted, reference-N loci are omitted, and ALT is
defined as total ACGT minus REF. If `--sample` is omitted, a unique `@RG SM`
is discovered from the input header; zero or multiple samples fail closed.
The count stage follows `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute
-> collect` with `RangePolicy`; the manifest records execution space,
batch/observation counts, prepare/execute seconds, output bytes and wall time.
HDF5/cloud/report extensions remain explicit fallback. See
`fastgatk-native/scripts/verify_collect_allelic_counts.py` and
`fastgatk-native/scripts/benchmark_collect_allelic_counts.py`.
`DepthOfCoverage` now adds the coverage-walker boundary on the same architecture:
HTSlib/Host traverses BAM/CRAM/SAM and applies GATK's default Wellformed,
mapped, non-duplicate, non-secondary-read filters; QC-fail and supplementary
reads remain countable, while a Wellformed CIGAR reference-skip (`N`) rejects
the complete read. A Kokkos atomic kernel counts per-locus
depth/base observations and emits the GATK locus CSV plus interval, sample,
histogram and cumulative sidecars. The default single-SM `COUNT_READS` path is
byte-identical to GATK 4.6.2.0 on the pinned chr17 fixture (all seven output
files). The manifest records strict determinism, Kokkos execution space, read
filter/observation counts, output bytes and wall time. The count stage now uses
the shared `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`
lifecycle and records batch/observation/prepare/execute telemetry (with
`RangePolicy` and the selected Kokkos execution space). The shared Host header
contract resolves `@RG ID` to `SM`, so multi-sample read-group partitions and
explicit `--sample` selection are deterministic across all summary/cumulative
outputs. Reference-N loci follow GATK semantics: they are excluded by default
and included with `--include-ref-n-sites`; this is covered by a synthetic
FAIDX/Java oracle in `verify_depth_of_coverage.py`. See
`verify_depth_of_coverage_multisample.py`. The native `COUNT_FRAGMENTS` mode
is an explicit read-name+locus de-duplication extension: pinned GATK 4.6.2.0
currently rejects fragment-based counting, so this mode is not advertised as
bit-identical. Genes, cloud feature traversal and unbounded whole-reference
materialization remain explicit fallback; `--omit-interval-statistics` and
`--omit-per-sample-statistics` suppress the same summary pairs as GATK (legacy
`--omit-intervals`/`--omit-sample-summary` aliases are accepted), and
`-XL/--exclude-intervals` is applied after included spans, with `-ip/-ixp`
padding and `--interval-set-rule UNION|INTERSECTION` handled by the shared
HtsReader; the post-subtraction locus set drives every summary denominator.
`--max-loci` is a memory guard. See
`fastgatk-native/scripts/verify_depth_of_coverage.py` and
`fastgatk-native/scripts/benchmark_depth_of_coverage.py` (the benchmark reports
both the pinned single-sample path, a two-sample `@RG→SM` workload, and the
native fragment-counting extension).
The advanced `--ignore-deletion-sites` boolean is also implemented: when
combined with `--include-deletions`, native suppresses deletion pileups from
effective depth and base-count output, matching GATK 4.6.2.0. The pinned
`verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py` compares
deletion-only loci and manifest effective settings on both CPU backends.
`ModelSegments` now closes the TSV CNV chain with a deterministic native
prototype: it validates denoised copy-ratio rows and optional allelic-count
loci, segments contiguous intervals with a bounded change-point threshold, and
emits `modelFinal.segments.tsv` plus OutputManifest telemetry. Repeating
`--denoised-copy-ratios`/`--allelic-counts` selects a bounded common-interval
joint segmentation and emits a Picard `.interval_list` for downstream runs.
`--mode`,
`--num-samples` and `--num-burn-in-iterations` select a bounded,
counter-driven deterministic random-walk Metropolis sampler in Kokkos
`RangePolicy`; before the chain, a deterministic Gibbs-conditional
responsibility pass uses midpoint-assigned raw points in Kokkos Views to update
point outlier weights, segment means, global variance and outlier probability.
The conditional means and global parameters initialize/scale the bounded chain.
Kokkos retains only a fixed sample block and Host maintains a fixed-capacity
reservoir per segment before aggregating posterior quantiles, so peak memory does
not grow linearly with `--num-samples`. Manifest telemetry identifies the
conditional model, global parameters, sampler, controls, block and reservoir
bounds, and repeatability is contract-tested. This is not a claim of Java
release-specific slice-sampler draws, latent-indicator MCMC, or report raw-bit
identity; complete
multi-sample fitting and HDF5/PoN/GC normalization remain explicit fallback.
When `--output-prefix` is supplied, the native path also materializes
`.modelBegin.seg`, `.modelFinal.seg`, `.cr.seg`, `.cr.igv.seg` and `.af.igv.seg`
sidecars and records their paths in the manifest for Picard/IGV/Nextflow hand-off.
GATK's advanced `--segments` Picard interval-list is honored as a fixed modeling
partition (midpoint point ownership, dictionary-order/non-overlap validation),
and skips kernel segmentation. The native writer additionally emits
`.modelBegin/.modelFinal.cr.param/.af.param` GATK-shaped posterior-decile tables;
these are deterministic bounded summaries rather than a claim of Java MCMC
bit identity. Standard GATK MCMC/smoothing option names are accepted and
recorded; Java's multi-round credible-interval smoothing is not executed by
the bounded native sampler and is explicitly reported as
`smoothing_applied=false`.
The manifest also marks strict determinism and reports modeled-output,
sidecar-output and wall-time telemetry for both single-sample and bounded
multi-sample interval-list paths. Heterozygous filtering, boundary detection,
kernel embedding and bounded posterior blocks use the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` lifecycle
and expose execution-space, policy, batch/observation and prepare/execute
telemetry, making resource use auditable under SLURM.
See `fastgatk-native/scripts/verify_model_segments.py`.
The reproducible file-boundary benchmark is
`fastgatk-native/scripts/benchmark_model_segments.py`.

Before segmentation, allelic counts now use a GATK-shaped
Beta-integrated heterozygous-site filter.  The native tool accepts
`--minimum-total-allele-count-case`, `--minimum-total-allele-count-normal`,
`--genotyping-homozygous-log-ratio-threshold`, and
`--genotyping-base-error-rate`; matched-normal mode intersects case loci with
normal heterozygous loci, while case-only multisample mode intersects all case
sets.  Position matching uses interval overlap, and single-sample runs emit
`.hets.tsv` (plus `.hets.normal.tsv` when a matched normal is supplied).  The
manifest records the retained locus count and all genotyping controls.  Full
Java MCMC/model-parameter/report parity remains explicit fallback.
`GatherTranches` adds the VQSR scatter/gather report boundary: it validates
plain or `.tranches.gz` version-6 inputs (zlib-decoded on Host),
version-6 VQSLOD CSV rows, merges scattered rows by `minVQSLod`, combines Ti/Tv
and truth-site counts with GATK's transition/transversion decomposition, and
implements GATK's stateful target-truth-sensitivity walk (including the
`callsAtTruthSites` stable output ordering and version-5 filter-name formatting).
It requires the same explicit `--mode` argument as GATK, accepts both
`--truth-sensitivity-tranche` and `-tranche`, and treats a header-only
version-6 scatter slice as a zero-row contributor.  The contract is compared
to a real GATK 4.6.2.0 oracle byte-for-byte for both a multi-shard fixture and
an empty-slice merge; VariantRecalibrator model fitting and full
recalibration-table/allele-specific ApplyVQSR annotation remain explicit
fallback. A `.gz` output suffix enables zlib streaming compression and is
recorded in the OutputManifest without changing the decompressed tranche
text. See `fastgatk-native/scripts/verify_gather_tranches.py` and
`fastgatk-native/scripts/benchmark_gather_tranches.py`.
`ApplyVQSR` provides the matching native VCF boundary when the recalibration
artifact is a VCF carrying INFO/VQSLOD: it joins records by contig/position/
alleles, applies a tranche or explicit LOD cutoff, preserves mode-specific
records, and emits BGZF/Tabix output. `--create-output-variant-index` accepts
GATK optional-boolean syntax (default true); false/0 leaves BGZF output without
`.tbi` and records the choice in OutputManifest. Recalibration-table and allele-specific
model paths remain explicit fallback, while an AS recal VCF carrying
`AS_VQSLOD` is consumed allele-by-allele and emits `AS_FilterStatus`; missing
AS fields fail closed. Tranche CSV inputs may be plain or `.tranches.gz` (decoded
by HTSlib on Host). See
`fastgatk-native/scripts/verify_apply_vqsr.py`; scalar output uses GATK's
lower-case `INFO/culprit` (legacy upper-case `CULPRIT` is accepted on input),
and the pinned Java oracle compares all three scalar record semantics plus
two single-ALT AS score/status records. The Java spelling
`--use-allele-specific-annotations true/false` (including `=true/false`) is
accepted alongside the native `--AS` alias. See also
`fastgatk-native/scripts/benchmark_apply_vqsr.py`.
The AS oracle also fixes the GATK recalibration-file shape: Java emits one
recalibration record per ALT with scalar `culprit`; native joins those records
by contig/POS/REF/ALT and projects the values into Number=A
`AS_VQSLOD`/`AS_FilterStatus`/`AS_culprit` vectors on the input multi-ALT
record. `verify_apply_vqsr_gatk_oracle.py` and the
`fastgatk-apply-vqsr-gatk-oracle` CTest pass on both OpenMP and Serial. This is
a bounded AS VCF boundary; complete VQSR model, recalibration-table, and
cross-run provenance semantics remain explicit fallback.
Repeated `-L/--intervals`/`--region` selectors reuse the shared HTSlib interval
parser and support `--interval-set-rule UNION|INTERSECTION`; the subset is applied
before mode/scoring, empty intersections fail closed, and rule/input/skip telemetry
is preserved in the summary and OutputManifest. OpenMP and Serial contracts cover
this boundary.
`VariantRecalibrator` now provides the preceding VQSR scoring boundary: it
accepts GATK-labelled training/truth/known resource VCFs, extracts numeric INFO
annotations, computes training-only normalization, applies each resource's
Phred `prior` as a log10-odds shift at the contrastive scoring boundary, and
retains records with missing annotation values using GATK-compatible seeded
normalization and 20-draw Gaussian marginalization (the complete rows remain
on the Kokkos missing-dimension path), and scores a deterministic diagonal Gaussian
model through Kokkos,
and optionally fits a diagonal GMM with GATK's default
`--max-gaussians 8` (overrideable down to one) and bounded EM iterations
(`--max-iterations` default 150, K-means default 100). GATK's separate
`--max-attempts` model-build retry budget defaults to 1; native preserves both
controls and records them independently in OutputManifest. The
`--maximum-training-variants` (alias `--max-num-training-data`, default
2,500,000) cap applies a deterministic Host-side bounded sampler before Kokkos
training to bound resident memory; the sampler consumes the same JavaRandom
LCG/Fisher–Yates order as GATK's `Collections.shuffle`. Negative-model rows are
selected after a positive-model
pass using `--bad-lod-score-cutoff` (alias `--bad-lod-cutoff`, default -5); if a
small fixture has no rows under the cutoff, the deterministic worst
`--minimum-bad-variants` rows (default 1,000) are used as a bounded fallback and
the counts are recorded in the manifest. The model can also be a
regularized full-covariance GMM with `--full-covariance`. E-step,
M-step and scoring batches use the shared `HostBatch -> KernelPlan.prepare ->
Kokkos Views -> execute -> collect` lifecycle (`RangePolicy`/`TeamPolicy`), with
aggregate policy/batch/observation/prepare/execute telemetry in the manifest. It
applies GATK's Normal-Wishart mean-prior covariance as one
`shrinkageFactor * deltaMu * deltaMu^T` term. The pinned two-Gaussian oracle
matches Java positive/negative PMix, means and covariances within `5.23e-9`
and six VQSLOD values within `5e-5`; the former squared-factor path differed by
`1.22744`. It
emits a VQSLOD recal VCF plus tranches and manifest. Compressed recal VCFs receive
.tbi; plain recal VCFs receive a shared HTSJDK/Tribble LinearIndex v3 .idx,
recorded as vcf_index in the manifest and verified by a GATK interval query.
GATKReport v1.1 models may
be read from or written to plain text or `.model.gz`; compressed models are
decoded/encoded on Host and reuse the same Kokkos scoring path. The manifest records the
effective prior log10-odds range. With `--AS`, annotations
declared as `Number=A` are scored per ALT and emitted as `AS_VQSLOD`; scalar
annotations fail closed. SNP mode includes MNP records and INDEL mode includes
mixed/symbolic records using GATK's `checkVariationClass` boundary. The
`--AS` and `--use-allele-specific-annotations` spellings both accept GATK's
optional boolean form (`flag`, `true`/`false`, or `1`/`0`); an explicit false
selects scalar scoring and invalid literals fail closed. Recal VCFs now also
carry the resource-derived
`POSITIVE_TRAIN_SITE`/`NEGATIVE_TRAIN_SITE` flags (stale input flags are
removed before regeneration); filtered/non-variant resource records are
excluded like GATK's `isValidVariant`, so the native output can be passed
directly to the native or GATK ApplyVQSR boundary. Missing-value draws use the
GATK random seed and deterministic model-pass ordering, so repeated runs are
byte-stable. Decode also supports GATK's zero/endpoint annotation jitter and
`--output-tranches-for-scatter`: this switches the tranche file to GATK's
version-6 `requestedVQSLOD` schema for scatter/gather, with repeatable
`--vqslod-tranche` overrides (the default -10..10 resolution is generated when
no overrides are supplied). Empty slices above the observed maximum retain the
GATK empty-tranche counts, and the manifest records `scatter_tranches` and the
number of requested slices. This makes the resulting tranche shards directly
consumable by `GatherTranches`; Java raw-bit VBEM/resource calibration and
full recalibration-table/provenance remain explicit fallback.
`--mq-cap`/`--mq-jitter` controls, with seed and effective values in the manifest;
Java raw-bit VBEM convergence, recalibration-table and exact AS model semantics
remain explicit fallback. See
`fastgatk-native/scripts/verify_variant_recalibrator.py` and
`fastgatk-native/scripts/benchmark_variant_recalibrator.py` (use
`--max-gaussians 2 --full-covariance` to benchmark the full-covariance path).
The bounded model-artifact boundary is also covered against pinned GATK 4.6.2.0:
the Host reads GATKReport v1.1 `AnnotationMeans`/`AnnotationStdevs`, positive and
negative Gaussian PMix, means, and covariance tables, reorders them by annotation
name, and reuses the Kokkos scoring path. The pinned Java oracle compares VQSLOD,
coordinates/record count, and stable positive-resource labels plus model-table
dimensions; GATK's model-side negative-site selection is not serialized and is
therefore reported as a non-bit-identical metadata boundary. Scalar
recalibration VCF output follows GATK's `N/<VQSR>` dummy allele, `END`, lower-case
`culprit`, missing QUAL/FILTER, and `--sites-only-vcf-output` eight-column shape.
The OpenMP and Serial model-oracle CTests pass. Java raw-bit VBEM convergence,
complete recalibration-table/provenance, and full AS model parity remain explicit
fallback; this is not a claim of complete VQSR equivalence.
VariantRecalibrator also follows GATK's input-validity boundary for model
normalization: a training annotation with standard deviation below `1e-5`
fails closed instead of being silently assigned unit scale. The pinned
constant QD/MQ Java/native oracle (`fastgatk-variant-recalibrator-zero-variance-gatk-oracle`)
passes on OpenMP and Serial; full VQSR raw-bit VBEM parity remains fallback.
`GenotypeGVCFs` now has a native HTSlib contract-compatible local path that materializes all concrete
non-`<NON_REF>` alleles from GVCFs, remaps multi-allelic GT/AD/PL with the VCF
triangular genotype-index formula, and writes BGZF/Tabix output. It accepts the
GATK `--heterozygosity`/`--indel-heterozygosity` parameters and implements the
`GenotypePriorCalculator.assumingHW` SNP/INDEL/OTHER formula. The optional
`--use-posteriors-to-calculate-qual` (`--gp-qual`) path evaluates those priors in
a reusable Kokkos posterior kernel and records the model/telemetry in
`OutputManifest`; default GT/GQ/PL assignment remains PL-based, as in GATK.
The cohort calculator switches `--use-new-qual-calculator`/`--new-qual` and
the posterior-QUAL aliases all use the shared GATK optional-boolean parser:
bare/`true`/`1` enables, `false`/`0` disables, and empty or unknown values fail
before the input is opened.
The GATK-compatible `--genotype-assignment-method` option additionally accepts
`PREFER_PLS` (default), `USE_PLS_TO_ASSIGN`, `SET_TO_NO_CALL`,
`SET_TO_NO_CALL_NO_ANNOTATIONS`, `BEST_MATCH_TO_ORIGINAL`,
`DO_NOT_ASSIGN_GENOTYPES`, and an explicit
`USE_POSTERIOR_PROBABILITIES` mode.  The latter consumes the same
`GenotypePriorCalculator.assumingHW` per-site prior as GATK and selects GT/GQ through a dedicated Kokkos posterior-assignment
kernel; its call count and execution-space timings are recorded separately.
The native output also writes GATK FORMAT/GP and FORMAT/PG phred-scaled
posterior/prior vectors from the same Kokkos result.
`USE_POSTERIORS_ANNOTATION` now preserves FORMAT/PP, remaps its Number=G
rows through the shared Kokkos path, and assigns GT/GQ from the existing
phred-scaled posterior vector; its assignment timing is recorded separately.
`SET_TO_NO_CALL_NO_ANNOTATIONS` now enumerates the output FORMAT dictionary
and removes every present genotype-level annotation (including GP/PG/PP,
FT/MIN_DP, and caller-specific fields), retaining only a no-call GT; the
removed-field count is recorded in OutputManifest telemetry.
The default path now also runs a Kokkos port of GATK's
`AlleleFrequencyCalculator` EM/Dirichlet model, writes MLEAC/MLEAF, and uses
the cohort posterior for QUAL; the oracle reports exact GT/GQ/PL and numerical
QUAL/AF agreement on the bundled GATK chr17 locus. For records carrying a
spanning-deletion `*` ALT, the shared cohort kernel matches GATK's non-variant
set by including REF-only and REF/`*` genotypes in P(no-variant); a concrete
deletion span retains `*` as an independently countable ALT, while an orphan
`*` is removed at the final compatibility output subset. The opt-in
`--use-posteriors-to-calculate-qual` path uses the same non-variant set, so
posterior QUAL does not regress to REF/REF-only semantics at spanning-deletion
sites. Its Kokkos call count and backend timing are recorded in the manifest,
and the contract exercises the default and posterior-QUAL paths on the same
spanning-deletion fixture. Indexed-reference mode now
splits a reference block that spans another sample's concrete variant and
 merges its point segment into the joint locus. A shared-prior Kokkos
cross-sample reference-confidence kernel now emits joint `RCQ`/`RCP` INFO
annotations and records prepare/execute telemetry; release-specific GATK
calibration, opaque TileDB GenomicsDB semantics, full-genome provenance/raw bit
identity and cloud workspace behavior still require explicit GATK fallback.
`--include-non-variant-sites` is supported for
reference-backed inputs: a pure `<NON_REF>` block is expanded to one REF-only
site per coordinate, using the FASTA base and `MIN_DP`/RGQ-compatible depth;
zero-quality sites retain GATK's `./.` no-call. The default mode continues to
drop pure reference-only groups; an indexed `-R` is required when a block must
be split around an internal cross-sample candidate.
The explicit `--gatk-compatible-annotations` (alias
`--strict-gatk-annotations`) profile removes native-only `RAW_MQandDP`, `RCQ`,
`RCP` and `SB` fields after deriving the standard GATK `MQ`, `QD`, `FS` and
`SOR` annotations. The bundled GATK oracle now covers the complete
`17:69000-70000` window (1,001 total sites: three variants plus 998 REF-only): it
asserts exact INFO/FORMAT key-set parity, semantic header-definition parity,
and raw data-row text identity under this profile. The writer applies GATK's
INFO order, htsjdk float formatting (including rank-sum scientific notation),
GVCF-band header cleanup, and uncovered-site `./.` semantics. Remaining full
file byte differences for GVCF remain in RCM block semantics and provenance;
ordinary VCF is byte-identical when command-line provenance is disabled.
The cohort oracle also runs a real two-sample
`HaplotypeCaller → RenameSampleInVcf → CombineGVCFs → GenotypeGVCFs` chain and
requires exact sample order and data-row text. A separate strict fixture carries
two concrete ALTs plus `<NON_REF>` with the complete 10-entry diploid PL vector;
it verifies concrete-ALT union/remapping, EM `MLEAC/MLEAF`, and GATK's
multi-allelic `ExcessHet` classification. A companion 20-sample oracle now
computes the full GATK `GenotypeUtils` multiallelic best-ALT projection and
requires `InbreedingCoeff=0.0022` plus exact row/header text. The strict writer
also canonicalizes `GT:AD:DP:GQ:PL` FORMAT/sample order when source headers use
a different declaration order. QualByDepth's high-QD branch uses the
release-seeded Java `Random.nextGaussian()` sequence, so this parity remains
deterministic across runs. The default path uses the shared Kokkos
`AlleleFrequencyCalculator` EM/Dirichlet kernel for cohort and HC candidate-site
MLEAC/MLEAF; full HC multi-allelic assembly provenance and release-complete
AFCalculationResult details remain explicit expansion work. The bounded
InbreedingCoeff reducer is kept at the Host annotation boundary; joint PL,
AF/EM and other high-volume likelihood/count work remain in the shared Kokkos
kernels.
`verify_hc_complex_multiallelic_oracle.py` additionally runs GATK's bundled
tetra-ploid/tetra-allelic assembly fixture and requires native locus/ALT
provenance, the complete 15-cell Number=G PL row, executed PairHMM, one-read-
one-allele FORMAT/AD ownership (`2,15,17:34`), and complete row text identity
including QUAL, INFO annotations and rank-sum values. OpenMP and Serial both
match the pinned GATK row byte-for-byte for this corpus.
For reblocked inputs that omit FORMAT/PL, the cohort path also matches GATK's
`GenotypeUtils.genotypeIsUsableForAFCalculation`: diploid hom-ref calls with
FORMAT/GQ use the deterministic approximate PL model (0 for hom-ref, GQ for
ref-containing heterozygotes, 10×GQ for hom-var rows). The approximation count
and Kokkos execution telemetry are recorded in `OutputManifest`; non-hom-ref
or unsupported-ploidy PL-less calls remain excluded rather than being guessed.
Repeated `-L`/`--intervals` selectors use GATK `--interval-set-rule UNION`
(default) or `INTERSECTION`; literal, Picard
interval-list and BED selectors share the Host interval reader and their parse
counts are recorded in the manifest.

For whole-cohort jobs, `--stream-by-locus` enables the bounded joint path. It
probes only concrete variant spans, keeps one HTSlib cursor/current record per
input, performs a stable Host k-way merge, and sends each merged locus through
the same Kokkos-backed `ThreeStagePipeline`. Reference-block segments are
generated lazily around cross-sample variants, so a long block does not create
an in-memory segment vector. The span index is proportional to concrete
variants and is reported as `stream_span_probe_records`; `streamed_loci`,
`stream_max_inflight_records`, and `streamed_peak_host_bytes` are included in
the manifest. With `-L` and a `.tbi`/`.csi` sidecar, each cursor now uses
HTSlib's indexed iterator; inputs without an index retain an explicit
sequential fallback. Dense `--include-non-variant-sites` expansion is also
lazy in stream mode: one-base REF-only records are emitted through the same
locus heap, without retaining the block's full base vector.

The linked de-Bruijn path has a real-corpus gate in
`scripts/verify_real_assembly_graph.py`. It runs the bundled
`NA12878.chr17_69k_70k.dictFix.bam` over `17:69000-70000` with adaptive
pruning, duplicate-reference-kmer handling and PairHMM enabled, and checks the
graph node/edge population, reference connectivity, haplotype materialization
and artificial-haplotype recovery. The same run with
`--disable-artificial-haplotype-recovery` must report zero recovery paths and
bases. This is an internal assembly/telemetry oracle, not a claim of Java GATK
VCF byte identity.

The GenomicsDBImport adapter also writes a deterministic
`fastgatk-inputs.tsv` workspace index. In explicit native mode the copied
files and `fastgatk-native-inputs.tsv` are preferred by the reader. Native
GenotypeGVCFs accepts
`-V gendb://<workspace>` for such adapter-produced workspaces. For a genuine
GATK TileDB/GenomicsDB workspace (with `callset.json`, `vidmap.json`,
`vcfheader.vcf`, and TileDB array metadata), the optional
`fastgatk-genomicsdb-export` child bridge uses the matching GATK
`libtiledbgenomicsdb.so` in its legacy C++ ABI, exports one VCF shard per array,
and maps canonical `--interval CONTIG:START-END` selectors to inclusive TileDB
column ranges before export. It feeds those shards into the same native Kokkos
genotyping path. The bridge
is isolated by a process boundary; configure a matching
`FASTGATK_GENOMICSDB_LIBRARY`/`FASTGATK_GENOMICSDB_JVM_LIBRARY` pair or let the
pinned GATK package be extracted automatically, and override the executable
with `FASTGATK_GENOMICSDB_BRIDGE` when needed. If the bridge is unavailable,
opaque workspaces still fail closed rather than being guessed at. When the native sidecar
`fastgatk-workspace.json` is present, the reader validates the schema,
`fastgatk-sparse-index` backend tag, input count, and local input-file
existence before expansion. Native import also writes
`fastgatk-record-index.tsv` with per-input contig/span/record counts; interval
queries match contig names (not header-local numeric RIDs) to skip unrelated
shard bodies while still loading every sample header, including disjoint shard
dictionaries. A selector absent from the first shard is retained by canonical
name, and both CSI and TBI indexed traversal map that name to each shard's
index sequence ID. Contigs are merged into one output dictionary and FORMAT arrays
are written only after global sample merge, so disjoint shard dictionaries do not
corrupt RID or sample-width encoding. Legacy index-only adapter workspaces remain accepted for backwards
compatibility. `benchmark_genomicsdb_import.py` reports both external-adapter
and native-sparse-index file-boundary timings.

`verify_genomicsdb_bridge.py` builds a real GATK GenomicsDB workspace, verifies
the isolated export manifest and an adjacent-record single-base pushdown, and
compares native direct-gendb GenotypeGVCFs
GT/AD/PL values with Java GATK. The bridge is intentionally not labeled
bit-identical: GenomicsDB storage/query compatibility is now covered, while
the downstream native annotation/QUAL contract remains the separate oracle.

When the OpenMP build is used, GenotypeGVCFs initialization honors
`OMP_NUM_THREADS` first and otherwise caps Kokkos to `SLURM_CPUS_PER_TASK`
through `ResourceSnapshot::effective_threads()`. This keeps direct dispatcher
replacement inside the scheduler allocation even though GATK exposes no
GenotypeGVCFs compute-thread flag.

The HaplotypeCaller RCM ref/het/non-ref pileup likelihood is covered by the
kernel API smoke, the mixed kernel benchmark (reference-confidence prepare and
kernel timings/checksum), and the bundled chr17 GATK oracle. The oracle requires
11 selected block-start GQ/PL triples and the terminal soft-clip POS→END
boundary now match on the bundled fixture. The JSON report exposes
`rcm_terminal_boundary_exact` plus the GATK/native POS→END maps; broader
realignment/RCM semantics and whole-file bit identity remain separate gates.
The arbitrary-ploidy RCM materializer now stages observations into stable
locus-local segments before the Kokkos `RangePolicy` reduction. It preserves
the original within-locus summation order (so strict PL/GQ bits are unchanged)
while avoiding a full observation-array scan for every locus/genotype; the
polyploid prepare/kernel p50 and checksum remain part of the benchmark contract.
It can merge disjoint sample headers into a joint VCF record, preserving missing
FORMAT values and compacting GT/AD/PL per sample; when shards at one locus carry
different concrete ALT subsets it first forms an allele-name union and remaps
GT/AD/PL/INFO-AD into the joint triangular layout. Diploid GT/GQ are re-derived
by the reusable Kokkos VCF-combination PL kernel (diploid triangular and
arbitrary fixed ploidy), while AC/AN/AF are produced by a
separate deterministic Kokkos allele-count kernel over the merged sample-major
GT vector; both kernel prepare/execute phases and execution spaces are recorded
in the manifest. Reference-only `<NON_REF>` blocks are retained through the
union, projected from REF/<NON_REF> likelihoods into REF/concrete-ALT PL rows,
and emitted as hom-ref samples at a concrete cohort site instead of being
dropped; this projection is structural and deterministic, not a claim of full
GATK cohort posterior equivalence. Repeated `-L`/`--intervals`/`--region`
selectors are applied with GATK `--interval-set-rule UNION` (default) or
`INTERSECTION` before materialization.
Literal, Picard interval-list and BED selectors use the shared Host interval
reader; parsed file counts are recorded in the OutputManifest.

The cross-shard concrete-ALT union and `<NON_REF>` projection use the same
shared Kokkos Number=G PL-remap API as ReblockGVCF and SelectVariants. A `-1`
target-to-source mapping produces missing PL rows for shard-local alleles, and
the kernel sorts mapped source tuples before ranking, so reversed ALT order is
deterministic. `verify_genotype_gvcf.py` covers reordered ALT, multi-sample and
triploid calls; `benchmark_genotype_gvcf.py` reports remap, GT/GQ and cohort
kernel telemetry at the file boundary. The benchmark accepts
`--records-per-shard N` and `--repetitions N`; its default 128-record workload
is deliberately a micro regression/telemetry smoke and emits
`speedup_claim_allowed=false`. It does not run a matched Java/GATK baseline,
so neither the 128-record p50 nor any other output from this script may be
presented as a Java-GATK speedup. Use at least 4096 records per shard, with a
matched Java command and the same input/output scope, before making a
throughput claim. FORMAT/INFO `AD` Number=R projections use the companion
Kokkos allele-field remap API and expose their own lifecycle telemetry.
duplicate sample names are deterministically de-duplicated. `GenomicsDBImport` uses the resource-aware
external-backend adapter described above; `--dry-run` keeps the complete
original argument vector visible.
`ReblockGVCF` has a native sample-major prototype for multi-sample GQ-band block merging,
low-quality hom-ref conversion, complete-PL GQ materialization, arbitrary-
ploidy allele/PL compaction (including triploid Number=G remapping),
PL[0]-based RGQ thresholding (including the GATK long alias
`--rgq-threshold-to-no-call` in separated or inline form, independently of
`--drop-low-quals`), TREE_SCORE no-call conversion,
GATK `HomRefBlock`-compatible element-wise PL minima, rounded per-sample median DP
and per-sample MIN_DP-preserving block merges, optional QUALapprox, and BGZF/Tabix output;
high-quality variant records also emit GATK's `RAW_GT_COUNT` (hom-ref/het/
hom-var) annotation for downstream ExcessHet/QD consumers;
low-quality deletion/span-deletion calls trim the reference block REF to its
leading base while preserving END;
when allele compaction drops an unused concrete ALT, high-quality variants
also apply GATK-compatible reverse common-suffix trimming (symbolic and
spanning-deletion alleles are preserved) and update the covered REF span;
variant FILTERs follow GATK defaults (cleared unless `--keep-site-filters`),
and `--add-site-filters-to-genotype` writes FORMAT/FT;
high-quality records zero `<NON_REF>` AD and subtract the removed count from
FORMAT/DP through a Kokkos cleanup kernel;
stale fixed and dynamic `GVCFBlock*` INFO annotations are removed from
records/header, while
`--annotations-to-keep` and `--annotations-to-remove` provide explicit
caller-field retention/removal;
repeated `-L`/`--intervals` selectors use GATK's `--interval-set-rule UNION`
(default) or `INTERSECTION` semantics and retain a gVCF block when its `END`
span overlaps the selected interval set;
literal, Picard interval-list and BED selectors are parsed by the shared Host
interval reader;
`--create-output-variant-index` follows GATK's optional-boolean syntax (default
true): a bare flag enables Tabix, while `false`/`0` suppresses the `.tbi`
sidecar and is recorded in the OutputManifest; invalid boolean literals fail
closed;
full annotation/posterior semantics and complex multi-sample
`addRefBlockIfNecessary` writer behavior still use explicit fallback; the
single-sample `-R` path is implemented as described above. The bundled chr17
HaplotypeCaller/GATK 4.6.2.0 oracle (`scripts/verify_reblock_gatk_oracle.py`)
matches all normalized record fields; the supported multi-input mode is also
checked by `scripts/verify_reblock_gatk_shards.py` using two non-overlapping
shards from one sample. GATK rejects multi-sample VCF inputs for this tool, so
native multi-sample handling remains an explicitly non-GATK-compatible
prototype boundary; only INFO-key serialization order is normalized by the
comparison.
`SelectVariants` has a native HTSlib streaming prototype for SNP/INDEL/MNP/MIXED
type selection, a deterministic JEXL subset (`QUAL`/INFO comparisons,
`vc.getAttribute`, common `VariantContext` methods (`isSNP`, `isIndel`,
`isMNP`, `isTransition`, `isTransversion`, `isMixed`, `isSymbolic`, `isVariant`,
`isBiallelic`, `isMultiallelic`, `isFiltered`, `isPass`, `hasAttribute`,
`isPolymorphicInSamples`, `isMonomorphicInSamples`,
`getStart`, `getEnd`, `getNAlleles`, `getNSamples`, `getContig`, `getID`,
`getAlleles().size`, `getFilters().size` and `getFilters().contains`,
reference/alternate base strings, `getReference()/getAlternateAllele(i)` allele
methods (`isSymbolic`, `isReference`, `isNoCall`, `isCalled`, `isNonReference`,
`isNonRefAllele`, `isBreakpoint`, `isSingleBreakend`, `length`),
`getType` string/enum comparisons, scalar arithmetic (`+`, `-`, `*`, `/`, `%`)
with Java precedence, and explicit boolean RHS comparisons for
site/allele methods,
`isNotFiltered` and `hasAlternateAllele`), string `contains`/`startsWith`/
`endsWith`/regex methods, INFO-null checks, and `&&`/`||`/`!`), sample subsetting, filtered-record exclusion
and unused ALT/GT/AD/PL compaction with BGZF/Tabix output. AC/AN/AF and GQ are
recomputed after sample/allele subsetting. Repeated site expressions are
OR-combined, `--select` also accepts explicit `vc.getGenotype("S").*`
receivers, and `--select-genotype` supports any-sample GQ/DP/AD/MIN_DP
comparisons plus indexed `getAD()[i]`/`getPL()[i]` and
`isHet`/`isHom`/`isHomRef`/`isHomVar`/`isCalled`/`isNoCall`/`isAvailable` predicates;
the numeric predicate spellings used by GATK (`isHet == 1`, `isHomRef == 1`,
and similar) are also lowered per sample;
site-level `hasGenotypes()`, `getCalledChrCount()`, `getNoCallCount()`,
`getHomRefCount()`, `getHetCount()`, and `getHomVarCount()` are evaluated from
the actual VCF GT vector (including arbitrary fixed ploidy);
`--concordance`/`--discordance` now provide deterministic contig/POS/REF/ALT
site-key filtering. Explicit `--concordance-genotypes`/
`--discordance-genotypes` now compare unphased GTs for shared samples using the
native any-shared-sample subset. Exact Java sample-set aggregation/edge cases
INFO vector operands accept GATK's `vc.getAttribute("TAG").get(i)` syntax and
the equivalent `TAG[i]` shorthand for Number=A/Number=R numeric arrays. Indexed
vectors remain direct comparison operands (SelectVariants/HTSJDK rejects vector
arithmetic); arithmetic on scalar site fields and per-sample GQ/DP is covered by
the native/GATK oracle. Exact Java sample-set aggregation/edge cases and the
remaining full JEXL surface remain explicit fallback.
Repeated `-L`/`--intervals` regions use GATK `--interval-set-rule UNION`
(default) or `INTERSECTION` before filtering and preserve
deterministic record order; literal selectors, Picard interval-list files and
zero-based BED files (plain or `.gz`) share the Host interval parser, with parsed file counts in
the OutputManifest.
`--exclude-non-variants` follows GATK `VariantContext.isVariant()` for
sample-bearing records (hom-ref/no-call-only sites are removed, called ALT and
sample-free ALT records are retained). `--min-indel-size` and
`--max-indel-size` filter indels by absolute REF/ALT length difference. A pure
site-selection pass preserves existing AC/AN/AF INFO; those annotations are
recomputed only after sample/allele subsetting, matching the GATK boundary and
covered by the native/GATK oracle.
The `--exclude-non-variants` and `--exclude-non-variant-sites` switches also
accept GATK/Barclay optional-boolean spellings: bare (true), `true`/`1`,
`false`/`0`, and inline `--option=value`. The focused
`scripts/verify_select_variants_gatk_oracle.py` test compares these forms with
GATK 4.6.2.0 and runs in both OpenMP and Serial CTest configurations; broader
JEXL and sample-set edge cases remain at the prototype/fallback boundary above.

SelectVariants also supports GATK 4.6.2.0's
`--sites-only-vcf-output` writer option. Selection, sample/allele remapping and
genotype predicates consume the complete FORMAT/sample payload; only the final
writer header and record are reduced to an 8-column site-only VCF. Bare and
separated true/false forms match GATK, while embedded `=` is rejected like
Barclay. `fastgatk-select-variants-sites-only-gatk-oracle` compares site rows
and the full-output false path against pinned GATK on OpenMP and Serial, and
checks the manifest compatibility/telemetry field.
`VariantsToTable` has a native HTSlib streaming prototype for the common GATK
table contract. Repeated `-F` site/INFO fields (including GATK wildcard fields
such as `A*`), `-GF` FORMAT fields, `-ASF` and
`-ASGF` allele-specific fields are rendered with deterministic sample-column
ordering; missing values are `NA`, filtered records are skipped by default,
`-SMA` emits one row per alternate allele, `-XL/--exclude-intervals` removes
overlapping records, and `--moltenize` emits the standard
`RecordID/Sample/Variable/Value` form. Repeated `-L`/`--intervals` selectors
(with GATK-compatible `-ip`/`-ixp` padding and `-XL` exclusions) use
the shared literal, Picard interval-list and plain/`.gz` BED reader, including
UCSC `track`/`browser` headers; selection is record-span overlap based, so a
gVCF `END` block is retained when it covers the requested locus. Java/GATK argument aliases are registered in
the dispatcher and an OutputManifest records field extraction, skipped rows and
compatibility limits. Standard VariantContext getters (`HET`, `HOM-REF`,
`HOM-VAR`, `NO-CALL`, `VAR`, `NSAMPLES`, `NCALLED`, `SAMPLE_NAME`, `TYPE`,
`EVENTLENGTH`, `TRANSITION`) and Java allele-valued GT rendering are covered by
a pinned GATK 4.6.2.0 table oracle; split Number=A/R fields follow the ASGF
rules while Number=G fields remain intact. QUAL missing values render as
GATK's `-10.0`, spanning-deletion `*` follows VariantContext type rules, and
unsupported namespace-prefixed/standard fields remain `NA` instead of acquiring
native-only semantics. Text-VCF INFO tokens are retained at the Host boundary,
so long decimals, exponent spelling, literal `.` and INFO flags match GATK;
typed FORMAT/QUAL values use the shared Java-style shortest serializer. Typed
BCF and every extreme Java floating-point/report edge are not promised to be
raw-bit-identical; complex expressions, cloud staging and release-specific
annotation expansion remain explicit fallback.
The focused `verify_variants_to_table_gatk_oracle.py` adds byte-level evidence
for GATK's split `-ASF` Number=R formatting (including its observable leading
spaces), biallelic `-SMA`, and `-EMD` presence semantics (a literal String
`NA` is not missing). It also preserves the release's unusual molten writer:
it sequentially consumes extracted values while emitting only `-F`/`-GF`
labels, so requested `-ASF`/`-ASGF` values shift later molten cells instead of
having dedicated molten rows.
The shared interval arguments `--interval-merging-rule ALL|OVERLAPPING_ONLY` and
`--variant-output-filtering STARTS_IN|ENDS_IN|OVERLAPS|CONTAINED|ANYWHERE` are
also accepted and recorded in the manifest. The latter follows the pinned GATK
table-writer behavior: after validating the mode and requiring `-L` for
non-`ANYWHERE`, traversal remains ordinary record-span `OVERLAPS` selection
because `VariantsToTable` writes a table directly rather than through
`IntervalFilteringVcfWriter`.
`VariantEval` has a native Host/HTSlib prototype for the high-value evaluator
subset: `CountVariants`, `TiTvVariantEvaluator`, `VariantAFEvaluator`,
`ThetaVariantEvaluator`, `MendelianViolationEvaluator`, `CompOverlap` and
`GenotypeConcordance`. It now emits deterministic `Contig`, `Filter`,
`VariantType`, `AlleleFrequency`, comparison-backed `Novelty`, and per-sample plus aggregate `Sample`
stratification tables; GATK `-ST/--stratification-module Sample` aliases are
accepted alongside `-S`. Per-sample `CountVariants` rows preserve hom-ref,
het, hom-var, no-call and singleton state and are checked against the pinned
GATK 4.6.2.0 oracle. It accepts repeatable `-eval`/`-comp` tracks and GATK-style
`-EV`/`-ST`/`-S` selectors, and records concordance/stratification telemetry in
OutputManifest; `-L` selection uses record-span/`END` overlap for gVCF blocks.
Repeated `-L`/`--intervals` selectors support GATK `UNION` (default) or
`INTERSECTION` via `-isr/--interval-set-rule`; eval and comparison tracks use
the same Host-side interval set and the selected rule is recorded in the
summary/manifest.
`CountVariants` follows GATK's genotype-aware ignore-AC0 locus semantics:
when `FORMAT/GT` is present, hom-ref and all-no-call ALT records are reference
loci, while sites-only records retain their site-level variant state. The
pinned GATK 4.6.2.0 oracle covers both genotyped and sites-only CountVariants
inputs, including singleton and insertion/deletion ratio fields.
`--keep-ac0` (and `-keep-ac0`) restores genotyped AC=0 SNP/INDEL records to the
native CountVariants aggregate path. A pinned GATK 4.6.2.0 oracle compares that
aggregate report byte-for-byte; native-only checks cover bare/true/false and the
short alias without claiming parity for other evaluator/stratifier tables.
`VariantAFEvaluator` implements GATK's default AC=0 exclusion and diploid
0/0.5/1.0 genotype AF accumulation; its five GATKReport fields are compared
directly against GATK 4.6.2.0. `ThetaVariantEvaluator` implements the
polymorphic-SNP heterozygosity, pairwise-difference and harmonic region
estimate fields, also with a direct GATKReport oracle. `-no-ev` now suppresses
standard evaluator tables unless they are explicitly selected.
`MendelianViolationEvaluator` accepts `--pedigree/-ped` and the default
`--mendelian-violation-qual-threshold 50`, emits all 31 GATK fields, and is
checked against a pinned trio oracle; missing pedigree relationships remain
explicitly skipped rather than guessed.
Explicit `-no-ev -EV ValidationReport -no-st` now implements GATK's
comparison-driven site-accuracy matrix. Host decode applies GT and INFO/AC
mono/poly classification plus comparison-genotype subsetting to the eval
sample set, then a Kokkos integer reduction accumulates the 4x4 status matrix.
The pinned GATK 4.6.2.0 oracle compares complete GATKReport bytes for genotyped
and sites-only inputs and fixes the default pre-evaluator removal of filtered
records. Multiple tracks, duplicate eval loci, and Filter-stratified
ValidationReport remain fail-closed fallback rather than inferred semantics.
`-ST/-S Family` now reuses the same PED membership to subset CountVariants
and emits `all` plus per-family rows, with numeric rate columns compared to
the pinned GATK Family-stratification oracle.
`-no-st -ST Novelty -EV CountVariants --gatk-report` emits the GATK-shaped
`all`/`known`/`novel` CountVariants rows using the same Host comparison matcher;
the manifest records known/novel record counts. This is a bounded comparison
track path, not a claim that every Java evaluator/stratifier is complete.
The release-pinned `verify_variant_eval_gatk_oracle.py` adds a narrower
byte-level gate: `CountVariants` aggregate and `Sample` stratification reports
match GATK 4.6.2.0 including GATKReport format metadata, Java numeric
serialization, and column padding. It also requires truncated VCF input to
fail closed. This byte-level promise is intentionally limited to the
`CountVariants` slice; other evaluator/stratifier tables and remote feature
inputs remain explicit fallback.
The native report remains `prototype`: the complete Java
plugin graph, exact formatting for every evaluator/stratifier and remote
feature inputs remain explicit fallback.
`ValidateVariants` has a native Host/HTSlib/FAIDX validation prototype. It
preserves the usual `-V/-R/--dbsnp/-L` surface, validates REF/reference
agreement, ALT uniqueness, AN/GT chromosome counts, INFO/FORMAT Number=A/R/G
and fixed-cardinality fields (including cross-sample GT ploidy), optional dbSNP IDs and
GVCF `<NON_REF>` ordering/overlap, and supports excluded validation types,
filtered-record skipping and warning mode. GVCF interval selection is
record-span based, so an `END` block covering a requested locus is validated
even when its POS is before the interval. Structured validation issues and
resource/interval telemetry are written to OutputManifest. Reference sequence
dictionary contig/length mismatches and reference-backed `--validate-gvcf`
per-contig/interval-union coverage (including gaps between concrete variant
records), leading/trailing coverage or `END` overrun fail closed and can only be bypassed
for the dictionary check with the explicit
`--disable-sequence-dictionary-validation` switch; full Tribble/htsjdk
validation and release-specific GVCF annotation semantics remain explicit
fallback. `-XL/--exclude-intervals` is applied after `-L/--intervals`, and
`-ip/-ixp` expand include/exclude intervals before traversal; repeated
`-L/--intervals` additionally support GATK `UNION` (default) or
`INTERSECTION` via `-isr/--interval-set-rule`; the behavior is
covered by a GATK command-level oracle and `benchmark_validate_variants.py`.
`GetPileupSummaries` now has a native HTSlib/Kokkos prototype for the common
contamination workflow. It reads SAM/BAM/CRAM through the shared Host reader,
selects AF-qualified biallelic SNP sites from `-V` (including a VCF `-L`
selector), applies the GATK default mapping/duplicate/primary/vendor/CIGAR
mask, and counts ref/ALT/other bases with a Kokkos kernel. When `-L` resolves
to coordinate or VCF-site selectors and the population VCF has a `.tbi`/`.csi` index, site
traversal uses HTSlib indexed iterators (including tabix sequence-name mapping
when VCF header `rid` values differ from TBI tids); unindexed resources retain
deterministic sequential fallback and the manifest records the traversal
mode/query count.
Repeated `-L`/`--intervals` selectors support GATK `UNION` (default) or
`INTERSECTION` via `-isr/--interval-set-rule`; the Host interval set is shared
by population-site selection and indexed traversal, and the rule is recorded
in the summary/manifest.
Supported GATK
`-RF/--read-filter` and `-DF/--disable-read-filter` names are validated explicitly,
`--disable-tool-default-read-filters` starts from an empty set, and
`--disable-sequence-dictionary-validation` is honored. The output table uses
`ReadLengthReadFilter` through typed `--min-read-length`/`--max-read-length`
bounds; filtered reads are removed before base projection and the effective
range is recorded in OutputManifest telemetry. The ref/ALT, observation, and
count Kokkos Views are persistent across read batches and grow only when
capacity is exceeded; allocation/reuse counts are included in the manifest so
prepare overhead is auditable. The output table uses the
GATK `SAMPLE` metadata and six-column schema. Text VCF AF tokens retain
GATK/htsjdk binary64 parsing and shortest round-trip output (typed BCF retains
HTSlib's float decode); the bundled GATK 4.6.2.0 oracle
(`scripts/verify_get_pileup_gatk_oracle.py`) is byte-exact on long-decimal,
scientific-notation and read-filter fixtures. Full LocusWalker pileup/filter,
cloud, dictionary edge cases and broader Java report semantics remain explicit
fallback.
`GoodCigarReadFilter` is evaluated per record through the shared Host
`ReadBatch::cigar_record_layout_valid()` API; reference-span/projection helpers
use the same record-local contract, so a malformed CIGAR cannot discard valid
sibling reads in one decoded batch. A zero-length-CIGAR mixed-batch regression
is covered by both OpenMP and Serial contracts.
When indexed region streaming reaches a tile with no alignments, the shared
filter kernel treats the empty batch as vacuously passing instead of requiring a
nonexistent CIGAR payload; this keeps all-contig `--stream-by-region` traversal
valid under the same default mask.
When AF-bearing population records are present but none pass the requested
AF/interval or biallelic-SNP selection, the native writer follows GATK and
successfully emits metadata plus the six-column header with zero data rows;
missing AF input still fails closed.
The reads header must contain an `@RG SM` sample tag, matching GATK's metadata
contract; a missing sample fails closed instead of emitting `SAMPLE=UNKNOWN`.
`CalculateContamination` consumes the same GATK-shaped pileup table and emits
the standard `sample/contamination/error` table. Its native model now follows
GATK's three-round MAF/contamination learning, genotype posterior and
hom-alt/hom-ref strategy cascade; segment likelihoods use the shared
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`/`RangePolicy`
lifecycle, and the standard-error binary search is reproduced.
Four bundled GATK 4.6.2.0 fixtures are checked against the Java oracle. The
changepoint stage now uses the seeded kernel approximation, deterministic
Commons-Math/JAMA-style bidiagonal SVD, persistence and backward-selection parameters used by Java
`KernelSegmenter`; Commons Math 3.5 elementwise SVD parity and cloud/report bit
format parity remain explicit follow-up gates, so the registry remains
`prototype`. Coverage filtering now matches Apache Commons Math
`Mean.evaluate(double[])`: Java's ordinary left-to-right sum followed by a
residual correction pass, with strict low/high thresholds. The manifest records
covered-site count, corrected mean, median, and both thresholds for auditability.
The manifest additionally records kernel execution space/policy, batch and
observation counts, and prepare/execute timing. The native watershed tie-break also follows Java `PersistenceOptimizer` by
retaining the lower-index minimum when equal minima meet. LikelihoodPlan Views
are persistent across the MAF/contamination rounds (and optional segmentation
re-learning), with allocation/reuse counters recorded in the manifest and
benchmark. Pileup tables may be plain or `.gz`; a `.gz` output or segmentation
suffix selects zlib streaming while preserving the decompressed table text.
The segmentation sidecar preserves Java `Collectors.groupingBy`/`HashMap`
contig bucket order and all floating-point fields use binary64 shortest-round-trip
serialization; final libm/optimizer ULP differences remain explicit in the
manifest rather than being presented as bit identity.
`GatherPileupSummaries` provides the scatter/gather half of that workflow. It
orders non-empty shards by the first record under a Picard `.dict`, `.fai`, or
FASTA contig order, then preserves each shard's internal row order, matching
GATK's `PileupSummaryFileComparator`; empty shards are ignored and duplicate
loci are retained by default. `--reject-overlaps` is an explicit stricter
extension, while `--allow-overlaps` remains a compatibility flag. The bundled
GATK 4.6.2.0 oracle (`scripts/verify_gather_pileup_gatk_oracle.py`) is byte-exact
on the multi-shard fixture. The gathered table and OutputManifest are readable
by the native contamination tools and retain the GATK six-column contract.
Contig lengths and site positions are checked against the selected sequence
dictionary; unknown contigs and out-of-range positions require the explicit
`--disable-sequence-dictionary-validation` override.
Plain and gzip-compressed (`.gz`) tables can be mixed as inputs or selected for
the output; the manifest records compressed input count, compressed-output mode,
and completion status. Cloud URI staging and full Java dictionary edge semantics
remain explicit fallback.
AF values are serialized with a shortest binary64 round-trip representation,
matching htsjdk/Java output for high-precision values rather than truncating to
fixed stream precision.

The reproducible file-boundary benchmark for these three tools is
`fastgatk-native/scripts/benchmark_pileup_validation.py`; it reports warmup,
p50/p95 wall time, output/manifest bytes and the Kokkos telemetry for each
case.

`ReblockGVCF` uses the shared `fastgatk-kernels` genotype API for Number=G PL
allele-subset remapping and the companion Number=R AD remap, including fixed
non-diploid ploidy; both APIs now consume all sample rows in one sample-major
batch, preserve per-sample `MIN_DP`/GQ bands during reference-block merges, and
materialize `RAW_GT_COUNT` across called samples. The native tool records remap
prepare/execute/execution-space telemetry and has a reproducible file-boundary
benchmark in `fastgatk-native/scripts/benchmark_reblock_gvcf.py`; the non-reference
AD cleanup kernel also uses the shared `HostBatch -> KernelPlan.prepare -> Kokkos
Views -> execute -> collect` `RangePolicy` lifecycle and reports its own
prepare/execute telemetry; annotation/posterior and
`QUALapprox`/`AS_QUALapprox` now use the same biallelic Kokkos PL projection;
`RAW_GT_COUNT` is materialized from the final called GT and its annotation
count is recorded in the OutputManifest. With `--do-qual-approx`, the native
path also emits GATK's `VarDP` and raw `AS_VarDP` depth annotations from called
variant genotypes;
high-quality ALT compaction now also performs reverse common-suffix trimming
with GATK's symbolic-allele rules; full annotation/posterior parity and the
reference-backed deletion-gap bridge now use the shared Kokkos PL projection;
overlap-aware writer semantics and full annotation/posterior parity remain
explicit fallback. GATK 4.6.2.0 accepts only non-overlapping shards from one
sample for this tool and rejects multi-sample VCF inputs; the native sample-major
multi-sample path is therefore a prototype-only extension, not a direct-replacement
claim.
`SelectVariants` uses the same Kokkos PL-remap and GT/GQ derivation path after
sample/ALT subsetting and Kokkos allele-field remap for AD; the file-boundary baseline is
`fastgatk-native/scripts/benchmark_select_variants.py` and its manifest records
both kernel lifecycles.
`GatherVcfs` has a native streaming prototype for coordinate-sorted shard
concatenation, later-header definition merge, sample-order validation and
BGZF/Tabix or plain VCF/Tribble `.idx` output. Default mode requires disjoint sorted intervals, including
gVCF `END` spans; `--allow-overlaps` is explicit. `-L`/`--intervals` now accepts
literal intervals, Picard `.interval_list`/`.intervals`/`.list` files and
zero-based half-open `.bed`/`.bed.gz` files (including UCSC `track`/`browser`
headers); parsed interval-file counts are recorded in
the manifest. Repeated `-L`/`--intervals` selectors support GATK
`UNION` (default) or `INTERSECTION` via `-isr/--interval-set-rule`, with the
selected rule recorded in summary/manifest. Record selection is span-overlap based, so a gVCF block whose POS
precedes `-L` is retained when its `END` covers the requested locus. Cloud URI
staging remains explicit fallback for the direct binary (the dispatcher stages
HTTP(S)/presigned inputs and commits outputs). `--reorder-input-by-first-variant`
(`-RI true`) performs a stable first-record sort using the first shard's
sequence dictionary, matching GATK's `REORDER_INPUT_BY_FIRST_VARIANT`; the
manifest and stdout telemetry record the mode. Contig names/order/lengths are
validated across shard headers; `--disable-sequence-dictionary-validation` is
an explicit escape hatch and is recorded in the manifest. Repeated regions are
applied before the coordinate contract. The manifest also records strict Host
execution, output/index completeness, output byte counts and wall time. GATK's
`-CO/--COMMENT` is emitted as `##GatherVcfs.comment=...` metadata; the native
`--CREATE_INDEX` alias controls `.tbi`/`.idx` creation and `--COMPRESSION_LEVEL` uses
the same 0–9 BGZF level range (default 2) as Picard. The repeatable
`benchmark_gather_vcfs.py` baseline measures shard/record throughput and checks
output/index completeness.
Native `.bcf` output is fail-closed: HTSlib emits BCF 2.2 while the pinned
GATK/htsjdk GatherVcfs reader requires the BCF 2.1 form produced by GATK, so use
the dispatcher's explicit `--fallback` path for BCF output rather than creating
an artifact that a downstream GATK stage cannot read.
`LeftAlignAndTrimVariants` has a reference-backed HTSlib prototype for common
prefix/suffix trimming and repeat left-shifting of primitive non-symbolic
alleles, plus deterministic multiallelic splitting with GT/AD/PL/GQ
compaction, split-record AC/AN/AF recomputation, and optional preservation of
AC/AF/AN as AC_Orig/AF_Orig/AN_Orig. Split PL and AD now use the shared Kokkos
genotype/allele-field remap kernels for arbitrary fixed ploidy (including
triploid), normalize PL by each sample minimum, project removed ALTs to REF as
GATK does, and record prepare/execute/execution-space telemetry. Repeated `-L`/`--intervals` selectors
use GATK `--interval-set-rule UNION` (default) or `INTERSECTION` before
normalization and record-span/`END` overlap; literal,
Picard interval-list and BED selectors use the shared Host interval reader. The
GATK's `--max-indel-length` and `--max-leading-bases` bounds are enforced
before/inside the repeat walk and recorded in the manifest, including
oversized-indel and left-shift-base telemetry. The file-boundary benchmark is
`fastgatk-native/scripts/benchmark_left_align.py`.
多等位记录含 `*`、`<...>`、BND 或 single-breakend symbolic ALT 时，split path
现在按 HTSJDK/GATK 的 fail-closed 语义逐 ALT 输出 no-call GT，清除不可安全
投影的 FORMAT/INFO，并在 OutputManifest 记录清除计数。复杂 symbolic 的
reference normalization 和完整 annotation/masking surface remain explicit
fallback。
`VariantFiltration` has a native deterministic JEXL subset for QUAL/INFO and
FORMAT comparisons, `vc.getAttribute(...)`, `vc.getGenotype(...).isHet()` /
`isHomRef()` / `isHom()` / `isHomVar()` / `isNoCall()`, `getDP()`/`getGQ()`
`getMIN_DP()` and indexed `getAD()[i]`/`getPL()[i]`, scalar arithmetic
(`+`, `-`, `*`, `/`, `%`) with Java precedence, common site
methods (`isSNP`, `isIndel`, `isMNP`, `isMixed`, `isSymbolic`, `isVariant`,
`isBiallelic`, `isMultiallelic`, `isFiltered`, `isPass`, `hasAttribute`,
`isPolymorphicInSamples`, `isMonomorphicInSamples`,
`isNotFiltered`, `isNoVariation`, `hasAlternateAllele`, `getStart`, `getEnd`,
`getNAlleles`, `getNSamples`, `getAlleles().size`, `getFilters().size`,
`getFilters().contains`, `getType`, contig/ID/reference/alternate string
accessors, string contains/startsWith/endsWith/regex methods and `=~`/`!~`
regex operators, enum comparisons,
`getReference()/getAlternateAllele(i)` allele methods
(`isSymbolic`, `isReference`, `isNoCall`, `isCalled`, `isNonReference`,
`isNonRefAllele`, `isBreakpoint`, `isSingleBreakend`, `length`), and
INFO/QUAL null checks), explicit boolean RHS comparisons (`== true/false`), and boolean
`&&`/`||`/`!`, lexical FILTER label ordering and BGZF/Tabix output. Explicit
`--genotype-filter-expression`/`--genotype-filter-name` rules now emit
genotype-level `FT`; site/genotype filter inversion is supported, while mask
VCF/BCF interval filtering (`--mask`, `--mask-name`, `--mask-extension`) is
supported, including `--filter-not-in-mask` reverse selection and the
`--mask-description` FILTER header override. Repeated `-L/--intervals`
selectors use GATK `--interval-set-rule`
`UNION` (default) or `INTERSECTION` and record-span/`END` overlap,
so a gVCF block covering a requested locus is not lost because its POS is
outside the interval. Deterministic clustered-event filtering is available through
`--cluster-size`/`--cluster-window-size`; the same genotype availability/count
methods as SelectVariants are supported; complete JEXL method and
allele-specific annotation coverage remains explicit fallback except for the
common numeric `AS_*` `Number=A` subset enabled by
`--apply-allele-specific-filters`, which writes per-ALT `AS_FilterStatus`.
Numeric INFO-vector operands additionally accept GATK's
`vc.getAttribute("TAG").get(i)` and the compatible `TAG[i]` shorthand, with
out-of-range elements treated as missing; vector arithmetic and all supported
genotype scalar arithmetic are covered by the GATK oracle.
Repeated
`-L`/`--intervals` regions are applied before filtering; literal selectors,
Picard interval-list files and BED files use the shared Host interval parser.
Multiple matching genotype rules are composed with stable semicolons, and
arbitrary VCF ploidy is handled by the genotype predicates.
The genotype expression parser also accepts GATK's compact per-sample forms
(`GQ < 20`, `DP == 0`, `AD[1] == 0`, `PL[2] < 10`, and numeric predicates such
as `isHet == 1`) and binds them to every sample in the record; explicit
`vc.getGenotype("SAMPLE")` expressions remain supported and continue to win
over the implicit sample context.
The focused `scripts/verify_variant_filtration_gatk_oracle.py` additionally
compares a pinned GATK 4.6.2.0 run for QUAL/INFO/FORMAT/GT/FT/FILTER records,
including site-filter inversion, on both OpenMP and Serial. The comparison
normalizes only writer-level INFO numeric spelling, FORMAT key order, trailing
missing FORMAT values, and FILTER label order; all semantic values remain
checked. Site filtering preserves pre-existing non-PASS FILTER labels when a
new rule is appended, matching chained GATK VariantFiltration stages. Full
JEXL/annotation-engine and cloud semantics remain explicit fallback.
The common collection accessors `vc.getGenotypes().size()`,
`vc.getGenotypes().isEmpty()`, and `vc.getFilters().isEmpty()` are also
implemented; the VCF spanning-deletion allele `*` follows HTSJDK symbolic
semantics in both SelectVariants and VariantFiltration.
All native optional booleans now use the shared GATK parser (HC, BQSR,
Mutect2, GenotypeGVCFs, ReblockGVCF, VQSR/ApplyVQSR, SortSam,
MarkDuplicates, GatherVcfs, DepthOfCoverage and the VCF tools):
bare/`true`/`1` enables a switch, `false`/`0` disables it, and an unknown or
empty literal fails before input staging. Output-index switches therefore
have one consistent `.tbi`/`.bai` sidecar contract; the registry and
`fastgatk-optional-boolean-contract` exercise both inline and next-token forms.
`SortSam` has a native external-memory prototype for coordinate/queryname/duplicate order,
bounded record staging, spill runs, deterministic k-way merge and BAM/CRAM
index output. The duplicate order follows the pinned Picard/GATK 4.6.2.0
`SAMRecordDuplicateComparator`: library, unclipped read and MC-derived mate
coordinates, orientation, mapped-end state, mapped reference-length score and
read-name/end tie breaks. Mapped pairs therefore require the Picard `MC` tag
contract. With explicit `--resume-spill`, a completed input pass is
checkpointed in scratch; an output-stage failure can be retried against the
same input/output/sort/budget contract and reuses spill runs only after input
size/mtime and run-completeness validation. The
checkpoint is removed only after a complete output manifest. `--add-pg-tag`
adds a configurable SAM `@PG` record (ID/name/version/command line) without
changing the default header. The pinned Picard/GATK 4.6.2.0 SortSam oracle
now covers HTSJDK coordinate/queryname tie-break ordering, the output `@HD SO`
declaration, and the Picard `CREATE_INDEX=false` default versus explicit
index creation. A separate pinned oracle covers duplicate ordering with bounded
multi-run spill/merge on both OpenMP and Serial builds. Cloud I/O, compression
identity and complete Picard provenance remain explicit fallback.
`MarkDuplicates`
has a native two-pass metadata prototype:
the first pass builds coordinate-sorted canonical duplicate-key/fragment
summaries and the second pass streams the original records, so mates separated
by other coordinate records still receive one quality-based consensus. It
supports `BAM_FDUP` marking or removal, read-group/library partitioning,
canonical paired-fragment consensus (both mates are marked as one duplicate
unit), optical duplicate distance detection, optional `DT:Z:SQ/LB` tagging, a
Picard-shaped metrics table and BAM/CRAM index output. With `--add-pg-tag`,
the output header also carries a configurable `@PG` ID/name/version/command
line. When `--tmp-dir` is
provided, duplicate-key summaries spill to sorted binary runs and are looked
up from disk during the second pass; large individual duplicate keys are
chunked across runs and merged by key, so the record metadata bound applies
to the worst single group as well. Spill files are cleaned after a complete
manifest. The first pass also publishes a versioned checkpoint containing an
input size/mtime fingerprint, effective memory limit, counters, library metrics,
and the complete spill-run index. If preemption or a resource failure
interrupts the second pass, rerun with the same options plus `--resume-spill`;
the checkpoint is validated before any run is opened and is removed only after
the output, metrics, index, and manifest all pass completeness checks. The
Picard-shaped metrics now include a deterministic estimated
library size for paired reads (with optical pairs excluded exactly as Picard
does), one deterministic row per read-group library (`LB`), and the overlaid
duplicate/optical/non-optical set-size histogram. PG-line/cloud and the
remaining release-specific Picard metrics semantics remain explicit fallback.
Mapped mates are paired by read name during the first pass and their actual
unclipped five-prime coordinates are used for one ReadEnds key; this covers
missing `MC` tags and terminal S/H clipping. Repaired pair keys are included in
the spill checkpoint so `--resume-spill` preserves the same duplicate decision.
The pinned four-record no-MC/clipped Picard oracle passes on both backends.
MarkDuplicates defaults to Picard's `ADD_PG_TAG_TO_READS=true` and
`CREATE_INDEX=false`; pass `--add-pg-tag=false` to suppress the native program
record, or explicitly pass `--create-output-bam-index=true` to create a BAM index.
Both tools parse `--add-pg-tag` and output-index switches as GATK optional
booleans (bare/true/1 or false/0); invalid literals fail before opening the
input or creating spill state.
The duplicate representative follows Picard's default
`SUM_OF_BASE_QUALITIES` exactly (MAPQ is not mixed into the score), and optical
duplicate classification is limited to paired fragments. The regression
compares duplicate flags and `DT` tags against the pinned Picard 4.6.2.0 oracle
when the local GATK jar is available. `--tagging-policy OpticalOnly` is also
supported: only optical duplicate records receive `DT:Z:SQ`, while ordinary
library duplicates remain untagged; a dedicated 14-record Java/native oracle
and CTest cover this boundary on both OpenMP and Serial builds.

```bash
fastgatk-native/dispatcher/fastgatk --version
fastgatk-native/dispatcher/fastgatk --list
fastgatk-native/dispatcher/fastgatk --dry-run HaplotypeCaller \
  -I reads.bam -R reference.fa -L 17:69000-69100 -O calls.vcf
python3 fastgatk-native/dispatcher/verify_dispatcher.py
python3 fastgatk-native/scripts/compute_progress_score.py
python3 fastgatk-native/scripts/verify_fixture_digests.py
python3 fastgatk-native/scripts/verify_bqsr.py
python3 fastgatk-native/scripts/verify_genotype_gvcf.py
python3 fastgatk-native/scripts/verify_gatk_genotype_gvcf.py
python3 fastgatk-native/scripts/verify_reblock_gvcf.py
python3 fastgatk-native/scripts/verify_select_variants.py
python3 fastgatk-native/scripts/verify_variants_to_table.py
python3 fastgatk-native/scripts/verify_variants_to_table_gatk_oracle.py
python3 fastgatk-native/scripts/benchmark_variants_to_table.py
python3 fastgatk-native/scripts/verify_gather_vcfs.py
python3 fastgatk-native/scripts/verify_left_align.py
python3 fastgatk-native/scripts/verify_variant_filtration.py
python3 fastgatk-native/scripts/verify_variant_filtration_gatk_oracle.py
python3 fastgatk-native/scripts/verify_sort_sam.py
python3 fastgatk-native/scripts/benchmark_sort_sam.py
python3 fastgatk-native/scripts/verify_mark_duplicates.py
python3 fastgatk-native/scripts/benchmark_mark_duplicates.py
python3 fastgatk-native/scripts/verify_genomicsdb_import.py
python3 fastgatk-native/scripts/verify_kokkos_backend_matrix.py
bash fastgatk-native/scripts/verify_all.sh
```

`verify_all.sh` is the reproducible local gate: it builds (unless
`FASTGATK_SKIP_BUILD=1`), runs the full OpenMP CTest and a lightweight focused
Serial CTest by default (including score and fixture-digest gates), executes the
GATK oracle when available (or fails when `FASTGATK_REQUIRE_GATK_ORACLE=1`), and
then runs the VQSR/CNV file-boundary verify+benchmark suites and the
Nextflow/SLURM local workflow contract. Set `FASTGATK_VERIFY_SERIAL_FULL=1` to
run the complete current Serial matrix. A real SLURM allocation and GPU backend
remain CI/environment-specific gates.

The registry currently marks the native HaplotypeCaller, CountBasesInReference,
CountReads, FlagStat, IndexFeatureFile, CompareReferences, CheckReferenceCompatibility,
AnnotateIntervals, FastaReferenceMaker, FastaAlternateReferenceMaker, ShiftFasta,
SplitIntervals, FilterIntervals, PreprocessIntervals, GenotypeGVCFs and GatherVcfs as `contract-compatible`: HC still has the limited
assembly/posterior boundary described above, while CountBasesInReference has
literal, BED, interval-list, gzip-selector, overlap-union, malformed-input,
and bundled GATK 4.6.2.0 oracle coverage on OpenMP and Serial. The registry
state may only be promoted after the per-tool oracle, output, resource, and
Nextflow/SLURM gates pass. The versioned weighted snapshot and its arithmetic
check are in [progress_score.json](</home/turing-agents/Documents/fast-gatk/progress_score.json>)
and `fastgatk-native/scripts/compute_progress_score.py`.

If HTSlib is not found, the target still builds as a core-only binary but exits with
`BACKEND_UNAVAILABLE` for BAM/CRAM instead of silently using a text approximation.

The latest bounded GATK 4.6.2.0 gates also cover VariantsToTable moltenize/allele-specific
field ordering and `Number=R` splitting, CalculateContamination's all-sites-filtered
`0.0/1.0` plus header-only segmentation and `-matched/-segments` aliases, and
GatherTranches' required `--mode`, `-tranche` alias, and version-6 empty-scatter
contributor behavior. These are scoped semantic slices; they do not imply complete
algorithmic or bit-identical replacement. The active CTest inventory is now 140 entries.

CollectAllelicCounts now resolves the pinned GATK 4.6.2.0 read-filter boundary: Wellformed
and Mapped from LocusWalker plus Mapped, NonZeroReferenceLength, NotDuplicate and
MappingQuality>=30 from the tool. Secondary/supplementary/QC-fail flags remain retained by
default. `--disable-read-filter/-DF`, `--read-filter/-RF`, and
`--disable-tool-default-read-filters` are parsed after argument collection with GATK's
default-removal then explicit-filter-addition ordering; Wellformed rejects CIGAR `N`.
The independent pinned oracle compares default, all-defaults-disabled, and explicit
NotDuplicate/MAPQ re-enable cases on OpenMP and Serial. CollectAllelicCounts remains a
prototype for multi-sample/cloud/full malformed-input and complete CNV workflow semantics;
the shared CTest inventory is now 140 entries.

The 2026-09-01 bounded compatibility increment adds a pinned GATK 4.6.2.0
BI/BD indel-quality oracle for BaseRecalibrator (all four report tables and
malformed-tag rejection), SelectVariants `--set-filtered-gt-to-nocall` with
post-compaction AC/AN/AF refresh, and CollectAllelicCounts default/read-filter
ordering (`-DF`, `-RF`, and `--disable-tool-default-read-filters`). All three
oracles pass on OpenMP and Serial; these checks do not promote the tools beyond
their registered prototype status or imply complete algorithmic parity.

The same increment adds an ApplyBQSR `--preserve-qscores-less-than` oracle:
17,748 bases below Q10 remain unchanged through native and pinned Java
`PrintReads`, while higher-quality bases continue through the Bayesian path.
The gate passes on both OpenMP and Serial and is limited to this quality-preserve
boundary rather than a claim of full BQSR model parity.

The 2026-09-01 CallCopyRatioSegments bounded slice preserves GATK 4.6.2.0
SimpleCopyRatioCaller IEEE-754 behavior for empty and singleton copy-neutral sets:
statistics remain NaN and resulting calls are neutral; the implementation no longer
turns singleton standard deviation into zero or falls back to the original set after
outlier filtering. `-O` now derives `<output stem>.igv.seg` automatically, while
`--legacy-output` overrides the path; non-finite manifest statistics are JSON null.
Pinned GATK called/legacy-byte oracles pass on OpenMP and Serial. Full ModelSegments,
cloud and release-specific semantics remain explicit fallback and the registry stays
prototype.

The same caller now preserves the GATK 4.6.2.0 TSV boundary for non-finite
`MEAN_LOG2_COPY_RATIO` values (`NaN`, `Infinity`, and `-Infinity`). These values flow
through the Kokkos `2^log2`/IEEE-754 statistics path, receive neutral calls, and retain
Java's textual representation in both called and `.igv.seg` outputs. The pinned
`fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle` compares main and legacy
output bytes plus JSON-null degenerate statistics on OpenMP and Serial. This is a
bounded input/output slice; complete ModelSegments probability semantics remain
explicit fallback.

SelectVariants also has a pinned GATK 4.6.2.0 FILTER-exclusion boundary: the optional
Boolean `--exclude-filtered`/`--exclude-filtered-variants` forms (bare/true/1 versus
false/0), PASS/dot handling, multi-label ordering, OutputManifest choice, and filtered
record counts are oracle-tested on OpenMP and Serial. This does not promote the
SelectVariants prototype to full ID/pedigree/random-fraction/cloud/JEXL parity.

GenotypeGVCFs now has a pinned dense-reference writer oracle for
`--include-non-variant-sites --gatk-compatible-annotations`.  For monomorphic
reference records the native compatibility profile matches GATK 4.6.2.0's
MIN_DP-to-DP recovery, GQ-to-RGQ move, and PL/MIN_DP removal.  Aggregate and
`--stream-by-locus` outputs both match 101 Java data rows exactly on OpenMP and
Serial; the default native diagnostic profile intentionally keeps RCQ/RCP and
GQ/PL for telemetry and regression compatibility.  This is a writer/profile
boundary, not a claim of complete joint-genotyping algorithm or genome-scale
bit identity.
HaplotypeCaller now also follows GATK 4.6.2.0's repeatable ReadThreadingAssembler
`--kmer-size` boundary: the CLI default is normalized to `[10,25]`, reordered repeated
requests are sorted/deduplicated at the Host boundary, and the legacy scalar telemetry
field retains the last explicit value. The shared Kokkos graph orchestration tries the
requested list before cyclic/non-unique reference expansion. The pinned
`fastgatk-hc-kmer-list-gatk-oracle` compares default and reordered `[25,10]` Java/native
VCFs byte-for-byte on OpenMP and Serial and checks normalized telemetry; full graph/path
and complex-assembly equivalence remain explicit follow-up scope.

The shared HC/Mutect2 ActivityProfile also consumes explicit CIGAR insertion/deletion
anchors as Kokkos activity seeds. This preserves GATK AssemblyRegion ownership for reads
whose aligned bases are reference-matching but whose CIGAR carries an indel, avoiding the
narrow fallback that can disqualify all PairHMM reads. The pinned
`fastgatk-hc-cigar-indel-activity-gatk-oracle` recreates 8x `60M2D58M` plus 8x `120M`
and matches GATK 4.6.2.0's normalized `chr1:250 GCA>G` record and genotype payload on
OpenMP and Serial; complete ReadThreadingAssembler/graph parity remains explicit fallback.

GetPileupSummaries now supports GATK-compatible `-XL/--exclude-intervals` with
`-ip`/`-ixp` padding and applies exclusions after `-L` (including VCF-site
selectors). Separate-token Boolean forms for default-filter and dictionary
switches are accepted; the pinned GATK 4.6.2.0 oracle covers these boundaries
on OpenMP and Serial.

`ModelSegments --segments` now accepts a Picard interval-list partition with
GATK-compatible dictionary-order/non-overlap validation. Supplied intervals
skip kernel segmentation, and copy-ratio observations are assigned by their
midpoints; boundaries, point counts, and `.cr.seg` rows are verified against
pinned GATK 4.6.2.0 on OpenMP and Serial by
`verify_model_segments_input_segments_gatk_oracle.py`. The native posterior
now includes the deterministic Gibbs-conditional copy-ratio prepass; its
midpoint partition/report shape and an outlier-responsibility fixture are
verified against pinned GATK 4.6.2.0 by
`verify_model_segments_copy_ratio_conditionals_gatk_oracle.py`. Java
slice-sampler draw/raw-bit parity remains bounded prototype scope.

The same probabilistic path now retains ref/alt counts aggregated from the
allelic-count input and runs a Kokkos `RangePolicy` allele-fraction conditional
pass. It updates responsibility-weighted segment minor fractions together with
mean-bias, bias-variance and outlier-probability summaries, then uses those
conditionals to initialize and scale the MAF chain. The pinned oracle also
checks finite MAF posterior deciles and the GATK-shaped `.af.param` report;
Java Gamma-bias marginalization and credible-interval smoothing remain
explicit fallback boundaries.

ApplyVQSR now also treats `--ignore-all-filters` and `--exclude-filtered` as
GATK optional Booleans. Its pinned GATK 4.6.2.0 oracle checks bare, `true`,
and `false`: pre-filtered records are rescored only when requested, and new
VQSR-filtered records are omitted only when requested. This improves the
writer/score boundary without claiming complete VQSR model provenance parity.

ModelSegments combined allele-fraction segmentation

For combined copy-ratio and allelic-count inputs, the Kokkos segmenter uses
the first coordinate-ordered heterozygous site's oriented ALT fraction in
each copy-ratio interval, matching GATK 4.6.2.0. The signal is kept separate
from folded MAF modeling evidence; missing sites are imputed to 0.5. The
`fastgatk-model-segments-first-alt-fraction-gatk-oracle` fixture has two
oppositely oriented sites per interval (average ALT fraction 0.5, MAF 0.1)
and pins the Java/native `1-1600/1601-3200` boundary on OpenMP and Serial.
Full per-locus Java MCMC/report parity remains explicit prototype scope.

For HDF5 SVD PoN input, DenoiseReadCounts now enforces GATK 4.6.2.0's
SimpleInterval identity contract: the case interval list must exactly match
the PoN original interval list by contig/start/end and order before filtered
panel intervals are selected. A pinned Java oracle verifies valid acceptance
and subset, superset, and coordinate-mismatch rejection on both OpenMP and
Serial. This avoids changing the fractional-coverage denominator while
silently accepting a non-GATK case; complete CNV model/Spark/cloud semantics
remain explicit fallback.

VariantFiltration now supports GATK's optional-Boolean
`--invalidate-previous-filters`: it clears prior site FILTER labels while
retaining genotype FT history before applying new predicates. The pinned
GATK 4.6.2.0 oracle passes on both OpenMP and Serial; unsupported full JEXL
and annotation-engine behavior remains explicit fallback.

Its mask boundary also accepts GATK's `--filter-not-in-mask` reverse mode:
normal mode filters records overlapping an indexed VCF/BCF mask, reverse mode
filters records outside it, and both preserve existing site FILTER labels.
`--mask-description` controls the emitted mask FILTER header description, with
GATK's normal/reverse default descriptions when omitted. The pinned Java oracle
indexes the mask through `IndexFeatureFile` and compares both directions,
custom description, chained history and manifest metadata on OpenMP and
Serial. Full JEXL/annotation-engine, cloud mask and complex allele-specific
semantics remain explicit fallback.

For the GATK-compatible GenotypeGVCFs writer, the inherited
`--genotype-assignment-method`/`--gam` enum follows GATK 4.6.2.0's effective
`PREFER_PLS` behavior (the Java engine forces it internally); no extra GP/PG
fields or posterior-GQ changes are emitted. The diagnostic native profile
still exposes richer assignment modes for kernel tests. The pinned alias/row
oracle passes on OpenMP and Serial.

GenotypeGVCFs now mirrors GATK's final output-allele subset boundary: the
Kokkos AF posterior removes unsupported sibling ALTs at the default confidence
30 threshold (`--standard-min-confidence-threshold-for-calling`/
`--stand-call-conf`) and remaps PL/AD/GT and Number=A annotations before the
writer. The pinned GATK 4.6.2.0 two-sample multi-ALT oracle passes on OpenMP
and Serial for both `PREFER_PLS` and `USE_PLS_TO_ASSIGN` spellings. This does
not claim full cohort calibration or raw bit identity.

For multi-sample `--include-non-variant-sites`, INFO/DP is recomputed from
the merged sample-major FORMAT/DP values, matching GATK's cohort depth instead
of retaining the first source block's depth. A pinned GATK 4.6.2.0 oracle
verifies 101 rows in aggregate and `--stream-by-locus` modes on OpenMP and
Serial.

`CombineGVCFs` also preserves PL-less reference-confidence inputs at the
writer boundary. If no merged sample contains FORMAT/PL, the native output
omits FORMAT/PL rather than creating `.,.,.` likelihoods; pure
`REF/<NON_REF>` blocks retain `END` without an INFO/DP promotion from sample
FORMAT/DP. The pinned GATK 4.6.2.0 two-sample
`fastgatk-combine-gvcfs-plless-gatk-oracle` compares rows exactly on OpenMP
and Serial. PL-bearing merge/remap and explicit genotype materialization
remain covered by the existing CombineGVCFs paths; full joint posterior
semantics are still explicit fallback.

`CombineGVCFs` interval selection also matches a GATK 4.6.2.0
reference-confidence boundary: if an interval ends inside a `REF/<NON_REF>`
block, Java warns and drops the block; if the block `END` is contained, Java
emits the original full block rather than clipping it. Aggregate and
`--stream-merge` use the same predicate, covered by
`fastgatk-combine-gvcfs-interval-refblock-gatk-oracle`.

GenotypeGVCFs now honors `--max-alternate-alleles` (default 6) before cohort
AF calculation. The Kokkos likelihood-score kernel follows GATK's
`AlleleSubsettingUtils` best-PL genotype rule, keeps stable input-order ties,
and the Kokkos remap path updates PL/AD/GT and Number=A fields. A genotype
whose ALT was dropped becomes `./.` with missing GQ and a zero-based projected
PL row. The pinned GATK 4.6.2.0 max-ALT oracle is exact in aggregate and
`--stream-by-locus` modes for OpenMP and Serial.

The same GenotypeGVCFs path supports GATK's optional
`--annotate-with-num-discovered-alleles` Boolean. With the switch enabled,
INFO/NDA records the number of concrete ALTs discovered before max-ALT
subsetting (excluding `<NON_REF>`); reference-only blocks do not receive NDA.
The pinned seven-ALT GATK 4.6.2.0 oracle asserts `NDA=7`, exact compatible
INFO ordering and manifest telemetry in aggregate and stream-by-locus modes
on OpenMP and Serial.

For upstream gVCFs that already contain FORMAT/GP and FORMAT/PG, native now
matches GATK 4.6.2.0's posterior-field lifecycle. `GenotypeGVCFsEngine` forces
`PREFER_PLS`; when the symbolic `<NON_REF>` allele is removed,
`AlleleSubsettingUtils` drops those pre-subset Number=G posterior vectors, and
`--gp-qual` therefore does not rewrite this PL-based result. The
GATK-compatible native writer clears GP/PG only after cohort and annotation
consumers finish, preventing stale vectors from leaking into final records;
the diagnostic profile keeps them for explicit posterior-assignment tests.
`fastgatk-genotype-gvcf-gp-input-gatk-oracle` uses a two-sample GP/PG fixture and
requires Java/native row identity in aggregate and `--stream-by-locus` modes;
OpenMP and Serial CTest #71 pass. This is a bounded posterior FORMAT cleanup
slice, not a claim of complete release-specific joint-posterior calibration.

GenotypeGVCFs also applies GATK's spanning-deletion ownership boundary. In the
GATK-compatible profile, an orphan `*` ALT with no concrete deletion span is
kept through the Kokkos cohort AF/QUAL calculation, then removed at the final
output-allele subset; samples whose genotype used that ALT become `./.` with
GQ=0 and normalized PL. A concrete deletion record covering the locus keeps
`*` unchanged. The pinned
`fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` compares both cases
exactly in aggregate and `--stream-by-locus` modes on OpenMP and Serial, and
the manifest reports `orphan_spanning_deletion_loci`. This bounded rule does
not claim full release-specific joint-posterior calibration.

ApplyVQSR default threshold boundary

When `--tranches-file` is present without an explicit
`--truth-sensitivity-filter-level`, pinned GATK 4.6.2.0 uses the effective
`VQSLOD >= 0.0` cutoff and writes `LOW_VQSLOD`; it does not select the first
tranche. Native follows this rule. The independent
`fastgatk-apply-vqsr-default-cutoff-gatk-oracle` passes on OpenMP and Serial.

ApplyVQSR also supports GATK's `--sites-only-vcf-output` writer boundary. The
native path retains full FORMAT/sample state during score and FILTER
processing, then publishes an exact eight-column VCF after all annotations are
complete. The pinned `fastgatk-apply-vqsr-sites-only-gatk-oracle` compares the
Java 4.6.2.0 rows and manifest metadata on OpenMP and Serial.

Mutect2 also supports GATK 4.6.2.0 `--force-active[=true|false]`. The shared
Kokkos ActivityProfile path keeps profile coordinates unchanged and, when
enabled, materializes inactive segments as AssemblyRegion graph/PairHMM work;
stats and manifest record the effective switch. The pinned
`fastgatk-mutect2-force-active-gatk-oracle` checks 22 Java/native IGV rows,
active-state changes, and candidate-site stability on both OpenMP and Serial.

ApplyVQSR also accepts GATK's repeatable `-XL/--exclude-intervals` and
`-ixp/--interval-exclusion-padding`. Exclusions are applied after `-L` include
selection and before score/filter work; OutputManifest records excluded rows,
padding, and interval-file statistics. The pinned
`fastgatk-apply-vqsr-exclude-intervals-gatk-oracle` covers include+exclude,
exclusion-only, and padded traversal on OpenMP and Serial. Full VQSR
model/recalibration-table provenance remains explicit fallback.
Mutect2 `--normal-lod` is supported at the matched-normal emission boundary
(default 2.2). Native applies the threshold to normal ALT log10 odds after
likelihood materialization and records the value in stats/manifest. NLOD uses
GATK's fixed diploid `hom-ref - ref/alt-het` likelihood ratio and is not derived
from the variational somatic evidence used for TLOD. The pinned matrix oracle
checks finite and zero-probability ALT rows, while the command-level
`fastgatk-mutect2-normal-lod-gatk-oracle` checks Java/native emit-versus-
suppress behavior on a deterministic tumor/normal fixture with both OpenMP
and Serial backends.

HaplotypeCaller now accepts its standard
`--min-base-quality-score` spelling and `-mbq` alias in addition to the
earlier native compatibility spelling.  Both select the same Host-side
pre-assembly quality threshold and are written to output telemetry, without
changing the Kokkos ActivityProfile/graph/PairHMM interfaces.  The pinned
CIGAR-deletion oracle confirms that GATK 4.6.2.0 and native both remove Q40
evidence at `-mbq 60` on OpenMP and Serial.

The same HC/gVCF boundary now honors
`--reference-model-deletion-quality` rather than baking its default Q30 into
reference-confidence CIGAR-deletion evidence.  The Host passes the configured
quality into the existing Kokkos Ref-vs-Any and arbitrary-ploidy PL/GQ path;
the indel-informative cache remains GATK's fixed Q45 model.  A dual Q10/Q60
BP_RESOLUTION CIGAR-deletion oracle compares every Java/native row and proves
that the two requested qualities produce distinct reference-confidence output
on both OpenMP and Serial.

HaplotypeCaller `--alleles` now follows GATK's GenotypeGivenAlleles boundary
for covered concrete SNP, multi-ALT SNP, and anchored-indel feature records.
The C++ Host reads indexed VCF/BCF features, excludes filtered records unless
`--force-call-filtered-alleles` (or GATK's
`--genotype-filtered-alleles` alias) is enabled, and injects accepted alleles
into the local haplotype set. A per-locus feature mask crosses the existing
Kokkos ActivityProfile kernel, so the forced coordinate receives exactly the
source activity state without enabling global `--force-active`. The same
Kokkos AF path supplies the complementary monomorphic-site confidence used by
GATK when a forced ALT is emitted as `0/0`. The pinned GATK 4.6.2.0 oracle
compares complete VCF rows for a reference-only SNP, a three-ALT record, and
an anchored deletion; it also checks the GVCF concrete-ALT plus `<NON_REF>`
record, filtered-feature policy, and the no-read region where both
implementations correctly emit no row. All cases run on OpenMP and Serial.

ModelSegments default kernel dispatch

When `--segments` is omitted, GATK uses `KernelSegmenter`; native now selects
the Kokkos kernel by default. The native `--change-point-threshold` remains an
explicit legacy threshold extension. Copy-ratio-only linear segmentation uses
the exact rank-one feature of GATK's SVD linear kernel (`x -> x`) instead of an
anchor-scaled outer-product embedding. The pinned
`fastgatk-model-segments-default-kernel-gatk-oracle` checks step, gradient and
blip profiles against GATK 4.6.2.0 on OpenMP/Serial. Gaussian/multidimensional
projection, full Java SVD raw-bit identity and MCMC posterior identity remain
explicit fallback boundaries.
