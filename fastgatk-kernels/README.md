# fastgatk-kernels

Reusable CMake library target for the first numerical kernels:

- fastgatk::pairhmm::compute_kokkos: the existing Kokkos 5.2.0 PairHMM API,
  using `Kokkos::Experimental::simd<double>` for scalar/AVX2/AVX-512 builds;
- `compute_kokkos_bucketed(..., PairHmmPrecision::Float32)` is an explicit
  approximate throughput mode using the same Kokkos `RangePolicy`, views and
  `Kokkos::Experimental::simd<float>` API.  It uses scale-invariant 2^120
  initialization to keep float DP finite, reports `precision="float32"`, and
  is never treated as a GATK raw-bit-compatible result.  The default remains
  `PairHmmPrecision::Float64` (`precision="double"`);
- the legacy `compute_scalar`/`compute_simd` compatibility entry points are
  implemented by `src/pairhmm_compat.cpp` as batch-level calls into the same
  Kokkos PairHMM path.  The production target contains no `_mm*` intrinsics;
  the standalone `pairhmm-demo/Makefile` remains a historical raw-SIMD oracle
  only.
- fastgatk::pairhmm::PersistentBucketPlan: a long-lived length/shape bucket
  plan that retains Kokkos Views and the three-row DP workspace across repeated
  region calls, while restoring request order and reporting cache hits;
- fastgatk::pairhmm::compute_kokkos_full_matrix: a ragged full
  read-by-haplotype matrix API. It creates no padded sequences, uses the same
  persistent length buckets, and returns row-major likelihoods/scaled sums;
- fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos and
  fastgatk::pairhmm::reduce_read_allele_best_kokkos: shared Kokkos
  max-marginalization plus GATK-shaped REF-first/ALT-order BestAllele owner
  selection. Marginalization now accepts an explicit concrete allele count,
  so the reusable kernel can materialize REF plus any number of ALTs without
  silently truncating a multi-ALT row; the two-argument overload remains the
  biallelic compatibility boundary. HC and Mutect2 use the owner result to
  count one informative read once at a multi-ALT locus; the API reports
  prepare/execute/backend telemetry and is covered by API smoke, kernel
  benchmark, and the chr17 HC oracle. Full cross-AssemblyRegion posterior
  remains a separate gate;
- fastgatk::pairhmm::compute_kokkos_flow: an independent four-step flow-space
  recurrence with `(read_flow_len,hap_flow_len)` buckets. Host
  `fastgatk::io::decode_flow_read` supplies GATK `tp`/`t0`-calibrated
  `[read_flow][256]` probabilities; no flow tag or string parsing occurs in
  the Kokkos kernel. GPU wavefront and release raw-bit corpus remain explicit
  follow-up gates;
- fastgatk::kernels::smith_waterman_score_reference and
  fastgatk::kernels::smith_waterman_align_reference: deterministic affine-gap
  score plus Host traceback/CIGAR API.  They implement GATK's
  `SOFTCLIP`, `INDEL`, `LEADING_INDEL` and `IGNORE` endpoint/overhang modes,
  Java tie priority and an explicit `alignment_offset`.
- fastgatk::kernels::smith_waterman_score_kokkos: variable-length batch score
  API using Kokkos Views and the selected execution space.  Traceback remains
  Host-owned because CIGAR strings are variable-length output. Rectangular
  batches additionally use `Kokkos::Experimental::simd<int>` across
  independent alignments; ragged batches retain the same Kokkos scalar
  recurrence. The API reports `simd_width` and `simd_groups`, and the mixed
  benchmark reports separate heterogeneous versus SIMD-friendly timings.
- fastgatk::kernels::build_kmer_graph_kokkos: deterministic k-mer encoding and
  stable node/edge construction with Kokkos occurrence/degree counting,
  optional reference-window scaffolding with `(tid,start,end)` path metadata,
  reference-connected component classification, bounded deterministic
  haplotype-path traversal/sequence materialization and min-count/dangling-
  branch pruning telemetry, including the GATK-shaped
  `min_dangling_matching_bases` exact-match gate for SW recovery. Each materialized path also reports its half-open
  reference end, read-supported alternate-edge flag, support count, and the
  read span of its first/last alternate edge. HC consumes coordinate-compatible
  paths as additional PairHMM haplotypes and can conservatively project a
  stable insertion from a path with exact terminal/internal soft-clip support.
  The graph also exposes SeqGraph-style non-branching-chain compression
  counters (`seqgraph_nodes/edges`) and applies a real Host sequence-level
  SeqGraph rewrite over materialized paths: shared-prefix/suffix diamonds,
  common-suffix tails (default 10-base GATK threshold) and linear-chain zips
  are split/merged and the equivalent path language is reconstructed. Rewrite
  counts (`seqgraph_*_merges`, `seqgraph_suffix_*`) are returned alongside the
  topology. Read-edge multiplicity still supplies a support-descending, stable
  outgoing-edge tie-break. GATK linear-chain, diamond/tail, SNP-rethreading
  fixtures are covered by the API smoke; the API also exposes per-sequence
  `non_unique_kmers` matching ReadThreadingGraph's repeat set and reference
  repeat rejection telemetry.
  `num_pruning_samples` provides the MultiSampleEdge sample gate; optional
  adaptive-pruning controls expose bounded error-rate/LOD/seeding-LOD and
  max-unpruned-variant policies. `linked_de_bruijn_graph` keeps the raw
  junction topology and, unless `disable_artificial_haplotype_recovery` is
  set, synthesizes uncovered pivotal-edge haplotypes from the best completed
  path. `disable_seqgraph_simplification` remains the explicit no-zip switch.
  Optional BAM flags reverse-complement bit-0x10 reads into reference
  orientation before graph construction.
- `fastgatk::kernels::seqgraph_constant_error_log_likelihood_ratio`: the
  scalar Mutect2/AdaptiveChainPruner constant-error likelihood helper. Its
  Commons Math/GATK digamma cutoff (49), recursive unwind order and
  `binomialCoefficientLog` exact/double/log-sum branches are covered by the
  GATK 4.6.2.0 numerical oracle in API smoke; graph traversal remains
  Host-owned while all backends use the same source and contract.
- fastgatk::kernels::correct_read_errors_kokkos and
  fastgatk::kernels::correct_reads_by_pileup_kokkos: GATK-shaped assembly-only
  correction primitives. The Nearby path builds a sparse-to-nearest-solid
  exact-k-mer map (up to two Hamming mismatches) and applies strict overlapping
  consensus; the pileup path reproduces Mutect2's flat-Beta log-odds threshold,
  CIGAR indel/soft-clip exclusion and three-mismatch guard. Both keep irregular
  Host maps/traceback separate from the final Kokkos edit mask, so Serial,
  OpenMP and accelerator builds share one deterministic API.
- fastgatk::kernels::compute_activity_profile_kokkos: Kokkos pileup activity
  signal plus Host halo-expanded, deterministically merged AssemblyRegions. A
  reference-projected batch uses the reference base to compute the
  non-reference fraction; an empty reference vector keeps the explicit
  no-reference diversity fallback and is reported in `reference_aware`.
- fastgatk::kernels::filter_reads_kokkos: fixed-order flag/MAPQ/read-length
  filtering with deterministic per-locus downsampling and ordinal-stable mask.
  Optional GATK sentinel predicates reject MAPQ=255 (unavailable) and MAPQ=0.
  The same Kokkos mask also implements GATK's
  `NonChimericOriginalAlignmentReadFilter` by comparing the OA contig prefix
  with XM; missing either tag is accepted exactly as in Java. OA/XM tag
  presence is carried separately from flat payload spans, so a present-but-
  empty Z tag is not confused with an absent tag.
- fastgatk::kernels::count_bqsr_quality_kokkos: deterministic integer quality/
  mismatch histogram over a flat observation batch using Kokkos `TeamPolicy`
  team-local 188-bin tables followed by a Host league-order merge; CIGAR,
  reference, known-sites and covariate-key semantics remain Host-owned.
- fastgatk::kernels::count_bqsr_covariates_kokkos: deterministic integer
  reduction over Host-assigned dense covariate IDs using a bounded Kokkos
  `TeamPolicy`-local histogram and deterministic Host league-order merge;
  batches whose dense key space would exceed the memory cap safely fall back
  to a `RangePolicy` atomic histogram. String/context/cycle key construction
  remains Host-owned, while the per-batch count and mismatch accumulation no
  longer mutates a Host map. The result reports the selected policy,
  histogram workspace bytes and whether the team-local path was used.
- fastgatk::kernels::apply_bqsr_quality_kokkos: deterministic quality
  `RangePolicy` transform with GATK `QualityUtils.boundQual` [1,93] clamping;
  preserve-low-Q, report hierarchy and optional quantization remain Host-owned.
- fastgatk::kernels::derive_genotype_gt_gq_kokkos (with the diploid wrapper
  derive_diploid_gt_gq_kokkos): deterministic VCF genotype-combination
  minimum/runner-up selection for fixed ploidy, including non-diploid
  combinations; HTSlib sentinel encoding and VCF writing remain Host-owned.
- fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos: explicit
  posterior-probability genotype assignment using a Host-built log10 cohort
  prior vector and the same Kokkos PL matrix.  GenotypeGVCFs keeps
  `PREFER_PLS` as its default and invokes this path only for
  `--genotype-assignment-method USE_POSTERIOR_PROBABILITIES`.
- fastgatk::kernels::count_alleles_kokkos: deterministic fixed-ploidy,
  sample-major allele-count/AN reduction used by GenotypeGVCFs AC/AN/AF
  annotation; missing/vector-end sentinel decoding remains Host-owned.
- fastgatk::kernels::calculate_cross_sample_reference_confidence_kokkos:
  shared-genotype-prior, sample-major hom-reference posterior and deterministic
  joint `RCQ`/`RCP` primitive used by GenotypeGVCFs; VCF INFO/header ownership
  remains Host-owned.

The library depends on parent-provided Kokkos::kokkos and fastgatk-core targets.
The PairHMM implementation now lives in `fastgatk-kernels/src/pairhmm_kokkos.cpp`;
both `pairhmm-demo` and `fastgatk-native` link this production source without
creating a second Kokkos configuration. The demo is therefore only a client,
oracle and benchmark harness, not a source dependency of the production target.
The transition-table boundary follows GATK `PairHMMModel`, including its
`log1p` low-quality clamp (quality-zero inputs no longer create negative
no-event probabilities). The source-tree build uses the pinned Java-generated
`pairhmm-demo/results/gatk-tables.hex` by default, so strict raw-bit behavior
does not depend on a manually exported environment variable. Packaged builds
install the same file under `share/fastgatk/gatk-tables.hex` and probe that
location relative to the executable, so `cmake --install --prefix <prefix>`
remains strict without an environment variable. A package can stage a
different file with `FASTGATK_PAIRHMM_TABLES_FILE`, or override it at runtime
with `FAST_GATK_PAIRHMM_TABLES`; if neither artifact is available, the
portable math path remains an explicit non-bit-identical mode.

The somatic likelihood result keeps GATK's two normal/tumor statistics
separate: `tlod` is Dirichlet variational all-alleles evidence, while
`normal_log10_odds` is `hom-ref - ref/alt-het` from the fixed 1/2 mixture used
by `SomaticGenotypingEngine.diploidAltLogOdds`. The pinned
`fastgatk-somatic-normal-lod-gatk-oracle` covers finite and `-Infinity` ALT
likelihoods on both configured execution spaces.

The old pairhmm-demo headers remain compatibility forwarders, so the existing
demo and Makefile continue to work. Enable the API smoke test from a parent build
with -DFASTGATK_KERNELS_BUILD_TESTS=ON, then run:

~~~text
ctest --test-dir <build> --output-on-failure -R fastgatk-kernels-api-smoke
~~~

The same test-enabled build also provides `fastgatk-kernels-benchmark`, a
512-request mixed-length k-mer graph/reference traversal + activity + SW +
PairHMM regular/flow + genotype PL/allele-count workload, plus a 64K-observation
BQSR TeamPolicy count and RangePolicy transform. It emits benchmark schema version 7: warmup/repeat counts,
backend/SIMD/concurrency/compiler metadata, deterministic checksums, and
p50/p95 plus raw samples for SW, regular/flow PairHMM, the explicit float32
PairHMM path, triangular-PL genotype, allele-count,
cohort-EM, posterior genotype-assignment, cross-sample reference-confidence,
somatic-likelihood, somatic-posterior and BQSR count/apply prepare/kernel phases;
BQSR execution policy is recorded explicitly. The
persistent PairHMM plan reports cached shapes and cache hits so allocation
amortization is visible rather than hidden in a single hot-kernel number:

~~~text
FASTGATK_BENCH_WARMUP=1 FASTGATK_BENCH_REPEATS=10 \
fastgatk-native/build/fastgatk-kernels/fastgatk-kernels-benchmark
~~~

The backend matrix helper runs the API smoke and the same benchmark on the
OpenMP and Serial artifacts (when both are present), requiring identical
deterministic checksums while reporting SIMD width differences:

~~~text
python3 fastgatk-native/scripts/verify_kokkos_backend_matrix.py
~~~

FASTGATK_BENCH_WARMUP defaults to 1 and FASTGATK_BENCH_REPEATS defaults to 5;
both are bounded and must be positive. The benchmark is a correctness and
regression harness, not a GATK speed claim: Java/GKL comparisons must use the
same input digest, affinity, thread count, and prepare/IO accounting.

On this Ryzen 9 7945HX, the current native build reports PairHMM
`simd_width=8` (double lanes) and SW `sw_simd_width=16` (int lanes); a portable
build without a host SIMD ABI reports width one. ISA selection is therefore a
Kokkos build/deployment choice, not a second kernel implementation, and the
same deterministic checksums are required in both cases.

This target is a kernel/API extraction, not a claim that PairHMM or
Smith-Waterman is already a complete HaplotypeCaller implementation.  The
Kokkos SW score path, the Java-GATK overhang/tie corpus (including the
original-default, indel, substring and long homopolymer vectors),
reference-connected graph classification and bounded path traversal are
covered by the API smoke test.  The native caller currently uses a bounded,
fail-closed graph projection: only length-changing `I/D` alignments with an
exact soft-clip support floor are promoted to graph-derived candidates; SNPs
remain owned by the pileup path.  Full haplotype enumeration, mixed/complex
CIGAR projection and GATK's complete pruning semantics remain pending.
