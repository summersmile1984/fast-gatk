# fastgatk-hc-call

**GATK equivalent:** `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller`
  https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCaller.java
**Source:** `fastgatk-native/src/hc_call.cpp`
**Status:** shipped (contract-compatible, deterministic prototype)
**Oracle(s):** `scripts/verify_hc_*.py` (60+ scripts in
`fastgatk-native/scripts/`)

## 1. Purpose

Single-sample germline short-variant caller.  Assembles reads into a
de-Bruijn graph, finds candidate haplotypes, runs Smith-Waterman
realignment against the reference, runs a per-read-haplotype PairHMM
likelihood (GKL float32 / logless float64 variants), and emits VCF or
GVCF.

## 2. CLI parity

The full per-flag parity table is in
[`docs/cli-alignment.md`](../cli-alignment.md) — see the HaplotypeCaller
section (rows for `fastgatk-hc-call`).  Highlights:

- **Shipped (73 `accepted`)**: standard HC flags including
  `--kmer-size`, `--pruning-lod-threshold`, `--max-haplotype-paths`,
  `--flow-assembly-collapse-hmer-size`, `--native-pair-hmm-use-double-precision`,
  `--create-output-variant-index`, `--sites-only-vcf-output`,
  `--floor-blocks`, `--add-output-vcf-command-line`,
  `--disable-cap-base-qualities-to-map-quality`, etc.
- **Tier-1b gaps (sprint-1)**:
  - `StrandBiasBySample` (FORMAT): header declared
    (`hc_call.cpp:3941`), value never populated.  Fix
    in Tier-1b.
  - `InbreedingCoeff` (INFO): header declared
    (`hc_call.cpp:3087`), value not computed.  Fix in Tier-1b.
- **Catalog-ignored (83)**: known GATK flags that are
  silently consumed by `fastgatk::cli::consume(...)` without action
  (debug / verbosity / domain-specific).
- **Unsupported (18)**: Tier-2a and beyond (e.g. `--flow-mode`,
  `--dragen-mode`, `--apply-bqd`, `--transform-dragen-mapping-quality`,
  `--use-flow-aligner-for-stepwise-hc-filtering`, `--smith-waterman-implementation WFA`).
  See [`PARITY_AUDIT.md`](../../../work/mutect2-priority/PARITY_AUDIT.md).

## 3. Algorithm(s) used

- [`docs/algorithms/activity-profile.md`](../algorithms/activity-profile.md)
- [`docs/algorithms/kmer-graph.md`](../algorithms/kmer-graph.md)
- [`docs/algorithms/smith-waterman.md`](../algorithms/smith-waterman.md)
- [`docs/algorithms/reference-confidence.md`](../algorithms/reference-confidence.md) (Tier-1b focus)
- [`docs/algorithms/genotype.md`](../algorithms/genotype.md)
- [`docs/algorithms/filter-reads.md`](../algorithms/filter-reads.md)
- [`docs/algorithms/pairhmm.md`](../algorithms/pairhmm.md) (Kokkos SIMD bucketed)

## 4. Output contract

- **VCF** / **VCF.gz** (tabixed by default) with sample column =
  one selected @RG SM sample.  Standard GATK 4.6.2.0 INFO/FORMAT
  schema (`##INFO=<ID=DP,Number=1,…>`, `##FORMAT=<ID=GT:…>`,
  `##FORMAT=<ID=AD,Number=R,…>`, etc.).  Native OutputManifest is
  valid JSON with `status="prototype"` or
  `"contract-compatible"`.
- **GVCF** / **GVCF.gz** via `-ERC GVCF` or `-ERC BP_RESOLUTION`.
  Standard GQ-band blocks; `<NON_REF>` rows for joint-genotyping
  candidates.

## 5. Oracle(s) and expected parity

Listed in [`docs/verifications/contracts-hc.md`](../verifications/contracts-hc.md).
Status from latest `evidence/<date>-rerun/`:

- 60 / 61 HC verify scripts pass on HEAD (`evidence/2026-09-25-rerun/`)
- 1 fail: `verify_hc_graph_gatk_oracle.py` (pre-existing
  `--graph-output` flag mismatch — Tier-2a, deferred)

## 6. Known gaps

Full inventory: [`work/mutect2-priority/PARITY_AUDIT.md`](../../../work/mutect2-priority/PARITY_AUDIT.md)
section "Tier 2 — HC specific gaps".  Sprint-1 focus on Tier 1b (this
doc's "Tier-1b gaps" section above).

## 7. Build / install

- CMake target: `fastgatk-hc-call`
- Link deps: `libfastgatk-calling.a`, `libfastgatk-kernels.a`,
  `libfastgatk-core.a`, `libfastgatk-runtime.a`, `libfastgatk-hts.a`,
  `libkokkoscore.a`, `libkokkossimd.a`, HTSlib
- Kokkos backends tested: OpenMP, Threads, Serial (default OpenMP)
- GPU backends (CUDA, HIP, SYCL): code compiles, **CI gate not
  validated** (`work/mutect2-priority/RESEARCH.md` notes this)
- GPU-safety policy: `fastgatk/kernels/gpu_safety.hpp` —
  `fastgatk-kernels/src/pairhmm_kokkos.cpp` and `smith_waterman.cpp`
  are GPU-portable via `WithoutInitializing` + `Kokkos::resize`.
  No `<immintrin.h>`, no `_mm_*`, no third-party host-only SIMD libs.

## 8. CLI quick reference

```bash
# Single-sample on NA12878 chr17 69k-70k (the canonical pinned fixture)
fastgatk-hc-call \
    --input gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam \
    --reference gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta \
    --sample-name NA12878 \
    -L 17:69000-70000 \
    -O out.vcf.gz

# GVCF output (reference-confidence)
fastgatk-hc-call \
    --input sample.bam \
    --reference ref.fasta \
    --sample-name S1 \
    -ERC GVCF \
    -O out.g.vcf.gz
```

Full option surface: see `cli-alignment.md` HaplotypeCaller section.

## 9. Cross-references

- **Upstream feeders**:
  - `fastgatk-mark-duplicates` (read coordinate sort + duplicate flag)
  - `fastgatk-bqsr` + `fastgatk-apply-bqsr` (BQSR tables consumed)
- **Downstream consumers**:
  - `fastgatk-genotype-gvcf` (joint genotyping)
  - `fastgatk-variant-filtration` (VCF filtering)
  - `fastgatk-variant-recalibrator` (VQSR training)
  - `fastgatk-collect-f1r2-counts` (Mutect2 orientation prior)
- **Parallel callers**: `fastgatk-mutect2` (shares calling_pipeline.cpp)
