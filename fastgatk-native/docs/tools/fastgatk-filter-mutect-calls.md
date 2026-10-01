# fastgatk-filter-mutect-calls

**GATK equivalent:** `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls`
  https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/filtering/FilterMutectCalls.java
**Source:** `fastgatk-native/src/filter_mutect_tool.cpp`
**Status:** shipped (contract-compatible; Tier-2b focus for sprint-1)
**Oracle(s):** `scripts/verify_filter_mutect_*.py`

## 1. Purpose

Post-filter Mutect2 output (VCF) by various biological priors:
strand bias (Fisher / strict strand artifact), contamination
(`--contamination-table`), germline (`--germline-resource`), read
orientation (`--orientation-model`), polymerase slippage, weak
tumor evidence, etc.  Emits a filtered VCF + matched tumor and
normal-sample tables.

## 2. CLI parity

The full per-flag parity table is in
[`docs/cli-alignment.md`](../cli-alignment.md) — see the
FilterMutectCalls section (rows for `fastgatk-filter-mutect-calls`).
Highlights:

- **Shipped (42 `accepted`)**: TLOD / NLOD thresholds, MLOD
  (mitochondria mode), strand bias filter (`--strance-bias-threshold` /
  `--strict-strand-bias`), orientation-model filter (`--orientation-model`),
  contamination filter, germline resource filter, weak-evidence
  filter, multiallelic / clustered / haplotype / max-alt-allele-count
  filters, MAPQ / MBQ / MFRL / N-ratios filters, slippage filter,
  ReadPosition filter, normal-artifact filter.
- **Tier-2b gaps (sprint-1)**:
  - `--artifact-mode` / `ArtifactReadFilter` / `ArtifactReadAttribute`
    — entire filter pathway missing.  Fix in Tier-2b.
  - `--filter-alignment-artifacts` /
    `FilterAlignmentArtifacts` — entire filter missing.  Fix in
    Tier-2b.
- **Catalog-ignored (491)**: known GATK flags that are silently
  consumed by `fastgatk::cli::consume(...)` without action.
- **Unsupported (212)**: Tier-3 and beyond.

## 3. Algorithm(s) used

- [`docs/algorithms/filter-reads.md`](../algorithms/filter-reads.md) — Tier-2b
  focus (ArtifactReadFilter, FilterAlignmentArtifacts)
- [`docs/algorithms/somatic-posterior.md`](../algorithms/somatic-posterior.md) —
  Beta-binomial / Dirichlet posterior (PSOMATIC / PGERMLINE / PARTIFACT / OBP)
- [`docs/algorithms/pileup-summary.md`](../algorithms/pileup-summary.md) —
  for contamination
- [`docs/algorithms/normal-artifact-orientation.md`](../algorithms/normal-artifact-orientation.md) —
  12-state beta-binomial EM for orientation-model filter
- [`docs/algorithms/slope-genotype.md`](../algorithms/slope-genotype.md)
  (TODO if applicable) — strand bias computation
- [`docs/algorithms/somatic-likelihood.md`](../algorithms/somatic-likelihood.md) —
  base evidence

## 4. Output contract

- **VCF** filtered from input, with new annotation columns
  (`##INFO=<ID=STRANDQ,…>`, `##INFO=<ID=GERMQ,…>`,
  `##INFO=<ID=ROQ,…>`, `##INFO=<ID=STRQ,…>`, `##INFO=<ID=CONTQ,…>`,
  `##INFO=<ID=NALOD,…>`, `##INFO=<ID=SEQQ,…>`).
- `--contamination-table` is optional; default uses population AF.

## 5. Oracle(s) and expected parity

Listed in [`docs/verifications/contracts-mutect.md`](../verifications/contracts-mutect.md).
Status from latest `evidence/<date>-rerun/`:

- 9 / 10 FilterMutectCalls verify scripts pass on HEAD
  (`evidence/2026-09-25-rerun/`)
- Tier-2b sprint-1 oracle scripts will be created

## 6. Known gaps

Full inventory: [`work/mutect2-priority/PARITY_AUDIT.md`](../../../work/mutect2-priority/PARITY_AUDIT.md)
section "Tier 3 — Mutect2 specific gaps".  Sprint-1 Tier-2b audit
closed-no-op per TDD (see `docs/algorithms/filter-reads.md` §6).

The `FASTGATK_FMC_JAVA_PASSES=1` env flag reaches byte-equal learning
for the joint contamination oracle
(`evidence/2026-09-25-rerun/`); the default INFO-learning path
matches GATK semantically but is not byte-equal.

## 7. Build / install

- CMake target: `fastgatk-filter-mutect-calls`
- Link deps: `libfastgatk-calling.a`, `libfastgatk-kernels.a`,
  `libfastgatk-core.a`, HTSlib
- Kokkos backends tested: OpenMP, Threads, Serial
- GPU-safety policy: `fastgatk/kernels/gpu_safety.hpp` enforced
  across all kernel files

## 8. CLI quick reference

```bash
# Default Tier-1 oracle (contamination joint path, byte-equal with
# FASTGATK_FMC_JAVA_PASSES=1; semantic-equivalent otherwise)
fastgatk-filter-mutect-calls \
    -V tumor.vcf.gz \
    -O filtered.vcf.gz \
    --contamination-table contamination.table \
    --segment-duplication-threshold 1.0 \
    --min-allele-fraction 0.0001

# Mito mode (default thresholds from M2FiltersArgumentCollection)
fastgatk-filter-mutect-calls \
    -V tumor.vcf.gz \
    -O filtered.mito.vcf.gz \
    --mitochondria-mode
```

Full option surface: see `cli-alignment.md` FilterMutectCalls section.

## 9. Cross-references

- **Upstream feeders**:
  - `fastgatk-mutect2` (VCF output)
  - `fastgatk-calculate-contamination` (`--contamination-table`)
  - `fastgatk-learn-read-orientation-model`
    (`--orientation-model`)
  - `fastgatk-get-pileup-summaries` (pileup-derived contamination)
- **Downstream consumers**:
  - downstream annotations pipelines (`Mutect2Plot` / `Mutect2FilterPanel` —
    not in scope)
  - `fastgatk-variant-filtration` (VCF post-filter)
- **Parallel callers**: `fastgatk-filter-mutect-calls` shares the
  somatic-posterior and StrandArtifact filters with `Mutect2Engine`'s
  embedded clustering (GATK behavior); the standalone
  FilterMutectCalls tool is a strict superset.
