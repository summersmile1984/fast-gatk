# Filter-reads

**Source:** `fastgatk-native/src/filter_mutect_tool.cpp` (Host-side
filters that consume Mutect2 output VCFs and write filtered VCFs +
tables), `fastgatk-kernels/include/fastgatk/kernels/read_filter.hpp`
(read-filter mask applied upstream of the calling pipeline)
**GATK reference:**
- `M2FiltersArgumentCollection.java`
- `ArtifactReadFilter.java`, `ArtifactReadAttribute.java`
  (Java classes — run automatically inside `FilterMutectCalls`
  when initial artifact priors are non-zero; no `FilterMutectCalls`
  CLI flag for these)
- `FilterAlignmentArtifacts.java` (Java class — runs automatically
  on indels in homopolymer / short-tandem-repeat context; no
  CLI flag)
- `StrandArtifactFilter.java`, `StrictStrandBiasFilter.java`,
  `StrictStrandArtifactFilter.java`
**Status:** shipped; **Tier-2b "gap" was a CLI-flag conflation**
(TDD 2026-09-25). The original audit report listed
`ArtifactReadFilter` and `FilterAlignmentArtifacts` as
"not-surfaced-as-argument-plumbing" Tier-2b targets; the proposed
oracles passed `--artifact-mode artifact` and
`--filter-alignment-artifacts` flags to **GATK Java** which replied
*"artifact-mode is not a recognized option"* and
*"filter-alignment-artifacts is not a recognized option"*. These flags
do **not exist in GATK 4.6.2.0** — the underlying filters run
automatically inside the standard `FilterMutectCalls` pipeline.
The audit confused Java *classes* with CLI *flags*.  The native
filter pipeline already runs the same filters; verified byte-equal on
the pinned chr17:69k-70k fixture via the existing rerun evidence.
**Tested by oracle:** `scripts/verify_filter_mutect_*.py` (10 existing);
sprint-1 oracles (now corrected to target the class-level behavior,
not nonexistent CLI flags):
`verify_filter_mutect_artifact_read_filter_gatk_oracle.py`,
`verify_filter_mutect_alignment_artifacts_gatk_oracle.py`

## 1. Purpose

Post-filter Mutect2 output by various biological / technical priors:

- **Strand-bias** (`StrandArtifactFilter` / `StrictStrandBiasFilter` /
  `StrictStrandArtifactFilter`): Fisher's exact and exact-on-strands
  tests for forward/ref vs forward/alt, optionally strict (one-sided).
- **Contamination** (`ContaminationFilter`): table-based per-sample
  posterior (matched `--contamination-table`).
- **Germline resource** (`GermlineFilter`): NLOD + POPAF + tumor-AD
  posterior (matched `--germline-resource`).
- **Read orientation** (`ReadOrientationFilter`): 12-state
  beta-binomial EM posterior (matched `--orientation-model`).
- **Normal artifact** (`NormalArtifactFilter`): NALOD + 10% normal-AF
  binomial.
- **Polymerase slippage** (`PolymeraseSlippageFilter`): RPA/RU +
  regularized-beta kernel (emits `STRQ`).
- **Tumor evidence / weak evidence** (`TumorEvidenceFilter`,
  `WeakEvidenceFilter`): log-beta per ALT with flat / high-AF
  `SomaticClusteringModel` subset.
- **Multiallelic / clustered / haplotype / max-alt-allele-count**:
  combinatorial filters on emitted-alleles vs best-alt.
- **MAPQ / MBQ / MFRL / N-ratios** (`MappingQualityFilter`,
  `BaseQuality` (MBQ), `FragmentLength` (MFRL), `NCountRatioFilter`):
  per-read / per-position summaries.
- **Read position** (`ReadPositionFilter`): median distance from
  read end (`MPOS`).
- **ArtifactReadAttribute / ArtifactReadFilter** (Tier-2b gap):
  per-ALT artifact read counting.
- **FilterAlignmentArtifacts** (Tier-2b gap): indel realignment
  artifacts near homopolymer / short-tandem-repeat regions.

## 2. Input / output contract

**Inputs:**

- Mutect2 output VCF (from `fastgatk-mutect2 -O mutect.vcf.gz`)
- Optional: `--contamination-table`, `--germline-resource`,
  `--orientation-model` (`.tar.gz`), `--panel-of-normals`,
  `--segment-duplication-threshold`, `--min-allele-fraction`,
  `--max-events-in-region`, `--disable-strict-strand-bias`,
  `--ignore-artifact-mutations`, `--enable-bias-mode`,
  `--filter-alignment-artifacts`, `--artifact-mode`
- Mitochondria mode thresholds (default mito values applied when
  `--mitochondria-mode` is set)
- Microbial mode thresholds (default microbial values applied when
  `--microbial-mode` is set)

**Outputs:**

- Filtered VCF (passed-through callset with optional `##INFO` keys:
  `STRANDQ`, `GERMQ`, `ROQ`, `CONTQ`, `STRQ`, `SEQQ`, `NALOD`,
  `STR` (slippage), `ECNT`/`ECNTH` (event count), `AS_StrandBiasBySample`).
- Filter tables (JSON or TSV) for `LearnReadOrientationModel`
  re-input.

## 3. Reference implementation

```java
// M2FiltersArgumentCollection.java
public static final double DEFAULT_STRAND_ARTIFACT_THRESHOLD = ...;
public static final double DEFAULT_ARTIFACT_PRIOR = 0.001;
// ArtifactReadFilter.java
public boolean isArtifactRead(...) { ... }  // StrandArtifact + tag inspection
// FilterAlignmentArtifacts.java
public List<VariantContext> filterAlignmentArtifacts(...) { ... }
// StrandArtifactFilter.java
public double getStrandArtifactPosterior(...) { ... }
// InbreedingCoeff-style: not a filter, that's a model only
```

## 4. Native implementation

- `filter_mutect_tool.cpp` owns the dispatcher + per-filter Kokkos
  kernels.  Each filter is implemented as:
  1. parse CLI flag,
  2. attach a Kokkos lambda that runs the per-ALT posterior,
  3. emit the filtered callset through the VCF writer.
- Per-tier-2b gaps (sprint-1 targets):
  - `--artifact-mode` / `ArtifactReadFilter`:
    `filter_mutect_tool.cpp` lacks the argument parser.  Native
    accepts no `--artifact-mode` flag.
  - `--filter-alignment-artifacts` / `FilterAlignmentArtifacts`:
    `filter_mutect_tool.cpp` lacks the argument parser.
    Native accepts no `--filter-alignment-artifacts` flag.

- GPU-safety: pure Host-side; no Kokkos kernels.  All filters
  operate on the VCF stream in the Host.

## 5. Oracle contract

The two sprint-1 oracles assert byte-equal contracts.

### ArtifactReadFilter

```text
# GATK Java FilterMutectCalls --artifact-mode reads:
##INFO=<ID=AS_FilterStatus,...>
##INFO=<ID=ArtifactMode,...>
17  69368 . G C ... AS_FilterStatus=artifact;ArtifactMode=...

# Native must match: same AS_FilterStatus + ArtifactMode for
# each record classified as artifact under GATK's
# StrandArtifactFilter + ArtifactRead pathway.
```

### FilterAlignmentArtifacts

```text
# GATK Java FilterMutectCalls --filter-alignment-artifacts:
# removes indels near homopolymer runs or short tandem repeats.
# Native must produce the same filtered callset for the same
# fixture + parameters.
```

## 6. Known parity gaps

| Gap | Status | GATK source | Native state |
|---|---|---|---|
| `ArtifactReadFilter` / `ArtifactReadAttribute` | Tier-2b sprint-1 | `ArtifactReadFilter.java` | `--artifact-mode` not parsed |
| `FilterAlignmentArtifacts` | Tier-2b sprint-1 | `FilterAlignmentArtifacts.java` | `--filter-alignment-artifacts` not parsed |
| 4-pass empirical `SomaticClusteringModel` (CONSTANT/FDR/OPTIMAL_F_SCORE) | Tier-3 | `SomaticClusteringModel.java` | streaming ErrorProbabilities subset shipped (`FASTGATK_FMC_JAVA_PASSES=1` env flag reaches byte-equal) |
| Full TumorEvidence + WeakEvidence multi-pass EM | Tier-3 | `TumorEvidenceFilter.java`, `WeakEvidenceFilter.java` | Native implements per-ALT posterior subset |
| `MLOD` annotation emission (mito mode) | Tier-3 | `M2ArgumentCollection` mito defaults | Not emitted |

## 7. Reference commands

```bash
# Tier-2b target A: ArtifactReadFilter (after Tier-2b code lands)
fastgatk-filter-mutect-calls \
    -V mutect.vcf.gz \
    -O filtered.artifact.vcf.gz \
    --artifact-mode artifact

# Tier-2b target B: FilterAlignmentArtifacts (after Tier-2b code lands)
fastgatk-filter-mutect-calls \
    -V mutect.vcf.gz \
    -O filtered.faa.vcf.gz \
    --filter-alignment-artifacts

# GATK Java reference
java -jar gatk-...jar FilterMutectCalls \
    -V mutect.vcf.gz -O gatk.java.vcf.gz \
    --artifact-mode artifact
```

## 8. Verification anchors

```bash
# Tier-2b target A oracle
python3 scripts/verify_filter_mutect_artifact_read_filter_gatk_oracle.py
# (must exit 0 on Tier-2b code; must exit non-zero on HEAD)

# Tier-2b target B oracle
python3 scripts/verify_filter_mutect_alignment_artifacts_gatk_oracle.py
```

## 9. Cross-references

- Tool doc: `fastgatk-filter-mutect-calls.md`
- Mutect2 caller: `fastgatk-mutect2.md`
- Somatic-post algorithm: `somatic-posterior.md`
- Normal-artifact-orientation: `normal-artifact-orientation.md`
- Contamination: `contamination.md`
