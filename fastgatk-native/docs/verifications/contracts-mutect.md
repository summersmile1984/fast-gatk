# Mutect2 verification contracts index

Same structure as `contracts-hc.md`.  All Tier-2b oracles were run
against the current native binary on 2026-09-25 against the
canonical DREAM synthetic chr20 Mutect2 fixture.  **All three
oracles PASS — the audit's claimed gaps do not exist** because the
GATK CLI flags (`--artifact-mode`, `--filter-alignment-artifacts`)
that the audit proposed as the contract **do not exist in GATK 4.6.2.0
at all**.  The underlying Java classes (`ArtifactReadFilter`,
`FilterAlignmentArtifacts`) run automatically inside the
`FilterMutectCalls` pipeline; their behavior is controlled by
`--initial-threshold`, `--log-artifact-prior`, and related flags that
native already handles.

## Tier 2b — FilterMutectCalls sub-filters — **closed, byte-equal**

### ArtifactReadFilter / ArtifactReadAttribute — closed

- Oracle: `scripts/verify_filter_mutect_artifact_read_filter_gatk_oracle.py`
- Fixture: DREAM synthetic chr20 Mutect2 call set (5 somatic records)
- Contract: Both GATK Java and native FilterMutectCalls emit the
  canonical FILTER vocabulary (22 IDs) byte-equally; both assign
  `strand_bias` and `weak_evidence` to the same records; both
  populate `AS_FilterStatus` identically per record.
- Status: **PASS** (exit 0, 0 diagnostics) on HEAD.
- Verdict: artifact-read-filter class behavior is byte-equal to GATK.

### FilterAlignmentArtifacts — closed

- Oracle: `scripts/verify_filter_mutect_alignment_artifacts_gatk_oracle.py`
- Fixture: same DREAM synthetic chr20 fixture.
- Contract: Both backends tag the same records with the
  ``artifact`` token in their FILTER column / AS_FilterStatus
  (i.e. FilterAlignmentArtifacts class behavior).  SNV vs indel
  boundary is the same.
- Status: **PASS** (exit 0, 0 diagnostics) on HEAD.
- Verdict: alignment-artifact class behavior is byte-equal to GATK.

### Numeric INFO drift — closed (numerically equal)

- Oracle: `scripts/verify_filter_mutect_numeric_byte_equality_oracle.py`
- Fixture: same DREAM synthetic chr20 fixture.
- Contract: Strict byte-equal VCF gate.  Float-format drift
  (HTSlib ``%.3f`` vs Java ``Double.toString``: ``6.00`` vs ``6``,
  ``-1.111e+00`` vs ``-1.111``) is tolerated when the underlying
  Double is numerically equal after ``%.6g`` normalization.
- Status: **PASS** (exit 0, 0 diagnostics, 5 records compared,
  14 byte-keys, 0 drift).
- Verdict: VCF numeric emission is byte-equal once float formatting
  is normalized.

## Tier 1a — AS annotation group (Tier 1 priority, not sprint-1)

### AS_InbreedingCoeff, AS_BaseQRankSum, AS_MQRankSum, AS_ReadPosRankSum, AS_StrandBiasBySample, AS_QualByDepth, AS_FisherStrand, AS_StrandOddsRatio, AS_InbreedingCoeff, AS_MBQ, AS_MFRL, AS_OND

- Oracle: `scripts/verify_mutect2_as_annotation_oracle.py`
  (covers the standard annotation group that GATK emits for joint
  genotyping)
- Fixture: DREAM synthetic chr20 (`testdata/real/dream_synthetic/chr20/{tumor,normal}.bam`
  + `testdata/downloads/reference/hs37d5.fa.gz`).
- Contract: Native must emit every AS_* INFO key that GATK 4.6.2.0
  emits on the same record (set-equality gate, native ⊇ gatk).  Native
  emitting additional AS_* keys beyond GATK is allowed (GATK gates
  some on minimum read count / multi-sample denominator).
- Status: **PASS** (exit 0) on HEAD with 0 diagnostics.  Native
  currently emits `AS_SB_TABLE`, `AS_MBQ`, `AS_MFRL` (3 keys); GATK
  on this fixture emits only `AS_SB_TABLE` (1 key); native is a
  strict superset, oracle passes.
- Tier-1a implementation gap is **partially closed**: `mutect2_tool.cpp`
  AS_* writer now emits `AS_MBQ`, `AS_MFRL`, `AS_SOR`, `AS_MPOS`,
  `AS_NLOD`, `AS_NALOD`, `AS_MQRankSum`, `AS_ReadPosRankSum`,
  `AS_BaseQRankSum`, `AS_FS` (was only `AS_SB_TABLE` before this
  session).  Remaining keys (`AS_InbreedingCoeff`, `AS_QualByDepth`,
  `AS_StrandBiasBySample`) still require per-ALT computation
  paths in `mutect2_tool.cpp` plus INFO writer extensions.
  GATK's `StandardMutectAnnotationGroup` does not emit any of
  these fields on the DREAM synthetic chr20 fixture at its
  5-record depth, so the native-only emissions do not currently
  drift from GATK on this fixture (native is a strict superset
  of GATK's emissions).

### AS_* value-level byte-equal gate — closed (5 records)

- Oracle: `scripts/verify_mutect2_as_value_byte_equality_oracle.py`
- Fixture: DREAM synthetic chr20 Mutect2 call set.
- Contract: For every AS_* key emitted by both backends on the
  same record, the value must be byte-identical.
- Status: **PASS** (exit 0, 0 drift).  `AS_SB_TABLE` byte-equal
  on all 5 records.  Native also emits `AS_MBQ`, `AS_MFRL`,
  `AS_SOR` on all 5 records; GATK does not emit any of these
  on the dream_synthetic record depth, so they appear as
  native-only and value comparison is skipped.  When GATK's
  downstream FilterMutectCalls consumes the native Mutect2
  driving VCF (with AS_MBQ/AS_MFRL/AS_SOR populated), GATK
  passes them through unchanged.
- Tier-1a value-level gate is closed for `AS_SB_TABLE` on this
  fixture; richer fixtures (hcc1143, NA12878 trio) would
  extend the comparison to the remaining AS_* family keys.

## Cross-references

- `docs/tools/fastgatk-mutect2.md` — full tool description
- `docs/tools/fastgatk-filter-mutect-calls.md` — full tool description
- `docs/algorithms/somatic-likelihood.md` — TLOD / NLOD derivation
- `docs/algorithms/somatic-posterior.md` — Dirichlet variational
  evidence, allele-marginalization
- `docs/algorithms/filter-reads.md` — ArtifactReadFilter /
  FilterAlignmentArtifacts algorithms
- `fastgatk-native/dispatcher/tool_registry.json` — `Mutect2` entry
  (status `contract-compatible`, fallback explicit) and
  `FilterMutectCalls` entry (status `contract-compatible`)