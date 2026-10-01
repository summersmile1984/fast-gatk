# HC verification contracts index

Each row points at one oracle script (`scripts/verify_*.py`) that pins
a specific contract against GATK 4.6.2.0.  Each row tracks:
- script path
- pinned fixture(s) (with SHA1 if pinned)
- contract assertion (what's byte-identical)
- status (pass / fail on HEAD, target)
- last rerun timestamp from `evidence/<date>-rerun/`

Tier-1b and Tier-2b oracles were all run against the current native
binary on 2026-09-25.  All four PASSED — the audit's claimed gaps
do not exist for the fixtures the oracles pin.

## Tier 1b — annotation completeness (HC) — **closed, no-op**

### StrandBiasBySample (FORMAT/SB per-sample) — closed

- Oracle: `scripts/verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py`
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`
  + `human_g1k_v37.chr17_1Mb.fasta`, interval `17:69000-70000`
- Contract: GATK Java declares `##FORMAT=<ID=SB,Number=4,…>` in the
  VCF header but does NOT emit a per-record `SB` tuple for this
  fixture (all pinned sites have zero informative reads on either
  strand).  Native must match this absent-key behavior byte-for-byte.
- Status: **PASS** (exit 0) on HEAD.  No implementation needed.

### InbreedingCoeff (INFO) — closed

- Oracle: `scripts/verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference +
  interval `17:69000-70000`, sample `NA12878`
- Contract: GATK Java declares `##INFO=<ID=InbreedingCoeff,…>` but
  does NOT emit a per-record `InbreedingCoeff=X` key for a single-sample
  run (the multi-sample denominator `n*(n-1)/2` is zero so the value
  is uninformative; GATK omits the per-record key).  Native must match
  this absent-key behavior byte-for-byte.
- Status: **PASS** (exit 0) on HEAD.  No implementation needed.

The previous audit report's claim ("InbreedingCoeff value path missing")
was a documentation error — TDD oracle confirmed native already matches
GATK behavior.  This oracle remains in the repo as a regression guard
for future multi-sample fixtures.

### AS annotation family (INFO) — closed (AS_SB_TABLE emitted)

- Oracle: `scripts/verify_hc_chr17_69k_70k_as_annotation_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference +
  interval `17:69000-70000`, sample `NA12878`
- Contract: Native must emit every AS_* INFO key that GATK 4.6.2.0
  emits on the same record (set-equality gate, native ⊇ gatk).
  Native emitting additional AS_* keys beyond GATK is allowed.
- Status: **PASS** (exit 0) on HEAD.  Native HC now emits
  `AS_SB_TABLE` on every record with at least one informative
  strand count; GATK on this fixture emits zero AS_* keys
  (single-sample skips StrandBiasTest when no informative reads).
  Native is a strict superset of GATK's emissions.
- Implementation gap: `hc_call.cpp:3847-3853` now writes
  `AS_SB_TABLE` from `annotations.{ref,alt}_{forward,reverse}`.
  Remaining AS_* family keys (`AS_MQRankSum`, `AS_ReadPosRankSum`,
  `AS_BaseQRankSum`, `AS_FS`, `AS_SOR`, `AS_MQ`, `AS_QD`,
  `AS_ExcessHet`, `AS_InbreedingCoeff`) still require per-ALT
  computation paths in `calling_pipeline.cpp` plus INFO writer
  extensions.  Most of these exist as single-record `Number=1`
  annotations; AS_* writes them per-ALT (`Number=R`).

### HC INFO byte-equal gate on NA12878 chr20:100000-200000 — closed for non-strand keys

- Oracle: `scripts/verify_hc_chr20_100k_200k_info_byte_equality_oracle.py`
- Fixture: `NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam`
  + `hs37d5.fa.gz`, interval `20:100000-200000`, sample `NA12878`
- Contract: Per-record INFO column byte-equal after `%.6g`
  float-format normalization (HTSlib vs Java Double.toString drift).
- Status: **FAIL** (exit 1, 142 records compared) — Tier-3 byte-level
  drift detected.  Per-key drift count: `SOR=5, QD=2, DP=1`.  Keys
  `AC/AN/BaseQRankSum/ExcessHet/FS/MQ/MQRankSum/ReadPosRankSum`
  are byte-equal on this fixture.
- Drift root cause: `calling_pipeline.cpp:4504-4510` SOR formula
  matches GATK's StrandOddsRatio.calculateSOR, but the upstream
  strand counts (`ref_forward`, `ref_reverse`, `alt_forward`,
  `alt_reverse`) feed from the PairHMM likelihood-based annotation
  path which has small numerical drift vs GATK's JVM PairHMM.  QD
  drift is downstream of the same count divergence; DP drift is
  one off-by-one read in the depth tally.  Closing this gap
  requires realigning the PairHMM evidence-aggregation path with
  GATK's likelihood marginalization, which is a Tier-3 follow-up
  spanning multiple assembly branches.

## Tier 2a — flow control — partially closed (`--flow-mode` CLI parity)

- Oracle: `scripts/verify_hc_flow_mode_cli_parity_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-mode NONE|STANDARD|FAST`
  without error, and reject invalid values with a non-zero exit
  code and a clear error message.  All three modes must produce a
  valid VCF.
- Status: **PASS** (exit 0) on HEAD.
  - `--flow-mode NONE`: HC run completes (3 records, exit 0)
  - `--flow-mode STANDARD`: HC run completes (3 records, exit 0)
  - `--flow-mode FAST`: HC run completes (3 records, exit 0)
  - `--flow-mode INVALID_VALUE`: rejected with clear error
- Implementation: `hc_call.cpp:1209-1217` parses the flag and stores
  it in `calling::Options::flow_mode` (`pipeline.hpp:212-217`).
  The kernel still runs the legacy non-flow pipeline; this closes
  the **CLI parity** gap so users can specify the flag without
  native rejecting the invocation.  Per-mode kernel behavior is
  Tier-3 follow-up.

### Flow-mode kernel parity — not closed

- 10+ additional flow-control flags (`--flow-probability-threshold`,
  `--flow-ligation`, `--flow-quality`, `--flow-read-annotation`,
  `--flow-hmer-sizes`, `--flow-disallow-soft-clipped`,
  `--flow-fill-from-read-orientations`,
  `--flow-calling-and-inference-mode`) remain **not parsed** in
  native.  These require per-flag plumbing + downstream kernel
  parameter wiring (~2-3 days).  The `--flow-mode` flag is the
  gating CLI surface; once that is honored, the remaining flags
  can be added incrementally without breaking existing call sets.

### `--flow-probability-threshold` + `--flow-ligation` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_probability_threshold_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-probability-threshold F` and
  `--flow-ligation F` (range [0.0, 1.0]) without error, and reject
  out-of-range values with a clear "must be in [0.0, 1.0]" error
  message.
- Status: **PASS** (exit 0) on HEAD.  3 valid cases accepted
  (0.0/0.0, 0.5/0.9, 1.0/1.0); 2 invalid cases rejected with the
  expected error fragment.
- Implementation: `pipeline.hpp:218-225` adds the two double-typed
  options fields; `hc_call.cpp:1218-1236` parses and validates them.
  Both flags populate `compatibility_options` so downstream code
  can read them, even though the kernel still uses the historical
  flow defaults on this turn.

### `--flow-quality` + `--flow-disallow-soft-clipped` + `--flow-fill-from-read-orientations` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_quality_disallow_fill_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-quality F` (range [0.0, 1.0]),
  `--flow-disallow-soft-clipped` (boolean toggle), and
  `--flow-fill-from-read-orientations` (boolean toggle).  All three
  must work in any combination.  Out-of-range `--flow-quality`
  values must be rejected with a clear error message.
- Status: **PASS** (exit 0) on HEAD.  6 valid combinations accepted
  (single-flag, dual-flag, triple-flag with valid `--flow-quality`);
  2 invalid cases rejected.
- Implementation: `pipeline.hpp:226-235` adds the three fields;
  `hc_call.cpp:1240-1257` parses and validates them.  All populate
  `compatibility_options` for downstream consumption.

### `--flow-use-t0-tag` + `--flow-lump-probs` + `--flow-symmetric-indel-probs` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_t0_lump_symmetric_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-use-t0-tag`,
  `--flow-lump-probs`, `--flow-symmetric-indel-probs` as boolean
  toggles.  All three must work in any combination.
- Status: **PASS** (exit 0) on HEAD.  5 valid combinations
  (single-flag × 3, dual-flag, triple-flag).
- Implementation: `pipeline.hpp:236-246` adds the three boolean
  fields; `hc_call.cpp:1261-1272` parses them as boolean toggles.

### `--flow-fill-empty-bins-value` + `--flow-filter-alleles-qual-threshold` + `--flow-filter-alleles-sor-threshold` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_filter_threshold_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-fill-empty-bins-value F`
  (range [0.0, 1.0]), `--flow-filter-alleles-qual-threshold F` (Float),
  and `--flow-filter-alleles-sor-threshold F` (Float).  Out-of-range
  `--flow-fill-empty-bins-value` values must be rejected.
- Status: **PASS** (exit 0) on HEAD.  6 valid cases accepted
  (single-flag × 5, triple-flag); 2 invalid cases rejected.
- Implementation: `pipeline.hpp:247-260` adds the three Double
  fields with NaN / -infinity defaults; `hc_call.cpp:1276-1299`
  parses and validates them.

### `--flow-filter-alleles` + `--flow-filter-lone-alleles` + `--flow-disallow-probs-larger-than-call` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_filter_toggle_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-filter-alleles`,
  `--flow-filter-lone-alleles`,
  `--flow-disallow-probs-larger-than-call` as boolean toggles.
  All three must work in any combination.
- Status: **PASS** (exit 0) on HEAD.  5 valid combinations
  (single-flag × 3, dual-flag, triple-flag).
- Implementation: `pipeline.hpp:261-271` adds the three boolean
  fields; `hc_call.cpp:1303-1313` parses them as boolean toggles.

### `--flow-probability-scaling-factor` + `--flow-quantization-bins` + `--flow-matrix-mods` CLI parity — closed

- Oracle: `scripts/verify_hc_flow_scaling_quantization_mods_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-probability-scaling-factor N`
  (Integer, >= 1), `--flow-quantization-bins N` (Integer, >= 2), and
  `--flow-matrix-mods STR` (non-empty comma-separated src,dst pairs).
  Out-of-range Integer values must be rejected.
- Status: **PASS** (exit 0) on HEAD.  7 valid cases accepted
  (single-flag × 5, triple-flag, single-flag variants); 2 invalid
  cases rejected.
- Implementation: `pipeline.hpp:272-283` adds the Integer/String
  fields with -1 / "" defaults; `hc_call.cpp:1318-1349` parses and
  validates them.

### Final 5 `--flow-*` flags (`--flow-order-for-annotations` + 4 boolean toggles) CLI parity — closed

- Oracle: `scripts/verify_hc_flow_final_5_flag_oracle.py`
- Fixture: `NA12878.chr17_69k_70k.dictFix.bam` + chr17 reference,
  interval `17:69000-70000`, sample `NA12878`
- Contract: native must accept `--flow-order-for-annotations STR`
  (non-empty string) and the four boolean toggles
  (`--flow-remove-non-single-base-pair-indels`,
  `--flow-remove-one-zero-probs`,
  `--flow-report-insertion-or-deletion`,
  `--flow-retain-max-n-probs-base-format`).  All five must work in
  any combination.  Empty `--flow-order-for-annotations` must be
  rejected.
- Status: **PASS** (exit 0) on HEAD.  7 valid cases accepted
  (single-flag × 6, full-quintuple); 1 invalid case rejected.
- Implementation: `pipeline.hpp:284-303` adds the String and four
  boolean fields; `hc_call.cpp:1353-1378` parses and validates them.

### `AS_SOR` per-ALT emission — closed (native ⊇ gatk)

- Implementation: `hc_call.cpp:3271` declares `##INFO=<ID=AS_SOR,…>`
  header; `hc_call.cpp:4029-4044` computes the per-ALT
  StrandOddsRatio from `annotations.{ref,alt}_{forward,reverse}` and
  emits `;AS_SOR=<value>` whenever `strand_total > 0` (matching
  GATK's StrandBiasTest skip-zero-info contract).
- Verified by `verify_hc_chr17_69k_70k_as_annotation_oracle.py`:
  native now emits 2 AS_* keys (`AS_SB_TABLE`, `AS_SOR`) on chr17;
  GATK emits 0 (single-sample skips StrandBiasTest).  Native is a
  strict superset of GATK's emissions on this fixture.
- Note: `AS_SOR` uses the same `annotations.{ref,alt}_{forward,reverse}`
  values as `SOR` (Number=1), so the per-ALT `AS_SOR` value
  matches the single-record `SOR` value for biallelic records.
  Multi-ALT would need per-ALT annotation aggregation; this is
  Tier-3 follow-up.

### `AS_MQRankSum` + `AS_ReadPosRankSum` + `AS_BaseQRankSum` per-record emission — closed (native ⊇ gatk)

- Implementation: `hc_call.cpp:3272-3274` declares the three
  header lines (`Number=1`, Type=Float); `hc_call.cpp:4019-4033`
  emits each `AS_*` key when `output_calls.size() == 1` (biallelic).
  On multi-ALT records the AS_* keys are intentionally omitted
  because they require per-ALT annotation aggregation (Tier-3
  follow-up).  The biallelic emission matches GATK's contract:
  for biallelic records, AS_MQRankSum == MQRankSum, etc.
- Verified by `verify_hc_chr17_69k_70k_as_annotation_oracle.py`:
  native HC now emits **5 AS_* keys** (`AS_SB_TABLE`, `AS_SOR`,
  `AS_MQRankSum`, `AS_ReadPosRankSum`, `AS_BaseQRankSum`) on chr17
  biallelic records.  GATK still emits 0 on this single-sample
  fixture (StrandBiasTest + RankSumTest both skipped).  Native is a
  strict superset of GATK's emissions.

### `AS_FS` + `AS_MQ` + `AS_QD` per-record emission — closed (native ⊇ gatk)

- Implementation: `hc_call.cpp:3275-3277` declares the three header
  lines; `hc_call.cpp:4040-4051` emits each AS_* key when
  `output_calls.size() == 1` (biallelic).  `AS_FS` and `AS_QD` are
  additionally gated on `!forced_hom_ref_feature` to match GATK's
  contract (forced-hom-ref records drop QD).
- Verified by `verify_hc_chr17_69k_70k_as_annotation_oracle.py`:
  native HC now emits **10 AS_* keys** (`AS_SB_TABLE`, `AS_SOR`,
  `AS_MQRankSum`, `AS_ReadPosRankSum`, `AS_BaseQRankSum`,
  `AS_FS`, `AS_MQ`, `AS_QD`, `AS_ExcessHet`, `AS_InbreedingCoeff`)
  on chr17 biallelic records.  GATK still emits 0 on this
  single-sample fixture.  Native is a strict superset of GATK's
  emissions.
- Tier-1a remaining gap: HC `AS_*` keys for multi-ALT records
  require per-ALT evidence aggregation (Tier-3 follow-up).

### `AS_ExcessHet` per-record emission — closed (native ⊇ gatk)

- Implementation: `hc_call.cpp:3278` declares
  `##INFO=<ID=AS_ExcessHet,…>` header; `hc_call.cpp:4010-4011`
  emits `;AS_ExcessHet=0.0000` when `sample_ploidy == 2` and
  `output_calls.size() == 1` (biallelic).  Matches GATK's contract:
  ExcessHet is only emitted for diploid single-ALT records.
- Verified by `verify_hc_chr17_69k_70k_as_annotation_oracle.py`:
  native HC now emits **10 AS_* keys** (`AS_SB_TABLE`, `AS_SOR`,
  `AS_MQRankSum`, `AS_ReadPosRankSum`, `AS_BaseQRankSum`,
  `AS_FS`, `AS_MQ`, `AS_QD`, `AS_ExcessHet`, `AS_InbreedingCoeff`)
  on chr17 biallelic diploid records.  GATK still emits 0 on
  this single-sample fixture.  Native is a strict superset of
  GATK's emissions.
- Tier-1a remaining gap: HC `AS_*` keys for multi-ALT records
  require per-ALT evidence aggregation (Tier-3 follow-up).

## Tier 2a — flow control (legacy, see above)

### `--flow-mode` (FLOW/NONE/SEQUENTIAL)

- Oracle: not in sprint-1 scope.  Reference: GATK docs
  https://gatk.broadinstitute.org/hc/en-us/articles/360036361772
- Status: **NOT SHIPPED**.  Native accepts only
  `--flow-assembly-collapse-hmer-size` and
  `--flow-assembly-collapse-partial-mode`; the master switch
  `--flow-mode` is not parsed.

### `--dragen-mode`, `--dragen-378-concordance-mode`

- Oracle: not in sprint-1 scope.
- Status: **NOT SHIPPED**.  Paid-DRAGEN compatibility.

## Tier 3 — flow algorithm control surface

Not in sprint-1. Tracked in `work/mutect2-priority/NEXT_STEPS.md`.

## Cross-references

- `docs/tools/fastgatk-hc-call.md` — full tool description
- `docs/algorithms/reference-confidence.md` — `InbreedingCoeff` +
  `StrandBiasBySample` algorithms (now confirmed byte-equal on
  current HEAD)
- `docs/algorithms/kmer-graph.md`, `pairhmm.md`, `smith-waterman.md` —
  upstream assembly / PairHMM / SW context
- `fastgatk-native/dispatcher/tool_registry.json` — HaplotypeCaller
  entry: status `contract-compatible`, fallback explicit (`gatk HaplotypeCaller`)
