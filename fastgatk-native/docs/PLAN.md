# Master plan — fastgatk-native ↔ GATK 4.6.2.0 algorithm & tool parity

**Date:** 2026-09-25
**Owner:** fastgatk-native maintainers
**Status:** In execution
**Inputs:** [`PARITY_AUDIT.md`](../../../work/mutect2-priority/PARITY_AUDIT.md), [`NEXT_STEPS.md`](../../../work/mutect2-priority/NEXT_STEPS.md)

## Mission

Bring every shipped `fastgatk-*` tool to **byte-identical-or-better parity** with GATK 4.6.2.0 on its supported oracle set, and produce a complete documentation set for every tool and algorithm — **before writing any new feature code**.

Two non-negotiable rules:

1. **TDD-first**: an oracle (verify script) runs against both GATK Java and the native binary, asserts the contract, and **must fail before the new code exists**.
2. **Doc-first for each gap**: for every algorithm and tool with an open parity gap, a `*.md` doc in `docs/algorithms/` or `docs/tools/` exists **before** the gap is closed. The doc is the spec.

## Repo layout (after this plan)

```
fastgatk-native/docs/
├── PLAN.md                          # this file
├── cli-alignment.md                 # existing CLI parity table (do not edit; regenerate via generate_cli_alignment_doc.py)
├── regression-evidence.md           # existing rerun protocol
├── tools/                           # NEW: per-tool doc — one md per shipped tool
│   ├── _TEMPLATE.md
│   ├── fastgatk-hc-call.md
│   ├── fastgatk-mutect2.md
│   ├── ...
├── algorithms/                      # NEW: per-algorithm doc — one md per kernel & pipeline stage
│   ├── _TEMPLATE.md
│   ├── pairhmm.md
│   ├── smith-waterman.md
│   ├── activity-profile.md
│   ├── kmer-graph.md
│   ├── haplotype-realn-cigar.md
│   ├── filter-reads.md
│   ├── pileup-error-correction.md
│   ├── reference-confidence.md
│   ├── bqsr.md
│   ├── vqsr.md
│   ├── cnv-denoise.md
│   ├── cnv-model-segments.md
│   ├── somatic-likelihood.md
│   ├── somatic-posterior.md
│   ├── fragment-aggregation.md
│   ├── allele-marginalization.md
│   ├── normal-artifact-orientation.md
│   ├── pileup-summary.md
│   ├── contamination.md
│   ├── gvcf-band-blocks.md
│   ├── combine-gvcfs.md
│   ├── genotype-gvcfs.md
│   ├── reblock-gvcfs.md
│   ├── select-variants.md
│   ├── validate-variants.md
│   ├── variant-eval.md
│   ├── variant-filtration.md
│   ├── variants-to-table.md
│   ├── mark-duplicates.md
│   ├── sort-sam.md
│   ├── gather-vcfs.md
│   ├── gather-bqsr.md
│   ├── analyze-covariates.md
│   ├── apply-bqsr.md
│   ├── f1r2-orientation.md
│   └── flow-hmer-collapse.md
└── verifications/                   # NEW: per-oracle-script index & contract registry
    ├── _README.md
    ├── contracts-hc.md
    ├── contracts-mutect2.md
    ├── contracts-bqsr-vqsr.md
    ├── contracts-picard.md
    ├── contracts-cnv.md
    └── contracts-rerun-protocol.md
```

## Master workflow — TDD-first for every gap

```
┌──────────────────────────────────────────────────────────────────────┐
│  1. ENUMERATE GAPS from docs/PARITY_AUDIT.md                         │
│     each gap is one Item: (tool|algo, gap, status)                  │
│                                                                      │
│  2. WRITE/UPDATE PER-ALGORITHM DOC at docs/algorithms/<name>.md     │
│     required sections per algorithm (see §1 doc template)            │
│     do NOT open code editor until doc merges                          │
│                                                                      │
│  3. WRITE/UPDATE PER-TOOL DOC at docs/tools/<binary>.md              │
│     required sections per tool (see §2 doc template)                  │
│                                                                      │
│  4. WRITE ORACLE at fastgatk-native/scripts/verify_*.py             │
│     contract: pin GATK 4.6.2.0 behavior with byte-identical VCF/  │
│     byte-identical INFO/BAM rows; or semantic equivalence for non- │
│     byte fields. Must exit 0 against current native binary (i.e.   │
│     FAILS today because the gap is un-shipped)                       │
│                                                                      │
│  5. RUN ORACLE → confirm it fails as expected                        │
│                                                                      │
│  6. INDEX in docs/verifications/contracts-<group>.md                 │
│     pin status: oracle existing, gap status, last run timestamp       │
│                                                                      │
│  7. OPEN CODE → minimum patch satisfying the oracle contract only   │
│     do NOT add features not asserted by the oracle                  │
│                                                                      │
│  8. RE-RUN ORACLE → confirm pass                                    │
│                                                                      │
│  9. UPDATE DOC with measured deltas, new flag surface, evidence       │
│                                                                      │
│ 10. RUN FULL RE-RUN via rerun_all_verify.py to catch regressions     │
└──────────────────────────────────────────────────────────────────────┘
```

Every gap MUST complete all 10 steps before the next gap starts.

## 1. Doc template — algorithm (`docs/algorithms/<name>.md`)

One md per algorithm (PairHMM, Smith-Waterman, ActivityProfile, kmer graph, haplotype SW, read filter, pileup correction, reference confidence, BQSR, VQSR, somatic likelihood, somatic posterior, fragment aggregation, allele marginalization, normal-artifact orientation, pileup summary, contamination, GVCF band blocks, combine GVCFs, genotype GVCFs, reblock GVCFs, select variants, validate variants, variant eval, variant filtration, variants-to-table, mark-duplicates, sort-sam, gather-vcfs, gather-bqsr, analyze-covariates, apply-bqsr, f1r2-orientation, flow-hmer-collapse).

```markdown
# <algorithm name>

**Source:** `fastgatk-kernels/<relative-path>.cpp` (or `fastgatk-native/src/...`)
**GATK reference:** `https://github.com/broadinstitute/gatk/blob/4.6.2.0/<gatk-path>.java`
**Status:** shipped | partial | not-shipped
**Tested by oracle:** `scripts/verify_*.py` (link to file in verifications/)

## 1. Purpose

What the algorithm computes.  GATK Java class name.  Math formula
in TeX if non-trivial.  One-paragraph summary.

## 2. Input / output contract

Inputs: types, shapes, NaN/null conventions, edge cases.
Outputs: same.

## 3. Reference implementation

Link to GATK source.  Excerpt of the Java method (or the math
reference, e.g. Wikipedia for Smith-Waterman).  Document the
**exact numerical policy** (log10 vs natural, FTZ, exact-zero
substitution, etc.).  Anything the oracle must verify.

## 4. Native implementation

Architecture summary: which Kokkos kernels, where on
Host/device.  Algorithmic differences from GATK (must be explicit).
GPU-safety analysis: any host-only intrinsics, any pure-host
codepaths.

## 5. Oracle contract

The byte-level (or semantic-equivalent) fields the oracle
verifies.  Cross-reference: `docs/verifications/contracts-*.md`.

## 6. Known parity gaps

For each open gap:
- gap description
- GATK source path
- native state
- blocking oracle (the next step 4)
- priority tier

## 7. Reference commands

CLI invocations for both GATK Java and native that exercise this
algorithm.  Output artifact paths.  Command-line flags relevant to
the algorithm.

## 8. Verification anchors

The shell + python one-liners that prove the contract.  Oracle
script path.  Manual reproduction snippet for offline debugging.

## 9. Cross-references

Tool docs that depend on this algorithm.  Other algorithm docs
that this one depends on.
```

## 2. Doc template — tool (`docs/tools/<binary>.md`)

One md per shipped `fastgatk-*` production binary (48 tools after
subtracting the 6 smoke / oracle / adapter binaries from `CMakeLists.txt`).

```markdown
# fastgatk-<tool-name>

**GATK equivalent:** `org.broadinstitute.hellbender.tools.<Tool>.java`
(URL)
**Source:** `fastgatk-native/src/<tool>_tool.cpp` (or wrapped file)
**Status:** shipped | partial | not-shipped
**Oracle(s):** `scripts/verify_<...>.py`

## 1. Purpose

One sentence.  When a user invokes this tool.

## 2. CLI parity

Table or link to `docs/cli-alignment.md` section.  Highlight only
the flags that differ from GATK 4.6.2.0.

## 3. Algorithm(s) used

Bullet list of `docs/algorithms/<name>.md` linked.

## 4. Output contract

What the tool writes: file format, compression, headers, sample
ordering, etc.

## 5. Oracle(s) and expected parity

For each oracle script that targets this tool:
- script path
- pinned fixture(s)
- what fields are byte-identical, what fields are semantic-equivalent
- last-verified timestamp from `evidence/<date>-rerun/`

## 6. Known gaps

Same shape as algorithm-doc §6.  Cross-reference to `docs/algorithms/`
and `PARITY_AUDIT.md`.

## 7. Build / install

CMake target.  Dependencies.  Kokkos backend requirement.

## 8. CLI quick reference

Example invocation with minimal flags.  Link to `cli-alignment.md`
for full option surface.

## 9. Cross-references

Pipeline-position: which upstream tools feed it, which downstream
tools consume its output.
```

## 3. Verifications index template (`docs/verifications/_README.md`)

```markdown
# Verifications index

Each tool/algorithm has a small number of oracle scripts.
This directory indexes them so a new contributor can find:
- "is there already an oracle for X?" — read the relevant
  `contracts-*.md`
- "where does GATK source come from for oracle X?" — link to
  `third_party/gatk-package/gatk-4.6.2.0/gatkdoc/<class>.json`
- "when did X last pass against native?" — `evidence/2026-09-25-rerun/`
  has the JSON rerun result

## How to add a new oracle

1. Create `fastgatk-native/scripts/verify_<thing>.py`.  Required
   skeleton:
   ```python
   #!/usr/bin/env python3
   """Pin <thing>: <one-line contract>.
   Source-of-truth: <GATK source path>.
   Pinned fixtures: <fixture list with SHA1 if pinned>.
   Contract: <what byte-identical means here>.
   TDD: this oracle MUST fail against the current binary (the gap is
   un-shipped) and pass after the corresponding Tier-N code change.
   """
   from __future__ import annotations
   import json, os, subprocess, tempfile
   from pathlib import Path
   import oracle_guard
   ROOT = Path(__file__).resolve().parents[2]
   # ... define fixtures, run GATK Java, run native, assert contract
   ```
2. Add entry to the appropriate `contracts-<group>.md` file.
3. Add a "this gap needs an oracle" entry to `NEXT_STEPS.md` if it
   doesn't have one.

## Re-running

`fastgatk-native/scripts/rerun_all_verify.py --tool <tool> --repo .`
```

## 4. Execution order (must respect TDD-first)

**Sprint-1 targets (decided 2026-09-25):**

- **Tier 1b** (1-2 days): HC `InbreedingCoeff` (header declared, value
  not computed) + HC `StrandBiasBySample` (FORMAT column; header
  declared at `hc_call.cpp:3941` but value never populated).
  Both shipped to `verify_hc_chr17_69k_70k` fixture's native VCF.

  - Algorithms: `reference-confidence.md` (§ ReferenceConfidenceModel
    + per-annotation extensions)
  - Tool: `fastgatk-hc-call.md`
  - Oracles (must FAIL today):
    - `scripts/verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py`
    - `scripts/verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py`

- **Tier 2b** (1-2 days each): Mutect2's `ArtifactReadFilter` /
  `FilterAlignmentArtifacts`.

  - Algorithms: `filter-reads.md`
  - Tool: `fastgatk-filter-mutect-calls.md`
  - Oracles (must FAIL today):
    - `scripts/verify_filter_mutect_artifact_read_filter_gatk_oracle.py`
    - `scripts/verify_filter_mutect_alignment_artifacts_gatk_oracle.py`

Per group, the order is:

1. **Documentation phase** (no code):
   - Update `docs/algorithms/<name>.md` for every algorithm in scope
   - Update `docs/tools/<binary>.md` for every tool in scope
   - Add `docs/verifications/contracts-<group>.md` index
2. **Oracle phase** (no new code):
   - Write `scripts/verify_*.py` for every gap; each must FAIL on
     current native binary
3. **Implementation phase**:
   - Pick the smallest patch that makes the oracle pass
   - No new features outside the oracle contract
4. **Verification phase**:
   - Re-run oracle → pass
   - Run `rerun_all_verify.py --tool <tool> --repo .` → no regressions
5. **Update phase**:
   - Doc §"Implementation" with measured deltas vs HEAD
   - Doc §"Known gaps" reflects remaining items only

## 5. Constraints (do NOT relax)

- **GPU-safety**: every kernel change must keep `fastgatk/kernels/gpu_safety.hpp` policy header happy.  No `<immintrin.h>`, no x86 intrinsics, no WFA/JNI.  Reviewed in `gpu_safety.hpp`.
- **Backend portability**: any change must build and run correctly on OpenMP, Threads, Serial, HPX, CUDA, HIP, SYCL.  Existing `rerun_all_verify.py` covers OpenMP+Serial; CUDA/HIP/SYCL CI is `KOKKOS_MODULE_VALIDATION.md`-tracked and not in scope here.
- **Backward compatibility**: oracle scripts MUST stay green for previously-shipped behavior.  No silent VCF schema changes.
- **Numerical contract**: any change to existing shipped algorithms must keep the per-algorithm oracle at byte-identical-or-better.  Drift >1e-8 requires Tier-2 review.
- **No new CLI flags** without an oracle asserting them.  Exception: `--disable-implementation-strict-mode` reserved for future Tier-3.

## 6. Sprint cadence

A sprint = 1 week (5 working days):

| Day | Activity |
|---|---|
| Mon | Doc updates (algorithms + tools + verifications index) |
| Tue | Oracle scripts for each gap |
| Wed-Thu | Implementation, smallest patch only |
| Thu EOD | Per-oracle re-run; commit doc + oracle + impl as one PR |
| Fri | Full `rerun_all_verify.py --tool <scope>`; if regressions, fix or document fallback |

## 7. Status

- [ ] Phase 0: All docs/algorithms/_TEMPLATE.md and docs/tools/_TEMPLATE.md
      + docs/verifications/_README.md written (this file is the plan)
- [ ] Phase 1: Tier 1b docs (reference-confidence, hc-call, contracts-hc)
- [ ] Phase 2: Tier 1b oracles (must fail on HEAD)
- [ ] Phase 3: Tier 1b implementation (smallest patch)
- [ ] Phase 4: Tier 1b verification (oracles pass, rerun green)
- [ ] Phase 5: Tier 2b docs
- [ ] Phase 6: Tier 2b oracles
- [ ] Phase 7: Tier 2b implementation
- [ ] Phase 8: Tier 2b verification
- [ ] Phase 9 (next sprint): Tier 1a AS annotation docs → oracles → impl
- [ ] Phase 10 (next+): Tier 3 full GATK families — scope decisions

## 8. Reviewers / owners

TBD; assign per-sprint.

## 9. References

- `PARITY_AUDIT.md` (this repo's `work/mutect2-priority/`) — full gap inventory
- `NEXT_STEPS.md` (this repo's `work/mutect2-priority/`) — ROI-ranked execution plan
- `RESEARCH.md` (this repo's `work/mutect2-priority/`) — algorithm literature
- `SUMMARY.md` (this repo's `work/mutect2-priority/`) — measured evidence
- `GATK_REMAINING_ALGORITHMS.md` — algorithm-gap P-class inventory
- `GATK_MODULE_INVENTORY.md` — GATK 4.6.2.0 command inventory
- `fastgatk-native/docs/cli-alignment.md` — CLI surface alignment matrix
- `fastgatk-native/dispatcher/tool_registry.json` — registry + fallback boundaries
- `evidence/<date>-rerun/` — most recent rerun snapshots
