# Verifications index

Each tool/algorithm has a small number of oracle scripts. This
directory indexes them so a new contributor can find:

- "is there already an oracle for X?" — read the relevant
  `contracts-*.md`
- "where does GATK source come from for oracle X?" — link to
  `third_party/gatk-package/gatk-4.6.2.0/gatkdoc/<class>.json`
- "when did X last pass against native?" — `evidence/2026-09-25-rerun/`
  has the JSON rerun result

## Subindex files

- `contracts-hc.md` — HaplotypeCaller oracles
- `contracts-mutect2.md` — Mutect2 + FilterMutectCalls + LearnReadOrientationModel
- `contracts-bqsr-vqsr.md` — BQSR + ApplyBQSR + VariantRecalibrator + ApplyVQSR
- `contracts-picard.md` — MarkDuplicates, SortSam, GatherVcfs
- `contracts-cnv.md` — CollectReadCounts, DenoiseReadCounts, CreateReadCountPanelOfNormals, ModelSegments, CallCopyRatioSegments
- `contracts-rerun-protocol.md` — how to re-run all oracles

## How to add a new oracle

1. Create `fastgatk-native/scripts/verify_<thing>.py`. Required
   skeleton (see `templates/verify_oracle_skeleton.py` or a real example):

   ```python
   #!/usr/bin/env python3
   """Pin <thing>: <one-line contract>.

   Source-of-truth: <GATK source path>.
   Pinned fixtures: <fixture list with SHA1 if pinned>.
   Contract: <what byte-identical means here>.

   TDD: this oracle MUST fail against the current native binary (the
   gap is un-shipped) and pass after the corresponding Tier-N code
   change.
   """
   from __future__ import annotations
   import json, os, subprocess, tempfile
   from pathlib import Path
   import oracle_guard
   ROOT = Path(__file__).resolve().parents[2]
   # ... define fixtures, run GATK Java, run native, assert contract
   ```

2. Add entry to the appropriate `contracts-<group>.md` file.

3. Add a "this gap needs an oracle" entry to
   `work/mutect2-priority/NEXT_STEPS.md` if it doesn't have one.

## Re-running

```bash
# Full rerun of one tool
python3 fastgatk-native/scripts/rerun_all_verify.py --tool <tool> --repo .

# Single oracle
python3 fastgatk-native/scripts/verify_<thing>.py
```

## TDD rule

Before opening the code editor for any gap, the corresponding oracle
script MUST exist AND MUST exit non-zero against the current binary.
That negative exit is the gate that proves the gap is real and
un-shipped.  The corresponding algorithm doc MUST also exist.

After the implementation lands, the same oracle must exit zero.  If
it does not, the implementation is wrong, not the oracle.
