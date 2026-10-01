# Verify-script re-run index — 2026-09-27

- Tools covered: **1** (excluding `genomicsdb-export` (native-only bridge, no verify script))
- Scripts total: **25**
- Passed: **25**    Failed: **0**    Skipped: **0**
- Total elapsed: **2091.244s**

Per-tool summaries (cross-ref with `fastgatk-native/docs/cli-alignment.md` *Bit-identical / bounded-parity evidence* table):

| Tool | Native binary | Scripts | Passed | Failed | Skipped | Elapsed (s) | Report |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `GenotypeGVCFs` | `fastgatk-genotype-gvcf` | 25 | 25 | 0 | 0 | 2091.244 | [genotype-gvcf-rerun-20260927.md](fastgatk-native/evidence/2026-09-27-rerun/genotype-gvcf/genotype-gvcf-rerun-20260927.md) |

## Bridge binary with no verify script

- `fastgatk-genomicsdb-export` — (native-only bridge, no verify script); cross-reference the docs.

## GATK-alignment status

Run `aggregate_alignment_status.py` after a re-run to cross-reference
this index with the *Bit-identical / bounded-parity evidence* table
in `fastgatk-native/docs/cli-alignment.md`.  It writes
`ALIGNMENT_STATUS.md` next to this file.

## Native-vs-Java resource comparison

Run `compare_native_vs_java.py` to classify every descendant process
of the rerun into `native` (comm starts with `fastgatk-`) vs `java`
(comm == `java`) buckets, and emit `NATIVE_VS_JAVA.md` next to this
file.  Anomalies (native slower than Java, or native RSS exceeding
Java's) are auto-flagged in that report.

## Workflow documentation

Full workflow, schema details, and an anomaly follow-up playbook
live in `fastgatk-native/docs/regression-evidence.md`.

## Reproduction

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --repo .
```

```bash
python3 fastgatk-native/scripts/verify_rerun_report.py --repo .
```

```bash
python3 fastgatk-native/scripts/aggregate_alignment_status.py --repo .
```

```bash
python3 fastgatk-native/scripts/compare_native_vs_java.py --repo .
```
