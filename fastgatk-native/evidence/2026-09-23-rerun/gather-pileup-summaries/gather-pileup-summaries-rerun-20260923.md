# GatherPileupSummaries rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 15.816s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gather_pileup_gatk_oracle.py` | 0 | 15.780s | 1339.54s | 4.8 MiB | 330.2 MiB | 15.75s | ✅ |
| `verify_gather_pileup_summaries.py` | 0 | 0.036s | 1330.73s | 3.2 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool gather-pileup-summaries
```

## JSON sidecar

See `gather-pileup-summaries-rerun-20260923.json` for full stdout/stderr tails.
