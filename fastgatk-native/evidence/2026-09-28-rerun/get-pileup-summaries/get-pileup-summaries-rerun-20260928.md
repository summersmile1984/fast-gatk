# GetPileupSummaries rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 45.021s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_get_pileup_gatk_oracle.py` | 0 | 44.904s | 3974.44s | 7.2 MiB | 368.5 MiB | 41.20s | ✅ |
| `verify_get_pileup_summaries.py` | 0 | 0.117s | 3941.94s | 6.7 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool get-pileup-summaries
```

## JSON sidecar

See `get-pileup-summaries-rerun-20260928.json` for full stdout/stderr tails.
