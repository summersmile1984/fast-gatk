# GetPileupSummaries rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 81.738s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_get_pileup_gatk_oracle.py` | 0 | 81.577s | 2845.91s | 7.1 MiB | 368.2 MiB | 76.40s | ✅ |
| `verify_get_pileup_summaries.py` | 0 | 0.161s | 2814.43s | 7.5 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool get-pileup-summaries
```

## JSON sidecar

See `get-pileup-summaries-rerun-20260923.json` for full stdout/stderr tails.
