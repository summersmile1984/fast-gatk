# CalculateContamination rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 34.792s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_calculate_contamination.py` | 0 | 0.053s | 542.57s | 5.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_calculate_contamination_gatk_oracle.py` | 0 | 34.739s | 571.60s | 11.8 MiB | 650.4 MiB | 32.64s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool calculate-contamination
```

## JSON sidecar

See `calculate-contamination-rerun-20260923.json` for full stdout/stderr tails.
