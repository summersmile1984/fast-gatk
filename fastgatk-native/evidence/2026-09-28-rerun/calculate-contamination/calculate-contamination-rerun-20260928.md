# CalculateContamination rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 31.306s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_calculate_contamination.py` | 0 | 0.069s | 606.98s | 5.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_calculate_contamination_gatk_oracle.py` | 0 | 31.237s | 637.29s | 11.7 MiB | 649.0 MiB | 29.07s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool calculate-contamination
```

## JSON sidecar

See `calculate-contamination-rerun-20260928.json` for full stdout/stderr tails.
