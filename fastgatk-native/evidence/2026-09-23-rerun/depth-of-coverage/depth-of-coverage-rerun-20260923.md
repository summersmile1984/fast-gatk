# DepthOfCoverage rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 131.323s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_depth_of_coverage.py` | 0 | 79.462s | 943.11s | 6.5 MiB | 377.8 MiB | 79.33s | ✅ |
| `verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py` | 0 | 35.776s | 910.96s | 6.8 MiB | 356.1 MiB | 17.13s | ✅ |
| `verify_depth_of_coverage_multisample.py` | 0 | 0.059s | 892.32s | 5.9 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_depth_of_coverage_read_filter_gatk_oracle.py` | 0 | 16.026s | 899.77s | 5.7 MiB | 345.9 MiB | 15.97s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool depth-of-coverage
```

## JSON sidecar

See `depth-of-coverage-rerun-20260923.json` for full stdout/stderr tails.
