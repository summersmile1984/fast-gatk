# DepthOfCoverage rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 87.494s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_depth_of_coverage.py` | 0 | 58.446s | 1008.35s | 7.0 MiB | 371.0 MiB | 58.32s | ✅ |
| `verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py` | 0 | 20.571s | 977.38s | 0.0 MiB | 399.0 MiB | 20.54s | ✅ |
| `verify_depth_of_coverage_multisample.py` | 0 | 0.046s | 959.12s | 6.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_depth_of_coverage_read_filter_gatk_oracle.py` | 0 | 8.431s | 966.36s | 0.0 MiB | 365.9 MiB | 8.38s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool depth-of-coverage
```

## JSON sidecar

See `depth-of-coverage-rerun-20260928.json` for full stdout/stderr tails.
