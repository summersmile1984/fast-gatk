# GatherVcfs rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 38.56s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gather_vcfs.py` | 0 | 34.548s | 1367.47s | 5.9 MiB | 364.1 MiB | 34.46s | ✅ |
| `verify_gather_vcfs_cli_boundary_gatk_oracle.py` | 0 | 4.012s | 1351.68s | 0.2 MiB | 289.2 MiB | 3.98s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool gather-vcfs
```

## JSON sidecar

See `gather-vcfs-rerun-20260923.json` for full stdout/stderr tails.
