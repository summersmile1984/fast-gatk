# GatherVcfs rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 21.317s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gather_vcfs.py` | 0 | 18.378s | 2154.99s | 5.7 MiB | 372.9 MiB | 18.30s | ✅ |
| `verify_gather_vcfs_cli_boundary_gatk_oracle.py` | 0 | 2.939s | 2141.84s | 0.0 MiB | 284.6 MiB | 2.91s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool gather-vcfs
```

## JSON sidecar

See `gather-vcfs-rerun-20260928.json` for full stdout/stderr tails.
