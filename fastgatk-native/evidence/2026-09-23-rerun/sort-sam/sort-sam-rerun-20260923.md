# SortSam rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 41.845s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_sort_sam.py` | 0 | 11.610s | 6807.35s | 3.7 MiB | 271.0 MiB | 11.56s | ✅ |
| `verify_sort_sam_cli_boundary_gatk_oracle.py` | 0 | 1.318s | 6795.81s | 0.0 MiB | 269.5 MiB | 1.28s | ✅ |
| `verify_sort_sam_duplicate_gatk_oracle.py` | 0 | 5.397s | 6799.94s | 0.8 MiB | 259.4 MiB | 5.37s | ✅ |
| `verify_sort_sam_gatk_oracle.py` | 0 | 23.520s | 6820.52s | 2.2 MiB | 269.4 MiB | 18.07s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool sort-sam
```

## JSON sidecar

See `sort-sam-rerun-20260923.json` for full stdout/stderr tails.
