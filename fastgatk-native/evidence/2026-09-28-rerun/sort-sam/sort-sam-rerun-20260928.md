# SortSam rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 27.479s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_sort_sam.py` | 0 | 7.583s | 7692.37s | 0.0 MiB | 272.9 MiB | 7.54s | ✅ |
| `verify_sort_sam_cli_boundary_gatk_oracle.py` | 0 | 1.278s | 7680.84s | 0.0 MiB | 250.5 MiB | 1.25s | ✅ |
| `verify_sort_sam_duplicate_gatk_oracle.py` | 0 | 3.899s | 7685.09s | 0.0 MiB | 270.2 MiB | 3.86s | ✅ |
| `verify_sort_sam_gatk_oracle.py` | 0 | 14.719s | 7705.12s | 1.9 MiB | 287.7 MiB | 14.68s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool sort-sam
```

## JSON sidecar

See `sort-sam-rerun-20260928.json` for full stdout/stderr tails.
