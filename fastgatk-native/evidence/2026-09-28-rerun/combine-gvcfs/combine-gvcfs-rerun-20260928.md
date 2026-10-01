# CombineGVCFs rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 87.367s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_combine_gvcfs.py` | 0 | 11.377s | 745.96s | 7.4 MiB | 349.0 MiB | 11.21s | ✅ |
| `verify_combine_gvcfs_gatk_oracle.py` | 0 | 35.797s | 793.08s | 7.0 MiB | 404.8 MiB | 35.76s | ✅ |
| `verify_combine_gvcfs_interval_refblock_gatk_oracle.py` | 0 | 22.452s | 773.76s | 7.1 MiB | 404.5 MiB | 22.40s | ✅ |
| `verify_combine_gvcfs_plless_oracle.py` | 0 | 17.741s | 757.74s | 5.6 MiB | 400.8 MiB | 17.70s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool combine-gvcfs
```

## JSON sidecar

See `combine-gvcfs-rerun-20260928.json` for full stdout/stderr tails.
