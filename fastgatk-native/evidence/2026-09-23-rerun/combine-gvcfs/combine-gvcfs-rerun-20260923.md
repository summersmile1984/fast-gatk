# CombineGVCFs rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 118.562s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_combine_gvcfs.py` | 0 | 9.050s | 679.34s | 7.4 MiB | 377.3 MiB | 8.85s | ✅ |
| `verify_combine_gvcfs_gatk_oracle.py` | 0 | 45.330s | 723.14s | 7.0 MiB | 393.6 MiB | 45.26s | ✅ |
| `verify_combine_gvcfs_interval_refblock_gatk_oracle.py` | 0 | 33.741s | 704.94s | 2.7 MiB | 395.3 MiB | 33.69s | ✅ |
| `verify_combine_gvcfs_plless_oracle.py` | 0 | 30.441s | 690.81s | 2.6 MiB | 376.6 MiB | 30.39s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool combine-gvcfs
```

## JSON sidecar

See `combine-gvcfs-rerun-20260923.json` for full stdout/stderr tails.
