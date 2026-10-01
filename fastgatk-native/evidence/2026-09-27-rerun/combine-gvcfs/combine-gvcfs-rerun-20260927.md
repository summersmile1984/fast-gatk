# CombineGVCFs rerun report — 2026-09-27

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 93.086s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_combine_gvcfs.py` | 0 | 11.314s | 8.81s | 6.9 MiB | 368.5 MiB | 11.12s | ✅ |
| `verify_combine_gvcfs_gatk_oracle.py` | 0 | 24.667s | 31.75s | 5.5 MiB | 392.5 MiB | 24.59s | ✅ |
| `verify_combine_gvcfs_interval_refblock_gatk_oracle.py` | 0 | 28.374s | 51.74s | 5.5 MiB | 387.1 MiB | 28.24s | ✅ |
| `verify_combine_gvcfs_plless_oracle.py` | 0 | 28.731s | 66.86s | 2.1 MiB | 384.9 MiB | 28.69s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool combine-gvcfs
```

## JSON sidecar

See `combine-gvcfs-rerun-20260927.json` for full stdout/stderr tails.
