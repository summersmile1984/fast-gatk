# ApplyVQSR rerun report — 2026-09-23

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 268.403s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_vqsr.py` | 0 | 55.895s | 160.91s | 5.2 MiB | 380.1 MiB | 55.81s | ✅ |
| `verify_apply_vqsr_default_cutoff_gatk_oracle.py` | 0 | 28.638s | 118.25s | 0.0 MiB | 354.4 MiB | 28.60s | ✅ |
| `verify_apply_vqsr_exclude_intervals_gatk_oracle.py` | 0 | 45.762s | 136.74s | 2.1 MiB | 358.0 MiB | 27.23s | ✅ |
| `verify_apply_vqsr_filter_booleans_gatk_oracle.py` | 0 | 85.718s | 196.23s | 5.0 MiB | 359.2 MiB | 85.66s | ✅ |
| `verify_apply_vqsr_gatk_oracle.py` | 0 | 23.897s | 90.72s | 4.9 MiB | 364.5 MiB | 23.86s | ✅ |
| `verify_apply_vqsr_sites_only_gatk_oracle.py` | 0 | 28.493s | 104.56s | 0.0 MiB | 357.5 MiB | 28.46s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-vqsr
```

## JSON sidecar

See `apply-vqsr-rerun-20260923.json` for full stdout/stderr tails.
