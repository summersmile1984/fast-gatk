# ApplyVQSR rerun report — 2026-09-25

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 255.566s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_vqsr.py` | 0 | 57.741s | 215.01s | 5.6 MiB | 395.1 MiB | 57.63s | ✅ |
| `verify_apply_vqsr_default_cutoff_gatk_oracle.py` | 0 | 22.749s | 117.67s | 1.8 MiB | 393.1 MiB | 22.70s | ✅ |
| `verify_apply_vqsr_exclude_intervals_gatk_oracle.py` | 0 | 51.294s | 181.42s | 5.4 MiB | 408.4 MiB | 51.24s | ✅ |
| `verify_apply_vqsr_filter_booleans_gatk_oracle.py` | 0 | 70.648s | 262.22s | 5.4 MiB | 383.8 MiB | 70.60s | ✅ |
| `verify_apply_vqsr_gatk_oracle.py` | 0 | 26.481s | 135.87s | 0.8 MiB | 396.7 MiB | 26.44s | ✅ |
| `verify_apply_vqsr_sites_only_gatk_oracle.py` | 0 | 26.653s | 154.10s | 5.0 MiB | 352.9 MiB | 13.73s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-vqsr
```

## JSON sidecar

See `apply-vqsr-rerun-20260925.json` for full stdout/stderr tails.
