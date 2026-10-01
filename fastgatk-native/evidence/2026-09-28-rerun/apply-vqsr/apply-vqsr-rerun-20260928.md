# ApplyVQSR rerun report — 2026-09-28

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 187.317s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_vqsr.py` | 0 | 35.550s | 179.65s | 5.2 MiB | 370.4 MiB | 35.47s | ✅ |
| `verify_apply_vqsr_default_cutoff_gatk_oracle.py` | 0 | 14.040s | 98.17s | 0.8 MiB | 384.2 MiB | 14.00s | ✅ |
| `verify_apply_vqsr_exclude_intervals_gatk_oracle.py` | 0 | 35.486s | 154.84s | 4.8 MiB | 373.0 MiB | 35.45s | ✅ |
| `verify_apply_vqsr_filter_booleans_gatk_oracle.py` | 0 | 69.301s | 219.13s | 5.4 MiB | 372.0 MiB | 69.25s | ✅ |
| `verify_apply_vqsr_gatk_oracle.py` | 0 | 14.500s | 115.12s | 0.0 MiB | 389.3 MiB | 14.46s | ✅ |
| `verify_apply_vqsr_sites_only_gatk_oracle.py` | 0 | 18.440s | 131.58s | 5.0 MiB | 374.0 MiB | 18.40s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-vqsr
```

## JSON sidecar

See `apply-vqsr-rerun-20260928.json` for full stdout/stderr tails.
