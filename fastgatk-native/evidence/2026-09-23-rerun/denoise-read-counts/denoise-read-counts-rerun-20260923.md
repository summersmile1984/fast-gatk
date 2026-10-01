# DenoiseReadCounts rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 61.985s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_denoise_read_counts.py` | 0 | 0.162s | 856.49s | 17.8 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py` | 0 | 5.148s | 860.62s | 15.7 MiB | 345.4 MiB | 5.06s | ✅ |
| `verify_denoise_read_counts_integer_input_gatk_oracle.py` | 0 | 28.769s | 892.12s | 12.5 MiB | 335.4 MiB | 28.67s | ✅ |
| `verify_denoise_read_counts_interval_identity_gatk_oracle.py` | 0 | 27.906s | 874.21s | 16.1 MiB | 337.4 MiB | 21.94s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool denoise-read-counts
```

## JSON sidecar

See `denoise-read-counts-rerun-20260923.json` for full stdout/stderr tails.
