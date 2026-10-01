# DenoiseReadCounts rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 48.196s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_denoise_read_counts.py` | 0 | 0.149s | 923.21s | 16.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py` | 0 | 3.950s | 927.19s | 13.0 MiB | 342.6 MiB | 3.88s | ✅ |
| `verify_denoise_read_counts_integer_input_gatk_oracle.py` | 0 | 24.337s | 958.86s | 12.5 MiB | 336.5 MiB | 24.25s | ✅ |
| `verify_denoise_read_counts_interval_identity_gatk_oracle.py` | 0 | 19.760s | 940.97s | 17.0 MiB | 352.0 MiB | 19.65s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool denoise-read-counts
```

## JSON sidecar

See `denoise-read-counts-rerun-20260928.json` for full stdout/stderr tails.
