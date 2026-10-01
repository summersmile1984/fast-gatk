# CallCopyRatioSegments rerun report — 2026-09-28

## Summary
- Total scripts: 5
- Passed: 5    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 40.602s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_call_copy_ratio_segments.py` | 0 | 0.050s | 637.34s | 4.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py` | 0 | 4.488s | 642.09s | 5.5 MiB | 334.8 MiB | 4.44s | ✅ |
| `verify_call_copy_ratio_segments_gatk_oracle.py` | 0 | 11.804s | 664.04s | 5.5 MiB | 367.5 MiB | 11.74s | ✅ |
| `verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py` | 0 | 7.982s | 650.98s | 5.5 MiB | 342.3 MiB | 7.92s | ✅ |
| `verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py` | 0 | 16.278s | 676.43s | 5.5 MiB | 347.4 MiB | 16.22s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool call-copy-ratio-segments
```

## JSON sidecar

See `call-copy-ratio-segments-rerun-20260928.json` for full stdout/stderr tails.
