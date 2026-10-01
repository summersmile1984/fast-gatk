# CallCopyRatioSegments rerun report — 2026-09-23

## Summary
- Total scripts: 5
- Passed: 5    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 53.955s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_call_copy_ratio_segments.py` | 0 | 0.070s | 571.66s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py` | 0 | 5.269s | 576.82s | 5.6 MiB | 354.1 MiB | 5.19s | ✅ |
| `verify_call_copy_ratio_segments_gatk_oracle.py` | 0 | 21.205s | 608.58s | 5.5 MiB | 345.2 MiB | 21.11s | ✅ |
| `verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py` | 0 | 10.563s | 585.28s | 5.7 MiB | 346.9 MiB | 10.48s | ✅ |
| `verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py` | 0 | 16.848s | 597.07s | 5.7 MiB | 370.1 MiB | 16.76s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool call-copy-ratio-segments
```

## JSON sidecar

See `call-copy-ratio-segments-rerun-20260923.json` for full stdout/stderr tails.
