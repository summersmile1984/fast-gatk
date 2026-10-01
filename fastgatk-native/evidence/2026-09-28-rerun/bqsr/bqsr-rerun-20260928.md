# BaseRecalibrator rerun report — 2026-09-28

## Summary
- Total scripts: 9
- Passed: 9    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 432.52s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_bqsr.py` | 0 | 180.609s | 606.39s | 26.2 MiB | 469.0 MiB | 132.30s | ✅ |
| `verify_bqsr_context_size_gatk_oracle.py` | 0 | 28.480s | 293.04s | 8.7 MiB | 414.6 MiB | 21.80s | ✅ |
| `verify_bqsr_cram_gatk_oracle.py` | 0 | 34.542s | 356.45s | 13.2 MiB | 447.1 MiB | 18.43s | ✅ |
| `verify_bqsr_gatk_oracle.py` | 0 | 51.991s | 454.33s | 9.1 MiB | 455.1 MiB | 34.40s | ✅ |
| `verify_bqsr_indel_gatk_oracle.py` | 0 | 27.631s | 267.08s | 33.5 MiB | 472.2 MiB | 19.10s | ✅ |
| `verify_bqsr_long_read_gatk_oracle.py` | 0 | 20.620s | 236.29s | 10.1 MiB | 493.7 MiB | 19.14s | ✅ |
| `verify_bqsr_preserve_gatk_oracle.py` | 0 | 32.823s | 326.39s | 9.0 MiB | 475.5 MiB | 25.46s | ✅ |
| `verify_bqsr_read_filter_gatk_oracle.py` | 0 | 37.341s | 394.75s | 7.9 MiB | 410.3 MiB | 37.08s | ✅ |
| `verify_bqsr_report_roundtrip_gatk_oracle.py` | 0 | 18.483s | 413.93s | 9.0 MiB | 376.7 MiB | 15.38s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool bqsr
```

## JSON sidecar

See `bqsr-rerun-20260928.json` for full stdout/stderr tails.
