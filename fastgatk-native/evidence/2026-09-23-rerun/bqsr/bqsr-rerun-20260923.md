# BaseRecalibrator rerun report — 2026-09-23

## Summary
- Total scripts: 9
- Passed: 9    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 576.811s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_bqsr.py` | 0 | 229.198s | 542.40s | 26.3 MiB | 505.8 MiB | 215.69s | ✅ |
| `verify_bqsr_context_size_gatk_oracle.py` | 0 | 42.085s | 319.96s | 8.8 MiB | 381.0 MiB | 41.89s | ✅ |
| `verify_bqsr_cram_gatk_oracle.py` | 0 | 38.557s | 268.41s | 15.4 MiB | 479.1 MiB | 35.91s | ✅ |
| `verify_bqsr_gatk_oracle.py` | 0 | 58.617s | 374.62s | 11.0 MiB | 469.1 MiB | 54.59s | ✅ |
| `verify_bqsr_indel_gatk_oracle.py` | 0 | 36.307s | 241.45s | 33.7 MiB | 495.2 MiB | 32.24s | ✅ |
| `verify_bqsr_long_read_gatk_oracle.py` | 0 | 22.141s | 214.73s | 10.1 MiB | 532.6 MiB | 21.39s | ✅ |
| `verify_bqsr_preserve_gatk_oracle.py` | 0 | 41.903s | 297.01s | 9.2 MiB | 505.3 MiB | 41.64s | ✅ |
| `verify_bqsr_read_filter_gatk_oracle.py` | 0 | 80.695s | 409.24s | 7.6 MiB | 402.9 MiB | 80.56s | ✅ |
| `verify_bqsr_report_roundtrip_gatk_oracle.py` | 0 | 27.308s | 335.67s | 10.3 MiB | 370.1 MiB | 26.17s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool bqsr
```

## JSON sidecar

See `bqsr-rerun-20260923.json` for full stdout/stderr tails.
