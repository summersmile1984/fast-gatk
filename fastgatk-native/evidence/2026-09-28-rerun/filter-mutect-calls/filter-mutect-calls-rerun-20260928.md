# FilterMutectCalls rerun report — 2026-09-28

## Summary
- Total scripts: 13
- Passed: 13    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 1192.591s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_filter_mutect_alignment_artifacts_gatk_oracle.py` | 0 | 254.899s | 1628.68s | 712.3 MiB | 291.9 MiB | 4.74s | ✅ |
| `verify_filter_mutect_artifact_read_filter_gatk_oracle.py` | 0 | 254.052s | 1440.62s | 711.5 MiB | 367.2 MiB | 4.67s | ✅ |
| `verify_filter_mutect_calls.py` | 0 | 8.294s | 1090.93s | 8.3 MiB | 375.1 MiB | 6.09s | ✅ |
| `verify_filter_mutect_calls_contamination_oracle.py` | 0 | 19.859s | 1175.19s | 7.9 MiB | 412.7 MiB | 14.14s | ✅ |
| `verify_filter_mutect_calls_germline_oracle.py` | 0 | 10.503s | 1115.00s | 7.9 MiB | 403.4 MiB | 9.35s | ✅ |
| `verify_filter_mutect_contamination_joint_oracle.py` | 0 | 10.647s | 1141.29s | 7.7 MiB | 395.6 MiB | 9.41s | ✅ |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | 0 | 269.719s | 2017.63s | 708.0 MiB | 1184.7 MiB | 59.88s | ✅ |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | 0 | 42.711s | 1229.66s | 8.6 MiB | 511.9 MiB | 42.22s | ✅ |
| `verify_filter_mutect_normal_artifact_oracle.py` | 0 | 13.347s | 1185.21s | 7.9 MiB | 393.4 MiB | 13.14s | ✅ |
| `verify_filter_mutect_numeric_byte_equality_oracle.py` | 0 | 251.743s | 1816.80s | 712.5 MiB | 294.3 MiB | 4.69s | ✅ |
| `verify_filter_mutect_orientation_gatk_oracle.py` | 0 | 9.105s | 1149.98s | 7.8 MiB | 400.6 MiB | 8.94s | ✅ |
| `verify_filter_mutect_orientation_joint_oracle.py` | 0 | 23.363s | 1238.72s | 7.8 MiB | 430.7 MiB | 23.28s | ✅ |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | 0 | 24.349s | 1253.40s | 7.8 MiB | 379.1 MiB | 24.25s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool filter-mutect-calls
```

## JSON sidecar

See `filter-mutect-calls-rerun-20260928.json` for full stdout/stderr tails.
