# FilterMutectCalls rerun report — 2026-09-25

## Summary
- Total scripts: 13
- Passed: 13    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 1030.599s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_filter_mutect_alignment_artifacts_gatk_oracle.py` | 0 | 201.254s | 540.31s | 697.4 MiB | 372.4 MiB | 8.19s | ✅ |
| `verify_filter_mutect_artifact_read_filter_gatk_oracle.py` | 0 | 197.339s | 396.94s | 698.0 MiB | 354.8 MiB | 8.32s | ✅ |
| `verify_filter_mutect_calls.py` | 0 | 12.817s | 64.90s | 8.3 MiB | 380.1 MiB | 8.91s | ✅ |
| `verify_filter_mutect_calls_contamination_oracle.py` | 0 | 36.667s | 134.41s | 8.2 MiB | 401.7 MiB | 36.37s | ✅ |
| `verify_filter_mutect_calls_germline_oracle.py` | 0 | 16.001s | 78.45s | 7.9 MiB | 394.9 MiB | 15.76s | ✅ |
| `verify_filter_mutect_contamination_joint_oracle.py` | 0 | 20.466s | 91.57s | 7.7 MiB | 407.4 MiB | 20.36s | ✅ |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | 0 | 220.872s | 836.95s | 696.3 MiB | 976.5 MiB | 71.00s | ✅ |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | 0 | 52.541s | 236.76s | 8.7 MiB | 516.7 MiB | 49.13s | ✅ |
| `verify_filter_mutect_normal_artifact_oracle.py` | 0 | 16.634s | 101.21s | 7.8 MiB | 383.7 MiB | 16.52s | ✅ |
| `verify_filter_mutect_numeric_byte_equality_oracle.py` | 0 | 197.206s | 682.14s | 697.3 MiB | 371.0 MiB | 8.41s | ✅ |
| `verify_filter_mutect_orientation_gatk_oracle.py` | 0 | 14.537s | 111.31s | 8.2 MiB | 417.1 MiB | 14.39s | ✅ |
| `verify_filter_mutect_orientation_joint_oracle.py` | 0 | 15.228s | 144.79s | 8.0 MiB | 426.9 MiB | 15.08s | ✅ |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | 0 | 29.037s | 260.83s | 7.9 MiB | 391.2 MiB | 28.38s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool filter-mutect-calls
```

## JSON sidecar

See `filter-mutect-calls-rerun-20260925.json` for full stdout/stderr tails.
