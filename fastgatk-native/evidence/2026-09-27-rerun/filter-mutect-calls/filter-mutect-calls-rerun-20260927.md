# FilterMutectCalls rerun report — 2026-09-27

## Summary
- Total scripts: 13
- Passed: 13    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 864.942s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_filter_mutect_alignment_artifacts_gatk_oracle.py` | 0 | 183.407s | 431.52s | 697.4 MiB | 340.5 MiB | 4.53s | ✅ |
| `verify_filter_mutect_artifact_read_filter_gatk_oracle.py` | 0 | 179.990s | 303.40s | 696.4 MiB | 353.8 MiB | 4.46s | ✅ |
| `verify_filter_mutect_calls.py` | 0 | 7.172s | 20.68s | 8.3 MiB | 402.0 MiB | 5.80s | ✅ |
| `verify_filter_mutect_calls_contamination_oracle.py` | 0 | 19.397s | 113.20s | 7.9 MiB | 421.0 MiB | 18.41s | ✅ |
| `verify_filter_mutect_calls_germline_oracle.py` | 0 | 9.850s | 40.13s | 7.7 MiB | 379.4 MiB | 8.90s | ✅ |
| `verify_filter_mutect_contamination_joint_oracle.py` | 0 | 10.047s | 64.66s | 7.8 MiB | 395.9 MiB | 8.84s | ✅ |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | 0 | 195.867s | 691.49s | 696.2 MiB | 1040.0 MiB | 60.15s | ✅ |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | 0 | 42.950s | 181.12s | 8.7 MiB | 476.9 MiB | 42.76s | ✅ |
| `verify_filter_mutect_normal_artifact_oracle.py` | 0 | 8.427s | 74.89s | 7.9 MiB | 389.4 MiB | 8.18s | ✅ |
| `verify_filter_mutect_numeric_byte_equality_oracle.py` | 0 | 176.962s | 552.28s | 696.3 MiB | 345.5 MiB | 4.53s | ✅ |
| `verify_filter_mutect_orientation_gatk_oracle.py` | 0 | 8.362s | 83.69s | 7.9 MiB | 417.4 MiB | 8.29s | ✅ |
| `verify_filter_mutect_orientation_joint_oracle.py` | 0 | 8.643s | 124.08s | 7.9 MiB | 433.6 MiB | 8.39s | ✅ |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | 0 | 13.868s | 140.25s | 7.7 MiB | 394.1 MiB | 13.70s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool filter-mutect-calls
```

## JSON sidecar

See `filter-mutect-calls-rerun-20260927.json` for full stdout/stderr tails.
