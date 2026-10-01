# FilterMutectCalls rerun report — 2026-09-23

## Summary
- Total scripts: 10
- Passed: 10    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 206.064s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_filter_mutect_calls.py` | 0 | 9.706s | 1022.49s | 8.5 MiB | 390.9 MiB | 8.64s | ✅ |
| `verify_filter_mutect_calls_contamination_oracle.py` | 0 | 35.644s | 1168.64s | 7.9 MiB | 426.5 MiB | 35.57s | ✅ |
| `verify_filter_mutect_calls_germline_oracle.py` | 0 | 13.600s | 1063.53s | 7.7 MiB | 402.3 MiB | 13.43s | ✅ |
| `verify_filter_mutect_contamination_joint_oracle.py` | 0 | 18.720s | 1075.23s | 7.8 MiB | 392.8 MiB | 18.63s | ✅ |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | 0 | 43.249s | 1237.21s | 727.9 MiB | 2020.3 MiB | 31.00s | ✅ |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | 0 | 27.678s | 1129.63s | 8.7 MiB | 1222.7 MiB | 27.45s | ✅ |
| `verify_filter_mutect_normal_artifact_oracle.py` | 0 | 13.218s | 1035.04s | 0.0 MiB | 397.7 MiB | 13.13s | ✅ |
| `verify_filter_mutect_orientation_gatk_oracle.py` | 0 | 13.263s | 1048.42s | 8.0 MiB | 428.1 MiB | 13.18s | ✅ |
| `verify_filter_mutect_orientation_joint_oracle.py` | 0 | 11.867s | 1083.54s | 8.0 MiB | 411.2 MiB | 11.81s | ✅ |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | 0 | 19.119s | 1148.53s | 7.9 MiB | 369.3 MiB | 18.89s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool filter-mutect-calls
```

## JSON sidecar

See `filter-mutect-calls-rerun-20260923.json` for full stdout/stderr tails.
