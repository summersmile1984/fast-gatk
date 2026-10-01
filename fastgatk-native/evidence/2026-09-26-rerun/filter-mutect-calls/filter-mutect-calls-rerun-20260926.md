# FilterMutectCalls rerun report — 2026-09-26

## Summary
- Total scripts: 13
- Passed: 13    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 975.701s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_filter_mutect_alignment_artifacts_gatk_oracle.py` | 0 | 199.220s | 427.13s | 696.5 MiB | 332.2 MiB | 5.20s | ✅ |
| `verify_filter_mutect_artifact_read_filter_gatk_oracle.py` | 0 | 200.212s | 563.29s | 696.8 MiB | 295.7 MiB | 6.05s | ✅ |
| `verify_filter_mutect_calls.py` | 0 | 13.511s | 134.47s | 8.4 MiB | 382.1 MiB | 7.21s | ✅ |
| `verify_filter_mutect_calls_contamination_oracle.py` | 0 | 23.538s | 179.97s | 7.9 MiB | 417.0 MiB | 22.60s | ✅ |
| `verify_filter_mutect_calls_germline_oracle.py` | 0 | 12.020s | 26.33s | 7.8 MiB | 390.0 MiB | 10.51s | ✅ |
| `verify_filter_mutect_contamination_joint_oracle.py` | 0 | 13.063s | 55.71s | 7.9 MiB | 384.3 MiB | 11.44s | ✅ |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | 0 | 217.632s | 853.13s | 697.5 MiB | 1364.6 MiB | 69.96s | ✅ |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | 0 | 51.326s | 287.28s | 8.7 MiB | 450.5 MiB | 50.66s | ✅ |
| `verify_filter_mutect_normal_artifact_oracle.py` | 0 | 10.515s | 146.02s | 7.9 MiB | 410.2 MiB | 5.37s | ✅ |
| `verify_filter_mutect_numeric_byte_equality_oracle.py` | 0 | 193.382s | 699.01s | 8.8 MiB | 333.5 MiB | 5.88s | ✅ |
| `verify_filter_mutect_orientation_gatk_oracle.py` | 0 | 10.316s | 193.71s | 8.0 MiB | 428.8 MiB | 10.03s | ✅ |
| `verify_filter_mutect_orientation_joint_oracle.py` | 0 | 9.883s | 205.97s | 8.0 MiB | 404.3 MiB | 9.68s | ✅ |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | 0 | 21.083s | 233.15s | 7.8 MiB | 431.9 MiB | 20.40s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool filter-mutect-calls
```

## JSON sidecar

See `filter-mutect-calls-rerun-20260926.json` for full stdout/stderr tails.
