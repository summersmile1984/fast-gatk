# CollectReadCounts rerun report — 2026-09-23

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 68.002s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_collect_read_counts.py` | 0 | 0.114s | 639.88s | 16.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_collect_read_counts_gatk_oracle.py` | 0 | 48.932s | 672.04s | 13.9 MiB | 367.6 MiB | 39.49s | ✅ |
| `verify_hdf5_simple_count_collection.py` | 0 | 18.956s | 649.65s | 13.2 MiB | 376.8 MiB | 18.90s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool collect-read-counts
```

## JSON sidecar

See `collect-read-counts-rerun-20260923.json` for full stdout/stderr tails.
