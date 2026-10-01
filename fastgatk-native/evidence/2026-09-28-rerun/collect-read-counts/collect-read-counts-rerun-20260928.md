# CollectReadCounts rerun report — 2026-09-28

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 45.712s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_collect_read_counts.py` | 0 | 0.093s | 706.84s | 13.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_collect_read_counts_gatk_oracle.py` | 0 | 33.307s | 739.75s | 13.1 MiB | 366.1 MiB | 33.22s | ✅ |
| `verify_hdf5_simple_count_collection.py` | 0 | 12.312s | 716.78s | 17.8 MiB | 376.6 MiB | 12.27s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool collect-read-counts
```

## JSON sidecar

See `collect-read-counts-rerun-20260928.json` for full stdout/stderr tails.
