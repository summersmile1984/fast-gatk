# CollectAllelicCounts rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 38.462s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_collect_allelic_counts.py` | 0 | 0.036s | 679.57s | 3.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_collect_allelic_counts_gatk_oracle.py` | 0 | 38.426s | 698.65s | 4.3 MiB | 376.9 MiB | 38.37s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool collect-allelic-counts
```

## JSON sidecar

See `collect-allelic-counts-rerun-20260928.json` for full stdout/stderr tails.
