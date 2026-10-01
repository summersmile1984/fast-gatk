# CollectAllelicCounts rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 46.938s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_collect_allelic_counts.py` | 0 | 0.060s | 611.84s | 4.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_collect_allelic_counts_gatk_oracle.py` | 0 | 46.878s | 630.88s | 5.6 MiB | 369.3 MiB | 46.82s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool collect-allelic-counts
```

## JSON sidecar

See `collect-allelic-counts-rerun-20260923.json` for full stdout/stderr tails.
