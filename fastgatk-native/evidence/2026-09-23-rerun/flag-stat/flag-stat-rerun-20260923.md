# FlagStat rerun report — 2026-09-23

## Summary
- Total scripts: 1
- Passed: 1    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 111.86s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_flag_stat.py` | 0 | 111.860s | 1300.59s | 6.6 MiB | 361.8 MiB | 111.65s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool flag-stat
```

## JSON sidecar

See `flag-stat-rerun-20260923.json` for full stdout/stderr tails.
