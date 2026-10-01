# MarkDuplicates rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 69.859s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_mark_duplicates.py` | 0 | 12.490s | 4413.27s | 4.7 MiB | 525.8 MiB | 12.44s | ✅ |
| `verify_mark_duplicates_gatk_oracle.py` | 0 | 45.458s | 4452.57s | 5.1 MiB | 7089.0 MiB | 45.42s | ✅ |
| `verify_mark_duplicates_pair_key_gatk_oracle.py` | 0 | 5.662s | 4399.40s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mark_duplicates_tagging_policy_gatk_oracle.py` | 0 | 6.249s | 4404.55s | 2.3 MiB | 520.2 MiB | 6.21s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool mark-duplicates
```

## JSON sidecar

See `mark-duplicates-rerun-20260923.json` for full stdout/stderr tails.
