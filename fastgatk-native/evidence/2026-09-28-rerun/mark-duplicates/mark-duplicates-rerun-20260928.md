# MarkDuplicates rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 53.373s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_mark_duplicates.py` | 0 | 7.951s | 5617.66s | 5.5 MiB | 537.5 MiB | 7.89s | ✅ |
| `verify_mark_duplicates_gatk_oracle.py` | 0 | 37.245s | 5655.92s | 0.0 MiB | 7090.6 MiB | 37.20s | ✅ |
| `verify_mark_duplicates_pair_key_gatk_oracle.py` | 0 | 4.076s | 5605.79s | 0.0 MiB | 529.9 MiB | 4.04s | ✅ |
| `verify_mark_duplicates_tagging_policy_gatk_oracle.py` | 0 | 4.101s | 5610.15s | 2.0 MiB | 521.0 MiB | 4.06s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool mark-duplicates
```

## JSON sidecar

See `mark-duplicates-rerun-20260928.json` for full stdout/stderr tails.
