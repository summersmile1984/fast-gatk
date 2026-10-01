# ApplyBQSR rerun report — 2026-09-25

## Summary
- Total scripts: 1
- Passed: 1    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 58.138s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_bqsr_alias_gatk_oracle.py` | 0 | 58.138s | 99.44s | 8.2 MiB | 442.8 MiB | 52.64s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-bqsr
```

## JSON sidecar

See `apply-bqsr-rerun-20260925.json` for full stdout/stderr tails.
