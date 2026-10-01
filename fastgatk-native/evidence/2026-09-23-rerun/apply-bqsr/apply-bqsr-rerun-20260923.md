# ApplyBQSR rerun report — 2026-09-23

## Summary
- Total scripts: 1
- Passed: 1    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 50.253s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_bqsr_alias_gatk_oracle.py` | 0 | 50.253s | 77.32s | 10.3 MiB | 433.0 MiB | 48.29s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-bqsr
```

## JSON sidecar

See `apply-bqsr-rerun-20260923.json` for full stdout/stderr tails.
