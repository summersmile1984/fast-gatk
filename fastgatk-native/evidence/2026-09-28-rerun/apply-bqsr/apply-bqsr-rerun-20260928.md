# ApplyBQSR rerun report — 2026-09-28

## Summary
- Total scripts: 1
- Passed: 1    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 46.614s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_apply_bqsr_alias_gatk_oracle.py` | 0 | 46.614s | 81.41s | 8.4 MiB | 451.2 MiB | 41.25s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool apply-bqsr
```

## JSON sidecar

See `apply-bqsr-rerun-20260928.json` for full stdout/stderr tails.
