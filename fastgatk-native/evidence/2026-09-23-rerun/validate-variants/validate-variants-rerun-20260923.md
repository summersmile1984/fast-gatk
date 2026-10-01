# ValidateVariants rerun report — 2026-09-23

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 141.143s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_validate_variants.py` | 0 | 20.487s | 6855.19s | 5.2 MiB | 308.2 MiB | 20.43s | ✅ |
| `verify_validate_variants_gatk_oracle.py` | 0 | 92.593s | 6915.68s | 5.2 MiB | 312.5 MiB | 92.54s | ✅ |
| `verify_validate_variants_symbolic_gatk_oracle.py` | 0 | 28.063s | 6868.29s | 4.7 MiB | 292.2 MiB | 28.03s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool validate-variants
```

## JSON sidecar

See `validate-variants-rerun-20260923.json` for full stdout/stderr tails.
