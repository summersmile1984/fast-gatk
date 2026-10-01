# ValidateVariants rerun report — 2026-09-28

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 111.394s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_validate_variants.py` | 0 | 14.245s | 7738.03s | 5.3 MiB | 315.3 MiB | 14.18s | ✅ |
| `verify_validate_variants_gatk_oracle.py` | 0 | 79.044s | 7797.65s | 5.1 MiB | 336.5 MiB | 78.99s | ✅ |
| `verify_validate_variants_symbolic_gatk_oracle.py` | 0 | 18.105s | 7751.60s | 0.0 MiB | 322.4 MiB | 18.08s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool validate-variants
```

## JSON sidecar

See `validate-variants-rerun-20260928.json` for full stdout/stderr tails.
