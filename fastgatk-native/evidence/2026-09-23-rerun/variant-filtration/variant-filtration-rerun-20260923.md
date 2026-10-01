# VariantFiltration rerun report — 2026-09-23

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 377.549s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_filtration.py` | 0 | 159.847s | 7209.10s | 6.2 MiB | 391.0 MiB | 159.62s | ✅ |
| `verify_variant_filtration_asfilterstatus_gatk_oracle.py` | 0 | 46.415s | 7101.48s | 5.2 MiB | 341.0 MiB | 46.35s | ✅ |
| `verify_variant_filtration_flag_only_gatk_oracle.py` | 0 | 40.061s | 7050.76s | 6.1 MiB | 359.2 MiB | 40.00s | ✅ |
| `verify_variant_filtration_gatk_oracle.py` | 0 | 55.632s | 7135.14s | 5.9 MiB | 385.3 MiB | 55.57s | ✅ |
| `verify_variant_filtration_missing_boolean_gatk_oracle.py` | 0 | 40.829s | 7080.32s | 2.1 MiB | 390.4 MiB | 40.78s | ✅ |
| `verify_variant_filtration_set_nocall_gatk_oracle.py` | 0 | 34.765s | 7023.44s | 4.2 MiB | 387.5 MiB | 34.72s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-filtration
```

## JSON sidecar

See `variant-filtration-rerun-20260923.json` for full stdout/stderr tails.
