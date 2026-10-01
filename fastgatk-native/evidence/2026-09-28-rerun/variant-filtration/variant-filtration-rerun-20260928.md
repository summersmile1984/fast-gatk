# VariantFiltration rerun report — 2026-09-28

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 324.97s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_filtration.py` | 0 | 126.711s | 202.45s | 6.2 MiB | 396.2 MiB | 126.52s | ✅ |
| `verify_variant_filtration_asfilterstatus_gatk_oracle.py` | 0 | 43.429s | 105.61s | 4.4 MiB | 350.5 MiB | 43.37s | ✅ |
| `verify_variant_filtration_flag_only_gatk_oracle.py` | 0 | 37.566s | 49.77s | 0.0 MiB | 322.3 MiB | 33.07s | ✅ |
| `verify_variant_filtration_gatk_oracle.py` | 0 | 38.431s | 82.25s | 5.9 MiB | 394.9 MiB | 38.38s | ✅ |
| `verify_variant_filtration_missing_boolean_gatk_oracle.py` | 0 | 43.885s | 132.95s | 3.7 MiB | 389.1 MiB | 43.84s | ✅ |
| `verify_variant_filtration_set_nocall_gatk_oracle.py` | 0 | 34.948s | 24.66s | 5.6 MiB | 386.3 MiB | 34.90s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-filtration
```

## JSON sidecar

See `variant-filtration-rerun-20260928.json` for full stdout/stderr tails.
