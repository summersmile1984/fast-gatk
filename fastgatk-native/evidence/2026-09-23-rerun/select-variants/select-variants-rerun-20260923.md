# SelectVariants rerun report — 2026-09-23

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 423.558s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_select_variants.py` | 0 | 132.562s | 6692.00s | 7.8 MiB | 383.7 MiB | 131.88s | ✅ |
| `verify_select_variants_filtered_nocall_oracle.py` | 0 | 27.655s | 6606.63s | 5.9 MiB | 371.9 MiB | 27.59s | ✅ |
| `verify_select_variants_filtered_oracle.py` | 0 | 14.466s | 6580.50s | 6.7 MiB | 384.4 MiB | 14.40s | ✅ |
| `verify_select_variants_gatk_oracle.py` | 0 | 14.671s | 6592.41s | 7.0 MiB | 395.7 MiB | 14.60s | ✅ |
| `verify_select_variants_refonly_gatk_oracle.py` | 0 | 206.071s | 6788.14s | 7.4 MiB | 354.2 MiB | 205.89s | ✅ |
| `verify_select_variants_sites_only_gatk_oracle.py` | 0 | 28.133s | 6624.82s | 5.7 MiB | 397.2 MiB | 28.07s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool select-variants
```

## JSON sidecar

See `select-variants-rerun-20260923.json` for full stdout/stderr tails.
