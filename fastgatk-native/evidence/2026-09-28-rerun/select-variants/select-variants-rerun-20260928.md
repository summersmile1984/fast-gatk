# SelectVariants rerun report — 2026-09-28

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 309.704s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_select_variants.py` | 0 | 91.134s | 7583.67s | 7.5 MiB | 389.7 MiB | 90.60s | ✅ |
| `verify_select_variants_filtered_nocall_oracle.py` | 0 | 14.783s | 7504.36s | 6.8 MiB | 389.8 MiB | 14.70s | ✅ |
| `verify_select_variants_filtered_oracle.py` | 0 | 10.068s | 7476.07s | 5.7 MiB | 381.9 MiB | 10.00s | ✅ |
| `verify_select_variants_gatk_oracle.py` | 0 | 10.278s | 7488.56s | 6.8 MiB | 387.3 MiB | 10.20s | ✅ |
| `verify_select_variants_refonly_gatk_oracle.py` | 0 | 154.454s | 7673.84s | 7.3 MiB | 348.5 MiB | 154.31s | ✅ |
| `verify_select_variants_sites_only_gatk_oracle.py` | 0 | 28.987s | 7520.91s | 5.6 MiB | 381.1 MiB | 28.94s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool select-variants
```

## JSON sidecar

See `select-variants-rerun-20260928.json` for full stdout/stderr tails.
