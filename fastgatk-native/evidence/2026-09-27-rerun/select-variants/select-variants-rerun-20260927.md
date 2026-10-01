# SelectVariants rerun report — 2026-09-27

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 384.689s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_select_variants.py` | 0 | 112.710s | 161.06s | 7.5 MiB | 385.2 MiB | 107.07s | ✅ |
| `verify_select_variants_filtered_nocall_oracle.py` | 0 | 16.487s | 52.26s | 5.7 MiB | 376.7 MiB | 11.21s | ✅ |
| `verify_select_variants_filtered_oracle.py` | 0 | 11.119s | 30.50s | 7.1 MiB | 379.3 MiB | 10.92s | ✅ |
| `verify_select_variants_gatk_oracle.py` | 0 | 11.098s | 15.13s | 6.9 MiB | 386.5 MiB | 5.11s | ✅ |
| `verify_select_variants_refonly_gatk_oracle.py` | 0 | 215.894s | 271.31s | 7.3 MiB | 378.0 MiB | 215.55s | ✅ |
| `verify_select_variants_sites_only_gatk_oracle.py` | 0 | 17.381s | 76.94s | 5.6 MiB | 401.4 MiB | 17.29s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool select-variants
```

## JSON sidecar

See `select-variants-rerun-20260927.json` for full stdout/stderr tails.
