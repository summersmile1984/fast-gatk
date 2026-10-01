# VariantEval rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 127.712s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_eval.py` | 0 | 34.537s | 6964.88s | 6.7 MiB | 471.6 MiB | 34.37s | ✅ |
| `verify_variant_eval_gatk_oracle.py` | 0 | 59.588s | 6999.23s | 6.6 MiB | 446.3 MiB | 52.22s | ✅ |
| `verify_variant_eval_keep_ac0_gatk_oracle.py` | 0 | 18.349s | 6939.60s | 5.7 MiB | 446.2 MiB | 18.25s | ✅ |
| `verify_variant_eval_validation_report_gatk_oracle.py` | 0 | 15.238s | 6928.67s | 6.8 MiB | 456.9 MiB | 15.12s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-eval
```

## JSON sidecar

See `variant-eval-rerun-20260923.json` for full stdout/stderr tails.
