# VariantEval rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 86.756s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_eval.py` | 0 | 28.860s | 7863.89s | 6.6 MiB | 463.3 MiB | 28.70s | ✅ |
| `verify_variant_eval_gatk_oracle.py` | 0 | 37.745s | 7898.29s | 5.7 MiB | 441.3 MiB | 37.66s | ✅ |
| `verify_variant_eval_keep_ac0_gatk_oracle.py` | 0 | 10.091s | 7835.70s | 5.7 MiB | 452.1 MiB | 9.94s | ✅ |
| `verify_variant_eval_validation_report_gatk_oracle.py` | 0 | 10.060s | 7820.82s | 6.6 MiB | 436.1 MiB | 9.93s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-eval
```

## JSON sidecar

See `variant-eval-rerun-20260928.json` for full stdout/stderr tails.
