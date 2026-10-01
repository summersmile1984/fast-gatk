# VariantRecalibrator rerun report — 2026-09-23

## Summary
- Total scripts: 10
- Passed: 10    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 238.426s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_recalibrator.py` | 0 | 8.965s | 7218.09s | 6.7 MiB | 389.0 MiB | 8.83s | ✅ |
| `verify_variant_recalibrator_annotation_order_gatk_oracle.py` | 0 | 19.842s | 7253.77s | 6.8 MiB | 366.8 MiB | 13.17s | ✅ |
| `verify_variant_recalibrator_attempt_iteration_gatk_oracle.py` | 0 | 20.860s | 7273.60s | 4.4 MiB | 385.4 MiB | 20.81s | ✅ |
| `verify_variant_recalibrator_culprit_gatk_oracle.py` | 0 | 75.453s | 7380.57s | 6.9 MiB | 375.4 MiB | 75.37s | ✅ |
| `verify_variant_recalibrator_gatk_model_oracle.py` | 0 | 30.790s | 7294.40s | 0.0 MiB | 362.5 MiB | 30.75s | ✅ |
| `verify_variant_recalibrator_sample_every_gatk_oracle.py` | 0 | 31.990s | 7337.34s | 6.0 MiB | 368.5 MiB | 31.94s | ✅ |
| `verify_variant_recalibrator_vbem_gatk_oracle.py` | 0 | 31.725s | 7315.50s | 6.8 MiB | 379.1 MiB | 31.68s | ✅ |
| `verify_variant_recalibrator_zero_variance_gatk_oracle.py` | 0 | 18.678s | 7235.85s | 0.0 MiB | 376.1 MiB | 18.65s | ✅ |
| `verify_vqsr_scatter_joint.py` | 0 | 0.057s | 7218.49s | 6.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_vqsr_scatter_joint_4shard_oracle.py` | 0 | 0.066s | 7218.80s | 6.7 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-recalibrator
```

## JSON sidecar

See `variant-recalibrator-rerun-20260923.json` for full stdout/stderr tails.
