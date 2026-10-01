# VariantRecalibrator rerun report — 2026-09-28

## Summary
- Total scripts: 10
- Passed: 10    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 161.204s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_recalibrator.py` | 0 | 6.428s | 7.95s | 7.0 MiB | 382.7 MiB | 6.31s | ✅ |
| `verify_variant_recalibrator_annotation_order_gatk_oracle.py` | 0 | 16.440s | 24.63s | 0.0 MiB | 396.6 MiB | 16.40s | ✅ |
| `verify_variant_recalibrator_attempt_iteration_gatk_oracle.py` | 0 | 18.544s | 73.19s | 3.1 MiB | 380.3 MiB | 18.49s | ✅ |
| `verify_variant_recalibrator_culprit_gatk_oracle.py` | 0 | 43.560s | 151.62s | 6.9 MiB | 348.1 MiB | 39.61s | ✅ |
| `verify_variant_recalibrator_gatk_model_oracle.py` | 0 | 16.988s | 42.88s | 4.3 MiB | 357.6 MiB | 16.95s | ✅ |
| `verify_variant_recalibrator_sample_every_gatk_oracle.py` | 0 | 20.883s | 93.08s | 2.2 MiB | 392.0 MiB | 20.84s | ✅ |
| `verify_variant_recalibrator_vbem_gatk_oracle.py` | 0 | 21.068s | 111.28s | 2.2 MiB | 368.5 MiB | 17.18s | ✅ |
| `verify_variant_recalibrator_zero_variance_gatk_oracle.py` | 0 | 17.193s | 58.32s | 0.9 MiB | 361.9 MiB | 11.47s | ✅ |
| `verify_vqsr_scatter_joint.py` | 0 | 0.042s | 8.10s | 6.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_vqsr_scatter_joint_4shard_oracle.py` | 0 | 0.058s | 8.39s | 7.1 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-recalibrator
```

## JSON sidecar

See `variant-recalibrator-rerun-20260928.json` for full stdout/stderr tails.
