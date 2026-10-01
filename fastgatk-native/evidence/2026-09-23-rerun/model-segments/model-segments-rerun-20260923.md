# ModelSegments rerun report — 2026-09-23

## Summary
- Total scripts: 10
- Passed: 10    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 95.668s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_model_segments.py` | 0 | 6.585s | 4478.00s | 6.1 MiB | 320.6 MiB | 6.35s | ✅ |
| `verify_model_segments_allele_fraction_initialization_gatk_oracle.py` | 0 | 7.578s | 4493.10s | 6.4 MiB | 339.1 MiB | 7.46s | ✅ |
| `verify_model_segments_allele_fraction_likelihood_gatk_oracle.py` | 0 | 7.847s | 4500.23s | 6.4 MiB | 349.9 MiB | 7.73s | ✅ |
| `verify_model_segments_copy_ratio_conditionals_gatk_oracle.py` | 0 | 6.358s | 4470.33s | 6.3 MiB | 340.8 MiB | 6.23s | ✅ |
| `verify_model_segments_default_kernel_gatk_oracle.py` | 0 | 24.460s | 4532.16s | 0.0 MiB | 347.7 MiB | 24.32s | ✅ |
| `verify_model_segments_first_alt_fraction_gatk_oracle.py` | 0 | 6.327s | 4462.85s | 0.0 MiB | 345.8 MiB | 6.30s | ✅ |
| `verify_model_segments_input_segments_gatk_oracle.py` | 0 | 7.553s | 4485.80s | 5.9 MiB | 377.6 MiB | 7.43s | ✅ |
| `verify_model_segments_mcmc_oracle.py` | 0 | 0.553s | 4455.53s | 17.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_model_segments_multisample_gatk_oracle.py` | 0 | 13.089s | 4510.33s | 5.6 MiB | 380.8 MiB | 13.03s | ✅ |
| `verify_model_segments_smoothing_gatk_oracle.py` | 0 | 15.318s | 4517.66s | 5.9 MiB | 340.3 MiB | 15.27s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool model-segments
```

## JSON sidecar

See `model-segments-rerun-20260923.json` for full stdout/stderr tails.
