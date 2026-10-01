# ModelSegments rerun report — 2026-09-28

## Summary
- Total scripts: 10
- Passed: 10    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 67.903s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_model_segments.py` | 0 | 4.903s | 5704.71s | 6.0 MiB | 369.7 MiB | 4.66s | ✅ |
| `verify_model_segments_allele_fraction_initialization_gatk_oracle.py` | 0 | 4.870s | 5688.45s | 5.7 MiB | 358.1 MiB | 4.79s | ✅ |
| `verify_model_segments_allele_fraction_likelihood_gatk_oracle.py` | 0 | 4.683s | 5673.06s | 6.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_model_segments_copy_ratio_conditionals_gatk_oracle.py` | 0 | 4.834s | 5680.83s | 6.2 MiB | 361.6 MiB | 4.75s | ✅ |
| `verify_model_segments_default_kernel_gatk_oracle.py` | 0 | 17.412s | 5735.80s | 5.7 MiB | 357.5 MiB | 17.33s | ✅ |
| `verify_model_segments_first_alt_fraction_gatk_oracle.py` | 0 | 4.662s | 5666.18s | 5.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_model_segments_input_segments_gatk_oracle.py` | 0 | 4.885s | 5695.40s | 5.7 MiB | 351.1 MiB | 4.82s | ✅ |
| `verify_model_segments_mcmc_oracle.py` | 0 | 0.493s | 5659.12s | 17.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_model_segments_multisample_gatk_oracle.py` | 0 | 8.448s | 5713.44s | 0.0 MiB | 339.6 MiB | 8.39s | ✅ |
| `verify_model_segments_smoothing_gatk_oracle.py` | 0 | 12.713s | 5720.98s | 5.7 MiB | 346.3 MiB | 12.66s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool model-segments
```

## JSON sidecar

See `model-segments-rerun-20260928.json` for full stdout/stderr tails.
