# GenotypeGVCFs rerun report — 2026-09-27

## Summary
- Total scripts: 25
- Passed: 25    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2091.244s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gatk_genotype_gvcf.py` | 0 | 36.692s | 179.85s | 12.2 MiB | 350.2 MiB | 26.41s | ✅ |
| `verify_gatk_genotype_gvcf_inbreeding.py` | 0 | 5.886s | 17.76s | 8.6 MiB | 282.4 MiB | 5.63s | ✅ |
| `verify_gatk_genotype_gvcf_legacy_qual.py` | 0 | 12.386s | 77.89s | 8.7 MiB | 363.2 MiB | 5.77s | ✅ |
| `verify_gatk_genotype_gvcf_multiallelic.py` | 0 | 5.802s | 8.99s | 8.4 MiB | 297.3 MiB | 5.57s | ✅ |
| `verify_gatk_genotype_gvcf_multisample.py` | 0 | 22.587s | 134.69s | 8.6 MiB | 331.2 MiB | 22.43s | ✅ |
| `verify_genotype_gvcf.py` | 0 | 8.173s | 56.63s | 13.2 MiB | 429.5 MiB | 5.40s | ✅ |
| `verify_genotype_gvcf_assignment_gatk_oracle.py` | 0 | 15.680s | 99.44s | 8.7 MiB | 424.7 MiB | 15.28s | ✅ |
| `verify_genotype_gvcf_confident_hom_alt_gatk_oracle.py` | 0 | 60.516s | 244.70s | 8.3 MiB | 342.6 MiB | 55.79s | ✅ |
| `verify_genotype_gvcf_dense_materialize_gatk_oracle.py` | 0 | 187.602s | 585.35s | 8.3 MiB | 331.0 MiB | 182.36s | ✅ |
| `verify_genotype_gvcf_depth_gate_gatk_oracle.py` | 0 | 81.927s | 335.26s | 8.2 MiB | 306.3 MiB | 75.36s | ✅ |
| `verify_genotype_gvcf_exclude_intervals_gatk_oracle.py` | 0 | 9.828s | 109.76s | 9.0 MiB | 294.4 MiB | 9.63s | ✅ |
| `verify_genotype_gvcf_gp_input_gatk_oracle.py` | 0 | 25.439s | 196.76s | 8.5 MiB | 352.1 MiB | 25.28s | ✅ |
| `verify_genotype_gvcf_header_order_gatk_oracle.py` | 0 | 71.186s | 290.07s | 8.3 MiB | 329.5 MiB | 70.99s | ✅ |
| `verify_genotype_gvcf_include_non_variant_gatk_oracle.py` | 0 | 7.142s | 143.14s | 13.2 MiB | 327.2 MiB | 5.12s | ✅ |
| `verify_genotype_gvcf_independent_corpus_gatk_oracle.py` | 0 | 236.014s | 1539.01s | 264.7 MiB | 878.3 MiB | 121.77s | ✅ |
| `verify_genotype_gvcf_lowqual_gatk_oracle.py` | 0 | 215.846s | 823.46s | 8.4 MiB | 328.7 MiB | 197.95s | ✅ |
| `verify_genotype_gvcf_malformed_gatk_oracle.py` | 0 | 218.768s | 923.29s | 8.2 MiB | 341.0 MiB | 212.04s | ✅ |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` | 0 | 23.994s | 254.76s | 8.6 MiB | 403.6 MiB | 23.80s | ✅ |
| `verify_genotype_gvcf_multisample_reference_confidence_oracle.py` | 0 | 38.127s | 357.28s | 9.3 MiB | 330.2 MiB | 36.72s | ✅ |
| `verify_genotype_gvcf_reverse_trim_gatk_oracle.py` | 0 | 126.356s | 471.73s | 8.3 MiB | 301.8 MiB | 105.98s | ✅ |
| `verify_genotype_gvcf_spandel_gatk_oracle.py` | 0 | 451.235s | 1827.60s | 8.4 MiB | 341.7 MiB | 429.59s | ✅ |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py` | 0 | 24.402s | 375.40s | 8.5 MiB | 359.7 MiB | 24.17s | ✅ |
| `verify_genotype_gvcf_spanning_source_gq_gatk_oracle.py` | 0 | 113.931s | 725.04s | 0.0 MiB | 868.6 MiB | 12.36s | ✅ |
| `verify_genotype_gvcf_star_only_locus_gatk_oracle.py` | 0 | 75.759s | 412.44s | 8.4 MiB | 314.0 MiB | 62.22s | ✅ |
| `verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py` | 0 | 15.966s | 596.65s | 8.6 MiB | 273.2 MiB | 15.82s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool genotype-gvcf
```

## JSON sidecar

See `genotype-gvcf-rerun-20260927.json` for full stdout/stderr tails.
