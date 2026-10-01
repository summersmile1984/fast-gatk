# GenotypeGVCFs rerun report — 2026-09-23

## Summary
- Total scripts: 25
- Passed: 25    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2138.171s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gatk_genotype_gvcf.py` | 0 | 29.241s | 62.34s | 12.1 MiB | 365.6 MiB | 28.54s | ✅ |
| `verify_gatk_genotype_gvcf_inbreeding.py` | 0 | 7.356s | 11.48s | 8.4 MiB | 268.4 MiB | 7.24s | ✅ |
| `verify_gatk_genotype_gvcf_legacy_qual.py` | 0 | 13.749s | 21.72s | 8.9 MiB | 305.9 MiB | 13.68s | ✅ |
| `verify_gatk_genotype_gvcf_multiallelic.py` | 0 | 6.401s | 5.30s | 8.4 MiB | 304.5 MiB | 6.34s | ✅ |
| `verify_gatk_genotype_gvcf_multisample.py` | 0 | 29.552s | 91.21s | 8.6 MiB | 319.3 MiB | 29.42s | ✅ |
| `verify_genotype_gvcf.py` | 0 | 8.705s | 41.70s | 12.4 MiB | 350.0 MiB | 7.06s | ✅ |
| `verify_genotype_gvcf_assignment_gatk_oracle.py` | 0 | 19.766s | 74.55s | 8.8 MiB | 397.1 MiB | 19.52s | ✅ |
| `verify_genotype_gvcf_confident_hom_alt_gatk_oracle.py` | 0 | 68.379s | 148.91s | 8.4 MiB | 316.9 MiB | 68.02s | ✅ |
| `verify_genotype_gvcf_dense_materialize_gatk_oracle.py` | 0 | 175.359s | 308.54s | 8.4 MiB | 320.9 MiB | 174.55s | ✅ |
| `verify_genotype_gvcf_depth_gate_gatk_oracle.py` | 0 | 66.306s | 186.90s | 8.6 MiB | 315.4 MiB | 66.11s | ✅ |
| `verify_genotype_gvcf_exclude_intervals_gatk_oracle.py` | 0 | 12.228s | 99.98s | 9.1 MiB | 360.1 MiB | 12.10s | ✅ |
| `verify_genotype_gvcf_gp_input_gatk_oracle.py` | 0 | 18.912s | 111.72s | 8.5 MiB | 344.5 MiB | 18.78s | ✅ |
| `verify_genotype_gvcf_header_order_gatk_oracle.py` | 0 | 52.049s | 214.74s | 8.2 MiB | 324.6 MiB | 44.78s | ✅ |
| `verify_genotype_gvcf_include_non_variant_gatk_oracle.py` | 0 | 7.297s | 153.77s | 12.5 MiB | 318.4 MiB | 6.99s | ✅ |
| `verify_genotype_gvcf_independent_corpus_gatk_oracle.py` | 0 | 255.567s | 854.08s | 271.6 MiB | 882.5 MiB | 111.40s | ✅ |
| `verify_genotype_gvcf_lowqual_gatk_oracle.py` | 0 | 206.790s | 433.37s | 8.3 MiB | 323.5 MiB | 206.30s | ✅ |
| `verify_genotype_gvcf_malformed_gatk_oracle.py` | 0 | 195.619s | 529.05s | 8.4 MiB | 319.1 MiB | 195.29s | ✅ |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` | 0 | 11.704s | 316.47s | 8.5 MiB | 392.4 MiB | 11.61s | ✅ |
| `verify_genotype_gvcf_multisample_reference_confidence_oracle.py` | 0 | 51.502s | 335.41s | 9.6 MiB | 319.6 MiB | 51.05s | ✅ |
| `verify_genotype_gvcf_reverse_trim_gatk_oracle.py` | 0 | 102.560s | 906.20s | 8.4 MiB | 323.5 MiB | 102.15s | ✅ |
| `verify_genotype_gvcf_spandel_gatk_oracle.py` | 0 | 570.321s | 1407.35s | 8.5 MiB | 327.2 MiB | 556.20s | ✅ |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py` | 0 | 29.514s | 543.51s | 8.4 MiB | 333.6 MiB | 29.39s | ✅ |
| `verify_genotype_gvcf_spanning_source_gq_gatk_oracle.py` | 0 | 104.053s | 1107.34s | 266.8 MiB | 860.7 MiB | 8.12s | ✅ |
| `verify_genotype_gvcf_star_only_locus_gatk_oracle.py` | 0 | 71.731s | 950.79s | 8.4 MiB | 327.0 MiB | 71.59s | ✅ |
| `verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py` | 0 | 23.510s | 917.16s | 8.5 MiB | 316.5 MiB | 23.40s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool genotype-gvcf
```

## JSON sidecar

See `genotype-gvcf-rerun-20260923.json` for full stdout/stderr tails.
