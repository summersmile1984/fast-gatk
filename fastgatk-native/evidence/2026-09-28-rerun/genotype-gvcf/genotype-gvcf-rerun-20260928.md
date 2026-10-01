# GenotypeGVCFs rerun report — 2026-09-28

## Summary
- Total scripts: 25
- Passed: 25    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2038.662s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gatk_genotype_gvcf.py` | 0 | 40.737s | 2353.95s | 12.0 MiB | 324.7 MiB | 39.00s | ✅ |
| `verify_gatk_genotype_gvcf_inbreeding.py` | 0 | 11.232s | 2259.17s | 8.5 MiB | 289.1 MiB | 11.03s | ✅ |
| `verify_gatk_genotype_gvcf_legacy_qual.py` | 0 | 15.847s | 2304.50s | 8.7 MiB | 370.0 MiB | 15.72s | ✅ |
| `verify_gatk_genotype_gvcf_multiallelic.py` | 0 | 6.024s | 2250.17s | 8.3 MiB | 280.6 MiB | 5.93s | ✅ |
| `verify_gatk_genotype_gvcf_multisample.py` | 0 | 44.389s | 2373.48s | 8.5 MiB | 335.7 MiB | 44.28s | ✅ |
| `verify_genotype_gvcf.py` | 0 | 12.802s | 2291.03s | 13.3 MiB | 401.4 MiB | 10.62s | ✅ |
| `verify_genotype_gvcf_assignment_gatk_oracle.py` | 0 | 16.182s | 2322.96s | 8.7 MiB | 404.4 MiB | 11.08s | ✅ |
| `verify_genotype_gvcf_confident_hom_alt_gatk_oracle.py` | 0 | 66.334s | 2439.73s | 8.2 MiB | 330.3 MiB | 60.70s | ✅ |
| `verify_genotype_gvcf_dense_materialize_gatk_oracle.py` | 0 | 202.187s | 3014.68s | 8.4 MiB | 338.1 MiB | 185.78s | ✅ |
| `verify_genotype_gvcf_depth_gate_gatk_oracle.py` | 0 | 69.041s | 2504.52s | 8.1 MiB | 319.7 MiB | 68.89s | ✅ |
| `verify_genotype_gvcf_exclude_intervals_gatk_oracle.py` | 0 | 23.905s | 2331.47s | 9.1 MiB | 385.2 MiB | 23.82s | ✅ |
| `verify_genotype_gvcf_gp_input_gatk_oracle.py` | 0 | 33.789s | 2391.68s | 8.4 MiB | 355.1 MiB | 33.66s | ✅ |
| `verify_genotype_gvcf_header_order_gatk_oracle.py` | 0 | 58.842s | 2467.55s | 8.3 MiB | 302.3 MiB | 58.73s | ✅ |
| `verify_genotype_gvcf_include_non_variant_gatk_oracle.py` | 0 | 10.283s | 2378.55s | 13.4 MiB | 323.2 MiB | 9.75s | ✅ |
| `verify_genotype_gvcf_independent_corpus_gatk_oracle.py` | 0 | 202.571s | 3648.77s | 272.1 MiB | 891.2 MiB | 80.08s | ✅ |
| `verify_genotype_gvcf_lowqual_gatk_oracle.py` | 0 | 165.215s | 3111.59s | 8.3 MiB | 316.2 MiB | 149.90s | ✅ |
| `verify_genotype_gvcf_malformed_gatk_oracle.py` | 0 | 183.096s | 3210.63s | 8.2 MiB | 319.4 MiB | 178.20s | ✅ |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` | 0 | 9.375s | 2400.43s | 8.6 MiB | 404.8 MiB | 9.10s | ✅ |
| `verify_genotype_gvcf_multisample_reference_confidence_oracle.py` | 0 | 35.295s | 2525.91s | 9.3 MiB | 369.7 MiB | 35.15s | ✅ |
| `verify_genotype_gvcf_reverse_trim_gatk_oracle.py` | 0 | 108.793s | 2600.03s | 8.4 MiB | 310.5 MiB | 107.99s | ✅ |
| `verify_genotype_gvcf_spandel_gatk_oracle.py` | 0 | 473.132s | 3940.95s | 8.6 MiB | 333.4 MiB | 472.21s | ✅ |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py` | 0 | 28.339s | 2543.81s | 8.6 MiB | 336.2 MiB | 27.84s | ✅ |
| `verify_genotype_gvcf_spanning_source_gq_gatk_oracle.py` | 0 | 105.378s | 2900.83s | 267.6 MiB | 861.7 MiB | 10.79s | ✅ |
| `verify_genotype_gvcf_star_only_locus_gatk_oracle.py` | 0 | 85.273s | 2636.95s | 8.4 MiB | 330.9 MiB | 84.95s | ✅ |
| `verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py` | 0 | 30.601s | 2913.28s | 8.6 MiB | 273.2 MiB | 30.35s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool genotype-gvcf
```

## JSON sidecar

See `genotype-gvcf-rerun-20260928.json` for full stdout/stderr tails.
