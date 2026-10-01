# HaplotypeCaller rerun report — 2026-09-23

## Summary
- Total scripts: 61
- Passed: 61    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2623.018s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_activity_region_controls.py` | 0 | 0.148s | 2846.19s | 15.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_bam_intervals.py` | 0 | 0.465s | 2847.44s | 14.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | 0 | 15.792s | 2874.27s | 158.4 MiB | 315.8 MiB | 12.53s | ✅ |
| `verify_hc_af_zero_format_gatk_oracle.py` | 0 | 199.358s | 3595.37s | 15.5 MiB | 369.4 MiB | 198.02s | ✅ |
| `verify_hc_alleles_deep_boundary.py` | 0 | 133.073s | 3228.59s | 15.5 MiB | 367.3 MiB | 131.18s | ✅ |
| `verify_hc_alleles_deep_limits.py` | 0 | 80.766s | 2996.80s | 16.0 MiB | 365.7 MiB | 79.44s | ✅ |
| `verify_hc_alleles_gatk_oracle.py` | 0 | 108.597s | 3147.66s | 12.4 MiB | 393.6 MiB | 95.30s | ✅ |
| `verify_hc_alleles_overlap_gate_oracle.py` | 0 | 56.554s | 2927.97s | 14.8 MiB | 314.0 MiB | 55.97s | ✅ |
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | 0 | 106.332s | 3102.98s | 14.6 MiB | 305.0 MiB | 105.57s | ✅ |
| `verify_hc_assembly_region_boundary_gatk_oracle.py` | 0 | 8.620s | 2854.59s | 18.0 MiB | 419.1 MiB | 8.46s | ✅ |
| `verify_hc_bp_resolution_gatk_oracle.py` | 0 | 5.948s | 2857.76s | 8.1 MiB | 248.1 MiB | 5.91s | ✅ |
| `verify_hc_broad_gatk_oracle.py` | 0 | 46.059s | 2957.43s | 24.0 MiB | 334.4 MiB | 37.70s | ✅ |
| `verify_hc_chr17_69k_70k_gatk_oracle.py` | 0 | 37.642s | 2899.05s | 22.0 MiB | 422.0 MiB | 23.42s | ✅ |
| `verify_hc_chr17_independent_regions.py` | 0 | 0.321s | 2899.74s | 18.8 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_chr20_100k_nocall_gatk_oracle.py` | 0 | 130.138s | 3425.37s | 758.6 MiB | 915.8 MiB | 8.49s | ✅ |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 0 | 46.988s | 3059.12s | 158.2 MiB | 381.5 MiB | 32.52s | ✅ |
| `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py` | 0 | 41.648s | 3028.95s | 158.2 MiB | 378.7 MiB | 37.03s | ✅ |
| `verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py` | 0 | 10.792s | 3003.86s | 158.1 MiB | 340.5 MiB | 6.58s | ✅ |
| `verify_hc_chr20_real_contract.py` | 0 | 31.307s | 3170.74s | 158.6 MiB | 374.5 MiB | 17.29s | ✅ |
| `verify_hc_cigar_indel_activity_gatk_oracle.py` | 0 | 83.864s | 3473.61s | 11.9 MiB | 366.5 MiB | 83.61s | ✅ |
| `verify_hc_complex_multiallelic_oracle.py` | 0 | 35.665s | 3261.97s | 127.3 MiB | 392.5 MiB | 27.92s | ✅ |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py` | 0 | 36.782s | 3283.53s | 157.2 MiB | 340.3 MiB | 34.01s | ✅ |
| `verify_hc_forced_alleles_emission_gate_oracle.py` | 0 | 57.264s | 3319.50s | 15.1 MiB | 404.4 MiB | 56.48s | ✅ |
| `verify_hc_genotype_priors.py` | 0 | 1.455s | 3174.52s | 11.4 MiB | 109.2 MiB | 0.32s | ✅ |
| `verify_hc_gq_bands_gatk_oracle.py` | 0 | 12.035s | 3236.11s | 11.0 MiB | 266.7 MiB | 11.96s | ✅ |
| `verify_hc_graph_gatk_oracle.py` | 0 | 5.681s | 3240.84s | 18.2 MiB | 274.8 MiB | 5.45s | ✅ |
| `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | 0 | 65.051s | 3628.28s | 16.0 MiB | 388.1 MiB | 42.05s | ✅ |
| `verify_hc_gvcf_indel_end_gatk_oracle.py` | 0 | 64.269s | 3656.49s | 15.5 MiB | 412.9 MiB | 63.37s | ✅ |
| `verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | 0 | 159.005s | 3924.29s | 15.4 MiB | 394.5 MiB | 158.19s | ✅ |
| `verify_hc_indel_zero_gatk_oracle.py` | 0 | 22.333s | 3293.78s | 19.3 MiB | 349.3 MiB | 21.82s | ✅ |
| `verify_hc_informative_overlap_margin_gatk_oracle.py` | 0 | 20.265s | 3441.77s | 18.0 MiB | 310.1 MiB | 19.80s | ✅ |
| `verify_hc_issue3845_gatk_oracle.py` | 0 | 24.993s | 3497.13s | 18.2 MiB | 404.4 MiB | 24.67s | ✅ |
| `verify_hc_kmer_list_gatk_oracle.py` | 0 | 12.246s | 3513.07s | 18.0 MiB | 302.4 MiB | 12.04s | ✅ |
| `verify_hc_likelihood_filter.py` | 0 | 7.003s | 3502.65s | 12.2 MiB | 402.3 MiB | 6.89s | ✅ |
| `verify_hc_min_base_quality_boundary_gatk_oracle.py` | 0 | 18.404s | 3669.57s | 18.8 MiB | 450.5 MiB | 18.04s | ✅ |
| `verify_hc_min_base_quality_gatk_oracle.py` | 0 | 21.294s | 3679.99s | 18.0 MiB | 286.1 MiB | 20.94s | ✅ |
| `verify_hc_min_pruning_gatk_oracle.py` | 0 | 33.810s | 3729.64s | 19.9 MiB | 371.9 MiB | 33.13s | ✅ |
| `verify_hc_multi_input.py` | 0 | 18.410s | 3695.54s | 32.9 MiB | 437.5 MiB | 6.57s | ✅ |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | 0 | 118.958s | 4032.90s | 13.4 MiB | 471.8 MiB | 105.37s | ✅ |
| `verify_hc_ploidy.py` | 0 | 0.237s | 3628.94s | 12.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_ploidy_window_invariance_gatk_oracle.py` | 0 | 64.845s | 3851.80s | 157.8 MiB | 362.3 MiB | 61.56s | ✅ |
| `verify_hc_polyploid_gatk_oracle.py` | 0 | 10.522s | 3684.51s | 12.1 MiB | 329.2 MiB | 10.46s | ✅ |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | 0 | 112.981s | 3979.78s | 14.9 MiB | 365.0 MiB | 101.64s | ✅ |
| `verify_hc_rcm_pl_range_gatk_oracle.py` | 0 | 10.169s | 3711.72s | 6.8 MiB | 254.5 MiB | 10.13s | ✅ |
| `verify_hc_rcm_realignment.py` | 0 | 0.263s | 3685.00s | 18.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_recover_all_gatk_oracle.py` | 0 | 38.090s | 3774.09s | 11.1 MiB | 315.7 MiB | 38.03s | ✅ |
| `verify_hc_region_streaming.py` | 0 | 6.350s | 3703.45s | 18.9 MiB | 296.3 MiB | 5.55s | ✅ |
| `verify_hc_reverse_strand_indel.py` | 0 | 0.189s | 3707.92s | 10.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_segmented_real_gatk_oracle.py` | 0 | 27.830s | 3757.33s | 14.5 MiB | 411.3 MiB | 27.68s | ✅ |
| `verify_hc_sites_only_gatk_oracle.py` | 0 | 36.301s | 3800.47s | 18.1 MiB | 376.1 MiB | 35.78s | ✅ |
| `verify_hc_softclip_contig_start_gatk_oracle.py` | 0 | 5.707s | 3734.20s | 12.5 MiB | 366.3 MiB | 5.63s | ✅ |
| `verify_hc_softclip_gatk_oracle.py` | 0 | 6.486s | 3739.66s | 18.6 MiB | 375.3 MiB | 6.34s | ✅ |
| `verify_hc_softclip_low_quality_gatk_oracle.py` | 0 | 22.421s | 3812.08s | 18.1 MiB | 322.4 MiB | 22.10s | ✅ |
| `verify_hc_span_del_qual_oracle.py` | 0 | 186.830s | 4230.46s | 16.1 MiB | 392.2 MiB | 185.65s | ✅ |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | 0 | 114.090s | 4148.08s | 14.9 MiB | 423.5 MiB | 113.58s | ✅ |
| `verify_hc_streaming.py` | 0 | 0.095s | 3800.80s | 11.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_variant_annotations.py` | 0 | 0.096s | 3801.04s | 18.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_window_invariance_gatk_oracle.py` | 0 | 99.818s | 4095.65s | 159.0 MiB | 396.2 MiB | 93.77s | ✅ |
| `verify_overlapping_quality_correction.py` | 0 | 0.042s | 3812.37s | 6.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_real_assembly_graph.py` | 0 | 0.279s | 3812.65s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_resource_limits.py` | 0 | 0.042s | 3812.78s | 7.2 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool hc-call
```

## JSON sidecar

See `hc-call-rerun-20260923.json` for full stdout/stderr tails.
