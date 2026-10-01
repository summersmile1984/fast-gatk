# HaplotypeCaller rerun report — 2026-09-27

## Summary
- Total scripts: 74
- Passed: 74    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2319.454s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_activity_region_controls.py` | 0 | 0.123s | 0.22s | 16.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_bam_intervals.py` | 0 | 0.475s | 1.70s | 13.9 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | 0 | 12.536s | 31.92s | 157.9 MiB | 332.3 MiB | 5.05s | ✅ |
| `verify_hc_af_zero_format_gatk_oracle.py` | 0 | 156.099s | 707.37s | 15.5 MiB | 377.7 MiB | 143.42s | ✅ |
| `verify_hc_alleles_deep_boundary.py` | 0 | 92.888s | 400.78s | 16.2 MiB | 362.8 MiB | 90.77s | ✅ |
| `verify_hc_alleles_deep_limits.py` | 0 | 71.408s | 198.12s | 17.6 MiB | 399.8 MiB | 69.20s | ✅ |
| `verify_hc_alleles_gatk_oracle.py` | 0 | 81.172s | 280.15s | 13.2 MiB | 385.8 MiB | 80.60s | ✅ |
| `verify_hc_alleles_overlap_gate_oracle.py` | 0 | 32.183s | 70.77s | 15.0 MiB | 445.8 MiB | 31.72s | ✅ |
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | 0 | 88.984s | 341.96s | 15.6 MiB | 323.6 MiB | 87.59s | ✅ |
| `verify_hc_assembly_region_boundary_gatk_oracle.py` | 0 | 5.777s | 10.83s | 18.4 MiB | 364.3 MiB | 5.42s | ✅ |
| `verify_hc_bp_resolution_gatk_oracle.py` | 0 | 3.865s | 14.65s | 6.4 MiB | 279.0 MiB | 3.82s | ✅ |
| `verify_hc_broad_gatk_oracle.py` | 0 | 45.065s | 119.36s | 19.0 MiB | 380.2 MiB | 42.80s | ✅ |
| `verify_hc_chr17_69k_70k_as_annotation_oracle.py` | 0 | 3.996s | 36.47s | 18.6 MiB | 350.0 MiB | 3.90s | ✅ |
| `verify_hc_chr17_69k_70k_gatk_oracle.py` | 0 | 39.140s | 148.42s | 22.2 MiB | 426.7 MiB | 22.00s | ✅ |
| `verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py` | 0 | 4.874s | 75.56s | 18.7 MiB | 364.1 MiB | 4.69s | ✅ |
| `verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py` | 0 | 4.674s | 81.80s | 18.6 MiB | 301.5 MiB | 4.51s | ✅ |
| `verify_hc_chr17_independent_regions.py` | 0 | 0.663s | 82.86s | 19.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_chr20_100k_200k_info_byte_equality_oracle.py` | 0 | 34.121s | 225.03s | 3254.2 MiB | 810.8 MiB | 12.19s | ✅ |
| `verify_hc_chr20_100k_nocall_gatk_oracle.py` | 0 | 160.719s | 1042.72s | 731.0 MiB | 865.8 MiB | 20.22s | ✅ |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 0 | 47.711s | 435.03s | 159.1 MiB | 421.1 MiB | 37.62s | ✅ |
| `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py` | 0 | 34.379s | 470.42s | 158.4 MiB | 333.0 MiB | 25.73s | ✅ |
| `verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py` | 0 | 11.545s | 287.93s | 255.9 MiB | 381.7 MiB | 5.31s | ✅ |
| `verify_hc_chr20_real_contract.py` | 0 | 33.384s | 512.52s | 159.2 MiB | 327.4 MiB | 9.65s | ✅ |
| `verify_hc_cigar_indel_activity_gatk_oracle.py` | 0 | 61.288s | 613.89s | 12.7 MiB | 354.3 MiB | 55.52s | ✅ |
| `verify_hc_complex_multiallelic_oracle.py` | 0 | 39.792s | 572.55s | 127.9 MiB | 375.9 MiB | 37.61s | ✅ |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py` | 0 | 30.274s | 540.83s | 154.6 MiB | 326.3 MiB | 27.58s | ✅ |
| `verify_hc_flow_filter_threshold_flag_oracle.py` | 0 | 2.109s | 436.81s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_filter_toggle_flag_oracle.py` | 0 | 1.538s | 471.79s | 18.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_final_5_flag_oracle.py` | 0 | 2.148s | 474.74s | 18.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_mode_flag_oracle.py` | 0 | 0.905s | 474.74s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_probability_threshold_flag_oracle.py` | 0 | 0.912s | 475.73s | 18.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_quality_disallow_fill_flag_oracle.py` | 0 | 2.019s | 477.64s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_scaling_quantization_mods_flag_oracle.py` | 0 | 2.443s | 479.70s | 18.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_t0_lump_symmetric_flag_oracle.py` | 0 | 1.518s | 480.86s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_forced_alleles_emission_gate_oracle.py` | 0 | 62.369s | 787.76s | 15.3 MiB | 426.6 MiB | 61.30s | ✅ |
| `verify_hc_genotype_priors.py` | 0 | 1.846s | 485.42s | 12.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_gq_bands_gatk_oracle.py` | 0 | 12.813s | 548.59s | 11.3 MiB | 279.7 MiB | 12.74s | ✅ |
| `verify_hc_graph_gatk_oracle.py` | 0 | 4.772s | 517.08s | 18.6 MiB | 273.8 MiB | 4.54s | ✅ |
| `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | 0 | 73.138s | 876.13s | 17.5 MiB | 322.5 MiB | 54.00s | ✅ |
| `verify_hc_gvcf_indel_end_gatk_oracle.py` | 0 | 52.981s | 815.64s | 16.6 MiB | 379.9 MiB | 51.75s | ✅ |
| `verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | 0 | 136.533s | 1239.46s | 15.4 MiB | 362.9 MiB | 120.66s | ✅ |
| `verify_hc_indel_zero_gatk_oracle.py` | 0 | 19.523s | 582.17s | 19.0 MiB | 366.8 MiB | 19.09s | ✅ |
| `verify_hc_informative_overlap_margin_gatk_oracle.py` | 0 | 19.260s | 733.61s | 18.5 MiB | 283.2 MiB | 18.47s | ✅ |
| `verify_hc_issue3845_gatk_oracle.py` | 0 | 19.221s | 756.36s | 20.1 MiB | 361.7 MiB | 18.73s | ✅ |
| `verify_hc_kmer_list_gatk_oracle.py` | 0 | 9.296s | 717.53s | 18.5 MiB | 298.4 MiB | 8.93s | ✅ |
| `verify_hc_likelihood_filter.py` | 0 | 4.658s | 760.84s | 13.6 MiB | 402.8 MiB | 4.50s | ✅ |
| `verify_hc_min_base_quality_boundary_gatk_oracle.py` | 0 | 11.154s | 838.91s | 33.9 MiB | 407.9 MiB | 9.76s | ✅ |
| `verify_hc_min_base_quality_gatk_oracle.py` | 0 | 9.888s | 826.14s | 18.4 MiB | 361.3 MiB | 9.25s | ✅ |
| `verify_hc_min_pruning_gatk_oracle.py` | 0 | 23.918s | 893.39s | 83.4 MiB | 361.6 MiB | 21.86s | ✅ |
| `verify_hc_multi_alt_as_family_pin_oracle.py` | 0 | 0.038s | 787.80s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_multi_input.py` | 0 | 20.448s | 905.31s | 102.9 MiB | 412.7 MiB | 16.75s | ✅ |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | 0 | 92.032s | 1291.96s | 14.6 MiB | 316.6 MiB | 90.73s | ✅ |
| `verify_hc_ploidy.py` | 0 | 0.331s | 826.63s | 13.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_ploidy_window_invariance_gatk_oracle.py` | 0 | 57.225s | 1128.02s | 157.4 MiB | 388.1 MiB | 51.22s | ✅ |
| `verify_hc_polyploid_gatk_oracle.py` | 0 | 5.361s | 843.26s | 13.2 MiB | 287.0 MiB | 5.28s | ✅ |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | 0 | 114.481s | 1346.72s | 15.5 MiB | 401.1 MiB | 100.46s | ✅ |
| `verify_hc_rcm_pl_range_gatk_oracle.py` | 0 | 9.094s | 910.01s | 7.5 MiB | 314.2 MiB | 9.04s | ✅ |
| `verify_hc_rcm_realignment.py` | 0 | 0.746s | 906.15s | 18.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_recover_all_gatk_oracle.py` | 0 | 28.492s | 1064.07s | 12.8 MiB | 384.9 MiB | 28.40s | ✅ |
| `verify_hc_region_streaming.py` | 0 | 11.414s | 933.60s | 19.0 MiB | 289.9 MiB | 9.32s | ✅ |
| `verify_hc_reverse_strand_indel.py` | 0 | 0.636s | 925.74s | 10.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_segmented_real_gatk_oracle.py` | 0 | 23.225s | 1089.15s | 15.6 MiB | 416.6 MiB | 22.94s | ✅ |
| `verify_hc_sites_only_gatk_oracle.py` | 0 | 46.273s | 1168.43s | 18.8 MiB | 339.0 MiB | 45.47s | ✅ |
| `verify_hc_softclip_contig_start_gatk_oracle.py` | 0 | 4.083s | 1047.42s | 13.6 MiB | 286.8 MiB | 3.96s | ✅ |
| `verify_hc_softclip_gatk_oracle.py` | 0 | 4.216s | 1069.48s | 19.4 MiB | 321.1 MiB | 4.10s | ✅ |
| `verify_hc_softclip_low_quality_gatk_oracle.py` | 0 | 17.818s | 1139.12s | 18.8 MiB | 390.4 MiB | 17.43s | ✅ |
| `verify_hc_span_del_qual_oracle.py` | 0 | 139.290s | 1531.82s | 16.6 MiB | 358.9 MiB | 138.37s | ✅ |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | 0 | 99.596s | 1455.24s | 15.4 MiB | 369.4 MiB | 99.01s | ✅ |
| `verify_hc_streaming.py` | 0 | 0.102s | 1128.50s | 12.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_variant_annotations.py` | 0 | 0.101s | 1128.74s | 18.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_window_invariance_gatk_oracle.py` | 0 | 71.760s | 1406.01s | 158.5 MiB | 359.2 MiB | 63.33s | ✅ |
| `verify_overlapping_quality_correction.py` | 0 | 0.045s | 1139.39s | 6.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_real_assembly_graph.py` | 0 | 0.505s | 1139.70s | 27.9 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_resource_limits.py` | 0 | 0.045s | 1139.81s | 7.1 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool hc-call
```

## JSON sidecar

See `hc-call-rerun-20260927.json` for full stdout/stderr tails.
