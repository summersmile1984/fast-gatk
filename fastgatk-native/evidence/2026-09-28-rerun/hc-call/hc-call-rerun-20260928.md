# HaplotypeCaller rerun report — 2026-09-28

## Summary
- Total scripts: 74
- Passed: 74    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2369.056s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_activity_region_controls.py` | 0 | 0.127s | 3974.70s | 16.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_bam_intervals.py` | 0 | 0.467s | 3975.94s | 14.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | 0 | 13.041s | 4005.53s | 157.3 MiB | 310.0 MiB | 9.57s | ✅ |
| `verify_hc_af_zero_format_gatk_oracle.py` | 0 | 165.204s | 4612.41s | 15.8 MiB | 374.2 MiB | 163.61s | ✅ |
| `verify_hc_alleles_deep_boundary.py` | 0 | 112.207s | 4350.75s | 17.4 MiB | 410.6 MiB | 109.61s | ✅ |
| `verify_hc_alleles_deep_limits.py` | 0 | 66.606s | 4132.67s | 17.4 MiB | 387.1 MiB | 64.58s | ✅ |
| `verify_hc_alleles_gatk_oracle.py` | 0 | 86.140s | 4258.95s | 13.0 MiB | 378.3 MiB | 85.58s | ✅ |
| `verify_hc_alleles_overlap_gate_oracle.py` | 0 | 36.887s | 4042.75s | 14.8 MiB | 406.2 MiB | 36.54s | ✅ |
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | 0 | 80.742s | 4210.03s | 15.6 MiB | 399.3 MiB | 79.73s | ✅ |
| `verify_hc_assembly_region_boundary_gatk_oracle.py` | 0 | 5.962s | 3984.79s | 18.6 MiB | 370.7 MiB | 5.64s | ✅ |
| `verify_hc_bp_resolution_gatk_oracle.py` | 0 | 4.124s | 3988.25s | 6.3 MiB | 273.8 MiB | 4.09s | ✅ |
| `verify_hc_broad_gatk_oracle.py` | 0 | 29.546s | 4076.35s | 18.9 MiB | 390.9 MiB | 24.05s | ✅ |
| `verify_hc_chr17_69k_70k_as_annotation_oracle.py` | 0 | 4.819s | 4010.39s | 18.3 MiB | 320.9 MiB | 4.72s | ✅ |
| `verify_hc_chr17_69k_70k_gatk_oracle.py` | 0 | 59.053s | 4159.35s | 22.2 MiB | 417.3 MiB | 31.30s | ✅ |
| `verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py` | 0 | 9.266s | 4087.62s | 18.6 MiB | 306.0 MiB | 9.09s | ✅ |
| `verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py` | 0 | 4.263s | 4081.56s | 18.5 MiB | 361.2 MiB | 4.15s | ✅ |
| `verify_hc_chr17_independent_regions.py` | 0 | 0.490s | 4082.76s | 19.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_chr20_100k_200k_info_byte_equality_oracle.py` | 0 | 51.739s | 4292.66s | 3255.4 MiB | 821.2 MiB | 12.18s | ✅ |
| `verify_hc_chr20_100k_nocall_gatk_oracle.py` | 0 | 171.563s | 4987.77s | 731.0 MiB | 778.2 MiB | 15.34s | ✅ |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 0 | 60.240s | 4468.59s | 159.1 MiB | 367.8 MiB | 49.39s | ✅ |
| `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py` | 0 | 49.660s | 4436.83s | 158.4 MiB | 335.5 MiB | 41.27s | ✅ |
| `verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py` | 0 | 15.544s | 4266.04s | 256.1 MiB | 361.1 MiB | 10.72s | ✅ |
| `verify_hc_chr20_real_contract.py` | 0 | 31.947s | 4376.15s | 159.3 MiB | 314.2 MiB | 14.35s | ✅ |
| `verify_hc_cigar_indel_activity_gatk_oracle.py` | 0 | 73.984s | 4645.48s | 12.7 MiB | 372.2 MiB | 73.67s | ✅ |
| `verify_hc_complex_multiallelic_oracle.py` | 0 | 26.429s | 4403.72s | 127.7 MiB | 365.9 MiB | 24.24s | ✅ |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py` | 0 | 36.968s | 4510.50s | 157.3 MiB | 351.1 MiB | 32.90s | ✅ |
| `verify_hc_flow_filter_threshold_flag_oracle.py` | 0 | 0.870s | 4377.71s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_filter_toggle_flag_oracle.py` | 0 | 0.696s | 4378.79s | 18.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_final_5_flag_oracle.py` | 0 | 0.998s | 4380.34s | 18.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_mode_flag_oracle.py` | 0 | 0.699s | 4381.06s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_probability_threshold_flag_oracle.py` | 0 | 1.065s | 4381.73s | 18.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_quality_disallow_fill_flag_oracle.py` | 0 | 2.542s | 4405.35s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_scaling_quantization_mods_flag_oracle.py` | 0 | 3.001s | 4406.93s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_flow_t0_lump_symmetric_flag_oracle.py` | 0 | 1.847s | 4438.16s | 18.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_forced_alleles_emission_gate_oracle.py` | 0 | 45.925s | 4671.48s | 15.4 MiB | 351.0 MiB | 45.08s | ✅ |
| `verify_hc_genotype_priors.py` | 0 | 1.652s | 4473.70s | 12.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_gq_bands_gatk_oracle.py` | 0 | 9.056s | 4486.90s | 11.5 MiB | 260.7 MiB | 8.93s | ✅ |
| `verify_hc_graph_gatk_oracle.py` | 0 | 5.385s | 4478.89s | 18.5 MiB | 277.5 MiB | 5.20s | ✅ |
| `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | 0 | 60.725s | 4796.49s | 17.6 MiB | 332.3 MiB | 59.62s | ✅ |
| `verify_hc_gvcf_indel_end_gatk_oracle.py` | 0 | 46.865s | 4703.74s | 17.1 MiB | 411.4 MiB | 45.61s | ✅ |
| `verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | 0 | 118.221s | 5173.97s | 15.5 MiB | 376.6 MiB | 116.71s | ✅ |
| `verify_hc_indel_zero_gatk_oracle.py` | 0 | 15.310s | 4522.38s | 19.1 MiB | 305.9 MiB | 14.65s | ✅ |
| `verify_hc_informative_overlap_margin_gatk_oracle.py` | 0 | 19.673s | 4730.56s | 18.5 MiB | 397.3 MiB | 18.79s | ✅ |
| `verify_hc_issue3845_gatk_oracle.py` | 0 | 20.566s | 4752.64s | 20.1 MiB | 492.1 MiB | 19.95s | ✅ |
| `verify_hc_kmer_list_gatk_oracle.py` | 0 | 10.285s | 4714.43s | 0.0 MiB | 344.5 MiB | 9.74s | ✅ |
| `verify_hc_likelihood_filter.py` | 0 | 4.434s | 4675.85s | 13.6 MiB | 403.6 MiB | 4.26s | ✅ |
| `verify_hc_min_base_quality_boundary_gatk_oracle.py` | 0 | 11.001s | 4764.87s | 34.1 MiB | 395.2 MiB | 10.19s | ✅ |
| `verify_hc_min_base_quality_gatk_oracle.py` | 0 | 10.321s | 4806.45s | 18.3 MiB | 290.5 MiB | 9.82s | ✅ |
| `verify_hc_min_pruning_gatk_oracle.py` | 0 | 16.310s | 4823.76s | 83.3 MiB | 378.1 MiB | 14.88s | ✅ |
| `verify_hc_multi_alt_as_family_pin_oracle.py` | 0 | 0.032s | 4730.59s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_multi_input.py` | 0 | 13.864s | 4835.94s | 103.0 MiB | 426.0 MiB | 9.74s | ✅ |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | 0 | 107.167s | 5226.45s | 14.5 MiB | 369.8 MiB | 106.15s | ✅ |
| `verify_hc_ploidy.py` | 0 | 0.204s | 4765.24s | 13.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_ploidy_window_invariance_gatk_oracle.py` | 0 | 49.090s | 5062.76s | 159.4 MiB | 358.0 MiB | 43.47s | ✅ |
| `verify_hc_polyploid_gatk_oracle.py` | 0 | 9.624s | 4840.32s | 13.8 MiB | 349.8 MiB | 9.50s | ✅ |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | 0 | 119.502s | 5276.82s | 15.4 MiB | 423.7 MiB | 118.84s | ✅ |
| `verify_hc_rcm_pl_range_gatk_oracle.py` | 0 | 4.423s | 4845.27s | 11.1 MiB | 283.2 MiB | 4.38s | ✅ |
| `verify_hc_rcm_realignment.py` | 0 | 0.863s | 4841.04s | 19.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_recover_all_gatk_oracle.py` | 0 | 29.405s | 5007.77s | 12.8 MiB | 316.1 MiB | 29.31s | ✅ |
| `verify_hc_region_streaming.py` | 0 | 7.429s | 4880.52s | 19.1 MiB | 274.6 MiB | 5.53s | ✅ |
| `verify_hc_reverse_strand_indel.py` | 0 | 1.111s | 4873.86s | 12.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_segmented_real_gatk_oracle.py` | 0 | 26.757s | 5025.45s | 15.5 MiB | 424.5 MiB | 26.49s | ✅ |
| `verify_hc_sites_only_gatk_oracle.py` | 0 | 32.782s | 5093.40s | 19.1 MiB | 311.9 MiB | 32.08s | ✅ |
| `verify_hc_softclip_contig_start_gatk_oracle.py` | 0 | 9.543s | 4992.10s | 13.5 MiB | 279.2 MiB | 9.47s | ✅ |
| `verify_hc_softclip_gatk_oracle.py` | 0 | 9.732s | 5068.01s | 19.3 MiB | 321.2 MiB | 9.52s | ✅ |
| `verify_hc_softclip_low_quality_gatk_oracle.py` | 0 | 19.247s | 5105.26s | 18.9 MiB | 345.0 MiB | 18.87s | ✅ |
| `verify_hc_span_del_qual_oracle.py` | 0 | 137.772s | 5466.94s | 16.2 MiB | 404.1 MiB | 126.01s | ✅ |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | 0 | 98.981s | 5388.76s | 15.4 MiB | 414.2 MiB | 98.47s | ✅ |
| `verify_hc_streaming.py` | 0 | 0.140s | 5068.47s | 12.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_variant_annotations.py` | 0 | 0.180s | 5068.79s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_window_invariance_gatk_oracle.py` | 0 | 86.306s | 5337.73s | 159.3 MiB | 373.5 MiB | 81.07s | ✅ |
| `verify_overlapping_quality_correction.py` | 0 | 0.043s | 5093.58s | 6.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_real_assembly_graph.py` | 0 | 0.356s | 5093.93s | 27.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_resource_limits.py` | 0 | 0.038s | 5094.16s | 7.0 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool hc-call
```

## JSON sidecar

See `hc-call-rerun-20260928.json` for full stdout/stderr tails.
