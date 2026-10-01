# HaplotypeCaller rerun report — 2026-09-24

## Summary
- Total scripts: 61
- Passed: 60    Failed: 1    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 2430.612s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_activity_region_controls.py` | 0 | 0.128s | 0.23s | 15.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_bam_intervals.py` | 0 | 0.561s | 1.78s | 14.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | 0 | 15.121s | 38.78s | 157.5 MiB | 353.5 MiB | 11.08s | ✅ |
| `verify_hc_af_zero_format_gatk_oracle.py` | 0 | 177.413s | 874.95s | 15.5 MiB | 419.5 MiB | 175.12s | ✅ |
| `verify_hc_alleles_deep_boundary.py` | 0 | 92.815s | 427.88s | 15.6 MiB | 412.7 MiB | 90.90s | ✅ |
| `verify_hc_alleles_deep_limits.py` | 0 | 78.018s | 210.28s | 16.2 MiB | 343.5 MiB | 76.36s | ✅ |
| `verify_hc_alleles_gatk_oracle.py` | 0 | 84.315s | 276.28s | 12.3 MiB | 374.0 MiB | 83.73s | ✅ |
| `verify_hc_alleles_overlap_gate_oracle.py` | 0 | 50.989s | 81.22s | 14.7 MiB | 346.0 MiB | 50.49s | ✅ |
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | 0 | 92.611s | 351.97s | 14.9 MiB | 354.1 MiB | 91.36s | ✅ |
| `verify_hc_assembly_region_boundary_gatk_oracle.py` | 0 | 7.141s | 12.55s | 17.8 MiB | 373.1 MiB | 6.79s | ✅ |
| `verify_hc_bp_resolution_gatk_oracle.py` | 0 | 4.616s | 17.31s | 10.4 MiB | 278.4 MiB | 4.56s | ✅ |
| `verify_hc_broad_gatk_oracle.py` | 0 | 42.205s | 120.62s | 18.2 MiB | 370.3 MiB | 40.69s | ✅ |
| `verify_hc_chr17_69k_70k_gatk_oracle.py` | 0 | 51.555s | 152.08s | 22.6 MiB | 432.3 MiB | 32.99s | ✅ |
| `verify_hc_chr17_independent_regions.py` | 0 | 0.656s | 82.53s | 18.8 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_chr20_100k_nocall_gatk_oracle.py` | 0 | 166.925s | 1250.52s | 730.8 MiB | 781.6 MiB | 16.00s | ✅ |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 0 | 55.235s | 523.56s | 159.2 MiB | 435.1 MiB | 44.27s | ✅ |
| `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py` | 0 | 38.013s | 481.50s | 157.9 MiB | 397.3 MiB | 25.46s | ✅ |
| `verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py` | 0 | 11.547s | 285.80s | 157.9 MiB | 328.5 MiB | 6.56s | ✅ |
| `verify_hc_chr20_real_contract.py` | 0 | 33.821s | 560.44s | 159.3 MiB | 335.8 MiB | 16.22s | ✅ |
| `verify_hc_cigar_indel_activity_gatk_oracle.py` | 0 | 75.579s | 734.88s | 11.6 MiB | 375.2 MiB | 66.30s | ✅ |
| `verify_hc_complex_multiallelic_oracle.py` | 0 | 47.211s | 643.30s | 126.8 MiB | 364.4 MiB | 45.31s | ✅ |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py` | 0 | 34.674s | 608.86s | 157.2 MiB | 408.2 MiB | 30.85s | ✅ |
| `verify_hc_forced_alleles_emission_gate_oracle.py` | 0 | 58.725s | 692.35s | 15.2 MiB | 408.4 MiB | 41.96s | ✅ |
| `verify_hc_genotype_priors.py` | 0 | 1.836s | 528.70s | 11.7 MiB | 109.1 MiB | 0.37s | ✅ |
| `verify_hc_gq_bands_gatk_oracle.py` | 0 | 14.508s | 570.15s | 10.9 MiB | 289.5 MiB | 14.40s | ✅ |
| `verify_hc_graph_gatk_oracle.py` | 1 | 10.171s | 614.33s | 0.0 MiB | 0.0 MiB | 0.00s | ❌ |
| `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | 0 | 54.460s | 952.07s | 16.1 MiB | 375.5 MiB | 53.45s | ✅ |
| `verify_hc_gvcf_indel_end_gatk_oracle.py` | 0 | 50.899s | 911.47s | 15.4 MiB | 321.1 MiB | 43.06s | ✅ |
| `verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | 0 | 128.944s | 1421.16s | 15.5 MiB | 393.3 MiB | 127.30s | ✅ |
| `verify_hc_indel_zero_gatk_oracle.py` | 0 | 10.858s | 656.95s | 19.1 MiB | 387.8 MiB | 10.34s | ✅ |
| `verify_hc_informative_overlap_margin_gatk_oracle.py` | 0 | 25.573s | 752.43s | 18.0 MiB | 336.6 MiB | 24.85s | ✅ |
| `verify_hc_issue3845_gatk_oracle.py` | 0 | 37.290s | 1043.06s | 18.8 MiB | 403.3 MiB | 36.67s | ✅ |
| `verify_hc_kmer_list_gatk_oracle.py` | 0 | 16.088s | 965.33s | 17.9 MiB | 284.9 MiB | 15.55s | ✅ |
| `verify_hc_likelihood_filter.py` | 0 | 10.051s | 971.88s | 12.6 MiB | 392.6 MiB | 9.84s | ✅ |
| `verify_hc_min_base_quality_boundary_gatk_oracle.py` | 0 | 16.027s | 986.94s | 18.6 MiB | 318.9 MiB | 15.08s | ✅ |
| `verify_hc_min_base_quality_gatk_oracle.py` | 0 | 17.151s | 1001.15s | 17.8 MiB | 307.0 MiB | 16.60s | ✅ |
| `verify_hc_min_pruning_gatk_oracle.py` | 0 | 22.768s | 1064.76s | 19.8 MiB | 380.8 MiB | 21.60s | ✅ |
| `verify_hc_multi_input.py` | 0 | 14.919s | 1015.11s | 32.2 MiB | 449.8 MiB | 10.82s | ✅ |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | 0 | 121.743s | 1639.43s | 13.6 MiB | 378.9 MiB | 120.73s | ✅ |
| `verify_hc_ploidy.py` | 0 | 0.258s | 988.09s | 12.8 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_ploidy_window_invariance_gatk_oracle.py` | 0 | 68.895s | 1504.58s | 158.3 MiB | 368.1 MiB | 63.75s | ✅ |
| `verify_hc_polyploid_gatk_oracle.py` | 0 | 10.333s | 1075.52s | 12.3 MiB | 320.8 MiB | 10.22s | ✅ |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | 0 | 104.999s | 1576.36s | 15.1 MiB | 377.4 MiB | 97.77s | ✅ |
| `verify_hc_rcm_pl_range_gatk_oracle.py` | 0 | 4.912s | 1070.20s | 7.9 MiB | 295.2 MiB | 4.85s | ✅ |
| `verify_hc_rcm_realignment.py` | 0 | 0.520s | 1065.95s | 18.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_recover_all_gatk_oracle.py` | 0 | 37.721s | 1313.43s | 12.0 MiB | 349.7 MiB | 37.59s | ✅ |
| `verify_hc_region_streaming.py` | 0 | 8.462s | 1134.12s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_reverse_strand_indel.py` | 0 | 2.603s | 1122.91s | 11.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_segmented_real_gatk_oracle.py` | 0 | 30.078s | 1289.89s | 14.1 MiB | 426.1 MiB | 29.81s | ✅ |
| `verify_hc_sites_only_gatk_oracle.py` | 0 | 46.247s | 1457.04s | 18.0 MiB | 332.7 MiB | 45.51s | ✅ |
| `verify_hc_softclip_contig_start_gatk_oracle.py` | 0 | 5.305s | 1257.89s | 0.0 MiB | 305.0 MiB | 5.18s | ✅ |
| `verify_hc_softclip_gatk_oracle.py` | 0 | 10.025s | 1264.32s | 18.7 MiB | 339.4 MiB | 9.90s | ✅ |
| `verify_hc_softclip_low_quality_gatk_oracle.py` | 0 | 10.168s | 1328.77s | 18.1 MiB | 362.7 MiB | 9.84s | ✅ |
| `verify_hc_span_del_qual_oracle.py` | 0 | 150.134s | 1842.79s | 15.7 MiB | 378.1 MiB | 143.34s | ✅ |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | 0 | 88.887s | 1761.71s | 15.1 MiB | 373.7 MiB | 88.28s | ✅ |
| `verify_hc_streaming.py` | 0 | 0.106s | 1329.53s | 11.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_variant_annotations.py` | 0 | 0.122s | 1329.79s | 17.9 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_hc_window_invariance_gatk_oracle.py` | 0 | 77.087s | 1706.01s | 158.3 MiB | 354.1 MiB | 67.14s | ✅ |
| `verify_overlapping_quality_correction.py` | 0 | 0.048s | 1421.37s | 6.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_real_assembly_graph.py` | 0 | 0.495s | 1421.68s | 18.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_resource_limits.py` | 0 | 0.047s | 1421.82s | 7.5 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool hc-call
```

## JSON sidecar

See `hc-call-rerun-20260924.json` for full stdout/stderr tails.
