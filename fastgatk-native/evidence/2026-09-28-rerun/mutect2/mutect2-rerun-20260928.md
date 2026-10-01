# Mutect2 rerun report — 2026-09-28

## Summary
- Total scripts: 40
- Passed: 40    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 1512.051s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_cnv_somatic_e2e_hcc1143_oracle.py` | 0 | 0.334s | 3.76s | 16.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_fragment_aggregation_gatk_oracle.py` | 0 | 1.106s | 6.96s | 6.0 MiB | 109.6 MiB | 0.28s | ✅ |
| `verify_mutect2.py` | 0 | 29.562s | 64.04s | 35.2 MiB | 406.9 MiB | 10.37s | ✅ |
| `verify_mutect2_as_annotation_oracle.py` | 0 | 97.942s | 182.89s | 714.6 MiB | 920.4 MiB | 25.17s | ✅ |
| `verify_mutect2_as_value_byte_equality_oracle.py` | 0 | 93.851s | 300.62s | 716.0 MiB | 941.0 MiB | 24.03s | ✅ |
| `verify_mutect2_bp_resolution_gatk_contract.py` | 0 | 3.866s | 186.48s | 7.6 MiB | 271.6 MiB | 3.83s | ✅ |
| `verify_mutect2_dream_low_bq_gatk_oracle.py` | 0 | 37.471s | 346.05s | 298.9 MiB | 453.4 MiB | 17.98s | ✅ |
| `verify_mutect2_dream_synthetic_oracle.py` | 0 | 52.414s | 436.58s | 716.0 MiB | 1819.8 MiB | 16.96s | ✅ |
| `verify_mutect2_feature_resource_gatk_oracle.py` | 0 | 43.497s | 466.35s | 15.9 MiB | 524.1 MiB | 42.74s | ✅ |
| `verify_mutect2_force_active_gatk_oracle.py` | 0 | 30.744s | 551.54s | 28.6 MiB | 441.9 MiB | 19.72s | ✅ |
| `verify_mutect2_gatk_oracle.py` | 0 | 13.750s | 483.91s | 28.5 MiB | 440.2 MiB | 12.68s | ✅ |
| `verify_mutect2_gvcf_eventmap_gatk_contract.py` | 0 | 5.652s | 491.15s | 81.7 MiB | 437.6 MiB | 4.83s | ✅ |
| `verify_mutect2_gvcf_matched_normal_gatk_contract.py` | 0 | 12.561s | 567.43s | 2.0 MiB | 341.1 MiB | 12.53s | ✅ |
| `verify_mutect2_gvcf_multitumor_gatk_contract.py` | 0 | 7.877s | 559.94s | 7.7 MiB | 400.4 MiB | 7.83s | ✅ |
| `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | 0 | 24.004s | 575.60s | 11.9 MiB | 399.5 MiB | 23.95s | ✅ |
| `verify_mutect2_hcc1143_chr20_oracle.py` | 0 | 130.227s | 983.87s | 3055.5 MiB | 1085.5 MiB | 19.59s | ✅ |
| `verify_mutect2_hcc1143_reference_boundary_oracle.py` | 0 | 4.223s | 579.33s | 12.2 MiB | 396.7 MiB | 4.16s | ✅ |
| `verify_mutect2_hcc1143_spot_gatk_oracle.py` | 0 | 10.680s | 590.49s | 34.4 MiB | 417.6 MiB | 8.94s | ✅ |
| `verify_mutect2_independent_mates_gatk_contract.py` | 0 | 1.159s | 593.00s | 24.7 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_issue3845_gatk_oracle.py` | 0 | 12.777s | 608.12s | 18.6 MiB | 423.6 MiB | 12.48s | ✅ |
| `verify_mutect2_itr_artifact_gatk_contract.py` | 0 | 0.222s | 608.67s | 12.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_kmer_list_gatk_oracle.py` | 0 | 15.794s | 628.07s | 158.6 MiB | 446.2 MiB | 12.89s | ✅ |
| `verify_mutect2_mismapping_rate_boundary.py` | 0 | 16.578s | 644.72s | 28.6 MiB | 411.7 MiB | 13.24s | ✅ |
| `verify_mutect2_mito_interval_halo_gatk_oracle.py` | 0 | 173.623s | 1082.22s | 3165.7 MiB | 918.7 MiB | 16.53s | ✅ |
| `verify_mutect2_mito_realign_gatk_oracle.py` | 0 | 170.649s | 1178.24s | 3165.7 MiB | 850.6 MiB | 15.88s | ✅ |
| `verify_mutect2_mitochondria_gatk_oracle.py` | 0 | 319.974s | 1656.60s | 3045.6 MiB | 686.6 MiB | 10.58s | ✅ |
| `verify_mutect2_multinormal_rg_namespace_gatk_oracle.py` | 0 | 6.237s | 1183.89s | 159.8 MiB | 418.1 MiB | 4.14s | ✅ |
| `verify_mutect2_multitumor_format_gatk_oracle.py` | 0 | 24.168s | 1207.50s | 159.4 MiB | 455.5 MiB | 13.47s | ✅ |
| `verify_mutect2_normal_lod_gatk_oracle.py` | 0 | 16.494s | 1221.53s | 13.6 MiB | 352.2 MiB | 16.33s | ✅ |
| `verify_mutect2_pcr_overlap_gatk_contract.py` | 0 | 0.047s | 1221.62s | 7.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_recheck_assembly_resultset_joint.py` | 0 | 14.521s | 1237.87s | 31.1 MiB | 378.7 MiB | 13.47s | ✅ |
| `verify_mutect2_recheck_normal_replay.py` | 0 | 82.740s | 1333.37s | 692.5 MiB | 850.3 MiB | 16.33s | ✅ |
| `verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py` | 0 | 36.771s | 1374.60s | 336.2 MiB | 423.5 MiB | 27.02s | ✅ |
| `verify_mutect2_tlod_formula.py` | 0 | 0.236s | 1375.22s | 27.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_zero_lod_indel_activity_gatk_oracle.py` | 0 | 11.204s | 1384.90s | 16.5 MiB | 443.4 MiB | 5.61s | ✅ |
| `verify_pairhmm_default_indel_quality_gatk_oracle.py` | 0 | 4.661s | 1390.76s | 28.6 MiB | 424.7 MiB | 4.21s | ✅ |
| `verify_pairhmm_results_oracle.py` | 0 | 4.348s | 1397.59s | 8.3 MiB | 424.5 MiB | 4.22s | ✅ |
| `verify_somatic_likelihood_oracle.py` | 0 | 0.037s | 1397.64s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_normal_lod_oracle.py` | 0 | 0.026s | 1397.69s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_posterior_normal_oracle.py` | 0 | 0.024s | 1397.73s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool mutect2
```

## JSON sidecar

See `mutect2-rerun-20260928.json` for full stdout/stderr tails.
