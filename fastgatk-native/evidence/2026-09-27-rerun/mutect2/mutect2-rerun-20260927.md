# Mutect2 rerun report — 2026-09-27

## Summary
- Total scripts: 40
- Passed: 39    Failed: 1    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 4084.136s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_cnv_somatic_e2e_hcc1143_oracle.py` | 0 | 0.998s | 2.96s | 15.6 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_fragment_aggregation_gatk_oracle.py` | 0 | 1.488s | 7.05s | 5.3 MiB | 107.6 MiB | 0.48s | ✅ |
| `verify_mutect2.py` | 0 | 131.694s | 272.39s | 35.6 MiB | 395.5 MiB | 10.62s | ✅ |
| `verify_mutect2_as_annotation_oracle.py` | 0 | 398.604s | 992.80s | 697.8 MiB | 832.5 MiB | 87.60s | ✅ |
| `verify_mutect2_as_value_byte_equality_oracle.py` | 0 | 392.949s | 839.31s | 696.6 MiB | 846.1 MiB | 84.41s | ✅ |
| `verify_mutect2_bp_resolution_gatk_contract.py` | 0 | 5.034s | 12.74s | 11.4 MiB | 313.7 MiB | 4.94s | ✅ |
| `verify_mutect2_dream_low_bq_gatk_oracle.py` | 0 | 154.692s | 485.92s | 298.5 MiB | 441.7 MiB | 24.83s | ✅ |
| `verify_mutect2_dream_synthetic_oracle.py` | 0 | 349.600s | 687.78s | 697.0 MiB | 1001.4 MiB | 83.70s | ✅ |
| `verify_mutect2_feature_resource_gatk_oracle.py` | 0 | 55.838s | 50.54s | 16.0 MiB | 416.9 MiB | 38.77s | ✅ |
| `verify_mutect2_force_active_gatk_oracle.py` | 0 | 148.723s | 397.96s | 28.5 MiB | 433.8 MiB | 26.39s | ✅ |
| `verify_mutect2_gatk_oracle.py` | 0 | 56.915s | 92.17s | 28.5 MiB | 427.0 MiB | 19.89s | ✅ |
| `verify_mutect2_gvcf_eventmap_gatk_contract.py` | 0 | 19.445s | 113.88s | 80.4 MiB | 430.8 MiB | 6.86s | ✅ |
| `verify_mutect2_gvcf_matched_normal_gatk_contract.py` | 0 | 10.780s | 100.02s | 4.3 MiB | 344.8 MiB | 10.74s | ✅ |
| `verify_mutect2_gvcf_multitumor_gatk_contract.py` | 0 | 11.483s | 123.20s | 7.4 MiB | 391.3 MiB | 11.29s | ✅ |
| `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | 0 | 10.061s | 131.54s | 12.6 MiB | 403.3 MiB | 9.83s | ✅ |
| `verify_mutect2_hcc1143_chr20_oracle.py` | 0 | 521.579s | 1709.38s | 3018.1 MiB | 460.3 MiB | 52.74s | ✅ |
| `verify_mutect2_hcc1143_reference_boundary_oracle.py` | 0 | 5.184s | 135.35s | 12.0 MiB | 401.3 MiB | 4.99s | ✅ |
| `verify_mutect2_hcc1143_spot_gatk_oracle.py` | 0 | 19.816s | 147.11s | 34.0 MiB | 423.8 MiB | 13.01s | ✅ |
| `verify_mutect2_independent_mates_gatk_contract.py` | 0 | 2.135s | 149.07s | 24.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_issue3845_gatk_oracle.py` | 0 | 16.712s | 163.86s | 18.5 MiB | 440.1 MiB | 15.71s | ✅ |
| `verify_mutect2_itr_artifact_gatk_contract.py` | 0 | 1.915s | 274.73s | 12.0 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_kmer_list_gatk_oracle.py` | 0 | 26.885s | 506.49s | 159.9 MiB | 432.9 MiB | 17.50s | ✅ |
| `verify_mutect2_mismapping_rate_boundary.py` | 0 | 14.943s | 291.48s | 28.3 MiB | 422.4 MiB | 11.08s | ✅ |
| `verify_mutect2_mito_interval_halo_gatk_oracle.py` | 0 | 464.718s | 1907.16s | 3165.6 MiB | 831.5 MiB | 47.48s | ✅ |
| `verify_mutect2_mito_realign_gatk_oracle.py` | 0 | 461.512s | 1807.82s | 3165.6 MiB | 577.0 MiB | 46.71s | ✅ |
| `verify_mutect2_mitochondria_gatk_oracle.py` | 1 | 328.909s | 1123.52s | 909.6 MiB | 624.7 MiB | 28.35s | ❌ |
| `verify_mutect2_multinormal_rg_namespace_gatk_oracle.py` | 0 | 14.400s | 512.04s | 159.0 MiB | 393.4 MiB | 5.83s | ✅ |
| `verify_mutect2_multitumor_format_gatk_oracle.py` | 0 | 64.739s | 536.01s | 160.9 MiB | 448.2 MiB | 20.18s | ✅ |
| `verify_mutect2_normal_lod_gatk_oracle.py` | 0 | 18.341s | 549.35s | 13.4 MiB | 307.0 MiB | 17.80s | ✅ |
| `verify_mutect2_pcr_overlap_gatk_contract.py` | 0 | 0.093s | 549.48s | 7.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_recheck_assembly_resultset_joint.py` | 0 | 23.881s | 567.46s | 28.4 MiB | 343.2 MiB | 17.85s | ✅ |
| `verify_mutect2_recheck_normal_replay.py` | 0 | 245.550s | 1224.40s | 0.0 MiB | 808.6 MiB | 54.79s | ✅ |
| `verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py` | 0 | 75.932s | 1058.77s | 338.6 MiB | 428.9 MiB | 37.99s | ✅ |
| `verify_mutect2_tlod_formula.py` | 0 | 0.895s | 840.01s | 28.4 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_zero_lod_indel_activity_gatk_oracle.py` | 0 | 13.825s | 1008.49s | 16.5 MiB | 409.6 MiB | 12.92s | ✅ |
| `verify_pairhmm_default_indel_quality_gatk_oracle.py` | 0 | 7.615s | 998.90s | 28.4 MiB | 412.0 MiB | 5.83s | ✅ |
| `verify_pairhmm_results_oracle.py` | 0 | 6.096s | 1015.59s | 8.2 MiB | 424.4 MiB | 5.83s | ✅ |
| `verify_somatic_likelihood_oracle.py` | 0 | 0.085s | 1009.00s | 5.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_normal_lod_oracle.py` | 0 | 0.038s | 1009.05s | 4.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_posterior_normal_oracle.py` | 0 | 0.034s | 1009.10s | 5.0 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool mutect2
```

## JSON sidecar

See `mutect2-rerun-20260927.json` for full stdout/stderr tails.
