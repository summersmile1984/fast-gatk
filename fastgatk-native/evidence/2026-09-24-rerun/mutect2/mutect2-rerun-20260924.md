# Mutect2 rerun report — 2026-09-24

## Summary
- Total scripts: 38
- Passed: 37    Failed: 1    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 3122.887s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_cnv_somatic_e2e_hcc1143_oracle.py` | 0 | 0.674s | 2.21s | 17.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_fragment_aggregation_gatk_oracle.py` | 0 | 2.168s | 8.58s | 0.0 MiB | 109.5 MiB | 0.77s | ✅ |
| `verify_mutect2.py` | 0 | 39.008s | 148.11s | 30.3 MiB | 410.1 MiB | 12.18s | ✅ |
| `verify_mutect2_bp_resolution_gatk_contract.py` | 0 | 6.178s | 15.10s | 11.5 MiB | 298.1 MiB | 6.07s | ✅ |
| `verify_mutect2_dream_low_bq_gatk_oracle.py` | 0 | 75.783s | 316.14s | 297.8 MiB | 442.7 MiB | 28.05s | ✅ |
| `verify_mutect2_dream_synthetic_oracle.py` | 0 | 367.899s | 664.76s | 690.2 MiB | 1050.7 MiB | 51.15s | ✅ |
| `verify_mutect2_feature_resource_gatk_oracle.py` | 0 | 52.268s | 200.78s | 15.6 MiB | 436.7 MiB | 49.14s | ✅ |
| `verify_mutect2_force_active_gatk_oracle.py` | 0 | 93.788s | 404.83s | 22.5 MiB | 450.2 MiB | 32.97s | ✅ |
| `verify_mutect2_gatk_oracle.py` | 0 | 28.626s | 68.18s | 18.4 MiB | 428.2 MiB | 19.36s | ✅ |
| `verify_mutect2_gvcf_eventmap_gatk_contract.py` | 0 | 8.700s | 22.87s | 20.0 MiB | 420.6 MiB | 7.00s | ✅ |
| `verify_mutect2_gvcf_matched_normal_gatk_contract.py` | 0 | 10.193s | 32.65s | 4.3 MiB | 351.5 MiB | 10.15s | ✅ |
| `verify_mutect2_gvcf_multitumor_gatk_contract.py` | 0 | 16.633s | 44.25s | 7.8 MiB | 389.2 MiB | 16.48s | ✅ |
| `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | 0 | 14.364s | 78.31s | 12.0 MiB | 418.1 MiB | 14.26s | ✅ |
| `verify_mutect2_hcc1143_chr20_oracle.py` | 0 | 488.346s | 1541.43s | 2931.7 MiB | 497.1 MiB | 35.15s | ✅ |
| `verify_mutect2_hcc1143_reference_boundary_oracle.py` | 0 | 5.014s | 84.48s | 12.1 MiB | 424.4 MiB | 4.87s | ✅ |
| `verify_mutect2_hcc1143_spot_gatk_oracle.py` | 0 | 17.311s | 163.04s | 31.1 MiB | 422.9 MiB | 12.62s | ✅ |
| `verify_mutect2_independent_mates_gatk_contract.py` | 0 | 1.478s | 85.93s | 18.9 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_issue3845_gatk_oracle.py` | 0 | 17.523s | 218.89s | 18.8 MiB | 426.9 MiB | 16.73s | ✅ |
| `verify_mutect2_itr_artifact_gatk_contract.py` | 0 | 0.343s | 149.58s | 11.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_kmer_list_gatk_oracle.py` | 0 | 33.821s | 260.81s | 160.2 MiB | 452.7 MiB | 23.10s | ✅ |
| `verify_mutect2_mismapping_rate_boundary.py` | 0 | 15.491s | 238.69s | 18.4 MiB | 419.2 MiB | 11.59s | ✅ |
| `verify_mutect2_mito_interval_halo_gatk_oracle.py` | 0 | 461.576s | 1055.40s | 3159.4 MiB | 856.1 MiB | 50.80s | ✅ |
| `verify_mutect2_mito_realign_gatk_oracle.py` | 0 | 460.154s | 954.87s | 3159.4 MiB | 565.5 MiB | 49.71s | ✅ |
| `verify_mutect2_mitochondria_gatk_oracle.py` | 1 | 329.468s | 735.64s | 884.9 MiB | 641.3 MiB | 29.23s | ❌ |
| `verify_mutect2_multinormal_rg_namespace_gatk_oracle.py` | 0 | 15.010s | 322.58s | 155.2 MiB | 424.2 MiB | 6.01s | ✅ |
| `verify_mutect2_multitumor_format_gatk_oracle.py` | 0 | 60.747s | 462.84s | 159.2 MiB | 431.8 MiB | 20.26s | ✅ |
| `verify_mutect2_normal_lod_gatk_oracle.py` | 0 | 24.212s | 419.73s | 12.9 MiB | 262.6 MiB | 23.74s | ✅ |
| `verify_mutect2_pcr_overlap_gatk_contract.py` | 0 | 0.081s | 405.21s | 7.8 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_recheck_assembly_resultset_joint.py` | 0 | 27.784s | 437.43s | 22.1 MiB | 375.6 MiB | 23.32s | ✅ |
| `verify_mutect2_recheck_normal_replay.py` | 0 | 328.435s | 853.14s | 668.9 MiB | 800.9 MiB | 52.12s | ✅ |
| `verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py` | 0 | 84.041s | 534.95s | 334.2 MiB | 430.1 MiB | 41.27s | ✅ |
| `verify_mutect2_tlod_formula.py` | 0 | 0.898s | 463.52s | 18.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_mutect2_zero_lod_indel_activity_gatk_oracle.py` | 0 | 19.431s | 474.16s | 15.8 MiB | 420.8 MiB | 18.31s | ✅ |
| `verify_pairhmm_default_indel_quality_gatk_oracle.py` | 0 | 8.139s | 480.60s | 18.4 MiB | 421.0 MiB | 6.36s | ✅ |
| `verify_pairhmm_results_oracle.py` | 0 | 7.158s | 487.67s | 8.1 MiB | 414.0 MiB | 6.81s | ✅ |
| `verify_somatic_likelihood_oracle.py` | 0 | 0.058s | 487.86s | 5.2 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_normal_lod_oracle.py` | 0 | 0.047s | 487.95s | 4.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_somatic_posterior_normal_oracle.py` | 0 | 0.039s | 488.00s | 4.3 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool mutect2
```

## JSON sidecar

See `mutect2-rerun-20260924.json` for full stdout/stderr tails.
