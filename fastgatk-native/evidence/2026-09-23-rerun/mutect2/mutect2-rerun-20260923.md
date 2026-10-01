# Mutect2 rerun report — historical status (2026-09-23)

## Evidence correction (2026-09-24)

The historical resource data are **unvalidated and invalid as performance
evidence**. Manually estimated values were mixed with sampled measurements;
CPU values used cumulative child-process accounting; sampler finalization
could omit tracked processes; and native/Java totals cover unlike workloads
and invocation counts. The former elapsed/CPU/RSS table and optimization
claims are withdrawn rather than reproduced as verified measurements.

The historical report recorded 38 scripts, 38 successful exits, no failures
and no skips. These are retained script-status records, not a fresh rerun
or proof of output parity. Exit 0 does not identify which outputs were
compared, establish full GATK parity, or validate resource measurements.
The JSON sidecar is retained unchanged for traceability.

## Historical script exit records

| Script | Recorded exit |
| --- | --- |
| `verify_cnv_somatic_e2e_hcc1143_oracle.py` | 0 |
| `verify_fragment_aggregation_gatk_oracle.py` | 0 |
| `verify_mutect2.py` | 0 |
| `verify_mutect2_bp_resolution_gatk_contract.py` | 0 |
| `verify_mutect2_dream_low_bq_gatk_oracle.py` | 0 |
| `verify_mutect2_dream_synthetic_oracle.py` | 0 |
| `verify_mutect2_feature_resource_gatk_oracle.py` | 0 |
| `verify_mutect2_force_active_gatk_oracle.py` | 0 |
| `verify_mutect2_gatk_oracle.py` | 0 |
| `verify_mutect2_gvcf_eventmap_gatk_contract.py` | 0 |
| `verify_mutect2_gvcf_matched_normal_gatk_contract.py` | 0 |
| `verify_mutect2_gvcf_multitumor_gatk_contract.py` | 0 |
| `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | 0 |
| `verify_mutect2_hcc1143_chr20_oracle.py` | 0 |
| `verify_mutect2_hcc1143_reference_boundary_oracle.py` | 0 |
| `verify_mutect2_hcc1143_spot_gatk_oracle.py` | 0 |
| `verify_mutect2_independent_mates_gatk_contract.py` | 0 |
| `verify_mutect2_issue3845_gatk_oracle.py` | 0 |
| `verify_mutect2_itr_artifact_gatk_contract.py` | 0 |
| `verify_mutect2_kmer_list_gatk_oracle.py` | 0 |
| `verify_mutect2_mismapping_rate_boundary.py` | 0 |
| `verify_mutect2_mito_interval_halo_gatk_oracle.py` | 0 |
| `verify_mutect2_mito_realign_gatk_oracle.py` | 0 |
| `verify_mutect2_mitochondria_gatk_oracle.py` | 0 |
| `verify_mutect2_multinormal_rg_namespace_gatk_oracle.py` | 0 |
| `verify_mutect2_multitumor_format_gatk_oracle.py` | 0 |
| `verify_mutect2_normal_lod_gatk_oracle.py` | 0 |
| `verify_mutect2_pcr_overlap_gatk_contract.py` | 0 |
| `verify_mutect2_recheck_assembly_resultset_joint.py` | 0 |
| `verify_mutect2_recheck_normal_replay.py` | 0 |
| `verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py` | 0 |
| `verify_mutect2_tlod_formula.py` | 0 |
| `verify_mutect2_zero_lod_indel_activity_gatk_oracle.py` | 0 |
| `verify_pairhmm_default_indel_quality_gatk_oracle.py` | 0 |
| `verify_pairhmm_results_oracle.py` | 0 |
| `verify_somatic_likelihood_oracle.py` | 0 |
| `verify_somatic_normal_lod_oracle.py` | 0 |
| `verify_somatic_posterior_normal_oracle.py` | 0 |

## Previous patch claims withdrawn

The old two-stage note described verifier harness changes to native thread
and assembly-region parameters and Java `-Xmx512m`. It did not establish a
native-source optimization. Its before/after timings, speedups, RSS
reductions and assertions that estimates were “real driver-sampled
measurements” are withdrawn. Neither unchanged nor reduced native memory
has been established by those records. The `read_all_many` root-cause
claim likewise requires profiling rather than inference from aggregate RSS.

## Fresh evidence pending

Fresh isolated measurements and output-parity checks are pending in
[work/mutect2-priority/](../../../../work/mutect2-priority/). No completed
optimization result is claimed. Comparisons must use identified binaries,
matching inputs and semantic parameters, raw process measurements and
independent output checks, not changed halo/downsampling settings.

See the [regression evidence workflow](../../../docs/regression-evidence.md)
and [resource report correction](../NATIVE_VS_JAVA.md).
