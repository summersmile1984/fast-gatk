#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
build_script="$root/fastgatk-native/scripts/build_native.sh"
build_dir="$root/fastgatk-native/build"
cmake_bin="$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake"
python_bin="${PYTHON:-python3}"

if [[ "${FASTGATK_SKIP_BUILD:-0}" != "1" ]]; then
    bash "$build_script"
fi

[[ -x "$build_dir/fastgatk-hc-call" ]] || {
    echo "missing native build: $build_dir/fastgatk-hc-call" >&2
    exit 2
}

"$cmake_bin" --build "$build_dir" --target fastgatk-runtime-smoke \
    fastgatk-runtime-pipeline-smoke --parallel >/dev/null
"$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
    --test-dir "$build_dir" --output-on-failure

# Keep the default sweep fast enough for developer/CI feedback while still
# exercising the second Kokkos execution space and the lightweight score,
# fixture-digest, and scheduler resource-limit gates on every run. Set
# FASTGATK_VERIFY_SERIAL_FULL=1 for the complete current Serial CTest matrix
# (the count is reported by ctest -N).
serial_build="${FASTGATK_SERIAL_BUILD:-$root/fastgatk-native/build-serial}"
if [[ -x "$serial_build/fastgatk-hc-call" && -x "$serial_build/fastgatk-genotype-gvcf" ]]; then
    serial_ctest_args=(--test-dir "$serial_build" --output-on-failure)
    if [[ "${FASTGATK_VERIFY_SERIAL_FULL:-0}" != "1" ]]; then
        serial_ctest_args+=( -R 'fastgatk-progress-score-contract|fastgatk-fixture-digest-contract|fastgatk-kokkos-api-boundary|fastgatk-bqsr-contract|fastgatk-bqsr-gatk-oracle|fastgatk-bqsr-context-size-gatk-oracle|fastgatk-apply-bqsr-alias-gatk-oracle|fastgatk-analyze-covariates-bqsr-alias-gatk-oracle|fastgatk-bqsr-indel-gatk-oracle|fastgatk-bqsr-read-filter-gatk-oracle|fastgatk-bqsr-preserve-gatk-oracle|fastgatk-bqsr-report-roundtrip-gatk-oracle|fastgatk-bqsr-cram-gatk-oracle|fastgatk-gather-bqsr-gatk-oracle|fastgatk-hc-streaming-contract|fastgatk-hc-region-streaming-contract|fastgatk-hc-reverse-strand-indel-oracle|fastgatk-hc-min-pruning-gatk-oracle|fastgatk-hc-kmer-list-gatk-oracle|fastgatk-hc-softclip-gatk-oracle|fastgatk-launcher-contract|fastgatk-genotype-gvcf-contract|fastgatk-genotype-gvcf-gatk-oracle|fastgatk-genotype-gvcf-include-non-variant-gatk-oracle|fastgatk-genotype-gvcf-assignment-gatk-oracle|fastgatk-genotype-gvcf-legacy-qual-gatk-oracle|fastgatk-genotype-gvcf-exclude-intervals-gatk-oracle|fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle|fastgatk-mutect2-contract|fastgatk-mutect2-gatk-oracle|fastgatk-mutect2-tlod-formula|fastgatk-mutect2-pairhmm-default-indel-quality-gatk-oracle|fastgatk-filter-mutect-calls-contract|fastgatk-filter-mutect-calls-contamination-oracle|fastgatk-filter-mutect-calls-germline-oracle|fastgatk-filter-mutect-calls-orientation-gatk-oracle|fastgatk-filter-mutect-calls-variant-index-alias-gatk-oracle|fastgatk-select-variants-contract|fastgatk-select-variants-gatk-oracle|fastgatk-select-variants-sites-only-gatk-oracle|fastgatk-select-variants-filtered-nocall-gatk-oracle|fastgatk-select-variants-filtered-oracle|fastgatk-variant-filtration-contract|fastgatk-variant-filtration-gatk-oracle|fastgatk-variant-filtration-set-nocall-gatk-oracle|fastgatk-left-align-trim-contract|fastgatk-collect-read-counts-contract|fastgatk-collect-read-counts-gatk-oracle|fastgatk-collect-allelic-counts-contract|fastgatk-collect-allelic-counts-gatk-oracle|fastgatk-call-copy-ratio-segments-contract|fastgatk-call-copy-ratio-segments-gatk-oracle|fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle|fastgatk-call-copy-ratio-segments-compensated-sum-gatk-oracle|fastgatk-runtime-output-smoke|fastgatk-depth-of-coverage-contract|fastgatk-depth-of-coverage-read-filter-gatk-oracle|fastgatk-resource-limits-contract|fastgatk-slurm-resource-wrapper-contract|fastgatk-scheduler-retry-contract|fastgatk-mark-duplicates-contract|fastgatk-mark-duplicates-optical-only-gatk-oracle|fastgatk-sort-sam-contract|fastgatk-sort-sam-gatk-oracle|fastgatk-get-pileup-summaries-contract|fastgatk-get-pileup-summaries-gatk-oracle|fastgatk-denoise-read-counts-interval-identity-gatk-oracle|fastgatk-denoise-read-counts-integer-input-gatk-oracle|fastgatk-denoise-read-counts-hdf5-metadata-gatk-oracle|fastgatk-variant-recalibrator-zero-variance-gatk-oracle|fastgatk-variant-recalibrator-attempt-iteration-gatk-oracle|fastgatk-genomicsdb-import-update-workspace-gatk-oracle|fastgatk-genomicsdb-import-native-interval-gatk-oracle|fastgatk-model-segments-multisample-gatk-oracle|fastgatk-hc-chr20-real-contract' )
    fi
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" "${serial_ctest_args[@]}"
    # These real chr20 gVCF byte oracles are intentionally always replayed on
    # Serial too: their AF, physical-phasing, and spanning-deletion decisions
    # are shared Kokkos/Host behavior and must not become OpenMP-only guarantees.
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-hc-chr20-(min-pruning|max-mnp(-(haploid|polyploid(4|8)?))?)-gvcf-gatk-oracle$'
    # The bounded slices added after the last continuous 140-test baseline
    # are always replayed on Serial as well; keeping them separate preserves
    # the fast default sweep while preventing new compatibility gates from
    # being OpenMP-only in the default verification path.
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R 'fastgatk-combine-gvcfs-plless-gatk-oracle|fastgatk-filter-mutect-calls-orientation-joint-oracle|fastgatk-filter-mutect-calls-normal-artifact-oracle|fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle|fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle|fastgatk-model-segments-smoothing-gatk-oracle|fastgatk-somatic-posterior-normal-count-gatk-oracle|fastgatk-somatic-normal-lod-gatk-oracle|fastgatk-hc-assembly-region-boundary-gatk-oracle|fastgatk-hc-cigar-indel-activity-gatk-oracle|fastgatk-hc-sites-only-gatk-oracle|fastgatk-mutect2-mismapping-rate-boundary|fastgatk-mutect2-pairhmm-default-indel-quality-gatk-oracle|fastgatk-apply-vqsr-default-cutoff-gatk-oracle|fastgatk-apply-vqsr-exclude-intervals-gatk-oracle|fastgatk-left-align-trim-sites-only-gatk-oracle|fastgatk-left-align-cli-boundary-gatk-oracle|fastgatk-gather-vcfs-cli-boundary-gatk-oracle|fastgatk-mutect2-force-active-gatk-oracle|fastgatk-mutect2-normal-lod-gatk-oracle|fastgatk-sort-sam-cli-boundary-gatk-oracle|fastgatk-validate-variants-symbolic-gatk-oracle|fastgatk-denoise-read-counts-hdf5-metadata-gatk-oracle|fastgatk-apply-bqsr-alias-gatk-oracle|fastgatk-analyze-covariates-bqsr-alias-gatk-oracle|fastgatk-fasta-alternate-iupac-hom-gatk-oracle|fastgatk-select-variants-sites-only-gatk-oracle|fastgatk-mark-duplicates-optical-only-gatk-oracle|fastgatk-variant-recalibrator-zero-variance-gatk-oracle|fastgatk-variant-recalibrator-attempt-iteration-gatk-oracle|fastgatk-depth-of-coverage-read-filter-gatk-oracle|fastgatk-variant-eval-keep-ac0-gatk-oracle|fastgatk-model-segments-multisample-gatk-oracle|fastgatk-genomicsdb-import-native-interval-gatk-oracle'
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-bqsr-long-read-gatk-oracle$'
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-mark-duplicates-pair-key-gatk-oracle$'
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-depth-of-coverage-ignore-deletion-sites-gatk-oracle$'
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-model-segments-allele-fraction-likelihood-gatk-oracle$'
    "$root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest" \
        --test-dir "$serial_build" --output-on-failure \
        -R '^fastgatk-genotype-gvcf-starts-in-intervals-gatk-oracle$'
else
    echo '{"status":"skip","suite":"fastgatk-native-serial-matrix","reason":"build-serial not present"}'
fi

"$python_bin" "$root/fastgatk-native/dispatcher/verify_dispatcher.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_native.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_multi_input.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bam_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_ploidy.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_real_assembly_graph.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_hc_ploidy.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_hc_likelihood_filter.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_activity_region_controls.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_count_bases_in_reference.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_count_bases_in_reference.py" \
    --bases 2000000 --threads 2
"$python_bin" "$root/fastgatk-native/scripts/verify_compare_references.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_compare_references.py" \
    --bases 2000000 --threads 2
"$python_bin" "$root/fastgatk-native/scripts/verify_check_reference_compatibility.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_check_reference_compatibility.py" \
    --contigs 128 --length 1000
"$python_bin" "$root/fastgatk-native/scripts/verify_fasta_reference_tools.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_fasta_alternate_iupac_hom_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_fasta_reference_tools.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_shift_fasta.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_shift_fasta.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_index_feature_file.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_index_feature_file.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_count_reads.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_flag_stat.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_read_metrics.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_annotate_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_annotate_intervals.py" \
    --intervals 10000 --interval-size 1000 --threads 2 --with-tracks
"$python_bin" "$root/fastgatk-native/scripts/verify_split_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_split_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_filter_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_preprocess_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_preprocess_intervals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_gather_bqsr_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_gather_bqsr_reports.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_indel_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_preserve_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_read_filter_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_context_size_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_long_read_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_bqsr_alias_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_report_roundtrip_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_bqsr_cram_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_bqsr.py" \
    --iterations 2 --threads 2 --batch-records 64
"$python_bin" "$root/fastgatk-native/scripts/verify_analyze_covariates.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_analyze_covariates_bqsr_alias_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_analyze_covariates.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genotype_gvcf.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genotype_gvcf_assignment_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genotype_gvcf_exclude_intervals_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genotype_gvcf_spanning_deletion_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_genotype_gvcf.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_genotype_gvcf.py" --stream-by-locus
"$python_bin" "$root/fastgatk-native/scripts/verify_combine_gvcfs.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_combine_gvcfs_interval_refblock_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_combine_gvcfs.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_reblock_gvcf.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_reblock_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_reblock_overlap_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_reblock_gatk_shards.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_reblock_gvcf.py"
# Always check the Java/native command-scope guardrail.  The measured mode is
# opt-in because it launches a JVM for every repetition; when enabled it uses
# the same pinned fixture, interval, output-index policy, GNU-time wall/RSS/
# filesystem-I/O schema, and refuses claims below the configured workload size.
"$python_bin" "$root/fastgatk-native/scripts/benchmark_reblock_gvcf_java.py" --dry-run
if [[ "${FASTGATK_RUN_JAVA_BENCHMARK:-0}" == "1" ]]; then
    "$python_bin" "$root/fastgatk-native/scripts/benchmark_reblock_gvcf_java.py" \
        --repetitions "${FASTGATK_JAVA_BENCHMARK_REPETITIONS:-3}" \
        --minimum-records "${FASTGATK_JAVA_BENCHMARK_MIN_RECORDS:-1024}"
fi
"$python_bin" "$root/fastgatk-native/scripts/verify_select_variants.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_select_variants_sites_only_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_select_variants_filtered_nocall_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_select_variants_filtered_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_select_variants.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_eval.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_eval_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_variant_eval_keep_ac0_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_variant_eval.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_validate_variants.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_validate_variants_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_validate_variants_symbolic_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_validate_variants.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variants_to_table.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variants_to_table_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_variants_to_table.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr_default_cutoff_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr_sites_only_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr_exclude_intervals_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_apply_vqsr_filter_booleans_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_apply_vqsr.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_gatk_model_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_vbem_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_annotation_order_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_sample_every_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_zero_variance_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_variant_recalibrator_attempt_iteration_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_variant_recalibrator.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gather_tranches.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_gather_tranches.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gather_vcfs.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_gather_vcfs_cli_boundary_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_gather_vcfs.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_left_align.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_left_align_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_left_align_sites_only_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_left_align_cli_boundary_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_left_align.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_filtration.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_filtration_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_variant_filtration_missing_boolean_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_variant_filtration_set_nocall_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_variant_filtration.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_sort_sam.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_sort_sam_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_sort_sam_cli_boundary_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_sort_sam.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mark_duplicates.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mark_duplicates_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_mark_duplicates_tagging_policy_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_mark_duplicates_pair_key_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_mark_duplicates.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_import.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_import_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_import_sample_map_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_import_update_workspace_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_import_native_interval_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_genomicsdb_import.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_genomicsdb_bridge.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mutect2.py"
FASTGATK_MUTECT2_BINARY="$build_dir/fastgatk-mutect2" "$python_bin" "$root/fastgatk-native/scripts/verify_mutect2_tlod_formula.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mutect2_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mutect2_mismapping_rate_boundary.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mutect2_force_active_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_mutect2_normal_lod_gatk_oracle.py"
FASTGATK_SOMATIC_ORACLE_BINARY="$build_dir/fastgatk-kernels/fastgatk-somatic-likelihood-oracle" "$python_bin" "$root/fastgatk-native/scripts/verify_somatic_likelihood_oracle.py"
FASTGATK_SOMATIC_ORACLE_BINARY="$build_dir/fastgatk-kernels/fastgatk-somatic-likelihood-oracle" "$python_bin" "$root/fastgatk-native/scripts/verify_somatic_posterior_normal_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_pairhmm_default_indel_quality_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_mutect2.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_collect_f1r2_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_collect_f1r2_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_learn_read_orientation_model.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_learn_read_orientation_model_multisample_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_learn_read_orientation_model.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_combine_gvcfs.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_combine_gvcfs_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_combine_gvcfs_interval_refblock_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_calls.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_calls_contamination_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_calls_germline_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_orientation_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_variant_index_alias_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_normal_artifact_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_contamination_joint_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_filter_mutect_orientation_joint_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_optional_boolean_contract.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_filter_mutect_calls.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_depth_of_coverage.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_depth_of_coverage_multisample.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_depth_of_coverage_read_filter_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_depth_of_coverage.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_collect_allelic_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_collect_allelic_counts_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_collect_allelic_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_collect_read_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_collect_read_counts_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_collect_read_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hdf5_simple_count_collection.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_denoise_read_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_denoise_read_counts_interval_identity_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_denoise_read_counts_integer_input_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_denoise_read_counts.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_create_read_count_panel_of_normals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_create_read_count_panel_of_normals.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_model_segments.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_input_segments_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_copy_ratio_conditionals_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_allele_fraction_initialization_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_allele_fraction_likelihood_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_smoothing_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_default_kernel_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_first_alt_fraction_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_model_segments_multisample_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_model_segments.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_model_segments.py" \
    --points 20000 --threads 2 --kernel --window-size 8 --window-size 16
"$python_bin" "$root/fastgatk-native/scripts/benchmark_model_segments.py" \
    --points 20000 --threads 2 --matched-normal
"$python_bin" "$root/fastgatk-native/scripts/verify_call_copy_ratio_segments.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_call_copy_ratio_segments_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_call_copy_ratio_segments.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_get_pileup_summaries.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_get_pileup_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_get_pileup_summaries.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gather_pileup_summaries.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gather_pileup_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_pileup_validation.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_calculate_contamination.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_calculate_contamination_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_calculate_contamination.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_resource_limits.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_scheduler_retry_contract.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_kernel_benchmark.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_kokkos_backend_matrix.py"

# The Kokkos PairHMM architecture builds are optional in minimal images, but
# when the pinned GATK/JDK oracle and all three builds are present this is a
# hard raw-bit gate for the 4096-value matrix. Set
# FASTGATK_REQUIRE_PAIRHMM_ORACLE=1 to make missing assets a CI failure.
pairhmm_java="$root/third_party/jdk17/bin/java"
pairhmm_jar="$root/third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
pairhmm_tables="$root/pairhmm-demo/results/gatk-tables.hex"
pairhmm_ready=1
for pairhmm_variant in scalar avx2 avx512; do
    [[ -x "$root/pairhmm-demo/build-kokkos-$pairhmm_variant/pairhmm-kokkos" ]] || pairhmm_ready=0
done
if [[ -x "$pairhmm_java" && -f "$pairhmm_jar" && -f "$pairhmm_tables" && "$pairhmm_ready" == "1" ]]; then
    "$python_bin" "$root/pairhmm-demo/scripts/verify_kokkos.py" \
        --records=512 --workload=matrix8 --variants=scalar,avx2,avx512
    if [[ "${FASTGATK_SKIP_GATK_BENCHMARK:-0}" != "1" ]]; then
        "$python_bin" "$root/pairhmm-demo/scripts/benchmark_same_cores.py" \
            --native-backend=auto --cores=1,2,4 --pairs=1024 --iterations=2 --repeats=1 \
            --output="$root/pairhmm-demo/results/benchmark-kokkos-same-core-20260825.json"
    else
        echo '{"status":"skip","suite":"pairhmm-kokkos-same-core-benchmark","reason":"FASTGATK_SKIP_GATK_BENCHMARK=1"}'
    fi
elif [[ "${FASTGATK_REQUIRE_PAIRHMM_ORACLE:-0}" == "1" ]]; then
    echo "missing pinned PairHMM oracle/build assets" >&2
    exit 2
else
    echo '{"status":"skip","suite":"pairhmm-kokkos-strict-oracle","reason":"oracle/build assets unavailable"}'
fi

# This script is safe in images without a local GATK jar: the oracle verifier
# reports a structured skip. Set FASTGATK_REQUIRE_GATK_ORACLE=1 to make that
# absence a hard CI failure.
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_broad_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_min_pruning_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_kmer_list_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_softclip_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_gq_bands_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_polyploid_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_complex_multiallelic_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/benchmark_hc_broad.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_bp_resolution_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_rcm_pl_range_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_rcm_realignment.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_hc_sites_only_gatk_oracle.py"
FASTGATK_REQUIRE_GATK_ORACLE=1 "$python_bin" "$root/fastgatk-native/scripts/verify_hc_cigar_indel_activity_gatk_oracle.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_genotype_gvcf.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_genotype_gvcf_legacy_qual.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_genotype_gvcf_multisample.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_genotype_gvcf_multiallelic.py"
"$python_bin" "$root/fastgatk-native/scripts/verify_gatk_genotype_gvcf_inbreeding.py"

# Nextflow/SLURM checks are kept separate because a real scheduler allocation
# is environment-dependent. The script always runs local equivalents and
# executes real Nextflow when FASTGATK_NEXTFLOW_BIN (or the pinned local jar)
# is available.
bash "$root/fastgatk-native/scripts/verify_pipeline_local.sh"

echo '{"status":"pass","suite":"fastgatk-native-all"}'
