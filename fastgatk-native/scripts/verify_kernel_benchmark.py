#!/usr/bin/env python3
"""Validate the reproducible Kokkos kernel benchmark contract."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-kernels/fastgatk-kernels-benchmark"
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"missing benchmark executable: {binary}")
    repeats = 3
    environment = os.environ.copy()
    environment["FASTGATK_BENCH_WARMUP"] = "1"
    environment["FASTGATK_BENCH_REPEATS"] = str(repeats)
    result = subprocess.run(
        [str(binary)], check=True, text=True, capture_output=True, env=environment
    )
    artifact = json.loads(result.stdout.strip().splitlines()[-1])
    assert artifact["schema_version"] == 7
    assert artifact["status"] == "pass"
    assert artifact["warmup"] == 1
    assert artifact["repeats"] == repeats
    assert artifact["requests"] > 0
    assert artifact["simd_width"] >= 1
    assert artifact["compiled_simd_width"] == artifact["simd_width"]
    assert artifact["compiled_simd_backend"] in {"scalar", "avx2", "avx512"}
    if artifact["compiled_simd_backend"] == "avx512":
        assert artifact["cpu_supports_avx512"] is True
    if artifact["compiled_simd_backend"] == "avx2":
        assert artifact["cpu_supports_avx2"] is True
    assert artifact["sw_simd_width"] >= 1
    assert artifact["sw_simd_groups"] > 0
    assert artifact["concurrency"] >= 1
    assert artifact["sw_checksum"] == 335360
    assert artifact["sw_uniform_checksum"] == 368640
    assert artifact["pairhmm_checksum"] == artifact["persistent_pairhmm_checksum"]
    assert artifact["float_pairhmm_simd_width"] >= 1
    assert artifact["float_pairhmm_checksum"] != 0.0
    assert artifact["float_pairhmm_max_abs_error"] >= 0.0
    # The float path is an explicitly approximate mode.  This gate catches
    # overflow/underflow or an accidentally unscaled recurrence without
    # pretending that float arithmetic can satisfy the strict raw-bit oracle.
    assert artifact["float_pairhmm_max_abs_error"] < 1.0e-2
    assert artifact["pairhmm_likelihood_normalization_checksum"] > 0
    assert artifact["pairhmm_allele_marginalization_checksum"] > 0
    assert artifact["pairhmm_read_allele_uncertainty_checksum"] > 0
    assert artifact["pairhmm_read_allele_best_checksum"] > 0
    assert artifact["flow_pairhmm_error_model"] == "flow"
    assert artifact["flow_pairhmm_execution_policy"] == "TeamPolicy"
    assert artifact["flow_pairhmm_checksum"] != 0.0
    assert artifact["genotype_checksum"] > 0
    assert artifact["joint_genotype_checksum"] > 0
    assert artifact["genotype_prior_checksum"] > 0
    assert artifact["posterior_assignment_checksum"] > 0
    assert artifact["allele_count_checksum"] == 1024
    assert artifact["posterior_checksum"] > 0
    assert artifact["cross_sample_reference_checksum"] > 0
    assert artifact["cohort_checksum"] > 0
    assert artifact["reference_confidence_checksum"] > 0
    assert artifact["reference_confidence_polyploid_checksum"] > 0
    assert artifact["somatic_checksum"] > 0
    assert artifact["somatic_posterior_checksum"] > 0
    assert artifact["somatic_multiallelic_checksum"] > 0
    assert artifact["bqsr_count_checksum"] > 0
    assert artifact["bqsr_apply_checksum"] > 0
    assert artifact["timings"]["genotype_execution_space"]
    assert artifact["timings"]["joint_genotype_execution_space"]
    assert artifact["timings"]["genotype_prior_execution_space"]
    assert artifact["timings"]["posterior_assignment_execution_space"]
    assert artifact["timings"]["allele_count_execution_space"]
    assert artifact["timings"]["posterior_execution_space"]
    assert artifact["timings"]["cross_sample_reference_execution_space"]
    assert artifact["timings"]["cohort_execution_space"]
    assert artifact["timings"]["reference_confidence_execution_space"]
    assert artifact["timings"]["reference_confidence_polyploid_execution_space"]
    assert artifact["timings"]["somatic_execution_space"]
    assert artifact["timings"]["somatic_posterior_execution_space"]
    assert artifact["timings"]["somatic_multiallelic_execution_space"]
    assert artifact["timings"]["bqsr_count_execution_space"]
    assert artifact["timings"]["bqsr_count_execution_policy"] == "TeamPolicy"
    assert artifact["timings"]["bqsr_covariate_execution_space"]
    assert artifact["timings"]["bqsr_covariate_execution_policy"] == "TeamPolicy"
    assert artifact["timings"]["bqsr_covariate_workspace_bytes"] > 0
    assert artifact["timings"]["bqsr_covariate_team_local_histogram"] is True
    assert artifact["timings"]["bqsr_apply_execution_space"]
    assert artifact["timings"]["bqsr_apply_execution_policy"] == "RangePolicy"
    assert artifact["timings"]["pairhmm_read_allele_uncertainty_execution_space"]
    assert artifact["timings"]["pairhmm_read_allele_best_execution_space"]
    assert artifact["persistent_pairhmm_cached_shapes"] > 0
    assert artifact["persistent_pairhmm_cache_hits"] >= repeats - 1
    assert 1 <= artifact["graph_kmer_size_selected"] <= 31
    assert artifact["graph_kmer_iterations"] >= 1
    assert isinstance(artifact["graph_has_non_reference_cycles"], bool)
    assert artifact["graph_dangling_branch_paths"] >= 0
    assert artifact["graph_dangling_branch_bases"] >= 0
    assert artifact["graph_dangling_recovered_paths"] >= 0
    assert artifact["graph_dangling_recovered_bases"] >= 0
    assert artifact["graph_artificial_haplotype_recovery_paths"] >= 0
    assert artifact["graph_artificial_haplotype_recovery_bases"] >= 0
    assert artifact["graph_seqgraph_nodes"] >= 0
    assert artifact["graph_seqgraph_edges"] >= 0
    assert artifact["graph_non_unique_kmers"] >= 0
    assert artifact["error_correction_solid_kmers"] >= 0
    assert artifact["error_correction_corrected_kmers"] >= 0
    assert artifact["error_correction_uncorrectable_kmers"] >= 0
    timings = artifact["timings"]["prepare_and_execute"]
    expected = {
        "sw_prepare",
        "sw_kernel",
        "sw_uniform_prepare",
        "sw_uniform_kernel",
        "pairhmm_prepare",
        "pairhmm_kernel",
        "float_pairhmm_prepare",
        "float_pairhmm_kernel",
        "pairhmm_likelihood_normalization_prepare",
        "pairhmm_likelihood_normalization_kernel",
        "pairhmm_allele_marginalization_prepare",
        "pairhmm_allele_marginalization_kernel",
        "pairhmm_read_allele_uncertainty_prepare",
        "pairhmm_read_allele_uncertainty_kernel",
        "pairhmm_read_allele_best_prepare",
        "pairhmm_read_allele_best_kernel",
        "flow_pairhmm_prepare",
        "flow_pairhmm_kernel",
        "persistent_pairhmm_prepare",
        "persistent_pairhmm_kernel",
        "genotype_prepare",
        "genotype_kernel",
        "joint_genotype_prepare",
        "joint_genotype_kernel",
        "genotype_prior_prepare",
        "genotype_prior_kernel",
        "genotype_posterior_assignment_prepare",
        "genotype_posterior_assignment_kernel",
        "genotype_allele_count_prepare",
        "genotype_allele_count_kernel",
        "genotype_posterior_prepare",
        "genotype_posterior_kernel",
        "cross_sample_reference_prepare",
        "cross_sample_reference_kernel",
        "genotype_cohort_prepare",
        "genotype_cohort_kernel",
        "reference_confidence_prepare",
        "reference_confidence_kernel",
        "reference_confidence_polyploid_prepare",
        "reference_confidence_polyploid_kernel",
        "somatic_prepare",
        "somatic_kernel",
        "somatic_posterior_prepare",
        "somatic_posterior_kernel",
        "somatic_multiallelic_prepare",
        "somatic_multiallelic_kernel",
        "bqsr_count_prepare",
        "bqsr_count_kernel",
        "bqsr_covariate_prepare",
        "bqsr_covariate_kernel",
        "bqsr_apply_prepare",
        "bqsr_apply_kernel",
    }
    assert set(timings) == expected
    for summary in timings.values():
        assert len(summary["samples"]) == repeats
        assert summary["p50_seconds"] > 0.0
        assert summary["p95_seconds"] >= summary["p50_seconds"]
    print(
        json.dumps(
            {
                "status": "pass",
                "suite": "fastgatk-kernel-benchmark",
                "schema_version": artifact["schema_version"],
                "backend": artifact["backend"],
                "simd_width": artifact["simd_width"],
                "compiled_simd_backend": artifact["compiled_simd_backend"],
                "compiled_simd_width": artifact["compiled_simd_width"],
                "pairhmm_kernel_p50_seconds": timings["pairhmm_kernel"]["p50_seconds"],
                "float_pairhmm_kernel_p50_seconds": timings["float_pairhmm_kernel"]["p50_seconds"],
                "float_pairhmm_simd_width": artifact["float_pairhmm_simd_width"],
                "float_pairhmm_max_abs_error": artifact["float_pairhmm_max_abs_error"],
                "pairhmm_likelihood_normalization_kernel_p50_seconds": timings[
                    "pairhmm_likelihood_normalization_kernel"]["p50_seconds"],
                "pairhmm_allele_marginalization_kernel_p50_seconds": timings[
                    "pairhmm_allele_marginalization_kernel"]["p50_seconds"],
                "pairhmm_read_allele_uncertainty_kernel_p50_seconds": timings[
                    "pairhmm_read_allele_uncertainty_kernel"]["p50_seconds"],
                "pairhmm_read_allele_best_kernel_p50_seconds": timings[
                    "pairhmm_read_allele_best_kernel"]["p50_seconds"],
                "flow_pairhmm_kernel_p50_seconds": timings["flow_pairhmm_kernel"]["p50_seconds"],
                "flow_pairhmm_execution_policy": artifact["flow_pairhmm_execution_policy"],
                "posterior_kernel_p50_seconds": timings["genotype_posterior_kernel"]["p50_seconds"],
                "joint_genotype_kernel_p50_seconds": timings["joint_genotype_kernel"]["p50_seconds"],
                "genotype_prior_kernel_p50_seconds": timings["genotype_prior_kernel"]["p50_seconds"],
                "reference_confidence_polyploid_kernel_p50_seconds": timings["reference_confidence_polyploid_kernel"]["p50_seconds"],
                "cohort_kernel_p50_seconds": timings["genotype_cohort_kernel"]["p50_seconds"],
                "cross_sample_reference_kernel_p50_seconds": timings["cross_sample_reference_kernel"]["p50_seconds"],
                "somatic_kernel_p50_seconds": timings["somatic_kernel"]["p50_seconds"],
                "somatic_posterior_kernel_p50_seconds": timings["somatic_posterior_kernel"]["p50_seconds"],
                "somatic_multiallelic_kernel_p50_seconds": timings["somatic_multiallelic_kernel"]["p50_seconds"],
                "bqsr_count_kernel_p50_seconds": timings["bqsr_count_kernel"]["p50_seconds"],
                "bqsr_count_execution_policy": artifact["timings"]["bqsr_count_execution_policy"],
                "bqsr_apply_kernel_p50_seconds": timings["bqsr_apply_kernel"]["p50_seconds"],
                "bqsr_apply_execution_policy": artifact["timings"]["bqsr_apply_execution_policy"],
                "bqsr_covariate_execution_policy": artifact["timings"][
                    "bqsr_covariate_execution_policy"],
                "bqsr_covariate_workspace_bytes": artifact["timings"][
                    "bqsr_covariate_workspace_bytes"],
                "bqsr_covariate_team_local_histogram": artifact["timings"][
                    "bqsr_covariate_team_local_histogram"],
                "persistent_cache_hits": artifact["persistent_pairhmm_cache_hits"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
