#!/usr/bin/env python3
"""Recompute and validate the versioned Fast-GATK progress snapshot."""

from __future__ import annotations

import collections
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCORE = ROOT / "progress_score.json"
REGISTRY = ROOT / "fastgatk-native/dispatcher/tool_registry.json"
FIXTURE_MANIFEST = ROOT / "fastgatk-native/tests/pinned_fixture_digests.sha256"

# Keep this map explicit instead of guessing from tool names: several GATK
# commands use short native aliases (HC, BQSR, PoN) and many have more than
# one bounded oracle.  A missing entry is a hard audit failure, so adding a
# registry tool cannot silently produce an un-evidenced progress report.
TOOL_EVIDENCE: dict[str, tuple[str, ...]] = {
    "Funcotator": (
        "verify_funcotator_vcf_gatk_oracle.py",
    ),
    "VariantAnnotator": (
        "verify_variant_annotator_resource_expression_gatk_oracle.py",
        "verify_variant_annotator_coverage_gatk_oracle.py",
    ),
    "HaplotypeCaller": (
        "verify_hc_broad_gatk_oracle.py", "verify_hc_streaming.py", "verify_hc_region_streaming.py",
        "verify_hc_assembly_region_boundary_gatk_oracle.py", "verify_hc_cigar_indel_activity_gatk_oracle.py",
        "verify_hc_complex_multiallelic_oracle.py",
        "verify_hc_multialt_owner_annotation_fixture_oracle.py", "verify_hc_polyploid_gatk_oracle.py",
        "verify_hc_gq_bands_gatk_oracle.py", "verify_hc_kmer_list_gatk_oracle.py",
        "verify_hc_min_pruning_gatk_oracle.py", "verify_hc_softclip_gatk_oracle.py",
        "verify_hc_bp_resolution_gatk_oracle.py", "verify_hc_rcm_pl_range_gatk_oracle.py",
        "verify_hc_variant_annotations.py", "verify_hc_sites_only_gatk_oracle.py",
        "verify_hc_chr17_69k_70k_gatk_oracle.py",
        "verify_flow_pairhmm_real_data_oracle.py",
    ),
    "BaseRecalibrator": ("verify_bqsr.py", "verify_bqsr_gatk_oracle.py", "verify_bqsr_indel_gatk_oracle.py", "verify_bqsr_read_filter_gatk_oracle.py", "verify_bqsr_long_read_gatk_oracle.py"),
    "ApplyBQSR": ("verify_bqsr.py", "verify_bqsr_gatk_oracle.py", "verify_bqsr_preserve_gatk_oracle.py", "verify_bqsr_cram_gatk_oracle.py", "verify_apply_bqsr_alias_gatk_oracle.py"),
    "GatherBQSRReports": ("verify_gather_bqsr_gatk_oracle.py",),
    "AnalyzeCovariates": ("verify_analyze_covariates.py", "verify_analyze_covariates_bqsr_alias_gatk_oracle.py"),
    "Mutect2": (
        "verify_mutect2.py", "verify_mutect2_gatk_oracle.py", "verify_somatic_likelihood_oracle.py",
        "verify_somatic_posterior_normal_oracle.py",
        "verify_mutect2_force_active_gatk_oracle.py", "verify_mutect2_normal_lod_gatk_oracle.py",
        "verify_mutect2_gvcf_reference_blocks_gatk_contract.py",
        "verify_mutect2_bp_resolution_gatk_contract.py",
        "verify_mutect2_gvcf_eventmap_gatk_contract.py",
        "verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py",
        "verify_mutect2_gvcf_matched_normal_gatk_contract.py",
        "verify_mutect2_gvcf_multitumor_gatk_contract.py",
        "verify_pairhmm_default_indel_quality_gatk_oracle.py",
        "verify_mutect2_dream_synthetic_oracle.py",
        "verify_mutect2_hcc1143_chr20_oracle.py",
    ),
    "CollectF1R2Counts": ("verify_collect_f1r2_counts.py",),
    "FilterMutectCalls": (
        "verify_filter_mutect_calls.py",
        "verify_filter_mutect_calls_contamination_oracle.py",
        "verify_filter_mutect_calls_germline_oracle.py",
        "verify_filter_mutect_orientation_gatk_oracle.py",
        "verify_filter_mutect_normal_artifact_oracle.py",
        "verify_filter_mutect_contamination_joint_oracle.py",
        "verify_filter_mutect_orientation_joint_oracle.py",
        "verify_filter_mutect_variant_index_alias_gatk_oracle.py",
        "verify_filter_mutect_dream_synthetic_joint_oracle.py",
        "verify_filter_mutect_hcc1143_joint_oracle.py",
    ),
    "LearnReadOrientationModel": (
        "verify_learn_read_orientation_model.py",
        "verify_learn_read_orientation_model_multisample_gatk_oracle.py",
    ),
    "AnnotateIntervals": ("verify_annotate_intervals.py",),
    "CountBasesInReference": ("verify_count_bases_in_reference.py",),
    "CompareReferences": ("verify_compare_references.py",),
    "CheckReferenceCompatibility": ("verify_check_reference_compatibility.py",),
    "FastaReferenceMaker": ("verify_fasta_reference_tools.py",),
    "FastaAlternateReferenceMaker": ("verify_fasta_reference_tools.py", "verify_fasta_alternate_iupac_hom_gatk_oracle.py"),
    "ShiftFasta": ("verify_shift_fasta.py",),
    "IndexFeatureFile": ("verify_index_feature_file.py",),
    "CountReads": ("verify_count_reads.py",),
    "FlagStat": ("verify_flag_stat.py",),
    "SplitIntervals": ("verify_split_intervals.py",),
    "FilterIntervals": ("verify_filter_intervals.py",),
    "PreprocessIntervals": ("verify_preprocess_intervals.py",),
    "CollectReadCounts": ("verify_collect_read_counts.py", "verify_collect_read_counts_gatk_oracle.py"),
    "DenoiseReadCounts": ("verify_denoise_read_counts.py", "verify_denoise_read_counts_interval_identity_gatk_oracle.py", "verify_denoise_read_counts_integer_input_gatk_oracle.py", "verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py"),
    "CreateReadCountPanelOfNormals": ("verify_create_read_count_panel_of_normals.py", "verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py", "verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py", "verify_pon_gatk_bundled_byte_structure_oracle.py"),
    "CallCopyRatioSegments": ("verify_call_copy_ratio_segments.py", "verify_call_copy_ratio_segments_gatk_oracle.py", "verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py", "verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py", "verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py", "verify_cnv_somatic_e2e_hcc1143_oracle.py"),
    "CollectAllelicCounts": ("verify_collect_allelic_counts.py", "verify_collect_allelic_counts_gatk_oracle.py"),
    "DepthOfCoverage": ("verify_depth_of_coverage.py", "verify_depth_of_coverage_multisample.py", "verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py"),
    "ModelSegments": ("verify_model_segments.py", "verify_model_segments_input_segments_gatk_oracle.py", "verify_model_segments_copy_ratio_conditionals_gatk_oracle.py", "verify_model_segments_allele_fraction_initialization_gatk_oracle.py", "verify_model_segments_allele_fraction_likelihood_gatk_oracle.py", "verify_model_segments_smoothing_gatk_oracle.py", "verify_model_segments_first_alt_fraction_gatk_oracle.py", "verify_model_segments_mcmc_oracle.py"),
    "GatherTranches": ("verify_gather_tranches.py",),
    "ApplyVQSR": (
        "verify_apply_vqsr.py", "verify_apply_vqsr_gatk_oracle.py",
        "verify_apply_vqsr_default_cutoff_gatk_oracle.py",
        "verify_apply_vqsr_sites_only_gatk_oracle.py",
        "verify_apply_vqsr_exclude_intervals_gatk_oracle.py",
        "verify_vqsr_scatter_joint_4shard_oracle.py",
    ),
    "VariantRecalibrator": (
        "verify_variant_recalibrator.py",
        "verify_variant_recalibrator_gatk_model_oracle.py",
        "verify_variant_recalibrator_annotation_order_gatk_oracle.py",
        "verify_variant_recalibrator_sample_every_gatk_oracle.py",
        "verify_variant_recalibrator_vbem_gatk_oracle.py",
    ),
    "GenotypeGVCFs": (
        "verify_genotype_gvcf.py",
        "verify_genotype_gvcf_assignment_gatk_oracle.py",
        "verify_genotype_gvcf_gp_input_gatk_oracle.py",
        "verify_genotype_gvcf_include_non_variant_gatk_oracle.py",
        "verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py",
        "verify_genotype_gvcf_multisample_reference_confidence_oracle.py",
        "verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py",
    ),
    "ReblockGVCF": (
        "verify_reblock_gvcf.py",
        "verify_reblock_gatk_oracle.py",
        "verify_reblock_gatk_shards.py",
        "verify_reblock_overlap_gatk_oracle.py",
    ),
    "SelectVariants": ("verify_select_variants.py", "verify_select_variants_gatk_oracle.py", "verify_select_variants_filtered_oracle.py", "verify_select_variants_sites_only_gatk_oracle.py"),
    "VariantsToTable": ("verify_variants_to_table.py", "verify_variants_to_table_gatk_oracle.py"),
    "VariantEval": ("verify_variant_eval.py", "verify_variant_eval_gatk_oracle.py", "verify_variant_eval_keep_ac0_gatk_oracle.py"),
    "GatherPileupSummaries": ("verify_gather_pileup_summaries.py", "verify_gather_pileup_gatk_oracle.py"),
    "CalculateContamination": ("verify_calculate_contamination.py", "verify_calculate_contamination_gatk_oracle.py"),
    "GetPileupSummaries": ("verify_get_pileup_summaries.py", "verify_get_pileup_gatk_oracle.py"),
    "ValidateVariants": ("verify_validate_variants.py", "verify_validate_variants_gatk_oracle.py", "verify_validate_variants_symbolic_gatk_oracle.py"),
    "GatherVcfs": ("verify_gather_vcfs.py", "verify_gather_vcfs_cli_boundary_gatk_oracle.py"),
    "LeftAlignAndTrimVariants": ("verify_left_align.py", "verify_left_align_gatk_oracle.py", "verify_left_align_sites_only_gatk_oracle.py", "verify_left_align_cli_boundary_gatk_oracle.py"),
    "VariantFiltration": ("verify_variant_filtration.py", "verify_variant_filtration_gatk_oracle.py", "verify_variant_filtration_missing_boolean_gatk_oracle.py", "verify_variant_filtration_set_nocall_gatk_oracle.py"),
    "SortSam": ("verify_sort_sam.py", "verify_sort_sam_gatk_oracle.py", "verify_sort_sam_cli_boundary_gatk_oracle.py"),
    "MarkDuplicates": ("verify_mark_duplicates.py", "verify_mark_duplicates_gatk_oracle.py", "verify_mark_duplicates_tagging_policy_gatk_oracle.py", "verify_mark_duplicates_pair_key_gatk_oracle.py"),
    "CombineGVCFs": (
        "verify_combine_gvcfs.py",
        "verify_combine_gvcfs_gatk_oracle.py",
        "verify_combine_gvcfs_plless_oracle.py",
        "verify_combine_gvcfs_interval_refblock_gatk_oracle.py",
    ),
    "GenomicsDBImport": ("verify_genomicsdb_import.py", "verify_genomicsdb_import_gatk_oracle.py", "verify_genomicsdb_import_sample_map_gatk_oracle.py", "verify_genomicsdb_import_native_interval_gatk_oracle.py", "verify_genomicsdb_bridge.py"),
}


def weighted(values: dict[str, float], weights: dict[str, float]) -> float:
    missing = set(weights) - set(values)
    if missing:
        raise SystemExit(f"missing score gates: {sorted(missing)}")
    return sum(values[name] * weights[name] for name in weights)


def _relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _fixture_digest(evidence: list[Path]) -> dict[str, object]:
    """Bind each audit row to the pinned fixture manifest and its evidence code.

    Most oracle fixtures are intentionally created in a temporary directory;
    hashing the manifest plus the exact verifier bytes gives the audit a
    stable, tool-specific provenance without checking generated temp files
    into the repository.
    """
    manifest_sha256 = _validate_fixture_manifest()
    hasher = hashlib.sha256()
    hasher.update(manifest_sha256.encode("ascii"))
    for path in evidence:
        hasher.update(_relative(path).encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(path.read_bytes())
    return {
        "algorithm": "sha256",
        "sha256": hasher.hexdigest(),
        "manifest": _relative(FIXTURE_MANIFEST),
        "manifest_sha256": manifest_sha256,
        "evidence_bytes": [_relative(path) for path in evidence],
    }


def _validate_fixture_manifest() -> str:
    """Validate every pinned fixture, not merely the manifest's existence."""
    if not FIXTURE_MANIFEST.is_file():
        raise SystemExit(f"missing fixture digest manifest: {FIXTURE_MANIFEST}")
    entries: list[tuple[str, str]] = []
    seen: set[str] = set()
    for line_number, raw in enumerate(FIXTURE_MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise SystemExit(f"invalid fixture digest manifest line {line_number}: {raw!r}")
        expected, relative = fields[0].lower(), fields[1].lstrip("*").strip()
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or relative in seen:
            raise SystemExit(f"invalid/duplicate fixture digest path at line {line_number}: {relative}")
        seen.add(relative)
        entries.append((expected, relative))
    if not entries:
        raise SystemExit("fixture digest manifest is empty")
    for expected, relative in entries:
        path = ROOT / relative
        if not path.is_file():
            raise SystemExit(f"missing pinned fixture: {relative}")
        hasher = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                hasher.update(chunk)
        if hasher.hexdigest() != expected:
            raise SystemExit(f"fixture digest mismatch: {relative}")
    return hashlib.sha256(FIXTURE_MANIFEST.read_bytes()).hexdigest()


def _fallback_reason(entry: dict[str, object]) -> str:
    notes = entry.get("notes", [])
    if not isinstance(notes, list):
        raise SystemExit("registry notes must be a list")
    for note in notes:
        if not isinstance(note, str):
            raise SystemExit("registry notes must contain strings")
        lowered = note.lower()
        if "fallback" in lowered or "remain" in lowered or "not fully" in lowered:
            return note
    status = entry.get("status", "unknown")
    return f"registry status={status}; no additional fallback boundary was recorded"


def _tool_audits(registry: dict[str, dict[str, object]], gate_weights: dict[str, float]) -> list[dict[str, object]]:
    registry_names = set(registry)
    mapped_names = set(TOOL_EVIDENCE)
    missing = sorted(registry_names - mapped_names)
    stale = sorted(mapped_names - registry_names)
    if missing or stale:
        raise SystemExit(f"tool evidence map mismatch: missing={missing}, stale={stale}")
    _validate_fixture_manifest()

    audits: list[dict[str, object]] = []
    for name, entry in registry.items():
        required = ("status", "aliases", "native_binary", "fallback_command", "fallback_policy", "backends")
        absent = [key for key in required if key not in entry]
        if absent:
            raise SystemExit(f"registry entry {name} missing audit fields: {absent}")
        if entry["status"] not in {"prototype", "contract-compatible", "adapter"}:
            raise SystemExit(f"registry entry {name} has invalid status: {entry['status']!r}")
        if not isinstance(entry["aliases"], list) or not entry["aliases"]:
            raise SystemExit(f"registry entry {name} has no aliases")
        if not isinstance(entry["native_binary"], str) or not entry["native_binary"]:
            raise SystemExit(f"registry entry {name} has no native_binary")
        if not isinstance(entry["fallback_command"], list) or not entry["fallback_command"]:
            raise SystemExit(f"registry entry {name} has no fallback_command")
        if entry["fallback_policy"] not in {"explicit", "required", "optional"}:
            raise SystemExit(f"registry entry {name} has invalid fallback_policy")
        if not isinstance(entry["backends"], list) or not entry["backends"]:
            raise SystemExit(f"registry entry {name} has no backends")
        evidence = [ROOT / "fastgatk-native/scripts" / relative for relative in TOOL_EVIDENCE[name]]
        missing_evidence = [_relative(path) for path in evidence if not path.is_file()]
        if missing_evidence:
            raise SystemExit(f"tool {name} missing evidence paths: {missing_evidence}")
        notes_text = " ".join(str(note) for note in entry.get("notes", []))
        has_oracle = any("oracle" in path.name or "gatk" in path.name for path in evidence)
        has_benchmark = "benchmark" in notes_text.lower() or any("benchmark" in path.name for path in evidence)
        has_kernel_telemetry = any(token in notes_text.lower() for token in ("kokkos", "telemetry", "resource", "hdf5"))
        has_output_contract = bool(entry.get("required_outputs") or entry.get("required_output_groups"))
        gate_scores = {
            "api_cli": 1.0,
            "oracle": 1.0 if has_oracle else 0.55,
            "format_sidecar": 1.0 if has_output_contract else 0.80,
            "resource_io": 0.75 if has_kernel_telemetry else 0.40,
            "e2e_perf": 0.75 if has_benchmark else 0.35,
        }
        for gate, value in gate_scores.items():
            if gate not in gate_weights or not 0.0 <= value <= 1.0:
                raise SystemExit(f"invalid tool gate score {name}.{gate}={value}")
        audits.append({
            "name": name,
            "status": entry["status"],
            "audit_status": "pass",
            "gate_scores": gate_scores,
            "score": round(weighted(gate_scores, gate_weights), 8),
            "evidence_paths": [_relative(path) for path in evidence],
            "fixture_digest": _fixture_digest(evidence),
            "fallback_reason": _fallback_reason(entry),
        })
    if len(audits) != len(TOOL_EVIDENCE):
        raise SystemExit(
            f"tool audit expected {len(TOOL_EVIDENCE)} entries, got {len(audits)}")
    return audits


def main() -> int:
    snapshot = json.loads(SCORE.read_text(encoding="utf-8"))
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))["tools"]
    gate_weights = snapshot["gate_weights"]
    workflow_weights = snapshot["workflow_weights"]
    tool_audits = _tool_audits(registry, gate_weights)
    if snapshot.get("tool_audit_schema_version") != 1:
        raise SystemExit("unsupported or missing tool audit schema version")
    if snapshot.get("tool_audit_expected_entries") != len(tool_audits):
        raise SystemExit("tool audit entry-count snapshot mismatch")
    computed_workflows: dict[str, float] = {}
    for name, workflow in snapshot["workflows"].items():
        score = weighted(workflow["gates"], gate_weights)
        computed_workflows[name] = score
        if abs(score - workflow["score"]) > 1.0e-12:
            raise SystemExit(f"workflow score mismatch for {name}: {score} != {workflow['score']}")
    global_score = weighted(computed_workflows, workflow_weights)
    if abs(global_score - snapshot["global_score"]) > 1.0e-12:
        raise SystemExit(f"global score mismatch: {global_score} != {snapshot['global_score']}")
    status = collections.Counter(entry.get("status") for entry in registry.values())
    expected = snapshot["registry_snapshot"]
    actual = {
        "entries": len(registry),
        "prototype": status.get("prototype", 0),
        "contract_compatible": status.get("contract-compatible", 0),
        "adapter": status.get("adapter", 0),
    }
    for key, value in actual.items():
        if value != expected[key]:
            raise SystemExit(f"registry snapshot mismatch for {key}: {value} != {expected[key]}")
    fraction = actual["contract_compatible"] / actual["entries"]
    if abs(fraction - expected["contract_compatible_fraction"]) > 1.0e-12:
        raise SystemExit("registry contract-compatible fraction mismatch")
    print(json.dumps({
        "status": "pass",
        "method": snapshot["method"],
        "global_score": round(global_score, 8),
        "global_percent_rounded": round(global_score * 100),
        "registry": actual,
        "tool_audit": {
            "schema_version": 1,
            "status": "pass",
            "entries": len(tool_audits),
            "fields": ["status", "gate_scores", "evidence_paths", "fixture_digest", "fallback_reason"],
            "mechanical_gate_scores": True,
        },
        "tool_audits": tool_audits,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
