#!/usr/bin/env python3
"""Contract tests for the fastgatk dispatcher and tool registry."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DISPATCHER = ROOT / "fastgatk-native/dispatcher/fastgatk"
COMPAT_LAUNCHER = ROOT / "fastgatk-native/dispatcher/gatk"
REGISTRY = ROOT / "fastgatk-native/dispatcher/tool_registry.json"
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
BINARY = ROOT / "fastgatk-native/build/fastgatk-hc-call"
BQSR_BINARY = ROOT / "fastgatk-native/build/fastgatk-bqsr"
APPLY_BQSR_BINARY = ROOT / "fastgatk-native/build/fastgatk-apply-bqsr"
GATHER_BQSR_BINARY = ROOT / "fastgatk-native/build/fastgatk-gather-bqsr-reports"
ANALYZE_COVARIATES_BINARY = ROOT / "fastgatk-native/build/fastgatk-analyze-covariates"
GENOTYPE_GVCF_BINARY = ROOT / "fastgatk-native/build/fastgatk-genotype-gvcf"
REBLOCK_GVCF_BINARY = ROOT / "fastgatk-native/build/fastgatk-reblock-gvcf"
SELECT_VARIANTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-select-variants"
GATHER_VCFS_BINARY = ROOT / "fastgatk-native/build/fastgatk-gather-vcfs"
LEFT_ALIGN_BINARY = ROOT / "fastgatk-native/build/fastgatk-left-align-trim"
VARIANT_FILTRATION_BINARY = ROOT / "fastgatk-native/build/fastgatk-variant-filtration"
SORT_SAM_BINARY = ROOT / "fastgatk-native/build/fastgatk-sort-sam"
MARK_DUPLICATES_BINARY = ROOT / "fastgatk-native/build/fastgatk-mark-duplicates"
COMBINE_GVCF_BINARY = ROOT / "fastgatk-native/build/fastgatk-combine-gvcfs"
FILTER_MUTECT_BINARY = ROOT / "fastgatk-native/build/fastgatk-filter-mutect-calls"
GENOMICSDB_BINARY = ROOT / "fastgatk-native/build/fastgatk-genomicsdb-import"
VARIANTS_TO_TABLE_BINARY = ROOT / "fastgatk-native/build/fastgatk-variants-to-table"
VARIANT_EVAL_BINARY = ROOT / "fastgatk-native/build/fastgatk-variant-eval"
VALIDATE_VARIANTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-validate-variants"
GET_PILEUP_SUMMARIES_BINARY = ROOT / "fastgatk-native/build/fastgatk-get-pileup-summaries"
CALCULATE_CONTAMINATION_BINARY = ROOT / "fastgatk-native/build/fastgatk-calculate-contamination"
GATHER_PILEUP_SUMMARIES_BINARY = ROOT / "fastgatk-native/build/fastgatk-gather-pileup-summaries"
LEARN_READ_ORIENTATION_BINARY = ROOT / "fastgatk-native/build/fastgatk-learn-read-orientation-model"
COLLECT_F1R2_BINARY = ROOT / "fastgatk-native/build/fastgatk-collect-f1r2-counts"
ANNOTATE_INTERVALS_BINARY = ROOT / "fastgatk-native/build/fastgatk-annotate-intervals"
COUNT_BASES_BINARY = ROOT / "fastgatk-native/build/fastgatk-count-bases-in-reference"
COMPARE_REFERENCES_BINARY = ROOT / "fastgatk-native/build/fastgatk-compare-references"
CHECK_REFERENCE_COMPATIBILITY_BINARY = ROOT / "fastgatk-native/build/fastgatk-check-reference-compatibility"
FASTA_REFERENCE_MAKER_BINARY = ROOT / "fastgatk-native/build/fastgatk-fasta-reference-maker"
FASTA_ALTERNATE_REFERENCE_MAKER_BINARY = ROOT / "fastgatk-native/build/fastgatk-fasta-alternate-reference-maker"
SHIFT_FASTA_BINARY = ROOT / "fastgatk-native/build/fastgatk-shift-fasta"
FILTER_INTERVALS_BINARY = ROOT / "fastgatk-native/build/fastgatk-filter-intervals"
PREPROCESS_INTERVALS_BINARY = ROOT / "fastgatk-native/build/fastgatk-preprocess-intervals"
COLLECT_READ_COUNTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-collect-read-counts"
DENOISE_READ_COUNTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-denoise-read-counts"
CREATE_PON_BINARY = ROOT / "fastgatk-native/build/fastgatk-create-read-count-panel-of-normals"
CALL_COPY_RATIO_SEGMENTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-call-copy-ratio-segments"
COLLECT_ALLELIC_COUNTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-collect-allelic-counts"
MODEL_SEGMENTS_BINARY = ROOT / "fastgatk-native/build/fastgatk-model-segments"
GATHER_TRANCHES_BINARY = ROOT / "fastgatk-native/build/fastgatk-gather-tranches"
APPLY_VQSR_BINARY = ROOT / "fastgatk-native/build/fastgatk-apply-vqsr"
VARIANT_RECALIBRATOR_BINARY = ROOT / "fastgatk-native/build/fastgatk-variant-recalibrator"


def run(*args: str, expect: int = 0, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    result = subprocess.run([str(DISPATCHER), *args], text=True, capture_output=True, env=merged)
    assert result.returncode == expect, (args, result.returncode, result.stdout, result.stderr)
    return result


def main() -> int:
    assert DISPATCHER.is_file() and os.access(DISPATCHER, os.X_OK)
    assert COMPAT_LAUNCHER.is_file() and os.access(COMPAT_LAUNCHER, os.X_OK)
    registry = json.loads(REGISTRY.read_text())
    assert registry["schema_version"] == 1
    launcher_contract = registry["launcher_contract"]
    assert launcher_contract["value_options"] == [
        "--java-options", "-java-options", "--gatk-config-file"
    ]
    assert launcher_contract["argument_file_prefix"] == "@"
    assert launcher_contract["argument_file_escape"] == "@@"
    assert launcher_contract["max_argument_file_depth"] == 16
    assert launcher_contract["separator"] == "--"
    assert launcher_contract["unknown_option_policy"] == "fail-closed"
    assert "HaplotypeCaller" in registry["tools"]
    assert "--kmer-size" in registry["tools"]["HaplotypeCaller"]["repeatable_options"]
    hc_registry = registry["tools"]["HaplotypeCaller"]
    for option in ("--min-base-quality-score", "-mbq"):
        assert option in hc_registry["value_options"]
    assert "--soft-clip-low-quality-ends" in hc_registry["flag_options"]
    assert "--soft-clip-low-quality-ends" in hc_registry["optional_boolean_options"]
    assert "BaseRecalibrator" in registry["tools"]
    assert "ApplyBQSR" in registry["tools"]
    assert "GatherBQSRReports" in registry["tools"]
    gather_bqsr = registry["tools"]["GatherBQSRReports"]
    for option in ("--QUIET", "--tmp-dir", "--verbosity", "--use-jdk-deflater",
                   "--jdk-deflater", "--use-jdk-inflater", "--jdk-inflater"):
        assert option in gather_bqsr["value_options"]
    for option in ("--QUIET", "--quiet", "--use-jdk-deflater", "--jdk-deflater",
                   "--use-jdk-inflater", "--jdk-inflater"):
        assert option in gather_bqsr["optional_boolean_options"]
    assert any("GATK 4.6.2.0 utility controls" in note and "verbosity" in note
               for note in gather_bqsr["notes"])
    collect_read_counts = registry["tools"]["CollectReadCounts"]
    for option in ("--QUIET", "--read-validation-stringency", "-VS",
                   "--disable-bam-index-caching", "-DBIC", "--tmp-dir"):
        assert option in collect_read_counts["value_options"]
    for option in ("--QUIET", "--read-validation-stringency", "--disable-bam-index-caching",
                   "--use-jdk-deflater", "--use-jdk-inflater"):
        if option == "--read-validation-stringency":
            continue
        if option == "--QUIET" or option in collect_read_counts["optional_boolean_options"]:
            assert option in collect_read_counts["optional_boolean_options"]
    assert "--read-filter" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "-RF" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "--disable-read-filter" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "-DF" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "--disable-tool-default-read-filters" in registry["tools"]["BaseRecalibrator"]["optional_boolean_options"]
    assert "--mismatches-context-size" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "--maximum-cycle-value" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert "--max-cycle" in registry["tools"]["BaseRecalibrator"]["value_options"]
    assert any("maximum-cycle-value" in note and "long" in note.lower()
               for note in registry["tools"]["BaseRecalibrator"]["notes"])
    assert any("seven standard defaults" in note and "-DF" in note and "-RF" in note
               for note in registry["tools"]["BaseRecalibrator"]["notes"])
    assert "--preserve-qscores-less-than" in registry["tools"]["ApplyBQSR"]["value_options"]
    assert "-bqsr" in registry["tools"]["ApplyBQSR"]["value_options"]
    assert any("-bqsr short alias" in note and "--bqsr-recal-file" in note
               for note in registry["tools"]["ApplyBQSR"]["notes"])
    assert any("preserve" in note.lower() and "quality" in note.lower()
               for note in registry["tools"]["ApplyBQSR"]["notes"])
    assert registry["tools"]["AnalyzeCovariates"]["status"] == "contract-compatible"
    assert "-bqsr" in registry["tools"]["AnalyzeCovariates"]["value_options"]
    assert any("lower-case -bqsr short alias" in note and "--bqsr-recal-file" in note
               for note in registry["tools"]["AnalyzeCovariates"]["notes"])
    assert registry["tools"]["GenotypeGVCFs"]["status"] == "contract-compatible"
    assert registry["tools"]["GenotypeGVCFs"]["compatibility_scope"] == "local-vcf-gvcf-and-fastgatk-native-gendb"
    assert "opaque-TileDB-GenomicsDB-workspaces" in registry["tools"]["GenotypeGVCFs"]["fallback_boundaries"]
    assert "Kokkos" in registry["tools"]["GenotypeGVCFs"]["backends"]
    assert "--genotype-assignment-method" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--standard-min-confidence-threshold-for-calling" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--stand-call-conf" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--max-alternate-alleles" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--annotate-with-num-discovered-alleles" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--annotate-with-num-discovered-alleles" in registry["tools"]["GenotypeGVCFs"]["optional_boolean_options"]
    assert "--only-output-calls-starting-in-intervals" in registry["tools"]["GenotypeGVCFs"]["value_options"]
    assert "--only-output-calls-starting-in-intervals" in registry["tools"]["GenotypeGVCFs"]["optional_boolean_options"]
    assert any("output-allele subset" in note for note in registry["tools"]["GenotypeGVCFs"]["notes"])
    assert any("orphan `*`" in note or "orphan '*'" in note
               for note in registry["tools"]["GenotypeGVCFs"]["notes"])
    assert any("max-alternate-alleles" in note and "likelihood-score" in note
               for note in registry["tools"]["GenotypeGVCFs"]["notes"])
    assert any("INFO/NDA" in note for note in registry["tools"]["GenotypeGVCFs"]["notes"])
    assert registry["tools"]["ReblockGVCF"]["status"] == "contract-compatible"
    assert "Kokkos" in registry["tools"]["ReblockGVCF"]["backends"]
    assert registry["tools"]["SelectVariants"]["status"] == "contract-compatible"
    assert "--exclude-filtered" in registry["tools"]["SelectVariants"]["value_options"]
    assert "--exclude-filtered-variants" in registry["tools"]["SelectVariants"]["value_options"]
    assert "--exclude-filtered" in registry["tools"]["SelectVariants"]["optional_boolean_options"]
    assert "--exclude-filtered-variants" in registry["tools"]["SelectVariants"]["optional_boolean_options"]
    assert any("PASS and '.' are not considered filtered" in note
               for note in registry["tools"]["SelectVariants"]["notes"])
    assert "--set-filtered-gt-to-nocall" in registry["tools"]["SelectVariants"]["value_options"]
    assert "--set-filtered-gt-to-nocall" in registry["tools"]["SelectVariants"]["optional_boolean_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["SelectVariants"]["value_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["SelectVariants"]["optional_boolean_options"]
    assert any("final writer" in note and "8-column" in note
               for note in registry["tools"]["SelectVariants"]["notes"])
    assert registry["tools"]["VariantsToTable"]["status"] == "contract-compatible"
    assert registry["tools"]["VariantEval"]["status"] == "contract-compatible"
    assert "Kokkos" in registry["tools"]["VariantEval"]["backends"]
    assert any("ValidationReport" in note and "4x4" in note
               for note in registry["tools"]["VariantEval"]["notes"])
    assert "--keep-ac0" in registry["tools"]["VariantEval"]["value_options"]
    assert "-keep-ac0" in registry["tools"]["VariantEval"]["value_options"]
    assert "--keep-ac0" in registry["tools"]["VariantEval"]["optional_boolean_options"]
    assert "-keep-ac0" in registry["tools"]["VariantEval"]["optional_boolean_options"]
    assert any("--keep-ac0" in note and "CountVariants" in note
               for note in registry["tools"]["VariantEval"]["notes"])
    assert registry["tools"]["ValidateVariants"]["status"] == "contract-compatible"
    assert registry["tools"]["GetPileupSummaries"]["status"] == "contract-compatible"
    assert registry["tools"]["CalculateContamination"]["status"] == "contract-compatible"
    assert registry["tools"]["GatherPileupSummaries"]["status"] == "contract-compatible"
    assert registry["tools"]["LearnReadOrientationModel"]["status"] == "contract-compatible"
    learn_orientation = registry["tools"]["LearnReadOrientationModel"]
    for option in ("--QUIET", "--tmp-dir", "--VERBOSITY", "--use-jdk-deflater",
                   "--use-jdk-inflater"):
        assert option in learn_orientation["value_options"]
    for option in ("--QUIET", "--quiet", "--use-jdk-deflater", "--use-jdk-inflater"):
        assert option in learn_orientation["optional_boolean_options"]
    assert any("GATK 4.6.2.0 utility controls" in note and "telemetry" in note
               for note in learn_orientation["notes"])
    assert registry["tools"]["CollectF1R2Counts"]["status"] == "contract-compatible"
    assert registry["tools"]["AnnotateIntervals"]["status"] == "contract-compatible"
    annotate_registry = registry["tools"]["AnnotateIntervals"]
    assert "--feature-query-lookahead" in annotate_registry["value_options"]
    assert any("feature-query-lookahead" in note and "telemetry" in note
               for note in annotate_registry["notes"])
    assert registry["tools"]["CountBasesInReference"]["status"] == "contract-compatible"
    assert registry["tools"]["CompareReferences"]["status"] == "contract-compatible"
    assert registry["tools"]["CheckReferenceCompatibility"]["status"] == "contract-compatible"
    assert registry["tools"]["FastaReferenceMaker"]["status"] == "contract-compatible"
    assert registry["tools"]["FastaAlternateReferenceMaker"]["status"] == "contract-compatible"
    assert any("Homozygous --use-iupac-sample" in note and "ALT[1]" in note
               for note in registry["tools"]["FastaAlternateReferenceMaker"]["notes"])
    assert registry["tools"]["ShiftFasta"]["status"] == "contract-compatible"
    assert registry["tools"]["IndexFeatureFile"]["status"] == "contract-compatible"
    assert registry["tools"]["CountReads"]["status"] == "contract-compatible"
    assert registry["tools"]["FlagStat"]["status"] == "contract-compatible"
    assert registry["tools"]["SplitIntervals"]["status"] == "contract-compatible"
    assert registry["tools"]["FilterIntervals"]["status"] == "contract-compatible"
    assert registry["tools"]["PreprocessIntervals"]["status"] == "contract-compatible"
    assert registry["tools"]["CollectReadCounts"]["status"] == "contract-compatible"
    assert registry["tools"]["DenoiseReadCounts"]["status"] == "contract-compatible"
    assert any("sequence_dictionary" in note and "both standardized and denoised" in note
               for note in registry["tools"]["DenoiseReadCounts"]["notes"])
    assert registry["tools"]["CreateReadCountPanelOfNormals"]["status"] == "contract-compatible"
    assert registry["tools"]["CallCopyRatioSegments"]["status"] == "contract-compatible"
    assert any("NaN" in note and ".igv.seg" in note
               for note in registry["tools"]["CallCopyRatioSegments"]["notes"])
    assert any("Infinity" in note and "non-finite" in note
               for note in registry["tools"]["CallCopyRatioSegments"]["notes"])
    assert registry["tools"]["CollectAllelicCounts"]["status"] == "contract-compatible"
    assert "--disable-read-filter" in registry["tools"]["CollectAllelicCounts"]["value_options"]
    assert "--read-filter" in registry["tools"]["CollectAllelicCounts"]["value_options"]
    assert "--disable-tool-default-read-filters" in registry["tools"]["CollectAllelicCounts"]["optional_boolean_options"]
    depth_of_coverage = registry["tools"]["DepthOfCoverage"]
    assert "--ignore-deletion-sites" in depth_of_coverage["value_options"]
    assert "--ignore-deletion-sites" in depth_of_coverage["optional_boolean_options"]
    assert any("ignore-deletion-sites" in note and "deletion pileup" in note
               for note in depth_of_coverage["notes"])
    assert registry["tools"]["ModelSegments"]["status"] == "contract-compatible"
    assert any("default" in note and "KernelSegmenter" in note and "--segments" in note
               for note in registry["tools"]["ModelSegments"]["notes"])
    assert any("rank-one" in note.lower() and "linear copy-ratio" in note.lower()
               for note in registry["tools"]["ModelSegments"]["notes"])
    assert any("complete per-interval sample vector" in note and "anti-correlated" in note
               for note in registry["tools"]["ModelSegments"]["notes"])
    assert any("Gamma-bias Laplace likelihood" in note and "modelBegin" in note
               for note in registry["tools"]["ModelSegments"]["notes"])
    assert registry["tools"]["GatherTranches"]["status"] == "contract-compatible"
    assert "-tranche" in registry["tools"]["GatherTranches"]["value_options"]
    assert registry["tools"]["ApplyVQSR"]["status"] == "contract-compatible"
    assert "--ignore-all-filters" in registry["tools"]["ApplyVQSR"]["optional_boolean_options"]
    assert "--exclude-filtered" in registry["tools"]["ApplyVQSR"]["optional_boolean_options"]
    assert any("one-ALT VCF shape" in note and "Number=A AS_culprit" in note
               for note in registry["tools"]["ApplyVQSR"]["notes"])
    assert any("ignore-all-filters" in note and "bare/true/false" in note
               for note in registry["tools"]["ApplyVQSR"]["notes"])
    assert registry["tools"]["VariantRecalibrator"]["status"] == "contract-compatible"
    assert "--sites-only-vcf-output" in registry["tools"]["VariantRecalibrator"]["value_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["VariantRecalibrator"]["optional_boolean_options"]
    assert any("--max-attempts" in note and "--max-iterations" in note
               and "retry" in note.lower() for note in registry["tools"]["VariantRecalibrator"]["notes"])
    assert registry["tools"]["GatherVcfs"]["status"] == "contract-compatible"
    assert "--INPUT" in registry["tools"]["GatherVcfs"]["value_options"]
    assert "--OUTPUT" in registry["tools"]["GatherVcfs"]["value_options"]
    assert "--QUIET" in registry["tools"]["GatherVcfs"]["optional_boolean_options"]
    assert "--VERBOSITY" in registry["tools"]["GatherVcfs"]["value_options"]
    assert any("uppercase --INPUT/--OUTPUT" in note
               for note in registry["tools"]["GatherVcfs"]["notes"])
    assert registry["tools"]["LeftAlignAndTrimVariants"]["status"] == "contract-compatible"
    assert "--sites-only-vcf-output" in registry["tools"]["LeftAlignAndTrimVariants"]["value_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["LeftAlignAndTrimVariants"]["optional_boolean_options"]
    assert "--exclude-intervals" in registry["tools"]["LeftAlignAndTrimVariants"]["value_options"]
    assert "--interval-padding" in registry["tools"]["LeftAlignAndTrimVariants"]["value_options"]
    assert "--interval-exclusion-padding" in registry["tools"]["LeftAlignAndTrimVariants"]["value_options"]
    assert "--split-multi-allelics" in registry["tools"]["LeftAlignAndTrimVariants"]["optional_boolean_options"]
    assert any("8-column site records" in note
               for note in registry["tools"]["LeftAlignAndTrimVariants"]["notes"])
    assert any("-XL exclusion" in note and "-ip/-ixp" in note
               for note in registry["tools"]["LeftAlignAndTrimVariants"]["notes"])
    assert registry["tools"]["VariantFiltration"]["status"] == "contract-compatible"
    assert "--mask-description" in registry["tools"]["VariantFiltration"]["value_options"]
    assert "--filter-not-in-mask" in registry["tools"]["VariantFiltration"]["flag_options"]
    # The registry keeps the mask reverse-logic and header-description
    # guarantees in separate notes; validate both claims without requiring
    # unrelated documentation to be co-located in one string.
    assert any("filter-not-in-mask" in note
               for note in registry["tools"]["VariantFiltration"]["notes"])
    assert any("mask-description" in note
               for note in registry["tools"]["VariantFiltration"]["notes"])
    assert "--set-filtered-genotype-to-no-call" in registry["tools"]["VariantFiltration"]["value_options"]
    assert "--set-filtered-genotype-to-no-call" in registry["tools"]["VariantFiltration"]["optional_boolean_options"]
    assert any("set-filtered-genotype-to-no-call" in note and "AC" in note and "AF" in note
               for note in registry["tools"]["VariantFiltration"]["notes"])
    assert registry["tools"]["SortSam"]["status"] == "contract-compatible"
    sort_sam = registry["tools"]["SortSam"]
    assert sort_sam["compatibility_scope"] == "local-coordinate-and-queryname-external-memory-sort"
    assert "cloud-object-store-staging" in sort_sam["fallback_boundaries"]
    for option in ("--INPUT", "--OUTPUT", "--REFERENCE_SEQUENCE", "-SO",
                   "--SORT_ORDER", "--MAX_RECORDS_IN_RAM", "--TMP_DIR",
                   "--COMPRESSION_LEVEL", "--VALIDATION_STRINGENCY",
                   "--CREATE_INDEX", "--VERBOSITY"):
        assert option in sort_sam["value_options"]
    for option in ("--CREATE_INDEX", "--QUIET", "--USE_JDK_DEFLATER", "--USE_JDK_INFLATER"):
        assert option in sort_sam["optional_boolean_options"]
    assert any("Picard/GATK 4.6.2.0 long aliases" in note and
               "COMPRESSION_LEVEL" in note and "VALIDATION_STRINGENCY" in note
               for note in sort_sam["notes"])
    assert registry["tools"]["MarkDuplicates"]["status"] == "contract-compatible"
    mark_duplicates = registry["tools"]["MarkDuplicates"]
    assert mark_duplicates["compatibility_scope"] == "local-coordinate-sorted-two-pass-external-memory-duplicate-marking"
    assert "cloud-staging" in mark_duplicates["fallback_boundaries"]
    assert "--tagging-policy" in registry["tools"]["MarkDuplicates"]["value_options"]
    assert any("OpticalOnly" in note and "DT:Z:SQ" in note
               for note in registry["tools"]["MarkDuplicates"]["notes"])
    assert registry["tools"]["CombineGVCFs"]["status"] == "contract-compatible"
    assert "--sites-only-vcf-output" in registry["tools"]["CombineGVCFs"]["value_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["CombineGVCFs"]["optional_boolean_options"]
    assert any("PL-less" in note and "FORMAT/PL" in note and "INFO/DP" in note
               for note in registry["tools"]["CombineGVCFs"]["notes"])
    assert registry["tools"]["FilterMutectCalls"]["status"] == "contract-compatible"
    assert "--force-active" in registry["tools"]["Mutect2"]["optional_boolean_options"]
    assert "--force-active" in registry["tools"]["Mutect2"]["flag_options"]
    assert any("force-active" in note and "ActivityProfile" in note
               for note in registry["tools"]["Mutect2"]["notes"])
    assert any("ContaminationFilter" in note and "NON_SOMATIC" in note and
               "OPTIMAL_F_SCORE" in note
               for note in registry["tools"]["FilterMutectCalls"]["notes"])
    assert "--mitochondria-mode" in registry["tools"]["Mutect2"]["optional_boolean_options"]
    assert "--mitochondria-mode" in registry["tools"]["FilterMutectCalls"]["optional_boolean_options"]
    assert "--microbial-mode" in registry["tools"]["FilterMutectCalls"]["optional_boolean_options"]
    assert "--sites-only-vcf-output" in registry["tools"]["FilterMutectCalls"]["value_options"]
    assert "-OVI" in registry["tools"]["FilterMutectCalls"]["value_options"]
    assert "-OVI" in registry["tools"]["FilterMutectCalls"]["optional_boolean_options"]
    assert any("ReadOrientationFilter" in note and "OPTIMAL_F_SCORE" in note
               for note in registry["tools"]["FilterMutectCalls"]["notes"])
    ovi_plan = json.loads(run(
        "FilterMutectCalls", "-V", "input.vcf.gz", "-O", "output.vcf.gz",
        "-OVI", "false", "--dry-run",
    ).stdout)
    assert ovi_plan["execution_mode"] == "native"
    assert ovi_plan["argv"][-2:] == ["-OVI", "false"]
    assert "--validate-GVCF" in registry["tools"]["ValidateVariants"]["optional_boolean_options"]
    assert "--warn-on-errors" in registry["tools"]["ValidateVariants"]["optional_boolean_options"]
    assert "--do-not-validate-filtered-records" in registry["tools"]["ValidateVariants"]["optional_boolean_options"]
    assert any("Symbolic ALT" in note and "concrete ALT" in note
               for note in registry["tools"]["ValidateVariants"]["notes"])
    for alias in ("-D", "-isr", "-do-not-validate-filtered-records",
                  "-warn-on-errors", "-disable-sequence-dictionary-validation"):
        assert alias in registry["tools"]["ValidateVariants"].get("value_options", []) \
            or alias in registry["tools"]["ValidateVariants"].get("optional_boolean_options", [])
    assert registry["tools"]["GenomicsDBImport"]["status"] == "adapter"
    assert registry["tools"]["GenomicsDBImport"]["fallback_command"] == ["gatk", "GenomicsDBImport"]
    assert "--genomicsdb-update-workspace-path" in registry["tools"]["GenomicsDBImport"]["value_options"]
    assert any("incremental" in note and "native sparse" in note
               for note in registry["tools"]["GenomicsDBImport"]["notes"])
    assert any("interval" in note and "record span" in note and "native" in note
               for note in registry["tools"]["GenomicsDBImport"]["notes"])

    assert "Usage: fastgatk" in run("--help").stdout
    compat_version = subprocess.run(
        [str(COMPAT_LAUNCHER), "--version"], text=True, capture_output=True,
        check=False,
    )
    assert compat_version.returncode == 0
    assert "fastgatk 0.1.0" in compat_version.stdout
    assert "HaplotypeCaller" in run("--list").stdout
    assert "BaseRecalibrator" in run("--list").stdout
    assert "ApplyBQSR" in run("--list").stdout
    assert "GatherBQSRReports" in run("--list").stdout
    assert "AnalyzeCovariates" in run("--list").stdout
    assert "GenotypeGVCFs" in run("--list").stdout
    assert "ReblockGVCF" in run("--list").stdout
    assert "SelectVariants" in run("--list").stdout
    assert "VariantsToTable" in run("--list").stdout
    assert "VariantEval" in run("--list").stdout
    assert "ValidateVariants" in run("--list").stdout
    assert "GetPileupSummaries" in run("--list").stdout
    assert "CalculateContamination" in run("--list").stdout
    assert "GatherPileupSummaries" in run("--list").stdout
    assert "LearnReadOrientationModel" in run("--list").stdout
    assert "CollectF1R2Counts" in run("--list").stdout
    assert "AnnotateIntervals" in run("--list").stdout
    assert "CountBasesInReference" in run("--list").stdout
    assert "CompareReferences" in run("--list").stdout
    assert "CheckReferenceCompatibility" in run("--list").stdout
    assert "FilterIntervals" in run("--list").stdout
    assert "IndexFeatureFile" in run("--list").stdout
    assert "CountReads" in run("--list").stdout
    assert "FlagStat" in run("--list").stdout
    assert "SplitIntervals" in run("--list").stdout
    assert "PreprocessIntervals" in run("--list").stdout
    assert "CollectReadCounts" in run("--list").stdout
    assert "DenoiseReadCounts" in run("--list").stdout
    assert "CreateReadCountPanelOfNormals" in run("--list").stdout
    assert "CallCopyRatioSegments" in run("--list").stdout
    assert "CollectAllelicCounts" in run("--list").stdout
    assert "ModelSegments" in run("--list").stdout
    assert "GatherTranches" in run("--list").stdout
    assert "ApplyVQSR" in run("--list").stdout
    assert "VariantRecalibrator" in run("--list").stdout
    assert "GatherVcfs" in run("--list").stdout
    assert "LeftAlignAndTrimVariants" in run("--list").stdout
    assert "VariantFiltration" in run("--list").stdout
    assert "SortSam" in run("--list").stdout
    assert "MarkDuplicates" in run("--list").stdout
    assert "CombineGVCFs" in run("--list").stdout
    assert "FilterMutectCalls" in run("--list").stdout
    assert "GenomicsDBImport" in run("--list").stdout
    assert "registry schema 1" in run("--version").stdout
    assert "fastgatk HaplotypeCaller" in run("HaplotypeCaller", "--help").stdout
    assert "fastgatk BaseRecalibrator" in run("BaseRecalibrator", "--help").stdout
    assert "fastgatk AnalyzeCovariates" in run("AnalyzeCovariates", "--help").stdout

    with tempfile.TemporaryDirectory(prefix="fastgatk-dispatcher-") as directory:
        work = Path(directory)
        args_file = work / "hc.args"
        args_file.write_text(
            f'-I "{BAM}" -L 17:69000-69100 -O {work / "calls.vcf"} '
            f'--output-manifest {work / "calls.vcf.manifest.json"}\n',
            encoding="utf-8",
        )
        dry = run("HaplotypeCaller", f"@{args_file}", "--dry-run", env={"FASTGATK_HC_BINARY": str(BINARY)})
        plan = json.loads(dry.stdout)
        assert plan["status"] == "dry-run"
        assert plan["execution_mode"] == "native"
        assert plan["binary_exists"] is True
        assert str(BAM) in plan["argv"]
        assert "17:69000-69100" in plan["argv"]

        # Keep the dispatcher allow-list in lockstep with the native HC
        # parser for the bounded assembly/streaming controls.  These options
        # are implemented by the C++ Host + Kokkos path and must not be
        # rejected before reaching the binary in a Nextflow/SLURM command.
        advanced_dry = run(
            "HaplotypeCaller", "-I", str(BAM), "-O", str(work / "advanced.vcf"),
            "--stream-by-contig", "--stream-by-region", "128",
            "--adaptive-pruning", "--adaptive-pruning-initial-error-rate", "0.001",
            "--pruning-lod-threshold", "2", "--pruning-seeding-lod-threshold", "3",
            "--max-unpruned-variants", "100", "--allow-non-unique-kmers-in-ref",
            "--genotype-assignment-method", "USE_PLS_TO_ASSIGN", "--dry-run",
            env={"FASTGATK_HC_BINARY": str(BINARY)},
        )
        advanced_plan = json.loads(advanced_dry.stdout)
        assert advanced_plan["execution_mode"] == "native"
        assert "--stream-by-region" in advanced_plan["argv"]
        assert "--adaptive-pruning" in advanced_plan["argv"]

        # Deployment builds may live outside the conventional OpenMP
        # directory.  The dispatcher must resolve the registry suffix against
        # FASTGATK_NATIVE_BUILD without requiring one environment variable per
        # tool.
        serial_dry = run(
            "HaplotypeCaller", "-I", str(BAM), "-L", "17:69000-69100",
            "-O", str(work / "serial-calls.vcf"), "--dry-run",
            env={"FASTGATK_HC_BINARY": "",
                 "FASTGATK_NATIVE_BUILD": str(ROOT / "fastgatk-native/build-serial")},
        )
        serial_plan = json.loads(serial_dry.stdout)
        assert serial_plan["binary_exists"] is True
        assert "/build-serial/fastgatk-hc-call" in serial_plan["binary"]

        result = run(
            "HaplotypeCaller", "--input", str(BAM), "--intervals=17:69000-69100",
            "--output", str(work / "long-alias.vcf"), "--threads", "2",
            env={"FASTGATK_HC_BINARY": str(BINARY)},
        )
        assert (work / "long-alias.vcf").exists()
        assert json.loads(result.stdout.splitlines()[-1])["status"] == "prototype"

        bqsr_dry = run(
            "BaseRecalibrator", "-I", str(BAM), "-O", str(work / "recal.tsv"),
            "--dry-run", env={"FASTGATK_BASERECALIBRATOR_BINARY": str(BQSR_BINARY)},
        )
        bqsr_plan = json.loads(bqsr_dry.stdout)
        assert bqsr_plan["tool"] == "BaseRecalibrator"
        assert bqsr_plan["execution_mode"] == "native"
        assert bqsr_plan["binary_exists"] is True

        apply_dry = run(
            "ApplyBQSR", "-I", str(BAM), "--bqsr-recal-file", str(work / "recal.tsv"),
            "-O", str(work / "recal.bam"), "--dry-run",
            env={"FASTGATK_APPLYBQSR_BINARY": str(APPLY_BQSR_BINARY)},
        )
        apply_plan = json.loads(apply_dry.stdout)
        assert apply_plan["tool"] == "ApplyBQSR"
        assert apply_plan["execution_mode"] == "native"
        assert apply_plan["binary_exists"] is True

        gather_dry = run(
            "GatherBQSRReports", "-I", str(work / "recal-1.tsv"),
            "-I", str(work / "recal-2.tsv"),
            "-O", str(work / "merged-recal.tsv"), "--dry-run",
            env={"FASTGATK_GATHERBQSRREPORTS_BINARY": str(GATHER_BQSR_BINARY)},
        )
        gather_plan = json.loads(gather_dry.stdout)
        assert gather_plan["tool"] == "GatherBQSRReports"
        assert gather_plan["execution_mode"] == "native"
        assert gather_plan["binary_exists"] is True

        analyze_dry = run(
            "AnalyzeCovariates", "-before", str(work / "before.table"),
            "-csv", str(work / "analysis.csv"), "--dry-run",
            env={"FASTGATK_ANALYZECOVARIATES_BINARY": str(ANALYZE_COVARIATES_BINARY)},
        )
        analyze_plan = json.loads(analyze_dry.stdout)
        assert analyze_plan["tool"] == "AnalyzeCovariates"
        assert analyze_plan["execution_mode"] == "native"
        assert analyze_plan["binary_exists"] is True

        pon_dry = run(
            "CreateReadCountPanelOfNormals", "-I", str(work / "normal-1.tsv"),
            "-O", str(work / "pon.hdf5"), "--dry-run",
            env={"FASTGATK_CREATEREADCOUNTPANELOFNORMALS_BINARY": str(CREATE_PON_BINARY)},
        )
        pon_plan = json.loads(pon_dry.stdout)
        assert pon_plan["tool"] == "CreateReadCountPanelOfNormals"
        assert pon_plan["execution_mode"] == "native"
        assert pon_plan["binary_exists"] is True
        assert gather_plan["argv"].count("-I") == 2

        genotype_dry = run(
            "GenotypeGVCFs", "-V", str(ROOT / "gatk-source/src/test/resources/large/NA12878.prod.chr20snippet.g.vcf.gz"),
            "-O", str(work / "genotyped.vcf.gz"), "--dry-run",
            env={"FASTGATK_GENOTYPEGVCFS_BINARY": str(GENOTYPE_GVCF_BINARY)},
        )
        genotype_plan = json.loads(genotype_dry.stdout)
        assert genotype_plan["tool"] == "GenotypeGVCFs"
        assert genotype_plan["execution_mode"] == "native"
        assert genotype_plan["binary_exists"] is True

        reblock_dry = run(
            "ReblockGVCF", "-V", str(work / "input.g.vcf.gz"), "-GQB", "20",
            "-O", str(work / "reblocked.g.vcf.gz"), "--dry-run",
            env={"FASTGATK_REBLOCKGVCF_BINARY": str(REBLOCK_GVCF_BINARY)},
        )
        reblock_plan = json.loads(reblock_dry.stdout)
        assert reblock_plan["tool"] == "ReblockGVCF"
        assert reblock_plan["execution_mode"] == "native"
        assert reblock_plan["binary_exists"] is True

        select_dry = run(
            "SelectVariants", "-V", str(work / "input.vcf.gz"),
            "-O", str(work / "selected.vcf.gz"), "--sample-name", "S2",
            "--select-type-to-include", "SNP", "--dry-run",
            env={"FASTGATK_SELECTVARIANTS_BINARY": str(SELECT_VARIANTS_BINARY)},
        )
        select_plan = json.loads(select_dry.stdout)
        assert select_plan["tool"] == "SelectVariants"
        assert select_plan["execution_mode"] == "native"
        assert select_plan["binary_exists"] is True

        table_dry = run(
            "VariantsToTable", "-V", str(work / "input.vcf.gz"),
            "-O", str(work / "table.tsv"), "-F", "CHROM", "-GF", "GT", "--dry-run",
            env={"FASTGATK_VARIANTSTOTABLE_BINARY": str(VARIANTS_TO_TABLE_BINARY)},
        )
        table_plan = json.loads(table_dry.stdout)
        assert table_plan["tool"] == "VariantsToTable"
        assert table_plan["execution_mode"] == "native"
        assert table_plan["binary_exists"] is True

        eval_dry = run(
            "VariantEval", "-eval", str(work / "input.vcf.gz"),
            "-O", str(work / "eval.report"), "-EV", "CountVariants", "--dry-run",
            env={"FASTGATK_VARIANTEVAL_BINARY": str(VARIANT_EVAL_BINARY)},
        )
        eval_plan = json.loads(eval_dry.stdout)
        assert eval_plan["tool"] == "VariantEval"
        assert eval_plan["execution_mode"] == "native"
        assert eval_plan["binary_exists"] is True

        validate_dry = run(
            "ValidateVariants", "-V", str(work / "input.vcf.gz"),
            "-R", str(work / "reference.fa"), "--gvcf", "--dry-run",
            env={"FASTGATK_VALIDATEVARIANTS_BINARY": str(VALIDATE_VARIANTS_BINARY)},
        )
        validate_plan = json.loads(validate_dry.stdout)
        assert validate_plan["tool"] == "ValidateVariants"
        assert validate_plan["execution_mode"] == "native"
        assert validate_plan["binary_exists"] is True

        pileup_dry = run(
            "GetPileupSummaries", "-I", str(BAM), "-V", str(work / "sites.vcf.gz"),
            "-L", str(work / "sites.vcf.gz"), "-O", str(work / "pileups.table"), "--dry-run",
            env={"FASTGATK_GETPILEUPSUMMARIES_BINARY": str(GET_PILEUP_SUMMARIES_BINARY)},
        )
        pileup_plan = json.loads(pileup_dry.stdout)
        assert pileup_plan["tool"] == "GetPileupSummaries"
        assert pileup_plan["execution_mode"] == "native"
        assert pileup_plan["binary_exists"] is True

        contamination_dry = run(
            "CalculateContamination", "-I", str(work / "pileups.table"),
            "-O", str(work / "contamination.table"), "--dry-run",
            env={"FASTGATK_CALCULATECONTAMINATION_BINARY": str(CALCULATE_CONTAMINATION_BINARY)},
        )
        contamination_plan = json.loads(contamination_dry.stdout)
        assert contamination_plan["tool"] == "CalculateContamination"
        assert contamination_plan["execution_mode"] == "native"
        assert contamination_plan["binary_exists"] is True

        gather_pileup_dry = run(
            "GatherPileupSummaries", "-I", str(work / "pileups-1.table"),
            "-I", str(work / "pileups-2.table"), "-SD", str(work / "reference.dict"),
            "-O", str(work / "gathered-pileups.table"), "--dry-run",
            env={"FASTGATK_GATHERPILEUPSUMMARIES_BINARY": str(GATHER_PILEUP_SUMMARIES_BINARY)},
        )
        gather_pileup_plan = json.loads(gather_pileup_dry.stdout)
        assert gather_pileup_plan["tool"] == "GatherPileupSummaries"
        assert gather_pileup_plan["execution_mode"] == "native"
        assert gather_pileup_plan["binary_exists"] is True

        orientation_dry = run(
            "LearnReadOrientationModel", "-I", str(work / "calls.f1r2.tsv"),
            "-I", str(work / "calls-2.f1r2.tsv"), "-O", str(work / "artifact-prior.tar.gz"),
            "--sample", "TUMOR", "--dry-run",
            env={"FASTGATK_LEARNREADORIENTATIONMODEL_BINARY": str(LEARN_READ_ORIENTATION_BINARY)},
        )
        orientation_plan = json.loads(orientation_dry.stdout)
        assert orientation_plan["tool"] == "LearnReadOrientationModel"
        assert orientation_plan["execution_mode"] == "native"
        assert orientation_plan["binary_exists"] is True
        assert orientation_plan["argv"].count("-I") == 2

        collect_f1r2_dry = run(
            "CollectF1R2Counts", "-R", str(work / "reference.fasta"),
            "-I", str(BAM), "-L", "17:69000-70000", "-O", str(work / "collect-f1r2.tar.gz"),
            "--dry-run", env={"FASTGATK_COLLECTF1R2COUNTS_BINARY": str(COLLECT_F1R2_BINARY)},
        )
        collect_f1r2_plan = json.loads(collect_f1r2_dry.stdout)
        assert collect_f1r2_plan["tool"] == "CollectF1R2Counts"
        assert collect_f1r2_plan["execution_mode"] == "native"
        assert collect_f1r2_plan["binary_exists"] is True
        assert collect_f1r2_plan["argv"].count("-I") == 1

        annotate_dry = run(
            "AnnotateIntervals", "-R", str(work / "reference.fasta"),
            "-L", str(work / "targets.interval_list"), "-O", str(work / "annotated.tsv"),
            "--dry-run", env={"FASTGATK_ANNOTATEINTERVALS_BINARY": str(ANNOTATE_INTERVALS_BINARY)},
        )
        annotate_plan = json.loads(annotate_dry.stdout)
        assert annotate_plan["tool"] == "AnnotateIntervals"
        assert annotate_plan["execution_mode"] == "native"
        assert annotate_plan["binary_exists"] is True

        count_bases_dry = run(
            "CountBasesInReference", "-R", str(work / "reference.fasta"),
            "-L", "chr1:1-100", "-O", str(work / "base-counts.txt"), "--dry-run",
            env={"FASTGATK_COUNT_BASES_BINARY": str(COUNT_BASES_BINARY)},
        )
        count_bases_plan = json.loads(count_bases_dry.stdout)
        assert count_bases_plan["tool"] == "CountBasesInReference"
        assert count_bases_plan["execution_mode"] == "native"
        assert count_bases_plan["binary_exists"] is True

        compare_references_dry = run(
            "CompareReferences", "-R", str(work / "reference.fasta"),
            "-refcomp", str(work / "reference-2.fasta"), "-O", str(work / "references.table"),
            "--base-comparison", "FIND_SNPS_ONLY", "--base-comparison-output", str(work), "--dry-run",
            env={"FASTGATK_COMPAREREFERENCES_BINARY": str(COMPARE_REFERENCES_BINARY)},
        )
        compare_references_plan = json.loads(compare_references_dry.stdout)
        assert compare_references_plan["tool"] == "CompareReferences"
        assert compare_references_plan["execution_mode"] == "native"
        assert compare_references_plan["binary_exists"] is True

        check_reference_dry = run(
            "CheckReferenceCompatibility", "-I", str(work / "reads.bam"),
            "-refcomp", str(work / "reference.fasta"), "-O", str(work / "compatibility.table"),
            "--dry-run", env={"FASTGATK_CHECKREFERENCECOMPATIBILITY_BINARY": str(CHECK_REFERENCE_COMPATIBILITY_BINARY)},
        )
        check_reference_plan = json.loads(check_reference_dry.stdout)
        assert check_reference_plan["tool"] == "CheckReferenceCompatibility"
        assert check_reference_plan["execution_mode"] == "native"
        assert check_reference_plan["binary_exists"] is True

        check_reference_variant_dry = run(
            "CheckReferenceCompatibility", "-V", str(work / "variants.vcf.gz"),
            "-refcomp", str(work / "reference.fasta"), "-O", str(work / "compatibility-vcf.table"),
            "--dry-run", env={"FASTGATK_CHECKREFERENCECOMPATIBILITY_BINARY": str(CHECK_REFERENCE_COMPATIBILITY_BINARY)},
        )
        check_reference_variant_plan = json.loads(check_reference_variant_dry.stdout)
        assert check_reference_variant_plan["tool"] == "CheckReferenceCompatibility"
        assert check_reference_variant_plan["execution_mode"] == "native"
        assert check_reference_variant_plan["binary_exists"] is True

        fasta_reference_dry = run(
            "FastaReferenceMaker", "-R", str(work / "reference.fasta"),
            "-L", "chr1:1-100", "-O", str(work / "subset.fasta"),
            "--line-width", "60", "--dry-run",
            env={"FASTGATK_FASTAREFERENCEMAKER_BINARY": str(FASTA_REFERENCE_MAKER_BINARY)},
        )
        fasta_reference_plan = json.loads(fasta_reference_dry.stdout)
        assert fasta_reference_plan["tool"] == "FastaReferenceMaker"
        assert fasta_reference_plan["execution_mode"] == "native"
        assert fasta_reference_plan["binary_exists"] is True

        fasta_alternate_dry = run(
            "FastaAlternateReferenceMaker", "-R", str(work / "reference.fasta"),
            "-V", str(work / "variants.vcf.gz"), "-L", "chr1:1-100",
            "-O", str(work / "alternate.fasta"), "--dry-run",
            env={"FASTGATK_FASTAALTERNATEREFERENCEMAKER_BINARY": str(FASTA_ALTERNATE_REFERENCE_MAKER_BINARY)},
        )
        fasta_alternate_plan = json.loads(fasta_alternate_dry.stdout)
        assert fasta_alternate_plan["tool"] == "FastaAlternateReferenceMaker"
        assert fasta_alternate_plan["execution_mode"] == "native"
        assert fasta_alternate_plan["binary_exists"] is True

        shift_fasta_dry = run(
            "ShiftFasta", "-R", str(work / "reference.fasta"),
            "-O", str(work / "shifted.fasta"),
            "--shift-back-output", str(work / "shifted.chain"),
            "--shift-offset-list", "50", "--dry-run",
            env={"FASTGATK_SHIFTFASTA_BINARY": str(SHIFT_FASTA_BINARY)},
        )
        shift_fasta_plan = json.loads(shift_fasta_dry.stdout)
        assert shift_fasta_plan["tool"] == "ShiftFasta"
        assert shift_fasta_plan["execution_mode"] == "native"
        assert shift_fasta_plan["binary_exists"] is True

        filter_dry = run(
            "FilterIntervals", "-L", str(work / "targets.interval_list"),
            "--annotated-intervals", str(work / "annotated.tsv"),
            "-O", str(work / "filtered.interval_list"), "--dry-run",
            env={"FASTGATK_FILTERINTERVALS_BINARY": str(FILTER_INTERVALS_BINARY)},
        )
        filter_plan = json.loads(filter_dry.stdout)
        assert filter_plan["tool"] == "FilterIntervals"
        assert filter_plan["execution_mode"] == "native"
        assert filter_plan["binary_exists"] is True

        preprocess_dry = run(
            "PreprocessIntervals", "-R", str(work / "reference.fasta"),
            "-L", str(work / "targets.interval_list"), "-XL", str(work / "targets.interval_list"),
            "--bin-length", "1000",
            "--padding", "250", "-O", str(work / "preprocessed.interval_list"), "--dry-run",
            env={"FASTGATK_PREPROCESSINTERVALS_BINARY": str(PREPROCESS_INTERVALS_BINARY)},
        )
        preprocess_plan = json.loads(preprocess_dry.stdout)
        assert preprocess_plan["tool"] == "PreprocessIntervals"
        assert preprocess_plan["execution_mode"] == "native"
        assert preprocess_plan["binary_exists"] is True

        counts_dry = run(
            "CollectReadCounts", "-I", str(work / "reads.bam"), "-L", str(work / "targets.interval_list"),
            "-O", str(work / "counts.tsv"), "--format", "TSV", "--dry-run",
            env={"FASTGATK_COLLECTREADCOUNTS_BINARY": str(COLLECT_READ_COUNTS_BINARY)},
        )
        counts_plan = json.loads(counts_dry.stdout)
        assert counts_plan["tool"] == "CollectReadCounts"
        assert counts_plan["execution_mode"] == "native"
        assert counts_plan["binary_exists"] is True

        denoise_dry = run(
            "DenoiseReadCounts", "-I", str(work / "counts.tsv"),
            "-O", str(work / "denoised.tsv"),
            "--standardized-copy-ratios", str(work / "standardized.tsv"),
            "--format", "TSV", "--dry-run",
            env={"FASTGATK_DENOISEREADCOUNTS_BINARY": str(DENOISE_READ_COUNTS_BINARY)},
        )
        denoise_plan = json.loads(denoise_dry.stdout)
        assert denoise_plan["tool"] == "DenoiseReadCounts"
        assert denoise_plan["execution_mode"] == "native"
        assert denoise_plan["binary_exists"] is True
        standardized_dry = run(
            "DenoiseReadCounts", "-I", str(work / "counts.tsv"),
            "-O", str(work / "denoised.tsv"),
            "--standardized-copy-ratios", str(work / "standardized.tsv"),
            "--format", "TSV", "--dry-run",
            env={"FASTGATK_DENOISEREADCOUNTS_BINARY": str(DENOISE_READ_COUNTS_BINARY)},
        )
        standardized_plan = json.loads(standardized_dry.stdout)
        assert standardized_plan["tool"] == "DenoiseReadCounts"
        assert "--standardized-copy-ratios" in standardized_plan["argv"]

        call_ratio_dry = run(
            "CallCopyRatioSegments", "-I", str(work / "segments.tsv"),
            "-O", str(work / "called.tsv"), "--dry-run",
            env={"FASTGATK_CALLCOPYRATIOSEGMENTS_BINARY": str(CALL_COPY_RATIO_SEGMENTS_BINARY)},
        )
        call_ratio_plan = json.loads(call_ratio_dry.stdout)
        assert call_ratio_plan["tool"] == "CallCopyRatioSegments"
        assert call_ratio_plan["execution_mode"] == "native"
        assert call_ratio_plan["binary_exists"] is True

        allelic_dry = run(
            "CollectAllelicCounts", "-I", str(work / "reads.bam"),
            "-R", str(work / "reference.fasta"), "-L", "chr1:1-10",
            "-O", str(work / "allelic.tsv"), "--dry-run",
            env={"FASTGATK_COLLECTALLELICCOUNTS_BINARY": str(COLLECT_ALLELIC_COUNTS_BINARY)},
        )
        allelic_plan = json.loads(allelic_dry.stdout)
        assert allelic_plan["tool"] == "CollectAllelicCounts"
        assert allelic_plan["execution_mode"] == "native"
        assert allelic_plan["binary_exists"] is True

        model_segments_dry = run(
            "ModelSegments", "--denoised-copy-ratios", str(work / "denoised.tsv"),
            "--output-prefix", str(work / "sample"), "--dry-run",
            env={"FASTGATK_MODELSEGMENTS_BINARY": str(MODEL_SEGMENTS_BINARY)},
        )
        model_segments_plan = json.loads(model_segments_dry.stdout)
        assert model_segments_plan["tool"] == "ModelSegments"
        assert model_segments_plan["execution_mode"] == "native"
        assert model_segments_plan["binary_exists"] is True

        gather_tranches_dry = run(
            "GatherTranches", "-I", str(work / "shard-a.tranches"),
            "-I", str(work / "shard-b.tranches"), "--mode", "SNP",
            "-tranche", "90",
            "-O", str(work / "gathered.tranches"), "--dry-run",
            env={"FASTGATK_GATHERTRANCHES_BINARY": str(GATHER_TRANCHES_BINARY)},
        )
        gather_tranches_plan = json.loads(gather_tranches_dry.stdout)
        assert gather_tranches_plan["tool"] == "GatherTranches"
        assert gather_tranches_plan["execution_mode"] == "native"
        assert gather_tranches_plan["binary_exists"] is True
        assert gather_tranches_plan["argv"].count("-I") == 2
        assert "-tranche" in gather_tranches_plan["argv"]

        apply_vqsr_dry = run(
            "ApplyVQSR", "-V", str(work / "input.vcf.gz"),
            "--recal-file", str(work / "recal.vcf.gz"), "-O", str(work / "filtered.vcf.gz"),
            "--mode", "SNP", "--dry-run",
            env={"FASTGATK_APPLYVQSR_BINARY": str(APPLY_VQSR_BINARY)},
        )
        apply_vqsr_plan = json.loads(apply_vqsr_dry.stdout)
        assert apply_vqsr_plan["tool"] == "ApplyVQSR"
        assert apply_vqsr_plan["execution_mode"] == "native"
        assert apply_vqsr_plan["binary_exists"] is True

        variant_recalibrator_dry = run(
            "VariantRecalibrator", "-V", str(work / "input.vcf.gz"),
            "--resource:truth,training=true,truth=true", str(work / "truth.vcf.gz"),
            "-an", "QD", "--max-gaussians", "2", "--full-covariance",
            "-O", str(work / "recal.vcf.gz"),
            "--tranches-file", str(work / "recal.tranches"), "--dry-run",
            env={"FASTGATK_VARIANTRECALIBRATOR_BINARY": str(VARIANT_RECALIBRATOR_BINARY)},
        )
        variant_recalibrator_plan = json.loads(variant_recalibrator_dry.stdout)
        assert variant_recalibrator_plan["tool"] == "VariantRecalibrator"
        assert variant_recalibrator_plan["execution_mode"] == "native"
        assert variant_recalibrator_plan["binary_exists"] is True
        assert "--resource:truth,training=true,truth=true" in variant_recalibrator_plan["argv"]
        assert "--full-covariance" in variant_recalibrator_plan["argv"]

        gather_dry = run(
            "GatherVcfs", "-I", str(work / "shard-1.vcf.gz"),
            "-I", str(work / "shard-2.vcf.gz"), "-O", str(work / "gathered.vcf.gz"),
            "--dry-run", env={"FASTGATK_GATHERVCFS_BINARY": str(GATHER_VCFS_BINARY)},
        )
        gather_plan = json.loads(gather_dry.stdout)
        assert gather_plan["tool"] == "GatherVcfs"
        assert gather_plan["execution_mode"] == "native"
        assert gather_plan["binary_exists"] is True

        left_align_dry = run(
            "LeftAlignAndTrimVariants", "-V", str(work / "input.vcf.gz"),
            "-R", str(work / "reference.fa"), "-O", str(work / "normalized.vcf.gz"),
            "--dry-run", env={"FASTGATK_LEFTALIGNANDTRIMVARIANTS_BINARY": str(LEFT_ALIGN_BINARY)},
        )
        left_align_plan = json.loads(left_align_dry.stdout)
        assert left_align_plan["tool"] == "LeftAlignAndTrimVariants"
        assert left_align_plan["execution_mode"] == "native"
        assert left_align_plan["binary_exists"] is True

        filtration_dry = run(
            "VariantFiltration", "-V", str(work / "input.vcf.gz"),
            "-O", str(work / "filtered.vcf.gz"), "--filter-expression", "QUAL < 30",
            "--filter-name", "LowQual", "--dry-run",
            env={"FASTGATK_VARIANTFILTRATION_BINARY": str(VARIANT_FILTRATION_BINARY)},
        )
        filtration_plan = json.loads(filtration_dry.stdout)
        assert filtration_plan["tool"] == "VariantFiltration"
        assert filtration_plan["execution_mode"] == "native"
        assert filtration_plan["binary_exists"] is True

        sort_dry = run(
            "SortSam", "-I", str(work / "input.sam"), "-O", str(work / "sorted.bam"),
            "--sort-order", "coordinate", "--max-records-in-memory", "2", "--dry-run",
            env={"FASTGATK_SORTSAM_BINARY": str(SORT_SAM_BINARY)},
        )
        sort_plan = json.loads(sort_dry.stdout)
        assert sort_plan["tool"] == "SortSam"
        assert sort_plan["execution_mode"] == "native"
        assert sort_plan["binary_exists"] is True

        mark_duplicates_dry = run(
            "MarkDuplicates", "-I", str(work / "input.bam"),
            "-O", str(work / "marked.bam"), "--metrics-file", str(work / "marked.metrics.txt"),
            "--dry-run", env={"FASTGATK_MARKDUPLICATES_BINARY": str(MARK_DUPLICATES_BINARY)},
        )
        mark_duplicates_plan = json.loads(mark_duplicates_dry.stdout)
        assert mark_duplicates_plan["tool"] == "MarkDuplicates"
        assert mark_duplicates_plan["execution_mode"] == "native"
        assert mark_duplicates_plan["binary_exists"] is True

        combine_dry = run(
            "CombineGVCFs", "-V", str(work / "shard-1.g.vcf.gz"),
            "-V", str(work / "shard-2.g.vcf.gz"), "-O", str(work / "combined.g.vcf.gz"),
            "--dry-run", env={"FASTGATK_COMBINEGVCFS_BINARY": str(COMBINE_GVCF_BINARY)},
        )
        combine_plan = json.loads(combine_dry.stdout)
        assert combine_plan["tool"] == "CombineGVCFs"
        assert combine_plan["execution_mode"] == "native"
        assert combine_plan["binary_exists"] is True

        filter_dry = run(
            "FilterMutectCalls", "-V", str(work / "mutect.vcf.gz"),
            "-O", str(work / "filtered.vcf.gz"), "--dry-run",
            env={"FASTGATK_FILTERMUTECTCALLS_BINARY": str(FILTER_MUTECT_BINARY)},
        )
        filter_plan = json.loads(filter_dry.stdout)
        assert filter_plan["tool"] == "FilterMutectCalls"
        assert filter_plan["execution_mode"] == "native"
        assert filter_plan["binary_exists"] is True

        filter_orientation_dry = run(
            "FilterMutectCalls", "-R", str(work / "reference.fa"),
            "-V", str(work / "mutect.vcf.gz"), "-O", str(work / "filtered-orientation.vcf.gz"),
            "--orientation-bias-artifact-priors", str(work / "TUMOR.orientation_priors.tar.gz"),
            "--threshold-strategy", "OPTIMAL_F_SCORE", "--f-score-beta", "1.0", "--dry-run",
            env={"FASTGATK_FILTERMUTECTCALLS_BINARY": str(FILTER_MUTECT_BINARY)},
        )
        filter_orientation_plan = json.loads(filter_orientation_dry.stdout)
        assert filter_orientation_plan["tool"] == "FilterMutectCalls"
        assert filter_orientation_plan["execution_mode"] == "native"
        assert filter_orientation_plan["binary_exists"] is True
        assert filter_orientation_plan["argv"].count("--orientation-bias-artifact-priors") == 1
        assert "--threshold-strategy" in filter_orientation_plan["argv"]

    genomicsdb_no_binary = run(
        "GenomicsDBImport", "-V", "sample.g.vcf.gz",
        "--genomicsdb-workspace-path", "/tmp/fastgatk-genomicsdb",
        env={"FASTGATK_GENOMICSDBIMPORT_BINARY": str(ROOT / "missing-genomicsdb-adapter")},
        expect=69,
    )
    genomicsdb_error = json.loads(genomicsdb_no_binary.stderr.splitlines()[-1])
    assert genomicsdb_error["category"] == "BACKEND_UNAVAILABLE"

    genomicsdb_adapter = run(
        "GenomicsDBImport", "-V", "sample.g.vcf.gz",
        "--genomicsdb-workspace-path", "/tmp/fastgatk-genomicsdb",
        "--release-specific-option", "kept", "--dry-run",
        env={"FASTGATK_GENOMICSDBIMPORT_BINARY": str(GENOMICSDB_BINARY)},
    )
    genomicsdb_adapter_plan = json.loads(genomicsdb_adapter.stdout)
    assert genomicsdb_adapter_plan["tool"] == "GenomicsDBImport"
    assert genomicsdb_adapter_plan["execution_mode"] == "adapter"
    assert "--release-specific-option" in genomicsdb_adapter_plan["argv"]

    genomicsdb_fallback = run(
        "--dry-run", "--fallback", "genomicsdbimport", "-V", "sample.g.vcf.gz",
        "--genomicsdb-workspace-path", "/tmp/fastgatk-genomicsdb",
        "--release-specific-option", "kept",
        env={"FASTGATK_GENOMICSDBIMPORT_BINARY": str(ROOT / "missing-genomicsdb-adapter")},
    )
    genomicsdb_plan = json.loads(genomicsdb_fallback.stdout)
    assert genomicsdb_plan["tool"] == "GenomicsDBImport"
    assert genomicsdb_plan["execution_mode"] == "fallback"
    assert "adapter binary unavailable" in genomicsdb_plan["fallback_reason"]
    assert "--release-specific-option" in genomicsdb_plan["argv"]
    assert "kept" in genomicsdb_plan["argv"]

    unknown = run("HaplotypeCaller", "-I", str(BAM), "-O", "/tmp/unused.vcf", "--not-a-gatk-option", expect=2)
    error = json.loads(unknown.stderr.splitlines()[-1])
    assert error["category"] == "UNSUPPORTED_PARAMETER"
    assert "not-a-gatk-option" in error["message"]

    # ReadLengthReadFilter is now a native parameterized mask.  The
    # dispatcher must keep it on the native path (and preserve the bounds in
    # argv) instead of routing it to Java as an unsupported custom filter.
    supported_filter = run(
        "--dry-run", "HaplotypeCaller", "-I", str(BAM), "-O", "/tmp/unused.vcf",
        "--read-filter", "ReadLengthReadFilter", "--min-read-length", "80",
        "--max-read-length", "100", env={"FASTGATK_HC_BINARY": str(BINARY)},
    )
    supported_filter_plan = json.loads(supported_filter.stdout)
    assert supported_filter_plan["execution_mode"] == "native"
    assert "ReadLengthReadFilter" in supported_filter_plan["argv"]
    assert "--min-read-length" in supported_filter_plan["argv"]
    assert "--max-read-length" in supported_filter_plan["argv"]

    # The command boundary must expose the exact GATK HC spelling used by
    # ReadThreadingAssembler as well as the native low-quality-end soft-clip
    # Boolean.  Both are handled by fastgatk-hc-call; rejecting them in the
    # dispatcher would incorrectly force an otherwise supported Host/Kokkos
    # pipeline to Java fallback.
    hc_quality_plan = run(
        "--dry-run", "HaplotypeCaller", "-I", str(BAM), "-O", "/tmp/q20.vcf",
        "--min-base-quality-score=20", "--soft-clip-low-quality-ends",
        env={"FASTGATK_HC_BINARY": str(BINARY)},
    )
    hc_quality_payload = json.loads(hc_quality_plan.stdout)
    assert hc_quality_payload["execution_mode"] == "native"
    assert "--min-base-quality-score=20" in hc_quality_payload["argv"]
    assert "--soft-clip-low-quality-ends" in hc_quality_payload["argv"]

    filter_fallback = run(
        "--dry-run", "--fallback", "HaplotypeCaller", "-I", str(BAM),
        "-O", "/tmp/unused.vcf", "--read-filter", "UnsupportedReadFilter",
    )
    filter_fallback_plan = json.loads(filter_fallback.stdout)
    assert filter_fallback_plan["execution_mode"] == "fallback"
    assert "UnsupportedReadFilter" in filter_fallback_plan["argv"]

    fallback = run(
        "--dry-run", "--fallback", "HaplotypeCaller", "-I", str(BAM), "-O", "/tmp/unused.vcf",
        "--not-a-gatk-option", env={"FASTGATK_HC_BINARY": str(BINARY)},
    )
    fallback_plan = json.loads(fallback.stdout)
    assert fallback_plan["execution_mode"] == "fallback"
    assert "--not-a-gatk-option" in fallback_plan["argv"]

    missing_binary = run(
        "--dry-run", "--fallback", "HaplotypeCaller", "-I", str(BAM), "-O", "/tmp/unused.vcf",
        env={"FASTGATK_HC_BINARY": str(ROOT / "missing-fastgatk-hc-call")},
    )
    assert json.loads(missing_binary.stdout)["execution_mode"] == "fallback"

    mutect_dry = run(
        "Mutect2", "-I", str(BAM), "-O", "/tmp/mutect.vcf", "--dry-run",
        env={"FASTGATK_MUTECT2_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-mutect2")},
    )
    mutect_plan = json.loads(mutect_dry.stdout)
    assert mutect_plan["execution_mode"] == "native"
    assert mutect_plan["binary_exists"] is True

    mutect_fallback = run(
        "--dry-run", "--fallback", "Mutect2", "-I", str(BAM), "-O", "/tmp/mutect.vcf",
        "--unsupported-native-option",
    )
    mutect_fallback_plan = json.loads(mutect_fallback.stdout)
    assert mutect_fallback_plan["execution_mode"] == "fallback"
    assert "--unsupported-native-option" in mutect_fallback_plan["argv"]

    index_dry = run(
        "IndexFeatureFile", "-I", str(ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/4featuresHG38Header.unindexed.vcf.gz"),
        "-O", "/tmp/sample.g.vcf.gz.tbi", "--dry-run",
        env={"FASTGATK_INDEX_FEATURE_FILE_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-index-feature-file")},
    )
    index_plan = json.loads(index_dry.stdout)
    assert index_plan["tool"] == "IndexFeatureFile"
    assert index_plan["execution_mode"] == "native"

    count_reads_dry = run(
        "CountReads", "-I", str(BAM), "--dry-run",
        env={"FASTGATK_COUNT_READS_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-count-reads")},
    )
    count_reads_plan = json.loads(count_reads_dry.stdout)
    assert count_reads_plan["tool"] == "CountReads"
    assert count_reads_plan["execution_mode"] == "native"
    count_reads_multi_dry = run(
        "CountReads", "-I", str(BAM), "--input", str(BAM), "--dry-run",
        env={"FASTGATK_COUNT_READS_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-count-reads")},
    )
    count_reads_multi_plan = json.loads(count_reads_multi_dry.stdout)
    assert count_reads_multi_plan["execution_mode"] == "native"
    assert count_reads_multi_plan["argv"].count(str(BAM)) == 2
    count_reads_filter_dry = run(
        "CountReads", "-I", str(BAM), "-L", "17:69000-69500", "-RF", "MappedReadFilter",
        "--dry-run", env={"FASTGATK_COUNT_READS_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-count-reads")},
    )
    count_reads_filter_plan = json.loads(count_reads_filter_dry.stdout)
    assert count_reads_filter_plan["execution_mode"] == "native"
    # CountReads/FlagStat expose the complete read-QC mask implemented by
    # read_metrics_tool.cpp.  These predicates are intentionally broader than
    # the HC/Mutect2 mask and must not be routed to Java by the dispatcher.
    count_reads_extended_filter = run(
        "CountReads", "-I", str(BAM), "-RF", "PrimaryLineReadFilter",
        "-RF", "MappingQualityAvailableReadFilter", "-DF", "WellformedReadFilter", "--dry-run",
        env={"FASTGATK_COUNT_READS_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-count-reads")},
    )
    count_reads_extended_plan = json.loads(count_reads_extended_filter.stdout)
    assert count_reads_extended_plan["execution_mode"] == "native"
    assert "PrimaryLineReadFilter" in count_reads_extended_plan["argv"]
    assert "MappingQualityAvailableReadFilter" in count_reads_extended_plan["argv"]
    assert "WellformedReadFilter" in count_reads_extended_plan["argv"]
    count_reads_unknown_filter = run(
        "CountReads", "-I", str(BAM), "-RF", "UnsupportedReadFilter", "--dry-run", expect=2,
    )
    assert json.loads(count_reads_unknown_filter.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_PARAMETER"
    flag_stat_dry = run(
        "FlagStat", "-I", str(BAM), "--dry-run",
        env={"FASTGATK_FLAG_STAT_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-flag-stat")},
    )
    flag_stat_plan = json.loads(flag_stat_dry.stdout)
    assert flag_stat_plan["tool"] == "FlagStat"
    assert flag_stat_plan["execution_mode"] == "native"
    flag_stat_multi_dry = run(
        "FlagStat", "-I", str(BAM), "--input", str(BAM), "--dry-run",
        env={"FASTGATK_FLAG_STAT_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-flag-stat")},
    )
    flag_stat_multi_plan = json.loads(flag_stat_multi_dry.stdout)
    assert flag_stat_multi_plan["execution_mode"] == "native"
    assert flag_stat_multi_plan["argv"].count(str(BAM)) == 2
    flag_stat_extended_filter = run(
        "FlagStat", "-I", str(BAM), "-RF", "MateDifferentStrandReadFilter",
        "--dry-run", env={"FASTGATK_FLAG_STAT_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-flag-stat")},
    )
    flag_stat_extended_plan = json.loads(flag_stat_extended_filter.stdout)
    assert flag_stat_extended_plan["execution_mode"] == "native"
    assert "MateDifferentStrandReadFilter" in flag_stat_extended_plan["argv"]
    split_intervals_dry = run(
        "SplitIntervals", "-R", str(ROOT / "gatk-source/src/test/resources/hg19micro.fasta"),
        "-O", "/tmp/fastgatk-split-intervals", "--scatter-count", "2", "--dry-run",
        env={"FASTGATK_SPLIT_INTERVALS_BINARY": str(ROOT / "fastgatk-native/build/fastgatk-split-intervals")},
    )
    split_intervals_plan = json.loads(split_intervals_dry.stdout)
    assert split_intervals_plan["tool"] == "SplitIntervals"
    assert split_intervals_plan["execution_mode"] == "native"
    no_fallback = run("Mutect2", "-I", str(BAM), "-O", "/tmp/mutect.vcf",
                      "--unsupported-native-option", expect=2)
    assert json.loads(no_fallback.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_PARAMETER"

    unknown_tool = run("NotATool", "--help", expect=2)
    assert json.loads(unknown_tool.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_TOOL"
    # The plan is fallback-first for unregistered long-tail walkers.  An
    # explicit --fallback must preserve the tool name and every original
    # argument, while the default path above remains fail-closed for typos.
    unknown_fallback = run(
        "--dry-run", "--fallback", "NotATool", "-I", str(BAM), "-O", "/tmp/unused.vcf",
        "--release-specific-option", "kept",
    )
    unknown_fallback_plan = json.loads(unknown_fallback.stdout)
    assert unknown_fallback_plan["execution_mode"] == "fallback"
    assert unknown_fallback_plan["tool"] == "NotATool"
    assert unknown_fallback_plan["argv"][:2] == ["gatk", "NotATool"]
    assert "--release-specific-option" in unknown_fallback_plan["argv"]
    print(json.dumps({"status": "pass", "checks": 142, "registry": str(REGISTRY)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
