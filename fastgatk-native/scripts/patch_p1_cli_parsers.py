#!/usr/bin/env python3
"""Insert shared GATK CLI consume() into native tool parsers (P1)."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "fastgatk-native/src"
INCLUDE = '#include "fastgatk/cli/gatk_common_cli.hpp"'

FILE_TOOL = {
    "hc_call.cpp": '"HaplotypeCaller"',
    "mutect2_tool.cpp": '"Mutect2"',
    "bqsr_tool.cpp": (
        "#if defined(FASTGATK_BQSR_APPLY)\n"
        '            "ApplyBQSR"\n'
        "#else\n"
        '            "BaseRecalibrator"\n'
        "#endif"
    ),
    "gather_bqsr_tool.cpp": '"GatherBQSRReports"',
    "analyze_covariates_tool.cpp": '"AnalyzeCovariates"',
    "collect_f1r2_counts_tool.cpp": '"CollectF1R2Counts"',
    "filter_mutect_tool.cpp": '"FilterMutectCalls"',
    "learn_read_orientation_model_tool.cpp": '"LearnReadOrientationModel"',
    "annotate_intervals_tool.cpp": '"AnnotateIntervals"',
    "count_bases_in_reference_tool.cpp": '"CountBasesInReference"',
    "compare_references_tool.cpp": '"CompareReferences"',
    "check_reference_compatibility_tool.cpp": '"CheckReferenceCompatibility"',
    "fasta_reference_tool.cpp": (
        "#ifdef FASTGATK_FASTA_ALTERNATE\n"
        '            "FastaAlternateReferenceMaker"\n'
        "#else\n"
        '            "FastaReferenceMaker"\n'
        "#endif"
    ),
    "shift_fasta_tool.cpp": '"ShiftFasta"',
    "index_feature_file_tool.cpp": '"IndexFeatureFile"',
    "read_metrics_tool.cpp": "kToolName",
    "split_intervals_tool.cpp": '"SplitIntervals"',
    "filter_intervals_tool.cpp": '"FilterIntervals"',
    "preprocess_intervals_tool.cpp": '"PreprocessIntervals"',
    "collect_read_counts_tool.cpp": '"CollectReadCounts"',
    "denoise_read_counts_tool.cpp": '"DenoiseReadCounts"',
    "create_read_count_panel_of_normals_tool.cpp": '"CreateReadCountPanelOfNormals"',
    "call_copy_ratio_segments_tool.cpp": '"CallCopyRatioSegments"',
    "collect_allelic_counts_tool.cpp": '"CollectAllelicCounts"',
    "depth_of_coverage_tool.cpp": '"DepthOfCoverage"',
    "model_segments_tool.cpp": '"ModelSegments"',
    "gather_tranches_tool.cpp": '"GatherTranches"',
    "apply_vqsr_tool.cpp": '"ApplyVQSR"',
    "variant_recalibrator_tool.cpp": '"VariantRecalibrator"',
    "genotype_gvcf_tool.cpp": '"GenotypeGVCFs"',
    "reblock_gvcf_tool.cpp": '"ReblockGVCF"',
    "select_variants_tool.cpp": '"SelectVariants"',
    "variants_to_table_tool.cpp": '"VariantsToTable"',
    "variant_eval_tool.cpp": '"VariantEval"',
    "gather_pileup_summaries_tool.cpp": '"GatherPileupSummaries"',
    "calculate_contamination_tool.cpp": '"CalculateContamination"',
    "get_pileup_summaries_tool.cpp": '"GetPileupSummaries"',
    "validate_variants_tool.cpp": '"ValidateVariants"',
    "gather_vcfs_tool.cpp": '"GatherVcfs"',
    "left_align_tool.cpp": '"LeftAlignAndTrimVariants"',
    "variant_filtration_tool.cpp": '"VariantFiltration"',
    "sort_sam_tool.cpp": '"SortSam"',
    "mark_duplicates_tool.cpp": '"MarkDuplicates"',
    "combine_gvcf_tool.cpp": '"CombineGVCFs"',
}

THROW_LINE = re.compile(
    r'^(\s*)throw std::(?:invalid_argument|runtime_error)\("unknown option: " \+ (argument|arg)\);\s*$'
)
ELSE_LINE = re.compile(r"^(\s*)\} else \{\s*$")
FOR_LINE = re.compile(r"for\s*\(\s*int\s+(\w+)\s*=")


def insert_include(text: str) -> str:
    if INCLUDE in text:
        return text
    lines = text.splitlines(keepends=True)
    last_include = 0
    for i, line in enumerate(lines):
        if line.startswith("#include"):
            last_include = i
    lines.insert(last_include + 1, INCLUDE + "\n")
    return "".join(lines)


def patch_file(path: Path, tool_expr: str) -> bool:
    original = path.read_text(encoding="utf-8")
    text = insert_include(original)
    lines = text.splitlines(keepends=True)
    index_var = "index"
    inserted = 0
    i = 0
    while i < len(lines):
        for_match = FOR_LINE.search(lines[i])
        if for_match:
            index_var = for_match.group(1)
        throw = THROW_LINE.match(lines[i].rstrip("\n"))
        if throw and i >= 1 and ELSE_LINE.match(lines[i - 1].rstrip("\n")):
            indent = ELSE_LINE.match(lines[i - 1].rstrip("\n")).group(1)
            var = throw.group(2)
            replacement = [
                f"{indent}}} else if (fastgatk::cli::consume({index_var}, argc, argv,\n",
                f"{indent}            {tool_expr})) {{\n",
                f"{indent}}} else {{\n",
                f"{indent}    throw std::invalid_argument(\"unknown option: \" + {var});\n",
            ]
            lines[i - 1:i + 1] = replacement
            inserted += 1
            i += len(replacement)
            continue
        i += 1
    if inserted == 0:
        return False
    path.write_text("".join(lines), encoding="utf-8")
    return True


def main() -> int:
    changed = []
    missing = []
    for name, tool in FILE_TOOL.items():
        path = SRC / name
        if not path.is_file():
            missing.append(name)
            continue
        if patch_file(path, tool):
            changed.append(name)
        elif "gatk_common_cli.hpp" not in path.read_text(encoding="utf-8"):
            missing.append(name)
    print({"changed": changed, "missing_or_unpatched": missing})
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
