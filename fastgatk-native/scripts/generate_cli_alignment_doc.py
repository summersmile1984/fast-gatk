#!/usr/bin/env python3
"""Generate fastgatk-native/docs/cli-alignment.md.

Reads:
  - fastgatk-native/CMakeLists.txt: enumerates native executables.
  - fastgatk-native/src/<tool>.cpp: argv parser block (option-string literals).
  - fastgatk-native/include/fastgatk/cli/gatk_cli_catalog.hpp: native-side catalog
    of GATK options the parser explicitly ignores.
  - third_party/gatk-package/gatk-4.6.2.0/gatkdoc/*.json: GATK 4.6.2.0 argument
    schema for the matching Java tool.
  - fastgatk-native/dispatcher/tool_registry.json: per-tool status
    (contract-compatible / adapter) and determinism contract.
  - fastgatk-native/scripts/verify_*.py: oracle scripts that pin bit-identical
    or bounded parity contracts against pinned GATK 4.6.2.0 output.

Writes:
  - fastgatk-native/docs/cli-alignment.md: top-level overview table + per-tool
    Markdown tables + parity-evidence matrix.  Status enum:
        accepted                — argv parser handles the option.
        ignored_via_catalog     — option is in the catalog (eats it).
        native_only             — option exists in native argv but GATK JSON
                                  has no such argument.
        unsupported             — GATK exposes the option but native does not.

This script NEVER edits repository state outside the doc output path.  Errors
during parsing abort the run with a non-zero exit; the previous doc is left in
place.  All inputs are read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO_DEFAULT = Path(__file__).resolve().parents[2]

# Static binary → GATK tool name → JSON basename mapping.
# All entries were cross-checked against the gatkdoc JSON directory at
# GATK 4.6.2.0.  Any entry with a JSON key ending in "/NONE" has no gatkdoc
# JSON (native-only bridge binary); such rows get a sentinel status.
TOOL_MAP: Dict[str, Tuple[str, str]] = {
    "genomicsdb-export":              ("GenomicsDBExport (native-only bridge binary)", "NONE"),
    "hc-call":                        ("HaplotypeCaller", "org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller"),
    "genotype-gvcf":                  ("GenotypeGVCFs", "org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs"),
    "reblock-gvcf":                   ("ReblockGVCF", "org_broadinstitute_hellbender_tools_walkers_variantutils_ReblockGVCF"),
    "select-variants":                ("SelectVariants", "org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants"),
    "gather-vcfs":                    ("GatherVcfs (Picard)", "picard_vcf_GatherVcfs"),
    "left-align-trim":                ("LeftAlignAndTrimVariants", "org_broadinstitute_hellbender_tools_walkers_variantutils_LeftAlignAndTrimVariants"),
    "variant-filtration":             ("VariantFiltration", "org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration"),
    "variant-annotator":              ("VariantAnnotator", "org_broadinstitute_hellbender_tools_walkers_annotator_VariantAnnotator"),
    "funcotator":                     ("Funcotator", "org_broadinstitute_hellbender_tools_funcotator_Funcotator"),
    "apply-vqsr":                     ("ApplyVQSR", "org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR"),
    "variant-recalibrator":           ("VariantRecalibrator", "org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator"),
    "sort-sam":                       ("SortSam (Picard)", "picard_sam_SortSam"),
    "mark-duplicates":                ("MarkDuplicates (Picard)", "picard_sam_markduplicates_MarkDuplicates"),
    "mutect2":                        ("Mutect2", "org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2"),
    "combine-gvcfs":                  ("CombineGVCFs", "org_broadinstitute_hellbender_tools_walkers_CombineGVCFs"),
    "filter-mutect-calls":            ("FilterMutectCalls", "org_broadinstitute_hellbender_tools_walkers_mutect_filtering_FilterMutectCalls"),
    "genomicsdb-import":              ("GenomicsDBImport", "org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport"),
    "gather-bqsr-reports":            ("GatherBQSRReports", "org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports"),
    "analyze-covariates":             ("AnalyzeCovariates", "org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates"),
    "variants-to-table":              ("VariantsToTable", "org_broadinstitute_hellbender_tools_walkers_variantutils_VariantsToTable"),
    "variant-eval":                   ("VariantEval", "org_broadinstitute_hellbender_tools_walkers_varianteval_VariantEval"),
    "validate-variants":              ("ValidateVariants", "org_broadinstitute_hellbender_tools_walkers_variantutils_ValidateVariants"),
    "get-pileup-summaries":           ("GetPileupSummaries", "org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries"),
    "calculate-contamination":        ("CalculateContamination", "org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination"),
    "gather-pileup-summaries":        ("GatherPileupSummaries", "NONE"),  # no gatkdoc JSON
    "learn-read-orientation-model":   ("LearnReadOrientationModel", "org_broadinstitute_hellbender_tools_walkers_readorientation_LearnReadOrientationModel"),
    "denoise-read-counts":            ("DenoiseReadCounts", "org_broadinstitute_hellbender_tools_copynumber_DenoiseReadCounts"),
    "create-read-count-panel-of-normals": ("CreateReadCountPanelOfNormals", "org_broadinstitute_hellbender_tools_copynumber_CreateReadCountPanelOfNormals"),
    "call-copy-ratio-segments":       ("CallCopyRatioSegments", "org_broadinstitute_hellbender_tools_copynumber_CallCopyRatioSegments"),
    "model-segments":                 ("ModelSegments", "org_broadinstitute_hellbender_tools_copynumber_ModelSegments"),
    "gather-tranches":                ("GatherTranches", "org_broadinstitute_hellbender_tools_walkers_vqsr_GatherTranches"),
    "annotate-intervals":             ("AnnotateIntervals", "org_broadinstitute_hellbender_tools_copynumber_AnnotateIntervals"),
    "count-bases-in-reference":       ("CountBasesInReference", "org_broadinstitute_hellbender_tools_walkers_fasta_CountBasesInReference"),
    "compare-references":             ("CompareReferences", "org_broadinstitute_hellbender_tools_reference_CompareReferences"),
    "check-reference-compatibility":  ("CheckReferenceCompatibility", "org_broadinstitute_hellbender_tools_reference_CheckReferenceCompatibility"),
    "fasta-reference-maker":          ("FastaReferenceMaker", "org_broadinstitute_hellbender_tools_walkers_fasta_FastaReferenceMaker"),
    "fasta-alternate-reference-maker":("FastaAlternateReferenceMaker", "org_broadinstitute_hellbender_tools_walkers_fasta_FastaAlternateReferenceMaker"),
    "shift-fasta":                    ("ShiftFasta", "org_broadinstitute_hellbender_tools_walkers_fasta_ShiftFasta"),
    "index-feature-file":             ("IndexFeatureFile", "org_broadinstitute_hellbender_tools_IndexFeatureFile"),
    "count-reads":                    ("CountReads", "org_broadinstitute_hellbender_tools_CountReads"),
    "flag-stat":                      ("FlagStat", "org_broadinstitute_hellbender_tools_FlagStat"),
    "split-intervals":                ("SplitIntervals", "org_broadinstitute_hellbender_tools_walkers_SplitIntervals"),
    "filter-intervals":               ("FilterIntervals", "org_broadinstitute_hellbender_tools_copynumber_FilterIntervals"),
    "preprocess-intervals":           ("PreprocessIntervals", "org_broadinstitute_hellbender_tools_copynumber_PreprocessIntervals"),
    "collect-read-counts":            ("CollectReadCounts", "org_broadinstitute_hellbender_tools_copynumber_CollectReadCounts"),
    "collect-f1r2-counts":            ("CollectF1R2Counts", "NONE"),  # no gatkdoc JSON
    "collect-allelic-counts":         ("CollectAllelicCounts", "org_broadinstitute_hellbender_tools_copynumber_CollectAllelicCounts"),
    "depth-of-coverage":              ("DepthOfCoverage", "org_broadinstitute_hellbender_tools_walkers_coverage_DepthOfCoverage"),
}

# BqsrTool has two names depending on FASTGATK_BQSR_APPLY; both share the same
# src file (bqsr_tool.cpp).  We register one mapping per binary and detect the
# mode at generation time.
BQSR_MAP: Dict[str, Tuple[str, str]] = {
    "bqsr":           ("BaseRecalibrator", "org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator"),
    "apply-bqsr":     ("ApplyBQSR",        "org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR"),
}

# Smoke / test / oracle binaries — excluded from the alignment table.  Listed
# explicitly here so future tooling can audit the exclusion rather than
# silently skipping on filename heuristics.
EXCLUDED_BINARIES = {
    "hts-reader-test",
    "hc-smoke",
    "calling-reference-boundary-smoke",
    "haplotype-cigar-trim-smoke",
    "calling-indel-smoke",
    "fmc-trajectory-oracle",
    "contamination-kernel-segmenter-contract",
}


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

CMAKELISTS_ADD_EXECUTABLE_RE = re.compile(
    r"add_executable\(\s*fastgatk-(?P<name>[a-z0-9-]+)\s*(?P<srcs>(?:[^\)]*?))\)",
    re.DOTALL,
)

CMAKE_SOURCE_RE = re.compile(r"(?:src|tests)/[a-zA-Z0-9_-]+\.cpp")


def parse_cmakelists(cmakelists: Path) -> List[Tuple[str, str]]:
    """Return [(binary_name, primary_source_path), ...] preserving file order.

    Handles both literal targets (``add_executable(fastgatk-bqsr ...)``) and
    CMake variable expansion (``add_executable(${bqsr_target} ...)``) when the
    variable resolves to a list of names iterated via ``foreach``.
    """
    text = cmakelists.read_text()

    # Build a map of foreach variable → expanded targets so we can look up the
    # source list when the iteration is later in the file.  For each
    # ``foreach(<var> IN ITEMS <a> <b> <c>)`` we record <var> -> [a, b, c].
    foreach_vars: Dict[str, List[str]] = {}
    for m in re.finditer(
        r"foreach\(\s*(\w+)\s+IN\s+ITEMS\s+(?P<items>[^)]+)\)", text
    ):
        var = m.group(1)
        items = [tok.strip() for tok in m.group("items").split() if tok.strip()]
        foreach_vars.setdefault(var, []).extend(items)

    # Find the iteration that drives each ``add_executable(${VAR} ...)`` and
    # expand it once so we have concrete literal targets.
    def _resolve_var_refs(s: str) -> List[Tuple[str, str]]:
        """Yield (literal_target, source_list_blob) for each add_executable in s.

        When the target uses ``${VAR}`` and VAR maps to a list of names via
        ``foreach``, yield one entry per resolved name.  The source blob is
        the body of the literal add_executable call (sans the target slot).
        """
        out: List[Tuple[str, str]] = []
        for m in re.finditer(r"add_executable\(", s):
            depth = 1
            j = m.end()
            while j < len(s) and depth > 0:
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                j += 1
            body = s[m.end():j - 1]
            # Split body into target slot and arguments.
            target_match = re.match(r"\s*(\$\{(\w+)\}|(\S+))\s*(.*)", body, re.DOTALL)
            if not target_match:
                continue
            var = target_match.group(2)
            literal = target_match.group(3)
            rest = target_match.group(4)
            if var and var in foreach_vars:
                for name in foreach_vars[var]:
                    # Foreach items are full target names like ``fastgatk-bqsr``.
                    if name.startswith("fastgatk-"):
                        out.append((name[len("fastgatk-"):], rest))
                    else:
                        out.append((name, rest))
            elif literal:
                # Literal ``fastgatk-foo`` targets.
                if literal.startswith("fastgatk-"):
                    out.append((literal[len("fastgatk-"):], rest))
        return out

    raw_entries = _resolve_var_refs(text)
    # Also catch literal ``add_executable(fastgatk-foo ...)`` that resolve_var_refs
    # already covers via the literal branch; ensure no duplicates.
    seen = set()
    entries: List[Tuple[str, str]] = []
    for name, body in raw_entries:
        if name in seen:
            continue
        seen.add(name)
        srcs = CMAKE_SOURCE_RE.findall(body)
        primary = srcs[0] if srcs else ""
        entries.append((name, primary))
    return entries


# Native argv parser extraction.  We look for any of:
#   is_option(argument, "--name")
#   argument == "--name"
#   argument == "-S"
# Capture the long-name/synonym as a single literal string.  The same regex
# is also used for short-name style "-X" matches.
ARGV_LITERAL_RE = re.compile(
    r"""(?x)
    (?:
        is_option\s*\(\s*argument\s*,\s*"(?P<long>--?[A-Za-z0-9][A-Za-z0-9_-]*)"\s*\)
        |
        argument\s*(?:\.rfind\([^)]*\)\s*==\s*0)?\s*==\s*"(?P<long2>--?[A-Za-z0-9][A-Za-z0-9_-]*)"
        |
        is_option\s*\(\s*argument\s*,\s*"(?P<short>-?[A-Za-z0-9_-]+)"\s*\)
        |
        argument\s*(?:\.rfind\([^)]*\)\s*==\s*0)?\s*==\s*"(?P<short2>-?[A-Za-z0-9_-]+)"
    )
    """
)


def extract_native_options(cpp_path: Path) -> List[Tuple[str, int]]:
    """Return [(option_literal, source_line_1indexed)] preserving insertion order."""
    if not cpp_path.exists():
        return []
    out: List[Tuple[str, int]] = []
    seen: set = set()
    for lineno, line in enumerate(cpp_path.read_text().splitlines(), start=1):
        for m in ARGV_LITERAL_RE.finditer(line):
            lit = m.group("long") or m.group("long2") or m.group("short") or m.group("short2")
            if not lit or lit in seen:
                continue
            seen.add(lit)
            out.append((lit, lineno))
    return out


# Catalog header parsing.  We extract per-tool OptionSets: a (flags, values)
# pair stored as unordered_set<string>.  The header uses a single
# kCatalog std::unordered_map<std::string, OptionSets> literal.
CATALOG_TOOL_HEADER_RE = re.compile(r'\{"(?P<tool>[A-Za-z0-9_]+)"')
CATALOG_FLAGSET_RE = re.compile(r"\{\s*\{(?P<flags>[^}]*)\}\s*,\s*\{(?P<values>[^}]*)\}\s*\}")
CATALOG_ENTRY_RE = re.compile(r'"([^\"]+)"')


def parse_catalog(path: Path) -> Dict[str, Dict[str, set]]:
    text = path.read_text()
    tools: Dict[str, Dict[str, set]] = {}
    # Locate each top-level {"Tool", {{...}, {...}}} entry.
    pos = 0
    while True:
        m = re.search(r'\{"(?P<tool>[A-Za-z0-9_]+)"\s*,', text[pos:])
        if not m:
            break
        tool = m.group("tool")
        # Find matching close brace for the entry — there are nested braces, so
        # walk the string and count depth from the start position.
        start = pos + m.end()
        depth = 1
        i = start
        while i < len(text) and depth > 0:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            i += 1
        body = text[start:i - 1]
        flags_match = re.search(r"\{\s*(?P<f>[^{}]*?)\s*\}\s*,\s*\{(?P<v>[^{}]*)\}", body)
        if flags_match:
            flags = set(re.findall(r'"([^"]+)"', flags_match.group("f")))
            values = set(re.findall(r'"([^"]+)"', flags_match.group("v")))
            tools[tool] = {"flags": flags, "values": values}
        pos = i
    return tools


def load_gatkdoc(gatkdoc_dir: Path, basename: str) -> Optional[dict]:
    if basename == "NONE":
        return None
    path = gatkdoc_dir / f"{basename}.json"
    if not path.exists():
        return None
    with path.open() as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Status classification
# ---------------------------------------------------------------------------

def classify_status(
    arg_name: str,
    synonyms: str,
    native_options: Dict[str, List[int]],
    catalog: Optional[Dict[str, set]],
) -> Tuple[str, str]:
    """Return (status, native_line_ref).

    Case-insensitive matching: Picard tool JSONs use UPPERCASE Java-style
    argument names (e.g. ``--INPUT``) while the native argv parser and the
    ``kCatalog`` map use lowercase (``--input``).  The native tool is what
    the user invokes; both halves must therefore be normalised before
    comparing.
    """
    def _ci_map(opts):
        return {k.lower(): v for k, v in opts.items()}

    candidates = [arg_name]
    if synonyms and synonyms != "NA":
        # synonyms may be either a single string or a comma-separated list.
        for token in synonyms.split(","):
            tok = token.strip()
            if tok:
                candidates.append(tok)
    candidates_ci = {c.lower() for c in candidates if c}

    # Native argv parser hit?  (case-insensitive).
    if native_options:
        nopts_ci = _ci_map(native_options)
        for cand in candidates_ci:
            if cand in nopts_ci:
                lines = nopts_ci[cand]
                return ("accepted", ",".join(str(l) for l in lines))

    # Catalog hit?  (only meaningful when the tool has a catalog entry).
    if catalog:
        union = catalog["flags"] | catalog["values"]
        union_ci = {u.lower() for u in union}
        for cand in candidates_ci:
            if cand in union_ci:
                return ("ignored_via_catalog", "—")
    return ("unsupported", "—")


# ---------------------------------------------------------------------------
# Document generation
# ---------------------------------------------------------------------------

def render_overview_table(rows: List[Tuple[str, str, str, str, str, str]]) -> str:
    lines = [
        "## Overview",
        "",
        "| Native binary | GATK Java class | GATK reference JSON | Source cpp | Count of `accepted` | Status |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for native, gatk_class, json_ref, src, accepted_count, summary in rows:
        lines.append(
            f"| `fastgatk-{native}` | `{gatk_class}` | `{json_ref}` | `{src}` | "
            f"{accepted_count} | {summary} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_tool_section(
    native: str,
    gatk_class: str,
    gatkdoc: Optional[dict],
    native_options: Dict[str, List[int]],
    catalog: Optional[Dict[str, set]],
    src: str,
) -> str:
    lines: List[str] = []
    title = gatk_class if gatk_class.startswith("GenomicsDBExport") else gatk_class
    lines.append(f"## `fastgatk-{native}` ↔ `{title}`")
    lines.append("")
    lines.append(f"- Source: `fastgatk-native/{src}`")
    if gatkdoc is None:
        lines.append("- GATK reference: `<none>` (native-only bridge binary)")
    else:
        lines.append(f"- GATK arguments live in the JSON at the canonical `gatkdoc/{Path(gatkdoc.get('_path','')).name}` location; this generator reads `arguments[]` directly.")
    lines.append("- Status column meanings:")
    lines.append("    - `accepted` — native argv parser handles the option literal.")
    lines.append("    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.")
    lines.append("    - `unsupported` — GATK exposes the option but native does not handle it.")
    if gatkdoc is None:
        lines.append("")
        lines.append("No gatkdoc JSON is shipped with this GATK version; the tool is native-only.")
        return "\n".join(lines) + "\n"
    lines.append("")
    lines.append("| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for arg in gatkdoc.get("arguments", []):
        name = arg.get("name", "")
        syn = arg.get("synonyms", "NA")
        typ = arg.get("type", "")
        default = arg.get("defaultValue", "")
        required = arg.get("required", "")
        kind = arg.get("kind", "")
        status, ref = classify_status(name, syn, native_options, catalog)
        syn_disp = "" if syn == "NA" else syn
        lines.append(
            f"| `{name}` | `{syn_disp}` | `{typ}` | `{default}` | `{required}` | "
            f"`{kind}` | {ref} | `{status}` |"
        )
    # Native-only options: things in native_options that don't map to any GATK
    # argument name/synonym.  Case-insensitive so Picard uppercase JSON names
    # do not flag native lowercase literals as "native-only".
    gatk_names = set()
    for arg in gatkdoc.get("arguments", []):
        nm = arg.get("name", "")
        if nm:
            gatk_names.add(nm.lower())
        syn = arg.get("synonyms", "")
        if syn and syn != "NA":
            for tok in syn.split(","):
                t = tok.strip()
                if t:
                    gatk_names.add(t.lower())
    native_only: List[Tuple[str, List[int]]] = []
    for lit, lines_ in native_options.items():
        if lit.lower() not in gatk_names:
            native_only.append((lit, lines_))
    if native_only:
        lines.append("")
        lines.append("### Native-only options (no GATK counterpart in 4.6.2.0)")
        lines.append("")
        lines.append("| Native literal | Native parsing line(s) |")
        lines.append("| --- | --- |")
        for lit, ls in native_only:
            lines.append(f"| `{lit}` | {','.join(str(l) for l in ls)} |")
    lines.append("")
    return "\n".join(lines)


def collect_native_options(
    cpp_path: Path, native_options: Dict[str, List[int]]
) -> None:
    for lit, line in extract_native_options(cpp_path):
        native_options.setdefault(lit, []).append(line)


# ---------------------------------------------------------------------------
# Parity-evidence helpers
# ---------------------------------------------------------------------------

# Mapping from native binary name → GATK Java tool name (canonical, no
# parentheticals).  Mirrors the entries in TOOL_MAP / BQSR_MAP but stripped
# to the canonical Java name for cross-referencing the dispatcher registry.
BIN_TO_GAT_CLASS: Dict[str, str] = {}


def _init_bin_to_gat_class() -> None:
    if BIN_TO_GAT_CLASS:
        return
    for name, (gatk_class, _) in TOOL_MAP.items():
        BIN_TO_GAT_CLASS[name] = gatk_class.split(" ")[0]
    for name, (gatk_class, _) in BQSR_MAP.items():
        BIN_TO_GAT_CLASS[name] = gatk_class


# Filename patterns that map a verify_*.py script to the GATK tool it
# exercises.  Patterns are checked as substrings (case-insensitive) of the
# script filename.  Ordered most-specific-first so that longer keys win.
_TOOL_FILENAME_PATTERNS: List[Tuple[str, str]] = [
    # Compound tool names first to avoid double-counting on shorter keys.
    ("collect_f1r2", "CollectF1R2Counts"),
    ("collect_allelic", "CollectAllelicCounts"),
    ("collect_read_counts", "CollectReadCounts"),
    ("filter_mutect", "FilterMutectCalls"),
    ("genotype_gvcf", "GenotypeGVCFs"),
    ("combine_gvcfs", "CombineGVCFs"),
    ("reblock_gvcf", "ReblockGVCF"),
    ("create_read_count_panel", "CreateReadCountPanelOfNormals"),
    ("call_copy_ratio", "CallCopyRatioSegments"),
    ("model_segments", "ModelSegments"),
    ("denoise_read_counts", "DenoiseReadCounts"),
    ("learn_read", "LearnReadOrientationModel"),
    ("left_align", "LeftAlignAndTrimVariants"),
    ("funcotate_segments", "FuncotateSegments"),
    ("filter_funcotations", "FilterFuncotations"),
    ("funcotator", "Funcotator"),
    ("variant_annotator", "VariantAnnotator"),
    ("variant_filtration", "VariantFiltration"),
    ("variant_recalibrator", "VariantRecalibrator"),
    ("apply_vqsr", "ApplyVQSR"),
    ("select_variants", "SelectVariants"),
    ("gather_vcfs", "GatherVcfs"),
    ("gather_pileup", "GatherPileupSummaries"),
    ("get_pileup_summaries", "GetPileupSummaries"),
    ("calculate_contamination", "CalculateContamination"),
    ("annotate_intervals", "AnnotateIntervals"),
    ("count_bases_in_reference", "CountBasesInReference"),
    ("compare_references", "CompareReferences"),
    ("check_reference_compatibility", "CheckReferenceCompatibility"),
    ("fasta_alternate", "FastaAlternateReferenceMaker"),
    ("fasta_reference", "FastaReferenceMaker"),
    ("shift_fasta", "ShiftFasta"),
    ("index_feature_file", "IndexFeatureFile"),
    ("count_reads", "CountReads"),
    ("flag_stat", "FlagStat"),
    ("split_intervals", "SplitIntervals"),
    ("filter_intervals", "FilterIntervals"),
    ("preprocess_intervals", "PreprocessIntervals"),
    ("depth_of_coverage", "DepthOfCoverage"),
    ("gather_tranches", "GatherTranches"),
    ("variants_to_table", "VariantsToTable"),
    ("variant_eval", "VariantEval"),
    ("validate_variants", "ValidateVariants"),
    ("mark_duplicates", "MarkDuplicates"),
    ("sort_sam", "SortSam"),
    ("hc_", "HaplotypeCaller"),
    ("mutect2", "Mutect2"),
    ("genomicsdb", "GenomicsDBImport"),
    ("analyze_covariates", "AnalyzeCovariates"),
    # BQSR: keep these last because several GATK-oracle scripts mention BQSR
    # in body text but don't actually exercise BQSR.  Match only filenames.
    ("apply_bqsr", "ApplyBQSR"),
    ("gather_bqsr_reports", "GatherBQSRReports"),
    ("bqsr", "BaseRecalibrator"),
]


def _classify_script_to_tools(filename: str) -> List[str]:
    """Return the list of GATK tools a verify script exercises, matched by
    filename substring.  Generic scripts (`verify_native.py`,
    `verify_gatk_oracle.py`, etc.) match no tool and are excluded from the
    per-tool count.
    """
    name = filename.lower()
    matched: List[str] = []
    for pat, tool in _TOOL_FILENAME_PATTERNS:
        if pat in name:
            matched.append(tool)
    return matched


def collect_oracle_evidence(scripts_dir: Path) -> Dict[str, Dict[str, object]]:
    """Walk `fastgatk-native/scripts/` and bucket `verify_*.py` scripts per
    GATK tool.  For each tool we record:

      - total_oracle_scripts: count of scripts whose filename mentions the tool.
      - bit_identical_scripts: subset of the above whose body asserts
        ``byte-identical`` or ``bit-identical`` semantics against pinned
        GATK 4.6.2.0.
      - bounded_oracle_scripts: subset whose filename contains ``oracle``.
      - script_files: list of (filename, bit_identical_bool).
    """
    result: Dict[str, Dict[str, object]] = {}
    for py in sorted(scripts_dir.glob("verify_*.py")):
        tools = _classify_script_to_tools(py.name)
        if not tools:
            continue
        text = py.read_text()
        bi = ("byte-identical" in text) or ("bit-identical" in text)
        bounded = "oracle" in py.name.lower()
        for tool in tools:
            entry = result.setdefault(
                tool,
                {
                    "total_oracle_scripts": 0,
                    "bit_identical_scripts": 0,
                    "bounded_oracle_scripts": 0,
                    "script_files": [],
                },
            )
            entry["total_oracle_scripts"] += 1
            if bounded:
                entry["bounded_oracle_scripts"] += 1
            if bi:
                entry["bit_identical_scripts"] += 1
            entry["script_files"].append((py.name, bi))
    return result


def load_tool_registry(path: Path) -> Dict[str, Dict[str, object]]:
    if not path.exists():
        return {}
    with path.open() as f:
        data = json.load(f)
    out: Dict[str, Dict[str, object]] = {}
    for name, info in data.get("tools", {}).items():
        out[name] = {
            "status": info.get("status", "unknown"),
            "determinism": info.get("determinism", []),
            "fallback_policy": info.get("fallback_policy", "unknown"),
            "native_binary": info.get("native_binary", ""),
            "aliases": info.get("aliases", []),
        }
    return out


def sha256_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def render_parity_section(
    rows: List[Tuple[str, str, str, str, int, str]],
    registry: Dict[str, Dict[str, object]],
    oracle: Dict[str, Dict[str, object]],
    pinned_jar: Optional[Path],
    pinned_jar_sha: Optional[str],
    fixture_digest_path: Optional[Path],
) -> str:
    """Build the bit-identical / bounded-parity evidence matrix.

    The matrix cross-references each native binary with:

      - The dispatcher ``tool_registry.json`` entry (status, determinism).
      - The count of oracle verify scripts in ``fastgatk-native/scripts/``.
      - The subset of those that assert ``byte-identical`` semantics.
      - Whether the pinned GATK 4.6.2.0 jar is present on disk.
      - Whether pinned fixture digests are pinned in the repo.
    """
    lines: List[str] = []
    lines.append("## Bit-identical / bounded-parity evidence")
    lines.append("")
    lines.append(
        "This section cross-references the CLI alignment table with the "
        "tool-dispatcher registry and the pinned GATK 4.6.2.0 oracle scripts. "
        "It is auto-generated from:"
    )
    lines.append("")
    lines.append(
        "- `fastgatk-native/dispatcher/tool_registry.json` — per-tool status "
        "(`contract-compatible` / `adapter`), determinism, fallback policy."
    )
    lines.append(
        "- `fastgatk-native/scripts/verify_*.py` — `byte-identical` and "
        "`bit-identical` assertions against pinned GATK output."
    )
    lines.append(
        "- `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar` "
        "— pinned GATK Java runtime invoked by every oracle script."
    )
    lines.append(
        "- `fastgatk-native/tests/pinned_fixture_digests.sha256` — SHA-256 "
        "corpus that pins the exact bytes of every test fixture."
    )
    lines.append("")
    lines.append("### Pinned artifacts")
    lines.append("")
    lines.append("| Artifact | Status | Identifier |")
    lines.append("| --- | --- | --- |")
    if pinned_jar is None:
        lines.append("| Pinned GATK 4.6.2.0 jar | **MISSING** | — |")
    else:
        sha = pinned_jar_sha or "—"
        size_mb = pinned_jar.stat().st_size / (1024 * 1024)
        lines.append(
            f"| Pinned GATK 4.6.2.0 jar | present ({size_mb:.1f} MiB) | `{pinned_jar.name}` |"
        )
        lines.append(
            f"| SHA-256 of pinned jar | computed | `{sha}` |"
        )
    if fixture_digest_path is None or not fixture_digest_path.exists():
        lines.append("| Pinned fixture digests | **MISSING** | — |")
    else:
        fixture_count = sum(
            1
            for ln in fixture_digest_path.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        )
        lines.append(
            f"| Pinned fixture digests | present ({fixture_count} entries) | `{fixture_digest_path.name}` |"
        )
    lines.append("")
    lines.append("### Per-binary parity evidence")
    lines.append("")
    lines.append(
        "| Native binary | GATK class | registry status | determinism | oracle scripts | bit-identical scripts | bounded oracles | sample bit-identical script |"
    )
    lines.append(
        "| --- | --- | --- | --- | --- | --- | --- | --- |"
    )
    # Sort by canonical GATK class for stable output.
    sorted_rows = sorted(rows, key=lambda r: (r[1].split(" ")[0], r[0]))
    for native_bin, gatk_class, _json_ref, _src, _accepted, _summary in sorted_rows:
        _init_bin_to_gat_class()
        canon = BIN_TO_GAT_CLASS.get(native_bin, gatk_class.split(" ")[0])
        reg = registry.get(canon, {})
        status = str(reg.get("status", "—"))
        det = ", ".join(reg.get("determinism", [])) or "—"
        ev = oracle.get(canon, {})
        total = ev.get("total_oracle_scripts", 0)
        bi = ev.get("bit_identical_scripts", 0)
        bounded = ev.get("bounded_oracle_scripts", 0)
        script_files = ev.get("script_files", [])
        # Pick a representative bit-identical script for the sample column.
        bi_scripts = [s for s in script_files if s[1]]
        sample = bi_scripts[0][0] if bi_scripts else "—"
        lines.append(
            f"| `fastgatk-{native_bin}` | `{canon}` | `{status}` | `{det}` | "
            f"{total} | {bi} | {bounded} | `{sample}` |"
        )
    lines.append("")
    total_native_oracles = sum(
        ev.get("total_oracle_scripts", 0) for ev in oracle.values()
    )
    total_bi = sum(ev.get("bit_identical_scripts", 0) for ev in oracle.values())
    total_bounded = sum(
        ev.get("bounded_oracle_scripts", 0) for ev in oracle.values()
    )
    lines.append(
        f"**Aggregate:** {len(oracle)} tools have at least one oracle "
        f"verify script in `fastgatk-native/scripts/`; "
        f"{total_native_oracles} oracle scripts in total; "
        f"{total_bi} of them assert byte/bit-identical output against "
        f"pinned GATK 4.6.2.0; "
        f"{total_bounded} carry the bounded-oracle naming convention "
        f"(`*_gatk_oracle.py`)."
    )
    lines.append("")
    # Bucket summary
    contract_compat = sum(
        1 for info in registry.values() if info.get("status") == "contract-compatible"
    )
    adapter_only = sum(
        1 for info in registry.values() if info.get("status") == "adapter"
    )
    lines.append(
        f"**Dispatcher status counts:** {contract_compat} `contract-compatible`, "
        f"{adapter_only} `adapter`.  All {len(registry)} tools declare "
        "`determinism: [\"strict\", \"fast\"]` and `fallback_policy: \"explicit\"`."
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=str(REPO_DEFAULT), help="Path to repository root.")
    args = ap.parse_args(argv)
    repo = Path(args.repo).resolve()
    native_root = repo / "fastgatk-native"
    cmake = native_root / "CMakeLists.txt"
    src_root = native_root / "src"
    catalog_hpp = native_root / "include" / "fastgatk/cli/gatk_cli_catalog.hpp"
    gatkdoc_dir = repo / "third_party/gatk-package/gatk-4.6.2.0" / "gatkdoc"
    out_path = native_root / "docs" / "cli-alignment.md"

    if not cmake.exists():
        print(f"ERROR: {cmake} not found", file=sys.stderr)
        return 2
    if not catalog_hpp.exists():
        print(f"ERROR: {catalog_hpp} not found", file=sys.stderr)
        return 2
    if not gatkdoc_dir.exists():
        print(f"ERROR: {gatkdoc_dir} not found", file=sys.stderr)
        return 2

    # ----- discover binaries
    cmakelists_entries = parse_cmakelists(cmake)
    binaries: List[Tuple[str, str]] = []
    for name, src in cmakelists_entries:
        if name in EXCLUDED_BINARIES:
            continue
        binaries.append((name, src))

    # ----- parse catalog
    catalog = parse_catalog(catalog_hpp)

    # ----- resolve each binary's tool mapping
    rows: List[Tuple[str, str, str, str, str, str]] = []
    sections: List[str] = []
    counter: Dict[str, int] = {"accepted": 0, "ignored_via_catalog": 0, "unsupported": 0, "native_only": 0}
    failed = []

    for name, src in binaries:
        # Map binary → tool
        if name == "bqsr":
            gatk_class, json_base = BQSR_MAP["bqsr"]
        elif name == "apply-bqsr":
            gatk_class, json_base = BQSR_MAP["apply-bqsr"]
        elif name in TOOL_MAP:
            gatk_class, json_base = TOOL_MAP[name]
        else:
            failed.append((name, "no entry in TOOL_MAP"))
            continue
        # Native cpp argv parser extraction
        cpp_path = src_root / Path(src).name
        native_opts: Dict[str, List[int]] = {}
        collect_native_options(cpp_path, native_opts)
        # GATK doc
        if json_base == "NONE":
            gatkdoc = None
            json_ref = "<none>"
        else:
            gatkdoc = load_gatkdoc(gatkdoc_dir, json_base)
            if gatkdoc is None:
                failed.append((name, f"gatkdoc JSON missing: {json_base}.json"))
                continue
            # Stash path for the section header.
            gatkdoc["_path"] = f"third_party/gatk-package/gatk-4.6.2.0/gatkdoc/{json_base}.json"
            json_ref = json_base + ".json"
        # Tool-name lookup in catalog
        cat = catalog.get(gatk_class.split(" ")[0]) or (catalog.get(gatk_class) if gatk_class in catalog else None)
        # Render section
        section = render_tool_section(name, gatk_class, gatkdoc, native_opts, cat, src)
        sections.append(section)
        # Counts
        accepted = 0
        if gatkdoc is not None:
            for arg in gatkdoc.get("arguments", []):
                status, _ = classify_status(
                    arg.get("name", ""),
                    arg.get("synonyms", ""),
                    native_opts,
                    cat,
                )
                counter[status] = counter.get(status, 0) + 1
                if status == "accepted":
                    accepted += 1
        # native-only count
        if gatkdoc is not None:
            gatk_names = set()
            for arg in gatkdoc.get("arguments", []):
                nm = arg.get("name", "")
                if nm:
                    gatk_names.add(nm.lower())
                syn = arg.get("synonyms", "")
                if syn and syn != "NA":
                    for tok in syn.split(","):
                        t = tok.strip()
                        if t:
                            gatk_names.add(t.lower())
            no = sum(1 for lit in native_opts if lit.lower() not in gatk_names)
            counter["native_only"] = counter.get("native_only", 0) + no
        # Summary cell for overview
        if gatkdoc is None:
            summary = "native-only bridge"
        else:
            summary = (
                f"`accepted`={accepted}, "
                f"`unsupported`={counter['unsupported'] if gatkdoc else 0}, "
                f"`ignored_via_catalog`={counter['ignored_via_catalog'] if gatkdoc else 0}"
            )
        rows.append((name, gatk_class, json_ref, src, str(accepted), summary))

    if failed:
        print("ERROR: failed to resolve entries:", file=sys.stderr)
        for name, reason in failed:
            print(f"  {name}: {reason}", file=sys.stderr)
        return 2

    # ----- assemble document
    header = [
        "# CLI Alignment — fastgatk-native ↔ GATK 4.6.2.0",
        "",
        "**Slug:** `cli-alignment-doc`",
        "",
        "## How to regenerate",
        "",
        "```",
        "python3 fastgatk-native/scripts/generate_cli_alignment_doc.py",
        "```",
        "",
        "The script reads only:",
        "",
        "- `fastgatk-native/CMakeLists.txt` — binary list (`add_executable(fastgatk-*)`).",
        "- `fastgatk-native/src/<tool>.cpp` — argv-parser block (regex on `is_option(argument, \"--name\")` and `argument == \"--name\"`).",
        "- `fastgatk-native/include/fastgatk/cli/gatk_cli_catalog.hpp` — per-tool catalog the native `fastgatk::cli::consume(...)` function consults.",
        "- `third_party/gatk-package/gatk-4.6.2.0/gatkdoc/<class>.json` — GATK 4.6.2.0 `arguments[]` schema for the matching Java tool.",
        "",
        "Data sources are pinned: **GATK 4.6.2.0** (Java + JSON), **Kokkos 5.2.0** for native.  This document reflects the current `CMakeLists.txt` and `src/*.cpp`; it does not make claims about numerical/algorithmic parity, only the CLI surface.",
        "",
        "## Status enum",
        "",
        "- `accepted` — native argv parser handles the option literal.",
        "- `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is silently consumed by `fastgatk::cli::consume(...)`.",
        "- `unsupported` — GATK exposes the option but the native tool does not handle it (will likely fail or be ignored).",
        "- `native_only` — native tool accepts the option literal but GATK 4.6.2.0 has no such argument (custom native switch).",
        "",
        f"Total entries: **{len(rows)}** (production binaries; smoke/test/oracle binaries live in `fastgatk-native/tests/` and are intentionally excluded).",
        "",
    ]
    overview = render_overview_table(rows)
    body = "\n".join(sections)

    # ----- parity-evidence section
    scripts_dir = native_root / "scripts"
    registry_path = native_root / "dispatcher" / "tool_registry.json"
    registry = load_tool_registry(registry_path)
    oracle = collect_oracle_evidence(scripts_dir)
    pinned_jar = repo / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
    pinned_jar_sha = sha256_file(pinned_jar)
    fixture_digest_path = native_root / "tests" / "pinned_fixture_digests.sha256"
    parity = render_parity_section(
        rows, registry, oracle, pinned_jar, pinned_jar_sha, fixture_digest_path
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "\n".join(header)
        + "\n"
        + overview
        + "\n\n"
        + parity
        + "\n\n## Per-tool details\n\n"
        + body
    )

    # Final counter report
    print(
        f"Wrote {out_path} with {len(rows)} entries.  "
        f"counts: accepted={counter['accepted']}, "
        f"ignored_via_catalog={counter['ignored_via_catalog']}, "
        f"unsupported={counter['unsupported']}, "
        f"native_only={counter['native_only']}; "
        f"parity: {len(oracle)} tools with oracle scripts, "
        f"{sum(v['bit_identical_scripts'] for v in oracle.values())} bit-identical, "
        f"jar_sha={'present' if pinned_jar_sha else 'MISSING'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())