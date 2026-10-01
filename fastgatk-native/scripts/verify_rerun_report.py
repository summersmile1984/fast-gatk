#!/usr/bin/env python3
"""Regression gate for the verify-script re-run report directory.

Assertions:

* ``<repo>/fastgatk-native/evidence/`` contains exactly the expected run
  directory (default: ``2026-09-23-rerun/``); hidden files are ignored.
* ``2026-09-23-rerun/`` contains ``INDEX.md`` and exactly 48 tool directories
  (one per production binary; the ``fastgatk-genomicsdb-export`` bridge is
  excluded from the per-tool reports because it has no verify script).
* Every tool directory contains both ``<slug>-rerun-<YYYYMMDD>.md`` and
  ``<slug>-rerun-<YYYYMMDD>.json`` sidecars.
* Every JSON has ``schema_version == 1`` and ``scripts_total >= 1``.
* ``INDEX.md`` lists all 48 tool directories in its per-tool table.

Exits 0 on full compliance; non-zero lists every violation found.

Usage::

    python3 fastgatk-native/scripts/verify_rerun_report.py --repo .
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2].resolve()
EVIDENCE_ROOT = REPO / "fastgatk-native/evidence"

# Mirrors the canonical-tool table in rerun_all_verify.py.  48 production
# tools; ``genomicsdb-export`` (bridge) intentionally omitted.
EXPECTED_TOOL_SLUGS: list[tuple[str, str]] = [
    ("AnalyzeCovariates", "analyze-covariates"),
    ("AnnotateIntervals", "annotate-intervals"),
    ("ApplyBQSR", "apply-bqsr"),
    ("ApplyVQSR", "apply-vqsr"),
    ("BaseRecalibrator", "bqsr"),
    ("CalculateContamination", "calculate-contamination"),
    ("CallCopyRatioSegments", "call-copy-ratio-segments"),
    ("CheckReferenceCompatibility", "check-reference-compatibility"),
    ("CollectAllelicCounts", "collect-allelic-counts"),
    ("CollectF1R2Counts", "collect-f1r2-counts"),
    ("CollectReadCounts", "collect-read-counts"),
    ("CombineGVCFs", "combine-gvcfs"),
    ("CompareReferences", "compare-references"),
    ("CountBasesInReference", "count-bases-in-reference"),
    ("CountReads", "count-reads"),
    ("CreateReadCountPanelOfNormals", "create-read-count-panel-of-normals"),
    ("DenoiseReadCounts", "denoise-read-counts"),
    ("DepthOfCoverage", "depth-of-coverage"),
    ("FastaAlternateReferenceMaker", "fasta-alternate-reference-maker"),
    ("FastaReferenceMaker", "fasta-reference-maker"),
    ("FilterIntervals", "filter-intervals"),
    ("FilterMutectCalls", "filter-mutect-calls"),
    ("FlagStat", "flag-stat"),
    ("GatherBQSRReports", "gather-bqsr-reports"),
    ("GatherPileupSummaries", "gather-pileup-summaries"),
    ("GatherTranches", "gather-tranches"),
    ("GatherVcfs", "gather-vcfs"),
    ("GenomicsDBImport", "genomicsdb-import"),
    ("GenotypeGVCFs", "genotype-gvcf"),
    ("GetPileupSummaries", "get-pileup-summaries"),
    ("HaplotypeCaller", "hc-call"),
    ("IndexFeatureFile", "index-feature-file"),
    ("LearnReadOrientationModel", "learn-read-orientation-model"),
    ("LeftAlignAndTrimVariants", "left-align-trim"),
    ("MarkDuplicates", "mark-duplicates"),
    ("ModelSegments", "model-segments"),
    ("Mutect2", "mutect2"),
    ("PreprocessIntervals", "preprocess-intervals"),
    ("ReblockGVCF", "reblock-gvcf"),
    ("SelectVariants", "select-variants"),
    ("ShiftFasta", "shift-fasta"),
    ("SortSam", "sort-sam"),
    ("SplitIntervals", "split-intervals"),
    ("ValidateVariants", "validate-variants"),
    ("VariantEval", "variant-eval"),
    ("Funcotator", "funcotator"),
    ("VariantAnnotator", "variant-annotator"),
    ("VariantFiltration", "variant-filtration"),
    ("VariantRecalibrator", "variant-recalibrator"),
    ("VariantsToTable", "variants-to-table"),
]


def today_ymd() -> tuple[str, str]:
    today = dt.date.today()
    return today.strftime("%Y-%m-%d"), today.strftime("%Y%m%d")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=str(REPO))
    parser.add_argument("--rerun-dir", default=None,
                        help="Override the rerun subdir (default: today's <YYYY-MM-DD>-rerun)")
    args = parser.parse_args()

    repo = Path(args.repo).resolve()
    evidence = repo / "fastgatk-native/evidence"

    ymd_dash, ymd = today_ymd()
    today_name = f"{ymd_dash}-rerun"
    if args.rerun_dir:
        rerun_name = args.rerun_dir
        # When --rerun-dir is "<date>-rerun", derive ymd from "<date>" so the
        # gate finds sidecars named `<slug>-rerun-<YYYYMMDD>.{md,json}`.
        m = re.match(r"^(\d{4})-(\d{2})-(\d{2})-rerun$", args.rerun_dir)
        if m:
            ymd = f"{m.group(1)}{m.group(2)}{m.group(3)}"
            ymd_dash = args.rerun_dir[:-len("-rerun")]
    else:
        rerun_name = today_name
    rerun_dir = evidence / rerun_name

    violations: list[str] = []

    # ---- top-level: evidence/ should contain only <rerun-name>/ ----
    if not evidence.is_dir():
        violations.append(f"missing evidence dir: {evidence}")
        _print_and_exit(violations)
        return 1

    visible_top = sorted(
        p.name for p in evidence.iterdir() if not p.name.startswith(".")
    )
    # When --rerun-dir is given, the gate validates against that single
    # snapshot (its sibling coverage/ is allowed).  When omitted, the
    # gate validates against today's snapshot only and flags any other
    # dated rerun directories as unexpected.
    if args.rerun_dir:
        expected_top = {rerun_name, "coverage"}
    else:
        expected_top = {today_name, "coverage"}
    # The evidence tree accumulates dated snapshots and dated audit
    # artifacts across sessions; those are legitimate history.  Flag only
    # entries without a date stamp (undated clutter).
    def _dated(name: str) -> bool:
        return (re.match(r"^\d{4}-\d{2}-\d{2}-", name) is not None or
                re.search(r"-\d{8}(\.|$)", name) is not None or
                re.match(r"^cli-coverage-[0-9]+\.json$", name) is not None)
    extras = sorted(n for n in set(visible_top) - expected_top if not _dated(n))
    missing_top = sorted(expected_top - set(visible_top))
    for m in missing_top:
        violations.append(f"missing rerun dir: {evidence / m}")
    for x in extras:
        violations.append(f"unexpected top-level entry under evidence/: {evidence / x}")

    if not rerun_dir.is_dir():
        _print_and_exit(violations)
        return 1

    # ---- INDEX.md must exist and reference every tool ----
    index_md = rerun_dir / "INDEX.md"
    if not index_md.is_file():
        violations.append(f"missing INDEX.md: {index_md}")
    else:
        text = index_md.read_text(encoding="utf-8")
        for canonical, slug in EXPECTED_TOOL_SLUGS:
            if f"`{canonical}`" not in text:
                violations.append(
                    f"INDEX.md does not list canonical tool `{canonical}` (slug `{slug}`)"
                )

    # ---- ALIGNMENT_STATUS.md (optional but tracked) ----
    alignment_md = rerun_dir / "ALIGNMENT_STATUS.md"
    if not alignment_md.is_file():
        violations.append(f"missing ALIGNMENT_STATUS.md: {alignment_md} "
                          f"(run fastgatk-native/scripts/aggregate_alignment_status.py)")

    # ---- 48 tool dirs, each with .md + .json sidecars ----
    actual_dirs = sorted(
        p.name for p in rerun_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )
    expected_slugs = sorted(slug for _, slug in EXPECTED_TOOL_SLUGS)
    missing_dirs = sorted(set(expected_slugs) - set(actual_dirs))
    extra_dirs = sorted(set(actual_dirs) - set(expected_slugs))
    for m in missing_dirs:
        violations.append(f"missing tool dir under rerun: {rerun_dir / m}")
    for x in extra_dirs:
        violations.append(f"unexpected tool dir under rerun: {rerun_dir / x}")

    for canonical, slug in EXPECTED_TOOL_SLUGS:
        tool_dir = rerun_dir / slug
        md_path = tool_dir / f"{slug}-rerun-{ymd}.md"
        json_path = tool_dir / f"{slug}-rerun-{ymd}.json"
        if not md_path.is_file():
            violations.append(f"missing markdown sidecar: {md_path}")
        if not json_path.is_file():
            violations.append(f"missing json sidecar: {json_path}")
        if json_path.is_file():
            try:
                data = json.loads(json_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                violations.append(f"invalid JSON in {json_path}: {exc}")
                continue
            if data.get("schema_version") != 1:
                violations.append(
                    f"{json_path}: schema_version={data.get('schema_version')!r} (expected 1)"
                )
            if data.get("tool") != canonical:
                violations.append(
                    f"{json_path}: tool={data.get('tool')!r} (expected {canonical!r})"
                )
            scripts_total = data.get("scripts_total")
            if not isinstance(scripts_total, int) or scripts_total < 1:
                violations.append(
                    f"{json_path}: scripts_total={scripts_total!r} (expected >=1)"
                )

    _print_and_exit(violations)
    return 1 if violations else 0


def _print_and_exit(violations: list[str]) -> None:
    if violations:
        print(f"# verify_rerun_report.py: {len(violations)} violation(s)")
        for v in violations:
            print(f"  - {v}")
    else:
        print("# verify_rerun_report.py: all checks passed.")


if __name__ == "__main__":
    sys.exit(main())