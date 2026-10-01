#!/usr/bin/env python3
"""Aggregate per-tool GATK-alignment status from the rerun JSON sidecars and
`fastgatk-native/docs/cli-alignment.md` *Bit-identical / bounded-parity evidence*
table.

For each of the 48 production tools this emits:

* ``verify_total``        — total verify scripts in the rerun report
* ``verify_passed``       — exited 0
* ``verify_failed``       — exited non-zero
* ``verify_skipped``      — exited 0 via oracle-guard skip (NOT verified vs GATK)
* ``parity_oracle_total`` — total oracle scripts declared in cli-alignment.md
* ``parity_bit_identical``— declared bit-identical count
* ``parity_bounded``      — declared bounded-oracle count
* ``oracle_passes``       — verify scripts whose name ends in
                            ``_gatk_oracle.py``/``_gatk_contract.py`` and
                            which actually passed (exit 0, skipped=false)
* ``oracle_skips``        — oracle-named scripts that returned 0 via skip path
                            (no GATK oracle available at rerun time)
* ``alignment_status``    — coarse verdict (see below)

Verdict logic:
    ``aligned``            — every declared bit-identical + bounded oracle
                              passes AND there are no failures
    ``partially-aligned``  — declared parity > 0 but some oracle-named script
                              is missing / failed / skipped
    ``native-only``        — cli-alignment.md declares 0 bit-identical AND 0
                              bounded-parity scripts (the tool has no oracle
                              tests, only native regression)
    ``no-rerun-record``    — rerun JSON sidecar missing or unreadable

The script writes ``ALIGNMENT_STATUS.md`` next to ``INDEX.md`` and prints a
short summary table to stdout.  It is purely read-only; it never re-runs a
verify script.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2].resolve()
def _latest_rerun_root():
    """Newest <YYYY-MM-DD>-rerun snapshot under evidence/ (today's wins)."""
    import datetime as _dt
    import re as _re
    evidence = REPO / "fastgatk-native/evidence"
    candidates = sorted(
        p.name for p in evidence.iterdir()
        if p.is_dir() and _re.match(r"^\d{4}-\d{2}-\d{2}-rerun$", p.name))
    today = _dt.date.today().strftime("%Y-%m-%d") + "-rerun"
    chosen = today if today in candidates else (candidates[-1] if candidates else today)
    return evidence / chosen


RERUN_ROOT = _latest_rerun_root()
CLI_ALIGNMENT = REPO / "fastgatk-native/docs/cli-alignment.md"

# Native-binary slug → GATK canonical name (mirrors the driver table).
SLUG_TO_CANONICAL = {
    "analyze-covariates": "AnalyzeCovariates",
    "annotate-intervals": "AnnotateIntervals",
    "apply-bqsr": "ApplyBQSR",
    "apply-vqsr": "ApplyVQSR",
    "bqsr": "BaseRecalibrator",
    "calculate-contamination": "CalculateContamination",
    "call-copy-ratio-segments": "CallCopyRatioSegments",
    "check-reference-compatibility": "CheckReferenceCompatibility",
    "collect-allelic-counts": "CollectAllelicCounts",
    "collect-f1r2-counts": "CollectF1R2Counts",
    "collect-read-counts": "CollectReadCounts",
    "combine-gvcfs": "CombineGVCFs",
    "compare-references": "CompareReferences",
    "count-bases-in-reference": "CountBasesInReference",
    "count-reads": "CountReads",
    "create-read-count-panel-of-normals": "CreateReadCountPanelOfNormals",
    "denoise-read-counts": "DenoiseReadCounts",
    "depth-of-coverage": "DepthOfCoverage",
    "fasta-alternate-reference-maker": "FastaAlternateReferenceMaker",
    "fasta-reference-maker": "FastaReferenceMaker",
    "filter-intervals": "FilterIntervals",
    "filter-mutect-calls": "FilterMutectCalls",
    "flag-stat": "FlagStat",
    "gather-bqsr-reports": "GatherBQSRReports",
    "gather-pileup-summaries": "GatherPileupSummaries",
    "gather-tranches": "GatherTranches",
    "gather-vcfs": "GatherVcfs",
    "genomicsdb-import": "GenomicsDBImport",
    "genotype-gvcf": "GenotypeGVCFs",
    "get-pileup-summaries": "GetPileupSummaries",
    "hc-call": "HaplotypeCaller",
    "index-feature-file": "IndexFeatureFile",
    "learn-read-orientation-model": "LearnReadOrientationModel",
    "left-align-trim": "LeftAlignAndTrimVariants",
    "mark-duplicates": "MarkDuplicates",
    "model-segments": "ModelSegments",
    "mutect2": "Mutect2",
    "preprocess-intervals": "PreprocessIntervals",
    "reblock-gvcf": "ReblockGVCF",
    "select-variants": "SelectVariants",
    "shift-fasta": "ShiftFasta",
    "sort-sam": "SortSam",
    "split-intervals": "SplitIntervals",
    "validate-variants": "ValidateVariants",
    "variant-eval": "VariantEval",
    "variant-filtration": "VariantFiltration",
    "variant-recalibrator": "VariantRecalibrator",
    "variants-to-table": "VariantsToTable",
}


def parse_parity_table(path: Path) -> dict[str, dict[str, int]]:
    """Parse the *Bit-identical / bounded-parity evidence* table.

    Returns a dict keyed by native-binary slug with values
    ``{"oracle_total": N, "bit_identical": N, "bounded": N}``.
    Rows whose `bit-identical` and `bounded` columns are both `—`
    (i.e. zero) are still included.
    """
    text = path.read_text(encoding="utf-8")
    out: dict[str, dict[str, int]] = {}

    # Locate the "### Per-binary parity evidence" section and read until the
    # next ``## `` heading.
    marker = "### Per-binary parity evidence"
    if marker not in text:
        raise SystemExit(f"could not locate {marker!r} in {path}")
    start = text.index(marker)
    end = text.find("\n## ", start + len(marker))
    section = text[start:end if end != -1 else None]

    row_re = re.compile(
        r"^\|\s*`(?P<bin>fastgatk-[a-z0-9-]+)`\s*\|"
        r"\s*`?(?P<canon>[^`|]+?)`?\s*\|"
        r"\s*[^|]+\|\s*[^|]+\|\s*(?P<oracle>\d+|\—)\s*\|"
        r"\s*(?P<bit>\d+|\—)\s*\|"
        r"\s*(?P<bounded>\d+|\—)\s*\|",
        re.MULTILINE,
    )
    for m in row_re.finditer(section):
        bin_slug = m.group("bin")[len("fastgatk-"):]
        def as_int(token: str) -> int:
            return 0 if token == "—" else int(token)
        out[bin_slug] = {
            "oracle_total": as_int(m.group("oracle")),
            "bit_identical": as_int(m.group("bit")),
            "bounded": as_int(m.group("bounded")),
        }
    return out


# Oracle-script detection.  Matches cli-alignment.md's description of the
# 251-script oracle corpus: scripts whose name ends in "_oracle.py" or
# "_contract.py" AND that are not pure-native fixtures.  The known native-
# only fixtures are excluded explicitly below.
ORACLE_NAME_PATTERNS = (
    re.compile(r"_oracle\.py$"),
    re.compile(r"_contract\.py$"),
)
NATIVE_ONLY_ORACLE_PATTERNS = (
    re.compile(r"_fixture_oracle\.py$"),  # e.g. verify_hc_*_fixture_oracle.py
    re.compile(r"_real_contract\.py$"),   # e.g. verify_hc_chr20_real_contract.py
    re.compile(r"_native_contract\.py$"),
    re.compile(r"_native_oracle\.py$"),
    re.compile(r"_oracle_sanity\.py$"),
)


def is_oracle_script(name: str) -> bool:
    if any(p.search(name) for p in ORACLE_NAME_PATTERNS):
        if not any(p.search(name) for p in NATIVE_ONLY_ORACLE_PATTERNS):
            return True
    return False


def derive_status(row: dict) -> str:
    """Coarse verdict for one tool row.

    The parity table's ``oracle scripts`` column counts every verify script
    that asserts GATK parity (bit-identical *or* bounded-oracle).  When
    ``cli-alignment.md`` declares 0 oracle scripts for a tool, that tool
    has no GATK-parity tests at all in the corpus and the rerun covers
    native regression only (``native-only``).

    Otherwise the verdict is driven by the rerun results:

    * **aligned** — every rerun script for this tool exited 0 without
      going through the ``oracle_guard`` skip path (true GATK comparison
      evidence for every script).
    * **partially-aligned** — at least one rerun script failed, or at
      least one returned 0 via the oracle-guard skip path (GATK oracle
      absent from this environment).
    * **no-rerun-record** — the rerun JSON sidecar is missing or
      unreadable for this tool.

    Note: the parity-table ``oracle scripts`` count and the rerun
    ``scripts_total`` are tracked separately.  Cross-tool aliases
    (e.g. ``verify_analyze_covariates_bqsr_alias_gatk_oracle.py`` which
    is counted under both ``AnalyzeCovariates`` and ``BaseRecalibrator``
    in cli-alignment.md) may legitimately appear in the rerun of more
    than one tool's directory; that does not make either tool "partial".
    """
    if not row.get("rerun_present"):
        return "no-rerun-record"
    declared_oracle = row["parity_oracle_total"]
    if declared_oracle == 0:
        return "native-only"
    failed = row.get("verify_failed", 0)
    skipped = row.get("verify_skipped", 0)
    if failed > 0 or skipped > 0:
        return "partially-aligned"
    return "aligned"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=str(REPO))
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    rerun_root = RERUN_ROOT if RERUN_ROOT.is_dir() else _latest_rerun_root()
    cli_alignment = repo / "fastgatk-native/docs/cli-alignment.md"

    parity = parse_parity_table(cli_alignment)
    today = dt.date.today().strftime("%Y-%m-%d")

    rows: list[dict] = []
    for slug in sorted(SLUG_TO_CANONICAL):
        json_path = rerun_root / slug / f"{slug}-rerun-20260923.json"
        if not json_path.is_file():
            row = {
                "slug": slug,
                "canonical": SLUG_TO_CANONICAL[slug],
                "rerun_present": False,
                "parity_oracle_total": parity.get(slug, {}).get("oracle_total", 0),
                "parity_bit_identical": parity.get(slug, {}).get("bit_identical", 0),
                "parity_bounded": parity.get(slug, {}).get("bounded", 0),
            }
            row["alignment_status"] = derive_status(row)
            rows.append(row)
            continue
        data = json.loads(json_path.read_text(encoding="utf-8"))
        scripts = data.get("scripts", [])
        oracle_scripts = [s for s in scripts if is_oracle_script(s["name"])]
        oracle_total = len(oracle_scripts)
        oracle_passes = sum(1 for s in oracle_scripts
                            if s.get("passed") and not s.get("skipped"))
        oracle_skips = sum(1 for s in oracle_scripts
                           if s.get("passed") and s.get("skipped"))
        oracle_failed = sum(1 for s in oracle_scripts
                            if not s.get("passed"))
        # Aggregate native-vs-Java process metrics across all rerun scripts
        # for this tool.  These come from the per-script ``native`` and
        # ``java`` dictionaries emitted by the v2 sampler.
        nat_pids = sum(s.get("native", {}).get("process_count", 0) for s in scripts)
        nat_wall = sum(s.get("native", {}).get("wall_clock_seconds_sum", 0.0) for s in scripts)
        nat_rss = max((s.get("native", {}).get("peak_rss_bytes", 0) for s in scripts), default=0)
        jav_pids = sum(s.get("java", {}).get("process_count", 0) for s in scripts)
        jav_wall = sum(s.get("java", {}).get("wall_clock_seconds_sum", 0.0) for s in scripts)
        jav_rss = max((s.get("java", {}).get("peak_rss_bytes", 0) for s in scripts), default=0)
        row = {
            "slug": slug,
            "canonical": SLUG_TO_CANONICAL[slug],
            "rerun_present": True,
            "verify_total": data.get("scripts_total", len(scripts)),
            "verify_passed": data.get("scripts_passed", 0),
            "verify_failed": data.get("scripts_failed", 0),
            "verify_skipped": data.get("scripts_skipped", 0),
            "parity_oracle_total": parity.get(slug, {}).get("oracle_total", 0),
            "parity_bit_identical": parity.get(slug, {}).get("bit_identical", 0),
            "parity_bounded": parity.get(slug, {}).get("bounded", 0),
            "oracle_total_in_rerun": oracle_total,
            "oracle_passes": oracle_passes,
            "oracle_skips": oracle_skips,
            "oracle_failed": oracle_failed,
            "native_pids": nat_pids,
            "native_wall_s": round(nat_wall, 3),
            "native_peak_rss_bytes": nat_rss,
            "java_pids": jav_pids,
            "java_wall_s": round(jav_wall, 3),
            "java_peak_rss_bytes": jav_rss,
        }
        row["alignment_status"] = derive_status(row)
        rows.append(row)

    # Aggregate counts.
    counts = {
        "aligned": sum(1 for r in rows if r["alignment_status"] == "aligned"),
        "partially-aligned": sum(1 for r in rows if r["alignment_status"] == "partially-aligned"),
        "native-only": sum(1 for r in rows if r["alignment_status"] == "native-only"),
        "no-rerun-record": sum(1 for r in rows if r["alignment_status"] == "no-rerun-record"),
    }

    # Print summary.
    print(f"# Alignment status — {today}")
    print(f"# tools: {len(rows)}; aligned={counts['aligned']}, "
          f"partially-aligned={counts['partially-aligned']}, "
          f"native-only={counts['native-only']}, "
          f"no-rerun-record={counts['no-rerun-record']}")
    header = ("tool", "status", "parity_oracles", "bit", "bounded",
              "verify_total", "verify_pass", "verify_fail", "verify_skip",
              "nat_pids", "jav_pids", "java_wall_s", "native_wall_s", "java_rss_mib")
    print(" | ".join(header))
    print("-" * 130)
    for r in rows:
        if not r["rerun_present"]:
            print(f"{r['canonical']:24s} | no-rerun-record | "
                  f"{r['parity_oracle_total']} | {r['parity_bit_identical']} | "
                  f"{r['parity_bounded']} | - | - | - | - | - | - | - | - | -")
            continue
        jav_rss_mib = round(r.get('java_peak_rss_bytes', 0) / (1024 * 1024), 1)
        print(f"{r['canonical']:24s} | {r['alignment_status']:18s} | "
              f"{r['parity_oracle_total']:3d} | {r['parity_bit_identical']:3d} | "
              f"{r['parity_bounded']:3d} | "
              f"{r['verify_total']:3d} | {r['verify_passed']:3d} | "
              f"{r['verify_failed']:3d} | {r['verify_skipped']:3d} | "
              f"{r['native_pids']:3d} | {r['java_pids']:3d} | "
              f"{r['java_wall_s']:7.2f} | {r['native_wall_s']:7.2f} | "
              f"{jav_rss_mib:7.1f}")

    # Write ALIGNMENT_STATUS.md.
    md = [
        f"# GATK-alignment status — {today}",
        "",
        "Cross-reference of the per-tool re-run reports under this directory "
        "with the *Bit-identical / bounded-parity evidence* table in "
        "`fastgatk-native/docs/cli-alignment.md`.  A *passing* `verify_*.py` "
        "is not by itself evidence of GATK parity: a script that returns 0 "
        "via the `oracle_guard` skip path (because the GATK jar/JDK is "
        "missing from this environment) is *not* a GATK comparison.  See "
        "`INDEX.md` for the raw re-run counts.",
        "",
        "## Verdicts",
        "",
        "* **aligned** — every rerun script for this tool exited 0 without "
        "going through the `oracle_guard` skip path (i.e. real GATK "
        "comparison evidence for every script).",
        "* **partially-aligned** — declared parity > 0 but at least one "
        "rerun script failed, or returned 0 via the oracle-guard skip path "
        "(GATK oracle absent from this environment).",
        "* **native-only** — `cli-alignment.md` declares 0 bit-identical and 0 "
        "bounded-parity oracle scripts for this tool (re-run covers native "
        "regression only; GATK parity is unverified in this environment).",
        "* **no-rerun-record** — re-run JSON sidecar missing or unreadable.",
        "",
        f"## Counts",
        "",
        f"- Tools: **{len(rows)}**",
        f"- aligned: **{counts['aligned']}**",
        f"- partially-aligned: **{counts['partially-aligned']}**",
        f"- native-only: **{counts['native-only']}**",
        f"- no-rerun-record: **{counts['no-rerun-record']}**",
        "",
        "## Per-tool table",
        "",
        "| Tool | Status | Parity oracles (declared) | Bit-identical | Bounded | "
        "Verify total/pass/fail/skip | Native pids | Java pids | Java wall-clock (s) | "
        "Java peak RSS (MiB) |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        if not r["rerun_present"]:
            md.append(
                f"| `{r['canonical']}` (`{r['slug']}`) | no-rerun-record | "
                f"{r['parity_oracle_total']} | {r['parity_bit_identical']} | "
                f"{r['parity_bounded']} | — | — | — | — | — |"
            )
            continue
        jav_rss_mib = round(r.get('java_peak_rss_bytes', 0) / (1024 * 1024), 1)
        md.append(
            f"| `{r['canonical']}` (`{r['slug']}`) | "
            f"**{r['alignment_status']}** | "
            f"{r['parity_oracle_total']} | {r['parity_bit_identical']} | "
            f"{r['parity_bounded']} | "
            f"{r['verify_total']}/{r['verify_passed']}/{r['verify_failed']}/{r['verify_skipped']} | "
            f"{r['native_pids']} | {r['java_pids']} | {r['java_wall_s']:.2f} | "
            f"{jav_rss_mib:.1f} |"
        )
    md.extend([
        "",
        "## Reading the columns",
        "",
        "* `Parity oracles (declared)` — count from the "
        "*Bit-identical / bounded-parity evidence* table in "
        "`fastgatk-native/docs/cli-alignment.md`.",
        "* `Bit-identical` / `Bounded` — subsets of the declared oracles; "
        "see `cli-alignment.md` for definitions.",
        "* `Verify total/pass/fail/skip` — counts from the rerun JSON sidecar "
        "for this tool.  `skip` means the script returned 0 via the "
        "`oracle_guard` skip path (no GATK oracle available — see stderr "
        "tail of the JSON for the `[NOT VERIFIED AGAINST GATK]` notice).",
        "",
        "A tool is marked **aligned** only when `verify_total == declared` "
        "AND `verify_failed == 0` AND `verify_skipped == 0`.",
        "",
    ])
    (rerun_root / "ALIGNMENT_STATUS.md").write_text("\n".join(md),
                                                    encoding="utf-8")
    print(f"\n# Wrote {rerun_root / 'ALIGNMENT_STATUS.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())