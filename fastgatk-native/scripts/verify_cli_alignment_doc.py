#!/usr/bin/env python3
"""Regression verifier for `fastgatk-native/docs/cli-alignment.md`.

Runs the generator, then parses the produced document and asserts:

  1. The doc file exists at the canonical location.
  2. The "Overview" table has exactly one row per native binary declared in
     `fastgatk-native/CMakeLists.txt` (excluding smoke/test/oracle binaries).
  3. Every row in the overview references a GATK JSON file that exists on
     disk under `third_party/gatk-package/gatk-4.6.2.0/gatkdoc/` (or
     explicitly carries the `<none>` marker for native-only bridge binaries).
  4. Every per-tool section contains at least one `accepted` row, or the
     section is a documented native-only bridge binary with no GATK JSON.
  5. The status enum counts (`accepted`/`ignored_via_catalog`/`unsupported`)
     are surfaced in the per-row summary string in the overview table.
  6. The "Bit-identical / bounded-parity evidence" section is present,
     cross-references every overview row, and shows pinned GATK jar +
     fixture-digest artefacts.

Exits 0 on success, non-zero on any violation.  All errors are reported with
the offending binary name and the violation.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

REPO_DEFAULT = Path(__file__).resolve().parents[2]

# Smoke / test / oracle binaries — same exclusion list as the generator.
EXCLUDED_BINARIES = {
    "hts-reader-test",
    "hc-smoke",
    "calling-reference-boundary-smoke",
    "haplotype-cigar-trim-smoke",
    "calling-indel-smoke",
    "fmc-trajectory-oracle",
    "contamination-kernel-segmenter-contract",
}

OVERVIEW_TABLE_RE = re.compile(
    r"^\|\s*`fastgatk-(?P<bin>[a-z0-9-]+)`\s*\|"
    r"\s*`(?P<gatk>[^`]+)`\s*\|"
    r"\s*`(?P<json>[^`]+)`\s*\|"
    r"\s*`(?P<src>[^`]+)`\s*\|"
    r"\s*(?P<accepted>\d+)\s*\|"
    r"\s*(?P<summary>[^|]+?)\s*\|",
    re.MULTILINE,
)

TOOL_SECTION_RE = re.compile(
    r"^##\s+`fastgatk-(?P<bin>[a-z0-9-]+)`",
    re.MULTILINE,
)

STATUS_CELL_RE = re.compile(r"\|\s*`(?P<status>accepted|ignored_via_catalog|native_only|unsupported)`\s*\|")

PARITY_SECTION_RE = re.compile(
    r"^##\s+Bit-identical\s*/\s*bounded-parity evidence\s*$",
    re.MULTILINE,
)
PARITY_TABLE_RE = re.compile(
    r"^\|\s*`fastgatk-(?P<bin>[a-z0-9-]+)`\s*\|"
    r"\s*`(?P<gatk>[^`]+)`\s*\|"
    r"\s*`(?P<status>[^`]+)`\s*\|"
    r"\s*`(?P<det>[^`]+)`\s*\|"
    r"\s*(?P<oracle>\d+)\s*\|"
    r"\s*(?P<bi>\d+)\s*\|"
    r"\s*(?P<bounded>\d+)\s*\|",
    re.MULTILINE,
)
PINNED_JAR_ROW_RE = re.compile(
    r"\|\s*Pinned GATK 4\.6\.2\.0 jar\s*\|\s*([^|]+?)\s*\|",
    re.MULTILINE,
)
FIXTURE_ROW_RE = re.compile(
    r"\|\s*Pinned fixture digests\s*\|\s*([^|]+?)\s*\|",
    re.MULTILINE,
)


def discover_binaries(cmakelists: Path) -> List[str]:
    """Delegate to the generator's parser so we see exactly the same set of
    native binaries (including those emitted from CMake ``foreach`` loops).
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import generate_cli_alignment_doc as gen  # noqa: E402

    raw = gen.parse_cmakelists(cmakelists)
    return [name for name, _ in raw if name not in EXCLUDED_BINARIES]


def parse_overview_table(doc: str) -> List[Tuple[str, str, str, str, int, str]]:
    """Parse ONLY the Overview table at the top of the doc.

    Both the Overview and the parity-evidence section use ``fastgatk-*`` as
    the first column; the regex would otherwise pick up parity rows.  We
    explicitly slice the doc to the Overview heading range.
    """
    overview_start = re.search(r"^##\s+Overview\s*$", doc, re.MULTILINE)
    overview_end = re.search(
        r"^##\s+(Bit-identical / bounded-parity evidence|Per-tool details)\s*$",
        doc,
        re.MULTILINE,
    )
    if not overview_start:
        return []
    start = overview_start.end()
    end = overview_end.start() if overview_end else len(doc)
    slice_ = doc[start:end]
    rows: List[Tuple[str, str, str, str, int, str]] = []
    for m in OVERVIEW_TABLE_RE.finditer(slice_):
        rows.append(
            (
                m.group("bin"),
                m.group("gatk"),
                m.group("json"),
                m.group("src"),
                int(m.group("accepted")),
                m.group("summary").strip(),
            )
        )
    return rows


def parse_per_tool_status(doc: str) -> Dict[str, Dict[str, int]]:
    sections: Dict[str, Dict[str, int]] = {}
    # Walk the doc by finding each "## fastgatk-..." heading and slicing the
    # text from that heading until the next heading.
    matches = list(TOOL_SECTION_RE.finditer(doc))
    for i, m in enumerate(matches):
        bin_name = m.group("bin")
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(doc)
        body = doc[start:end]
        counts = {k: 0 for k in ("accepted", "ignored_via_catalog", "native_only", "unsupported")}
        for s in STATUS_CELL_RE.finditer(body):
            counts[s.group("status")] += 1
        sections[bin_name] = counts
    return sections


def parse_parity_table(doc: str) -> List[Tuple[str, str, str, str, int, int, int, str]]:
    """Parse the per-binary parity evidence table."""
    rows: List[Tuple[str, str, str, str, int, int, int, str]] = []
    for m in PARITY_TABLE_RE.finditer(doc):
        rows.append(
            (
                m.group("bin"),
                m.group("gatk"),
                m.group("status"),
                m.group("det"),
                int(m.group("oracle")),
                int(m.group("bi")),
                int(m.group("bounded")),
                # "sample bit-identical script" is the last cell; we don't
                # extract it here because the table allows `—`.
            )
        )
    return rows


def parse_pinned_artefacts(doc: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (pinned_jar_status, fixture_digest_status) from the parity
    section's "Pinned artifacts" table.
    """
    jar_match = PINNED_JAR_ROW_RE.search(doc)
    fix_match = FIXTURE_ROW_RE.search(doc)
    return (
        jar_match.group(1).strip() if jar_match else None,
        fix_match.group(1).strip() if fix_match else None,
    )


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=str(REPO_DEFAULT), help="Path to repository root.")
    ap.add_argument(
        "--skip-regenerate",
        action="store_true",
        help="Skip re-running the generator and verify the existing doc.",
    )
    args = ap.parse_args(argv)
    repo = Path(args.repo).resolve()
    gen_script = repo / "fastgatk-native" / "scripts" / "generate_cli_alignment_doc.py"
    doc_path = repo / "fastgatk-native" / "docs" / "cli-alignment.md"
    cmakelists = repo / "fastgatk-native" / "CMakeLists.txt"
    gatkdoc_dir = repo / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatkdoc"

    if not gen_script.exists():
        print(f"ERROR: generator script not found at {gen_script}", file=sys.stderr)
        return 2
    if not cmakelists.exists():
        print(f"ERROR: CMakeLists.txt not found at {cmakelists}", file=sys.stderr)
        return 2
    if not gatkdoc_dir.exists():
        print(f"ERROR: gatkdoc dir not found at {gatkdoc_dir}", file=sys.stderr)
        return 2

    # ----- Regenerate the doc (idempotent).
    if not args.skip_regenerate:
        rc = subprocess.run(
            [sys.executable, str(gen_script), "--repo", str(repo)],
            check=False,
        )
        if rc.returncode != 0:
            print("ERROR: generator failed", file=sys.stderr)
            return rc.returncode

    if not doc_path.exists():
        print(f"ERROR: doc file not found at {doc_path}", file=sys.stderr)
        return 2
    doc = doc_path.read_text()

    # ----- Assertion 1: doc exists (already checked).

    # ----- Assertion 2: overview table covers every native binary.
    expected_bins = discover_binaries(cmakelists)
    rows = parse_overview_table(doc)
    seen_bins = {r[0] for r in rows}
    errors: List[str] = []
    missing = [b for b in expected_bins if b not in seen_bins]
    extra = sorted(seen_bins - set(expected_bins))
    if missing:
        errors.append(f"missing overview rows for binaries: {missing}")
    if extra:
        errors.append(f"extra overview rows for binaries: {extra}")
    if len(rows) != len(expected_bins):
        errors.append(
            f"overview row count mismatch: doc has {len(rows)} rows, "
            f"expected {len(expected_bins)} (one per native binary in CMakeLists.txt, "
            f"excluding smoke/test/oracle binaries)"
        )

    # ----- Assertion 3: GATK JSON files exist on disk.
    for bin_name, gatk, json_ref, src, accepted, summary in rows:
        if json_ref == "<none>":
            # native-only bridge binary — must be flagged as such in the
            # per-row summary.
            if "native-only bridge" not in summary.lower():
                errors.append(
                    f"{bin_name}: json_ref is <none> but summary lacks 'native-only bridge': {summary!r}"
                )
            continue
        path = gatkdoc_dir / json_ref
        if not path.exists():
            errors.append(f"{bin_name}: gatkdoc JSON missing on disk: {path}")

    # ----- Assertion 4: every section has >= 1 accepted row, or it is a
    # native-only bridge with no gatkdoc JSON.
    sections = parse_per_tool_status(doc)
    for bin_name, _, json_ref, _, accepted, _ in rows:
        sec = sections.get(bin_name)
        if sec is None:
            errors.append(f"{bin_name}: missing per-tool section")
            continue
        if json_ref == "<none>":
            continue
        if accepted < 1:
            errors.append(f"{bin_name}: overview reports 0 accepted but JSON exists ({json_ref})")

    # ----- Assertion 5: per-row summary contains the status counts for non-bridge rows.
    for bin_name, _, json_ref, _, accepted, summary in rows:
        if json_ref == "<none>":
            continue
        if f"`accepted`={accepted}" not in summary:
            errors.append(
                f"{bin_name}: overview summary does not report `accepted`={accepted}: {summary!r}"
            )

    # ----- Assertion 6: parity section present + cross-references every row.
    if not PARITY_SECTION_RE.search(doc):
        errors.append(
            "missing '## Bit-identical / bounded-parity evidence' section"
        )
    parity_rows = parse_parity_table(doc)
    parity_bins = {r[0] for r in parity_rows}
    missing_in_parity = [
        bin_name for bin_name, _, _, _, _, _ in rows if bin_name not in parity_bins
    ]
    extra_in_parity = sorted(parity_bins - {b for b, *_ in rows})
    if missing_in_parity:
        errors.append(
            f"parity table is missing rows for binaries: {missing_in_parity}"
        )
    if extra_in_parity:
        errors.append(
            f"parity table has extra rows for binaries: {extra_in_parity}"
        )
    # Pinned artefacts section must report both GATK jar + fixture digest.
    jar_status, fixture_status = parse_pinned_artefacts(doc)
    if jar_status is None:
        errors.append("parity section: 'Pinned GATK 4.6.2.0 jar' row missing")
    elif "present" not in jar_status.lower() and "missing" not in jar_status.lower():
        errors.append(
            f"parity section: GATK jar status unclear: {jar_status!r}"
        )
    elif "missing" in jar_status.lower():
        errors.append(f"parity section: GATK jar missing: {jar_status!r}")
    if fixture_status is None:
        errors.append("parity section: 'Pinned fixture digests' row missing")
    elif "missing" in fixture_status.lower():
        errors.append(
            f"parity section: pinned fixture digests missing: {fixture_status!r}"
        )
    # Cross-check that on-disk dispatcher registry is consistent with what
    # the parity table says (no `contract-compatible` tool may appear as
    # `—` in the parity table unless the registry doesn't know about it).
    registry_path = repo / "fastgatk-native" / "dispatcher" / "tool_registry.json"
    if registry_path.exists():
        with registry_path.open() as f:
            reg = json.load(f).get("tools", {})
        for bin_name, gatk_class, _, _, _, _ in rows:
            canon = gatk_class.split(" ")[0]
            parity_row = next((r for r in parity_rows if r[0] == bin_name), None)
            if parity_row is None:
                continue
            parity_status = parity_row[2]
            reg_status = reg.get(canon, {}).get("status")
            if reg_status and reg_status not in parity_status:
                errors.append(
                    f"{bin_name}: parity table status '{parity_status}' does not "
                    f"match dispatcher registry '{reg_status}' for {canon}"
                )

    # ----- Print per-binary counts table for visibility.
    print("== Per-binary status counts (from doc) ==")
    print(f"{'binary':48} {'accepted':>9} {'ignored_via_catalog':>19} {'unsupported':>12} {'native_only':>12}")
    for bin_name in sorted(sections):
        c = sections[bin_name]
        print(
            f"fastgatk-{bin_name:42} {c['accepted']:>9} {c['ignored_via_catalog']:>19} "
            f"{c['unsupported']:>12} {c['native_only']:>12}"
        )

    if errors:
        print("\nERRORS:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1
    print("\nOK")
    return 0


if __name__ == "__main__":
    sys.exit(main())