#!/usr/bin/env python3
"""Strict oracle for VariantFiltration ``--apply-allele-specific-filters``
AS_FilterStatus encoding (pinned GATK 4.6.2.0).

Observation under test
----------------------
``fastgatk-native/scripts/verify_variant_filtration.py`` used to assert
``AS_FilterStatus=LowASQD,PASS`` / ``AS_FilterStatus=PASS,LowASQD`` (lines
504-505) for the fixture below.  That encoding is native-only: pinned GATK
emits ``AS_FilterStatus=SITE|SITE``, joins alleles with ``|`` and applies
**no** allele filter at all.  This script re-measures both tools on identical
input with identical arguments and gates the FILTER column and the
``AS_FilterStatus`` INFO key.

GATK's rule (read from ``gatk-source/``, then confirmed by measurement)
----------------------------------------------------------------------
1. ``VariantFiltration.apply`` (``VariantFiltration.java:359-365``) is taken
   only when ``--apply-allele-specific-filters`` is set.  It calls
   ``splitMultiAllelics`` (``:371-378``), which builds each single-ALT
   ``VariantContext`` from

       new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(),
                                 vc.getEnd(), Arrays.asList(vc.getReference(), Allele.NO_CALL))

   That builder constructor carries **only** contig/start/stop/alleles: the
   split context has *no INFO attributes and no genotypes*.  Measured with the
   pinned jar: ``vc.hasAttribute("AS_QD")`` is false and
   ``vc.getGenotype("S1")`` is null (``JexlEngine - attempting to call method
   on null``) inside this mode, while the same expressions match without the
   flag.
2. Each split context is filtered by ``filter`` (``:400-438``) and each
   expression is evaluated against that bare context, so a rule that
   dereferences an INFO/Number=A annotation sees JEXL ``null``; JEXL compares
   null as false, so the rule never adds a filter label.  This also holds under
   ``--missing-values-evaluate-as-failing``: the annotation is *present but
   null*, not an undefined property, so the missing-value policy does not
   apply.
3. The per-allele filter sets are then assembled by
   ``AlleleFilterUtils.addAlleleAndSiteFilters`` (``:94-122``).  The record's
   own ``AS_FilterStatus`` is decoded first (``decodeASFilters``, ``:24-28``)
   and only labels for alleles whose rule fired are merged in
   (``addAlleleFilters``, ``:69-82``); with no pre-existing value the slots are
   seeded with the placeholder ``GATKVCFConstants.SITE_LEVEL_FILTERS`` =
   ``"SITE"`` (``:104-106``, constant at ``GATKVCFConstants.java:201``).  A
   pre-existing vector with a different arity leaves the record untouched
   (``:99-102``) and ``--invalidate-previous-filters`` reseeds it (``:104-106``,
   ``:112-114``).
4. Alleles are joined by ``AnnotationUtils.encodeAnyASListWithRawDelim``
   (``ALLELE_SPECIFIC_RAW_DELIM`` = ``"|"``, ``AnnotationUtils.java:21``) and
   the labels of one allele by ``LIST_DELIMITER`` = ``","`` — **not** a comma
   between alleles — and an unfiltered allele keeps the literal ``SITE`` rather
   than ``PASS``.
5. The site FILTER column is the **intersection** of the allele filter sets
   (``:115-117``: ``retainAll`` across alleles) and falls back to PASS when the
   intersection is empty (``:118-120``).  With every allele unfiltered the
   intersection is empty, hence ``FILTER=PASS`` while ``AS_FilterStatus`` is
   still written.

Scope
-----
Every case byte-compares the data-row identity columns ``CHROM..FILTER`` and the
``AS_FilterStatus`` value, and additionally pins the literal GATK expectation
measured for that case so fixture/CLI drift cannot pass silently.  Writer-only
differences outside this subject -- INFO key order and the ``AS_QD`` float
spelling (GATK ``1.0`` vs htslib ``1``) -- are reported, not gated.

Case ``invert-filter-expression`` is *reported only*: there an allele filter
really does fire, so GATK's site FILTER (step 5) comes into play and native
still writes PASS.  That is a separate, previously-scoped divergence in the site
FILTER column, recorded here with its measured rows rather than hidden.

Exit status: 0 when every gated check passes, non-zero otherwise.
``--expect-divergence`` turns the run into a pure diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

BASE_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
    "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=Allele quality by depth>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
)
STATUS_HEADER = (
    "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=Existing per-allele status>\n"
    "##FILTER=<ID=foo,Description=An input filter label>\n"
    "##FILTER=<ID=bar,Description=An input filter label>\n"
)
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

AS_QD_FILTER_ARGS = [
    "--filter-expression", 'vc.getAttribute("AS_QD") < 2',
    "--filter-name", "LowASQD",
    "--apply-allele-specific-filters",
]

TWO_RECORDS = (
    "chr1\t1\t.\tA\tC,G\t50\tPASS\tAS_QD=1.0,3.0\tGT\t0/1\n"
    "chr1\t2\t.\tC\tT,G\t50\tPASS\tAS_QD=3.0,1.0\tGT\t0/1\n"
)
ONE_RECORD = "chr1\t1\t.\tA\tC,G\t50\tPASS\tAS_QD=1.0,3.0\tGT\t0/1\n"

# ``expect`` is the measured (FILTER, AS_FilterStatus) pair per record.
CASES = [
    {
        "case": "as-qd-number-a",
        "why": "the fixture of verify_variant_filtration.py:496-514: a Number=A "
               "INFO rule in allele mode sees the split context's missing INFO "
               "and fires for no allele",
        "body": TWO_RECORDS,
        "extra_header": "",
        "args": AS_QD_FILTER_ARGS,
        "gated": True,
        "expect": [("PASS", "SITE|SITE"), ("PASS", "SITE|SITE")],
    },
    {
        "case": "preexisting-as-filter-status",
        "why": "AlleleFilterUtils.java:98-108 starts from the decoded input "
               "AS_FilterStatus and only adds labels for alleles that fired",
        "body": "chr1\t1\t.\tA\tC,G\t50\tPASS\tAS_QD=1.0,3.0;AS_FilterStatus=foo|bar\tGT\t0/1\n",
        "extra_header": STATUS_HEADER,
        "args": AS_QD_FILTER_ARGS,
        "gated": True,
        "expect": [("PASS", "foo|bar")],
    },
    {
        "case": "invalidate-previous-filters",
        "why": "AlleleFilterUtils.java:104-106 reseeds the SITE placeholder when "
               "--invalidate-previous-filters is set",
        "body": "chr1\t1\t.\tA\tC,G\t50\tPASS\tAS_QD=1.0,3.0;AS_FilterStatus=foo|bar\tGT\t0/1\n",
        "extra_header": STATUS_HEADER,
        "args": AS_QD_FILTER_ARGS + ["--invalidate-previous-filters"],
        "gated": True,
        "expect": [(".", "SITE|SITE")],
    },
    {
        "case": "missing-values-evaluate-as-failing",
        "why": "the split context holds a present-but-null INFO attribute, so "
               "--missing-values-evaluate-as-failing does not turn it into a match",
        "body": ONE_RECORD,
        "extra_header": "",
        "args": AS_QD_FILTER_ARGS + ["--missing-values-evaluate-as-failing"],
        "gated": True,
        "expect": [("PASS", "SITE|SITE")],
    },
    {
        "case": "invert-filter-expression",
        "why": "inverting the predicate makes every allele match; AS_FilterStatus "
               "then carries the real label, but the site FILTER column (GATK's "
               "intersection, AlleleFilterUtils.java:115-120) is a separately "
               "scoped divergence, so this case is reported rather than gated",
        "body": ONE_RECORD,
        "extra_header": "",
        "args": AS_QD_FILTER_ARGS + ["--invert-filter-expression"],
        "gated": False,
        "expect": [("LowASQD", "LowASQD|LowASQD")],
    },
]


def data_rows(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                rows.append(line.rstrip("\n"))
    return rows


def info_map(row: str) -> dict[str, str | None]:
    field = row.split("\t")[7]
    items: dict[str, str | None] = {}
    if field in ("", "."):
        return items
    for item in field.split(";"):
        key, separator, value = item.partition("=")
        items[key] = value if separator else None
    return items


def observed(row: str) -> tuple[str, str | None]:
    return row.split("\t")[6], info_map(row).get("AS_FilterStatus")


def run(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-4000:])
    return result


def run_case(case: dict, work: pathlib.Path, java: pathlib.Path, jar: pathlib.Path,
             native: pathlib.Path, timeout: int) -> dict:
    source = work / f"{case['case']}.vcf"
    source.write_text(BASE_HEADER + case["extra_header"] + COLUMNS + case["body"],
                      encoding="utf-8")
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf.gz"
    gatk_result = run([str(java), "-Xmx1g", "-jar", str(jar), "VariantFiltration",
                       "-V", str(source), "-O", str(gatk_out), *case["args"]],
                      f"GATK VariantFiltration [{case['case']}]", timeout)
    native_result = run([str(native), "-V", str(source), "-O", str(native_out), *case["args"]],
                        f"native VariantFiltration [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    gatk_observed = [observed(row) for row in gatk_rows]
    native_observed = [observed(row) for row in native_rows]

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "args": case["args"],
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_observed": gatk_observed,
        "native_observed": native_observed,
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "violations": [],
    }
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    if gatk_observed != case["expect"]:
        result["violations"].append(
            f"GATK (FILTER, AS_FilterStatus) moved away from the measured truth: "
            f"GATK={gatk_observed} expected={case['expect']}")
    if case["gated"]:
        if len(gatk_observed) != len(native_observed):
            result["violations"].append(
                f"record count differs: GATK={len(gatk_observed)} native={len(native_observed)}")
        for index, (gatk_pair, native_pair) in enumerate(zip(gatk_observed, native_observed)):
            if gatk_pair[0] != native_pair[0]:
                result["violations"].append(
                    f"row {index}: FILTER differs: GATK={gatk_pair[0]!r} NATIVE={native_pair[0]!r}")
            if gatk_pair[1] != native_pair[1]:
                result["violations"].append(
                    f"row {index}: AS_FilterStatus differs: "
                    f"GATK={gatk_pair[1]!r} NATIVE={native_pair[1]!r}")
    for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows)):
        if info_map(gatk_row) != info_map(native_row):
            result.setdefault("reported_info_differences", []).append(
                [index, info_map(gatk_row), info_map(native_row)])
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the VariantFiltration AS_FilterStatus "
                    "encoding against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_VARIANT_FILTRATION_BINARY"),
        help="native VariantFiltration binary (default: "
             "$FASTGATK_VARIANT_FILTRATION_BINARY or "
             "fastgatk-native/build/fastgatk-variant-filtration)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the measured GATK/native rows and exit 0 "
             "instead of gating parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-variant-filtration")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_variant_filtration_asfilterstatus_gatk_oracle.py", java, jar)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    # A plain (unindexed) VCF keeps both tools on byte-identical input: GATK
    # refuses block-compressed input without an index, native does not.
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-variant-filtration-asfilterstatus-oracle-") as directory:
        work = pathlib.Path(directory)
        for case in selected:
            results.append(run_case(case, work, java, jar, native, arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: VariantFiltration.splitMultiAllelics "
          "(VariantFiltration.java:371-378) evaluates each expression per ALT "
          "against a context built without INFO attributes or genotypes, so an "
          "AS_* INFO rule sees JEXL null and never fires; "
          "AlleleFilterUtils.addAlleleAndSiteFilters (:94-122) then keeps the "
          "'SITE' placeholder (GATKVCFConstants.java:201), joins alleles with '|' "
          "(AnnotationUtils.encodeAnyASListWithRawDelim) and derives the site "
          "FILTER from the intersection of the allele filter sets (:115-120).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args'])}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    expected (FILTER, AS_FilterStatus): {result['expect']}")
        print(f"    GATK   (FILTER, AS_FilterStatus): {result['gatk_observed']}")
        print(f"    NATIVE (FILTER, AS_FilterStatus): {result['native_observed']}")
        for index, gatk_info, native_info in result.get("reported_info_differences", []):
            print(f"    row {index}: reported, not gated, INFO differ: "
                  f"GATK={gatk_info} NATIVE={native_info}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "gatk_version": "4.6.2.0",
        "mode": "strict" if strict_mode else "expect-divergence",
        "cases": [
            {key: value for key, value in result.items() if key != "gated"}
            for result in results
        ],
        "violations": violations,
        "headers_and_index_bytes_compared": False,
    }
    print(json.dumps(payload, sort_keys=True))

    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
