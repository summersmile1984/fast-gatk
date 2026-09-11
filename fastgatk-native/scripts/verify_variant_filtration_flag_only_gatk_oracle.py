#!/usr/bin/env python3
"""Strict oracle for VariantFiltration ``--apply-allele-specific-filters`` when
**no** ``AS_*`` expression is in play (pinned GATK 4.6.2.0).

The adjacent gap this gate pins
-------------------------------
The previous round fixed the ``AS_FilterStatus`` encoding for the case
"flag set *and* an ``AS_*`` rule present" and deliberately reported, but did not
fix, the adjacent case "flag set, no ``AS_*`` rule".  Re-measuring that case
against the pinned jar shows the divergence is **not** confined to a missing
``AS_FilterStatus``: with the flag set GATK does not evaluate site rules against
the record at all.

GATK's rule (read from ``gatk-source/``, then confirmed by measurement)
----------------------------------------------------------------------
1. ``VariantFiltration.apply`` (``VariantFiltration.java:359-365``) takes the
   allele path whenever ``--apply-allele-specific-filters`` is set -- *not* only
   when an ``AS_*`` expression is present.  Every site expression, the mask and
   the cluster test are evaluated once per split ALT.
2. ``splitMultiAllelics`` (``:371-378``) rebuilds each ALT as a fresh
   ``VariantContext`` from
   ``new VariantContextBuilder("SimpleSplit", contig, start, end, [ref, NO_CALL])``.
   That 5-argument htsjdk constructor carries **only** contig/start/stop/alleles:
   the split context has no INFO attributes and no genotypes, and its
   ``log10PError`` keeps the builder default ``VariantContext.NO_LOG10_PERROR``
   = ``1.0d`` (verified with ``javap -constants``), so the JEXL ``QUAL``
   attribute -- ``VariantJEXLContext`` computes it as
   ``-10.0 * vc.getLog10PError()`` -- is the constant ``-10.0`` for every split
   context, whatever the record's real QUAL is.
3. Consequences, all measured below:
   * an INFO rule such as ``DP > 5`` never fires (the attribute is unavailable);
   * a QUAL rule such as ``QUAL < 60`` fires for **every** ALT of **every**
     record, including records whose real QUAL is above the threshold;
   * ``--mask`` is applied per split ALT, so a masked record's ALTs all carry
     the mask name.
4. ``AlleleFilterUtils.addAlleleAndSiteFilters`` (``:94-122``) then merges the
   per-ALT filter sets onto the decoded ``AS_FilterStatus`` (``SITE`` placeholder
   for an ALT no filter selected, ``GATKVCFConstants.java:201``), joins alleles
   with ``|`` (``AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM``) and sets the site
   FILTER to the **intersection** of the new per-ALT sets (``:115-117``), added
   to the record's existing filters, PASS when the result is empty and
   unfiltered under ``--invalidate-previous-filters`` (``:118-120``).

So ``AS_FilterStatus`` appears whenever the flag is set (``LowQual|LowQual``,
``Mask|Mask``, ``SITE|SITE``), *and* the FILTER column follows the split-context
evaluation rather than the record's own QUAL/INFO.

Scope
-----
Each gated case runs pinned GATK and native with identical arguments and
compares the data rows: ``CHROM..FILTER`` and the whole INFO map (key order is
not part of the contract, and no floating-point INFO is used, so an exact map
comparison is safe here).  The literal GATK expectation measured for each case
is asserted too, so fixture/CLI drift cannot pass silently.

Case ``flag-without-expression`` is *reported only*: GATK accepts
``--apply-allele-specific-filters`` with no expression and no mask, while native
rejects the invocation (``at least one site/genotype filter expression or --mask
is required``).  That is a CLI-validation difference, not an ``AS_FilterStatus``
difference, and it is out of scope for this gate; its measured rows are printed
rather than gated.

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
    "##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
    "##FILTER=<ID=foo,Description=A pre-existing input filter>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
)
MASK_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
)
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

# Record 1: multi-ALT, QUAL 50 (below the threshold).
# Record 2: multi-ALT, QUAL 80 (above the threshold) and already filtered.
# Record 3: biallelic, QUAL 30 (below the threshold), so the per-ALT vector
#           is one element long.
BODY = (
    "chr1\t1\t.\tA\tC,G\t50\tPASS\tDP=10\tGT\t0/1\n"
    "chr1\t2\t.\tC\tT,G\t80\tfoo\tDP=10\tGT\t0/1\n"
    "chr1\t3\t.\tG\tA\t30\tPASS\tDP=10\tGT\t1/1\n"
)

FLAG = "--apply-allele-specific-filters"
QUAL_RULE = ["--filter-expression", "QUAL < 60", "--filter-name", "LowQual"]
INFO_RULE = ["--filter-expression", "DP > 5", "--filter-name", "HighDP"]

# Record 1 overlaps the mask, records 2 and 3 do not.
MASK_BODY = "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n"

CASES = [
    {
        "case": "flag-qual-expression",
        "why": "VariantFiltration.java:359-365 always takes the allele path when "
               "the flag is set; the split context's QUAL is the constant -10.0 "
               "(VariantContext.NO_LOG10_PERROR, VariantJEXLContext), so "
               "'QUAL < 60' fires for every ALT of every record -- including "
               "record 2, whose real QUAL is 80 -- and the site FILTER becomes "
               "the intersection of the per-ALT sets (:115-120)",
        "body": BODY,
        "args": QUAL_RULE + [FLAG],
        "gated": True,
        "expect": [
            ("LowQual", "LowQual|LowQual"),
            ("LowQual;foo", "LowQual|LowQual"),
            ("LowQual", "LowQual"),
        ],
    },
    {
        "case": "flag-info-expression",
        "why": "the split context has no INFO attributes at all "
               "(VariantFiltration.java:371-378), so an INFO rule such as "
               "'DP > 5' never fires even though the input record has DP=10; "
               "AS_FilterStatus is still written with the SITE placeholder",
        "body": BODY,
        "args": INFO_RULE + [FLAG],
        "gated": True,
        "expect": [
            ("PASS", "SITE|SITE"),
            ("foo", "SITE|SITE"),
            ("PASS", "SITE"),
        ],
    },
    {
        "case": "flag-mask",
        "why": "the mask is applied inside filter() to each split ALT "
               "(VariantFiltration.java:379-392, :400-403), so a masked record "
               "gets the mask name in every per-ALT slot and in the site "
               "FILTER, while an unmasked record keeps SITE|SITE",
        "body": BODY,
        "args": ["--mask", "@mask", "--mask-name", "Mask", FLAG],
        "gated": True,
        "expect": [
            ("Mask", "Mask|Mask"),
            ("foo", "SITE|SITE"),
            ("PASS", "SITE"),
        ],
    },
    {
        "case": "no-flag-control",
        "why": "control: without the flag GATK filters the record itself "
               "(VariantFiltration.java:367) and writes no AS_FilterStatus at "
               "all -- this is the behaviour native must keep",
        "body": BODY,
        "args": QUAL_RULE,
        "gated": True,
        "expect": [
            ("LowQual", None),
            ("foo", None),
            ("LowQual", None),
        ],
    },
    {
        "case": "flag-without-expression",
        "why": "reported only: with the flag and no expression and no mask "
               "(=every per-ALT filter set is empty) GATK still writes the "
               "SITE placeholder vector; native rejects the invocation during "
               "argument validation, a separate CLI-contract difference",
        "body": BODY,
        "args": [FLAG],
        "gated": False,
        "expect": [
            ("PASS", "SITE|SITE"),
            ("foo", "SITE|SITE"),
            ("PASS", "SITE"),
        ],
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


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-4000:])
    return result


def run_case(case: dict, work: pathlib.Path, java: pathlib.Path, jar: pathlib.Path,
             native: pathlib.Path, timeout: int) -> dict:
    source = work / f"{case['case']}.vcf"
    source.write_text(BASE_HEADER + COLUMNS + case["body"], encoding="utf-8")
    if any(item == "@mask" for item in case["args"]):
        # GATK's --mask is a FeatureInput and needs a random-access index
        # (measured: "Input mask.vcf must support random access to enable
        # queries by interval"); native does not need one.
        mask = work / "mask.vcf"
        mask.write_text(MASK_HEADER + COLUMNS + MASK_BODY, encoding="utf-8")
        invoke([str(java), "-jar", str(jar), "IndexFeatureFile", "-I", str(mask)],
               "GATK IndexFeatureFile for the VariantFiltration mask", timeout)
    args = [str(mask) if item == "@mask" else item for item in case["args"]]
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf.gz"
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar), "VariantFiltration",
                          "-V", str(source), "-O", str(gatk_out), *args],
                         f"GATK VariantFiltration [{case['case']}]", timeout)
    native_result = invoke([str(native), "-V", str(source), "-O", str(native_out), *args],
                           f"native VariantFiltration [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    gatk_observed = [observed(row) for row in gatk_rows]
    native_observed = [observed(row) for row in native_rows]

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "args": args,
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
            "GATK (FILTER, AS_FilterStatus) moved away from the measured truth: "
            f"GATK={gatk_observed} expected={case['expect']}")
    if case["gated"]:
        if len(gatk_observed) != len(native_observed):
            result["violations"].append(
                f"record count differs: GATK={len(gatk_observed)} native={len(native_observed)}")
        for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows)):
            if gatk_row.split("\t")[:7] != native_row.split("\t")[:7]:
                result["violations"].append(
                    f"row {index}: CHROM..FILTER differ: "
                    f"GATK={gatk_row.split(chr(9))[:7]} NATIVE={native_row.split(chr(9))[:7]}")
            if info_map(gatk_row) != info_map(native_row):
                result["violations"].append(
                    f"row {index}: INFO differ: "
                    f"GATK={info_map(gatk_row)} NATIVE={info_map(native_row)}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for VariantFiltration "
                    "--apply-allele-specific-filters with no AS_* expression, "
                    "against pinned GATK 4.6.2.0.")
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
            "verify_variant_filtration_flag_only_gatk_oracle.py", java, jar)
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
            prefix="fastgatk-variant-filtration-flag-only-oracle-") as directory:
        work = pathlib.Path(directory)
        for case in selected:
            results.append(run_case(case, work, java, jar, native, arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: --apply-allele-specific-filters always takes the "
          "allele path (VariantFiltration.java:359-365); each ALT is filtered as a "
          "fresh context built without INFO or genotypes and with the default "
          "log10PError (NO_LOG10_PERROR = 1.0, so JEXL QUAL == -10.0) "
          "(:371-378); the mask and cluster tests run per ALT (:400-424); "
          "AlleleFilterUtils.addAlleleAndSiteFilters (:94-122) then writes the "
          "per-ALT AS_FilterStatus with the SITE placeholder and sets the site "
          "FILTER to the intersection of the per-ALT filter sets.")
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
        for index, (gatk_row, native_row) in enumerate(zip(result["gatk_rows"],
                                                           result["native_rows"])):
            print(f"    row {index}: GATK={gatk_row}")
            print(f"    row {index}: NATIVE={native_row}")
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
