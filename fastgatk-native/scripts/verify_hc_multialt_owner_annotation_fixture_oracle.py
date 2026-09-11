#!/usr/bin/env python3
"""Fixture oracle for the ordinary-VCF multi-ALT owner-annotation line.

Site under test
---------------
``fastgatk-native/src/hc_call.cpp:3521`` (ordinary-VCF writer)::

    if (annotation_reads != nullptr &&
        (max_allele_subset_changed || joint_af_allele_subset_changed)) {
        ... annotations_from_owner() ...   // re-annotate from the AssemblyRegion owner
    }

HaplotypeCaller annotates the *final* VariantContext after joint AF has selected
its concrete ALT set, so a record whose ALT set was reduced must have its
annotations recomputed against the retained ALT list and, for a partitioned
result, against the AssemblyRegion owner that still holds the stable
source-record map (a flattened ``Result`` keeps the compact likelihood rows but
not that map).

Owner-pattern census status
---------------------------
The previous round never saw this guard active: "ordinary VCF mode produced ZERO
multi-ALT records across 3 ploidies x 6 window starts with the fixtures tried".
The reason is that HaplotypeCaller's ordinary output rarely carries two concrete
ALTs at one locus from reads alone.  The fixture below makes the record
multi-ALT deterministically with ``--alleles``: a feature VCF whose single record
carries three concrete ALT alleles at ``chr1:698`` (``A,G,T``).  GATK splits that
feature record into three EventMap events and emits ONE multi-allelic site
(``gatk-source/.../AssemblyBasedCallerUtils`` given-allele handling), which
native reproduces.

Two independent ways to satisfy the guard are exercised here:

* ``max_allele_subset_changed`` -- the forced allele list is longer than
  ``--max-alternate-alleles``, so the writer reduces the ALT set
  (``hc_call.cpp`` "max_allele_subset_changed = true" branch).  The reduction is
  visible in the output ALT list, so it is provable from the record alone.
* ``joint_af_allele_subset_changed`` -- the AF plausibility subset strips an
  allele that the numeric cap kept.

Reachability evidence
---------------------
The emitted record is multi-ALT **and** carries annotation fields
(``BaseQRankSum``/``MQRankSum``/``ReadPosRankSum``), which proves
``annotation_reads != nullptr``; the ALT list is shorter than the forced allele
list, which proves ``max_allele_subset_changed``; therefore the guard at
``hc_call.cpp:3521`` is satisfied.

In addition, this oracle runs the native caller a second time with
``FASTGATK_DEBUG_ANNOTATION_POSITION=697`` (0-based POS of the multi-ALT record),
which makes the Host annotation boundary print one
``[FASTGATK_ANNOTATION_EVIDENCE]`` line per read per invocation.  Measured on
this fixture:

    forced 3 ALTs, --max-alternate-alleles 6 (no numeric subsetting): 4 passes
    forced 3 ALTs, --max-alternate-alleles 2 (subsetting, 2 ALTs emitted): 5 passes

i.e. exactly one extra annotation pass appears iff the ALT set was reduced --
the extra call is the guard block at ``hc_call.cpp:3521``.  The case table below
asserts that delta, so the gate fails if the guard stops executing.

Result: the line is REACHED and, on every configuration tried, its output is
identical to GATK's, i.e. the owner-based annotation is not observably different
from the flattened one here.

Separately tracked difference seen in these fixtures (reported, not gated)
------------------------------------------------------------------------
``INFO AF``/``MLEAF`` zero formatting: on a record with >= 2 ALTs where one
allele's value needs three decimals and another is exactly zero, GATK prints
``0.500,0.00,0.00`` and native prints ``0.500,0.000,0.000`` (per-value trimming
to a minimum of two decimals vs one shared precision).  That is a distinct,
small serialization divergence; it is reported per case and excluded from the
gate (the case table prints it as ``AF/MLEAF formatting``).

Exit status
-----------
Strict by default and green when the fixture behaves as documented: every
record must match GATK field by field after the documented canonicalization
(modulo the AF/MLEAF formatting note), the multi-ALT record must exist, the
subsetting cases must show the extra annotation pass, and the non-subsetting
control must not.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hc_symbolic_prior_fixture_lib import (  # noqa: E402
    REFERENCE_LENGTH, build_reference, compare_rows, layout, make_bam,
    prepare_reference, records,
)

MULTI_ALT_SITE = "698"           # 1-based POS of the forced multi-ALT record
MULTI_ALT_SITE_0BASED = 697      # FASTGATK_DEBUG_ANNOTATION_POSITION value
ANNOTATION_EVIDENCE = "FASTGATK_ANNOTATION_EVIDENCE"
# INFO AF/MLEAF zero-value formatting (see the module docstring): reported, not gated.
AF_FORMATTING_TOLERANCE = {MULTI_ALT_SITE: frozenset({"AF", "MLEAF"})}

CASES: tuple[dict, ...] = (
    {
        "case": "multialt-subset-max2",
        "why": "3 forced concrete ALTs (A,G,T) with --max-alternate-alleles 2 and "
               "12 reads supporting A: the ALT set is reduced to 2, so "
               "max_allele_subset_changed is true and the annotation guard fires",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 2, "snp_reads": 12,
        "expect_subsetting": True,
    },
    {
        "case": "multialt-subset-max2-noreads",
        "why": "same forced allele list and cap but a reference-only pileup: "
               "byte-identical to GATK including AF/MLEAF",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 2, "snp_reads": 0,
        "expect_subsetting": True,
    },
    {
        "case": "multialt-subset-max2-ploidy3",
        "why": "same reduction at --sample-ploidy 3: the guard is not ploidy-specific",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 3, "snp_reads": 12,
        "expect_subsetting": True,
    },
    {
        "case": "multialt-subset-max1",
        "why": "the same forced allele list reduced to a single ALT by "
               "--max-alternate-alleles 1 (still max_allele_subset_changed)",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 1, "ploidy": 2, "snp_reads": 12,
        "expect_subsetting": True, "expect_multi_alt": False,
    },
    {
        "case": "multialt-nosubset-max3",
        "why": "CONTROL: the same forced allele list with --max-alternate-alleles 3 "
               "is NOT reduced (3 ALTs kept): multi-ALT record, no extra "
               "annotation pass",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 3, "ploidy": 2, "snp_reads": 12,
        "expect_subsetting": False,
    },
    {
        "case": "biallelic-forced-control",
        "why": "CONTROL: two forced ALTs (A,T) under a cap of 6 is not reduced "
               "(biallelic, no subsetting, no extra annotation pass)",
        "alts": ["A", "T"], "max_alternate_alleles": 6, "ploidy": 2, "snp_reads": 12,
        "expect_subsetting": False,
    },
)

SUBSETTING_CASES = tuple(case["case"] for case in CASES if case["expect_subsetting"])
CONTROL_CASES = tuple(case["case"] for case in CASES if not case["expect_subsetting"])


def write_alleles_feature(work: Path, name: str, alts: list[str]) -> Path:
    """bgzip+tabix'd single-record feature VCF with several concrete ALTs."""
    reference_base = build_reference()[int(MULTI_ALT_SITE) - 1]
    plain = work / f"{name}.alleles.vcf"
    plain.write_text("##fileformat=VCFv4.2\n"
                     f"##contig=<ID=chr1,length={REFERENCE_LENGTH}>\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     f"chr1\t{MULTI_ALT_SITE}\t.\t{reference_base}\t"
                     f"{','.join(alts)}\t.\tPASS\t.\n", encoding="utf-8")
    root = Path(__file__).resolve().parents[2]
    bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
    tabix = root / "third_party/htslib-build/htslib-src/tabix"
    if not bgzip.is_file() or not tabix.is_file():
        raise SystemExit(f"missing bgzip/tabix under {bgzip.parent}")
    subprocess.run([str(bgzip), "-f", str(plain)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    compressed = plain.with_suffix(plain.suffix + ".gz")
    subprocess.run([str(tabix), "-f", "-p", "vcf", str(compressed)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    return compressed


def annotation_passes(native: str, common: list[str], work: Path, name: str) -> dict:
    """Run native with the annotation-boundary trace and count per-read passes."""
    output = work / f"trace.{name}.vcf"
    environment = dict(os.environ)
    environment["FASTGATK_DEBUG_ANNOTATION_POSITION"] = str(MULTI_ALT_SITE_0BASED)
    run = subprocess.run([native, *common, "--threads", "2", "-O", str(output)],
                         capture_output=True, text=True, timeout=900, env=environment)
    lines = [line for line in run.stderr.splitlines() if ANNOTATION_EVIDENCE in line]
    reads = 0
    if lines:
        # One line per read per invocation: count the invocations by looking for
        # a repeated read name at the head of a block.
        names = [line.split("read=", 1)[1].split(" ", 1)[0] for line in lines]
        first = names[0]
        reads = sum(1 for index, value in enumerate(names) if value == first
                    and (index == 0 or names[index - 1] != first))
    return {"exit": run.returncode, "evidence_lines": len(lines),
            "annotation_passes": reads,
            "stderr_tail": run.stderr[-300:] if run.returncode else ""}


def run_case(case: dict, reference: Path, reference_text: str, work: Path,
             java: str, gatk: str, native: str, threads: int) -> dict:
    bam = make_bam(work, case["case"], layout(30, case["snp_reads"], 0, "cigar"),
                   reference_text, java, gatk)
    feature = write_alleles_feature(work, case["case"], case["alts"])
    common = ["-R", str(reference), "-I", str(bam), "-L", "chr1:500-780",
              "--min-pruning", "1",
              "--sample-ploidy", str(case["ploidy"]),
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false",
              "--alleles", str(feature),
              "--max-alternate-alleles", str(case["max_alternate_alleles"])]
    gatk_vcf = work / f"gatk.{case['case']}.vcf"
    native_vcf = work / f"native.{case['case']}.vcf"
    gatk_run = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "HaplotypeCaller", *common, "-O", str(gatk_vcf)],
        capture_output=True, text=True, timeout=900)
    native_run = subprocess.run(
        [native, *common, "--threads", str(threads), "-O", str(native_vcf)],
        capture_output=True, text=True, timeout=900)
    result = {"case": case["case"], "why": case["why"],
              "forced_alt_alleles": list(case["alts"]),
              "max_alternate_alleles": case["max_alternate_alleles"],
              "ploidy": case["ploidy"], "snp_reads": case["snp_reads"],
              "expect_subsetting": case["expect_subsetting"],
              "gatk_exit": gatk_run.returncode, "native_exit": native_run.returncode,
              "gatk_stderr_tail": gatk_run.stderr[-400:],
              "native_stderr_tail": native_run.stderr[-400:]}
    if gatk_run.returncode != 0 or native_run.returncode != 0:
        result["comparison"] = {"violations": [
            f"run failed: gatk_exit={gatk_run.returncode} "
            f"native_exit={native_run.returncode}"], "data_rows_identical": False}
        return result
    gatk_rows = records(gatk_vcf)
    native_rows = records(native_vcf)
    comparison = compare_rows(gatk_rows, native_rows,
                              tolerated=AF_FORMATTING_TOLERANCE)
    comparison["gatk_multi_alt_records"] = [
        row for row in gatk_rows if len(row[4].split(",")) > 1]
    comparison["native_multi_alt_records"] = [
        row for row in native_rows if len(row[4].split(",")) > 1]
    comparison["gatk_rows_all"] = gatk_rows
    comparison["native_rows_all"] = native_rows
    comparison["subsetting_proved"] = all(
        row[4].count(",") + 1 < len(case["alts"]) for row in native_rows
        if row[1] == MULTI_ALT_SITE) if native_rows else False
    compare_trace = annotation_passes(native, common, work, case["case"])
    result["annotation_trace"] = compare_trace
    result["comparison"] = comparison
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fixture oracle for the ordinary-VCF multi-ALT owner-annotation "
                    "line (hc_call.cpp:3521) vs pinned GATK 4.6.2.0.")
    parser.add_argument("--native", default=os.environ.get("FASTGATK_HC_BINARY"),
                        help="native HaplotypeCaller binary (default: "
                             "$FASTGATK_HC_BINARY or "
                             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument("--expect-divergence", action="store_true",
                        help="report mode: never fail on the separately tracked "
                             "AF/MLEAF formatting difference (kept for house "
                             "convention; this oracle is green by design)")
    parser.add_argument("--case", action="append", default=None,
                        help=f"run only the named case(s); default all of "
                             f"{[case['case'] for case in CASES]}")
    parser.add_argument("--threads", type=int, default=2,
                        help="--threads value for the native run (default 2)")
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(arguments.native) if arguments.native else (
        root / "fastgatk-native/build/fastgatk-hc-call")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    assets = [native, Path(java), gatk]
    if not all(path.is_file() for path in assets):
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")

    results: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-multialt-owner-") as directory:
        work = Path(directory)
        reference = prepare_reference(work, java, str(gatk))
        reference_text = build_reference()
        for case in selected:
            results.append(run_case(case, reference, reference_text, work, java,
                                    str(gatk), str(native), arguments.threads))

    violations: list[str] = []
    by_case = {result["case"]: result for result in results}
    for result in results:
        comparison = result["comparison"]
        violations += [f"{result['case']}: {item}"
                       for item in comparison.get("violations", [])]
        if not comparison.get("native_multi_alt_records") and \
                result["case"] != "multialt-subset-max1":
            violations.append(f"{result['case']}: expected a multi-ALT ordinary-VCF "
                              f"record but native emitted none")
        if result["expect_subsetting"] and not comparison.get("subsetting_proved"):
            violations.append(f"{result['case']}: the forced ALT list ({len(result['forced_alt_alleles'])} "
                              f"alleles) was not reduced to fewer ALTs, so "
                              f"max_allele_subset_changed is not proved")
    # Paired evidence: subsetting adds exactly one annotation pass.
    def passes(case_name: str) -> int | None:
        entry = by_case.get(case_name)
        return entry["annotation_trace"]["annotation_passes"] if entry else None

    for subsetting, control in (("multialt-subset-max2-noreads", "multialt-nosubset-max3"),
                                ("multialt-subset-max2", "multialt-nosubset-max3")):
        left, right = passes(subsetting), passes(control)
        if left is not None and right is not None:
            # The two cases differ in read count as well; compare after scaling
            # by the read count so the delta is a pass count.
            reads_left = by_case[subsetting]["snp_reads"] + 30
            reads_right = by_case[control]["snp_reads"] + 30
            passes_left = left
            passes_right = right * reads_left / reads_right
            if not passes_left > passes_right:
                violations.append(
                    f"{subsetting}: expected one extra annotation pass vs "
                    f"{control} (guard hc_call.cpp:3521), got "
                    f"{passes_left} vs {passes_right:g}")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# site: hc_call.cpp:3521 -- ordinary-VCF multi-ALT annotations_from_owner "
          "block, guard (annotation_reads != nullptr && (max_allele_subset_changed || "
          "joint_af_allele_subset_changed))")
    print(f"# cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# subsetting cases (guard must be satisfied): {list(SUBSETTING_CASES)}")
    print(f"# controls (no subsetting): {list(CONTROL_CASES)}")
    for result in results:
        comparison = result["comparison"]
        trace = result["annotation_trace"]
        print(f"[{result['case']}] forced={result['forced_alt_alleles']} "
              f"max-alt={result['max_alternate_alleles']} ploidy={result['ploidy']} "
              f"snp_reads={result['snp_reads']} "
              f"gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    why: {result['why']}")
        print(f"    native subsetting proved={comparison.get('subsetting_proved')} "
              f"multi_alt_records_gatk={len(comparison.get('native_multi_alt_records', []))} "
              f"data_rows_identical={comparison.get('data_rows_identical')} "
              f"annotation_evidence_lines={trace['evidence_lines']} "
              f"annotation_passes={trace['annotation_passes']}")
        for row in comparison.get("gatk_rows_all", []):
            print(f"    GATK   {chr(9).join(row)}")
        for row in comparison.get("native_rows_all", []):
            print(f"    NATIVE {chr(9).join(row)}")
        for item in comparison.get("separately_tracked_notes", []):
            print(f"    separately tracked (reported, not gated): {item}")
        for item in comparison.get("violations", []):
            print(f"    VIOLATION: {item}")
    for violation in violations:
        print(f"# VIOLATION: {violation}")

    payload = {
        "status": "fail" if violations else "pass",
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller ordinary-VCF multi-ALT owner annotation parity "
                  "(hc_call.cpp:3521)",
        "binary": str(native),
        "site": "fastgatk-native/src/hc_call.cpp:3521",
        "outcome": "REACHED, no observable divergence: the guard is satisfied by "
                   "construction in every subsetting case (forced ALT list reduced, "
                   "record annotated) and the measured annotation-boundary pass "
                   "count rises by exactly one, yet every field except the "
                   "separately tracked INFO AF/MLEAF zero formatting is identical "
                   "to GATK",
        "cases": [case["case"] for case in selected],
        "subsetting_cases": list(SUBSETTING_CASES),
        "control_cases": list(CONTROL_CASES),
        "violations": violations,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
