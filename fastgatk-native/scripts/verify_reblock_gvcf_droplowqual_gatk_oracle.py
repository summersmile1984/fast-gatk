#!/usr/bin/env python3
"""Strict oracle for ReblockGVCF ``--drop-low-quals`` + ``--rgq-threshold``
(``--drop-low-quals`` wins; pinned GATK 4.6.2.0).

The gap this gate pins
----------------------
``fastgatk-native/scripts/verify_reblock_gvcf.py`` runs ReblockGVCF with
``-GQB 20 -GQB 100 --drop-low-quals --rgq-threshold 10 --floor-blocks`` over a
GVCF whose fourth record is a low-quality concrete call::

    chr1 20 . C T,<NON_REF> . PASS DP=5 GT:DP:AD:PL:GQ 0/1:5:4,1,0:2,0,40,99,99,99:99

and asserted that the record survives as a third, GQ0 hom-ref **reference
block** at POS 20 (``END=20``).  Pinned GATK emits **two** rows for the whole
file and drops that site entirely; the native-only third row is the divergence
this gate pins.

GATK's rule (read from ``gatk-source/``, then confirmed by measurement)
----------------------------------------------------------------------
1. ``ReblockGVCF.regenotypeVC`` handles a concrete variant in two ordered
   steps.  First, with ``--drop-low-quals`` set and a non-zero site depth,
   GATK re-genotypes the site with ``standardConfidenceForCalling``
   (``ReblockGVCF.java:415-424``; the confidence is armed only in drop mode at
   ``:327``).  A site the genotyping engine cannot call is **dropped** -- the
   method returns without writing anything.
2. Only afterwards does ``shouldBeReblocked`` (``ReblockGVCF.java:427``,
   defined at ``:514-535``) run, whose ``pls[0] < rgqThreshold`` clause
   (``:528``) is the ``--rgq-threshold`` conversion.  So ``--rgq-threshold``
   can never resurrect a site that drop-mode genotyping already removed: the
   PL[0] rule is a *later* filter, not an earlier one.
3. ``lowQualVariantToGQ0HomRef`` (``:537-563``) additionally returns ``null``
   -- dropping the record -- when drop mode is on and the variant is not a
   monomorphic hom-ref call with concrete ALTs (``:542-545``).

Native ran the ``--rgq-threshold`` conversion *before* its drop-mode
re-genotyping emulation, so the low-quality site at POS 20 was converted into a
GQ0 reference block instead of being dropped.

Scope and comparison contract
-----------------------------
Each gated case runs pinned GATK and native with identical arguments on the
same plain (unindexed) VCF and asserts the **data rows are byte-identical**.
Reference-block rows are used as the surviving control because both
implementations encode them identically (verified by measurement), so a
whole-row byte comparison is meaningful and a record-count difference -- the
pinned defect -- cannot hide behind a projection.

The exact fixture from ``verify_reblock_gvcf.py`` (four records, including the
multiallelic site whose QUAL/annotation FORMAT ordering native has not yet been
reconciled with GATK) is also run and **reported only**: its surviving variant
row differs in QUAL ("." vs the re-genotyped 42.64) and in INFO/FORMAT ordering,
which are unrelated open divergences and not what this gate is about.

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

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
    "##INFO=<ID=TREE_SCORE,Number=1,Type=Float,Description=Tree score>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
    "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum DP within a block>\n"
)
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

# The reference block is written in the FORMAT subset (GT:DP:GQ) that GATK and
# native encode byte-identically, so the surviving row is a real control.
SURVIVING_BLOCK = "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:GQ\t0/0:10:20\n"
LOW_QUAL_SITE = (
    "chr1\t20\t.\tC\tT,<NON_REF>\t.\t.\tDP=5\t"
    "GT:DP:AD:PL:GQ\t0/1:5:4,1,0:2,0,40,99,99,99:99\n"
)

# The literal fixture of verify_reblock_gvcf.py (reported only, see docstring).
AUDIT_FIXTURE = (
    "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT:DP:AD:PL:GQ\t0/0:10:10,0:0,10,200:10\n"
    "chr1\t6\t.\tC\t<NON_REF>\t.\tPASS\tEND=10\tGT:DP:AD:PL:GQ\t0/0:8:8,0:0,15,200:15\n"
    "chr1\t12\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL:GQ\t0/1:20:12,8,0,0:50,0,80,99,99,99,99,99,99,20:.\n"
    "chr1\t20\t.\tC\tT,<NON_REF>\t.\tPASS\tDP=5\t"
    "GT:DP:AD:PL:GQ\t0/1:5:4,1,0:2,0,40,99,99,99:99\n"
)

# The exact option vector of the main case of verify_reblock_gvcf.py.
DROP_FLAGS = ["-GQB", "20", "-GQB", "100", "--drop-low-quals",
              "--rgq-threshold", "10", "--floor-blocks"]

CASES = [
    {
        "case": "drop-low-quals-with-rgq-threshold",
        "why": "ReblockGVCF.java:415 re-genotypes (and drops) a concrete variant "
               "under --drop-low-quals before shouldBeReblocked (:427) can apply "
               "the PL[0] < rgqThreshold clause (:528), so --rgq-threshold must "
               "not turn a site that drop mode removes into a GQ0 block",
        "fixture": "reduced",
        "body": SURVIVING_BLOCK + LOW_QUAL_SITE,
        "args": DROP_FLAGS,
        "gated": True,
        # Only the untouched reference block survives; the POS 20 site is gone.
        # --floor-blocks drops MIN_DP/PL from the block, so the surviving row
        # carries GT:DP:GQ only (measured).
        "expect": ["chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:GQ\t0/0:10:20"],
    },
    {
        "case": "drop-low-quals-alone",
        "why": "control: with drop mode but no --rgq-threshold native already "
               "matched GATK here, so the gate is not vacuous and the surviving "
               "block row is genuinely byte-identical",
        "fixture": "reduced",
        "body": SURVIVING_BLOCK + LOW_QUAL_SITE,
        "args": ["--drop-low-quals"],
        "gated": True,
        "expect": ["chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:GQ:MIN_DP\t0/0:10:20:10"],
    },
    {
        "case": "audit-fixture-diagnostic",
        "why": "reported only: the literal verify_reblock_gvcf.py fixture.  Its "
               "surviving multiallelic row still differs from GATK in QUAL (native "
               "leaves '.', GATK writes the re-genotyped value) and in INFO/FORMAT "
               "ordering, which are unrelated open divergences; the record set "
               "(2 rows, no POS 20 row) is the part this gate pins",
        "fixture": "audit",
        "body": AUDIT_FIXTURE,
        "args": DROP_FLAGS,
        "gated": False,
        "expect": None,
    },
]


def data_rows(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows: list[str] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                rows.append(line.rstrip("\n"))
    return rows


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
    bases = ["A"] * 100
    bases[5] = "C"   # POS 6  REF=C
    bases[19] = "C"  # POS 20 REF=C
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "".join(bases) + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    # GATK resolves the sequence dictionary as <stem>.dict (measured: it
    # reports "Fasta dict file ... reference.dict ... does not exist" for a
    # "reference.fa.dict" sibling).
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-4000:])
    return result


def run_case(case: dict, work: pathlib.Path, reference: pathlib.Path,
             java: pathlib.Path, jar: pathlib.Path, native: pathlib.Path,
             timeout: int) -> dict:
    source = work / f"{case['case']}.vcf"
    source.write_text(HEADER + COLUMNS + case["body"], encoding="utf-8")
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    # Identical arguments on both sides; only the output path differs.  A plain
    # (unindexed) VCF keeps the input bytes identical: GATK refuses
    # block-compressed input without an index, native does not.
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar), "ReblockGVCF",
                          *common, "-O", str(gatk_out)],
                         f"GATK ReblockGVCF [{case['case']}]", timeout)
    native_result = invoke([str(native), *common, "-O", str(native_out)],
                           f"native ReblockGVCF [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "violations": [],
    }
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    if case["expect"] is not None and gatk_rows != case["expect"]:
        result["violations"].append(
            "GATK moved away from the measured truth: "
            f"GATK={gatk_rows} expected={case['expect']}")
    if case["gated"]:
        if len(gatk_rows) != len(native_rows):
            result["violations"].append(
                f"record count differs: GATK={len(gatk_rows)} native={len(native_rows)}")
        for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows)):
            if gatk_row != native_row:
                result["violations"].append(
                    f"row {index} is not byte-identical: GATK={gatk_row!r} "
                    f"NATIVE={native_row!r}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for ReblockGVCF --drop-low-quals with "
                    "--rgq-threshold, against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_REBLOCK_BINARY"),
        help="native ReblockGVCF binary (default: $FASTGATK_REBLOCK_BINARY or "
             "fastgatk-native/build/fastgatk-reblock-gvcf)")
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
        / "fastgatk-reblock-gvcf")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_reblock_gvcf_droplowqual_gatk_oracle.py", java, jar)
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
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-reblock-droplowqual-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: regenotypeVC re-genotypes concrete variants "
          "under --drop-low-quals first (ReblockGVCF.java:415-424, confidence "
          "armed at :327) and returns without writing when the site cannot be "
          "called; only then does shouldBeReblocked (:427, defined :514-535) "
          "apply the PL[0] < rgqThreshold clause (:528) and "
          "lowQualVariantToGQ0HomRef nulls out non-monomorphic calls (:542-545).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args'])}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        if result["expect"] is not None:
            print(f"    expected rows: {result['expect']}")
        print(f"    GATK   rows ({len(result['gatk_rows'])}):")
        for row in result["gatk_rows"]:
            print(f"        {row}")
        print(f"    NATIVE rows ({len(result['native_rows'])}):")
        for row in result["native_rows"]:
            print(f"        {row}")
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
        "rows_compared": "data rows byte-identical (CHROM..sample columns)",
    }
    print(json.dumps(payload, sort_keys=True))

    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
