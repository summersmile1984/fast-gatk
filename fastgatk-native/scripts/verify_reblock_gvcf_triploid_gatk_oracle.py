#!/usr/bin/env python3
"""Strict oracle for ReblockGVCF polyploid reblocking driven by ``PL``
(pinned GATK 4.6.2.0).

The gap this gate pins
----------------------
``fastgatk-native/scripts/verify_reblock_gvcf.py`` reblocks a triploid record::

    chr1 30 . A C,G,<NON_REF> . PASS DP=17 \
        GT:DP:AD:PL:GQ 0/1/1:17:12,5,0,0:0,10,20,30,...,190:99

and asserted that it stays a variant -- ``GT 0/1/1``, ALT ``C,<NON_REF>`` and a
PL array compacted to the 10 triploid entries of the three remaining alleles.
Pinned GATK instead emits a **reference block**: ``GT 0/0/0``, ALT
``<NON_REF>``, ``END=30`` and 4 PLs.  The native-only variant encoding is the
divergence this gate pins.

GATK's rule (read from ``gatk-source/``, then confirmed by measurement)
----------------------------------------------------------------------
``ReblockGVCF.shouldBeReblocked`` (``ReblockGVCF.java:514-535``) does **not**
look at the called GT.  It takes the likelihood vector, finds the *minimum*
entry (``:524`` ``MathUtils.minElementIndex(pls)``), decodes that index into a
genotype at the sample's ploidy (``:525``
``GenotypesCache.get(genotype.getPloidy(), minLikelihoodIndex)``), maps it onto
the record's alleles (``:527`` ``alleleCounts.asAlleleList(...)``) and reblocks
whenever

* ``pls[0] < rgqThreshold`` (``:528``), or
* that genotype contains no concrete ALT allele (``:529``
  ``!GATKVariantContextUtils.genotypeHasConcreteAlt(finalAlleles)``), or
* it contains ``<NON_REF>`` (``:530``), or
* the genotype has neither PL nor GQ (``:531``), or
* ``TREE_SCORE`` is below the threshold (``:532``).

For the fixture above the minimum PL is ``PL[0] = 0``, whose triploid index 0
decodes to ``A/A/A`` -- reference only -- so clause ``:529`` fires and
``lowQualVariantToGQ0HomRef`` (``:537-563``) rewrites the call:
``changeCallToHomRefVersusNonRef`` (``:576-631``) subsets the likelihoods to the
reference plus the single most likely ALT
(``AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, ploidy, 1, true)``,
``:607-612``), so the output PL has ``numLikelihoods(2, ploidy)`` entries, the
genotype becomes the reference repeated at the sample ploidy (``:629``) and
``END`` is set to the record's end (``:627``).

Native classified "variant versus reference block" from the record's ALT list
alone, so a concrete-ALT record whose likeliest call is hom-ref was never
reblocked.

Scope and comparison contract
-----------------------------
Each gated case runs pinned GATK and native with identical arguments on the same
plain (unindexed) VCF and asserts the **data rows are byte-identical**.  The
first case also carries an untouched reference block as a surviving control, so
the gate cannot be satisfied by deleting records.

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

SURVIVING_BLOCK = "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:GQ\t0/0:10:20\n"

# The literal triploid record of verify_reblock_gvcf.py.  PL is the complete
# Number=G vector for 4 alleles at ploidy 3: C(4+3-1, 3) = 20 entries.
TRIPLOID = (
    "chr1\t30\t.\tA\tC,G,<NON_REF>\t.\t.\tDP=17\t"
    "GT:DP:AD:PL:GQ\t0/1/1:17:12,5,0,0:"
    "0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190:99\n"
)

# GATK's reblocked encoding: the reference repeated at the sample ploidy, the
# likelihoods subset to REF plus the most likely ALT (4 = numLikelihoods(2, 3)),
# GQ recomputed from that vector as 10 - 0, and END retained from the record.
TRIPLOID_REBLOCKED = (
    "chr1\t30\t.\tA\t<NON_REF>\t.\t.\tEND=30\t"
    "GT:DP:GQ:MIN_DP:PL\t0/0/0:17:10:17:0,10,20,30"
)

CASES = [
    {
        "case": "triploid-with-surviving-block",
        "why": "shouldBeReblocked uses the minimum-likelihood PL genotype "
               "(ReblockGVCF.java:524-527), not the called GT: PL[0]=0 decodes to "
               "0/0/0 at ploidy 3, so genotypeHasConcreteAlt is false (:529) and "
               "the record is reblocked instead of kept as a variant; the "
               "preceding reference block is the control row",
        "body": SURVIVING_BLOCK + TRIPLOID,
        "args": [],
        "gated": True,
        "expect": [
            "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:GQ:MIN_DP\t0/0:10:20:10",
            TRIPLOID_REBLOCKED,
        ],
    },
    {
        "case": "triploid-alone",
        "why": "the literal verify_reblock_gvcf.py fixture: a single triploid "
               "record becomes one reference block with END=30 and 4 PLs",
        "body": TRIPLOID,
        "args": [],
        "gated": True,
        "expect": [TRIPLOID_REBLOCKED],
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
    """A 100 bp chr1 whose bases match every fixture REF allele (A)."""
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + ("A" * 100) + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    # GATK resolves the sequence dictionary as <stem>.dict.
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
    if gatk_rows != case["expect"]:
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
        description="Strict oracle for ReblockGVCF polyploid reblocking driven "
                    "by the minimum-likelihood PL genotype, against pinned GATK "
                    "4.6.2.0.")
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
            "verify_reblock_gvcf_triploid_gatk_oracle.py", java, jar)
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
            prefix="fastgatk-reblock-triploid-oracle-") as directory:
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
    print("# GATK rule under test: shouldBeReblocked (ReblockGVCF.java:514-535) "
          "decodes the minimum-likelihood PL index into a genotype at the "
          "sample ploidy (:524-527) and reblocks when that genotype has no "
          "concrete ALT (:529) or contains <NON_REF> (:530); "
          "lowQualVariantToGQ0HomRef/changeCallToHomRefVersusNonRef (:537-631) "
          "then subsets PLs to REF plus the most likely ALT "
          "(AlleleSubsettingUtils.calculateMostLikelyAlleles), sets GT to the "
          "ploidy-fold reference and keeps END.")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args'])}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
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
