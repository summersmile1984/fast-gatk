#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the orphan spanning-deletion gVCF fixture used
by ``fastgatk-native/scripts/verify_genotype_gvcf.py``.

The gap this gate pins
----------------------
``verify_genotype_gvcf.py`` builds, at lines 373-381, a single-sample gVCF
record whose ALT list is ``*,G,<NON_REF>`` and asserts at line 391 that native
GenotypeGVCFs emits exactly one record whose ALT contains the symbolic
spanning deletion ``*``.

Pinned GATK 4.6.2.0 emits **zero** records for that fixture by default and, with
``--include-non-variant-sites`` (short name ``-all-sites``, both spellings bound
to the same ``@Argument`` at ``GenotypeGVCFs.java:116-117,126-127``), exactly one
record that carries ``ALT='.'`` and ``GT='./.'``.  The old assertion therefore
pinned native-only output rather than the GATK contract.

GATK's rule (read from ``gatk-source/`` and then confirmed by measurement)
-------------------------------------------------------------------------
1. ``GenotypingEngine.calculateOutputAlleleSubset`` only publishes an ALT that
   is individually plausible at ``--standard-min-confidence-threshold-for-calling``
   (default 30), and it drops a symbolic ``*`` outright when no emitted deletion
   covers the locus (``GenotypingEngine.java:311-327``).  For this fixture both
   ``*`` and ``G`` are pruned, so the surviving ALT set is empty.
2. With an empty ALT set the site stays "monomorphic" and
   ``passesEmitThreshold()`` is false because the best guess is the reference
   (``GenotypingEngine.java:425-427``), so ``calculateGenotypes`` returns
   ``null`` and the default traversal writes nothing
   (``GenotypingEngine.java:167-169``).
3. ``--include-non-variant-sites`` sets ``OutputMode.EMIT_ALL_ACTIVE_SITES``
   (``GenotypeGVCFsEngine.java:381-383``), which makes ``emitAllActiveSites()``
   true (``GenotypingEngine.java:421-423``) and skips that ``null`` return.  The
   rebuilt VariantContext is REF-only: ``subsetToRefOnly()`` supplies the
   genotypes (``GenotypingEngine.java:190``) and
   ``cleanupGenotypeAnnotations(result, true, false)`` then replaces each of
   them with ``./.`` and drops every FORMAT field but GT
   (``GenotypeGVCFsEngine.java:191-194`` and ``:479-491``).  ``INFO/DP`` survives
   from the merged record and the Number=A ``MLEAC``/``MLEAF`` vectors are
   written as missing values, giving the literal row

       chr1 2 . A . 127.78 . DP=20;MLEAC=.;MLEAF=. GT ./.

   ``GenotypeGVCFs.apply()`` writes it because the record is not
   "spanning-deletion only" -- it has no ALT at all
   (``GenotypeGVCFs.java:324-330``, ``GATKVariantContextUtils.java:2089-2091``).

The orphan-``*``-with-surviving-ALT cases
-----------------------------------------
``surviving-concrete-alt`` was carried as REPORTED ONLY until the round
documented in ``.diag/round-prefer-pls.md`` promoted it.  It pins the other half
of ``GATKVariantContextUtils.makeGenotypeCall()``'s PREFER_PLS split: when the
likelihood row that ``AlleleSubsettingUtils.subsetAlleles()`` projected onto the
KEPT alleles is not informative (``isInformative`` is ``sum(log10GL) < -0.1``,
``GATKVariantContextUtils.java:54-58``), the row is ignored and the call is
projected from the source genotype with ``bestMatchToOriginalGT()``
(``:333-338``, definition ``:397-403``) -- a surviving source allele keeps its
identity, every pruned allele becomes the reference, the source copy order and
phase bit are preserved, and no GQ is assigned.  The companion cases pin the
position-independent projection (``...-star-not-first``), the copy order
(``...-reversed-source-gt``) and the phase bit (``...-phased-source-gt``), which
together distinguish "project the source call" from "derive the call from the
PLs" (whose argmax over this fixture is the ``(*,G)`` cell).

Scope and comparison contract
-----------------------------
Pinned GATK and native run with identical arguments on the same plain
(unindexed .vcf.gz-refusing) VCF input and the same generated reference; the
data rows must be byte-identical.  Native additionally receives
``--gatk-compatible-annotations`` because the default native profile is an
explicitly non-GATK diagnostic output (it publishes RCQ/RCP and skips GATK's
AF-based output-allele pruning); every registered genotype-gvcf oracle uses the
same convention.

Exit status: 0 when every gated case passes, non-zero otherwise.
``--expect-divergence`` turns the run into a diagnostic that always exits 0.
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

# Copied verbatim from verify_genotype_gvcf.py:328-337 so the gate pins the
# same fixture bytes as the registered contract test.
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""

# verify_genotype_gvcf.py:378-381, byte for byte.
STAR_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100\n"
)

# The concrete ALT G owns the best genotype (PL index 4 = the (*,G) cell), so
# only the orphan '*' is pruned.  GATK keeps ALT='G' and publishes 0/1.
G_PLAUSIBLE_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Same locus with the orphan '*' in the middle of the ALT list, so the pruned
# allele is neither first nor last.  bestMatchToOriginalGT() projects the source
# call by allele identity, not by index, so the published row is unchanged.
G_STAR_SECOND_RECORD = (
    "chr1\t2\t.\tA\tG,*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Reported only: a concrete DELETION owns the best genotype.  GATK prunes the
# orphan '*' here too, but native keeps it in the ALT list, which is a separate
# deletion-ownership divergence (isVcCoveredByDeletion) rather than the
# PREFER_PLS subsetting rule this oracle pins.
G_DELETION_ALT_RECORD = (
    "chr1\t2\t.\tAA\t*,A,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# The same locus with the source call written in the opposite copy order.  The
# measured GATK row keeps that order, so the fallback must project the source
# genotype rather than any PL-argmax or sorted genotype.
G_REVERSED_GT_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t2/0:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# ... and the phased spelling of the same call, which GATK also preserves.
G_PHASED_GT_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t2|0:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

GATK_DEFAULT_ROW = "chr1\t2\t.\tA\t.\t127.78\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./."

# Measured with pinned GATK 4.6.2.0 (see .diag/round-prefer-pls.md).
GATK_G_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
              "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
              "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
GATK_G_REVERSED_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
                       "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                       "GT:AD:DP:PL\t1/0:0,20:20:0,0,0")
GATK_G_PHASED_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
                     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                     "GT:AD:DP:PL\t1|0:0,20:20:0,0,0")

CASES = [
    {
        "case": "default",
        "why": "the exact invocation of verify_genotype_gvcf.py:383-386 (no "
               "option flags): GATK's GenotypingEngine returns null for this "
               "monomorphic-after-pruning locus and writes no record, so the "
               "stale assertion at line 391 pinned native-only output",
        "body": STAR_RECORD,
        "args": [],
        "gated": True,
        "expect": [],
    },
    {
        "case": "include-non-variant-sites",
        "why": "same fixture with OutputMode.EMIT_ALL_ACTIVE_SITES: the "
               "REF-only no-call row survives and is the measured GATK contract "
               "for the audit's '-all-sites emits ALT=. GT=./.' claim",
        "body": STAR_RECORD,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEFAULT_ROW],
    },
    {
        "case": "stand-call-conf-50-include-non-variant-sites",
        "why": "the standard-confidence option cannot change the outcome here: "
               "the ALT set is already empty, so the confidence value is "
               "irrelevant to allele pruning and the row is unchanged",
        "body": STAR_RECORD,
        "args": ["--standard-min-confidence-threshold-for-calling", "50",
                 "--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEFAULT_ROW],
    },
    {
        "case": "surviving-concrete-alt",
        "why": "GATK keeps ALT='G' and calls 0/1.  The subset projection onto the "
               "kept alleles [A,G] collapses the source PL row [100,100,100] to a "
               "constant, so GATKVariantContextUtils.makeGenotypeCall() takes its "
               "PREFER_PLS fallback: isInformative is false "
               "(GATKVariantContextUtils.java:54-58; :333-338) and the call is "
               "projected from the SOURCE genotype with bestMatchToOriginalGT() "
               "(:397-403), which keeps every source allele that survived the "
               "subset and replaces a pruned one by the reference.  No GQ is "
               "assigned on that branch, and the published PL row is the "
               "min-shifted projection (:90-95, GenotypeLikelihoods.GLsToPLs) = "
               "0,0,0.  AC/AF/AN then follow from that genotype through the "
               "StandardAnnotation pass: ChromosomeCounts counts the called "
               "alleles of the published genotypes "
               "(ChromosomeCounts.java:43-53 -> "
               "VariantContextUtils.calculateChromosomeCounts) after "
               "GenotypeGVCFsEngine.regenotypeVC() has finished them "
               "(GenotypeGVCFsEngine.java:187-190)",
        "body": G_PLAUSIBLE_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_ROW],
    },
    {
        "case": "surviving-concrete-alt-star-not-first",
        "why": "the same PREFER_PLS fallback with the pruned '*' in the middle of "
               "the ALT list: the projection is by allele identity, so the row is "
               "byte-identical to surviving-concrete-alt and pins that GATK does "
               "not depend on the position of the dropped allele",
        "body": G_STAR_SECOND_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_ROW],
    },
    {
        "case": "surviving-concrete-alt-reversed-source-gt",
        "why": "the measured contract for the copy order: bestMatchToOriginalGT() "
               "maps the source allele LIST position by position "
               "(GATKVariantContextUtils.java:397-403), so a source call written "
               "2/0 yields 1/0 and not a canonicalized 0/1.  This case separates "
               "'project the source genotype' from 'derive the genotype from the "
               "PLs' (whose argmax genotype would be (*,G) here)",
        "body": G_REVERSED_GT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_REVERSED_ROW],
    },
    {
        "case": "surviving-concrete-alt-phased-source-gt",
        "why": "GenotypeBuilder retains the source genotype's phase flag "
               "(GATKVariantContextUtils.java:397-403 rebuilds the allele list "
               "over the same builder), so a phased source call stays phased "
               "after the projection",
        "body": G_PHASED_GT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_PHASED_ROW],
    },
    {
        "case": "surviving-deletion-alt",
        "why": "reported only: with a concrete deletion owning the best genotype "
               "GATK also prunes the orphan '*' and publishes ALT='A' with 0/1, "
               "but native keeps '*' in the ALT list.  That is a separate "
               "deletion-ownership (isVcCoveredByDeletion) divergence, not the "
               "PREFER_PLS subsetting rule this oracle gates, so it is carried "
               "without gating",
        "body": G_DELETION_ALT_RECORD,
        "args": [],
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
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
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
    source = work / f"{case['case']}.g.vcf"
    source.write_text(HEADER + case["body"], encoding="utf-8")
    # GATK's --include-non-variant-sites switches to a group-by-locus traversal
    # that requires an index; a plain .vcf gets a tribble .idx.  Creating it for
    # every case keeps the two tools' inputs byte-identical.
    index_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                           "IndexFeatureFile", "-I", str(source)],
                          f"GATK IndexFeatureFile [{case['case']}]", timeout)
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                          "GenotypeGVCFs", *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    # Native's default profile is a non-GATK diagnostic output, so the
    # GATK-compatibility writer is the only byte-comparable surface here; the
    # semantic arguments above are identical on both sides.
    native_result = invoke([str(native), *common,
                            "--gatk-compatible-annotations",
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

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
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
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
        description="Strict oracle for GenotypeGVCFs on a gVCF record whose ALT "
                    "list carries the symbolic spanning deletion '*', against "
                    "pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"),
        help="native GenotypeGVCFs binary (default: $FASTGATK_GENOTYPE_BINARY "
             "or fastgatk-native/build/fastgatk-genotype-gvcf)")
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
        / "fastgatk-genotype-gvcf")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_genotype_gvcf_spandel_gatk_oracle.py", java, jar)
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
            prefix="fastgatk-genotype-spandel-oracle-") as directory:
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
    print("# GATK rule under test: an ALT survives only when it passes "
          "standardConfidenceForCalling and, for '*', when an emitted deletion "
          "covers the locus (GenotypingEngine.java:311-327); with an empty ALT "
          "set the default traversal writes nothing (:167-169) while "
          "--include-non-variant-sites forces the REF-only no-call row "
          "(GenotypeGVCFsEngine.java:191-194, :381-383; cleanup at :479-491).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args']) if result['args'] else '(none)'}")
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
