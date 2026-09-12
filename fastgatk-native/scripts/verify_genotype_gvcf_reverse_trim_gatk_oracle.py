#!/usr/bin/env python3
"""Strict pinned-GATK oracle for GenotypeGVCFs' reverse allele trim.

The rule under test
-------------------
``GenotypeGVCFsEngine`` trims the emitted alleles back to their minimal
representation for every locus it regenotypes, in BOTH traversal modes::

    if (originalVC.isVariant() && originalVC.getAttributeAsInt(DP_KEY, 0) > 0
            && (isProperlyPolymorphic(regenotypedVC) || includeNonVariants)) {
        ...
        finalizeAnnotations(...);                      // :165
        regenotypedVC = GATKVariantContextUtils.reverseTrimAlleles(regenotypedVC);  // :167
        ...
        annotateContext(regenotypedVC, ...);           // :189
    }

(``gatk-source/.../walkers/GenotypeGVCFsEngine.java:160-169``, reached through
``GATKVariantContextUtils.reverseTrimAlleles()`` at ``:1443-1445`` =
``trimAlleles(vc, trimForward=false, trimReverse=true)``).

Four properties of that rule are pinned here, because each is separately
falsifiable and a plausible implementation can get each of them wrong:

1. **The trim is a trailing-base rewrite of REF and every concrete ALT, and it
   never moves the record start.**  ``trimAlleles(vc, -1, revTrim)`` clips the
   last ``revTrim`` bases of each non-symbolic, non-``*`` allele
   (``:1504``) and rebuilds the record with the *same* start and
   ``stop = start + REF.length() - 1`` (``:1515-1518``).  Measured on
   ``2 AAAA AACA,<NON_REF>``: GATK emits ``chr1 2 . AAA AAC``, i.e. one base
   clipped from both alleles at an unchanged POS -- not a moved or
   normalised record.
2. **The clip length is the common trailing run, capped by the shortest
   candidate allele.**  ``AlignmentUtils.normalizeAlleles(..., maxShift=0,
   trim=true)`` (``AlignmentUtils.java:818-845``) trims while the last base of
   every candidate is equal *and* no candidate is empty yet, so the run is
   capped at the shortest allele's length; one base is then restored when that
   cap was reached and nothing was clipped from the front
   (``GATKVariantContextUtils.java:1469-1475``), so an allele can never be
   emptied.  The two corners are separately measured: ``ACGTACGT/ACGT`` has a
   common run of four but clips three (``2 ACGTA A``), while ``AAAA/AACA``
   clips one (``2 AAA AAC``).
3. **The guard is on the emitted allele list, after the output-allele subset.**
   ``:1458`` returns the record untouched when it has one allele or when ANY
   allele has ``length() == 1`` and is not ``*`` -- which is why SNP and
   HaplotypeCaller-minimal indel records are never trimmed.  Because
   ``calculateOutputAlleleSubset()`` runs earlier in ``calculateGenotypes()``
   (``GenotypingEngine.java:155``), the guard sees the *surviving* alleles: a
   merged locus whose one-base ALT was pruned by the standard-confidence
   threshold is trimmed even though the merged input carried that ALT
   (``trim-runs-after-output-allele-subset``).
4. **Symbolic alleles and ``*`` are copied untouched** (``:1497-1501``), so a
   ``*`` record with a multi-base REF keeps its ``*`` while its REF is clipped
   (``single-concrete-candidate-clipped-to-one-base``, dense mode, compared at
   the record's own start position only -- see the case's ``compare_positions``
   note).

Reach, measured
---------------
The gate's cases are default-mode loci that do not depend on the covered-locus
materialization gap (item (1) of the dense round), so a green result here is
evidence about the trim alone and not about the traversal.  The one dense case
compares only the positions that both tools publish from a record start.

Usage / exit status
-------------------
``--native`` / ``$FASTGATK_GENOTYPE_BINARY`` select the native binary
(``$FASTGATK_NATIVE_BUILD``/``fastgatk-native/build`` supplies the default), the
pinned GATK jar and the bundled JDK17 are read from ``third_party/``, all
scratch lives in a ``tempfile.TemporaryDirectory``, and the run ends with a
single-line JSON status payload.  Exit status: 0 when every gated case passes,
non-zero otherwise.  ``--expect-divergence`` turns the run into a diagnostic
that always exits 0 (used while the defect is unfixed).
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from oracle_guard import oracle_not_verified  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""

# Two-sample header for the case that pins the ordering of the trim against the
# output-allele subset.
HEADER_TWO_SAMPLES = HEADER.replace(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR\tSTAR2")

DENSE = "--include-non-variant-sites"

# Strong heterozygous evidence: PL[het] = 0 with a fully supporting AD column.
GT_FIELDS = "GT:DP:AD:PL"
HET_STRONG = "0/1:20:0,20,0:100,0,100,100,100,100"
HET_STRONG_STAR = "0/2:20:0,0,20:100,100,100,100,0,100,100,100,100,100"
HOM_REF_STRONG = "0/0:20:20,0,0:0,100,100,100,100,100"

# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

# Property 1 and 2 (one-base clip).  Both alleles end in 'A' and the second
# base from the end differs, so exactly one base is clipped from REF and ALT.
SUFFIX_SUBSTITUTION = (
    f"chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n")

# Property 2 (no clip): the two alleles share no trailing base, so
# lastBaseOnRightIsSame() fails at once and the record is returned unchanged
# (:1492-1493).
NO_SHARED_TRAILING_BASE = (
    f"chr1\t2\t.\tAAAA\tAACC,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n")

# Property 3 (guard on a one-base ALT): HaplotypeCaller's minimal
# representation of this deletion is AAA>A, whose ALT has length() == 1, so
# :1458 returns the record untouched.
ONE_BASE_ALT_GUARD = (
    f"chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n")

# Property 2 (the emptiness guard): on an ACGT reference, REF = ACGTACGT and
# ALT = ACGT share a four-base trailing run, so normalizeAlleles() reaches
# minSize == 0 and endShift == 4; restoreOneBaseAtEnd then clips 4-1 = 3 bases,
# leaving REF = ACGTA and ALT = A rather than an empty ALT.
EMPTINESS_GUARD = (
    f"chr1\t2\t.\tACGTACGT\tACGT,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n")

# Property 2 (clip length is the run length, not a fixed one base): a longer
# allele pair with the same one-base trailing run.
LONGER_ALLELES_ONE_BASE_RUN = (
    f"chr1\t2\t.\tAAAAA\tAAACA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n")

# Property 3 (the trim runs on the SURVIVING allele list): sample STAR carries
# AACA, sample STAR2 carries a one-base ALT 'A' that its own hom-ref likelihoods
# prune from the merged locus, so the emitted allele list has no one-base
# concrete allele left and the guard at :1458 does not fire.
TWO_SAMPLE_PRUNED_ONE_BASE_ALT = (
    f"chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\t{HET_STRONG}\n"
    f"chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HOM_REF_STRONG}\t{HOM_REF_STRONG}\n")

# Property 4 (a '*' survives and the REF is still clipped): one emitted
# deletion at 2 spanning 2-4, plus the shape the previous round measured -- a
# '*' record with a two-base REF at 4.  Its only non-symbolic, non-'*' allele
# is the REF itself, so the common run equals the REF length, the emptiness
# guard restores one base and REF is clipped from 'AA' to 'A' while the '*' is
# copied untouched (:1497-1501).
DENSE_STAR_TWO_BASE_REF = (
    f"chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t{HET_STRONG}\n"
    f"chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\t{GT_FIELDS}\t1/1:20:0,20:100,100,0\n")

# ---------------------------------------------------------------------------
# measured GATK rows (this gate's fixtures, pinned from a pinned-GATK run)
# ---------------------------------------------------------------------------

ROW_SUFFIX_SUBSTITUTION = (
    "chr1\t2\t.\tAAA\tAAC\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_NO_SHARED_TRAILING_BASE = (
    "chr1\t2\t.\tAAAA\tAACC\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_ONE_BASE_ALT = (
    "chr1\t2\t.\tAAA\tA\t92.60\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_EMPTINESS_GUARD = (
    "chr1\t2\t.\tACGTA\tA\t92.60\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_LONGER_ALLELES = (
    "chr1\t2\t.\tAAAA\tAAAC\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_TWO_SAMPLE = (
    "chr1\t2\t.\tAAA\tAAC\t190.50\t.\t"
    "AC=2;AF=0.500;AN=4;DP=20;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=4.76\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100\t0/1:0,20:20:99:100,0,100")
# The dense case is compared at the two positions where both tools publish a row
# from a record start; position 3 is the round-4 materialization gap and is out
# of this gate's scope on purpose.
ROW_DENSE_STAR_AT_2 = (
    "chr1\t2\t.\tAAA\tA\t92.60\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_DENSE_STAR_AT_4 = (
    "chr1\t4\t.\tA\t*\t0\tLowQual\t"
    "AC=2;AF=1.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=2;MLEAF=1.00;QD=0.00\t"
    "GT:AD:DP:GQ:PL\t1/1:0,20:20:99:100,100,0")

CASES = [
    {
        "case": "suffix-substitution-default",
        "why": "property 1+2: GATK clips one trailing base from REF and from "
               "every concrete ALT and leaves POS alone "
               "(GenotypeGVCFsEngine.java:167 -> "
               "GATKVariantContextUtils.java:1443-1445, :1489-1521)",
        "body": SUFFIX_SUBSTITUTION,
        "reference": "A",
        "args": [],
        "expect": [ROW_SUFFIX_SUBSTITUTION],
    },
    {
        "case": "no-shared-trailing-base-unchanged",
        "why": "property 2 negative control: no common trailing base, so "
               "normalizeAlleles() returns endShift == 0 and trimAlleles() "
               "returns the input untouched (:1492-1493)",
        "body": NO_SHARED_TRAILING_BASE,
        "reference": "A",
        "args": [],
        "expect": [ROW_NO_SHARED_TRAILING_BASE],
    },
    {
        "case": "one-base-alt-guard",
        "why": "property 3: any allele with length() == 1 that is not '*' "
               "disables the trim (:1458), which is why the registered "
               "HaplotypeCaller-minimal corpus never sees it",
        "body": ONE_BASE_ALT_GUARD,
        "reference": "A",
        "args": [],
        "expect": [ROW_ONE_BASE_ALT],
    },
    {
        "case": "emptiness-guard-clips-one-less",
        "why": "property 2 corner: the common trailing run (4) would empty the "
               "ALT, so restoreOneBaseAtEnd clips one base fewer "
               "(:1469-1475); POS stays 2 -- the trim is not a left "
               "normalisation",
        "body": EMPTINESS_GUARD,
        "reference": "ACGT",
        "args": [],
        "expect": [ROW_EMPTINESS_GUARD],
    },
    {
        "case": "longer-alleles-same-one-base-run",
        "why": "property 2: the clip length is the shared run and not a fixed "
               "one base for byte-length-5 alleles either",
        "body": LONGER_ALLELES_ONE_BASE_RUN,
        "reference": "A",
        "args": [],
        "expect": [ROW_LONGER_ALLELES],
    },
    {
        "case": "trim-runs-after-output-allele-subset",
        "why": "property 3 ordering: the one-base ALT of the second sample is "
               "pruned by the standard-confidence subset first "
               "(GenotypingEngine.java:155), so the guard at :1458 does not "
               "fire and GATK trims the emitted alleles",
        "body": TWO_SAMPLE_PRUNED_ONE_BASE_ALT,
        "header": HEADER_TWO_SAMPLES,
        "reference": "A",
        "args": [],
        "expect": [ROW_TWO_SAMPLE],
    },
    {
        "case": "single-concrete-candidate-clipped-to-one-base",
        "why": "property 4: with only the REF as a non-symbolic, non-'*' "
               "allele the common run equals the REF length, so the emptiness "
               "guard clips it to a single base and the '*' is copied "
               "untouched (:1497-1501)",
        "body": DENSE_STAR_TWO_BASE_REF,
        "reference": "A",
        "args": [DENSE],
        "compare_positions": [2, 4],
        "expect_at_positions": {2: ROW_DENSE_STAR_AT_2, 4: ROW_DENSE_STAR_AT_4},
    },
]


def write_reference(work: pathlib.Path, kind: str) -> pathlib.Path:
    """A 100 bp chr1 of the requested repeat, with its index and dictionary."""
    reference = work / f"reference-{kind}.fa"
    reference.write_text(">chr1\n" + (kind * 25 if kind != "A" else "A" * 100) + "\n",
                         encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def read_lines(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8").splitlines()


def data_rows(path: pathlib.Path) -> list[str]:
    return [line for line in read_lines(path) if line and not line.startswith("#")]


def rows_at_positions(rows: list[str], positions: list[int]) -> dict[int, str]:
    wanted = {str(position) for position in positions}
    selected: dict[int, str] = {}
    for row in rows:
        fields = row.split("\t")
        if fields[1] in wanted:
            selected[int(fields[1])] = row
    return selected


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-2000:])
    return result


def run_case(case: dict, work: pathlib.Path, java: pathlib.Path, jar: pathlib.Path,
             native: pathlib.Path, timeout: int) -> dict:
    reference = write_reference(work, case["reference"])
    source = work / f"{case['case']}.g.vcf"
    source.write_text(case.get("header", HEADER) + case["body"], encoding="utf-8")
    # --include-non-variant-sites switches GATK to a group-by-locus traversal
    # that requires an index; indexing both sides' identical input keeps the two
    # invocations comparable.
    index_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                           "IndexFeatureFile", "-I", str(source)],
                          f"GATK IndexFeatureFile [{case['case']}]", timeout)
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
                          *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    # Native's default profile is a non-GATK diagnostic output, so the
    # GATK-compatibility writer is the only byte-comparable surface; the
    # semantic arguments are identical on both sides.
    native_result = invoke([str(native), *common, "--gatk-compatible-annotations",
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    positions = case.get("compare_positions")
    if positions:
        gatk_compared = rows_at_positions(gatk_rows, positions)
        native_compared = rows_at_positions(native_rows, positions)
        expect = case["expect_at_positions"]
    else:
        gatk_compared = {index: row for index, row in enumerate(gatk_rows)}
        native_compared = {index: row for index, row in enumerate(native_rows)}
        expect = {index: row for index, row in enumerate(case["expect"])}

    result = {
        "case": case["case"],
        "why": case["why"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "compared_positions": positions,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    # GATK's own behaviour is pinned from the measurement, not re-derived.
    if gatk_compared != expect:
        result["violations"].append(
            f"GATK's rows moved away from the measured contract: "
            f"expected={expect} measured={gatk_compared}")
    if native_compared != gatk_compared:
        result["violations"].append(
            f"native rows differ from GATK: GATK={gatk_compared} "
            f"NATIVE={native_compared}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypeGVCFs' reverse allele trim "
                    "(GenotypeGVCFsEngine.java:160-169, "
                    "GATKVariantContextUtils.java:1443-1521) against pinned "
                    "GATK 4.6.2.0.")
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
        oracle_not_verified("verify_genotype_gvcf_reverse_trim_gatk_oracle.py", java, jar)
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-revtrim-oracle-") as directory:
        work = pathlib.Path(directory)
        for case in selected:
            results.append(run_case(case, work, java, jar, native, arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: GenotypeGVCFsEngine.java:167 calls "
          "GATKVariantContextUtils.reverseTrimAlleles(), which clips the common "
          "trailing run of the non-symbolic, non-'*' alleles, restores one base "
          "when that would empty an allele (:1469-1475), and rewrites the record "
          "at an unchanged start (:1515-1518).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        print(f"[{result['case']}] {result['args']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        if result["compared_positions"]:
            print(f"    compared positions: {result['compared_positions']}")
        print(f"    GATK   rows: {result['gatk_rows']}")
        print(f"    native rows: {result['native_rows']}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "mode": "strict" if strict_mode else "expect-divergence",
        "gatk_version": "4.6.2.0",
        "cases": results,
        "violations": violations,
        "rule": "GenotypeGVCFsEngine.java:160-169 (reverseTrimAlleles at :167, "
                "both traversal modes), GATKVariantContextUtils.java:1443-1445 "
                "(reverseTrimAlleles), :1458 (one-allele / one-base-allele "
                "guard), :1462-1467+AlignmentUtils.java:818-845 "
                "(normalizeAlleles endShift, capped by the shortest allele), "
                ":1469-1475 (restore one base), :1489-1521 (clip and rewrite at "
                "an unchanged start), :1497-1501 (symbolic and '*' untouched)",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
