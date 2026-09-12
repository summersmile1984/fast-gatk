#!/usr/bin/env python3
"""Round-49 adversarial harness: pinned GATK 4.6.2.0 vs native fastgatk-genotype-gvcf.

Independent of the registered gate (verify_genotype_gvcf_reverse_trim_gatk_oracle.py):
own fixtures, own reference sequences, own comparison.  Scratch persists under
.diag/round49/work so every command line and row set can be re-inspected.

Usage:  python3 probe.py [fixture-name ...]        (default: all)
        python3 probe.py --list
"""
from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path("/home/turing-agents/Documents/fast-gatk")
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = ROOT / "fastgatk-native/build/fastgatk-genotype-gvcf"
WORK = pathlib.Path(__file__).resolve().parent / "work"
DENSE = "--include-non-variant-sites"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=400>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""

GT = "GT:DP:AD:PL"
HET2 = "0/1:20:0,20,0:100,0,100,100,100,100"          # 2 alleles
HET3 = "0/1:20:0,20,0,0:100,0,100,100,100,100,100,100,100,100"
HET4 = "0/1:20:0,20,0,0,0:100,0,100,100,100,100,100,100,100,100,100,100,100,100,100"
HOMREF2 = "0/0:20:20,0,0:0,100,100,100,100,100"
HOMREF3 = "0/0:20:20,0,0,0:0,100,100,100,100,100,100,100,100,100"
STAR2_HET = "0/1:20:0,20:100,0,100"                    # '*' + <NON_REF>, diploid
HOM_STAR = "1/1:20:0,20:100,100,0"
AAA = "A" * 400

FIXTURES: list[dict] = []


def fixture(name, why, body, reference=AAA, args=None, header=HEADER, native_args=None):
    FIXTURES.append(dict(name=name, why=why, body=body, reference=reference,
                         args=args or [], header=header,
                         native_args=native_args or []))


# --- axis (b): multi-base REF with a symbolic-only ALT list -------------------
fixture(
    "b-block-multibase-ref-default",
    "axis b: REF=AAAA with ALT=<NON_REF> only (a reference block); GATK's "
    "regenotypeVC takes the `else result = originalVC` arm because the record "
    "has no non-symbolic non-ref allele, so no trim should be applied",
    "chr1\t2\t.\tAAAA\t<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n")

fixture(
    "b-block-multibase-ref-dense",
    "axis b+f: same block in dense mode",
    "chr1\t2\t.\tAAAA\t<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n",
    args=[DENSE])

fixture(
    "b-block-multibase-ref-end",
    "axis a+b: block with INFO/END (ACGT-repeat reference so REF matches)",
    "chr1\t2\t.\tACGT\t<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n",
    reference="ACGT" * 100)

fixture(
    "b-block-plus-symbolic-del",
    "axis b: REF=ACGT with ALTs <NON_REF>,<DEL> -- still no concrete ALT, so "
    "still a non-symbolic-free allele list (Type=SYMBOLIC)",
    "chr1\t2\t.\tACGT\t<NON_REF>,<DEL>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0,0:0,100,100,100,100,100,100\n",
    reference="ACGT" * 100)

# --- axis (a): INFO/END with a multi-base REF and a concrete ALT -------------
fixture(
    "a-end-with-concrete-alt",
    "axis a: multi-base REF, concrete ALT, and INFO/END present; does the trim "
    "rewrite END or leave the attribute stale in either tool?",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t" + HET2 + "\n")

fixture(
    "a-end-with-concrete-alt-dense",
    "axis a+f: same with dense mode",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t" + HET2 + "\n",
    args=[DENSE])

# --- reach: the DP>0 precondition of GenotypeGVCFsEngine --------------------
fixture(
    "r-no-info-dp",
    "reach: GATK regenotypes only when originalVC.getAttributeAsInt(DP,0) > 0, "
    "so a record without INFO/DP never reaches :167; both tools should agree on "
    "whatever they emit",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\t.\tGT:DP:AD:PL\t" + HET2 + "\n")

fixture(
    "r-info-dp-zero-dense",
    "reach: INFO/DP=0 in dense mode -- GATK's `result = originalVC` arm, so the "
    "alleles must come out UNTRIMMED",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=0\tGT:DP:AD:PL\t" + HET2 + "\n",
    args=[DENSE])

fixture(
    "r-info-dp-zero-default",
    "reach: INFO/DP=0 in default mode",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=0\tGT:DP:AD:PL\t" + HET2 + "\n")

# --- axis (c): spanning deletion with a multi-base REF -----------------------
fixture(
    "c-star-multibase-ref-no-upstream",
    "axis c: a '*' record with REF=AA and no covering upstream deletion",
    "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n")

fixture(
    "c-star-multibase-ref-with-upstream",
    "axis c: upstream concrete deletion AAAA>AA at 2 spanning 2-5, then a '*' "
    "record with REF=AA at 4",
    "chr1\t2\t.\tAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
    "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n")

# --- axis (g): placement of the trim vs EmittedDeletions::record -------------
fixture(
    "g-placement-deletion-size",
    "axis g: AAAA>AA at 2 emits a deletion of size 2 recorded BEFORE the trim "
    "(GenotypingEngine :178-179).  Trimming first would emit AAA>A, which is a "
    "deletion of size 2 as well, but the recorded span differs in length: "
    "untrimmed 2..4 vs trimmed 2..4 -- the later '*' at 4 must stay owned.",
    "chr1\t2\t.\tAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
    "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n",
    args=[DENSE])

fixture(
    "g-placement-long-ref-star-4",
    "axis g: AAAAAA>AA at 2 (deletion size 4, span 2..6) then a '*' record with "
    "REF=AAAA at 4; the trim of the emitted list is REF-clipping, so the "
    "recorded span must stay the untrimmed one",
    "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
    "chr1\t4\t.\tAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n",
    args=[DENSE])

fixture(
    "g-star-at-emitted-span-end",
    "axis g: deletion AAAAAA>AA at 2 spans 2..6 untrimmed; a '*' record with "
    "REF=AAAAA at 6 sits exactly on the last covered base",
    "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
    "chr1\t6\t.\tAAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n",
    args=[DENSE])

# --- axis (d): merged loci with a surviving one-base ALT ---------------------
HEADER2 = HEADER.replace("\tSTAR\n", "\tSTAR\tSTAR2\n")
HEADER3 = HEADER.replace("\tSTAR\n", "\tSTAR\tSTAR2\tSTAR3\n")

fixture(
    "d-two-sample-one-base-alt-survives",
    "axis d: sample STAR2 carries a strongly-supported one-base ALT A at 2 "
    "(AAAA>A); it survives the output-allele subset, so the guard at :1458 must "
    "protect the whole record from the trim",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\t" + GT + "\t" + HET2 + "\t" + HET2 + "\n"
    "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t" + GT + "\t" + HET2 + "\t" + HET2 + "\n",
    header=HEADER2)

fixture(
    "d-three-sample-one-base-alt-pruned",
    "axis d: three samples; the one-base ALT carried by STAR2 must be pruned by "
    "the subset before the guard is consulted, so the record IS trimmed",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\t" + GT + "\t" + HET2 + "\t" + HET2 + "\t" + HET2 + "\n"
    "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t" + GT + "\t" + HOMREF2 + "\t" + HOMREF2 + "\t" + HOMREF2 + "\n",
    header=HEADER3)

# --- axis (e): polyploid Number=G ------------------------------------------
def fixture_ploidy(name, ploidy, gt_fields, why):
    fixture(name, why,
            "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + gt_fields + "\n",
            args=["--sample-ploidy", str(ploidy)])


fixture_ploidy("e-ploidy3-trim", 3, HET3,
               "axis e: ploidy 3, Number=G width 6, trailing-run clip APPLIES")
fixture_ploidy("e-ploidy4-trim", 4, HET4,
               "axis e: ploidy 4, Number=G width 10, trailing-run clip APPLIES")
fixture_ploidy("e-ploidy3-guard", 3, "0/1:20:0,20,0:100,0,100,100,100,100".replace(":0,20,0:", ":0,20,0:"),
               "axis e: ploidy 3 with a one-base ALT that prunes (guard must fire after subset)")

fixture(
    "e-ploidy3-one-base-alt",
    "axis e: ploidy 3 with a surviving one-base ALT A at 2",
    "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/1/1:20:0,20,0:100,0,100,100,100,100\n",
    args=["--sample-ploidy", "3"])

# --- axis (h): streaming path ----------------------------------------------
fixture(
    "h-stream-suffix-substitution",
    "axis h: the same locus through --stream-by-locus; the helper is reached "
    "from the streaming compute stage as well",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n",
    native_args=["--stream-by-locus"])

fixture(
    "h-stream-star-locus-group",
    "axis h: a two-record locus group (upstream deletion + downstream '*') "
    "through --stream-by-locus; the bounded k-way merge + shared helper",
    "chr1\t2\t.\tAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
    "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOM_STAR + "\n",
    native_args=["--stream-by-locus"])

fixture(
    "h-stream-dense",
    "axis h+f: dense mode through --stream-by-locus",
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n",
    args=[DENSE], native_args=["--stream-by-locus"])


def write_reference(path: pathlib.Path, sequence: str) -> None:
    path.write_text(">chr1\n" + sequence + "\n", encoding="utf-8")
    line_bases = 400
    path.with_name(path.name + ".fai").write_text(
        f"chr1\t{len(sequence)}\t6\t{line_bases}\t{line_bases + 1}\n", encoding="utf-8")
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:{len(sequence)}\n", encoding="utf-8")


def data_rows(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def run(command: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(command, text=True, capture_output=True, check=False, timeout=600)


def run_fixture(case: dict) -> bool:
    WORK.mkdir(parents=True, exist_ok=True)
    name = case["name"]
    reference = WORK / f"{name}.fa"
    write_reference(reference, case["reference"])
    source = WORK / f"{name}.g.vcf"
    source.write_text(case["header"] + case["body"], encoding="utf-8")

    gatk_out = WORK / f"{name}.gatk.vcf"
    native_out = WORK / f"{name}.native.vcf"
    for stale in (gatk_out, native_out):
        if stale.exists():
            stale.unlink()

    index = run([str(JAVA), "-Xmx1g", "-jar", str(JAR), "IndexFeatureFile", "-I", str(source)])
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_cmd = [str(JAVA), "-Xmx1g", "-jar", str(JAR), "GenotypeGVCFs", *common,
                "-O", str(gatk_out), "--create-output-variant-index", "false"]
    native_cmd = [str(NATIVE), *common, "--gatk-compatible-annotations",
                  *case["native_args"], "-O", str(native_out)]
    gatk = run(gatk_cmd)
    native = run(native_cmd)

    gatk_rows = data_rows(gatk_out)
    native_rows = data_rows(native_out)
    print("=" * 100)
    print(f"# {name}\n# why: {case['why']}")
    print(f"# gatk  : {' '.join(gatk_cmd)}   [exit {gatk.returncode}]")
    print(f"# native: {' '.join(native_cmd)}   [exit {native.returncode}]")
    if index.returncode != 0:
        print(f"# IndexFeatureFile exit {index.returncode}: {index.stderr[-500:]}")
    if gatk.returncode != 0:
        print(f"# GATK stderr tail: {gatk.stderr[-1200:]}")
    if native.returncode != 0:
        print(f"# native stderr tail: {native.stderr[-1200:]}")
    for label, rows in (("GATK", gatk_rows), ("NATIVE", native_rows)):
        if not rows:
            print(f"{label:>6} rows: <none>")
        for row in rows:
            print(f"{label:>6}      {row}")
    verdict = "EQUAL" if gatk_rows == native_rows else "DIFFER"
    print(f"# VERDICT: {verdict}")
    return verdict == "EQUAL"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="*")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.list:
        for case in FIXTURES:
            print(case["name"])
        return 0
    selected = [case for case in FIXTURES if not args.names or case["name"] in args.names]
    if not selected:
        raise SystemExit(f"no such fixture: {args.names}")
    results = {case["name"]: run_fixture(case) for case in selected}
    print("=" * 100)
    for name, ok in results.items():
        print(f"{'EQUAL ' if ok else 'DIFFER'}  {name}")
    print(f"# {sum(results.values())}/{len(results)} fixtures byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
