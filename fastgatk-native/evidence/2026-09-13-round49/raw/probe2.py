#!/usr/bin/env python3
"""Round-49 batch 2: reach, INFO/END, polyploidy and the deletion-bookkeeping
placement constraint, with per-position comparison for dense-mode fixtures.

Independent of the registered gate.  Reuses probe.py's helpers.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import probe  # noqa: E402

DENSE = probe.DENSE
HEADER = probe.HEADER
GT = probe.GT
HET2 = probe.HET2
HOMREF2 = probe.HOMREF2
HOME_STAR = probe.HOM_STAR
AAA = probe.AAA

CASES: list[dict] = []


def case(name, why, body, reference=AAA, args=None, compare_positions=None,
         header=HEADER, native_args=None):
    CASES.append(dict(name=name, why=why, body=body, reference=reference,
                      args=args or [], header=header,
                      native_args=native_args or [],
                      compare_positions=compare_positions))


# --- INFO/END -----------------------------------------------------------------
case("end-no-trim-control",
     "control for the END axis: REF/ALT share no trailing base, so GATK never "
     "trims; isolates END pass-through from the trim",
     "chr1\t2\t.\tAAAA\tAACC,<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t" + HET2 + "\n")

case("end-trim-consistent",
     "END=4 equals the POST-trim span (2..4): the trim cannot make this END "
     "stale, yet GATK drops END in the merge and native appears to keep it",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20;END=4\tGT:DP:AD:PL\t" + HET2 + "\n")

case("end-guard-no-trim",
     "the guard fires (one-base ALT), so no trim; END=4 is consistent with "
     "the untrimmed REF=AAA span 2..4 either way",
     "chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20;END=4\tGT:DP:AD:PL\t" + HET2 + "\n")

case("end-block-dense",
     "a multi-base-REF block with END in dense mode: if the trim ran on an "
     "END-bearing record, htsjdk's validateStop would reject it "
     "(measured in TrimProbe/P2), so this checks that the merger really "
     "removes END before :167",
     "chr1\t2\t.\tACGT\t<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n",
     reference="ACGT" * 100, args=[DENSE], compare_positions=[2])

# --- reach --------------------------------------------------------------------
case("prune-refonly-dense",
     "every concrete ALT is pruned (hom-ref likelihoods), so GATK's emitted "
     "list is the REF alone and :1458's getNAlleles()<=1 arm returns it "
     "untrimmed; native keeps <NON_REF> in the merged allele list, which "
     "would make the list two alleles and trigger a clip",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n",
     args=[DENSE], compare_positions=[2])

case("prune-refonly-default",
     "same input in default mode (GATK drops it: not properly polymorphic)",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n")

case("dp1-dense",
     "the DP>0 gate boundary: INFO/DP=1 in dense mode",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=1\tGT:DP:AD:PL\t" + HET2 + "\n",
     args=[DENSE], compare_positions=[2])

case("dp0-dense",
     "INFO/DP=0 makes GATK skip regenotyping altogether, so the merged record "
     "is emitted with UNTRIMMED alleles",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=0\tGT:DP:AD:PL\t" + HET2 + "\n",
     args=[DENSE], compare_positions=[2])

case("star-only-nonref-default",
     "a '*' record whose only ALT is '*' plus <NON_REF>: isProperlyPolymorphic() "
     "is false, so default mode drops it",
     "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n")

case("star-only-nonref-dense",
     "the same '*' record in dense mode, where :167 is reached",
     "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[4])

# --- polyploid ----------------------------------------------------------------
case("p3-trim",
     "axis e: --sample-ploidy 3, Number=G width 10; trailing-run clip applies "
     "and the wider PL row must survive",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/0/1:20:0,20,0:100,0,100,100,100,100,100,100,100,100\n",
     args=["--sample-ploidy", "3"])

case("p4-trim",
     "axis e: --sample-ploidy 4, Number=G width 15",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/0/0/1:20:0,20,0:" + ",".join(["100", "0"] + ["100"] * 13) + "\n")

case("p3-guard",
     "axis e: ploidy 3 with a surviving one-base ALT (guard protects the record)",
     "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/0/1:20:0,20,0:100,0,100,100,100,100,100,100,100,100\n",
     args=["--sample-ploidy", "3"])

# --- axis (g): the placement constraint --------------------------------------
case("g-star-chain-dense",
     "axis g: AAAAAA>AA at 2 emits a deletion of size 4 (span 2..6 recorded "
     "BEFORE the trim); the '*' record at 4 with REF=AAAAA is owned by it and "
     "its own emitted deletion size (5-1=4, span 4..8) is what keeps the '*' "
     "record at 7 non-spurious.  Trimming before recordDeletions would give "
     "record 2 a REF of one base, hence deletion size 0, and the record at 7 "
     "would be pruned as an orphan spanning deletion.",
     "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t4\t.\tAAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n"
     "chr1\t7\t.\tAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[2, 4, 7])

case("g-longref-phantom",
     "axis g adjacency check: a 12-base REF with an 11-base concrete ALT "
     "(deletion size 1, span 2..3) plus <NON_REF>.  htsjdk gives a symbolic "
     "allele length()==0, so if <NON_REF> reached recordDeletions() either "
     "tool would record a spurious span; native would compute 12-9=3.  The "
     "'*' at 4 is covered by neither real span (2..3) and must therefore be "
     "treated as a spurious spanning deletion by BOTH tools.",
     "chr1\t2\t.\tAAAAAAAAAAAA\tAAAAAAAAAAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t4\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[2, 4])

case("g-longref-default",
     "the 12-base REF / 11-base ALT locus alone, default mode: the emptiness "
     "rule with a long shared run (clip 10 of 11, REF->AA, ALT->A)",
     "chr1\t2\t.\tAAAAAAAAAAAA\tAAAAAAAAAAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n")

# --- axis (h): streaming ------------------------------------------------------
case("h-stream-end",
     "axis h: the END fixture through --stream-by-locus",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20;END=5\tGT:DP:AD:PL\t" + HET2 + "\n",
     native_args=["--stream-by-locus"])

case("h-stream-placement-chain",
     "axis h+g: the three-record '*' chain through --stream-by-locus, where "
     "the locus group is built by the bounded k-way merge",
     "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t4\t.\tAAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n"
     "chr1\t7\t.\tAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], native_args=["--stream-by-locus"], compare_positions=[2, 4, 7])

case("h-stream-p3",
     "axis e+h: ploidy 3 through the streaming path",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/0/1:20:0,20,0:100,0,100,100,100,100,100,100,100,100\n",
     args=["--sample-ploidy", "3"], native_args=["--stream-by-locus"])


def rows_at(path: pathlib.Path, positions) -> list[str]:
    out = []
    for row in probe.data_rows(path):
        if int(row.split("\t")[1]) in positions:
            out.append(row)
    return out


def run_case(case: dict) -> bool:
    probe.WORK.mkdir(parents=True, exist_ok=True)
    name = case["name"]
    reference = probe.WORK / f"{name}.fa"
    probe.write_reference(reference, case["reference"])
    source = probe.WORK / f"{name}.g.vcf"
    source.write_text(case["header"] + case["body"], encoding="utf-8")
    gatk_out = probe.WORK / f"{name}.gatk.vcf"
    native_out = probe.WORK / f"{name}.native.vcf"
    for stale in (gatk_out, native_out):
        if stale.exists():
            stale.unlink()
    index = probe.run([str(probe.JAVA), "-Xmx1g", "-jar", str(probe.JAR),
                       "IndexFeatureFile", "-I", str(source)])
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_cmd = [str(probe.JAVA), "-Xmx1g", "-jar", str(probe.JAR), "GenotypeGVCFs",
                *common, "-O", str(gatk_out), "--create-output-variant-index", "false"]
    native_cmd = [str(probe.NATIVE), *common, "--gatk-compatible-annotations",
                  *case["native_args"], "-O", str(native_out)]
    gatk = probe.run(gatk_cmd)
    native = probe.run(native_cmd)
    positions = case["compare_positions"]
    if positions:
        gatk_rows = rows_at(gatk_out, set(positions))
        native_rows = rows_at(native_out, set(positions))
    else:
        gatk_rows = probe.data_rows(gatk_out)
        native_rows = probe.data_rows(native_out)
    print("=" * 100)
    print(f"# {name}\n# why: {case['why']}")
    print(f"# gatk  : {' '.join(gatk_cmd)}   [exit {gatk.returncode}]")
    print(f"# native: {' '.join(native_cmd)}   [exit {native.returncode}]")
    if index.returncode != 0:
        print(f"# IndexFeatureFile exit {index.returncode}: {index.stderr[-400:]}")
    if gatk.returncode != 0:
        print(f"# GATK stderr tail: {gatk.stderr[-1500:]}")
    if native.returncode != 0:
        print(f"# native stderr tail: {native.stderr[-1500:]}")
    if positions:
        print(f"# compared positions: {positions} (dense-mode row shape/materialization "
              f"differences at other positions are out of scope)")
        print(f"# GATK rows at all positions: {[r.split(chr(9))[1] for r in probe.data_rows(gatk_out)]}")
        print(f"# native rows at all positions: {[r.split(chr(9))[1] for r in probe.data_rows(native_out)]}")
    for label, rows in (("GATK", gatk_rows), ("NATIVE", native_rows)):
        if not rows:
            print(f"{label:>6} rows: <none>")
        for row in rows:
            print(f"{label:>6}      {row}")
    verdict = "EQUAL" if gatk_rows == native_rows else "DIFFER"
    print(f"# VERDICT: {verdict}")
    return verdict == "EQUAL"


def main() -> int:
    selected = [c for c in CASES if not sys.argv[1:] or c["name"] in sys.argv[1:]]
    results = {c["name"]: run_case(c) for c in selected}
    print("=" * 100)
    for name, ok in results.items():
        print(f"{'EQUAL ' if ok else 'DIFFER'}  {name}")
    print(f"# {sum(results.values())}/{len(results)} fixtures byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
