#!/usr/bin/env python3
"""Round-49 batch 3: polyploid without --sample-ploidy, the merged-locus guard
case in dense mode, and the dense-block REF-source attribution probe."""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import probe  # noqa: E402

DENSE = probe.DENSE
CASES: list[dict] = []


def case(name, why, body, reference=probe.AAA, args=None, compare_positions=None,
         header=probe.HEADER, native_args=None):
    CASES.append(dict(name=name, why=why, body=body, reference=reference,
                      args=args or [], header=header, native_args=native_args or [],
                      compare_positions=compare_positions))


# --- axis (e): polyploid samples, no --sample-ploidy --------------------------
case("p3-nop-triploid-gt",
     "axis e: the native binary rejects --sample-ploidy, so the polyploid record "
     "is compared with the ploidy carried by the GT field on BOTH sides; "
     "Number=G width 10 for a triploid 3-allele record",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/1/1:20:0,20,0:100,0,100,100,100,100,100,100,100,100\n")

case("p4-nop-tetraploid-gt",
     "axis e: tetraploid GT, Number=G width 15",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/1/1/1:20:0,20,0:" + ",".join(["100", "0"] + ["100"] * 13) + "\n")

case("p3-nop-guard",
     "axis e: triploid with a surviving one-base ALT (guard)",
     "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/1/1:20:0,20,0:100,0,100,100,100,100,100,100,100,100\n")

case("e-wide-diploid-3alt",
     "axis e: a wider Number=G array at diploid ploidy -- two concrete ALTs, "
     "so the pre-trim PL row is 10 values wide",
     "chr1\t2\t.\tAAAA\tAACA,AA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t"
     "0/1:20:0,20,0,0:100,0,100,100,100,100,100,100,100,100\n")

# --- axis (d): merged locus where the one-base ALT survives -------------------
HEADER2 = probe.HEADER.replace("\tSTAR\n", "\tSTAR\tSTAR2\n")

case("d2-merged-two-het-dense",
     "axis d in DENSE mode, where GATK really does merge same-position records: "
     "both records are het for their own ALT, so the merged emitted list "
     "[AAAA, AACA, A] contains a one-base concrete ALT and the guard at :1458 "
     "must protect the whole record from the trim",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + probe.HET2 + "\t" + probe.HET2 + "\n"
     "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + probe.HET2 + "\t" + probe.HET2 + "\n",
     header=HEADER2, args=[DENSE], compare_positions=[2])

case("d2-merged-two-het-default",
     "axis d in DEFAULT mode: GATK's by-variant traversal regenotypes each input "
     "record separately (VariantLocusWalker.traverse + "
     "changeTraversalModeToByVariant), so the two records are never merged and "
     "the one-base-ALT record is returned untrimmed by the guard",
     "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + probe.HET2 + "\t" + probe.HET2 + "\n"
     "chr1\t2\t.\tAAAA\tA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + probe.HET2 + "\t" + probe.HET2 + "\n",
     header=HEADER2)

# --- attribution probe for the dense-block REF --------------------------------
case("x-attrib-block-ref-base",
     "attribution probe (input REF deliberately disagrees with the reference at "
     "this locus: reference 2-5 is CGTA but the record says ACGT).  The block's "
     "dense-mode REF therefore discriminates the three candidate sources: "
     "FASTA base at POS ('C'), the reverse trim of REF ('A'), or no rewrite "
     "('ACGT')",
     "chr1\t2\t.\tACGT\t<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/0:20:20,0:0,100,100\n",
     reference="ACGT" * 100, args=[DENSE], compare_positions=[2])


def rows_at(path, positions):
    return [r for r in probe.data_rows(path) if int(r.split("\t")[1]) in positions]


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
    probe.run([str(probe.JAVA), "-Xmx1g", "-jar", str(probe.JAR),
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
        gatk_rows, native_rows = rows_at(gatk_out, set(positions)), rows_at(native_out, set(positions))
    else:
        gatk_rows, native_rows = probe.data_rows(gatk_out), probe.data_rows(native_out)
    print("=" * 100)
    print(f"# {name}\n# why: {case['why']}")
    print(f"# gatk  : {' '.join(gatk_cmd)}   [exit {gatk.returncode}]")
    print(f"# native: {' '.join(native_cmd)}   [exit {native.returncode}]")
    if gatk.returncode != 0:
        print(f"# GATK stderr tail: {gatk.stderr[-1200:]}")
    if native.returncode != 0:
        print(f"# native stderr tail: {native.stderr[-1200:]}")
    if positions:
        print(f"# GATK positions: {[r.split(chr(9))[1] for r in probe.data_rows(gatk_out)]}")
        print(f"# native positions: {[r.split(chr(9))[1] for r in probe.data_rows(native_out)]}")
        print(f"# compared positions: {positions}")
    for label, rs in (("GATK", gatk_rows), ("NATIVE", native_rows)):
        if not rs:
            print(f"{label:>6} rows: <none>")
        for r in rs:
            print(f"{label:>6}      {r}")
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
