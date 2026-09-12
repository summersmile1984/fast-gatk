#!/usr/bin/env python3
"""Wide-band differential probe: many allele shapes in ONE GATK run.

The registered gate pins seven hand-picked shapes.  This probe instead puts
~50 loci in a single gVCF so one GATK startup covers all of them, and compares
three things per locus:

  * GATK's emitted REF/ALT  (the pinned 4.6.2.0 oracle),
  * native's emitted REF/ALT,
  * a Python transcription of the native helper apply_gatk_reverse_trim()
    applied to the input allele list, so a disagreement localises to the
    helper's algorithm rather than to the surrounding pipeline.

Usage: python3 wide.py [--stage A|B]
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

ROOT = pathlib.Path("/home/turing-agents/Documents/fast-gatk")
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = ROOT / "fastgatk-native/build/fastgatk-genotype-gvcf"
WORK = pathlib.Path(__file__).resolve().parent / "wide"

# A deterministic ACGT sequence with no long homopolymers, so every REF allele
# below is exactly the reference's own bases at that locus.
_BASES = "ACGTTGCACGATCGTTAGCACTGATCGTTGCAACGATTCGAGCTAGCTTCGATCGTAGCTAGCATCG"
REFERENCE = (_BASES * 30)[:1600]

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1600>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""


# --- transcription of apply_gatk_reverse_trim() ------------------------------
def symbolic(allele: str) -> bool:
    return allele.startswith("<")


def model_trim(alleles: list[str]) -> list[str]:
    if len(alleles) <= 1:
        return list(alleles)
    for allele in alleles:                                       # GATKVariantContextUtils:1458
        if not symbolic(allele) and allele != "*" and len(allele) == 1:
            return list(alleles)
    candidates = [a for a in alleles if not symbolic(a) and a != "*"]
    if not candidates:
        return list(alleles)
    shortest = min(len(a) for a in candidates)
    if shortest == 0:
        return list(alleles)
    end_trim = 0
    while end_trim < shortest and len({a[len(a) - 1 - end_trim] for a in candidates}) == 1:
        end_trim += 1
    rev_trim = end_trim - 1 if end_trim == shortest else end_trim  # :1469-1475
    if rev_trim == 0:
        return list(alleles)
    return [a if (symbolic(a) or a == "*") else a[:len(a) - rev_trim] for a in alleles]


# --- shape generators --------------------------------------------------------
def other(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def shapes(reference: str, start: int, length: int) -> list[tuple[str, str, list[str]]]:
    ref = reference[start - 1:start - 1 + length]
    assert len(ref) == length, (start, length)
    out: list[tuple[str, str, list[str]]] = []
    for k in (0, 1, 2, 3):                       # common trailing run of exactly k
        if length - k - 1 < 0:
            continue
        alt = ref[:length - k - 1] + other(ref[length - k - 1]) + ref[length - k:]
        out.append((f"run{k}", ref, [alt]))
    for d in (1, 2, 3):                          # trailing deletions
        if length - d >= 1:
            out.append((f"del{d}", ref, [ref[:length - d]]))
    out.append(("snp", ref, [ref[0]]))           # one-base ALT: the :1458 guard
    out.append(("ins-repeat", ref, [ref + ref[-1]]))
    out.append(("ins-diff", ref, [ref + other(ref[-1])]))
    out.append(("ins-shared", ref, [ref + ref[-2:]]))
    return out


def build_stage_a() -> list[tuple[str, int, str, list[str]]]:
    """Every record carries the symbolic <NON_REF>; GATK rejects gVCF input
    without it ("The list of input alleles must contain <NON_REF>")."""
    loci: list[tuple[str, int, str, list[str]]] = []
    position = 12
    for length in (2, 3, 4, 5, 6, 7):
        for label, ref, alts in shapes(REFERENCE, position, length):
            if position + length + 3 >= len(REFERENCE):
                break
            loci.append((f"L{length}-{label}", position, ref, [*alts, "<NON_REF>"]))
            position += 12                       # > any REF or deletion span
    return loci


def build_stage_b() -> list[tuple[str, int, str, list[str]]]:
    """Two concrete ALT alleles per record plus <NON_REF>."""
    loci: list[tuple[str, int, str, list[str]]] = []
    position = 12
    for length in (4, 5, 6, 8):
        for label, ref, alts in shapes(REFERENCE, position, length):
            if position + length + 3 >= len(REFERENCE):
                break
            second = None
            for index in (1, 2, 3, 0):
                if index >= length:
                    continue
                candidate = ref[:index] + other(ref[index]) + ref[index + 1:]
                if candidate != ref and candidate not in alts:
                    second = candidate
                    break
            if second is None:
                continue
            loci.append((f"L{length}-{label}", position, ref, [*alts, second, "<NON_REF>"]))
            position += 12
    return loci


def pl_vector(n_alleles: int, zero_index: int = 1) -> str:
    width = n_alleles * (n_alleles + 1) // 2
    values = [100] * width
    values[zero_index] = 0
    return ",".join(str(v) for v in values)


def write_inputs(loci, name):
    WORK.mkdir(parents=True, exist_ok=True)
    reference = WORK / f"{name}.fa"
    reference.write_text(">chr1\n" + REFERENCE + "\n", encoding="utf-8")
    length = len(REFERENCE)
    reference.with_name(reference.name + ".fai").write_text(
        f"chr1\t{length}\t6\t{length}\t{length + 1}\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:{length}\n", encoding="utf-8")
    source = WORK / f"{name}.g.vcf"
    lines = []
    for _label, position, ref, alts in loci:
        alleles = [ref, *alts]
        if len(alleles) >= 4:
            # REF, two concrete ALTs, <NON_REF>: call 1/2 so BOTH concrete
            # ALTs are supported and neither is pruned by the subset.
            gt = "1/2"
            ad = "0,10,10" + ",0" * (len(alleles) - 3)
            pl = pl_vector(len(alleles), zero_index=4)
        else:
            gt = "0/1"
            ad = ",".join(["0", "20"] + ["0"] * (len(alleles) - 2))
            pl = pl_vector(len(alleles), zero_index=1)
        lines.append(f"chr1\t{position}\t.\t{ref}\t{','.join(alts)}\t.\tPASS\tDP=20\t"
                     f"GT:DP:AD:PL\t{gt}:20:{ad}:{pl}")
    source.write_text(HEADER + "\n".join(lines) + "\n", encoding="utf-8")
    return reference, source


def run(command):
    return subprocess.run(command, text=True, capture_output=True, check=False, timeout=2400)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", default="A")
    args = parser.parse_args()
    loci = build_stage_a() if args.stage == "A" else build_stage_b()
    name = f"wide{args.stage}"
    reference, source = write_inputs(loci, name)
    run([str(JAVA), "-Xmx1g", "-jar", str(JAR), "IndexFeatureFile", "-I", str(source)])
    gatk_out = WORK / f"{name}.gatk.vcf"
    native_out = WORK / f"{name}.native.vcf"
    gatk = run([str(JAVA), "-Xmx1g", "-jar", str(JAR), "GenotypeGVCFs",
                "-R", str(reference), "-V", str(source), "-O", str(gatk_out),
                "--create-output-variant-index", "false"])
    native = run([str(NATIVE), "-R", str(reference), "-V", str(source),
                  "--gatk-compatible-annotations", "-O", str(native_out)])
    print(f"# stage {args.stage}: {len(loci)} loci; GATK exit={gatk.returncode} "
          f"native exit={native.returncode}")
    if gatk.returncode != 0:
        print("# GATK stderr tail:", gatk.stderr[-1500:])
    if native.returncode != 0:
        print("# native stderr tail:", native.stderr[-1500:])

    def rows(path):
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        table = {}
        for line in text.splitlines():
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            table[int(fields[1])] = line
        return table

    gatk_rows, native_rows = rows(gatk_out), rows(native_out)
    differ = 0
    model_mismatch = 0
    print(f"  {'locus':>18} {'POS':>4} {'input':<26} {'GATK':<26} {'NATIVE':<26} "
          f"{'model':<26} verdict")
    for label, position, ref, alts in loci:
        g = gatk_rows.get(position)
        n = native_rows.get(position)
        concrete = [a for a in alts if a != "<NON_REF>"]
        predicted = model_trim([ref, *concrete])
        prediction = "/".join(predicted)
        if g is None and n is None:
            print(f"  {label:>18} {position:>4} {ref}/{','.join(alts):<24} "
                  f"{'(none)':<26} {'(none)':<26} {prediction:<26} both-empty")
            continue
        gatk_alleles = f"{g.split(chr(9))[3]}/{g.split(chr(9))[4]}" if g else "(none)"
        native_alleles = f"{n.split(chr(9))[3]}/{n.split(chr(9))[4]}" if n else "(none)"
        g_parts = g.split("\t") if g else None
        expected = f"{predicted[0]}/" + (",".join(predicted[1:]) or ".")
        if g_parts and g_parts[4] != "." and f"{g_parts[3]}/{g_parts[4]}" != expected:
            model_mismatch += 1
            verdict = "MODEL!=GATK"
        elif g is None or n is None or g != n:
            verdict = "DIFFER"
        else:
            verdict = "EQUAL"
        if verdict != "EQUAL":
            differ += 1
        print(f"  {label:>18} {position:>4} {ref}/{','.join(alts):<24} {gatk_alleles:<26} "
              f"{native_alleles:<26} {prediction:<26} {verdict}")
        if verdict != "EQUAL":
            print(f"        GATK   row: {g}")
            print(f"        native row: {n}")
    print(f"# stage {args.stage}: {len(loci)} loci, {differ} not EQUAL, "
          f"{model_mismatch} where the model disagrees with GATK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
