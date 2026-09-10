#!/usr/bin/env python3
"""Deep oracle for the HaplotypeCaller --alleles / GenotypeGivenAlleles boundary.

Complements verify_hc_alleles_gatk_oracle.py (one central forced SNP) with the
structural cases where AssemblyResultSet.addGivenAlleles() must agree exactly:

  * indels inside a homopolymer run (leftmost / rightmost / deletion form)
  * indels inside a tandem repeat
  * a given allele whose REF span crosses the -L interval end
  * two adjacent given alleles (MNP coalescing window)
  * two mutually overlapping given alleles
  * a given allele overlapping an *assembled* non-reference haplotype event
    (base-haplotype selection must skip that base haplotype, not the allele)

One shared synthetic reference/BAM keeps the run short; only the feature VCF
and the command line change per case.  Every case runs the pinned GATK 4.6.2.0
jar and the native binary and compares VCF data rows field-by-field.

Exit code is non-zero when any case diverges (a divergence is a real finding,
not a skip).  Use --case NAME to run a single case.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path

REFERENCE_LENGTH = 1500
HOMOPOLYMER_START = 600          # 1-based, 10 x 'A' -> 600..609
TANDEM_START = 900               # 1-based, 8 x 'CAG' -> 900..923
ASSEMBLED_DELETION_POS = 700     # 1-based anchor of an assembled 1bp deletion
HOMOPOLYMER_INTERVAL = "chr1:500-780"
TANDEM_INTERVAL = "chr1:800-1080"


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def build_reference() -> str:
    random.seed(11)
    bases = list("".join(random.choice("ACGT") for _ in range(REFERENCE_LENGTH)))
    bases[HOMOPOLYMER_START - 1:HOMOPOLYMER_START - 1 + 10] = list("A" * 10)
    bases[TANDEM_START - 1:TANDEM_START - 1 + 24] = list("CAG" * 8)
    return "".join(bases)


def write_inputs(work: Path) -> tuple[Path, Path, str]:
    reference_text = build_reference()
    reference = work / "ref.fa"
    reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
    (work / "ref.fa.fai").write_text(
        f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
        encoding="utf-8")

    header = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}",
              "@RG\tID:rg\tSM:S1"]
    placed: list[tuple[int, str]] = []
    index = 0

    def add_reads(first_start: int, count: int, length: int, step: int) -> None:
        nonlocal index
        for offset in range(count):
            start = first_start + offset * step
            sequence = reference_text[start - 1:start - 1 + length]
            placed.append((start, f"r{index}\t0\tchr1\t{start}\t60\t{length}M\t*\t0\t0\t"
                                  f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg"))
            index += 1

    add_reads(440, 30, 300, 3)      # covers the 500-780 homopolymer window
    add_reads(800, 30, 300, 3)      # covers the 800-1080 tandem window

    # A small group carrying an assembled 1bp deletion at ASSEMBLED_DELETION_POS
    # so the base-haplotype population has one real non-reference haplotype.
    for offset in range(6):
        start = 560 + offset * 4
        sequence = reference_text[start - 1:start - 1 + 300]
        ref_offset = ASSEMBLED_DELETION_POS - start
        sequence = sequence[:ref_offset] + sequence[ref_offset + 1:]
        placed.append((start, f"d{offset}\t0\tchr1\t{start}\t60\t299M\t*\t0\t0\t{sequence}\t"
                              f"{'I' * len(sequence)}\tRG:Z:rg"))
    lines = header + [line for _, line in sorted(placed, key=lambda item: item[0])]
    sam = work / "input.sam"
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return reference, sam, reference_text


def write_feature(work: Path, name: str, rows: list[tuple[str, int, str, str]],
                  root: Path) -> Path:
    """rows: (chrom, pos, ref, alt) tuples; one VCF record each."""
    path = work / f"{name}.vcf"
    body = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={REFERENCE_LENGTH}>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, position, ref, alt in rows:
        body.append(f"{chrom}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.")
    path.write_text("\n".join(body) + "\n", encoding="utf-8")
    bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
    tabix = root / "third_party/htslib-build/htslib-src/tabix"
    subprocess.run([str(bgzip), "-f", str(path)], check=True, stdout=subprocess.DEVNULL)
    compressed = path.with_suffix(path.suffix + ".gz")
    subprocess.run([str(tabix), "-f", "-p", "vcf", str(compressed)], check=True,
                   stdout=subprocess.DEVNULL)
    return compressed


def case_definitions(reference_text: str) -> list[dict]:
    homopolymer_ref = reference_text[HOMOPOLYMER_START - 1:HOMOPOLYMER_START + 9]
    tandem_ref = reference_text[TANDEM_START - 1:TANDEM_START + 11]
    assembled_ref = reference_text[ASSEMBLED_DELETION_POS - 1:ASSEMBLED_DELETION_POS + 1]
    return [
        {
            "name": "homopolymer-insertion-left",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", HOMOPOLYMER_START, "A", "AA")],
            "why": "one-base insertion written at the leftmost repeat base",
        },
        {
            "name": "homopolymer-insertion-right",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", HOMOPOLYMER_START + 9, "A", "AA")],
            "why": "same insertion written at the rightmost repeat base",
        },
        {
            "name": "homopolymer-deletion",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", HOMOPOLYMER_START, "AA", "A")],
            "why": "one-base deletion inside the homopolymer run",
        },
        {
            "name": "tandem-repeat-deletion-anchored",
            "interval": TANDEM_INTERVAL,
            "rows": [("chr1", TANDEM_START + 2, reference_text[TANDEM_START + 1:TANDEM_START + 5],
                      reference_text[TANDEM_START + 1])],
            "why": "anchor-base form of a one-CAG-unit deletion (no suffix trim in Event)",
        },
        {
            "name": "tandem-repeat-one-unit-literal",
            "interval": TANDEM_INTERVAL,
            "rows": [("chr1", TANDEM_START + 3, "CAGCAG", "CAG")],
            "why": "ref=2 repeat units / alt=1 unit: Event trims to an empty ALT (GATK aborts)",
        },
        {
            "name": "assembled-event-overlap",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", ASSEMBLED_DELETION_POS, assembled_ref,
                      assembled_ref[0])],
            "why": "given deletion re-stated where an assembled non-ref hap already has it",
        },
        {
            "name": "adjacent-snps",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", HOMOPOLYMER_START + 1, "A", "G"),
                     ("chr1", HOMOPOLYMER_START + 2, "A", "C")],
            "why": "two adjacent forced SNPs inside the MNP distance window",
        },
        {
            "name": "overlapping-features",
            "interval": TANDEM_INTERVAL,
            "rows": [("chr1", TANDEM_START, tandem_ref[:6], tandem_ref[0]),
                     ("chr1", TANDEM_START + 4, tandem_ref[4], "T")],
            "why": "two forced events whose reference spans overlap",
        },
        {
            "name": "region-boundary-span",
            "interval": "chr1:500-640",
            "rows": [("chr1", 635, reference_text[634:654], reference_text[634])],
            "why": "forced deletion REF span crosses the -L interval end",
        },
    ]


def compare(java_rows: list[list[str]], native_rows: list[list[str]]) -> dict:
    detail: dict = {"java_rows": len(java_rows), "native_rows": len(native_rows)}
    for index in range(max(len(java_rows), len(native_rows))):
        java_row = java_rows[index] if index < len(java_rows) else None
        native_row = native_rows[index] if index < len(native_rows) else None
        if java_row is None or native_row is None:
            detail.update({"first_diff_row": index, "first_diff_field": "ROW_COUNT",
                           "java_row": java_row, "native_row": native_row})
            return detail
        for field in range(max(len(java_row), len(native_row))):
            java_value = java_row[field] if field < len(java_row) else None
            native_value = native_row[field] if field < len(native_row) else None
            if java_value != native_value:
                detail.update({"first_diff_row": index, "first_diff_field": field,
                               "java_row": java_row, "native_row": native_row})
                return detail
    return detail


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--native", default=os.environ.get("FASTGATK_HC_BINARY"))
    parser.add_argument("--drop-alleles", action="store_true",
                        help="control run: omit --alleles entirely, so a case that "
                             "still diverges is NOT an injection-boundary divergence")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(args.native) if args.native else root / "fastgatk-native/build/fastgatk-hc-call"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    if not native.exists() or not gatk.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit("missing HC binary or pinned GATK 4.6.2.0 jar")
        print(json.dumps({"status": "skipped", "reason": "GATK jar or native binary absent"}))
        return 0

    results = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-alleles-deep-") as directory:
        work = Path(directory)
        reference, sam, reference_text = write_inputs(work)
        bam = work / "input.bam"
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam),
                        "-O", str(bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)

        for case in case_definitions(reference_text):
            if args.case and case["name"] not in args.case:
                continue
            feature = write_feature(work, case["name"], case["rows"], root)
            allele_args = [] if args.drop_alleles else ["--alleles", str(feature)]
            common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
                      *allele_args, "--min-pruning", "1",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false"]
            java_vcf = work / f"gatk.{case['name']}.vcf"
            native_vcf = work / f"native.{case['name']}.vcf"
            java_run = subprocess.run(
                [java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                 "-O", str(java_vcf)],
                capture_output=True, text=True, timeout=300)
            native_run = subprocess.run(
                [str(native), *common, "--threads", "2", "-O", str(native_vcf)],
                capture_output=True, text=True, timeout=300)
            entry = {"case": case["name"], "why": case["why"],
                     "java_exit": java_run.returncode, "native_exit": native_run.returncode}
            if java_run.returncode != 0 or native_run.returncode != 0:
                entry.update({"status": "error",
                              "java_stderr": java_run.stderr[-400:],
                              "native_stderr": native_run.stderr[-400:]})
                results.append(entry)
                continue
            java_rows, native_rows = records(java_vcf), records(native_vcf)
            entry.update(compare(java_rows, native_rows))
            entry["match"] = "first_diff_field" not in entry
            entry["status"] = "match" if entry["match"] else "diverge"
            if not entry["match"]:
                entry["java_rows_all"] = java_rows
                entry["native_rows_all"] = native_rows
            results.append(entry)

    diverged = [entry["case"] for entry in results if entry["status"] != "match"]
    print(json.dumps({
        "status": "pass" if not diverged else "diverge",
        "release": "GATK 4.6.2.0",
        "binary": str(native),
        "cases": results,
        "diverged_cases": diverged,
    }, indent=2, sort_keys=True))
    return 0 if not diverged else 1


if __name__ == "__main__":
    raise SystemExit(main())
