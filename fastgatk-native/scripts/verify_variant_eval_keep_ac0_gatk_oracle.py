#!/usr/bin/env python3
"""Pinned GATK oracle for VariantEval's --keep-ac0 CountVariants boundary.

GATK defaults to deriving polymorphism from genotypes, so genotyped AC=0
SNP/INDEL records are reference loci for CountVariants.  --keep-ac0 restores
their site-level type.  This fixture keeps both hom-ref and all-no-call AC=0
records, and compares the aggregate CountVariants GATKReport byte-for-byte.
Native-only alias checks cover the remaining optional-boolean spellings so
the pinned Java portion remains bounded and does not claim parity for the
other evaluator/stratifier tables.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


VCF = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t1\thomref-snp\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0
17\t2\tnocall-snp\tC\tT\t50\tPASS\t.\tGT\t./.\t./.
17\t3\thet-snp\tG\tA\t50\tPASS\t.\tGT\t0/1\t./.
17\t4\thomvar-snp\tT\tC\t50\tPASS\t.\tGT\t1/1\t1/1
17\t5\thomref-ins\tA\tAT\t50\tPASS\t.\tGT\t0/0\t0/0
17\t6\tnocall-ins\tC\tCT\t50\tPASS\t.\tGT\t./.\t./.
"""


def invoke(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def run_java(java: Path, jar: Path, reference: Path, source: Path,
             output: Path, extra: list[str]) -> None:
    result = invoke([
        str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
        "-eval", str(source), "-L", "17:1-6", "-O", str(output), "-no-ev", "-EV",
        "CountVariants", "-no-st", *extra,
    ])
    if result.returncode != 0:
        raise RuntimeError(f"GATK VariantEval failed: {result.stderr[-3000:]}")


def run_native(binary: Path, source: Path, output: Path, extra: list[str]) -> None:
    result = invoke([
        str(binary), "-eval", str(source), "-L", "17:1-6", "-O", str(output),
        "-no-ev", "-EV", "CountVariants", "-no-st", "--gatk-report",
        *extra,
    ])
    if result.returncode != 0:
        raise RuntimeError(f"native VariantEval failed: {result.stderr[-3000:]}")


def compare_pair(binary: Path, java: Path, jar: Path, reference: Path,
                 source: Path, work: Path, label: str, extra: list[str]) -> None:
    native = work / f"native.{label}.report"
    oracle = work / f"java.{label}.report"
    run_native(binary, source, native, extra)
    run_java(java, jar, reference, source, oracle, extra)
    if native.read_bytes() != oracle.read_bytes():
        raise AssertionError(f"VariantEval report differs for {label}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_EVAL_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-eval"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
    tabix = root / "third_party/htslib-build/htslib-src/tabix"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, java, jar, reference, bgzip, tabix)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantEval keep-ac0 oracle inputs are required")
        print('{"status":"skip","reason":"GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-eval-keep-ac0-") as directory:
        work = Path(directory)
        plain_source = work / "input.vcf"
        plain_source.write_text(VCF, encoding="utf-8")
        compressed = invoke([str(bgzip), "-f", str(plain_source)])
        if compressed.returncode != 0:
            raise RuntimeError(f"fixture bgzip failed: {compressed.stderr[-3000:]}")
        source = work / "input.vcf.gz"
        indexed = invoke([str(tabix), "-f", "-p", "vcf", str(source)])
        if indexed.returncode != 0:
            raise RuntimeError(f"fixture tabix failed: {indexed.stderr[-3000:]}")
        # Bare and explicit true are equivalent in Barclay; false restores the
        # default AC=0 filtering.  Keep short and long spellings covered.
        # Compare the default and bare keep modes to Java.  The remaining
        # spellings are checked as native aliases against those two reports,
        # keeping the pinned Java portion bounded while covering all parser
        # forms accepted by GATK.
        compare_pair(binary, java, jar, reference, source, work, "default", [])
        compare_pair(binary, java, jar, reference, source, work, "bare", ["--keep-ac0"])
        native_default = work / "native.default.report"
        native_bare = work / "native.bare.report"
        for label, extra, expected in (
            ("true", ["--keep-ac0", "true"], native_bare),
            ("short_true", ["-keep-ac0", "true"], native_bare),
            ("false", ["--keep-ac0", "false"], native_default),
        ):
            native_alias = work / f"native.{label}.report"
            result = invoke([
                str(binary), "-eval", str(source), "-L", "17:1-6",
                "-O", str(native_alias), "-no-ev", "-EV", "CountVariants",
                "-no-st", "--gatk-report", *extra,
            ])
            if result.returncode != 0:
                raise RuntimeError(f"native keep-ac0 {label} failed: {result.stderr[-3000:]}")
            if native_alias.read_bytes() != expected.read_bytes():
                raise AssertionError(f"native keep-ac0 alias differs for {label}")

    print('{"status":"pass","gatk_version":"4.6.2.0",'
          '"keep_ac0_countvariants_aggregate_report":true,'
          '"keep_ac0_native_aliases":true}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
