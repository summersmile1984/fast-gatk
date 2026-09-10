#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for the bounded VariantEval report slice.

The broad contract intentionally checks native evaluator counters.  This
oracle is narrower: CountVariants (aggregate and Sample stratification) is
compared byte-for-byte, including GATKReport format metadata, Java numeric
formatting, and column padding.  A malformed/truncated VCF is also required
to fail closed instead of being mistaken for a clean EOF.
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
17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
17\t2\t.\tC\tT\t50\tPASS\t.\tGT\t./.\t1/1
"""


def invoke(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def gatk_table_rows(report: Path, table: str) -> tuple[list[str], list[list[str]]]:
    """Read one whitespace-delimited GATKReport table."""
    lines = report.read_text(encoding="utf-8").splitlines()
    marker = next(index for index, line in enumerate(lines)
                  if line.startswith("#:GATKTable:") and f":{table}:" in line)
    header = lines[marker + 1].split()
    rows: list[list[str]] = []
    index = marker + 2
    while index < len(lines) and lines[index] and not lines[index].startswith("#:"):
        rows.append(lines[index].split())
        index += 1
    return header, rows


def compare_count_fields(native: Path, java_report: Path, fields: tuple[str, ...]) -> None:
    native_header, native_rows = gatk_table_rows(native, "CountVariants")
    java_header, java_rows = gatk_table_rows(java_report, "CountVariants")
    if native_header != java_header or len(native_rows) != 1 or len(java_rows) != 1:
        raise AssertionError("CountVariants genotype-state report shape differs")
    for field in fields:
        left = native_rows[0][native_header.index(field)]
        right = java_rows[0][java_header.index(field)]
        if left != right:
            raise AssertionError(f"CountVariants {field} differs: native={left} java={right}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_EVAL_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-eval"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, java, jar, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantEval Java oracle inputs are required")
        print('{"status":"skip","reason":"GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-eval-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")
        common = ["-eval", str(source), "-O"]

        native_aggregate = work / "native.aggregate.report"
        java_aggregate = work / "java.aggregate.report"
        native = invoke([str(binary), *common, str(native_aggregate),
                         "-no-ev", "-EV", "CountVariants", "-no-st", "--gatk-report"])
        if native.returncode != 0:
            raise RuntimeError(f"native aggregate failed: {native.stderr[-3000:]}")
        oracle = invoke([str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
                         "-eval", str(source), "-O", str(java_aggregate),
                         "-no-ev", "-EV", "CountVariants", "-no-st"])
        if oracle.returncode != 0:
            raise RuntimeError(f"GATK aggregate failed: {oracle.stderr[-3000:]}")
        if native_aggregate.read_bytes() != java_aggregate.read_bytes():
            raise AssertionError("aggregate CountVariants GATKReport differs byte-for-byte")

        native_sample = work / "native.sample.report"
        java_sample = work / "java.sample.report"
        native = invoke([str(binary), *common, str(native_sample),
                         "-no-ev", "-EV", "CountVariants", "-no-st", "-ST", "Sample",
                         "--gatk-report"])
        if native.returncode != 0:
            raise RuntimeError(f"native Sample stratification failed: {native.stderr[-3000:]}")
        oracle = invoke([str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
                         "-eval", str(source), "-O", str(java_sample),
                         "-no-ev", "-EV", "CountVariants", "-no-st", "-ST", "Sample"])
        if oracle.returncode != 0:
            raise RuntimeError(f"GATK Sample stratification failed: {oracle.stderr[-3000:]}")
        if native_sample.read_bytes() != java_sample.read_bytes():
            raise AssertionError("Sample CountVariants GATKReport differs byte-for-byte")

        # Keep the evaluator side of this bounded slice honest as well.  The
        # AF/Theta/TiTv tables are small, stable GATK evaluators; compare their
        # complete report bytes (metadata, fields, formatting, and values).
        for evaluator in ("VariantAFEvaluator", "ThetaVariantEvaluator",
                           "TiTvVariantEvaluator"):
            native_path = work / f"native.{evaluator}.report"
            java_path = work / f"java.{evaluator}.report"
            native = invoke([str(binary), "-eval", str(source), "-O", str(native_path),
                             "-no-ev", "-EV", evaluator, "-no-st", "--gatk-report"])
            if native.returncode != 0:
                raise RuntimeError(f"native {evaluator} failed: {native.stderr[-3000:]}")
            oracle = invoke([str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
                             "-eval", str(source), "-O", str(java_path),
                             "-no-ev", "-EV", evaluator, "-no-st"])
            if oracle.returncode != 0:
                raise RuntimeError(f"GATK {evaluator} failed: {oracle.stderr[-3000:]}")
            if native_path.read_bytes() != java_path.read_bytes():
                raise AssertionError(f"{evaluator} GATKReport differs byte-for-byte")

        # CountVariants must derive the locus state from genotypes when GT is
        # present: hom-ref and all-no-call ALT records are reference loci,
        # while a called ALT genotype makes the locus variant.  Sites-only
        # records retain their site-level variant state and still count as
        # called loci.  This catches the easy-to-miss distinction between
        # VariantContext ALT alleles and GATK's ignore-AC0 behavior.
        genotype_states = work / "genotype-states.vcf"
        genotype_states.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t1\tallhomref\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0
17\t2\tallnocall\tC\tT\t50\tPASS\t.\tGT\t./.\t./.
17\t3\thet\tG\tA\t50\tPASS\t.\tGT\t0/1\t./.
17\t4\thomvar\tT\tC\t50\tPASS\t.\tGT\t1/1\t1/1
17\t5\tindelhomref\tA\tAT\t50\tPASS\t.\tGT\t0/0\t0/0
17\t6\tindelhet\tC\tCT\t50\tPASS\t.\tGT\t0/1\t./.
""",
            encoding="utf-8",
        )
        native_states = work / "native.genotype-states.report"
        java_states = work / "java.genotype-states.report"
        native = invoke([str(binary), "-eval", str(genotype_states),
                         "-O", str(native_states), "-no-ev", "-EV", "CountVariants",
                         "-no-st", "--gatk-report"])
        if native.returncode != 0:
            raise RuntimeError(f"native genotype-state CountVariants failed: {native.stderr[-3000:]}")
        oracle = invoke([str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
                         "-eval", str(genotype_states), "-O", str(java_states),
                         "-no-ev", "-EV", "CountVariants", "-no-st"])
        if oracle.returncode != 0:
            raise RuntimeError(f"GATK genotype-state CountVariants failed: {oracle.stderr[-3000:]}")
        compare_count_fields(
            native_states, java_states,
            ("nProcessedLoci", "nCalledLoci", "nRefLoci", "nVariantLoci",
             "nSNPs", "nInsertions", "nNoCalls", "nHets", "nHomRef",
             "nHomVar", "nSingletons", "variantRate", "variantRatePerBp",
             "heterozygosity", "heterozygosityPerBp", "hetHomRatio",
             "indelRate", "indelRatePerBp", "insertionDeletionRatio"),
        )

        sites_only = work / "sites-only.vcf"
        sites_only.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t1\tsite-snp\tA\tG\t50\tPASS\t.
17\t2\tsite-indel\tC\tCT\t50\tPASS\t.
17\t3\tsite-ref\tG\t.\t50\tPASS\t.
""",
            encoding="utf-8",
        )
        native_sites = work / "native.sites-only.report"
        java_sites = work / "java.sites-only.report"
        native = invoke([str(binary), "-eval", str(sites_only), "-O", str(native_sites),
                         "-no-ev", "-EV", "CountVariants", "-no-st", "--gatk-report"])
        if native.returncode != 0:
            raise RuntimeError(f"native sites-only CountVariants failed: {native.stderr[-3000:]}")
        oracle = invoke([str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
                         "-eval", str(sites_only), "-O", str(java_sites),
                         "-no-ev", "-EV", "CountVariants", "-no-st"])
        if oracle.returncode != 0:
            raise RuntimeError(f"GATK sites-only CountVariants failed: {oracle.stderr[-3000:]}")
        compare_count_fields(
            native_sites, java_sites,
            ("nProcessedLoci", "nCalledLoci", "nRefLoci", "nVariantLoci",
             "nSNPs", "nInsertions", "nNoCalls", "nHets", "nHomRef",
             "nHomVar", "nSingletons", "variantRate", "variantRatePerBp",
             "heterozygosity", "heterozygosityPerBp", "hetHomRatio",
             "indelRate", "indelRatePerBp", "insertionDeletionRatio"),
        )

        malformed = work / "truncated.vcf"
        malformed.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=10>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t1\tbroken\tA\n", encoding="utf-8")
        bad_output = work / "must-not-be-written.report"
        result = invoke([str(binary), "-eval", str(malformed), "-O", str(bad_output)])
        if result.returncode == 0 or bad_output.exists():
            raise AssertionError("truncated VCF was not rejected fail-closed")

    print('{"status":"pass","gatk_version":"4.6.2.0",'
          '"count_variants_report_bytes":true,"sample_stratification_bytes":true,'
          '"evaluator_report_bytes":{"VariantAFEvaluator":true,'
          '"ThetaVariantEvaluator":true,"TiTvVariantEvaluator":true},'
          '"genotype_aware_count_variants":true,"sites_only_count_variants":true,'
          '"truncated_input_fail_closed":true}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
