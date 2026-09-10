#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for VariantEval ValidationReport.

The fixtures pin NO_CALL/MONO/POLY status, default upstream exclusion of
filtered records, comparison-only traversal, INFO/AC sites-only
classification, and the GATK rule that subsets a genotyped comparison track
to evaluation samples.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


GENOTYPE_EVAL = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FILTER=<ID=LowQual,Description=Low quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
17\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/0
17\t3\t.\tG\tA\t50\tLowQual\t.\tGT\t0/1
17\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/1
17\t6\t.\tC\tG\t50\tPASS\t.\tGT\t0/0
17\t8\t.\tG\tT\t50\tPASS\t.\tGT\t0/1
"""

GENOTYPE_COMP = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FILTER=<ID=LowQual,Description=Low quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tEXTRA
17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
17\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0
17\t3\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\t0/0
17\t4\t.\tT\tC\t50\tPASS\t.\tGT\t0/1\t0/0
17\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\t0/1
17\t6\t.\tC\tG\t50\tPASS\t.\tGT\t0/0\t0/1
17\t7\t.\tT\tG\t50\tPASS\t.\tGT\t0/0\t0/1
17\t8\t.\tG\tT\t50\tLowQual\t.\tGT\t0/1\t0/1
"""

SITES_EVAL = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FILTER=<ID=LowQual,Description=Low quality>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t1\t.\tA\tG\t50\tPASS\tAC=1
17\t2\t.\tC\tT\t50\tPASS\tAC=0
17\t3\t.\tG\tA\t50\tLowQual\tAC=1
17\t5\t.\tA\tC\t50\tPASS\tAC=1
17\t6\t.\tC\tG\t50\tPASS\tAC=0
17\t8\t.\tG\tT\t50\tPASS\tAC=1
"""

SITES_COMP = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FILTER=<ID=LowQual,Description=Low quality>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t1\t.\tA\tG\t50\tPASS\tAC=1
17\t2\t.\tC\tT\t50\tPASS\tAC=1
17\t3\t.\tG\tA\t50\tPASS\tAC=1
17\t4\t.\tT\tC\t50\tPASS\tAC=1
17\t5\t.\tA\tC\t50\tPASS\tAC=0
17\t6\t.\tC\tG\t50\tPASS\tAC=0
17\t7\t.\tT\tG\t50\tPASS\tAC=0
17\t8\t.\tG\tT\t50\tLowQual\tAC=1
"""


def invoke(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def run_case(
    *, name: str, eval_text: str, comp_text: str, work: Path,
    binary: Path, java: Path, jar: Path, reference: Path,
) -> dict[str, object]:
    eval_path = work / f"{name}.eval.vcf"
    comp_path = work / f"{name}.comp.vcf"
    native_report = work / f"{name}.native.report"
    java_report = work / f"{name}.java.report"
    manifest = work / f"{name}.manifest.json"
    eval_path.write_text(eval_text, encoding="utf-8")
    comp_path.write_text(comp_text, encoding="utf-8")

    native = invoke([
        str(binary), "-eval", str(eval_path), "-comp", str(comp_path),
        "-O", str(native_report), "-no-st", "-no-ev", "-EV", "ValidationReport",
        "--gatk-report", "--output-manifest", str(manifest),
    ])
    if native.returncode != 0:
        raise RuntimeError(f"native {name} failed: {native.stderr[-3000:]}")
    oracle = invoke([
        str(java), "-jar", str(jar), "VariantEval", "-R", str(reference),
        "-eval", str(eval_path), "-comp", str(comp_path), "-O", str(java_report),
        "-no-st", "-no-ev", "-EV", "ValidationReport",
    ])
    if oracle.returncode != 0:
        raise RuntimeError(f"GATK {name} failed: {oracle.stderr[-3000:]}")
    if native_report.read_bytes() != java_report.read_bytes():
        raise AssertionError(
            f"{name} ValidationReport differs byte-for-byte\n"
            f"native:\n{native_report.read_text(encoding='utf-8')}\n"
            f"java:\n{java_report.read_text(encoding='utf-8')}")

    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    assert metadata["compatibility"]["validation_report"] is True
    assert metadata["fallback"]["unsupported_modules"] is False
    assert metadata["telemetry"]["validation_n_comp"] == 7
    assert metadata["telemetry"]["validation_tp"] == 1
    assert metadata["telemetry"]["validation_fp"] == 1
    assert metadata["telemetry"]["validation_fn"] == 3
    assert metadata["telemetry"]["validation_tn"] == 2
    assert metadata["telemetry"]["validation_report_execution_space"] in {
        "OpenMP", "Serial", "Cuda", "HIP", "SYCL"
    }
    return {
        "report_bytes": True,
        "execution_space": metadata["telemetry"]["validation_report_execution_space"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_EVAL_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-eval"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.is_file() for path in (binary, java, jar, reference)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantEval ValidationReport oracle inputs are required")
        print('{"status":"skip","reason":"GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-validation-report-oracle-") as directory:
        work = Path(directory)
        cases = {
            "genotypes": run_case(
                name="genotypes", eval_text=GENOTYPE_EVAL, comp_text=GENOTYPE_COMP,
                work=work, binary=binary, java=java, jar=jar, reference=reference),
            "sites_ac": run_case(
                name="sites-ac", eval_text=SITES_EVAL, comp_text=SITES_COMP,
                work=work, binary=binary, java=java, jar=jar, reference=reference),
        }
        unsafe_report = work / "must-not-write-stratified.report"
        unsafe = invoke([
            str(binary), "-eval", str(work / "genotypes.eval.vcf"),
            "-comp", str(work / "genotypes.comp.vcf"), "-O", str(unsafe_report),
            "-no-ev", "-EV", "ValidationReport", "--gatk-report",
        ])
        if unsafe.returncode == 0 or unsafe_report.exists() or "requires -no-st" not in unsafe.stderr:
            raise AssertionError("ValidationReport did not fail closed on standard stratifications")

        duplicate_eval = work / "duplicate.eval.vcf"
        duplicate_eval.write_text(
            GENOTYPE_EVAL.replace(
                "17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
                "17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
                "17\t1\tduplicate\tA\tG\t50\tPASS\t.\tGT\t0/1\n"),
            encoding="utf-8",
        )
        duplicate_report = work / "must-not-write-duplicate.report"
        duplicate = invoke([
            str(binary), "-eval", str(duplicate_eval),
            "-comp", str(work / "genotypes.comp.vcf"), "-O", str(duplicate_report),
            "-no-st", "-no-ev", "-EV", "ValidationReport", "--gatk-report",
        ])
        if (duplicate.returncode == 0 or duplicate_report.exists() or
                "at most one eval record per locus" not in duplicate.stderr):
            raise AssertionError("ValidationReport did not fail closed on duplicate eval loci")

    print(json.dumps({
        "status": "pass", "gatk_version": "4.6.2.0",
        "validation_report": cases,
        "site_states": ["NO_CALL", "MONO", "POLY"],
        "filtered_records_excluded_upstream": True,
        "comparison_sample_subset": True,
        "unsafe_stratification_fail_closed": True,
        "duplicate_eval_locus_fail_closed": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
