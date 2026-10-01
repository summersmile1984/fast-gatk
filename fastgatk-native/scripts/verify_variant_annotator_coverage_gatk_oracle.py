#!/usr/bin/env python3
"""VariantAnnotator oracle: -A Coverage over real reads (INFO/DP).

Runs GATK 4.6.2.0's VariantAnnotator and the native fastgatk-variant-
annotator on the same real chr17 fixture (NA12878 dictFix BAM + the b37
chr17 1Mb reference) with a small hand-written input VCF, and requires
byte-identical output (headers aside from execution provenance, data rows
verbatim).

Covered semantics:
  * -A Coverage reports likelihoods.evidenceCount() as INFO/DP — one
    evidence row per pileup element at the variant's start position
    (VariantAnnotator.makeLikelihoods builds one row per ReadPileup
    element, and Coverage.annotate emits the row count).
  * the -R header contract: contig lines rebuilt from the sequence
    dictionary with assembly = the reference file name, plus the
    ##reference=file:// URI line (VcfUtils.updateHeaderContigLines).
  * records keep their input ID/INFO verbatim otherwise.

The fixture BAM/reference live in gatk-source/src/test/resources (both
materialized in this checkout); the input VCF is written by this script.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import oracle_guard


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_VARIANT_ANNOTATOR_BINARY",
    str(ROOT / "fastgatk-native/build/fastgatk-variant-annotator")))
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"

INPUT_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878
17\t69000\t.\tT\tC\t.\t.\t.\tGT\t0/1
17\t69010\t.\tA\tG\t.\t.\t.\tGT\t0/1
17\t69040\t.\tA\tT\t.\t.\t.\tGT\t1/1
"""


def without_provenance(path: Path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith(("##GATKCommandLine=", "##fileDate="))]


def run(command: list[str], timeout: int = 600) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, capture_output=True, text=True, check=False,
                          timeout=timeout)


def main() -> int:
    if not all(path.is_file() for path in (JAVA, GATK, NATIVE, REFERENCE, BAM)):
        oracle_guard.oracle_not_verified(
            'verify_variant_annotator_coverage_gatk_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK VariantAnnotator inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK, native VariantAnnotator or chr17 fixture unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-annotator-cov-") as directory:
        work = Path(directory)
        (work / "input.vcf").write_text(INPUT_VCF, encoding="utf-8")
        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf"

        gatk = run([str(JAVA), "-Xmx512m", "-jar", str(GATK), "VariantAnnotator",
                    "-R", str(REFERENCE), "-I", str(BAM),
                    "-V", str(work / "input.vcf"), "-O", str(gatk_output),
                    "-A", "Coverage"])
        if gatk.returncode != 0:
            print(json.dumps({"status": "fail", "stage": "gatk",
                              "stderr": gatk.stderr[-2000:]}))
            return 1
        native = run([str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
                      "-V", str(work / "input.vcf"), "-O", str(native_output),
                      "-A", "Coverage"])
        if native.returncode != 0:
            print(json.dumps({"status": "fail", "stage": "native",
                              "stderr": native.stderr[-2000:]}))
            return 1

        expected = without_provenance(gatk_output)
        observed = without_provenance(native_output)
        rows = [line for line in expected if not line.startswith("#")]
        depths = [line.split("\t")[7] for line in rows]
        mismatches = [
            {"line": index, "gatk": line, "native": observed[index]
                if index < len(observed) else None}
            for index, line in enumerate(expected)
            if index >= len(observed) or observed[index] != line
        ]
        payload = {
            "status": "pass" if not mismatches and len(observed) == len(expected)
                      else "fail",
            "rows": len(rows),
            "info_dp": depths,
            "mismatches": mismatches,
            "oracle": "GATK-4.6.2.0",
            "covered": ["coverage-evidence-count", "reference-header-contract"],
        }
        print(json.dumps(payload, sort_keys=True))
        return 0 if payload["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
