#!/usr/bin/env python3
"""Pinned GATK oracle for ApplyVQSR's no-target default cutoff.

ApplyVQSR accepts a tranches file, but it only enters tranche walking when
``--truth-sensitivity-filter-level`` is explicit.  With a tranches file and no
target GATK still applies its default VQSLOD >= 0 cutoff and emits
``LOW_VQSLOD``.  This small fixture prevents the native implementation from
silently selecting the first tranche as a requested target.
"""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_APPLY_VQSR_BINARY",
    ROOT / "fastgatk-native" / "build" / "fastgatk-apply-vqsr",
))
JAVA = pathlib.Path(os.environ.get("JAVA", ROOT / "third_party" / "jdk17" / "bin" / "java"))
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def records(path: pathlib.Path) -> list[tuple[str, str, str, float]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    result: list[tuple[str, str, str, float]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(
                item.split("=", 1) for item in fields[7].split(";")
                if "=" in item
            ) if fields[7] != "." else {}
            result.append((fields[0], fields[1], fields[6], float(info["VQSLOD"])))
    return result


def index(path: pathlib.Path) -> None:
    result = run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)])
    if result.returncode != 0:
        raise AssertionError(f"GATK IndexFeatureFile failed for {path}:\n{result.stderr}")


def main() -> int:
    required = (BINARY, JAVA, GATK)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ApplyVQSR default-cutoff GATK oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-default-cutoff-") as temporary:
        work = pathlib.Path(temporary)
        input_vcf = work / "input.vcf"
        recal_vcf = work / "recal.vcf"
        tranches = work / "tranches.csv"
        body = (
            "chr1\t1\t.\tA\tG\t50\tPASS\t.\n"
            "chr1\t2\t.\tC\tT\t50\tPASS\t.\n"
            "chr1\t3\t.\tG\tA\t50\tPASS\t.\n"
        )
        input_vcf.write_text(HEADER + body, encoding="utf-8")
        recal_vcf.write_text(
            HEADER
            + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.0\n"
            + "chr1\t2\t.\tC\tT\t.\tPASS\tVQSLOD=0.5\n"
            + "chr1\t3\t.\tG\tA\t.\tPASS\tVQSLOD=-1.0\n",
            encoding="utf-8",
        )
        tranches.write_text(
            "# Variant quality score tranches file\n# Version number 5\n"
            "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
            "90.00,1,1,2.0,1.0,1.0000,VQSRTrancheSNP90.00to99.00,SNP,100,90,0.9000\n"
            "99.00,1,1,2.0,1.0,0.0000,VQSRTrancheSNP99.00to100.00,SNP,100,99,0.9900\n",
            encoding="utf-8",
        )
        index(input_vcf)
        index(recal_vcf)
        common = [
            "-V", str(input_vcf), "--recal-file", str(recal_vcf),
            "--tranches-file", str(tranches), "--mode", "SNP",
            "--create-output-variant-index", "false",
        ]
        java_output = work / "java.vcf.gz"
        java_result = run([
            str(JAVA), "-jar", str(GATK), "ApplyVQSR", *common,
            "-O", str(java_output),
        ])
        if java_result.returncode != 0:
            raise AssertionError(f"GATK ApplyVQSR failed:\n{java_result.stderr}")
        native_output = work / "native.vcf.gz"
        native_result = run([str(BINARY), *common, "-O", str(native_output)])
        if native_result.returncode != 0:
            raise AssertionError(f"native ApplyVQSR failed:\n{native_result.stderr}")
        expected = records(java_output)
        actual = records(native_output)
        if actual != expected:
            raise AssertionError({"expected": expected, "actual": actual})
        expected_filters = [row[2] for row in expected]
        if expected_filters != ["PASS", "PASS", "LOW_VQSLOD"]:
            raise AssertionError(f"unexpected pinned GATK default filters: {expected_filters}")
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "tranches_file_supplied": True,
            "truth_sensitivity_target_supplied": False,
            "effective_cutoff": 0.0,
            "filter_vector_exact": True,
            "default_cutoff_semantics_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
