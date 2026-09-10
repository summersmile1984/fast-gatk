#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for ApplyVQSR's sites-only writer boundary."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
BINARY = Path(os.environ.get(
    "FASTGATK_APPLY_VQSR_BINARY",
    str(ROOT / "fastgatk-native/build/fastgatk-apply-vqsr"),
))

INPUT_HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)
RECAL_HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def run(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"{label} failed ({result.returncode}):\n{result.stderr}")
    return result


def header_and_rows(path: Path) -> tuple[list[str], list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    headers: list[str] = []
    rows: list[str] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                headers.append(line.rstrip("\r\n"))
            elif line and not line.startswith("#"):
                rows.append(line.rstrip("\r\n"))
    return headers, rows


def semantic_rows(rows: list[str]) -> list[tuple[object, ...]]:
    result: list[tuple[object, ...]] = []
    for row in rows:
        fields = row.split("\t")
        info = {}
        if fields[7] != ".":
            for token in fields[7].split(";"):
                key, separator, value = token.partition("=")
                if separator:
                    info[key] = value
        result.append((*fields[:7], float(info["VQSLOD"]), info.get("culprit", "")))
    return result


def main() -> int:
    required = (JAVA, GATK, BINARY)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ApplyVQSR sites-only oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-sites-only-oracle-") as directory:
        work = Path(directory)
        input_vcf = work / "input.vcf"
        recal_vcf = work / "recal.vcf"
        input_vcf.write_text(
            INPUT_HEADER
            + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:20\n"
            + "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT:DP\t0/0:18\n",
            encoding="utf-8",
        )
        recal_vcf.write_text(
            RECAL_HEADER
            + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=1.0\n"
            + "chr1\t2\t.\tC\tT\t.\tPASS\tVQSLOD=-1.0\n",
            encoding="utf-8",
        )
        for path in (input_vcf, recal_vcf):
            run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)],
                "GATK IndexFeatureFile")
        common = [
            "-V", str(input_vcf), "--recal-file", str(recal_vcf),
            "--lod-score-cutoff", "0", "--mode", "SNP",
            "--sites-only-vcf-output", "true",
            "--create-output-variant-index", "false",
        ]
        java_output = work / "java.vcf.gz"
        native_output = work / "native.vcf.gz"
        native_manifest = work / "native.manifest.json"
        run([str(JAVA), "-jar", str(GATK), "ApplyVQSR", *common,
             "-O", str(java_output)], "GATK ApplyVQSR")
        native_result = run([str(BINARY), *common, "-O", str(native_output),
                             "--output-manifest", str(native_manifest)],
                            "native ApplyVQSR")
        java_header, java_rows = header_and_rows(java_output)
        native_header, native_rows = header_and_rows(native_output)
        expected_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        assert java_header == native_header == [expected_header], (java_header, native_header)
        assert semantic_rows(java_rows) == semantic_rows(native_rows) and len(native_rows) == 2, {
            "java": java_rows, "native": native_rows,
        }
        assert all(len(row.split("\t")) == 8 for row in native_rows)
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert manifest["sites_only_vcf_output"] is True
        payload = json.loads(native_result.stdout.splitlines()[-1])
        assert payload["sites_only_vcf_output"] is True
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(native_rows),
            "eight_column_shape_exact": True,
            "filter_and_vqslod_rows_exact": True,
            "manifest_telemetry_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
