#!/usr/bin/env python3
"""Pinned GATK oracle for LeftAlignAndTrimVariants sites-only VCF output.

The normalization path must retain FORMAT/sample evidence while trimming and
left-aligning records, then remove it only at the final writer boundary when
``--sites-only-vcf-output`` is enabled.  This guard compares the emitted site
columns and checks the exact 8-column VCF shape against GATK 4.6.2.0.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_vcf(path: Path) -> tuple[list[str], list[list[str]]]:
    opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
    with opener as handle:
        lines = [line.rstrip("\n") for line in handle]
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line.split("\t") for line in lines if line and not line.startswith("#")]
    return header.split("\t"), records


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-left-align-trim"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not all(path.is_file() for path in (binary, java, gatk)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("LeftAlignAndTrimVariants sites-only oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "native or pinned GATK assets unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-left-align-sites-only-oracle-") as directory:
        work = Path(directory)
        reference = work / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 12 + "\n", encoding="utf-8")
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\n", encoding="utf-8")
        reference.with_suffix(".fa.fai").write_text(
            "chr1\t12\t6\t12\t13\n", encoding="utf-8")
        source = work / "input.vcf"
        source.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=12>\n"
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t4\tindel\tAA\tA\t50\tPASS\tDP=10\tGT:AD:DP\t0/1:5,5:10\n"
            "chr1\t9\tsnp\tA\tC\t60\tPASS\tDP=12\tGT:AD:DP\t0/1:6,6:12\n",
            encoding="utf-8")

        results: dict[str, dict[str, object]] = {}
        for label, value_args in {
            "true-separated": ["--sites-only-vcf-output", "true"],
            "true-bare": ["--sites-only-vcf-output"],
            "false-separated": ["--sites-only-vcf-output", "false"],
        }.items():
            native_output = work / f"{label}.native.vcf.gz"
            gatk_output = work / f"{label}.gatk.vcf.gz"
            native_manifest = work / f"{label}.manifest.json"
            native_run = run([
                str(binary), "-V", str(source), "-R", str(reference),
                "-O", str(native_output), "--output-manifest", str(native_manifest),
                *value_args])
            assert native_run.returncode == 0, f"{label}: native failed: {native_run.stderr}"
            gatk_run = run([
                str(java), "-Xmx1g", "-jar", str(gatk), "LeftAlignAndTrimVariants",
                "-V", str(source), "-R", str(reference), "-O", str(gatk_output),
                *value_args])
            assert gatk_run.returncode == 0, f"{label}: GATK failed: {gatk_run.stderr}"

            native_header, native_records = read_vcf(native_output)
            gatk_header, gatk_records = read_vcf(gatk_output)
            expected_sites_only = value_args != ["--sites-only-vcf-output", "false"]
            expected_columns = 8 if expected_sites_only else 10
            assert len(native_header) == expected_columns, (label, native_header)
            assert len(gatk_header) == expected_columns, (label, gatk_header)
            assert len(native_records) == len(gatk_records)
            for index, (native_record, gatk_record) in enumerate(zip(native_records, gatk_records)):
                assert len(native_record) == expected_columns, (label, index, native_record)
                assert len(gatk_record) == expected_columns, (label, index, gatk_record)
                assert native_record[:8] == gatk_record[:8], (
                    f"{label}/{index}: site columns differ: {native_record[:8]} != {gatk_record[:8]}")
                if not expected_sites_only:
                    assert native_record[8:] == gatk_record[8:], (
                        f"{label}/{index}: FORMAT/sample differs: {native_record[8:]} != {gatk_record[8:]}")
            manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
            assert manifest["compatibility"]["sites_only_vcf_output"] is expected_sites_only
            results[label] = {
                "columns": expected_columns,
                "records": len(native_records),
                "shape_exact": True,
            }

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "optional_boolean_forms": results,
            "sites_only_shape_exact": True,
            "format_payload_preserved_when_false": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
