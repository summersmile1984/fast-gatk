#!/usr/bin/env python3
"""Pinned GATK oracle for ApplyVQSR's filter-control optional booleans.

GATK 4.6.2.0 binds both ``--ignore-all-filters`` and
``--exclude-filtered`` as optional Booleans.  These switches are materially
observable: the former decides whether an already-filtered input is rescored,
and the latter decides whether a newly VQSR-filtered record is published.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
JAVA = Path(os.environ.get("JAVA", ROOT / "third_party/jdk17/bin/java"))
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
BINARY = Path(os.environ.get(
    "FASTGATK_APPLY_VQSR_BINARY",
    ROOT / "fastgatk-native/build/fastgatk-apply-vqsr",
))
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##FILTER=<ID=LowQual,Description=low quality>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def run(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"{label} failed ({result.returncode}):\n{result.stderr}")
    return result


def index(path: Path) -> None:
    run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)],
        f"IndexFeatureFile {path.name}")


def rows(path: Path) -> list[tuple[str, str, str, float | None]]:
    opener = gzip.open if path.suffix == ".gz" else open
    result: list[tuple[str, str, str, float | None]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            info = dict(token.split("=", 1) for token in fields[7].split(";")
                        if "=" in token)
            result.append((fields[0], fields[1], fields[6],
                           float(info["VQSLOD"]) if "VQSLOD" in info else None))
    return result


def invoke(executable: list[str], input_vcf: Path, recal_vcf: Path,
           output: Path, option: str, value: str | None) -> list[tuple[str, str, str, float | None]]:
    command = [*executable, "ApplyVQSR"] if executable[0] == str(JAVA) else list(executable)
    command += ["-V", str(input_vcf), "--recal-file", str(recal_vcf),
                "--lod-score-cutoff", "0", "--mode", "SNP",
                "--create-output-variant-index", "false", option]
    if value is not None:
        command.append(value)
    command += ["-O", str(output)]
    run(command, f"{'GATK' if executable[0] == str(JAVA) else 'native'} {option} {value}")
    return rows(output)


def main() -> int:
    required = (JAVA, GATK, BINARY)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ApplyVQSR filter boolean oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    java = [str(JAVA), "-jar", str(GATK)]
    native = [str(BINARY)]
    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-filter-bool-") as directory:
        work = Path(directory)
        filtered_input = work / "filtered.input.vcf"
        filtered_recal = work / "filtered.recal.vcf"
        filtered_input.write_text(HEADER + "chr1\t10\t.\tA\tG\t50\tLowQual\t.\n", encoding="utf-8")
        filtered_recal.write_text(HEADER + "chr1\t10\t.\tA\tG\t.\tPASS\tVQSLOD=2.0\n", encoding="utf-8")
        excluded_input = work / "excluded.input.vcf"
        excluded_recal = work / "excluded.recal.vcf"
        excluded_input.write_text(HEADER + "chr1\t20\t.\tC\tT\t50\tPASS\t.\n", encoding="utf-8")
        excluded_recal.write_text(HEADER + "chr1\t20\t.\tC\tT\t.\tPASS\tVQSLOD=-1.0\n", encoding="utf-8")
        for path in (filtered_input, filtered_recal, excluded_input, excluded_recal):
            index(path)

        for value, expected in ((None, [("chr1", "10", "PASS", 2.0)]),
                                ("true", [("chr1", "10", "PASS", 2.0)]),
                                ("false", [("chr1", "10", "LowQual", None)])):
            java_rows = invoke(java, filtered_input, filtered_recal,
                               work / f"java.ignore.{value or 'bare'}.vcf.gz",
                               "--ignore-all-filters", value)
            native_rows = invoke(native, filtered_input, filtered_recal,
                                 work / f"native.ignore.{value or 'bare'}.vcf.gz",
                                 "--ignore-all-filters", value)
            if java_rows != native_rows or java_rows != expected:
                raise AssertionError({"option": "ignore-all-filters", "value": value,
                                      "java": java_rows, "native": native_rows,
                                      "expected": expected})

        for value, expected in ((None, []), ("true", []),
                                ("false", [("chr1", "20", "LOW_VQSLOD", -1.0)])):
            java_rows = invoke(java, excluded_input, excluded_recal,
                               work / f"java.exclude.{value or 'bare'}.vcf.gz",
                               "--exclude-filtered", value)
            native_rows = invoke(native, excluded_input, excluded_recal,
                                 work / f"native.exclude.{value or 'bare'}.vcf.gz",
                                 "--exclude-filtered", value)
            if java_rows != native_rows or java_rows != expected:
                raise AssertionError({"option": "exclude-filtered", "value": value,
                                      "java": java_rows, "native": native_rows,
                                      "expected": expected})

        print(json.dumps({
            "status": "pass", "gatk_version": "4.6.2.0",
            "ignore_all_filters_bare_true_false_exact": True,
            "exclude_filtered_bare_true_false_exact": True,
            "filter_and_vqslod_rows_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
