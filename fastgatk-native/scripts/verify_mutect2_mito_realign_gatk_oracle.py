#!/usr/bin/env python3
"""Pinned GATK oracle for Mutect2 best-haplotype realignment in a repeat.

The chrM:8877 T>C,TC locus is a high-depth homopolymer-adjacent assembly
region.  One first-mate read is within GATK's 0.2-log10 BestAllele tie window.
The source priority then selects its indel-bearing haplotype and
AlignmentUtils.createReadAlignedToRef() moves its endpoint before the T>C
site.  A wrong log-domain conversion or direct base-map projection leaves one
extra F2R1 alternate observation (416 rather than the GATK value 415).
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / \
    "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "real" / "mitomode" / "mito_shifted_8000.fasta"
READS = ROOT / "testdata" / "real" / "mitomode" / "mito.bam"
INTERVAL = "chrM:8860-8890"
TARGET = ("chrM", "8877", "T", "C,TC")
TLOD_TARGET = ("chrM", "8868", "A", "AC")
EXPECTED_F2R1 = "6,184,415"
EXPECTED_TLOD = "1.02"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def target_row(path: Path, target: tuple[str, str, str, str] = TARGET) -> list[str]:
    rows = [line.rstrip("\n").split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]
    matches = [row for row in rows if tuple(row[index] for index in (0, 1, 3, 4)) == target]
    assert len(matches) == 1, (path, target, rows)
    return matches[0]


def format_values(row: list[str]) -> dict[str, str]:
    keys = row[8].split(":")
    values = row[9].split(":")
    assert len(keys) == len(values), (keys, values)
    return dict(zip(keys, values, strict=True))


def info_values(row: list[str]) -> dict[str, str]:
    return dict(field.split("=", 1) for field in row[7].split(";") if "=" in field)


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, READS)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError(f"missing mitochondrial realignment oracle inputs: {missing}")
        print(json.dumps({"status": "skip", "reason": "mitochondrial oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-mito-realign-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        common = [
            "-R", str(REFERENCE), "-I", str(READS), "-L", INTERVAL,
            "--mitochondria-mode", "--max-reads-per-alignment-start", "75",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ]
        gatk_run = run([
            str(JAVA), "-jar", str(GATK), "Mutect2", *common,
            "--native-pair-hmm-threads", "1", "-O", str(gatk_vcf),
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
        native_run = run([
            str(NATIVE), *common, "--threads", "1", "-O", str(native_vcf),
        ])
        assert native_run.returncode == 0, native_run.stderr[-4000:]

        gatk = format_values(target_row(gatk_vcf))
        native = format_values(target_row(native_vcf))
        for key in ("AD", "DP", "F1R2", "F2R1", "FAD", "SB"):
            assert native[key] == gatk[key], (key, native[key], gatk[key])
        assert gatk["F2R1"] == EXPECTED_F2R1, gatk["F2R1"]

        gatk_tlod = info_values(target_row(gatk_vcf, TLOD_TARGET))["TLOD"]
        native_tlod = info_values(target_row(native_vcf, TLOD_TARGET))["TLOD"]
        assert native_tlod == gatk_tlod, (native_tlod, gatk_tlod)
        assert gatk_tlod == EXPECTED_TLOD, gatk_tlod

    print(json.dumps({
        "status": "pass",
        "fixture": "mito_shifted_8000_high_depth_realign",
        "target": ":".join(TARGET[:2]),
        "f2r1": EXPECTED_F2R1,
        "tlod_target": ":".join(TLOD_TARGET[:2]),
        "tlod": EXPECTED_TLOD,
        "fields_exact": ["AD", "DP", "F1R2", "F2R1", "FAD", "SB", "TLOD"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
