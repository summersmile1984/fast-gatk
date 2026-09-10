#!/usr/bin/env python3
"""GATK oracle for Mutect2's zero-LOD reference/indel ActivityProfile edge."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import pysam


REFERENCE = (
    "TTCTCTGCCAATTGAGATCACTGTGATCAGATCCTCACGCTCCAGGAAAAATGAAGTGATTGAGTAGAACGCGAAAATGTATCTGCCCCCTGGTGATTAATGTATGCGCATGGACAAAAAGCAGGTCTGCTCTTGCGTAAGCCGCCCAAGTGCGTCCTGAGTGTCACGCCCTGCTATACGTTTCATCTCTAAGGCGACCGGCACTGGAGAAAAACGAAACCTTAGTAGCAACGAATTTCCTTCCAGCCTCCAGGGTCTTTGTTGGTTCCGAGCCTCCAAAGGTTCACCCAACCATCGGTAGGCATAGTCTAAGGTTCTCCGCTCGAAACTAAAAACCGACAGTCTGAGTTCTGAGACCGTTATTTCTTACATCGTTCCTTTGTATTTAACCGTTTATTCTTAGACTACGTTAAGTGACGCCAAACCGTAGTCACTCGGGTAATAGTCTACGCTGTACGGTTTCCATAGTTCACGGCCTACGGAAGGCCGTTCCCGATGAATCCGAAATCCTATGTCTACTGGGCTTGCCCACGAAGATACTCTATTAGCCTAATGGTCGGCATTGAAACCACCTGGATGCCAAGCGATTGTACAGCATAG"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def vcf_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def igv_rows(path: Path) -> list[tuple[str, int, int, str, str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("Chromosome\t"):
            continue
        fields = line.split("\t")
        assert len(fields) == 5, fields
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]))
    return rows


def write_fixture(work: Path) -> tuple[Path, Path]:
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + REFERENCE + "\n", encoding="utf-8")
    pysam.faidx(str(reference))
    (work / "reference.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:600\n", encoding="utf-8")

    bam = work / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(REFERENCE)}],
        "RG": [{"ID": "rg1", "SM": "sample", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as output:
        # One CIGAR insertion occurs at each distinct anchor.  No anchor has
        # enough support for a positive PileupQualBuffer LOD, but the old
        # native nAlt=0 shortcut made every intervening reference pileup
        # active at --initial-tumor-lod 0 and allowed graph-only calls.
        for start in range(50, 100):
            ref = pysam.AlignedSegment()
            ref.query_name = f"ref{start}"
            ref.flag = 0
            ref.reference_id = 0
            ref.reference_start = start
            ref.mapping_quality = 60
            ref.cigarstring = "180M"
            ref.query_sequence = REFERENCE[start:start + 180]
            ref.query_qualities = pysam.qualitystring_to_array("I" * 180)
            ref.set_tag("RG", "rg1")
            output.write(ref)
            if start >= 78:
                continue
            insertion = pysam.AlignedSegment()
            insertion.query_name = f"ins{start}"
            insertion.flag = 0
            insertion.reference_id = 0
            insertion.reference_start = start
            insertion.mapping_quality = 60
            insertion.cigarstring = "70M1I110M"
            insertion.query_sequence = (
                REFERENCE[start:start + 70] + "C" + REFERENCE[start + 70:start + 180]
            )
            insertion.query_qualities = pysam.qualitystring_to_array("I" * 181)
            insertion.set_tag("RG", "rg1")
            output.write(insertion)
    pysam.index(str(bam))
    return reference, bam


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    if not all(path.is_file() for path in (java, gatk, native)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK zero-LOD ActivityProfile oracle is required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-zero-lod-indel-") as directory:
        work = Path(directory)
        reference, bam = write_fixture(work)
        gatk_default = work / "gatk-default.vcf"
        gatk_default_igv = work / "gatk-default.igv"
        native_default = work / "native-default.vcf"
        native_default_igv = work / "native-default.igv"
        common = [
            "-R", str(reference), "-I", str(bam), "--tumor-sample", "sample",
            "-L", "chr1:1-600", "--initial-tumor-lod", "0", "--tumor-lod-to-emit", "0",
            "--min-pruning", "1", "--create-output-variant-index", "false",
        ]
        gatk_run = run([
            str(java), "-jar", str(gatk), "Mutect2", *common,
            "-O", str(gatk_default), "--assembly-region-out", str(gatk_default_igv),
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = run([
            str(native), *common,
            "-O", str(native_default), "--assembly-region-out", str(native_default_igv),
        ])
        assert native_run.returncode == 0, native_run.stderr
        assert vcf_records(gatk_default) == []
        assert vcf_records(native_default) == []
        assert igv_rows(gatk_default_igv) == igv_rows(native_default_igv)
        assert igv_rows(gatk_default_igv) == [
            ("chr1", 0, 1, "end-marker", "0.00000"),
            ("chr1", 0, 300, "size=300", "-1.00000"),
            ("chr1", 300, 301, "end-marker", "0.00000"),
            ("chr1", 300, 600, "size=300", "-1.00000"),
        ]

        # --force-active bypasses only the activity gate.  The source and
        # native paths must therefore recover the same graph calls, including
        # concrete GT phasing and the associated PGT/PID/PS annotations.
        gatk_forced = work / "gatk-forced.vcf"
        gatk_forced_igv = work / "gatk-forced.igv"
        native_forced = work / "native-forced.vcf"
        native_forced_igv = work / "native-forced.igv"
        gatk_forced_run = run([
            str(java), "-jar", str(gatk), "Mutect2", *common,
            "--force-active", "true", "-O", str(gatk_forced),
            "--assembly-region-out", str(gatk_forced_igv),
        ])
        assert gatk_forced_run.returncode == 0, gatk_forced_run.stderr
        native_forced_run = run([
            str(native), *common, "--force-active=true", "-O", str(native_forced),
            "--assembly-region-out", str(native_forced_igv),
        ])
        assert native_forced_run.returncode == 0, native_forced_run.stderr
        gatk_forced_records = vcf_records(gatk_forced)
        native_forced_records = vcf_records(native_forced)
        assert gatk_forced_records == native_forced_records, {
            "gatk": gatk_forced_records,
            "native": native_forced_records,
        }
        assert len(gatk_forced_records) == 2
        assert igv_rows(gatk_forced_igv) == igv_rows(native_forced_igv)
        print(json.dumps({
            "status": "pass",
            "zero_lod_reference_profile_exact": True,
            "zero_lod_default_records": 0,
            "force_active_records_exact": 2,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
