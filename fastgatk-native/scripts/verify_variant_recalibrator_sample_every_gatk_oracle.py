#!/usr/bin/env python3
"""Pinned GATK oracle for VariantRecalibrator sample-every output semantics.

GATK's hidden ``--sample-every-Nth-variant`` is a traversal-level option: the
first record is retained (zero-based counter), and the recalibration table is
itself downsampled.  This catches the two easy native mistakes: selecting the
second record instead of the first, or scoring a subset while still emitting
all input records.
"""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import re
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_VARIANT_RECALIBRATOR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-variant-recalibrator")))
GATK_JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"command failed ({result.returncode}):\n{' '.join(map(str, command))}\n{result.stderr}")
    return result


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def records(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines() if line and not line.startswith("#")]


def info_map(record: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for token in record.split("\t")[7].split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            values[key] = value
        elif token:
            values[token] = ""
    return values


def main() -> int:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not GATK_JAVA.exists() or not GATK_JAR.exists():
        raise SystemExit("VariantRecalibrator sample-every oracle inputs are required")

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-sample-oracle-") as temporary:
        work = pathlib.Path(temporary)
        rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\tGT\t0/1"
            for index, (qd, mq) in enumerate([(30, 60), (28, 58), (25, 55), (4, 25), (3, 20), (2, 15)], 1)
        ]
        input_vcf, training_vcf, known_vcf = (work / name for name in ("input.vcf", "training.vcf", "known.vcf"))
        write_vcf(input_vcf, rows)
        write_vcf(training_vcf, rows[:3])
        write_vcf(known_vcf, rows[:2])
        for path in (input_vcf, training_vcf, known_vcf):
            run([str(GATK_JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(path)])

        resources = [
            "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
            "--resource:known,training=false,truth=false,known=true", str(known_vcf),
        ]
        java_common = [
            "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
            "--max-iterations", "20", "--max-attempts", "1", "--k-means-iterations", "20",
            "--bad-lod-score-cutoff", "100.0", "--dont-run-rscript", "true",
            "--create-output-variant-index", "false", "--sample-every-Nth-variant", "2",
            "--sites-only-vcf-output", "true",
        ]
        java_output = work / "java.recal.vcf"
        run([
            str(GATK_JAVA), "-jar", str(GATK_JAR), "VariantRecalibrator", "-V", str(input_vcf),
            *resources, *java_common, "-O", str(java_output),
            "--tranches-file", str(work / "java.tranches"), "--truth-sensitivity-tranche", "100",
        ])

        native_output = work / "native.recal.vcf.gz"
        native_manifest = work / "native.manifest.json"
        native_common = [
            "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
            "--max-iterations", "20", "--k-means-iterations", "20",
            "--bad-lod-score-cutoff", "100.0", "--create-output-variant-index", "false",
            "--sample-every-Nth-variant", "2", "--sites-only-vcf-output", "true",
        ]
        native_result = run([
            str(BINARY), "-V", str(input_vcf), *resources, *native_common,
            "-O", str(native_output), "--tranches-file", str(work / "native.tranches"),
            "--output-manifest", str(native_manifest),
        ])

        java_records = records(java_output)
        native_records = records(native_output)
        assert len(java_records) == len(native_records) == 3, (java_records, native_records)
        java_positions = [int(record.split("\t")[1]) for record in java_records]
        native_positions = [int(record.split("\t")[1]) for record in native_records]
        assert java_positions == native_positions == [1, 3, 5], (java_positions, native_positions)
        assert all(len(record.split("\t")) == 8 for record in java_records), java_records
        assert all(len(record.split("\t")) == 8 for record in native_records), native_records

        java_scores = []
        native_scores = []
        for java_record, native_record in zip(java_records, native_records):
            java_info = info_map(java_record)
            native_info = info_map(native_record)
            java_scores.append(float(java_info["VQSLOD"]))
            native_scores.append(float(native_info["VQSLOD"]))
            assert abs(java_scores[-1] - native_scores[-1]) < 1.0e-3, (java_record, native_record)
            for label in ("POSITIVE_TRAIN_SITE", "NEGATIVE_TRAIN_SITE"):
                assert (label in java_info) == (label in native_info), (java_record, native_record)

        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert manifest["sample_every_nth_variant"] == 2, manifest
        assert manifest["input_records"] == 6, manifest
        assert manifest["scoreable_records"] == 3, manifest
        assert manifest["scored_records"] == 3, manifest
        assert manifest["training_records"] == 2, manifest
        payload = json.loads(native_result.stdout)
        assert payload["scored_records"] == 3, payload
        print(json.dumps({
            "gatk_version": "4.6.2.0",
            "sample_every_nth_variant": 2,
            "positions": native_positions,
            "records": len(native_records),
            "max_abs_vqslod_delta": max(abs(left - right) for left, right in zip(java_scores, native_scores)),
            "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
