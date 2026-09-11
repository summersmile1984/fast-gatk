#!/usr/bin/env python3
"""Pinned GATK oracle for VariantRecalibrator annotation-dimension ordering.

GATK reorders dimensions after normalization using the training-vs-nontraining
mean shift.  The order is part of model initialization (including JavaRandom
covariance draws), so preserving the CLI order changes multi-Gaussian models and
the serialized model schema.  This fixture makes the expected MQ,QD order
observable while keeping a one-Gaussian score comparison stable.
"""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import re
import subprocess
import tempfile
import oracle_guard

ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_VARIANT_RECALIBRATOR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-variant-recalibrator")))
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        raise AssertionError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}")
    return result


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def model_annotation_order(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    lines = opener(path, "rt", encoding="utf-8").read().splitlines()
    marker = "#:GATKTable:AnnotationMeans:"
    start = next(index for index, line in enumerate(lines) if line.startswith(marker))
    names: list[str] = []
    for line in lines[start + 2:]:
        if not line.strip():
            break
        fields = line.split()
        if fields and fields[0] not in {"Annotation", "Mean"}:
            names.append(fields[0])
    return names


def scores(path: pathlib.Path) -> list[float]:
    opener = gzip.open if path.name.endswith(".gz") else open
    values: list[float] = []
    for line in opener(path, "rt", encoding="utf-8").read().splitlines():
        if not line or line.startswith("#"):
            continue
        match = re.search(r"(?:^|;)VQSLOD=([^;]+)", line.split("\t")[7])
        if match:
            values.append(float(match.group(1)))
    return values


def main() -> int:
    if not BINARY.exists() or not JAVA.exists() or not GATK.exists():
        oracle_guard.oracle_not_verified('verify_variant_recalibrator_annotation_order_gatk_oracle.py', JAVA, GATK)
        raise SystemExit("VariantRecalibrator annotation-order oracle inputs are required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-order-oracle-") as temporary:
        work = pathlib.Path(temporary)
        rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}"
            for index, (qd, mq) in enumerate([(30, 60), (28, 58), (25, 55),
                                               (4, 25), (3, 20), (2, 15)], 1)
        ]
        input_vcf, training_vcf = work / "input.vcf", work / "training.vcf"
        write_vcf(input_vcf, rows)
        write_vcf(training_vcf, rows[:3])
        for path in (input_vcf, training_vcf):
            run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)])
        resources = ["--resource:truth,training=true,truth=true,prior=15.0", str(training_vcf)]
        common = [
            "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
            "--max-attempts", "20", "--max-iterations", "20", "--k-means-iterations", "20",
            "--bad-lod-score-cutoff", "100.0", "--dont-run-rscript", "true",
            "--create-output-variant-index", "false", "--sites-only-vcf-output", "true",
            "--truth-sensitivity-tranche", "100",
        ]
        java_output, java_model = work / "java.vcf", work / "java.model"
        run([str(JAVA), "-jar", str(GATK), "VariantRecalibrator", "-V", str(input_vcf),
             *resources, *common, "-O", str(java_output), "--tranches-file", str(work / "java.tranches"),
             "--output-model", str(java_model)])
        native_output, native_model = work / "native.vcf.gz", work / "native.model"
        native_manifest = work / "native.manifest.json"
        native = run([str(BINARY), "-V", str(input_vcf), *resources, *common,
                      "-O", str(native_output), "--tranches-file", str(work / "native.tranches"),
                      "--output-model", str(native_model), "--output-manifest", str(native_manifest)])
        java_order = model_annotation_order(java_model)
        native_order = model_annotation_order(native_model)
        assert java_order == native_order == ["MQ", "QD"], (java_order, native_order)
        java_scores, native_scores = scores(java_output), scores(native_output)
        assert len(java_scores) == len(native_scores) == len(rows)
        max_delta = max(abs(left - right) for left, right in zip(java_scores, native_scores))
        assert max_delta < 1.0e-3, (java_scores, native_scores, max_delta)
        payload = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert payload["annotations"] == ["MQ", "QD"], payload
        print(json.dumps({
            "gatk_version": "4.6.2.0", "annotation_order": native_order,
            "records": len(rows), "max_abs_vqslod_delta": max_delta, "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
