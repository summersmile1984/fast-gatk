#!/usr/bin/env python3
"""Pinned GATK oracle for the Java multi-Gaussian VBEM model path.

The two-component fixture exercises the Normal-Wishart shrinkage outer product,
random initialization, convergence, and final-model evaluation.  Component IDs
are not stable model identities, so serialized parameters are matched by their
mean vectors before comparison.  GATKReport whitespace/number formatting is a
separate boundary; numeric model parameters and emitted VQSLOD are strict.
"""

from __future__ import annotations

import json
import gzip
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


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def vqslod_texts(path: pathlib.Path) -> list[str]:
    values: list[str] = []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    for line in lines:
        if not line or line.startswith("#"):
            continue
        info = line.split("\t")[7]
        match = re.search(r"(?:^|;)VQSLOD=([^;]+)", info)
        if match:
            values.append(match.group(1))
    return values


def report_table(path: pathlib.Path, name: str) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    marker = next(index for index, line in enumerate(lines)
                  if line.startswith(f"#:GATKTable:{name}:"))
    header: list[str] = []
    rows: list[list[str]] = []
    for line in lines[marker + 1:]:
        if line.startswith("#:"):
            break
        fields = line.split()
        if not fields:
            if rows:
                break
            continue
        if not header:
            header = fields
        else:
            rows.append(fields)
    if not header or not rows:
        raise AssertionError(f"missing GATKReport table {name}: {path}")
    return header, rows


def model_components(path: pathlib.Path, positive: bool) -> list[tuple[float, ...]]:
    prefix = "Positive" if positive else "Negative"
    pmix_name = "GoodGaussianPMix" if positive else "BadGaussianPMix"
    mean_header, mean_rows = report_table(path, f"{prefix}ModelMeans")
    covariance_header, covariance_rows = report_table(path, f"{prefix}ModelCovariances")
    _, pmix_rows = report_table(path, pmix_name)
    dimensions = mean_header[1:]
    means = {int(row[0]): tuple(float(value) for value in row[1:]) for row in mean_rows}
    pmix = {int(row[0]): float(row[1]) for row in pmix_rows}
    covariance: dict[int, dict[str, tuple[float, ...]]] = {}
    for row in covariance_rows:
        covariance.setdefault(int(row[0]), {})[row[1]] = tuple(float(value) for value in row[2:])
    components: list[tuple[float, ...]] = []
    for gaussian, mean in means.items():
        matrix = tuple(value for annotation in dimensions
                       for value in covariance[gaussian][annotation])
        components.append((*mean, *matrix, pmix[gaussian]))
    return sorted(components, key=lambda component: component[:len(dimensions)])


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"command failed ({result.returncode}):\n{' '.join(command)}\n{result.stderr}")
    return result


def main() -> int:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not GATK_JAVA.exists() or not GATK_JAR.exists():
        raise SystemExit("VariantRecalibrator VBEM oracle inputs are required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-vbem-oracle-") as temporary:
        work = pathlib.Path(temporary)
        rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\tGT\t0/1"
            for index, (qd, mq) in enumerate([
                (30, 60), (28, 58), (25, 55), (24, 54), (23, 53), (22, 52),
            ], 1)
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
        common = [
            "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "2",
            "--max-attempts", "20", "--k-means-iterations", "20",
            "--bad-lod-score-cutoff", "100.0", "--create-output-variant-index", "false",
        ]
        gatk_output = work / "gatk.recal.vcf"
        run([
            str(GATK_JAVA), "-jar", str(GATK_JAR), "VariantRecalibrator", "-V", str(input_vcf),
            *resources, *common, "--dont-run-rscript", "true", "--sites-only-vcf-output", "true",
            "-O", str(gatk_output), "--tranches-file", str(work / "gatk.tranches"),
            "--truth-sensitivity-tranche", "100", "--output-model", str(work / "gatk.model"),
        ])
        native_output = work / "native.recal.vcf.gz"
        native = run([
            str(BINARY), "-V", str(input_vcf), *resources, *common, "-O", str(native_output),
            "--tranches-file", str(work / "native.tranches"), "--output-model", str(work / "native.model"),
        ])
        java_tokens = vqslod_texts(gatk_output)
        native_tokens = vqslod_texts(native_output)
        assert len(java_tokens) == len(rows), java_tokens
        assert len(native_tokens) == len(rows), native_tokens
        # Raw-bit gate through serialization: GATK writes the VQSLOD INFO
        # attribute as String.format("%.4f", datum.lod) of the full-precision
        # double (VariantDataManager); native now emits the identical %.4f
        # text, so exact token equality is required.  This replaces the old
        # 1e-3 bound that the historical squared-shrinkage bug also tripped
        # (maximum delta 1.22744 on this fixture).
        assert java_tokens == native_tokens, {
            "java": java_tokens, "native": native_tokens,
        }
        deltas = [abs(float(left) - float(right))
                  for left, right in zip(java_tokens, native_tokens)]
        model_deltas: list[float] = []
        for positive in (True, False):
            java_components = model_components(work / "gatk.model", positive)
            native_components = model_components(work / "native.model", positive)
            assert len(java_components) == len(native_components) == 2, (
                java_components, native_components,
            )
            for java_component, native_component in zip(java_components, native_components):
                assert len(java_component) == len(native_component), (
                    java_component, native_component,
                )
                model_deltas.extend(abs(left - right)
                                    for left, right in zip(java_component, native_component))
        assert max(model_deltas) < 1.0e-6, max(model_deltas)
        payload = json.loads(native.stdout)
        assert payload["components"] == 2, payload
        assert payload["variational_normal_wishart"] is True, payload
        print(json.dumps({
            "gatk_version": "4.6.2.0",
            "records": len(rows),
            "max_abs_vqslod_delta": max(deltas),
            "max_abs_model_parameter_delta": max(model_deltas),
            "gaussians": payload["components"],
            "native_model": payload["model"],
            "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
