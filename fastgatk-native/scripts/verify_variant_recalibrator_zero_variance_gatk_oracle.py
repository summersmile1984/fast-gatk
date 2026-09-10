#!/usr/bin/env python3
"""Pinned GATK oracle for VariantRecalibrator's zero-variance guard.

GATK's VariantDataManager rejects a training annotation whose standard
deviation is below 1e-5.  The native path previously replaced that scale with
1.0 and continued, producing a model from an input for which Java fails.  The
fixture keeps both annotations constant and requires both implementations to
fail at the same normalization boundary.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import tempfile


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


def main() -> int:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not JAVA.exists() or not GATK.exists():
        raise SystemExit("VariantRecalibrator zero-variance oracle inputs are required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-zero-variance-") as temporary:
        work = pathlib.Path(temporary)
        rows = "".join(
            f"chr1\t{position}\t.\tA\tG\t50\tPASS\tQD=10;MQ=60\n"
            for position in (1, 2, 3)
        )
        input_vcf = work / "input.vcf"
        training_vcf = work / "training.vcf"
        input_vcf.write_text(HEADER + rows, encoding="utf-8")
        training_vcf.write_text(HEADER + rows, encoding="utf-8")
        for path in (input_vcf, training_vcf):
            indexed = subprocess.run(
                [str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)],
                text=True, capture_output=True, check=False,
            )
            if indexed.returncode != 0:
                raise AssertionError(f"GATK IndexFeatureFile failed:\n{indexed.stderr}")
        common = [
            "-V", str(input_vcf),
            "--resource:truth,training=true,truth=true", str(training_vcf),
            "-an", "QD", "-an", "MQ", "--mode", "SNP",
            "--max-gaussians", "1", "--max-attempts", "5",
            "--k-means-iterations", "2", "--create-output-variant-index", "false",
        ]
        java_result = subprocess.run(
            [str(JAVA), "-jar", str(GATK), "VariantRecalibrator", *common,
             "--dont-run-rscript", "true", "-O", str(work / "java.vcf"),
             "--tranches-file", str(work / "java.tranches")],
            text=True, capture_output=True, check=False,
        )
        if java_result.returncode == 0 or "zero variance" not in java_result.stderr.lower():
            raise AssertionError(
                f"GATK did not expose the zero-variance guard (rc={java_result.returncode}):\n"
                f"{java_result.stderr}")
        native_result = subprocess.run(
            [str(BINARY), *common, "-O", str(work / "native.vcf"),
             "--tranches-file", str(work / "native.tranches")],
            text=True, capture_output=True, check=False,
        )
        if native_result.returncode == 0 or "zero variance" not in native_result.stderr.lower():
            raise AssertionError(
                f"native did not expose the zero-variance guard (rc={native_result.returncode}):\n"
                f"{native_result.stderr}")
        print('{"gatk_version":"4.6.2.0","java_returncode":%d,"native_returncode":%d,'
              '"zero_variance_guard":true,"status":"pass"}' %
              (java_result.returncode, native_result.returncode))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
