#!/usr/bin/env python3
"""Pinned GATK oracle for VariantRecalibrator's retry/iteration controls.

GATK exposes two independent controls: ``--max-attempts`` is the number of
model-build retries (default 1), while ``--max-iterations`` is the VBEM loop
bound (default 150).  This oracle checks the pinned 4.6.2.0 help/runtime
surface and verifies that native telemetry preserves that distinction.
"""

from __future__ import annotations

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
        / "fastgatk-variant-recalibrator"),
))
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def index(path: pathlib.Path) -> None:
    result = run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)])
    if result.returncode != 0:
        raise AssertionError(f"GATK IndexFeatureFile failed for {path}:\n{result.stderr}")


def main() -> int:
    if not BINARY.exists() or not JAVA.exists() or not GATK.exists():
        oracle_guard.oracle_not_verified('verify_variant_recalibrator_attempt_iteration_gatk_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantRecalibrator attempt/iteration oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    java_help = run([str(JAVA), "-jar", str(GATK), "VariantRecalibrator", "--help"])
    if java_help.returncode != 0:
        raise AssertionError(f"GATK VariantRecalibrator help failed:\n{java_help.stderr}")
    java_help_text = java_help.stdout + java_help.stderr
    assert re.search(r"--max-attempts.*Default value:\s*1", java_help_text), java_help_text
    assert re.search(r"--max-iterations.*Default value:\s*150", java_help_text), java_help_text
    native_help = run([str(BINARY), "--help"])
    assert native_help.returncode == 0, native_help.stderr
    assert "model-build retries (default 1)" in native_help.stdout, native_help.stdout
    assert "VBEM iterations (default 150)" in native_help.stdout, native_help.stdout

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-attempts-") as temporary:
        work = pathlib.Path(temporary)
        rows = [
            f"chr1\t{position}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\n"
            for position, (qd, mq) in enumerate(((18, 45), (20, 50), (22, 55)), 1)
        ]
        input_vcf = work / "input.vcf"
        training_vcf = work / "training.vcf"
        input_vcf.write_text(HEADER + "".join(rows), encoding="utf-8")
        training_vcf.write_text(HEADER + "".join(rows), encoding="utf-8")
        index(input_vcf)
        index(training_vcf)
        common = [
            "-V", str(input_vcf),
            "--resource:truth,training=true,truth=true", str(training_vcf),
            "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
            "--max-iterations", "1", "--k-means-iterations", "1",
            "--bad-lod-score-cutoff", "100", "--dont-run-rscript", "true",
            "--create-output-variant-index", "false", "--sites-only-vcf-output", "true",
            "--truth-sensitivity-tranche", "100",
        ]
        java_result = run([
            str(JAVA), "-jar", str(GATK), "VariantRecalibrator", *common,
            "-O", str(work / "java.vcf"), "--tranches-file", str(work / "java.tranches"),
        ])
        assert java_result.returncode == 0, java_result.stderr

        def native_run(attempts: str, label: str) -> dict[str, object]:
            manifest = work / f"{label}.manifest.json"
            result = run([
                str(BINARY), *common, "--max-attempts", attempts,
                "-O", str(work / f"{label}.vcf.gz"),
                "--tranches-file", str(work / f"{label}.tranches"),
                "--output-manifest", str(manifest),
            ])
            assert result.returncode == 0, result.stderr
            return json.loads(manifest.read_text(encoding="utf-8"))

        one = native_run("1", "one")
        two = native_run("2", "two")
        assert one["max_attempts"] == 1 and one["em_iterations"] == 1, one
        assert two["max_attempts"] == 2 and two["em_iterations"] == 1, two
        print(json.dumps({
            "gatk_version": "4.6.2.0",
            "java_help_default_max_attempts": 1,
            "java_help_default_max_iterations": 150,
            "native_max_attempts": two["max_attempts"],
            "native_em_iterations": two["em_iterations"],
            "retry_iteration_controls_distinct": True,
            "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
