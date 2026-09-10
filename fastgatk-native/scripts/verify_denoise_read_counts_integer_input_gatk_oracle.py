#!/usr/bin/env python3
"""Pinned GATK oracle for DenoiseReadCounts integer COUNT decoding.

SimpleCountCollection's TSV codec uses ``DataLine.getInt`` for COUNT.  The
native reader must therefore reject decimal/scientific tokens and values
outside the Java signed-int range instead of accepting them as doubles and
silently changing fractional-coverage normalization.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_DENOISE_READ_COUNTS_BINARY", str(BUILD / "fastgatk-denoise-read-counts")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def write_counts(path: pathlib.Path, count_tokens: list[str]) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:1000",
        "@RG\tID:S\tSM:S",
        "CONTIG\tSTART\tEND\tCOUNT",
    ]
    for index, token in enumerate(count_tokens):
        start = index * 100 + 1
        lines.append(f"chr1\t{start}\t{start + 99}\t{token}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def native_command(input_path: pathlib.Path, prefix: pathlib.Path) -> list[str]:
    return [
        str(NATIVE), "-I", str(input_path), "-O", str(prefix) + ".native.tsv",
        "--standardized-copy-ratios", str(prefix) + ".native.standardized.tsv",
        "--output-manifest", str(prefix) + ".native.json",
    ]


def java_command(input_path: pathlib.Path, prefix: pathlib.Path) -> list[str]:
    return [
        str(JAVA), "-jar", str(JAR), "DenoiseReadCounts", "-I", str(input_path),
        "--standardized-copy-ratios", str(prefix) + ".java.standardized.tsv",
        "--denoised-copy-ratios", str(prefix) + ".java.denoised.tsv",
    ]


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, env=os.environ.copy())


def output_values(path: pathlib.Path) -> list[float]:
    return [
        float(line.split("\t")[-1])
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
                raise SystemExit(f"missing DenoiseReadCounts integer oracle input: {path}")
            print(json.dumps({"status": "skip", "reason": f"missing {path}"}))
            return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-denoise-integer-") as directory:
        work = pathlib.Path(directory)
        valid = work / "valid.tsv"
        write_counts(valid, ["10", "20", "30"])
        # DenoiseReadCounts has two required writer products in GATK:
        # accepting only one makes a dispatcher replacement silently omit the
        # standardized copy-ratio artifact.
        native_missing = run([
            str(NATIVE), "-I", str(valid), "-O", str(work / "missing.native.tsv")])
        java_missing = run([
            str(JAVA), "-jar", str(JAR), "DenoiseReadCounts", "-I", str(valid),
            "--denoised-copy-ratios", str(work / "missing.java.tsv")])
        if native_missing.returncode == 0 or java_missing.returncode == 0:
            raise AssertionError({
                "missing_standardized_output": {
                    "native": native_missing.returncode,
                    "java": java_missing.returncode,
                },
            })
        valid_prefix = work / "valid"
        native_valid = run(native_command(valid, valid_prefix))
        java_valid = run(java_command(valid, valid_prefix))
        if native_valid.returncode != 0:
            raise AssertionError(f"native valid integer counts failed:\n{native_valid.stderr}")
        if java_valid.returncode != 0:
            raise AssertionError(f"GATK valid integer counts failed:\n{java_valid.stderr[-4000:]}")
        native_values = output_values(valid_prefix.with_suffix(".native.tsv"))
        java_values = output_values(valid_prefix.with_suffix(".java.denoised.tsv"))
        if len(native_values) != len(java_values) or any(
                abs(left - right) > 2.0e-6 for left, right in zip(native_values, java_values)):
            raise AssertionError({"native": native_values, "java": java_values})

        rejected: dict[str, dict[str, bool]] = {}
        for name, tokens in {
            "decimal": ["10.5", "20", "30"],
            "scientific": ["1e1", "20", "30"],
            "int_overflow": ["2147483648", "20", "30"],
            "negative": ["-1", "20", "30"],
        }.items():
            input_path = work / f"{name}.tsv"
            write_counts(input_path, tokens)
            prefix = work / name
            native = run(native_command(input_path, prefix))
            java = run(java_command(input_path, prefix))
            rejected[name] = {
                "native": native.returncode != 0,
                "java": java.returncode != 0,
            }
            if native.returncode == 0 or java.returncode == 0:
                raise AssertionError({
                    "case": name,
                    "native_stderr": native.stderr,
                    "java_stderr": java.stderr[-4000:],
                })

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "valid_integer_acceptance": True,
            "native_values_match_java": True,
            "both_output_files_required": True,
            "invalid_count_rejection_exact": rejected,
            "count_semantics": "signed-int-nonnegative",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
