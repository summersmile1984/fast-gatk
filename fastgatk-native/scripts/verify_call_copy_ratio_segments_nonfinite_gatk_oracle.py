#!/usr/bin/env python3
"""Pinned GATK oracle for non-finite copy-ratio segment values.

``CopyRatioSegment`` accepts NaN and +/-Infinity in the TSV decoder.  The
SimpleCopyRatioCaller then carries those IEEE-754 values through ``pow`` and
the degenerate statistics path, yielding a neutral call while preserving the
original value in both output tables.  This is an observable input/output
boundary, not a relaxed assertion: native and GATK output bytes must match.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_CALL_COPY_RATIO_BINARY",
    str(BUILD / "fastgatk-call-copy-ratio-segments"),
))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


HEADER = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
    "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n"
)


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, text=True, capture_output=True, env=os.environ.copy())
    if result.returncode != 0:
        raise AssertionError({"command": cmd, "stdout": result.stdout, "stderr": result.stderr[-4000:]})
    return result


def main() -> int:
    if not NATIVE.exists():
        raise SystemExit(f"missing native binary: {NATIVE}")
    if not JAVA.exists() or not JAR.exists():
        oracle_guard.oracle_not_verified('verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py', JAVA, JAR)
        raise SystemExit("missing pinned GATK 4.6.2.0 runtime")

    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-nonfinite-") as directory:
        work = pathlib.Path(directory)
        cases = ("NaN", "Infinity", "-Infinity")
        for index, value in enumerate(cases):
            input_path = work / f"case-{index}.cr.seg"
            native_output = work / f"case-{index}.native.seg"
            java_output = work / f"case-{index}.java.seg"
            input_path.write_text(
                HEADER
                + f"chr1\t1\t100\t10\t{value}\n"
                + "chr1\t101\t200\t10\t0.000000\n",
                encoding="ascii",
            )
            native_manifest = native_output.with_suffix(".json")
            run([
                str(NATIVE), "-I", str(input_path), "-O", str(native_output),
                "--output-manifest", str(native_manifest),
            ])
            run([
                str(JAVA), "-jar", str(JAR), "CallCopyRatioSegments",
                "-I", str(input_path), "-O", str(java_output),
            ])
            native_legacy = native_output.with_name(native_output.stem + ".igv.seg")
            java_legacy = java_output.with_name(java_output.stem + ".igv.seg")
            if native_output.read_bytes() != java_output.read_bytes():
                raise AssertionError({
                    "case": value,
                    "native": native_output.read_text(encoding="ascii"),
                    "java": java_output.read_text(encoding="ascii"),
                })
            if native_legacy.read_bytes() != java_legacy.read_bytes():
                raise AssertionError({
                    "case": value,
                    "native_legacy": native_legacy.read_text(encoding="ascii"),
                    "java_legacy": java_legacy.read_text(encoding="ascii"),
                })
            manifest = json.loads(native_manifest.read_text(encoding="ascii"))
            if manifest["statistics_mean"] is not None or manifest["statistics_sd"] is not None:
                raise AssertionError({"case": value, "manifest": manifest})

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "cases": list(cases),
        "main_output_exact": True,
        "legacy_output_exact": True,
        "ieee754_values_preserved": True,
        "degenerate_calls_neutral": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
