#!/usr/bin/env python3
"""Pinned GATK oracle for compensated CallCopyRatioSegments statistics.

Java's ``DoubleStream.sum`` returns the high-order accumulator minus its
negated low-order compensation word.  Returning the high word alone changes
the final standard deviation by one ulp in this fixture and flips the final
segment exactly at the calling-z boundary.  Main and IGV outputs are compared
byte-for-byte with GATK 4.6.2.0.
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
    "FASTGATK_CALL_COPY_RATIO_BINARY",
    str(BUILD / "fastgatk-call-copy-ratio-segments"),
))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"

LENGTHS = (1941171830, 320091646, 1171646238, 1)
MEAN_LOG2_COPY_RATIOS = (-2, -2, 8, 11)
CALLING_Z = "13.200364270596019"


def input_table() -> str:
    header = "@HD\tVN:1.6\n"
    header += "".join(
        f"@SQ\tSN:chr{index + 1}\tLN:{length}\n"
        for index, length in enumerate(LENGTHS)
    )
    header += (
        "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
        "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n"
    )
    return header + "".join(
        f"chr{index + 1}\t1\t{length}\t10\t{mean}\n"
        for index, (length, mean) in enumerate(zip(LENGTHS, MEAN_LOG2_COPY_RATIOS))
    )


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True, env=os.environ.copy())
    if result.returncode != 0:
        raise AssertionError({
            "command": command,
            "stdout": result.stdout,
            "stderr": result.stderr[-5000:],
        })


def calls(path: pathlib.Path) -> list[str]:
    return [
        line.split("\t")[-1]
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def main() -> int:
    if not NATIVE.exists():
        raise SystemExit(f"missing native binary: {NATIVE}")
    if not JAVA.exists() or not JAR.exists():
        raise SystemExit("missing pinned GATK 4.6.2.0 runtime")

    common_args = [
        "--neutral-segment-copy-ratio-lower-bound", "0.2",
        "--neutral-segment-copy-ratio-upper-bound", "300",
        "--outlier-neutral-segment-copy-ratio-z-score-threshold", "100",
        "--calling-copy-ratio-z-score-threshold", CALLING_Z,
    ]
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-compensated-") as directory:
        work = pathlib.Path(directory)
        input_path = work / "compensated.cr.seg"
        native_output = work / "native.called.seg"
        java_output = work / "java.called.seg"
        manifest = work / "native.json"
        input_path.write_text(input_table(), encoding="ascii")

        run([
            str(NATIVE), "-I", str(input_path), "-O", str(native_output),
            "--output-manifest", str(manifest), *common_args,
        ])
        run([
            str(JAVA), "-jar", str(JAR), "CallCopyRatioSegments",
            "-I", str(input_path), "-O", str(java_output), *common_args,
        ])

        native_legacy = native_output.with_name(native_output.stem + ".igv.seg")
        java_legacy = java_output.with_name(java_output.stem + ".igv.seg")
        if native_output.read_bytes() != java_output.read_bytes():
            raise AssertionError({
                "native": native_output.read_text(encoding="ascii"),
                "java": java_output.read_text(encoding="ascii"),
            })
        if native_legacy.read_bytes() != java_legacy.read_bytes():
            raise AssertionError({
                "native_legacy": native_legacy.read_text(encoding="ascii"),
                "java_legacy": java_legacy.read_text(encoding="ascii"),
            })
        if calls(native_output) != ["0", "0", "0", "0"]:
            raise AssertionError({"calls": calls(native_output)})

        native_manifest = json.loads(manifest.read_text(encoding="ascii"))
        if not isinstance(native_manifest["statistics_mean"], float):
            raise AssertionError(native_manifest)
        if not isinstance(native_manifest["statistics_sd"], float):
            raise AssertionError(native_manifest)

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "java_double_stream_compensation_exact": True,
        "calling_z_boundary_call": "0",
        "main_output_exact": True,
        "legacy_output_exact": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
