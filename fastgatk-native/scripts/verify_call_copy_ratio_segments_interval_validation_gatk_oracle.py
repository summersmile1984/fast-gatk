#!/usr/bin/env python3
"""Pinned GATK oracle for CallCopyRatioSegments interval validation."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_CALL_COPY_RATIO_BINARY",
    str(ROOT / "fastgatk-native/build/fastgatk-call-copy-ratio-segments"),
))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
HEADER = ("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
          "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
          "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n")


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True)


def main() -> int:
    if not NATIVE.exists() or not JAVA.exists() or not JAR.exists():
        raise SystemExit("missing native binary or pinned GATK runtime")
    cases = {
        "overlap": "chr1\t1\t100\t10\t0\nchr1\t100\t200\t10\t1\n",
        "unsorted": "chr1\t101\t200\t10\t0\nchr1\t1\t100\t10\t1\n",
    }
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-interval-") as directory:
        work = pathlib.Path(directory)
        for label, rows in cases.items():
            input_path = work / f"{label}.tsv"
            input_path.write_text(HEADER + rows, encoding="ascii")
            native = run([str(NATIVE), "-I", str(input_path), "-O", str(work / f"{label}.native.tsv")])
            java = run([str(JAVA), "-jar", str(JAR), "CallCopyRatioSegments",
                        "-I", str(input_path), "-O", str(work / f"{label}.java.tsv")])
            if native.returncode == 0 or java.returncode == 0:
                raise AssertionError({"case": label, "native": native.returncode, "java": java.returncode})
    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "overlapping_and_unsorted_partitions_rejected": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
