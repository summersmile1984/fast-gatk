#!/usr/bin/env python3
"""Pinned GATK oracle for the degenerate SVD PoN failure boundary.

GATK 4.6.2.0 refuses to write a multi-sample PoN when SVD was requested but
the standardized panel has no singular value greater than ``EPSILON``.  This
is common for accidentally duplicated/identical normals: treating the empty
SVD as a valid HDF5 panel silently disables denoising in the next stage.  The
native implementation must reject the input before creating an output file.
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
    "FASTGATK_CREATE_PON_BINARY",
    str(BUILD / "fastgatk-create-read-count-panel-of-normals")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def write_counts(path: pathlib.Path, sample: str) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:10000",
        f"@RG\tID:{sample}\tSM:{sample}",
        "CONTIG\tSTART\tEND\tCOUNT",
    ]
    # Every normal is identical.  Fractional coverage and all subsequent
    # preprocessing therefore produce a zero matrix after centering.
    for index in range(12):
        start = index * 100 + 1
        lines.append(f"chr1\t{start}\t{start + 99}\t100")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing degenerate-PoN oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-pon-degenerate-") as directory:
        work = pathlib.Path(directory)
        inputs = []
        for index in range(3):
            path = work / f"normal-{index + 1}.tsv"
            write_counts(path, f"NORMAL_{index + 1}")
            inputs.append(path)

        common = [
            *sum((["-I", str(path)] for path in inputs), []),
            "--number-of-eigensamples", "1",
            "--minimum-interval-median-percentile", "0",
            "--maximum-zeros-in-sample-percentage", "100",
            "--maximum-zeros-in-interval-percentage", "100",
            "--extreme-sample-median-percentile", "0",
            "--extreme-outlier-truncation-percentile", "0",
        ]
        native_output = work / "native.pon.hdf5"
        java_output = work / "java.pon.hdf5"
        native = subprocess.run(
            [str(NATIVE), *common, "-O", str(native_output)],
            text=True, capture_output=True,
        )
        java = subprocess.run(
            [str(JAVA), "-jar", str(JAR), "CreateReadCountPanelOfNormals",
             *common, "-O", str(java_output)],
            text=True, capture_output=True,
        )
        assert native.returncode != 0, native.stdout
        assert java.returncode != 0, java.stdout
        native_text = native.stdout + native.stderr
        java_text = java.stdout + java.stderr
        assert "singular" in native_text.lower() or "svd" in native_text.lower(), native_text[-4000:]
        assert "singular" in java_text.lower() or "svd" in java_text.lower(), java_text[-4000:]
        assert not native_output.exists(), native_output
        assert not java_output.exists(), java_output

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "input_samples": 3,
            "requested_eigensamples": 1,
            "native_rejected": True,
            "java_rejected": True,
            "output_absent_after_rejection": True,
            "boundary": "multi-sample zero-singular-value SVD",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
