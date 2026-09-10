#!/usr/bin/env python3
"""Pinned GATK oracle for CreateReadCountPanelOfNormals sample-file metadata.

The v7 HDF5 PoN contains both the original and filtered panel sample filename
arrays.  GATK writes the absolute input count-file paths (not the ``@RG SM``
sample identifier) to both arrays.  This is a small metadata boundary, but it
matters when a downstream workflow uses the PoN as a reusable artifact and
needs to identify/align its normal rows.  ``strings`` is sufficient here: the
vlen HDF5 string datasets are emitted twice per input (original + panel), and
the command line contributes one additional occurrence.
"""

from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_CREATE_PON_BINARY",
    str(BUILD / "fastgatk-create-read-count-panel-of-normals")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def write_counts(path: pathlib.Path, sample: str, values: list[int]) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:10000",
        f"@RG\tID:{sample}\tSM:{sample}",
        "CONTIG\tSTART\tEND\tCOUNT",
    ]
    for index, value in enumerate(values):
        start = index * 100 + 1
        lines.append(f"chr1\t{start}\t{start + 99}\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def hdf5_strings(path: pathlib.Path) -> list[str]:
    strings = shutil.which("strings")
    if strings is None:
        raise RuntimeError("coreutils strings is required for the HDF5 sample metadata oracle")
    return subprocess.run([strings, str(path)], check=True, text=True,
                          capture_output=True).stdout.splitlines()


def assert_sample_paths(path: pathlib.Path, inputs: list[pathlib.Path], samples: list[str]) -> None:
    values = hdf5_strings(path)
    for input_path in inputs:
        expected = str(input_path.absolute())
        # original_data/sample_filenames + panel/sample_filenames.  The
        # command-line dataset is intentionally not counted: HDF5 string
        # extraction can split long command lines at an implementation-owned
        # boundary, while both metadata arrays are short standalone values.
        if values.count(expected) < 2:
            raise AssertionError({"file": str(path), "expected": expected,
                                  "occurrences": values.count(expected)})
    for sample in samples:
        wrong = str(path.parent / sample)
        if wrong in values:
            raise AssertionError({"file": str(path), "sample_identifier_written_as_path": wrong})


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing sample-metadata oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-pon-sample-metadata-") as directory:
        work = pathlib.Path(directory)
        samples = ["RG_SAMPLE_ALPHA", "RG_SAMPLE_BETA", "RG_SAMPLE_GAMMA"]
        inputs = []
        for index, sample in enumerate(samples):
            path = work / f"normal-{index + 1}.tsv"
            write_counts(path, sample, [100 + index, 102 + index, 98 + index, 101 + index])
            inputs.append(path)

        common = [
            "--number-of-eigensamples", "0",
            "--minimum-interval-median-percentile", "0",
            "--maximum-zeros-in-sample-percentage", "100",
            "--maximum-zeros-in-interval-percentage", "100",
            "--extreme-sample-median-percentile", "0",
            "--extreme-outlier-truncation-percentile", "0",
        ]
        native_output = work / "native.pon.hdf5"
        native_command = [str(NATIVE)]
        native_command += sum((["-I", str(path)] for path in inputs), [])
        native_command += ["-O", str(native_output), *common]
        native = subprocess.run(native_command, text=True, capture_output=True)
        if native.returncode != 0:
            raise AssertionError("native CreateReadCountPanelOfNormals failed:\n" + native.stderr[-4000:])
        assert_sample_paths(native_output, inputs, samples)

        java_output = work / "java.pon.hdf5"
        java_command = [str(JAVA), "-jar", str(JAR), "CreateReadCountPanelOfNormals"]
        java_command += sum((["-I", str(path)] for path in inputs), [])
        java_command += ["-O", str(java_output), *common]
        java = subprocess.run(java_command, text=True, capture_output=True)
        if java.returncode != 0:
            raise AssertionError("Java GATK CreateReadCountPanelOfNormals failed:\n" + java.stderr[-4000:])
        assert_sample_paths(java_output, inputs, samples)

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "input_samples": len(inputs),
            "requested_eigensamples": 0,
            "native_sample_filename_alignment": True,
            "java_sample_filename_alignment": True,
            "metadata": "absolute-input-paths-in-original-and-panel-arrays",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
