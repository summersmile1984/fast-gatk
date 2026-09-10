#!/usr/bin/env python3
"""Pinned GATK oracle for DenoiseReadCounts HDF5 input metadata propagation.

SimpleCountCollection stores its SAM sequence dictionary in HDF5 metadata,
not in the interval/count matrices.  GATK copies that dictionary to both
CopyRatioCollection outputs.  This check catches a native path that computes
the right rows but emits only ``@HD``/``@RG`` after reading HDF5, which breaks
downstream interval ordering and metadata validation.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
DENOISE = pathlib.Path(os.environ.get(
    "FASTGATK_DENOISE_READ_COUNTS_BINARY", str(BUILD / "fastgatk-denoise-read-counts")))
COLLECT = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_READ_COUNTS_BINARY", str(BUILD / "fastgatk-collect-read-counts")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
BAM = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/collect-read-counts-NA12878.bam"
INTERVALS = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/collect-read-counts-test.interval_list"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True)


def header(path: pathlib.Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("@")]


def data(path: pathlib.Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@") and not line.startswith("CONTIG")]


def main() -> int:
    required = (DENOISE, COLLECT, JAVA, JAR, BAM, INTERVALS)
    if any(not path.exists() for path in required):
        missing = next(path for path in required if not path.exists())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit(f"missing DenoiseReadCounts HDF5 metadata oracle input: {missing}")
        print(json.dumps({"status": "skip", "reason": f"missing {missing}"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-denoise-hdf5-metadata-") as directory:
        work = pathlib.Path(directory)
        counts = work / "counts.hdf5"
        collect = run([
            str(COLLECT), "-I", str(BAM), "-L", str(INTERVALS), "-O", str(counts),
            "--format", "HDF5", "--sample", "NA12878",
        ])
        if collect.returncode != 0:
            raise AssertionError(f"native CollectReadCounts failed:\n{collect.stderr}")

        native_std = work / "native.standardized.tsv"
        native_den = work / "native.denoised.tsv"
        native = run([
            str(DENOISE), "-I", str(counts),
            "--standardized-copy-ratios", str(native_std),
            "--denoised-copy-ratios", str(native_den),
        ])
        if native.returncode != 0:
            raise AssertionError(f"native DenoiseReadCounts failed:\n{native.stderr}")

        java_std = work / "java.standardized.tsv"
        java_den = work / "java.denoised.tsv"
        java = run([
            str(JAVA), "-jar", str(JAR), "DenoiseReadCounts", "-I", str(counts),
            "--standardized-copy-ratios", str(java_std),
            "--denoised-copy-ratios", str(java_den),
        ])
        if java.returncode != 0:
            raise AssertionError(f"GATK DenoiseReadCounts failed:\n{java.stderr[-4000:]}")

        native_std_header = header(native_std)
        java_std_header = header(java_std)
        native_den_header = header(native_den)
        java_den_header = header(java_den)
        if native_std_header != java_std_header or native_den_header != java_den_header:
            raise AssertionError({
                "native_standardized_header": native_std_header[:3],
                "java_standardized_header": java_std_header[:3],
                "native_denoised_header": native_den_header[:3],
                "java_denoised_header": java_den_header[:3],
            })
        if not any(line.startswith("@SQ\t") and "\tAS:" in line and "\tM5:" in line
                   for line in native_std_header):
            raise AssertionError("HDF5 sequence dictionary was not propagated to DenoiseReadCounts output")
        if data(native_std) != data(java_std) or data(native_den) != data(java_den):
            raise AssertionError("native and GATK output rows differ after HDF5 metadata decoding")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "hdf5_input": True,
            "sequence_dictionary_propagated": True,
            "standardized_header_exact": True,
            "denoised_header_exact": True,
            "rows_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
