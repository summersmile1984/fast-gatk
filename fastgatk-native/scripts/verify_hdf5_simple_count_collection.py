#!/usr/bin/env python3
"""Audit the native SimpleCountCollection HDF5 metadata against GATK 4.6.2.0.

This is intentionally a narrow HDF5 boundary test.  It does not claim that
TSV or this small subset is a complete CNV implementation: the oracle checks
that sample metadata, the full SAM sequence dictionary, interval metadata and
counts survive a Java reader, while PoN/SVD model state remains covered by the
separate CreateReadCountPanelOfNormals contract.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_READ_COUNTS_BINARY", str(BUILD / "fastgatk-collect-read-counts")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
BAM = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/collect-read-counts-NA12878.bam"
INTERVALS = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/collect-read-counts-test.interval_list"


def run(command: list[str], env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True, env=env)


def hdf5_dictionary(path: pathlib.Path) -> list[str]:
    """Read the variable-length dictionary dataset without requiring h5py."""
    strings = shutil.which("strings")
    if strings is None:
        raise RuntimeError("coreutils strings is required for the HDF5 metadata audit")
    text = subprocess.run([strings, str(path)], check=True, text=True,
                          capture_output=True).stdout.splitlines()
    return [line for line in text if line.startswith("@HD\t") or line.startswith("@SQ\t")]


def denoise_with_java(source: pathlib.Path, output_prefix: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    standardized = output_prefix.with_suffix(".standardized.tsv")
    denoised = output_prefix.with_suffix(".denoised.tsv")
    run([str(JAVA), "-jar", str(JAR), "DenoiseReadCounts", "-I", str(source),
         "--standardized-copy-ratios", str(standardized),
         "--denoised-copy-ratios", str(denoised)])
    if not standardized.is_file() or not denoised.is_file():
        raise AssertionError(f"Java DenoiseReadCounts did not write both outputs for {source}")
    return standardized, denoised


def main() -> int:
    for path in (NATIVE, JAVA, JAR, BAM, INTERVALS):
        if not path.exists():
            raise SystemExit(f"missing HDF5 metadata audit input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hdf5-simple-count-collection-") as directory:
        work = pathlib.Path(directory)
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")

        native_hdf5 = work / "native.hdf5"
        run([str(NATIVE), "-I", str(BAM), "-L", str(INTERVALS), "-O", str(native_hdf5),
             "--format", "HDF5", "--sample", "NA12878"], env)
        if not native_hdf5.is_file() or native_hdf5.stat().st_size == 0:
            raise AssertionError("native SimpleCountCollection HDF5 output is empty")

        # Produce a release-pinned HDF5 file from the same BAM/intervals.  The
        # binary payload need not be byte-identical (HDF5 allocation order is
        # not a GATK contract), but metadata and Java-decoded values must be.
        java_hdf5 = work / "java.hdf5"
        run([str(JAVA), "-jar", str(JAR), "CollectReadCounts", "-I", str(BAM),
             "-L", str(INTERVALS), "-O", str(java_hdf5), "--format", "HDF5",
             "--interval-merging-rule", "OVERLAPPING_ONLY"], env)

        native_dictionary = hdf5_dictionary(native_hdf5)
        java_dictionary = hdf5_dictionary(java_hdf5)
        if native_dictionary != java_dictionary:
            raise AssertionError({"native_dictionary": native_dictionary[:3],
                                   "java_dictionary": java_dictionary[:3],
                                   "native_count": len(native_dictionary),
                                   "java_count": len(java_dictionary)})
        if native_dictionary[0] != "@HD\tVN:1.6" or not any("\tAS:" in line and "\tM5:" in line
                                                                 for line in native_dictionary[1:]):
            raise AssertionError("native HDF5 sequence dictionary lost HTSJDK @SQ metadata")

        # A Java reader is the oracle for the sample/dictionary/interval/count
        # dataset layout.  Decoding both HDF5 files through the same release
        # also checks that the native interval matrix and count row are not
        # merely self-consistent in the native reader.
        native_standardized, native_denoised = denoise_with_java(native_hdf5, work / "native")
        java_standardized, java_denoised = denoise_with_java(java_hdf5, work / "java")
        if native_standardized.read_bytes() != java_standardized.read_bytes():
            raise AssertionError("Java decoded native HDF5 counts differ from Java HDF5 counts")
        if native_denoised.read_bytes() != java_denoised.read_bytes():
            raise AssertionError("Java decoded native HDF5 denoised output differs")

        raw_strings = subprocess.run([shutil.which("strings") or "strings", str(native_hdf5)],
                                     check=True, text=True, capture_output=True).stdout
        if "sample_metadata" not in raw_strings or "locatable_metadata" not in raw_strings or "NA12878" not in raw_strings:
            raise AssertionError("native HDF5 omitted sample/locatable metadata paths")

        print(
            '{"status":"pass","gatk_version":"4.6.2.0",'
            '"schema":"HDF5SimpleCountCollection",'
            '"sample_dictionary_interval_counts_java_roundtrip":true,'
            '"sequence_dictionary_sq_tags":true,"pon_svd_scope":"separate-contract"}'
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
