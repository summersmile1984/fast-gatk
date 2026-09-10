#!/usr/bin/env python3
"""Pinned GATK oracle for credible-interval segment smoothing.

GATK's ``MultidimensionalModeller`` merges adjacent segments when the
posterior median difference is less than the configured threshold times the
10--90% credible width of either segment.  This verifier exercises both sides
of that decision with an explicitly supplied two-segment partition and checks
the native final boundaries against GATK 4.6.2.0.
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
    "FASTGATK_MODEL_SEGMENTS_BINARY", str(BUILD / "fastgatk-model-segments")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError("command failed: " + " ".join(command) +
                             "\nstdout:\n" + result.stdout[-2000:] +
                             "\nstderr:\n" + result.stderr[-4000:])
    return result


def write_copy_ratios(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
        "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
        "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
        + "\n".join(
            f"chr1\t{index * 10 + 1}\t{index * 10 + 10}\t{value:.6f}"
            for index, value in enumerate((0.000, 0.004, 0.006, 0.008, 0.010, 0.012))
        ) + "\n", encoding="utf-8")


def write_segments(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
        "chr1\t1\t30\t+\tseg1\nchr1\t31\t60\t+\tseg2\n",
        encoding="utf-8")


def data_rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def run_case(work: pathlib.Path, copy_ratios: pathlib.Path, segments: pathlib.Path,
             threshold: float, native_name: str, java_name: str) -> dict[str, object]:
    native_dir = work / native_name
    java_dir = work / java_name
    native_dir.mkdir()
    java_dir.mkdir()
    native_manifest = native_dir / "manifest.json"
    native_result = run([
        str(NATIVE), "--denoised-copy-ratios", str(copy_ratios), "--segments", str(segments),
        "--output-prefix", "sample", "-O", str(native_dir),
        "--num-burn-in-iterations", "4", "--num-samples", "16",
        "--maximum-number-of-smoothing-iterations", "1",
        "--number-of-smoothing-iterations-per-fit", "0",
        "--smoothing-credible-interval-threshold-copy-ratio", str(threshold),
        "--smoothing-credible-interval-threshold-allele-fraction", str(threshold),
        "--output-manifest", str(native_manifest),
    ])
    run([
        str(JAVA), "-jar", str(JAR), "ModelSegments",
        "--denoised-copy-ratios", str(copy_ratios), "--segments", str(segments),
        "--output-prefix", "sample", "-O", str(java_dir),
        "--number-of-samples-copy-ratio", "16",
        "--number-of-burn-in-samples-copy-ratio", "4",
        "--number-of-samples-allele-fraction", "16",
        "--number-of-burn-in-samples-allele-fraction", "4",
        "--maximum-number-of-smoothing-iterations", "1",
        "--number-of-smoothing-iterations-per-fit", "0",
        "--smoothing-credible-interval-threshold-copy-ratio", str(threshold),
        "--smoothing-credible-interval-threshold-allele-fraction", str(threshold),
    ])
    native_rows = data_rows(native_dir / "sample.modelFinal.seg")
    java_rows = data_rows(java_dir / "sample.modelFinal.seg")
    native_bounds = [row[:3] for row in native_rows]
    java_bounds = [row[:3] for row in java_rows]
    assert native_bounds == java_bounds, {"native": native_rows, "java": java_rows}
    manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
    return {
        "native_segments": len(native_rows),
        "java_segments": len(java_rows),
        "boundaries_exact": True,
        "smoothing_applied": manifest["smoothing_applied"],
        "native_summary": json.loads(native_result.stdout.splitlines()[-1]),
    }


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing ModelSegments smoothing oracle input: {path}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-smoothing-") as directory:
        work = pathlib.Path(directory)
        copy_ratios = work / "copy-ratios.tsv"
        segments = work / "segments.interval_list"
        write_copy_ratios(copy_ratios)
        write_segments(segments)
        merged = run_case(work, copy_ratios, segments, 100.0, "native-merged", "java-merged")
        separate = run_case(work, copy_ratios, segments, 0.0, "native-separate", "java-separate")
        assert merged["native_segments"] == 1 and merged["java_segments"] == 1, merged
        assert separate["native_segments"] == 2 and separate["java_segments"] == 2, separate
        assert merged["smoothing_applied"] is True, merged
        assert separate["smoothing_applied"] is False, separate
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "merge_threshold": 100.0,
            "merge_case": merged,
            "non_merge_threshold": 0.0,
            "non_merge_case": separate,
            "credible_interval_smoothing_boundaries_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
