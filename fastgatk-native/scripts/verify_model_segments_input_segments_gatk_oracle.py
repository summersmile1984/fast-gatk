#!/usr/bin/env python3
"""Pinned GATK oracle for ModelSegments ``--segments`` partition input.

GATK skips kernel segmentation when an input Picard interval list is supplied
and models exactly those intervals.  This verifier checks that the native
path accepts the same interval-list shape, preserves the supplied boundaries,
uses midpoint ownership for copy-ratio points, and agrees with the pinned
Java writer on copy-ratio segment rows.
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


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=check)


def write_copy_ratios(path: pathlib.Path) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:1000",
        "@RG\tID:GATKCopyNumber\tSM:SAMPLE",
        "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO",
    ]
    for index, value in enumerate((0.00, 0.02, 0.03, 0.80, 0.82, 0.81)):
        lines.append(f"chr1\t{index * 10 + 1}\t{index * 10 + 10}\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_segments(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
        "chr1\t1\t30\t+\tseg1\nchr1\t31\t60\t+\tseg2\n",
        encoding="utf-8",
    )


def data_rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing ModelSegments --segments oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-input-") as directory:
        work = pathlib.Path(directory)
        copy_ratios = work / "copy-ratios.tsv"
        segments = work / "segments.interval_list"
        write_copy_ratios(copy_ratios)
        write_segments(segments)

        native_dir = work / "native"
        java_dir = work / "java"
        native_dir.mkdir()
        java_dir.mkdir()
        native_manifest = work / "native.json"
        native_result = run([
            str(NATIVE), "-I", str(copy_ratios), "--segments", str(segments),
            "--output-prefix", "sample", "-O", str(native_dir),
            "--output-manifest", str(native_manifest),
        ])
        java_result = run([
            str(JAVA), "-jar", str(JAR), "ModelSegments",
            "--denoised-copy-ratios", str(copy_ratios), "--segments", str(segments),
            "--output-prefix", "sample", "-O", str(java_dir),
            "--number-of-samples-copy-ratio", "2",
            "--number-of-burn-in-samples-copy-ratio", "1",
            "--number-of-samples-allele-fraction", "2",
            "--number-of-burn-in-samples-allele-fraction", "1",
            "--maximum-number-of-smoothing-iterations", "0",
            "--number-of-smoothing-iterations-per-fit", "0",
        ])

        native_prefix = native_dir / "sample"
        java_prefix = java_dir / "sample"
        native_final = data_rows(native_prefix.with_suffix(".modelFinal.seg"))
        java_final = data_rows(java_prefix.with_suffix(".modelFinal.seg"))
        expected_intervals = [["chr1", "1", "30"], ["chr1", "31", "60"]]
        assert [row[:3] for row in native_final] == expected_intervals, native_final
        assert [row[:3] for row in java_final] == expected_intervals, java_final
        # Both implementations assign each copy-ratio point by interval
        # midpoint; the supplied partition contains 3 points per segment.
        assert [int(row[3]) for row in native_final] == [3, 3], native_final
        assert [int(row[3]) for row in java_final] == [3, 3], java_final

        native_cr = data_rows(native_prefix.with_suffix(".cr.seg"))
        java_cr = data_rows(java_prefix.with_suffix(".cr.seg"))
        assert native_cr == java_cr, {"native": native_cr, "java": java_cr}
        for suffix in (".modelBegin.cr.param", ".modelBegin.af.param",
                       ".modelFinal.cr.param", ".modelFinal.af.param"):
            parameter_path = native_prefix.with_suffix(suffix)
            assert parameter_path.is_file(), parameter_path
            parameter_lines = parameter_path.read_text(encoding="utf-8").splitlines()
            assert parameter_lines[2].startswith("PARAMETER_NAME\tPOSTERIOR_10"), parameter_lines
        assert native_prefix.with_suffix(".modelFinal.cr.param").read_text(encoding="utf-8").splitlines()[3].startswith("VARIANCE\t")
        assert native_prefix.with_suffix(".modelFinal.af.param").read_text(encoding="utf-8").splitlines()[3].startswith("MEAN_BIAS\t")
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert manifest["segmentation_method"] == "provided", manifest
        assert manifest["segments_input"] == str(segments), manifest
        assert manifest["segments"] == 2, manifest
        assert json.loads(native_result.stdout.splitlines()[-1])["segmentation_method"] == "provided"

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "provided_segments": 2,
            "copy_ratio_points_per_segment": [3, 3],
            "boundaries_exact": True,
            "copy_ratio_rows_exact": True,
            "segmentation_skipped": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
