#!/usr/bin/env python3
"""Pinned GATK oracle for multisample vector-kernel segmentation.

GATK's ``MultisampleMultidimensionalKernelSegmenter`` evaluates a kernel on
the complete sample vector at each interval and sums the per-sample kernels.
Reducing those vectors to a per-interval arithmetic mean is not equivalent:
anti-correlated sample events cancel under the mean while they remain a clear
change point in GATK's linear kernel.  This fixture pins both the
anti-correlated and concordant two-sample boundaries to GATK 4.6.2.0.
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


def write_copy_ratios(path: pathlib.Path, sample: str, values: list[float]) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:10000",
        f"@RG\tID:GATKCopyNumber{sample}\tSM:{sample}",
        "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO",
    ]
    for index, value in enumerate(values):
        lines.append(f"chr1\t{index * 100 + 1}\t{index * 100 + 100}\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def interval_rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")[:3]
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("@")
    ]


def run_case(work: pathlib.Path, label: str, first: list[float], second: list[float]) -> dict[str, object]:
    first_path = work / f"{label}.s1.tsv"
    second_path = work / f"{label}.s2.tsv"
    write_copy_ratios(first_path, "S1", first)
    write_copy_ratios(second_path, "S2", second)
    native_dir = work / f"{label}.native"
    java_dir = work / f"{label}.java"
    native_dir.mkdir()
    java_dir.mkdir()
    manifest = native_dir / "manifest.json"
    run([
        str(NATIVE), "--denoised-copy-ratios", str(first_path),
        "--denoised-copy-ratios", str(second_path),
        "--output-prefix", "sample", "-O", str(native_dir),
        "--output-manifest", str(manifest),
    ])
    run([
        str(JAVA), "-jar", str(JAR), "ModelSegments",
        "--denoised-copy-ratios", str(first_path),
        "--denoised-copy-ratios", str(second_path),
        "--output-prefix", "sample", "-O", str(java_dir),
    ])
    native_rows = interval_rows(native_dir / "sample.interval_list")
    java_rows = interval_rows(java_dir / "sample.interval_list")
    if native_rows != java_rows:
        raise AssertionError({"label": label, "native": native_rows, "java": java_rows})
    native_manifest = json.loads(manifest.read_text(encoding="ascii"))
    if native_manifest["sample_count"] != 2 or native_manifest["segmentation_method"] != "KernelSegmenter":
        raise AssertionError({"label": label, "manifest": native_manifest})
    return {
        "segments": len(native_rows),
        "boundaries": native_rows,
        "boundaries_exact": True,
        "joint_feature_observations": native_manifest["telemetry"]["kernel_observations"],
    }


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing ModelSegments multisample oracle input: {path}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-multisample-") as directory:
        work = pathlib.Path(directory)
        first = [0.0] * 16 + [1.0] * 16
        anti = run_case(work, "anti-correlated", first, list(reversed(first)))
        concordant = run_case(work, "concordant", first, first)
        expected = [["chr1", "1", "1600"], ["chr1", "1601", "3200"]]
        assert anti["boundaries"] == expected, anti
        assert concordant["boundaries"] == expected, concordant
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "sample_count": 2,
            "anti_correlated_vector_boundary_exact": True,
            "concordant_vector_boundary_exact": True,
            "mean_reduction_would_cancel_anti_correlated_event": True,
            "cases": {"anti_correlated": anti, "concordant": concordant},
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
