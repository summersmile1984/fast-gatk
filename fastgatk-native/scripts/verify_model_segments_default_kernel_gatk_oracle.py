#!/usr/bin/env python3
"""Pinned GATK oracle for ModelSegments' default KernelSegmenter dispatch.

GATK uses KernelSegmenter whenever ``--segments`` is omitted.  The native
adapter historically used its threshold extension unless a kernel option was
spelled explicitly.  This oracle keeps that workflow boundary observable on a
linear copy-ratio fixture (the rank-one SVD path) while avoiding Java MCMC
posterior byte identity.
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


HEADER = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:chr1\tLN:100000\n"
    "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
    "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
)


def write_copy_ratios(path: pathlib.Path, values: list[float]) -> None:
    path.write_text(
        HEADER + "".join(
            f"chr1\t{index * 100 + 1}\t{index * 100 + 100}\t{value}\n"
            for index, value in enumerate(values)
        ), encoding="ascii")


def bounds(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")[:3]
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def run_case(work: pathlib.Path, label: str, values: list[float]) -> dict[str, object]:
    input_path = work / f"{label}.tsv"
    native_dir = work / f"{label}.native"
    java_dir = work / f"{label}.java"
    native_dir.mkdir()
    java_dir.mkdir()
    write_copy_ratios(input_path, values)
    native_manifest = native_dir / "manifest.json"
    native = subprocess.run(
        [str(NATIVE), "--denoised-copy-ratios", str(input_path),
         "--output-prefix", "sample", "-O", str(native_dir),
         "--maximum-number-of-smoothing-iterations", "0",
         "--output-manifest", str(native_manifest)],
        text=True, capture_output=True, check=False)
    if native.returncode != 0:
        raise AssertionError(f"native failed:\n{native.stdout}\n{native.stderr}")
    java = subprocess.run(
        [str(JAVA), "-jar", str(JAR), "ModelSegments",
         "--denoised-copy-ratios", str(input_path),
         "--output-prefix", "sample", "-O", str(java_dir),
         "--maximum-number-of-smoothing-iterations", "0"],
        text=True, capture_output=True, check=False)
    if java.returncode != 0:
        raise AssertionError(f"GATK failed:\n{java.stdout}\n{java.stderr[-5000:]}")
    native_bounds = bounds(native_dir / "sample.modelFinal.seg")
    java_bounds = bounds(java_dir / "sample.modelFinal.seg")
    if native_bounds != java_bounds:
        raise AssertionError({"label": label, "native": native_bounds, "java": java_bounds})
    manifest = json.loads(native_manifest.read_text(encoding="ascii"))
    if manifest["segmentation_method"] != "KernelSegmenter":
        raise AssertionError({"label": label, "manifest": manifest})
    return {
        "native_segments": len(native_bounds),
        "java_segments": len(java_bounds),
        "boundaries_exact": True,
        "segmentation_method": manifest["segmentation_method"],
    }


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing ModelSegments default-kernel oracle input: {path}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-default-kernel-") as directory:
        work = pathlib.Path(directory)
        cases = {
            "step": [0.0] * 16 + [1.0] * 16,
            "gradient": [index * 0.02 for index in range(32)],
            "blip": [0.0] * 8 + [1.0] * 8 + [0.0] * 8 + [1.0] * 8,
        }
        results = {label: run_case(work, label, values) for label, values in cases.items()}
        assert results["step"]["native_segments"] == 2, results
        assert results["gradient"]["native_segments"] == 1, results
        assert results["blip"]["native_segments"] == 1, results
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "default_without_segments_uses_kernel_segmenter": True,
            "linear_rank_one_kernel_feature": True,
            "cases": results,
            "boundaries_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
