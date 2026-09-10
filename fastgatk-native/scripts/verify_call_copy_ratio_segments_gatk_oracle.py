#!/usr/bin/env python3
"""Pinned GATK oracle for the finite-statistics CallCopyRatioSegments edge.

The SimpleCopyRatioCaller in GATK 4.6.2.0 intentionally lets an empty or
singleton copy-neutral set flow through IEEE-754 statistics.  A native caller
must not turn those valid cases into an exception or a zero standard deviation:
the resulting NaN makes every non-neutral segment call neutral.  This oracle
also checks GATK's automatically derived ``.igv.seg`` sidecar.
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


HEADER = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
    "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n"
)


def table(values: list[tuple[int, int, float]]) -> str:
    return HEADER + "".join(
        f"chr1\t{start}\t{end}\t10\t{mean}\n"
        for start, end, mean in values
    )


def run_native(input_path: pathlib.Path, output_path: pathlib.Path,
               extra_args: list[str] | None = None) -> dict:
    manifest = output_path.with_suffix(".json")
    result = subprocess.run(
        [str(NATIVE), "-I", str(input_path), "-O", str(output_path),
         "--output-manifest", str(manifest), *(extra_args or [])],
        text=True, capture_output=True, env=os.environ.copy(),
    )
    if result.returncode != 0:
        raise AssertionError(f"native failed:\n{result.stdout}\n{result.stderr}")
    return json.loads(manifest.read_text(encoding="ascii"))


def run_java(input_path: pathlib.Path, output_path: pathlib.Path,
             extra_args: list[str] | None = None) -> None:
    result = subprocess.run(
        [str(JAVA), "-jar", str(JAR), "CallCopyRatioSegments",
         "-I", str(input_path), "-O", str(output_path), *(extra_args or [])],
        text=True, capture_output=True, env=os.environ.copy(),
    )
    if result.returncode != 0:
        raise AssertionError(f"GATK failed:\n{result.stdout}\n{result.stderr[-5000:]}")


def compare_case(label: str, values: list[tuple[int, int, float]], work: pathlib.Path,
                 *, java_args: list[str] | None = None,
                 native_args: list[str] | None = None) -> None:
    input_path = work / f"{label}.cr.seg"
    native_output = work / f"{label}.called.seg"
    java_output = work / f"{label}.java.called.seg"
    input_path.write_text(table(values), encoding="ascii")
    manifest = run_native(input_path, native_output, native_args)
    run_java(input_path, java_output, java_args)
    if native_output.read_bytes() != java_output.read_bytes():
        raise AssertionError({
            "label": label,
            "native": native_output.read_text(encoding="ascii"),
            "java": java_output.read_text(encoding="ascii"),
        })
    native_legacy = native_output.with_name(native_output.stem + ".igv.seg")
    java_legacy = java_output.with_name(java_output.stem + ".igv.seg")
    if not native_legacy.exists() or not java_legacy.exists():
        raise AssertionError({"label": label, "native_legacy": str(native_legacy),
                              "java_legacy": str(java_legacy)})
    if native_legacy.read_bytes() != java_legacy.read_bytes():
        raise AssertionError({
            "label": label,
            "native_legacy": native_legacy.read_text(encoding="ascii"),
            "java_legacy": java_legacy.read_text(encoding="ascii"),
        })
    if manifest["statistics_mean"] is not None or manifest["statistics_sd"] is not None:
        raise AssertionError({"label": label, "manifest": manifest})


def main() -> int:
    if not NATIVE.exists():
        raise SystemExit(f"missing native binary: {NATIVE}")
    if not JAVA.exists() or not JAR.exists():
        raise SystemExit("missing pinned GATK 4.6.2.0 runtime")
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-gatk-") as directory:
        work = pathlib.Path(directory)
        compare_case("one-neutral", [
            (1, 100, 0.0), (101, 200, 1.5), (201, 300, -1.5),
        ], work)
        compare_case("none-neutral", [
            (1, 100, 1.5), (101, 200, -1.5),
        ], work)
        compare_case("full-z-score-aliases", [
            (1, 100, 0.0), (101, 200, 1.5), (201, 300, -1.5),
        ], work,
                     java_args=[
                         "--outlier-neutral-segment-copy-ratio-z-score-threshold", "2.0",
                         "--calling-copy-ratio-z-score-threshold", "2.0",
                     ],
                     native_args=[
                         "--outlier-neutral-segment-copy-ratio-z-score-threshold", "2.0",
                         "--calling-copy-ratio-z-score-threshold", "2.0",
                     ])
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "singleton_neutral_nan_statistics_exact": True,
            "empty_neutral_nan_statistics_exact": True,
            "all_calls_neutral_in_degenerate_cases": True,
            "derived_igv_sidecar_exact": True,
            "full_gatk_z_score_aliases": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
