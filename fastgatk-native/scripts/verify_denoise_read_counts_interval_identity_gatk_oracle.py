#!/usr/bin/env python3
"""Pinned GATK oracle for DenoiseReadCounts PoN interval identity.

GATK 4.6.2.0 requires a case SimpleCountCollection's interval list to be
identical to the original interval list used to build an SVD PoN.  A native
implementation that only looks for the post-filter panel intervals can
silently accept subsets or supersets; that changes the fractional-coverage
denominator and is not a direct replacement.  This oracle creates a native
v7 PoN, proves a valid case is accepted by both engines, and proves subset,
superset, and coordinate-mismatch cases are rejected by both engines.
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
CREATE_PON = pathlib.Path(os.environ.get(
    "FASTGATK_CREATE_PON_BINARY", str(BUILD / "fastgatk-create-read-count-panel-of-normals")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def write_count_table(path: pathlib.Path, sample: str, intervals: list[tuple[str, int, int, int]]) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:1000",
        f"@RG\tID:{sample}\tSM:{sample}",
        "CONTIG\tSTART\tEND\tCOUNT",
    ]
    lines.extend("\t".join(str(value) for value in row) for row in intervals)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=check)


def run_denoise(command_binary: pathlib.Path, case: pathlib.Path, pon: pathlib.Path,
                output_prefix: pathlib.Path) -> subprocess.CompletedProcess[str]:
    command = ([str(JAVA), "-jar", str(JAR), "DenoiseReadCounts"]
               if command_binary == JAVA else [str(command_binary)])
    return run(command + [
        "-I", str(case), "--count-panel-of-normals", str(pon),
        "--standardized-copy-ratios", str(output_prefix) + ".standardized.tsv",
        "--denoised-copy-ratios", str(output_prefix) + ".denoised.tsv",
        "--number-of-eigensamples", "1",
    ], check=False)


def main() -> int:
    for path in (DENOISE, CREATE_PON, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing DenoiseReadCounts interval oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-denoise-interval-identity-") as directory:
        work = pathlib.Path(directory)
        original = [("chr1", index * 100 + 1, index * 100 + 100, 100 + index)
                    for index in range(4)]
        normals = []
        for sample_index, delta in enumerate((0, -2, 3), 1):
            normal = work / f"normal-{sample_index}.tsv"
            write_count_table(
                normal, f"NORMAL_{sample_index}",
                [(contig, start, end, count + delta)
                 for contig, start, end, count in original],
            )
            normals.append(normal)

        pon = work / "panel.pon.hdf5"
        pon_result = run([
            str(CREATE_PON),
            *sum((["-I", str(path)] for path in normals), []),
            "-O", str(pon), "--number-of-eigensamples", "1",
            "--minimum-interval-median-percentile", "0",
            "--maximum-zeros-in-sample-percentage", "100",
            "--maximum-zeros-in-interval-percentage", "100",
            "--extreme-sample-median-percentile", "0",
        ])
        assert pon.is_file() and pon.stat().st_size > 0, pon_result.stderr

        valid = work / "valid.tsv"
        write_count_table(valid, "CASE", original)
        native_valid = run_denoise(DENOISE, valid, pon, work / "native-valid")
        assert native_valid.returncode == 0, native_valid.stderr
        java_valid = run_denoise(JAVA, valid, pon, work / "java-valid")
        assert java_valid.returncode == 0, java_valid.stderr
        assert (work / "native-valid.standardized.tsv").is_file()
        assert (work / "java-valid.standardized.tsv").is_file()

        cases = {
            # GATK validates equality before it subsets filtered panel intervals.
            "subset": original[:-1],
            "superset": original + [("chr1", 401, 500, 104)],
            "coordinate_mismatch": original[:2] + [("chr1", 201, 299, 102)] + original[3:],
        }
        native_rejected: dict[str, bool] = {}
        java_rejected: dict[str, bool] = {}
        for name, intervals in cases.items():
            case = work / f"{name}.tsv"
            write_count_table(case, "CASE", intervals)
            native = run_denoise(DENOISE, case, pon, work / f"native-{name}")
            java = run_denoise(JAVA, case, pon, work / f"java-{name}")
            native_rejected[name] = native.returncode != 0
            java_rejected[name] = java.returncode != 0
            assert native_rejected[name], {"case": name, "stderr": native.stderr}
            assert java_rejected[name], {"case": name, "stderr": java.stderr}
            assert "identical" in native.stderr or "interval" in native.stderr
            assert "identical" in java.stderr or "interval" in java.stderr

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "valid_case_accepted": True,
            "interval_identity": True,
            "native_rejected": native_rejected,
            "java_rejected": java_rejected,
            "poN_schema": "HDF5-SVD-ReadCountPanelOfNormals-v7",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
