#!/usr/bin/env python3
"""Compare the native ContaminationModel against the bundled GATK oracle."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


FIXTURES = (
    "NA12891_0.01_NA12892_0.99.table",
    "NA12891_0.03_NA12892_0.97.table",
    "NA12891_0.05_NA12892_0.95.table",
    "NA12891_0.08_NA12892_0.92.table",
)


def read_result(path: Path) -> tuple[float, float]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    fields = lines[-1].split("\t")
    if len(fields) != 3:
        raise AssertionError({"path": str(path), "lines": lines[-4:]})
    return float(fields[1]), float(fields[2])


def read_segments(path: Path) -> list[tuple[str, int, int, float]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("contig"):
            continue
        contig, start, end, maf = line.split("\t")
        result.append((contig, int(start), int(end), float(maf)))
    return result


def write_empty_coverage_table(path: Path, sample: str) -> None:
    """Write a valid table whose rows all fail GATK's strict depth > 10 gate."""
    path.write_text(
        f"#<METADATA>SAMPLE={sample}\n"
        "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        "chr1\t1\t10\t0\t0\t0.1\n"
        "chr1\t2\t9\t1\t0\t0.2\n",
        encoding="utf-8",
    )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-calculate-contamination"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk"
    java = root / "third_party/jdk17/bin"
    fixture_root = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/calculatecontamination"
    if not native.is_file() or not gatk.is_file() or not java.joinpath("java").is_file():
        oracle_guard.oracle_not_verified('verify_calculate_contamination_gatk_oracle.py', gatk, java.joinpath('java'))
        raise AssertionError("native/GATK/JDK oracle artifacts are required")

    report = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-contamination-oracle-") as directory:
        work = Path(directory)
        for fixture in FIXTURES:
            source = fixture_root / fixture
            native_out = work / f"native-{fixture}"
            java_out = work / f"java-{fixture}"
            native_segments = work / f"native-{fixture}.segments"
            java_segments = work / f"java-{fixture}.segments"
            native_run = subprocess.run(
                [str(native), "-I", str(source), "-O", str(native_out),
                 "--tumor-segmentation", str(native_segments)],
                text=True, capture_output=True, check=False,
            )
            if native_run.returncode != 0:
                raise AssertionError({"fixture": fixture, "stderr": native_run.stderr})
            environment = os.environ.copy()
            environment["PATH"] = f"{java}:{environment.get('PATH', '')}"
            java_run = subprocess.run(
                ["python3", str(gatk), "--java-options", "-Xmx1g",
                 "CalculateContamination", "-I", str(source), "-O", str(java_out),
                 "--tumor-segmentation", str(java_segments)],
                text=True, capture_output=True, check=False, env=environment,
            )
            if java_run.returncode != 0:
                raise AssertionError({"fixture": fixture, "stderr": java_run.stderr[-4000:]})
            native_value = read_result(native_out)
            java_value = read_result(java_out)
            contamination_delta = abs(native_value[0] - java_value[0])
            error_delta = abs(native_value[1] - java_value[1])
            contamination_bit_exact = native_value[0].hex() == java_value[0].hex()
            error_bit_exact = native_value[1].hex() == java_value[1].hex()
            native_segment_rows = read_segments(native_segments)
            java_segment_rows = read_segments(java_segments)
            native_sorted = sorted(native_segment_rows)
            java_sorted = sorted(java_segment_rows)
            segment_coordinates_exact = [row[:3] for row in native_sorted] == [row[:3] for row in java_sorted]
            segment_file_order_exact = [row[:3] for row in native_segment_rows] == [row[:3] for row in java_segment_rows]
            if not segment_coordinates_exact or len(native_segment_rows) != len(java_segment_rows):
                raise AssertionError({
                    "fixture": fixture,
                    "native_segments": len(native_segment_rows),
                    "gatk_segments": len(java_segment_rows),
                    "segment_coordinates_exact": segment_coordinates_exact,
                })
            if not segment_file_order_exact:
                raise AssertionError({
                    "fixture": fixture,
                    "segment_file_order_exact": False,
                    "native_first": [row[:3] for row in native_segment_rows[:5]],
                    "gatk_first": [row[:3] for row in java_segment_rows[:5]],
                })
            maf_delta = max(
                (abs(native_row[3] - java_row[3])
                 for native_row, java_row in zip(native_segment_rows, java_segment_rows)),
                default=0.0,
            )
            # The native path uses the same Commons-Math Brent tolerances and
            # Java-seeded kernel segmenter as GATK.  Keep this oracle strict:
            # the remaining differences are only platform libm/VCF text
            # round-trip noise, not a biologically meaningful tolerance.
            if contamination_delta > 2.0e-12 or error_delta > 2.0e-12 or maf_delta > 1.0e-10:
                raise AssertionError({
                    "fixture": fixture,
                    "native": native_value,
                    "java": java_value,
                    "contamination_delta": contamination_delta,
                    "error_delta": error_delta,
                    "segment_maf_delta": maf_delta,
                })
            if not contamination_bit_exact or not error_bit_exact:
                raise AssertionError({
                    "fixture": fixture,
                    "contamination_bit_exact": contamination_bit_exact,
                    "error_bit_exact": error_bit_exact,
                    "native": native_value,
                    "gatk": java_value,
                })
            report.append({
                "fixture": fixture,
                "native_contamination": native_value[0],
                "gatk_contamination": java_value[0],
                "contamination_delta": contamination_delta,
                "error_delta": error_delta,
                "contamination_bit_exact": contamination_bit_exact,
                "error_bit_exact": error_bit_exact,
                "segments": len(native_segment_rows),
                "segment_coordinates_exact": segment_coordinates_exact,
                "segment_file_order_exact": segment_file_order_exact,
                "segment_maf_delta": maf_delta,
            })
        # GATK does not reject a panel whose sites are all removed by the
        # strict MIN_COVERAGE (>10) filter.  It emits 0.0 +/- 1.0 and a
        # header-only tumor MAF segmentation sidecar.  Exercise the documented
        # Barclay aliases too, because established pipelines use
        # -matched/-segments rather than only their long forms.
        empty_tumor = work / "empty-coverage-tumor.table"
        empty_normal = work / "empty-coverage-normal.table"
        write_empty_coverage_table(empty_tumor, "EMPTY_TUMOR")
        write_empty_coverage_table(empty_normal, "EMPTY_NORMAL")
        native_empty_out = work / "native-empty-coverage.table"
        java_empty_out = work / "java-empty-coverage.table"
        native_empty_segments = work / "native-empty-coverage.segments"
        java_empty_segments = work / "java-empty-coverage.segments"
        native_empty_run = subprocess.run(
            [str(native), "-I", str(empty_tumor), "-matched", str(empty_normal),
             "-O", str(native_empty_out), "-segments", str(native_empty_segments)],
            text=True, capture_output=True, check=False,
        )
        if native_empty_run.returncode != 0:
            raise AssertionError({"empty_coverage_native_stderr": native_empty_run.stderr})
        environment = os.environ.copy()
        environment["PATH"] = f"{java}:{environment.get('PATH', '')}"
        java_empty_run = subprocess.run(
            ["python3", str(gatk), "--java-options", "-Xmx1g", "CalculateContamination",
             "-I", str(empty_tumor), "-matched", str(empty_normal),
             "-O", str(java_empty_out), "-segments", str(java_empty_segments)],
            text=True, capture_output=True, check=False, env=environment,
        )
        if java_empty_run.returncode != 0:
            raise AssertionError({"empty_coverage_gatk_stderr": java_empty_run.stderr[-4000:]})
        native_empty_text = native_empty_out.read_text(encoding="utf-8")
        java_empty_text = java_empty_out.read_text(encoding="utf-8")
        if native_empty_text != java_empty_text:
            raise AssertionError({
                "empty_coverage_output_exact": False,
                "native": native_empty_text,
                "gatk": java_empty_text,
            })
        native_empty_segmentation = native_empty_segments.read_text(encoding="utf-8")
        java_empty_segmentation = java_empty_segments.read_text(encoding="utf-8")
        if native_empty_segmentation != java_empty_segmentation:
            raise AssertionError({
                "empty_coverage_segmentation_exact": False,
                "native": native_empty_segmentation,
                "gatk": java_empty_segmentation,
            })
        if read_result(native_empty_out) != (0.0, 1.0):
            raise AssertionError({"empty_coverage_native_result": read_result(native_empty_out)})
        report.append({
            "fixture": "all-sites-filtered-by-coverage-with-matched-short-alias",
            "contamination_bit_exact": True,
            "error_bit_exact": True,
            "segment_coordinates_exact": True,
            "segment_file_order_exact": True,
            "header_only_segmentation": True,
        })
    print(json.dumps({"status": "pass", "fixtures": report}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
