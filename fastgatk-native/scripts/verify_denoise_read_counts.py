#!/usr/bin/env python3
"""Contract checks for the TSV DenoiseReadCounts native paths."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_DENOISE_READ_COUNTS_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-denoise-read-counts")))
COLLECT_BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_READ_COUNTS_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-collect-read-counts")))


def write_counts(path: pathlib.Path, counts: list[float]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\n")
        handle.write("@SQ\tSN:chr1\tLN:1000\n")
        handle.write("@RG\tID:fastgatk\tSM:SAMPLE\n")
        handle.write("CONTIG\tSTART\tEND\tCOUNT\n")
        for index, count in enumerate(counts):
            start = index * 100 + 1
            handle.write(f"chr1\t{start}\t{start + 99}\t{count}\n")


def rows(path: pathlib.Path) -> list[float]:
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@") or line.startswith("CONTIG"):
            continue
        values.append(float(line.split("\t")[3]))
    return values


def write_panel_matrix(path: pathlib.Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("CONTIG\tSTART\tEND\tN1\tN2\tN3\n")
        matrix = [
            (10, 12, 8),
            (20, 18, 22),
            (30, 31, 29),
            (40, 39, 41),
            (50, 52, 48),
        ]
        for index, values in enumerate(matrix):
            start = index * 100 + 1
            handle.write(f"chr1\t{start}\t{start + 99}\t{values[0]}\t{values[1]}\t{values[2]}\n")


def write_annotated_intervals(path: pathlib.Path, gc_values: list[float]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("CONTIG\tSTART\tEND\tGC_CONTENT\n")
        for index, gc in enumerate(gc_values):
            start = index * 100 + 1
            handle.write(f"chr1\t{start}\t{start + 99}\t{gc}\n")


def run(input_path: pathlib.Path, output_path: pathlib.Path, manifest: pathlib.Path) -> dict:
    standardized = output_path.with_name(output_path.stem + ".standardized.tsv")
    result = subprocess.run(
        [str(BINARY), "-I", str(input_path), "-O", str(output_path),
         "--standardized-copy-ratios", str(standardized),
         "--output-manifest", str(manifest), "--threads", "2"],
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(result.stdout)


def run_with_annotations(input_path: pathlib.Path, output_path: pathlib.Path,
                         manifest: pathlib.Path, annotated: pathlib.Path,
                         panel: pathlib.Path) -> dict:
    standardized = output_path.with_name(output_path.stem + ".standardized.tsv")
    result = subprocess.run(
        [str(BINARY), "-I", str(input_path), "-O", str(output_path),
         "--standardized-copy-ratios", str(standardized),
         "--panel-of-normals", str(panel), "--annotated-intervals", str(annotated),
         "--output-manifest", str(manifest), "--threads", "2"],
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(result.stdout)


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-denoise-") as temporary:
        work = pathlib.Path(temporary)
        plain = work / "plain.tsv"
        doubled = work / "doubled.tsv"
        write_counts(plain, [10, 20, 30, 40, 50])
        write_counts(doubled, [20, 40, 60, 80, 100])
        plain_out = work / "plain.den.tsv"
        doubled_out = work / "doubled.den.tsv"
        plain_manifest = work / "plain.json"
        doubled_manifest = work / "doubled.json"
        summary = run(plain, plain_out, plain_manifest)
        run(doubled, doubled_out, doubled_manifest)
        plain_values = rows(plain_out)
        doubled_values = rows(doubled_out)
        assert len(plain_values) == 5
        assert plain_values == doubled_values, (plain_values, doubled_values)
        assert abs(plain_values[2]) < 1e-6
        assert plain_values[0] < 0 < plain_values[-1]
        manifest = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert manifest["tool"] == "DenoiseReadCounts"
        assert manifest["format"] == "TSV"
        assert manifest["execution_space"] in {"OpenMP", "Serial"}
        assert manifest["determinism"] == "strict"
        assert manifest["panel_of_normals"] is False
        assert manifest["svd"] is False
        telemetry = manifest["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "MDRangePolicy+RangePolicy"
        assert telemetry["kernel_batches"] >= 1
        assert telemetry["kernel_observations"] >= 1
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        assert manifest["telemetry"]["output_bytes"] == plain_out.stat().st_size
        assert manifest["telemetry"]["standardized_output_bytes"] > 0
        assert manifest["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in manifest["outputs"])

        pon = work / "pon.tsv"
        write_counts(pon, [10, 10, 30, 40, 50])
        pon_out = work / "pon.den.tsv"
        pon_manifest = work / "pon.json"
        pon_summary = subprocess.run(
            [str(BINARY), "-I", str(plain), "-O", str(pon_out),
             "--standardized-copy-ratios", str(work / "pon.standardized.tsv"),
             "--panel-of-normals", str(pon), "--output-manifest", str(pon_manifest)],
            text=True, capture_output=True, check=True,
        )
        pon_values = rows(pon_out)
        assert len(pon_values) == 5
        assert pon_values[1] > pon_values[0]  # the PoN has a lower baseline at interval 2
        pon_payload = json.loads(pon_manifest.read_text(encoding="utf-8"))
        assert pon_payload["panel_of_normals"] is True
        assert json.loads(pon_summary.stdout)["panel_of_normals"] is True

        # GATK ignores --annotated-intervals when a panel is supplied.  Use
        # deliberately extreme GC values to make an accidental override
        # observable, and require byte-identical output plus telemetry that
        # reports GC correction as unused for this TSV panel path.
        annotated = work / "annotated.tsv"
        write_annotated_intervals(annotated, [0.0, 1.0, 0.0, 1.0, 0.0])
        pon_annotated_out = work / "pon.annotated.den.tsv"
        pon_annotated_manifest = work / "pon.annotated.json"
        run_with_annotations(plain, pon_annotated_out, pon_annotated_manifest, annotated, pon)
        assert pon_annotated_out.read_bytes() == pon_out.read_bytes()
        pon_annotated_payload = json.loads(pon_annotated_manifest.read_text(encoding="utf-8"))
        assert pon_annotated_payload["gc_content"] is False

        matrix = work / "pon-matrix.tsv"
        write_panel_matrix(matrix)
        matrix_out = work / "matrix.den.tsv"
        matrix_manifest = work / "matrix.json"
        matrix_summary = subprocess.run(
            [str(BINARY), "-I", str(plain), "-O", str(matrix_out),
             "--standardized-copy-ratios", str(work / "matrix.standardized.tsv"),
             "--panel-of-normals", str(matrix), "--number-of-eigensamples", "1",
             "--output-manifest", str(matrix_manifest)],
            text=True, capture_output=True, check=True,
        )
        matrix_values = rows(matrix_out)
        matrix_payload = json.loads(matrix_manifest.read_text(encoding="utf-8"))
        assert len(matrix_values) == 5 and all(map(lambda value: abs(value) < 100, matrix_values))
        assert matrix_payload["panel_samples"] == 3
        assert matrix_payload["svd"] is True
        assert matrix_payload["svd_rank"] == 1
        assert json.loads(matrix_summary.stdout)["svd"] is True

        zero = work / "zero.tsv"
        zero_out = work / "zero.den.tsv"
        write_counts(zero, [0, 0, 0])
        failed = subprocess.run([str(BINARY), "-I", str(zero), "-O", str(zero_out),
                                 "--standardized-copy-ratios", str(work / "zero.standardized.tsv")],
                                text=True, capture_output=True)
        assert failed.returncode != 0
        assert "positive sample median" in failed.stderr
        sam = work / "hdf5.sam"
        intervals = work / "hdf5.interval_list"
        sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:500\n"
                       "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                       "r2\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                       "r3\t0\tchr1\t201\t60\t4M\t*\t0\t0\tACGT\tIIII\n", encoding="utf-8")
        intervals.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:500\n"
                             "chr1\t1\t100\t+\ta\nchr1\t101\t200\t+\tb\nchr1\t201\t300\t+\tc\n", encoding="utf-8")
        hdf5_input = work / "counts.h5"
        subprocess.run([str(COLLECT_BINARY), "-I", str(sam), "-L", str(intervals),
                        "-O", str(hdf5_input), "--format", "HDF5"],
                       text=True, capture_output=True, check=True)
        hdf5_out = work / "hdf5.den.tsv"
        hdf5_result = subprocess.run(
            [str(BINARY), "-I", str(hdf5_input), "-O", str(hdf5_out),
             "--standardized-copy-ratios", str(work / "hdf5.standardized.tsv"),
             "--format", "HDF5"],
            text=True, capture_output=True, check=True,
        )
        assert json.loads(hdf5_result.stdout)["format"] == "HDF5"
        # CollectReadCounts defaults to GATK OVERLAPPING_ONLY: adjacent
        # half-open intervals remain distinct unless --interval-merging-rule
        # ALL is requested explicitly.
        assert len(rows(hdf5_out)) == 3
        hdf5_panel_out = work / "hdf5.panel.den.tsv"
        hdf5_panel_result = subprocess.run(
            [str(BINARY), "-I", str(hdf5_input), "-O", str(hdf5_panel_out), "--format", "HDF5",
             "--standardized-copy-ratios", str(work / "hdf5.panel.standardized.tsv"),
             "--panel-of-normals", str(hdf5_input)],
            text=True, capture_output=True, check=True,
        )
        assert json.loads(hdf5_panel_result.stdout)["panel_of_normals"] is True
        assert len(rows(hdf5_panel_out)) == 3
        print(json.dumps({
            "status": "pass",
            "intervals": len(plain_values),
            "scale_free": True,
            "zero_row_rejected": True,
            "hdf5": True,
            "panel_of_normals_tsv": True,
            "panel_matrix_svd": True,
            "panel_ignores_annotated_intervals": True,
            "hdf5_panel": True,
            "summary": summary,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
