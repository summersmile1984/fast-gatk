#!/usr/bin/env python3
"""Verify the native CollectReadCounts TSV/interval contract."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True, env=env)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_COLLECT_READ_COUNTS_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-collect-read-counts"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-read-counts-") as directory:
        work = Path(directory)
        sam = work / "reads.sam"
        intervals = work / "targets.interval_list"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            "r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            "r3\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            "r4\t1024\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            "r5\t0\tchr1\t1\t20\t4M\t*\t0\t0\tACGT\tIIII\n",
            encoding="utf-8",
        )
        intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
            "chr1\t1\t5\t+\ttarget-1\nchr1\t11\t20\t+\ttarget-2\n",
            encoding="utf-8",
        )
        output = work / "counts.tsv"
        manifest = work / "counts.manifest.json"
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        result = json.loads(run([
            str(native), "-I", str(sam), "-L", str(intervals), "-O", str(output),
            "--format", "TSV", "--sample", "SAMPLE1", "--output-manifest", str(manifest)
        ], env).stdout.splitlines()[-1])
        assert result["tool"] == "CollectReadCounts"
        assert result["status"] == "prototype"
        assert result["format"] == "TSV"
        assert result["intervals"] == 2
        assert result["interval_merging_rule"] == "OVERLAPPING_ONLY"
        assert result["reads_seen"] == 5 and result["reads_used"] == 3
        assert result["reads_filtered"] == 2
        assert result["output_bytes"] == output.stat().st_size
        assert result["hdf5_roundtrip"] is False
        assert result["wall_seconds"] >= 0.0
        rows = [line.split("\t") for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("@")] 
        assert rows == [["CONTIG", "START", "END", "COUNT"],
                        ["chr1", "1", "5", "2"], ["chr1", "11", "20", "1"]]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["tsv"] is True
        assert metadata["compatibility"]["hdf5"] is False
        assert metadata["compatibility"]["exclude_intervals"] is True
        assert metadata["compatibility"]["zero_interval_padding_required"] is True
        assert metadata["execution_space"] in {"OpenMP", "Serial"}
        assert metadata["determinism"] == "strict"
        assert metadata["telemetry"]["output_bytes"] == output.stat().st_size
        assert metadata["telemetry"]["hdf5_roundtrip"] is False
        assert metadata["telemetry"]["wall_seconds"] >= 0.0
        assert metadata["telemetry"]["interval_file_inputs"] == 1
        assert metadata["telemetry"]["interval_file_records"] == 2
        telemetry = metadata["telemetry"]
        assert telemetry["count_kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["count_kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["count_kernel_execution_policy"] == "RangePolicy"
        assert telemetry["count_kernel_batches"] >= 1
        assert telemetry["count_kernel_records"] == telemetry["reads_seen"]
        assert telemetry["count_kernel_prepare_seconds"] >= 0.0
        assert telemetry["count_kernel_execute_seconds"] >= 0.0
        assert telemetry["pipeline_lifecycle"] == "decode->compute->encode->sink"
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] == telemetry["count_kernel_batches"]
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0

        # GATK's command default is HDF5; TSV is intentionally explicit so a
        # direct dispatcher replacement cannot emit text under a .hdf5 path.
        default_hdf5 = work / "default.hdf5"
        default_result = json.loads(run([
            str(native), "-I", str(sam), "-L", str(intervals), "-O", str(default_hdf5),
            "--sample", "SAMPLE1",
        ], env).stdout.splitlines()[-1])
        assert default_result["format"] == "HDF5"
        assert default_result["hdf5_roundtrip"] is True
        assert default_hdf5.is_file() and default_hdf5.stat().st_size > 0

        # GATK's OVERLAPPING_ONLY interval-merging rule must union overlapping
        # targets before counting so a read start is never double-counted.
        overlap_output = work / "counts-overlap.tsv"
        overlap_manifest = work / "counts-overlap.manifest.json"
        overlap_result = json.loads(run([
            str(native), "-I", str(sam), "-L", "chr1:1-5", "-L", "chr1:4-8",
            "-O", str(overlap_output), "--format", "TSV", "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--output-manifest", str(overlap_manifest)
        ], env).stdout.splitlines()[-1])
        assert overlap_result["intervals"] == 1
        assert overlap_result["requested_intervals"] == 2
        assert overlap_result["interval_merges"] == 1
        assert overlap_result["reads_used"] == 2
        overlap_rows = [line.split("\t") for line in overlap_output.read_text(encoding="utf-8").splitlines()
                        if line and not line.startswith("@")] 
        assert overlap_rows[-1] == ["chr1", "1", "8", "2"]
        overlap_metadata = json.loads(overlap_manifest.read_text(encoding="utf-8"))
        assert overlap_metadata["compatibility"]["interval_merging_rule"] == "OVERLAPPING_ONLY"
        assert overlap_metadata["telemetry"]["requested_intervals"] == 2
        assert overlap_metadata["telemetry"]["interval_merges"] == 1

        # GATK's CopyNumberArgumentValidationUtils rejects ALL even though
        # the generic interval option advertises it.  Keep this native
        # failure explicit rather than silently changing target geometry.
        all_rule = subprocess.run([
            str(native), "-I", str(sam), "-L", "chr1:1-5", "-L", "chr1:6-8",
            "-O", str(work / "counts-adjacent.tsv"), "--format", "TSV",
            "--interval-merging-rule", "ALL"
        ], text=True, capture_output=True, env=env)
        assert all_rule.returncode != 0
        assert "Interval merging rule must be set to OVERLAPPING_ONLY" in all_rule.stderr

        excluded_output = work / "counts-excluded.tsv"
        excluded_manifest = work / "counts-excluded.manifest.json"
        excluded_result = json.loads(run([
            str(native), "-I", str(sam), "-L", "chr1:1-12", "-XL", "chr1:5-8",
            "-O", str(excluded_output), "--format", "TSV",
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--output-manifest", str(excluded_manifest),
        ], env).stdout.splitlines()[-1])
        assert excluded_result["excluded_intervals"] == 1
        excluded_rows = [line.split("\t") for line in excluded_output.read_text(encoding="utf-8").splitlines()
                         if line and not line.startswith("@")]
        assert excluded_rows[-2:] == [["chr1", "1", "4", "1"], ["chr1", "9", "12", "1"]]
        assert json.loads(excluded_manifest.read_text(encoding="utf-8"))["telemetry"]["excluded_intervals"] == 1

        include_output = work / "counts-with-duplicates.tsv"
        include_result = json.loads(run([
            str(native), "-I", str(sam), "-L", str(intervals), "-O", str(include_output),
            "--format", "TSV", "--sample", "SAMPLE1", "--include-duplicates"
        ], env).stdout.splitlines()[-1])
        assert include_result["reads_used"] == 4
        include_rows = [line.split("\t") for line in include_output.read_text(encoding="utf-8").splitlines()
                        if line and not line.startswith("@")] 
        assert include_rows[-2:] == [["chr1", "1", "5", "3"], ["chr1", "11", "20", "1"]]
        hdf5_output = work / "counts.h5"
        hdf5_manifest = work / "counts.h5.manifest.json"
        hdf5_result = json.loads(run([
            str(native), "-I", str(sam), "-L", str(intervals), "-O", str(hdf5_output),
            "--format", "HDF5", "--sample", "SAMPLE1", "--output-manifest", str(hdf5_manifest)
        ], env).stdout.splitlines()[-1])
        assert hdf5_result["format"] == "HDF5"
        assert hdf5_result["hdf5_roundtrip"] is True
        assert hdf5_result["output_bytes"] == hdf5_output.stat().st_size
        assert hdf5_output.is_file() and hdf5_output.stat().st_size > 0
        hdf5_metadata = json.loads(hdf5_manifest.read_text(encoding="utf-8"))
        assert hdf5_metadata["compatibility"]["hdf5"] is True
        assert hdf5_metadata["determinism"] == "strict"
        assert hdf5_metadata["telemetry"]["hdf5_roundtrip"] is True
        assert hdf5_metadata["telemetry"]["output_bytes"] == hdf5_output.stat().st_size
        assert hdf5_metadata["outputs"][0]["kind"] == "hdf5"
        print(json.dumps({"status": "pass", "intervals": 2, "reads_seen": 5,
                          "reads_used": 3, "include_duplicates_reads_used": 4,
                          "hdf5": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
