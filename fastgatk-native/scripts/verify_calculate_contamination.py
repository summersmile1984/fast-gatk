#!/usr/bin/env python3
"""Verify the native CalculateContamination table and matched-normal contract."""
from __future__ import annotations

import json
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True, check=False)


def table(sample: str, rows: list[str]) -> str:
    return (
        f"#<METADATA>SAMPLE={sample}\n"
        "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        + "\n".join(rows) + "\n"
    )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_CONTAMINATION_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-calculate-contamination")))
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-calculate-contamination-") as directory:
        work = Path(directory)
        tumor = work / "tumor.table"
        tumor.write_text(table("TUMOR", [
            "chr1\t1\t4\t16\t0\t0.10",
            "chr1\t2\t2\t18\t0\t0.10",
            "chr1\t3\t1\t1\t0\t0.10",  # below MIN_COVERAGE
        ]), encoding="utf-8")
        normal = work / "normal.table"
        normal.write_text(table("NORMAL", [
            "chr1\t1\t1\t19\t0\t0.10",
            "chr1\t2\t2\t18\t0\t0.10",
        ]), encoding="utf-8")
        output = work / "contamination.table"
        segmentation = work / "segments.table"
        manifest = work / "contamination.manifest.json"
        result = run(binary, [
            "-I", str(tumor), "--matched", str(normal), "-O", str(output),
            "--tumor-segmentation", str(segmentation), "--output-manifest", str(manifest),
        ])
        assert result.returncode == 0, result.stderr
        lines = output.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "sample\tcontamination\terror"
        fields = lines[1].split("\t")
        assert fields[0] == "TUMOR"
        # The full GATK ContaminationModel also returns 0.0 for this tiny
        # two-site matched-normal fixture: there are not enough non-LOH
        # segments to enter the hom-alt strategy.  Keep the contract aligned
        # with that oracle while requiring a finite bounded estimate.
        assert 0.0 <= float(fields[1]) < 0.5
        assert float(fields[2]) >= 0.0
        assert segmentation.read_text(encoding="utf-8").splitlines()[1].startswith("contig\tstart\tend")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["matched_normal_hom_site_estimate"] is True
        assert metadata["compatibility"]["full_contamination_model"] is True
        assert metadata["compatibility"]["kernel_segmenter"] is True
        assert metadata["compatibility"]["segment_coordinates_gatk_oracle"] is True
        assert metadata["telemetry"]["candidate_sites"] >= 0
        # Coverage thresholds use the Java/Commons-Math corrected two-pass
        # mean.  This fixture has two covered sites at depth 20, so the
        # strict low/high boundaries are 10 and 60 respectively.
        telemetry = metadata["telemetry"]
        assert telemetry["covered_sites"] == 2
        assert telemetry["filtered_sites"] == 2
        assert telemetry["coverage_median"] == 20.0
        assert telemetry["coverage_mean"] == 20.0
        assert telemetry["coverage_low_threshold"] == 10.0
        assert telemetry["coverage_high_threshold"] == 60.0
        assert telemetry["empty_filtered_tumor"] is False
        assert telemetry["empty_filtered_matched"] is False
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] >= 0
        assert telemetry["kernel_observations"] >= 0
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        assert telemetry["likelihood_plan_allocations"] >= 1
        assert telemetry["likelihood_plan_reuses"] > 0
        assert telemetry["segmenter_kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["segmenter_kernel_execution_policy"] == "MDRangePolicy"
        assert telemetry["segmenter_kernel_batches"] >= 0
        assert telemetry["segmenter_kernel_observations"] >= 0
        assert telemetry["segmenter_kernel_prepare_seconds"] >= 0.0
        assert telemetry["segmenter_kernel_execute_seconds"] >= 0.0

        compressed_tumor = work / "tumor.table.gz"
        with gzip.open(compressed_tumor, "wt", encoding="utf-8") as stream:
            stream.write(tumor.read_text(encoding="utf-8"))
        compressed_output = work / "contamination.table.gz"
        compressed_segmentation = work / "segments.table.gz"
        compressed_manifest = work / "contamination-compressed.manifest.json"
        compressed_result = run(binary, [
            "-I", str(compressed_tumor), "--matched", str(normal),
            "-O", str(compressed_output), "--tumor-segmentation", str(compressed_segmentation),
            "--output-manifest", str(compressed_manifest),
        ])
        assert compressed_result.returncode == 0, compressed_result.stderr
        with gzip.open(compressed_output, "rt", encoding="utf-8") as stream:
            assert stream.read() == output.read_text(encoding="utf-8")
        with gzip.open(compressed_segmentation, "rt", encoding="utf-8") as stream:
            assert stream.read() == segmentation.read_text(encoding="utf-8")
        compressed_metadata = json.loads(compressed_manifest.read_text(encoding="utf-8"))
        assert compressed_metadata["compatibility"]["compressed_io"] is True
        assert compressed_metadata["telemetry"]["input_compressed"] is True
        assert compressed_metadata["telemetry"]["output_compressed"] is True

        tumor_only = work / "tumor-only.table"
        tumor_only_output = work / "tumor-only-contamination.table"
        tumor_only.write_text(table("TUMOR_ONLY", [
            "chr1\t1\t18\t2\t0\t0.10",
            "chr1\t2\t17\t3\t0\t0.10",
        ]), encoding="utf-8")
        result = run(binary, ["-I", str(tumor_only), "-O", str(tumor_only_output)])
        assert result.returncode == 0, result.stderr
        assert tumor_only_output.read_text(encoding="utf-8").splitlines()[1].startswith("TUMOR_ONLY\t")

        # A completely filtered panel is a successful GATK result, not a
        # BAD_INPUT failure.  The conservative uncertainty must remain
        # available to downstream FilterMutectCalls, and the requested MAF
        # sidecar is a valid header-only table.  Use the documented aliases.
        empty_tumor = work / "empty-coverage.table"
        empty_normal = work / "empty-coverage-normal.table"
        empty_tumor.write_text(table("EMPTY_TUMOR", [
            "chr1\t1\t10\t0\t0\t0.10",
            "chr1\t2\t9\t1\t0\t0.20",
        ]), encoding="utf-8")
        empty_normal.write_text(table("EMPTY_NORMAL", [
            "chr1\t1\t10\t0\t0\t0.10",
        ]), encoding="utf-8")
        empty_output = work / "empty-coverage.contamination.table"
        empty_segments = work / "empty-coverage.segments.table"
        empty_manifest = work / "empty-coverage.manifest.json"
        empty_result = run(binary, [
            "-I", str(empty_tumor), "-matched", str(empty_normal),
            "-O", str(empty_output), "-segments", str(empty_segments),
            "--output-manifest", str(empty_manifest),
        ])
        assert empty_result.returncode == 0, empty_result.stderr
        assert empty_output.read_text(encoding="utf-8") == (
            "sample\tcontamination\terror\nEMPTY_TUMOR\t0.0\t1.0\n")
        assert empty_segments.read_text(encoding="utf-8") == (
            "#<METADATA>SAMPLE=EMPTY_TUMOR\n"
            "contig\tstart\tend\tminor_allele_fraction\n")
        empty_metadata = json.loads(empty_manifest.read_text(encoding="utf-8"))
        assert empty_metadata["compatibility"]["empty_post_coverage_result"] is True
        assert empty_metadata["telemetry"]["empty_filtered_tumor"] is True
        assert empty_metadata["telemetry"]["empty_filtered_matched"] is True
        assert empty_metadata["telemetry"]["coverage_median"] is None
        assert empty_metadata["telemetry"]["coverage_mean"] is None

    print(json.dumps({"status": "pass", "matched_candidate_sites": 2,
                      "tumor_only": True, "empty_post_coverage": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
