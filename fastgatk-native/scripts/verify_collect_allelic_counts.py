#!/usr/bin/env python3
"""Contract checks for CollectAllelicCounts' HTSlib/Kokkos path."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_ALLELIC_COUNTS_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-collect-allelic-counts"),
))


def write_reference(path: pathlib.Path) -> None:
    # chr1 positions 5 and 6 are N and therefore intentionally omitted from output.
    sequence = "ACGTNNGTAC" + "ACGT" * 10
    path.write_text(f">chr1\n{sequence}\n", encoding="ascii")
    # FASTA: >chr1\n is six bytes; the sequence is one unwrapped line.
    path.with_suffix(path.suffix + ".fai").write_text(
        f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n",
        encoding="ascii",
    )


def write_sam(path: pathlib.Path) -> None:
    bases = "ACGTACGTAC"
    qualities = "?" * len(bases)
    low_first = "4" + "?" * (len(bases) - 1)  # Phred 19 then Phred 30.
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:50",
        "@RG\tID:GATKCopyNumber\tSM:NA1",
        f"r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\t{bases}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\t{bases}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r3\t0\tchr1\t1\t60\t10M\t*\t0\t0\t{'A' * 10}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r4\t0\tchr1\t1\t60\t10M\t*\t0\t0\t{'C' * 10}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r5\t0\tchr1\t1\t20\t10M\t*\t0\t0\t{'G' * 10}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r6\t1024\tchr1\t1\t60\t10M\t*\t0\t0\t{bases}\t{qualities}\tRG:Z:GATKCopyNumber",
        f"r7\t0\tchr1\t1\t60\t10M\t*\t0\t0\t{bases}\t{low_first}\tRG:Z:GATKCopyNumber",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def run(sam: pathlib.Path, reference: pathlib.Path, output: pathlib.Path,
        manifest: pathlib.Path, include_duplicates: bool = False,
        sample: str | None = "NA1") -> tuple[dict, list[list[str]]]:
    command = [str(BINARY), "-I", str(sam), "-R", str(reference), "-L", "chr1:1-10",
               "-O", str(output), "--output-manifest", str(manifest),
               "--batch-records", "2"]
    if sample is not None:
        command.extend(["--sample", sample])
    if include_duplicates:
        command.append("--include-duplicates")
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    rows = [line.split("\t") for line in output.read_text(encoding="ascii").splitlines()
            if line and not line.startswith("@") and not line.startswith("CONTIG")]
    return json.loads(result.stdout), rows


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-allelic-") as temporary:
        work = pathlib.Path(temporary)
        reference = work / "reference.fasta"
        sam = work / "reads.sam"
        output = work / "counts.tsv"
        manifest = work / "counts.json"
        write_reference(reference)
        write_sam(sam)
        summary, rows = run(sam, reference, output, manifest)
        assert len(rows) == 8, rows  # positions 5 and 6 have reference N.
        assert rows[0][:6] == ["chr1", "1", "3", "1", "A", "C"], rows[0]
        assert rows[1][:6] == ["chr1", "2", "4", "1", "C", "A"], rows[1]
        assert all(row[0] == "chr1" for row in rows)
        assert summary["reads_seen"] == 7
        assert summary["reads_used"] == 5  # r1-r4 and r7; r5 MAPQ and r6 duplicate filtered.
        auto_output = work / "counts-auto.tsv"
        auto_manifest = work / "counts-auto.json"
        auto_summary, _ = run(sam, reference, auto_output, auto_manifest, sample=None)
        assert auto_summary["sample"] == "NA1"
        assert auto_summary["sample_source"] == "read-group-header"
        assert json.loads(auto_manifest.read_text(encoding="ascii"))["sample_source"] == "read-group-header"
        metadata = json.loads(manifest.read_text(encoding="ascii"))
        assert metadata["tool"] == "CollectAllelicCounts"
        assert metadata["execution_space"] in {"OpenMP", "Serial"}
        assert metadata["determinism"] == "strict"
        assert metadata["minimum_base_quality"] == 20
        assert metadata["reads_filtered"] == 2
        assert metadata["telemetry"]["output_bytes"] == output.stat().st_size
        assert metadata["telemetry"]["wall_seconds"] >= 0.0
        telemetry = metadata["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] >= 1
        assert telemetry["kernel_observations"] == summary["projected_bases"]
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        assert telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] >= 1
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0
        assert all(item["complete"] for item in metadata["outputs"])
        duplicate_summary, duplicate_rows = run(
            sam, reference, work / "counts-with-duplicates.tsv", work / "counts-with-duplicates.json", True)
        assert duplicate_summary["reads_used"] == 6
        assert int(duplicate_rows[0][2]) == 4
        print(json.dumps({
            "status": "pass",
            "rows": len(rows),
            "reads_seen": summary["reads_seen"],
            "reads_used": summary["reads_used"],
            "include_duplicates_reads_used": duplicate_summary["reads_used"],
            "reference_n_rows_omitted": 2,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
