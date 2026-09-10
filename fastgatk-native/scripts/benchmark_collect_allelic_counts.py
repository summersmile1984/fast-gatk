#!/usr/bin/env python3
"""File-boundary benchmark for CollectAllelicCounts."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_ALLELIC_COUNTS_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build")) /
        "fastgatk-collect-allelic-counts"),
))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", type=int, default=10_000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.reads < 1:
        raise SystemExit("--reads must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-allelic-bench-") as temporary:
        work = pathlib.Path(temporary)
        reference = work / "reference.fasta"
        sam = work / "reads.sam"
        output = work / "allelic.tsv"
        manifest = work / "allelic.manifest.json"
        sequence = "ACGT" * 2500
        reference.write_text(f">chr1\n{sequence}\n", encoding="ascii")
        (work / "reference.fasta.fai").write_text(
            f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n", encoding="ascii")
        bases = "ACGT" * 25
        quality = "?" * len(bases)
        with sam.open("w", encoding="ascii") as handle:
            handle.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len(sequence)}\n")
            handle.write("@RG\tID:GATKCopyNumber\tSM:BENCH\n")
            for index in range(args.reads):
                start = index % 900 + 1
                handle.write(f"r{index}\t0\tchr1\t{start}\t60\t100M\t*\t0\t0\t{bases}\t{quality}\n")
        begin = time.perf_counter()
        result = subprocess.run(
            [str(BINARY), "-I", str(sam), "-R", str(reference), "-L", "chr1:1-1000",
             "-O", str(output), "--sample", "BENCH", "--threads", str(args.threads),
             "--output-manifest", str(manifest)],
            # Keep the lifecycle telemetry in the benchmark payload.
            # The output-manifest is a normal GATK-compatible sidecar option.
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        rows = sum(1 for line in output.read_text(encoding="ascii").splitlines()
                   if line and not line.startswith("@") and not line.startswith("CONTIG"))
        telemetry = json.loads(manifest.read_text(encoding="ascii"))["telemetry"]
        print(json.dumps({
            "status": "pass",
            "tool": "CollectAllelicCounts",
            "reads": args.reads,
            "rows": rows,
            "wall_seconds": elapsed,
            "reads_per_second": args.reads / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "kernel": {
                "execution_space": telemetry["kernel_execution_space"],
                "batches": telemetry["kernel_batches"],
                "observations": telemetry["kernel_observations"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["kernel_execute_seconds"],
            },
            "pipeline": {
                "decoded_items": telemetry["pipeline_decoded_items"],
                "computed_items": telemetry["pipeline_computed_items"],
                "encoded_items": telemetry["pipeline_encoded_items"],
                "peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
                "peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
                "peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
            },
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
