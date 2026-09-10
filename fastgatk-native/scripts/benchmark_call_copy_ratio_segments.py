#!/usr/bin/env python3
"""File-boundary benchmark for CallCopyRatioSegments."""

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
    "FASTGATK_CALL_COPY_RATIO_BINARY",
    ROOT / "fastgatk-native" / "build" / "fastgatk-call-copy-ratio-segments",
))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--segments", type=int, default=100_000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.segments < 1:
        raise SystemExit("--segments must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-bench-") as temporary:
        work = pathlib.Path(temporary)
        input_path = work / "segments.tsv"
        output_path = work / "called.tsv"
        manifest_path = work / "called.manifest.json"
        with input_path.open("w", encoding="utf-8") as handle:
            handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000000\n@RG\tID:GATKCopyNumber\tSM:BENCH\n")
            handle.write("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n")
            for index in range(args.segments):
                start = index * 100 + 1
                value = (index % 200 - 100) / 2000.0
                handle.write(f"chr1\t{start}\t{start + 99}\t10\t{value}\n")
        begin = time.perf_counter()
        result = subprocess.run(
            [str(BINARY), "-I", str(input_path), "-O", str(output_path), "--threads", str(args.threads),
             "--output-manifest", str(manifest_path)],
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in output_path.read_text(encoding="utf-8").splitlines()
                      if line and not line.startswith("@") and not line.startswith("CONTIG"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        telemetry = manifest["telemetry"]
        print(json.dumps({
            "status": "pass",
            "tool": "CallCopyRatioSegments",
            "segments": args.segments,
            "records": records,
            "wall_seconds": elapsed,
            "segments_per_second": records / elapsed if elapsed else 0.0,
            "output_bytes": output_path.stat().st_size,
            "execution_space": manifest["execution_space"],
            "simd_width": telemetry["copy_ratio_kernel_simd_width"],
            "simd_groups": telemetry["copy_ratio_kernel_simd_groups"],
            "kernel": {
                "lifecycle": telemetry["kernel_lifecycle"],
                "execution_policy": telemetry["kernel_execution_policy"],
                "batches": telemetry["kernel_batches"],
                "observations": telemetry["kernel_observations"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["copy_ratio_kernel_execute_seconds"],
            },
            "copy_ratio_kernel_execute_seconds": telemetry["copy_ratio_kernel_execute_seconds"],
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
