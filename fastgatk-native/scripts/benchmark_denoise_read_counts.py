#!/usr/bin/env python3
"""Small file-boundary benchmark for the TSV DenoiseReadCounts path."""

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
    "FASTGATK_DENOISE_READ_COUNTS_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-denoise-read-counts")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", type=int, default=100_000)
    parser.add_argument("--panel-samples", type=int, default=0,
                        help="wide TSV PoN samples; 0 benchmarks the scalar path")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.intervals < 3:
        raise SystemExit("--intervals must be at least 3")
    if args.panel_samples < 0:
        raise SystemExit("--panel-samples must be non-negative")
    with tempfile.TemporaryDirectory(prefix="fastgatk-denoise-bench-") as temporary:
        work = pathlib.Path(temporary)
        input_path = work / "counts.tsv"
        output_path = work / "copy-ratio.tsv"
        standardized_path = work / "standardized-copy-ratio.tsv"
        panel_path = work / "panel.tsv"
        manifest_path = work / "denoise.manifest.json"
        with input_path.open("w", encoding="utf-8") as handle:
            handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000000\n")
            handle.write("@RG\tID:fastgatk\tSM:BENCH\nCONTIG\tSTART\tEND\tCOUNT\n")
            for index in range(args.intervals):
                start = index * 100 + 1
                handle.write(f"chr1\t{start}\t{start + 99}\t{10 + (index % 97)}\n")
        if args.panel_samples > 1:
            with panel_path.open("w", encoding="utf-8") as handle:
                names = [f"PANEL_{index + 1}" for index in range(args.panel_samples)]
                handle.write("CONTIG\tSTART\tEND\t" + "\t".join(names) + "\n")
                for index in range(args.intervals):
                    start = index * 100 + 1
                    values = [str(10 + (index % 97) + ((index + sample) % 5))
                              for sample in range(args.panel_samples)]
                    handle.write(f"chr1\t{start}\t{start + 99}\t" + "\t".join(values) + "\n")
        command = [str(BINARY), "-I", str(input_path), "-O", str(output_path),
                   "--standardized-copy-ratios", str(standardized_path),
                   "--threads", str(args.threads), "--output-manifest", str(manifest_path)]
        if args.panel_samples > 1:
            command += ["--panel-of-normals", str(panel_path), "--number-of-eigensamples",
                        str(min(10, args.panel_samples))]
        begin = time.perf_counter()
        result = subprocess.run(
            command,
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in output_path.read_text(encoding="utf-8").splitlines()
                      if line and not line.startswith("@") and not line.startswith("CONTIG"))
        telemetry = json.loads(manifest_path.read_text(encoding="utf-8"))["telemetry"]
        print(json.dumps({
            "status": "pass",
            "tool": "DenoiseReadCounts",
            "intervals": args.intervals,
            "panel_samples": args.panel_samples,
            "records": records,
            "wall_seconds": elapsed,
            "intervals_per_second": records / elapsed if elapsed else 0.0,
            "output_bytes": output_path.stat().st_size,
            "kernel": {
                "execution_space": telemetry["kernel_execution_space"],
                "policy": telemetry["kernel_execution_policy"],
                "batches": telemetry["kernel_batches"],
                "observations": telemetry["kernel_observations"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["kernel_execute_seconds"],
            },
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
