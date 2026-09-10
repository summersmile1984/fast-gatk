#!/usr/bin/env python3
"""File-boundary benchmark for native CreateReadCountPanelOfNormals."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--intervals", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.samples < 1 or args.intervals < 2 or args.threads < 1:
        raise SystemExit("samples/intervals/threads must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_CREATE_PON_BINARY",
        str(root / "fastgatk-native/build/fastgatk-create-read-count-panel-of-normals")))
    with tempfile.TemporaryDirectory(prefix="fastgatk-pon-benchmark-") as directory:
        work = pathlib.Path(directory)
        inputs = []
        for sample in range(args.samples):
            path = work / f"normal-{sample}.tsv"
            lines = ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:100000000",
                     f"@RG\tID:N{sample}\tSM:N{sample}", "CONTIG\tSTART\tEND\tCOUNT"]
            for interval in range(args.intervals):
                start = interval * 100 + 1
                lines.append(f"chr1\t{start}\t{start + 99}\t{100 + ((interval + sample * 7) % 101)}")
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            inputs.append(path)
        output = work / "panel.hdf5"
        manifest = work / "panel.manifest.json"
        command = [str(binary)]
        for path in inputs:
            command += ["-I", str(path)]
        command += ["-O", str(output), "--minimum-interval-median-percentile", "0",
                    "--maximum-zeros-in-sample-percentage", "100",
                    "--maximum-zeros-in-interval-percentage", "100",
                    "--extreme-sample-median-percentile", "0",
                    "--number-of-eigensamples", str(min(20, args.samples)),
                    "--output-manifest", str(manifest),
                    "--threads", str(args.threads)]
        begin = time.perf_counter()
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(result.stdout.splitlines()[-1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        print(json.dumps({"status": "pass", "tool": "CreateReadCountPanelOfNormals",
                          "samples": args.samples, "intervals": args.intervals,
                          "eigensamples": summary["eigensamples"],
                          "wall_seconds": elapsed,
                          "intervals_per_second": args.intervals / elapsed if elapsed else 0.0,
                          "output_bytes": output.stat().st_size,
                          "kernel": {
                              "execution_space": telemetry["kernel_execution_space"],
                              "policy": telemetry["kernel_execution_policy"],
                              "batches": telemetry["kernel_batches"],
                              "observations": telemetry["kernel_observations"],
                              "prepare_seconds": telemetry["kernel_prepare_seconds"],
                              "execute_seconds": telemetry["kernel_execute_seconds"],
                          }}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
