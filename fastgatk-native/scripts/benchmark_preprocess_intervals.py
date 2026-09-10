#!/usr/bin/env python3
"""Reproducible file-boundary benchmark for native PreprocessIntervals."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def write_reference(path: pathlib.Path, length: int) -> None:
    sequence = ("ACGT" * ((length + 3) // 4))[:length]
    header = f">chr1\n".encode()
    body = (sequence + "\n").encode()
    with path.open("wb") as output:
        output.write(header)
        output.write(body)
    path.with_suffix(path.suffix + ".fai").write_text(
        f"chr1\t{length}\t{len(header)}\t{length}\t{length + 1}\n", encoding="utf-8")
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:{length}\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", type=int, default=10000)
    parser.add_argument("--bin-length", type=int, default=1000)
    parser.add_argument("--padding", type=int, default=250)
    parser.add_argument("--excluded-intervals", type=int, default=0,
                        help="number of deterministic -XL intervals to subtract")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if (args.intervals < 1 or args.bin_length < 0 or args.padding < 0 or
            args.excluded_intervals < 0 or args.excluded_intervals > args.intervals or
            args.threads < 1):
        raise SystemExit("intervals/bin-length/padding/excluded-intervals must be non-negative, excluded-intervals <= intervals, threads positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_PREPROCESS_INTERVALS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-preprocess-intervals"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-preprocess-benchmark-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        reference_length = max(1000, args.intervals * max(1, args.bin_length // 2 + 1))
        write_reference(reference, reference_length)
        selectors = work / "targets.interval_list"
        spacing = max(1, reference_length // args.intervals)
        with selectors.open("w", encoding="utf-8") as output:
            output.write(f"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:{reference_length}\n")
            for index in range(args.intervals):
                start = index * spacing + 1
                end = min(reference_length, start + max(1, spacing // 3) - 1)
                output.write(f"chr1\t{start}\t{end}\t+\tbin{index}\n")
        exclude_path = None
        if args.excluded_intervals:
            exclude_path = work / "excluded.interval_list"
            with exclude_path.open("w", encoding="utf-8") as output:
                output.write(f"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:{reference_length}\n")
                for ordinal in range(args.excluded_intervals):
                    index = (ordinal * args.intervals) // args.excluded_intervals
                    start = index * spacing + 1
                    end = min(reference_length, start + max(1, spacing // 3) - 1)
                    output.write(f"chr1\t{start}\t{end}\t+\tmask{index}\n")
        destination = work / "preprocessed.interval_list"
        manifest = work / "preprocessed.manifest.json"
        command = [str(native), "-R", str(reference), "-L", str(selectors),
                   "--bin-length", str(args.bin_length), "--padding", str(args.padding),
                   "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(destination),
                   "--threads", str(args.threads), "--output-manifest", str(manifest)]
        if exclude_path is not None:
            command.extend(["-XL", str(exclude_path)])
        started = time.perf_counter()
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - started
        summary = json.loads(result.stdout.splitlines()[-1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        rows = sum(1 for line in destination.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("@"))
        print(json.dumps({"status": "pass", "tool": "PreprocessIntervals",
                          "input_intervals": args.intervals,
                          "unfiltered_bins": summary["unfiltered_bins"],
                          "excluded_intervals": summary.get("excluded_intervals", 0),
                          "intervals_after_exclusion": summary.get("intervals_after_exclusion"),
                          "output_intervals": rows,
                          "intervals_per_second": args.intervals / elapsed,
                          "bin_length": args.bin_length, "padding": args.padding,
                          "requested_excluded_intervals": args.excluded_intervals,
                          "threads": args.threads, "wall_seconds": elapsed,
                          "output_bytes": destination.stat().st_size,
                          "kernel": {
                              "execution_space": telemetry["kernel_execution_space"],
                              "batches": telemetry["kernel_batches"],
                              "records": telemetry["kernel_records"],
                              "prepare_seconds": telemetry["kernel_prepare_seconds"],
                              "execute_seconds": telemetry["kernel_execute_seconds"],
                          }}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
