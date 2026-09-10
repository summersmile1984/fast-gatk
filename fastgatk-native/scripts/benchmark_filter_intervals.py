#!/usr/bin/env python3
"""Reproducible file-boundary benchmark for annotation-based FilterIntervals."""

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
    parser.add_argument("--intervals", type=int, default=10000)
    parser.add_argument("--count-samples", type=int, default=0,
                        help="use TSV count filtering with this many synthetic samples instead of annotations")
    # FilterIntervals itself validates UNION only in GATK 4.6.2.0; keep the
    # benchmark surface aligned with that command contract.
    parser.add_argument("--interval-set-rule", choices=("UNION",), default="UNION")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.intervals < 2 or args.threads < 1 or args.count_samples < 0:
        raise SystemExit("intervals/threads must be positive and count-samples non-negative")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_FILTER_INTERVALS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-intervals"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-intervals-benchmark-") as directory:
        work = pathlib.Path(directory)
        interval_list = work / "targets.interval_list"
        annotated = work / "annotated.tsv"
        with interval_list.open("w", encoding="ascii") as interval_stream, annotated.open("w", encoding="ascii") as annotation_stream:
            interval_stream.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:" + str(args.intervals * 100) + "\n")
            annotation_stream.write("@HD\tVN:1.5\n@SQ\tSN:chr1\tLN:" + str(args.intervals * 100) + "\n")
            annotation_stream.write("CONTIG\tSTART\tEND\tGC_CONTENT\n")
            for index in range(args.intervals):
                start = index * 100 + 1
                end = start + 99
                interval_stream.write(f"chr1\t{start}\t{end}\t+\tI{index}\n")
                annotation_stream.write(f"chr1\t{start}\t{end}\t0.500000\n")
        output = work / "filtered.interval_list"
        manifest = work / "filtered.manifest.json"
        if args.count_samples:
            count_inputs = []
            for sample in range(args.count_samples):
                count_path = work / f"counts-{sample}.tsv"
                with count_path.open("w", encoding="ascii") as stream:
                    stream.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:" + str(args.intervals * 100) +
                                 f"\n@RG\tID:S{sample}\tSM:S{sample}\nCONTIG\tSTART\tEND\tCOUNT\n")
                    for index in range(args.intervals):
                        start = index * 100 + 1
                        stream.write(f"chr1\t{start}\t{start + 99}\t{100 + ((index + sample) % 17)}\n")
                count_inputs.append(count_path)
            command = [str(binary), "-L", str(interval_list), "--interval-merging-rule", "OVERLAPPING_ONLY",
                       "--interval-set-rule", args.interval_set_rule]
            for count_path in count_inputs:
                command += ["-I", str(count_path)]
            command += ["--low-count-filter-count-threshold", "0", "--low-count-filter-percentage-of-samples", "50",
                        "--extreme-count-filter-minimum-percentile", "1", "--extreme-count-filter-maximum-percentile", "99",
                        "--extreme-count-filter-percentage-of-samples", "90", "-O", str(output),
                        "--output-manifest", str(manifest),
                        "--threads", str(args.threads)]
            mode = "tsv-count"
        else:
            command = [str(binary), "-L", str(interval_list), "--annotated-intervals", str(annotated),
                       "--interval-merging-rule", "OVERLAPPING_ONLY", "--interval-set-rule", args.interval_set_rule,
                       "-O", str(output), "--output-manifest", str(manifest),
                       "--threads", str(args.threads)]
            mode = "annotation"
        begin = time.perf_counter()
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(completed.stdout.splitlines()[-1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        print(json.dumps({
            "status": "pass",
            "tool": "FilterIntervals",
            "mode": mode,
            "intervals": args.intervals,
            "count_samples": args.count_samples,
            "interval_set_rule": args.interval_set_rule,
            "threads": args.threads,
            "wall_seconds": elapsed,
            "intervals_per_second": args.intervals / elapsed if elapsed else 0.0,
            "output_intervals": summary["output_intervals"],
            "output_bytes": output.stat().st_size,
            "kernel": {
                "execution_space": telemetry["kernel_execution_space"],
                "records": telemetry["kernel_records"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["kernel_execute_seconds"],
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
