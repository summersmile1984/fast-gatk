#!/usr/bin/env python3
"""File-boundary benchmark for the deterministic ModelSegments prototype."""

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
    "FASTGATK_MODEL_SEGMENTS_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-model-segments")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=100_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--kernel", action="store_true",
                        help="benchmark the deterministic Kokkos KernelSegmenter path")
    parser.add_argument("--max-segments", type=int, default=8)
    parser.add_argument("--window-size", type=int, action="append", default=[])
    parser.add_argument("--probabilistic", action="store_true",
                        help="benchmark the bounded deterministic Kokkos posterior sampler")
    parser.add_argument("--num-samples", type=int, default=256)
    parser.add_argument("--num-burn-in-iterations", type=int, default=32)
    parser.add_argument("--multisample", action="store_true",
                        help="benchmark the bounded repeated-input joint segmentation path")
    parser.add_argument("--matched-normal", action="store_true",
                        help="benchmark matched-normal heterozygous-site genotyping/filtering")
    args = parser.parse_args()
    if args.points < 1 or args.threads < 1:
        raise SystemExit("--points and --threads must be positive")
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-bench-") as temporary:
        work = pathlib.Path(temporary)
        input_path = work / "denoised.tsv"
        output_path = work / ("sample.interval_list" if args.multisample else "sample.modelFinal.segments.tsv")
        manifest_path = work / "sample.modelFinal.segments.json"
        with input_path.open("w", encoding="utf-8") as handle:
            handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000000\n")
            handle.write("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n")
            for index in range(args.points):
                start = index * 10 + 1
                value = ((index % 200) - 100) / 10000.0
                handle.write(f"chr1\t{start}\t{start + 9}\t1\t{value}\n")
        case_allelic = work / "case.allelic.tsv"
        normal_allelic = work / "normal.allelic.tsv"
        if args.matched_normal:
            for path, is_normal in ((case_allelic, False), (normal_allelic, True)):
                with path.open("w", encoding="utf-8") as handle:
                    handle.write("@HD\tVN:1.6\n@RG\tID:GATKCopyNumber\tSM:{}\n".format(
                        "NORMAL" if is_normal else "CASE"))
                    handle.write("CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n")
                    for index in range(args.points):
                        start = index * 10 + 1
                        # Periodic homozygous loci exercise the normal filter;
                        # remaining rows are balanced enough for the -10 log
                        # ratio threshold to classify them as heterozygous.
                        alt = 0 if is_normal and index % 17 == 0 else (8 if index % 23 == 0 else 10)
                        handle.write(f"chr1\t{start}\t10\t{alt}\tA\tC\n")
        begin = time.perf_counter()
        command = [str(BINARY), "--denoised-copy-ratios", str(input_path),
                   "--output-prefix", str(work / "sample"), "--output-manifest", str(manifest_path),
                   "--threads", str(args.threads)]
        if args.multisample:
            command.extend(["--denoised-copy-ratios", str(input_path)])
        if args.matched_normal:
            command.extend(["--allelic-counts", str(case_allelic),
                            "--normal-allelic-counts", str(normal_allelic),
                            "--minimum-total-allele-count-normal", "10"])
        if args.kernel:
            command.extend(["--maximum-number-of-segments-per-chromosome", str(args.max_segments),
                            "--kernel-approximation-dimension", "100",
                            "--number-of-changepoints-penalty-factor", "0"])
            for window in args.window_size or [8, 16, 32]:
                command.extend(["--window-size", str(window)])
        if args.probabilistic:
            if args.num_samples < 1 or args.num_burn_in_iterations < 0:
                raise SystemExit("--num-samples must be positive and --num-burn-in-iterations non-negative")
            command.extend(["--mode", "COPY_RATIO",
                            "--num-samples", str(args.num_samples),
                            "--num-burn-in-iterations", str(args.num_burn_in_iterations)])
        result = subprocess.run(
            command,
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        telemetry = manifest_data["telemetry"]
        records = sum(1 for line in output_path.read_text(encoding="utf-8").splitlines()
                      if line and not line.startswith("@") and not line.startswith("CONTIG"))
        print(json.dumps({
            "status": "pass",
            "tool": "ModelSegments",
            "segmentation_method": "KernelSegmenter" if args.kernel else "threshold",
            "multisample": args.multisample,
            "matched_normal": args.matched_normal,
            "points": args.points,
            "segments": records,
            "wall_seconds": elapsed,
            "points_per_second": args.points / elapsed if elapsed else 0.0,
            "output_bytes": output_path.stat().st_size,
            "probabilistic": args.probabilistic,
            "num_samples": args.num_samples if args.probabilistic else 0,
            "num_burn_in_iterations": args.num_burn_in_iterations if args.probabilistic else 0,
            "kernel": {
                "execution_space": telemetry.get("execution_space", ""),
                "policy": telemetry.get("kernel_execution_policy", ""),
                "batches": telemetry.get("kernel_batches", 0),
                "observations": telemetry.get("kernel_observations", 0),
                "prepare_seconds": telemetry.get("kernel_prepare_seconds", 0.0),
                "execute_seconds": telemetry.get("kernel_execute_seconds", 0.0),
            },
            "manifest": manifest_data,
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
