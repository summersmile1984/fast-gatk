#!/usr/bin/env python3
"""File-boundary benchmark for the Kokkos DepthOfCoverage counting path."""

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
    "FASTGATK_DEPTH_OF_COVERAGE_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build")) /
        "fastgatk-depth-of-coverage"),
))
BAM = ROOT / "gatk-source" / "src" / "test" / "resources" / "NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source" / "src" / "test" / "resources" / "human_g1k_v37.chr17_1Mb.fasta"


def write_multisample_sam(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:17\tLN:1000000\n"
        "@RG\tID:rg1\tSM:S1\n"
        "@RG\tID:rg2\tSM:S2\n"
        "s1a\t0\t17\t69000\t60\t6M\t*\t0\t0\tAAAAAA\tIIIIII\tRG:Z:rg1\n"
        "s1b\t0\t17\t69000\t60\t6M\t*\t0\t0\tCCCCCC\tIIIIII\tRG:Z:rg1\n"
        "s2a\t0\t17\t69000\t60\t6M\t*\t0\t0\tGGGGGG\tIIIIII\tRG:Z:rg2\n"
        "s2b\t0\t17\t69001\t60\t6M\t*\t0\t0\tTTTTTT\tIIIIII\tRG:Z:rg2\n",
        encoding="utf-8",
    )


def write_fragment_sam(path: pathlib.Path) -> None:
    """Two paired fragments; mates are intentionally split across batches."""
    path.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:17\tLN:1000000\n"
        "@RG\tID:rgF\tSM:F_SAMPLE\n"
        "pairA\t99\t17\t69000\t60\t6M\t=\t69000\t6\tAAAAAA\tIIIIII\tRG:Z:rgF\n"
        "pairA\t147\t17\t69000\t60\t6M\t=\t69000\t-6\tCCCCCC\tIIIIII\tRG:Z:rgF\n"
        "pairB\t99\t17\t69000\t60\t6M\t=\t69000\t6\tGGGGGG\tIIIIII\tRG:Z:rgF\n"
        "pairB\t147\t17\t69000\t60\t6M\t=\t69000\t-6\tTTTTTT\tIIIIII\tRG:Z:rgF\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--batch-records", type=int, default=4096)
    args = parser.parse_args()
    if args.repeats < 1 or args.threads < 1 or args.batch_records < 1:
        raise SystemExit("repeats, threads and batch-records must be positive")
    if not BINARY.exists() or not BAM.exists() or not REFERENCE.exists():
        raise SystemExit("missing DepthOfCoverage benchmark assets")
    with tempfile.TemporaryDirectory(prefix="fastgatk-depth-of-coverage-bench-") as temporary:
        work = pathlib.Path(temporary)
        elapsed = []
        output_bytes = 0
        kernel_execute_seconds = []
        kernel_prepare_seconds = []
        kernel_batches = []
        kernel_observations = []
        for index in range(args.repeats):
            output = work / f"out.{index}"
            manifest = work / f"out.{index}.manifest.json"
            begin = time.perf_counter()
            subprocess.run(
                [str(BINARY), "-R", str(REFERENCE), "-I", str(BAM), "-L", "17:69000-70000",
                 "-O", str(output), "--threads", str(args.threads),
                 "--batch-records", str(args.batch_records), "--output-manifest", str(manifest)],
                check=True, text=True, capture_output=True,
            )
            elapsed.append(time.perf_counter() - begin)
            output_bytes += output.stat().st_size
            telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
            kernel_execute_seconds.append(telemetry["count_kernel_execute_seconds"])
            kernel_prepare_seconds.append(telemetry["count_kernel_prepare_seconds"])
            kernel_batches.append(telemetry["count_kernel_batches"])
            kernel_observations.append(telemetry["count_kernel_observations"])
        loci = 10001
        multi_sam = work / "multisample.sam"
        write_multisample_sam(multi_sam)
        multi_elapsed = []
        multi_output_bytes = 0
        multi_kernel_execute_seconds = []
        multi_kernel_prepare_seconds = []
        multi_kernel_batches = []
        multi_kernel_observations = []
        for index in range(args.repeats):
            output = work / f"multisample.{index}.out"
            manifest = work / f"multisample.{index}.manifest.json"
            begin = time.perf_counter()
            subprocess.run(
                [str(BINARY), "-R", str(REFERENCE), "-I", str(multi_sam), "-L", "17:69000-69005",
                 "-O", str(output), "--threads", str(args.threads),
                 "--batch-records", str(args.batch_records), "--output-manifest", str(manifest)],
                check=True, text=True, capture_output=True,
            )
            multi_elapsed.append(time.perf_counter() - begin)
            multi_output_bytes += output.stat().st_size
            telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
            multi_kernel_execute_seconds.append(telemetry["count_kernel_execute_seconds"])
            multi_kernel_prepare_seconds.append(telemetry["count_kernel_prepare_seconds"])
            multi_kernel_batches.append(telemetry["count_kernel_batches"])
            multi_kernel_observations.append(telemetry["count_kernel_observations"])
        fragment_sam = work / "fragments.sam"
        write_fragment_sam(fragment_sam)
        fragment_elapsed = []
        fragment_output_bytes = 0
        fragment_kernel_execute_seconds = []
        fragment_kernel_prepare_seconds = []
        fragment_kernel_batches = []
        fragment_kernel_observations = []
        fragment_keys = []
        for index in range(args.repeats):
            output = work / f"fragments.{index}.out"
            manifest = work / f"fragments.{index}.manifest.json"
            begin = time.perf_counter()
            subprocess.run(
                [str(BINARY), "-R", str(REFERENCE), "-I", str(fragment_sam), "-L", "17:69000-69005",
                 "-O", str(output), "--count-type", "COUNT_FRAGMENTS", "--threads", str(args.threads),
                 "--batch-records", "1", "--output-manifest", str(manifest)],
                check=True, text=True, capture_output=True,
            )
            fragment_elapsed.append(time.perf_counter() - begin)
            fragment_output_bytes += output.stat().st_size
            telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
            fragment_kernel_execute_seconds.append(telemetry["count_kernel_execute_seconds"])
            fragment_kernel_prepare_seconds.append(telemetry["count_kernel_prepare_seconds"])
            fragment_kernel_batches.append(telemetry["count_kernel_batches"])
            fragment_kernel_observations.append(telemetry["count_kernel_observations"])
            fragment_keys.append(telemetry["fragment_keys"])
        fragment_mean_wall = sum(fragment_elapsed) / len(fragment_elapsed)
        print(json.dumps({
            "status": "pass",
            "tool": "DepthOfCoverage",
            "repeats": args.repeats,
            "threads": args.threads,
            "batch_records": args.batch_records,
            "loci_per_run": loci,
            "mean_wall_seconds": sum(elapsed) / len(elapsed),
            "loci_per_second": loci / (sum(elapsed) / len(elapsed)),
            "output_bytes_total": output_bytes,
            "kernel": {
                "execution_space": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["count_kernel_execution_space"],
                "batches": kernel_batches[0],
                "observations": kernel_observations[0],
                "mean_prepare_seconds": sum(kernel_prepare_seconds) / len(kernel_prepare_seconds),
                "mean_execute_seconds": sum(kernel_execute_seconds) / len(kernel_execute_seconds),
            },
            "pipeline": {
                "decoded_items": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_decoded_items"],
                "computed_items": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_computed_items"],
                "encoded_items": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_encoded_items"],
                "peak_decoded_bytes": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_peak_decoded_bytes"],
                "peak_computed_bytes": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_peak_computed_bytes"],
                "peak_encoded_bytes": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["pipeline_peak_encoded_bytes"],
            },
            "multisample": {
                "sample_count": 2,
                "loci_per_run": 6,
                "mean_wall_seconds": sum(multi_elapsed) / len(multi_elapsed),
                "loci_per_second": 6 / (sum(multi_elapsed) / len(multi_elapsed)),
                "output_bytes_total": multi_output_bytes,
                "partition": "@RG-to-SM",
                "kernel": {
                    "batches": multi_kernel_batches[0],
                    "observations": multi_kernel_observations[0],
                    "mean_prepare_seconds": sum(multi_kernel_prepare_seconds) / len(multi_kernel_prepare_seconds),
                    "mean_execute_seconds": sum(multi_kernel_execute_seconds) / len(multi_kernel_execute_seconds),
                },
            },
            "count_fragments_extension": {
                "count_type": "COUNT_FRAGMENTS",
                "gatk_4_6_2_0_supported": False,
                "loci_per_run": 6,
                "mean_wall_seconds": fragment_mean_wall,
                "loci_per_second": 6 / fragment_mean_wall,
                "output_bytes_total": fragment_output_bytes,
                "fragment_keys": fragment_keys[0],
                "kernel": {
                    "batches": fragment_kernel_batches[0],
                    "observations": fragment_kernel_observations[0],
                    "mean_prepare_seconds": sum(fragment_kernel_prepare_seconds) / len(fragment_kernel_prepare_seconds),
                    "mean_execute_seconds": sum(fragment_kernel_execute_seconds) / len(fragment_kernel_execute_seconds),
                },
            },
            "execution_space": json.loads((work / "out.0.manifest.json").read_text(encoding="utf-8"))["telemetry"]["count_kernel_execution_space"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
