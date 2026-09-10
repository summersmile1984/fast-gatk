#!/usr/bin/env python3
"""File-boundary benchmark for the native BaseRecalibrator/ApplyBQSR loop.

The timings include HTSlib decode, Kokkos prepare/execute, BAM writing and
index creation.  This deliberately does not report a kernel-only speedup.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"


def run(command: list[str], env: dict[str, str]) -> tuple[float, dict]:
    started = time.perf_counter()
    result = subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    elapsed = time.perf_counter() - started
    return elapsed, json.loads(result.stdout.splitlines()[-1])


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    index = min(len(ordered) - 1, int(len(ordered) * fraction))
    return ordered[index]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--batch-records", type=int, default=64)
    args = parser.parse_args()
    if args.iterations < 1 or args.threads < 1 or args.batch_records < 1:
        raise SystemExit("iterations/threads/batch-records must be positive")
    # Keep the per-tool overrides for bisecting mixed builds, but make the
    # common FASTGATK_NATIVE_BUILD selector work like the other benchmark
    # drivers.  Without this fallback a requested Serial run silently used
    # the default OpenMP binaries, making backend comparisons misleading.
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
    bqsr = Path(os.environ.get("FASTGATK_BQSR_BINARY", build / "fastgatk-bqsr"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY", build / "fastgatk-apply-bqsr"))
    if not BAM.is_file() or not REFERENCE.is_file() or not bqsr.is_file() or not apply.is_file():
        print(json.dumps({"status": "skip", "suite": "bqsr-benchmark"}, sort_keys=True))
        return 0

    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    env["OMP_NUM_THREADS"] = str(args.threads)
    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-benchmark-") as directory:
        work = Path(directory)
        seed_report = work / "seed-recal.table"
        seed_manifest = work / "seed-recal.manifest.json"
        _, seed_summary = run([
            str(bqsr), "-I", str(BAM), "-R", str(REFERENCE), "-L", "17:69000-69100",
            "--batch-records", str(args.batch_records), "-O", str(seed_report),
            "--output-manifest", str(seed_manifest),
        ], env)

        phases = []
        recal_times: list[float] = []
        apply_times: list[float] = []
        last_recal_summary = seed_summary
        last_apply_summary: dict = {}
        for iteration in range(args.iterations):
            report = work / f"recal-{iteration}.table"
            manifest = work / f"recal-{iteration}.manifest.json"
            elapsed, summary = run([
                str(bqsr), "-I", str(BAM), "-R", str(REFERENCE), "-L", "17:69000-69100",
                "--batch-records", str(args.batch_records), "-O", str(report),
                "--output-manifest", str(manifest),
            ], env)
            recal_times.append(elapsed)
            last_recal_summary = summary
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            phases.append({
                "phase": "BaseRecalibrator", "iteration": iteration,
                "wall_seconds": elapsed, "output_bytes": report.stat().st_size,
                "kernel_execution_space": metadata["telemetry"]["kernel_execution_space"],
                "kernel_execution_policy": metadata["telemetry"]["kernel_execution_policy"],
                "kernel_prepare_seconds": metadata["telemetry"]["kernel_prepare_seconds"],
                "kernel_execute_seconds": metadata["telemetry"]["kernel_execute_seconds"],
                "pipeline_lifecycle": metadata["telemetry"]["pipeline_lifecycle"],
                "pipeline_decoded_items": metadata["telemetry"]["pipeline_decoded_items"],
                "pipeline_computed_items": metadata["telemetry"]["pipeline_computed_items"],
                "pipeline_encoded_items": metadata["telemetry"]["pipeline_encoded_items"],
                "pipeline_decoded_bytes": metadata["telemetry"]["pipeline_decoded_bytes"],
                "pipeline_computed_bytes": metadata["telemetry"]["pipeline_computed_bytes"],
                "pipeline_encoded_bytes": metadata["telemetry"]["pipeline_encoded_bytes"],
                "pipeline_peak_decoded_bytes": metadata["telemetry"]["pipeline_peak_decoded_bytes"],
                "pipeline_peak_computed_bytes": metadata["telemetry"]["pipeline_peak_computed_bytes"],
                "pipeline_peak_encoded_bytes": metadata["telemetry"]["pipeline_peak_encoded_bytes"],
            })

            output = work / f"recal-{iteration}.bam"
            output_manifest = work / f"recal-{iteration}.bam.manifest.json"
            elapsed, summary = run([
                str(apply), "-I", str(BAM), "-R", str(REFERENCE),
                "--bqsr-recal-file", str(report), "--batch-records", str(args.batch_records),
                "-O", str(output),
                "--output-manifest", str(output_manifest),
            ], env)
            apply_times.append(elapsed)
            last_apply_summary = summary
            metadata = json.loads(output_manifest.read_text(encoding="utf-8"))
            phases.append({
                "phase": "ApplyBQSR", "iteration": iteration,
                "wall_seconds": elapsed, "output_bytes": output.stat().st_size,
                "index_bytes": Path(str(output) + ".bai").stat().st_size,
                "kernel_execution_space": metadata["telemetry"]["kernel_execution_space"],
                "kernel_execution_policy": metadata["telemetry"]["kernel_execution_policy"],
                "kernel_prepare_seconds": metadata["telemetry"]["kernel_prepare_seconds"],
                "kernel_execute_seconds": metadata["telemetry"]["kernel_execute_seconds"],
                "pipeline_lifecycle": metadata["telemetry"]["pipeline_lifecycle"],
                "pipeline_decoded_items": metadata["telemetry"]["pipeline_decoded_items"],
                "pipeline_computed_items": metadata["telemetry"]["pipeline_computed_items"],
                "pipeline_encoded_items": metadata["telemetry"]["pipeline_encoded_items"],
                "pipeline_decoded_bytes": metadata["telemetry"]["pipeline_decoded_bytes"],
                "pipeline_computed_bytes": metadata["telemetry"]["pipeline_computed_bytes"],
                "pipeline_encoded_bytes": metadata["telemetry"]["pipeline_encoded_bytes"],
                "pipeline_peak_decoded_bytes": metadata["telemetry"]["pipeline_peak_decoded_bytes"],
                "pipeline_peak_computed_bytes": metadata["telemetry"]["pipeline_peak_computed_bytes"],
                "pipeline_peak_encoded_bytes": metadata["telemetry"]["pipeline_peak_encoded_bytes"],
            })

        print(json.dumps({
            "status": "pass", "suite": "bqsr-benchmark", "threads": args.threads,
            "batch_records": args.batch_records, "iterations": args.iterations,
            "input_records": last_apply_summary.get("records", 0),
            "input_bytes": BAM.stat().st_size,
            "summary": {
                "BaseRecalibrator": {
                    "p50_seconds": percentile(recal_times, 0.50),
                    "p95_seconds": percentile(recal_times, 0.95),
                    "reports_per_second": len(recal_times) / sum(recal_times),
                },
                "ApplyBQSR": {
                    "p50_seconds": percentile(apply_times, 0.50),
                    "p95_seconds": percentile(apply_times, 0.95),
                    "bam_per_second": len(apply_times) / sum(apply_times),
                },
            },
            "phases": phases,
            "last_recal_summary": last_recal_summary,
            "last_apply_summary": last_apply_summary,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
