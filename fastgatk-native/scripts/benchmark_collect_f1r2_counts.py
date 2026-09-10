#!/usr/bin/env python3
"""Small file-boundary benchmark for native CollectF1R2Counts."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_COLLECT_F1R2_COUNTS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-collect-f1r2-counts")))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    repeats = max(1, int(os.environ.get("FASTGATK_F1R2_BENCH_REPEATS", "3")))
    batch_records = max(1, int(os.environ.get("FASTGATK_F1R2_BATCH_RECORDS", "4096")))
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-f1r2-bench-") as directory:
        work = Path(directory)
        summaries = []
        samples = []
        unindexed_bam = work / "unindexed.bam"
        shutil.copyfile(bam, unindexed_bam)
        sequential_samples = []
        for index in range(repeats):
            output = work / f"native-{index}.tar.gz"
            manifest = work / f"native-{index}.manifest.json"
            start = time.perf_counter()
            result = subprocess.run([
                str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
                "-O", str(output), "--output-manifest", str(manifest), "--threads", "2",
                "--batch-records", str(batch_records),
            ], env=env, text=True, capture_output=True, check=True)
            elapsed = time.perf_counter() - start
            summary = json.loads(result.stdout.splitlines()[-1])
            samples.append(elapsed)
            summaries.append(summary)
            sequential_output = work / f"sequential-{index}.tar.gz"
            start = time.perf_counter()
            sequential_result = subprocess.run([
                str(native), "-R", str(reference), "-I", str(unindexed_bam), "-L", "17:69000-70000",
                "-O", str(sequential_output), "--threads", "2", "--batch-records", str(batch_records),
            ], env=env, text=True, capture_output=True, check=True)
            sequential_elapsed = time.perf_counter() - start
            sequential_summary = json.loads(sequential_result.stdout.splitlines()[-1])
            assert sequential_summary["indexed_inputs"] == 0
            assert sequential_summary["sequential_inputs"] == 1
            sequential_samples.append(sequential_elapsed)
        summary = summaries[-1]
        telemetry = json.loads((work / f"native-{repeats - 1}.manifest.json").read_text())["telemetry"]
        print(json.dumps({
            "status": "pass", "tool": "CollectF1R2Counts", "backend": "Kokkos",
            "repeats": repeats, "records_seen": summary["records_seen"],
            "observations": summary["observations"], "loci": summary["loci"],
            "accepted_loci": summary["accepted_loci"],
            "batch_records": batch_records,
            "wall_seconds": min(samples), "p50_wall_seconds": sorted(samples)[len(samples) // 2],
            "records_per_second": summary["records_seen"] / min(samples),
            "indexed_traversal": summary["indexed_inputs"] == 1,
            "indexed_p50_wall_seconds": sorted(samples)[len(samples) // 2],
            "sequential_p50_wall_seconds": sorted(sequential_samples)[len(sequential_samples) // 2],
            "kernel": {
                "execution_space": telemetry["kernel_execution_space"],
                "batches": telemetry["kernel_batches"],
                "observations": telemetry["kernel_observations"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["kernel_execute_seconds"],
                "pipeline_decoded_items": telemetry["pipeline_decoded_items"],
                "pipeline_computed_items": telemetry["pipeline_computed_items"],
                "pipeline_encoded_items": telemetry["pipeline_encoded_items"],
                "pipeline_peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
                "pipeline_peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
                "pipeline_peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
            },
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
