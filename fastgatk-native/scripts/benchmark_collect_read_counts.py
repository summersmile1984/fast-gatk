#!/usr/bin/env python3
"""File-boundary benchmark for the native CollectReadCounts TSV path."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_COLLECT_READ_COUNTS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-collect-read-counts"),
    ))
    records = int(os.environ.get("FASTGATK_READ_COUNTS_BENCH_RECORDS", "100000"))
    if records < 1:
        raise SystemExit("FASTGATK_READ_COUNTS_BENCH_RECORDS must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-read-counts-benchmark-") as directory:
        work = Path(directory)
        sam = work / "reads.sam"
        intervals = work / "targets.interval_list"
        output = work / "counts.tsv"
        manifest = work / "counts.tsv.manifest.json"
        with sam.open("w", encoding="utf-8") as stream:
            stream.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n")
            for index in range(records):
                start = (index % 100000) + 1
                stream.write(f"r{index}\t0\tchr1\t{start}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")
        intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000000\nchr1\t1\t100000\t+\tall\n",
            encoding="utf-8",
        )
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        start = time.perf_counter()
        completed = subprocess.run([
            str(native), "-I", str(sam), "-L", str(intervals), "-O", str(output),
            "--format", "TSV", "--sample", "BENCH", "--threads", "4", "--output-manifest", str(manifest)
        ], check=True, text=True, capture_output=True, env=env)
        wall_seconds = time.perf_counter() - start
        summary = json.loads(completed.stdout.splitlines()[-1])
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        if manifest_payload["telemetry"]["output_bytes"] != output.stat().st_size:
            raise AssertionError(manifest_payload)
        if manifest_payload["telemetry"]["reads_seen"] != records:
            raise AssertionError(manifest_payload)
        telemetry = manifest_payload["telemetry"]
        print(json.dumps({
            "status": "pass",
            "tool": "CollectReadCounts",
            "records": records,
            "wall_seconds": wall_seconds,
            "records_per_second": records / wall_seconds,
            "output_bytes": output.stat().st_size,
            "manifest_bytes": manifest.stat().st_size,
            "manifest_contract": True,
            "reads_seen": summary["reads_seen"],
            "reads_used": summary["reads_used"],
            "kernel": {
                "execution_space": telemetry["count_kernel_execution_space"],
                "batches": telemetry["count_kernel_batches"],
                "records": telemetry["count_kernel_records"],
                "prepare_seconds": telemetry["count_kernel_prepare_seconds"],
                "execute_seconds": telemetry["count_kernel_execute_seconds"],
            },
            "pipeline": {
                "lifecycle": telemetry["pipeline_lifecycle"],
                "decoded_items": telemetry["pipeline_decoded_items"],
                "computed_items": telemetry["pipeline_computed_items"],
                "encoded_items": telemetry["pipeline_encoded_items"],
                "peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
                "peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
                "peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
