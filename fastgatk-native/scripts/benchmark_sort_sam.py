#!/usr/bin/env python3
"""File-boundary benchmark for bounded external-memory SortSam."""

from __future__ import annotations

import json
import os
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def make_sam(path: Path, records: int = 512) -> None:
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:100000"]
    # Reverse the input order so the benchmark exercises spill sorting and
    # k-way merge rather than only the already-sorted fast path.
    for index in reversed(range(records)):
        position = 10 + (index % 256) * 10
        lines.append(
            f"read-{index:06d}\t0\tchr1\t{position}\t60\t4M\t*\t0\t0\tACGT\tIIII"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(binary: Path, source: Path, work: Path) -> tuple[float, dict]:
    work.mkdir(parents=True, exist_ok=True)
    output = work / "output.sam"
    manifest = work / "output.manifest.json"
    scratch = work / "scratch"
    start = time.perf_counter()
    result = subprocess.run(
        [str(binary), "-I", str(source), "-O", str(output),
         "--sort-order", "coordinate", "--max-records-in-memory", "32",
         "--tmp-dir", str(scratch), "--create-output-bam-index=false",
         "--output-manifest", str(manifest)],
        text=True, capture_output=True, check=False,
        env={**os.environ, "OMP_PROC_BIND": "true", "OMP_PLACES": "threads"},
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return time.perf_counter() - start, json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-sort-sam"
    if not binary.is_file():
        raise SystemExit(f"missing binary: {binary}")
    repetitions = int(os.environ.get("FASTGATK_BENCH_REPEATS", "5"))
    if repetitions < 2 or repetitions > 100:
        raise SystemExit("FASTGATK_BENCH_REPEATS must be in [2,100]")
    with tempfile.TemporaryDirectory(prefix="fastgatk-sort-sam-bench-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        make_sam(source)
        warmup_seconds, warmup = run(binary, source, work / "warmup")
        samples = []
        last = warmup
        for index in range(repetitions):
            elapsed, last = run(binary, source, work / f"run-{index}")
            samples.append(elapsed)
    ordered = sorted(samples)
    p95_index = min(len(ordered) - 1, int(round(0.95 * (len(ordered) - 1))))
    print(json.dumps({
        "benchmark": "sort-sam-file-boundary",
        "tool": "SortSam",
        "status": "pass",
        "schema_version": 1,
        "records": last["output_records"],
        "spill_runs": last["spill_runs"],
        "warmup_seconds": warmup_seconds,
        "p50_seconds": statistics.median(samples),
        "p95_seconds": ordered[p95_index],
        "repetitions": repetitions,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
