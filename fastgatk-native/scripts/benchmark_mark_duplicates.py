#!/usr/bin/env python3
"""File-boundary benchmark for bounded MarkDuplicates metadata/spill paths."""

from __future__ import annotations

import json
import os
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def make_sam(path: Path, records: int = 512) -> None:
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:100000",
        "@RG\tID:rg1\tSM:S1\tLB:lib1",
    ]
    for index in range(records):
        position = 10 + (index // 4) * 10
        quality = 40 if index % 4 == 0 else 20
        name = f"INST:1:1101:{index + 100}:{index + 100}"
        lines.append(
            f"{name}\t0\tchr1\t{position}\t{quality}\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(binary: Path, source: Path, work: Path) -> tuple[float, dict]:
    work.mkdir(parents=True, exist_ok=True)
    output = work / "output.sam"
    metrics = work / "output.metrics"
    manifest = work / "output.manifest.json"
    scratch = work / "scratch"
    start = time.perf_counter()
    result = subprocess.run(
        [str(binary), "-I", str(source), "-O", str(output),
         "--metrics-file", str(metrics), "--output-manifest", str(manifest),
         "--tagging-policy", "All", "--max-records-in-memory", "32",
         "--tmp-dir", str(scratch), "--create-output-bam-index=false"],
        text=True, capture_output=True, check=False,
        env={**os.environ, "OMP_PROC_BIND": "true", "OMP_PLACES": "threads"},
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return time.perf_counter() - start, json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-mark-duplicates"
    if not binary.is_file():
        raise SystemExit(f"missing binary: {binary}")
    repetitions = int(os.environ.get("FASTGATK_BENCH_REPEATS", "5"))
    if repetitions < 2 or repetitions > 100:
        raise SystemExit("FASTGATK_BENCH_REPEATS must be in [2,100]")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mark-duplicates-bench-") as directory:
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
        "benchmark": "mark-duplicates-file-boundary",
        "tool": "MarkDuplicates",
        "status": "pass",
        "schema_version": 1,
        "records": last["input_records"],
        "duplicate_records": last["duplicate_records"],
        "spill_run_count": last["spill_run_count"],
        "warmup_seconds": warmup_seconds,
        "p50_seconds": statistics.median(samples),
        "p95_seconds": ordered[p95_index],
        "repetitions": repetitions,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
