#!/usr/bin/env python3
"""File-boundary benchmark for native SplitIntervals."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "gatk-source/src/test/resources/hg19micro.fasta"
BINARY = Path(os.environ.get("FASTGATK_SPLIT_INTERVALS_BINARY", ROOT / "fastgatk-native/build/fastgatk-split-intervals"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--scatter-count", type=int, default=8)
    args = parser.parse_args()
    if not REFERENCE.is_file() or not BINARY.is_file():
        print(json.dumps({"status": "skip", "suite": "split-intervals-benchmark"}, sort_keys=True))
        return 0
    timings: list[float] = []
    output_bytes: list[int] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-split-intervals-bench-") as directory:
        root = Path(directory)
        for iteration in range(max(1, args.iterations)):
            outdir = root / f"run-{iteration}"
            start = time.perf_counter()
            result = subprocess.run(
                [str(BINARY), "-R", str(REFERENCE), "-O", str(outdir),
                 "--scatter-count", str(args.scatter_count)],
                text=True, capture_output=True,
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr)
            timings.append(time.perf_counter() - start)
            output_bytes.append(sum(path.stat().st_size for path in outdir.glob("*.interval_list")))
    timings.sort()
    p50 = timings[len(timings) // 2]
    p95 = timings[min(len(timings) - 1, int(len(timings) * 0.95))]
    print(json.dumps({
        "status": "pass", "suite": "split-intervals-benchmark", "backend": "native-host",
        "iterations": len(timings), "scatter_count": args.scatter_count,
        "p50_seconds": p50, "p95_seconds": p95,
        "output_bytes": output_bytes[len(output_bytes) // 2],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
