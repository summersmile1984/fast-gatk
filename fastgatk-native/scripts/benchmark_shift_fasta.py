#!/usr/bin/env python3
"""Small reproducible throughput benchmark for native ShiftFasta."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "gatk-source/src/test/resources/exampleFASTA.fasta"
BIN = ROOT / "fastgatk-native/build/fastgatk-shift-fasta"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repetitions", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.repetitions < 1 or args.threads < 1:
        raise SystemExit("repetitions and threads must be positive")
    assert BIN.exists() and REF.exists()
    bases = sum(len(line.strip()) for line in REF.read_text().splitlines()
                if not line.startswith(">"))
    elapsed: list[float] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-bench-shift-") as directory:
        work = Path(directory)
        for index in range(args.repetitions):
            output = work / f"shifted-{index}.fasta"
            chain = work / f"shifted-{index}.chain"
            start = time.perf_counter()
            subprocess.run([str(BIN), "-R", str(REF), "-O", str(output),
                            "--shift-back-output", str(chain), "--threads", str(args.threads)],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            elapsed.append(time.perf_counter() - start)
    p50 = statistics.median(elapsed)
    p95 = sorted(elapsed)[max(0, int(len(elapsed) * 0.95 + 0.999) - 1)]
    print(json.dumps({"tool": "ShiftFasta", "repetitions": args.repetitions,
                      "threads": args.threads, "bases": bases,
                      "p50_seconds": p50, "p95_seconds": p95,
                      "p50_bases_per_second": bases / p50 if p50 else 0.0,
                      "p95_bases_per_second": bases / p95 if p95 else 0.0}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
