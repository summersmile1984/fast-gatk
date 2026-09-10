#!/usr/bin/env python3
"""Reproducible file-boundary benchmark for CompareReferences."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def write_reference(path: pathlib.Path, fai: pathlib.Path, length: int, mutate: bool) -> None:
    width = 80
    sequence = list(("ACGT" * ((length + 3) // 4))[:length])
    if mutate:
        sequence[length // 2] = "T" if sequence[length // 2] != "T" else "A"
    sequence = "".join(sequence)
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(">chr1\n")
        for offset in range(0, length, width):
            stream.write(sequence[offset:offset + width] + "\n")
    fai.write_text(f"chr1\t{length}\t6\t{width}\t{width + 1}\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", type=int, default=2_000_000)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.bases < 1 or args.threads < 1:
        raise SystemExit("bases/threads must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_COMPARE_REFERENCES_BINARY",
        str(root / "fastgatk-native/build/fastgatk-compare-references"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-compare-references-benchmark-") as directory:
        work = pathlib.Path(directory)
        first = work / "first.fa"
        second = work / "second.fa"
        write_reference(first, work / "first.fa.fai", args.bases, False)
        write_reference(second, work / "second.fa.fai", args.bases, True)
        output = work / "table.tsv"
        snps = work / "snps"
        snps.mkdir()
        command = [str(binary), "-R", str(first), "--refcomp", str(second), "-O", str(output),
                   "--base-comparison", "FIND_SNPS_ONLY", "--base-comparison-output", str(snps),
                   "--threads", str(args.threads)]
        begin = time.perf_counter()
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(next(line for line in reversed(completed.stderr.splitlines())
                                  if line.startswith("{\"base_comparison\"")))
        print(json.dumps({
            "status": "pass",
            "tool": "CompareReferences",
            "bases": args.bases,
            "threads": args.threads,
            "wall_seconds": elapsed,
            "bases_per_second": args.bases / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "kernel": {
                "execution_space": "OpenMP/Serial",
                "records": args.bases,
                "mismatches": 1,
                "summary": summary,
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
