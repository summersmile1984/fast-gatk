#!/usr/bin/env python3
"""Reproducible file-boundary benchmark for CountBasesInReference."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def write_reference(path: pathlib.Path, fai: pathlib.Path, length: int) -> None:
    width = 80
    sequence = ("ACGTN" * ((length + 4) // 5))[:length]
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(">chr1\n")
        for offset in range(0, length, width):
            stream.write(sequence[offset:offset + width] + "\n")
    fai.write_text(f"chr1\t{length}\t6\t{width}\t{width + 1}\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", type=int, default=10_000_000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.bases < 1 or args.threads < 1:
        raise SystemExit("bases/threads must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_COUNT_BASES_BINARY",
        str(root / "fastgatk-native/build/fastgatk-count-bases-in-reference"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-count-bases-benchmark-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fa"
        write_reference(reference, work / "reference.fa.fai", args.bases)
        output = work / "counts.txt"
        command = [str(binary), "-R", str(reference), "-O", str(output),
                   "--threads", str(args.threads)]
        begin = time.perf_counter()
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(completed.stderr.splitlines()[-1])
        print(json.dumps({
            "status": "pass",
            "tool": "CountBasesInReference",
            "bases": args.bases,
            "threads": args.threads,
            "wall_seconds": elapsed,
            "bases_per_second": args.bases / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "kernel": {
                "execution_space": summary["kernel_execution_space"],
                "records": summary["kernel_records"],
                "execute_seconds": summary["kernel_execute_seconds"],
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
