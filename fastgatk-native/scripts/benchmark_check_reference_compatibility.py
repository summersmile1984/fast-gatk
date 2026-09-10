#!/usr/bin/env python3
"""Reproducible dictionary-boundary benchmark for CheckReferenceCompatibility."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import subprocess
import tempfile
import time


def write_inputs(reference: pathlib.Path, query: pathlib.Path, contigs: int, length: int) -> None:
    width = 80
    records: list[tuple[str, int, int, str]] = []
    offset = 0
    with reference.open("wb") as fasta:
        for index in range(contigs):
            name = f"chr{index + 1}"
            sequence = ("ACGT" * ((length + 3) // 4))[:length]
            digest = hashlib.md5(sequence.encode("ascii")).hexdigest()
            header = f">{name}\n".encode("ascii")
            fasta.write(header)
            offset += len(header)
            records.append((name, length, offset, digest))
            for begin in range(0, length, width):
                line = (sequence[begin:begin + width] + "\n").encode("ascii")
                fasta.write(line)
                offset += len(line)
    with reference.with_suffix(reference.suffix + ".fai").open("w", encoding="ascii") as fai:
        # The benchmark uses fixed 80-base lines, so each sequence offset is
        # recovered from the FASTA layout above by a second deterministic pass.
        offset = 0
        for name, seq_length, _, _ in records:
            offset += len(f">{name}\n".encode("ascii"))
            fai.write(f"{name}\t{seq_length}\t{offset}\t{width}\t{width + 1}\n")
            offset += ((seq_length + width - 1) // width) * (width + 1)
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n" + "".join(
            f"@SQ\tSN:{name}\tLN:{length}\tM5:{digest}\n" for name, length, _, digest in records
        ), encoding="ascii")
    query.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
            f"@SQ\tSN:{name}\tLN:{length}\tM5:{digest}\n" for name, length, _, digest in records
        ), encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contigs", type=int, default=128)
    parser.add_argument("--length", type=int, default=1000)
    args = parser.parse_args()
    if args.contigs < 1 or args.length < 1:
        raise SystemExit("contigs/length must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_CHECK_REFERENCE_COMPATIBILITY_BINARY",
        str(root / "fastgatk-native/build/fastgatk-check-reference-compatibility"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-check-reference-benchmark-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fa"
        query = work / "query.sam"
        write_inputs(reference, query, args.contigs, args.length)
        output = work / "compatibility.table"
        begin = time.perf_counter()
        subprocess.run([str(binary), "-I", str(query), "-refcomp", str(reference), "-O", str(output)],
                       check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        print(json.dumps({
            "status": "pass",
            "tool": "CheckReferenceCompatibility",
            "contigs": args.contigs,
            "bases": args.contigs * args.length,
            "wall_seconds": elapsed,
            "contigs_per_second": args.contigs / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "mode": "Host dictionary IO",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
