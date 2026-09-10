#!/usr/bin/env python3
"""File-boundary throughput benchmark for GatherVcfs."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import tempfile
import time


HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length={length}>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--records-per-shard", type=int, default=10000)
    args = parser.parse_args()
    if args.shards < 1 or args.records_per_shard < 1:
        raise SystemExit("--shards and --records-per-shard must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_GATHER_VCFS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-gather-vcfs"),
    ))
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"missing executable: {binary}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-vcfs-benchmark-") as directory:
        work = pathlib.Path(directory)
        inputs: list[pathlib.Path] = []
        total_records = args.shards * args.records_per_shard
        for shard in range(args.shards):
            path = work / f"shard-{shard}.vcf.gz"
            with gzip.open(path, "wt", encoding="ascii") as stream:
                stream.write(HEADER.format(length=total_records + 1))
                first = shard * args.records_per_shard + 1
                for position in range(first, first + args.records_per_shard):
                    stream.write(f"chr1\t{position}\t.\tA\tG\t50\tPASS\t.\n")
            inputs.append(path)
        output = work / "gathered.vcf.gz"
        manifest = work / "gathered.manifest.json"
        command = [str(binary)]
        for path in inputs:
            command.extend(["-I", str(path)])
        command.extend(["--CREATE_INDEX", "true", "-O", str(output),
                        "--output-manifest", str(manifest)])
        begin = time.perf_counter()
        completed = subprocess.run(command, text=True, capture_output=True, check=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(completed.stdout.splitlines()[-1])
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert summary["input_records"] == total_records
        assert metadata["telemetry"]["input_records"] == total_records
        assert output.is_file() and pathlib.Path(f"{output}.tbi").is_file()
        print(json.dumps({
            "status": "pass",
            "tool": "GatherVcfs",
            "shards": args.shards,
            "records_per_shard": args.records_per_shard,
            "records": total_records,
            "wall_seconds": elapsed,
            "records_per_second": total_records / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "index_bytes": pathlib.Path(f"{output}.tbi").stat().st_size,
            "manifest_contract": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
