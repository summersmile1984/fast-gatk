#!/usr/bin/env python3
"""File-boundary benchmark for GatherPileupSummaries shard merging."""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = ROOT / "fastgatk-native" / "build" / "fastgatk-gather-pileup-summaries"
HEADER = "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"


def table(sample: str, rows: list[str]) -> str:
    return f"#<METADATA>SAMPLE={sample}\n" + HEADER + "".join(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=32)
    parser.add_argument("--rows", type=int, default=100)
    args = parser.parse_args()
    if args.shards < 1 or args.rows < 1:
        raise SystemExit("--shards and --rows must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-pileup-summaries-bench-") as temporary:
        work = pathlib.Path(temporary)
        dictionary = work / "ref.dict"
        dictionary.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n", encoding="utf-8",
        )
        inputs: list[pathlib.Path] = []
        for shard in range(args.shards):
            path = work / f"shard-{shard}.pileup"
            rows = [
                f"chr1\t{shard * args.rows + index + 1}\t10\t2\t0\t0.05\n"
                for index in range(args.rows)
            ]
            path.write_text(table("TUMOR", rows), encoding="utf-8")
            inputs.append(path)
        output = work / "gathered.pileup"
        manifest = work / "gathered.manifest.json"
        begin = time.perf_counter()
        result = subprocess.run(
            [str(BINARY), *sum((["-I", str(path)] for path in inputs), []),
             "-SD", str(dictionary), "-O", str(output),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in output.read_text(encoding="utf-8").splitlines()
                      if line and not line.startswith("#") and not line.startswith("contig"))
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        if not manifest_payload["outputs"][0]["complete"]:
            raise AssertionError(manifest_payload)
        if records != args.shards * args.rows:
            raise AssertionError({"expected_rows": args.shards * args.rows,
                                  "actual_rows": records})
        print(json.dumps({
            "status": "pass",
            "tool": "GatherPileupSummaries",
            "shards": args.shards,
            "input_rows_per_shard": args.rows,
            "gathered_rows": records,
            "wall_seconds": elapsed,
            "rows_per_second": args.shards * args.rows / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "manifest_bytes": manifest.stat().st_size,
            "manifest_contract": True,
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
