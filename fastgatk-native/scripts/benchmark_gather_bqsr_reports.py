#!/usr/bin/env python3
"""File-boundary benchmark for GatherBQSRReports GATKReport merging.

Generates synthetic but schema-valid GATKReport v1.1 shards (Arguments,
Quantized, RecalTable0/1/2 with the GATK table headers emitted by the native
writer), then measures warmup/p50/p95 of the native merge across a bounded
shard set.  This is the K3 file-boundary layer; Java comparison belongs to
verify_gather_bqsr_gatk_oracle.py.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = ROOT / "fastgatk-native" / "build" / "fastgatk-gather-bqsr-reports"


def write_shard(path: pathlib.Path, shard: int, rows: int) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#:GATKReport.v1.1:5\n")
        handle.write("#:GATKTable:2:3:%s:%s:;\n")
        handle.write("#:GATKTable:Arguments:Recalibration argument collection values used in this run\n")
        handle.write("Argument\tValue\n")
        handle.write("binary_tag_name\tnull\n")
        handle.write("covariate\tReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate\n")
        handle.write("quantizing_levels\t16\n")
        handle.write("\n")
        handle.write("#:GATKTable:3:94:%d:%d:%d:;\n")
        handle.write("#:GATKTable:Quantized:Quality quantization map\n")
        handle.write("QualityScore\tCount\tQuantizedScore\n")
        for quality in range(94):
            handle.write(f"{quality}\t{100 + shard}\t{min(quality, 93)}\n")
        handle.write("\n")
        handle.write(f"#:GATKTable:6:{rows}:%s:%s:%.4f:%.4f:%d:%.2f:;\n")
        handle.write("#:GATKTable:RecalTable0:\n")
        handle.write("ReadGroup\tEventType\tEmpiricalQuality\tEstimatedQReported\tObservations\tErrors\n")
        for row in range(rows):
            handle.write(f"RG{shard}_{row}\tM\t{20.0 + (row % 20):.4f}\t{20.0 + (row % 20):.4f}\t{1000 + row}\t{1 + (row % 7):.2f}\n")
        handle.write("\n")
        handle.write(f"#:GATKTable:6:{rows}:%s:%s:%s:%.4f:%d:%.2f:;\n")
        handle.write("#:GATKTable:RecalTable1:\n")
        handle.write("ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations\tErrors\n")
        for row in range(rows):
            handle.write(f"RG{shard}_{row % 4}\t{10 + row % 30}\tM\t{10.0 + row % 30:.4f}\t{500 + row}\t{row % 5:.2f}\n")
        handle.write("\n")
        handle.write(f"#:GATKTable:8:{rows}:%s:%s:%s:%s:%s:%.4f:%d:%.2f:;\n")
        handle.write("#:GATKTable:RecalTable2:\n")
        handle.write("ReadGroup\tQualityScore\tCovariateValue\tCovariateName\tEventType\tEmpiricalQuality\tObservations\tErrors\n")
        for row in range(rows):
            handle.write(f"RG{shard}_{row % 4}\t{10 + row % 30}\tAC\tContext\tM\t{10.0 + row % 30:.4f}\t{300 + row}\t{row % 3:.2f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=16)
    parser.add_argument("--rows", type=int, default=64)
    parser.add_argument("--repetitions", type=int, default=5)
    args = parser.parse_args()
    if args.shards < 1 or args.rows < 1 or args.repetitions < 1:
        raise SystemExit("--shards, --rows and --repetitions must be positive")
    if not BINARY.is_file():
        raise SystemExit(f"missing native binary: {BINARY}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-bqsr-bench-") as temporary:
        work = pathlib.Path(temporary)
        inputs: list[pathlib.Path] = []
        for shard in range(args.shards):
            path = work / f"shard-{shard}.table"
            write_shard(path, shard, args.rows)
            inputs.append(path)

        def one_run(tag: str) -> float:
            output = work / f"gathered-{tag}.table"
            manifest = work / f"gathered-{tag}.manifest.json"
            begin = time.perf_counter()
            subprocess.run(
                [str(BINARY), *sum((["-I", str(path)] for path in inputs), []),
                 "-O", str(output), "--output-manifest", str(manifest)],
                text=True, capture_output=True, check=True)
            elapsed = time.perf_counter() - begin
            manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
            if manifest_payload["primary_output"] != str(output):
                raise AssertionError(manifest_payload)
            if not manifest_payload["outputs"][0]["complete"]:
                raise AssertionError(manifest_payload)
            return elapsed

        warmup = one_run("warmup")
        timings = sorted(one_run(f"repeat-{index}") for index in range(args.repetitions))
        p50 = timings[len(timings) // 2]
        p95 = timings[min(len(timings) - 1, int(len(timings) * 0.95 + 0.999999) - 1)]
        total_rows = args.shards * args.rows
        print(json.dumps({
            "status": "pass",
            "benchmark": "gather-bqsr-reports-file-boundary",
            "tool": "GatherBQSRReports",
            "schema_version": 1,
            "shards": args.shards,
            "rows_per_shard": args.rows,
            "total_rows": total_rows,
            "repetitions": args.repetitions,
            "warmup_seconds": warmup,
            "p50_seconds": p50,
            "p95_seconds": p95,
            "rows_per_second_p50": total_rows / p50 if p50 else 0.0,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
