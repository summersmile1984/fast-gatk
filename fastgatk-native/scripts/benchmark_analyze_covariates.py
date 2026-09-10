#!/usr/bin/env python3
"""File-boundary benchmark for the native AnalyzeCovariates report path."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def write_report(path: Path, rows: int) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write("#:GATKTable:Arguments:\nArgument Value\nquantizing_levels 16\n")
        stream.write("#:GATKTable:RecalTable1:\n")
        stream.write("ReadGroup QualityScore EventType EmpiricalQuality Observations Errors\n")
        for index in range(rows):
            stream.write(f"rg{index % 8} {20 + index % 20} M {20 + index % 20:.2f} 100 1\n")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = root / "fastgatk-native/build/fastgatk-analyze-covariates"
    rows = int(os.environ.get("FASTGATK_ANALYZE_COVARIATES_BENCH_ROWS", "10000"))
    if rows < 1:
        raise SystemExit("FASTGATK_ANALYZE_COVARIATES_BENCH_ROWS must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-analyze-covariates-benchmark-") as directory:
        work = Path(directory)
        report = work / "report.table"
        csv = work / "analysis.csv"
        manifest = work / "analysis.manifest.json"
        write_report(report, rows)
        begin = time.perf_counter()
        result = subprocess.run(
            [str(native), "-BQSR", str(report), "-csv", str(csv),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=True,
        )
        wall_seconds = time.perf_counter() - begin
        summary = json.loads(result.stdout.splitlines()[-1])
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        output_rows = metadata["telemetry"]["rows"]
        if output_rows < 1 or summary["rows"] != output_rows:
            raise AssertionError(metadata)
        if metadata["telemetry"]["csv_bytes"] != csv.stat().st_size:
            raise AssertionError(metadata)
        print(json.dumps({
            "status": "pass",
            "tool": "AnalyzeCovariates",
            "input_rows": rows,
            "rows": output_rows,
            "output_bytes": csv.stat().st_size,
            "manifest_bytes": manifest.stat().st_size,
            "manifest_contract": True,
            "wall_seconds": wall_seconds,
            "rows_per_second": rows / wall_seconds if wall_seconds else 0.0,
            "summary": summary,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
