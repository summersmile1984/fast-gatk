#!/usr/bin/env python3
"""File-boundary benchmark for the native VariantsToTable extraction path."""

from __future__ import annotations

import json
import os
import pathlib
import statistics
import subprocess
import tempfile
import time


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-variants-to-table"
    if not binary.is_file():
        raise SystemExit("native VariantsToTable binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variants-to-table-benchmark-") as directory:
        work = pathlib.Path(directory)
        source = work / "input.vcf"
        records = 10000
        lines = [HEADER]
        for index in range(records):
            pos = index + 1
            if index % 11 == 0:
                alt, ac, ad = "G,T", "2,3", "8,4,2"
            else:
                alt, ac, ad = "G", "1", "8,4"
            lines.append(
                f"chr1\t{pos}\tvar{index}\tA\t{alt}\t50\tPASS\tDP=12;AC={ac}\tGT:AD\t0/1:{ad}\n"
            )
        source.write_text("".join(lines), encoding="ascii")
        output = work / "table.tsv"
        manifest = work / "table.manifest.json"
        command = [
            str(binary), "-V", str(source), "-O", str(output),
            "-F", "CHROM", "-F", "POS", "-F", "ALT", "-F", "TYPE", "-F", "AC",
            "-GF", "GT", "-GF", "AD", "--output-manifest", str(manifest),
        ]
        samples: list[float] = []
        for iteration in range(6):
            started = time.perf_counter()
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            elapsed = time.perf_counter() - started
            if completed.returncode != 0:
                raise SystemExit(completed.stderr or completed.stdout)
            if iteration:
                samples.append(elapsed)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        ordered = sorted(samples)
        p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]
        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-variants-to-table-file-boundary",
            "tool": "VariantsToTable",
            "backend": "HTSlib+Host",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "records": records,
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": source.stat().st_size,
            "output_bytes": output.stat().st_size,
            "output_rows": metadata["telemetry"]["output_rows"],
            "records_per_second": records / statistics.median(samples),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
