#!/usr/bin/env python3
"""File-boundary benchmark for VariantFiltration genotype FT handling."""

from __future__ import annotations

import gzip
import json
import pathlib
import statistics
import subprocess
import tempfile
import time


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=20000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-variant-filtration"
    if not binary.is_file():
        raise SystemExit("native VariantFiltration binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-filtration-benchmark-") as directory:
        work = pathlib.Path(directory)
        source = work / "input.vcf.gz"
        records = 10000
        lines = [HEADER]
        for index in range(records):
            position = index + 1
            depth = 5 if index % 10 == 0 else 20
            existing = "Prior" if index % 3 == 0 else "PASS"
            lines.append(
                f"chr1\t{position}\tvar{index}\tA\tG\t50\tPASS\t.\tGT:DP:FT\t"
                f"0/1:{depth}:{existing}\t0/0:20:Existing2\n"
            )
        with gzip.open(source, "wt", encoding="ascii") as stream:
            stream.write("".join(lines))
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        command = [
            str(binary), "-V", str(source), "-O", str(output),
            "--genotype-filter-expression", "DP < 10",
            "--genotype-filter-name", "LowDP",
            "--output-manifest", str(manifest),
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
            "suite": "fastgatk-variant-filtration-file-boundary",
            "tool": "VariantFiltration",
            "backend": "HTSlib+Host",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "records": records,
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": source.stat().st_size,
            "output_bytes": output.stat().st_size,
            "genotype_filtered": metadata["telemetry"]["genotype_filtered"],
            "records_per_second": records / statistics.median(samples),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
