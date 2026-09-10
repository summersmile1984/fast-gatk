#!/usr/bin/env python3
"""File-boundary benchmark for ApplyVQSR recal-VCF lookup and filtering."""

from __future__ import annotations

import argparse
import gzip
import json
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = ROOT / "fastgatk-native" / "build" / "fastgatk-apply-vqsr"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=10_000)
    args = parser.parse_args()
    if args.records < 1:
        raise SystemExit("--records must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-bench-") as temporary:
        work = pathlib.Path(temporary)
        input_vcf = work / "input.vcf"
        recal_vcf = work / "recal.vcf"
        output = work / "output.vcf.gz"
        with input_vcf.open("w", encoding="utf-8") as data, recal_vcf.open("w", encoding="utf-8") as recal:
            data.write(HEADER)
            recal.write(HEADER)
            for index in range(args.records):
                position = index * 10 + 1
                score = 1.0 if index % 10 else -1.0
                data.write(f"chr1\t{position}\t.\tA\tG\t50\tPASS\t.\n")
                recal.write(f"chr1\t{position}\t.\tA\tG\t.\tPASS\tVQSLOD={score}\n")
        begin = time.perf_counter()
        result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "-O", str(output)],
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in gzip.open(output, "rt", encoding="utf-8")
                      if line and not line.startswith("#"))
        print(json.dumps({
            "status": "pass",
            "tool": "ApplyVQSR",
            "input_records": args.records,
            "output_records": records,
            "wall_seconds": elapsed,
            "records_per_second": args.records / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
