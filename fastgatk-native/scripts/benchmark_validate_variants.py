#!/usr/bin/env python3
"""File-boundary throughput benchmark for ValidateVariants.

The benchmark keeps reference/dbSNP I/O out of the timing path while still
exercising HTSlib decode, allele checks, chromosome counts and FORMAT
cardinality validation on a deterministic synthetic VCF.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=50000)
    args = parser.parse_args()
    if args.records < 1:
        raise SystemExit("--records must be positive")
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_VALIDATE_VARIANTS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-validate-variants"),
    ))
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"missing executable: {binary}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-validate-variants-benchmark-") as directory:
        work = pathlib.Path(directory)
        variant = work / "input.vcf"
        report = work / "validation.tsv"
        manifest = work / "validation.manifest.json"
        with variant.open("w", encoding="ascii") as stream:
            stream.write(
                "##fileformat=VCFv4.2\n"
                f"##contig=<ID=chr1,length={args.records + 1}>\n"
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            )
            for position in range(1, args.records + 1):
                ref = "A" if position % 2 else "C"
                alt = "G" if position % 2 else "T"
                stream.write(
                    f"chr1\t{position}\trs{position}\t{ref}\t{alt}\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n"
                )
        begin = time.perf_counter()
        completed = subprocess.run(
            [str(binary), "-V", str(variant), "-O", str(report),
             "--output-manifest", str(manifest),
             "--validation-type-to-exclude", "REF",
             "--validation-type-to-exclude", "IDS"],
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        summary = json.loads(completed.stdout.splitlines()[-1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["validated_records"] == args.records
        print(json.dumps({
            "status": "pass",
            "tool": "ValidateVariants",
            "records": args.records,
            "wall_seconds": elapsed,
            "records_per_second": args.records / elapsed if elapsed else 0.0,
            "output_bytes": report.stat().st_size,
            "manifest_bytes": manifest.stat().st_size,
            "summary": summary,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
