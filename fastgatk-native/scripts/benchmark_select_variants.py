#!/usr/bin/env python3
"""File-boundary benchmark for SelectVariants sample/ALT PL compaction."""

from __future__ import annotations

import gzip
import json
import pathlib
import statistics
import subprocess
import tempfile
import time


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AS_INT,Number=A,Type=Integer,Description=Per-ALT score>
##INFO=<ID=AS_FLOAT,Number=A,Type=Float,Description=Per-ALT frequency>
##INFO=<ID=ADINFO,Number=R,Type=Integer,Description=Per-allele score>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-select-variants"
    if not binary.is_file():
        raise SystemExit("native SelectVariants binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-select-benchmark-") as directory:
        work = pathlib.Path(directory)
        source = work / "input.vcf.gz"
        lines = [HEADER]
        for index in range(512):
            position = index * 100 + 1
            lines.append(
                f"chr1\t{position}\trs{index}\tA\tC,G\t50\tPASS\tDP=20;AS_INT=11,22;AS_FLOAT=0.1,0.2;ADINFO=100,11,22\t"
                "GT:AD:PL:GQ\t0/0:20,0,0:0,50,60,70,80,90:50\t"
                "0/1:12,8,0:50,0,60,70,80,90:50\n"
            )
        with gzip.open(source, "wt", encoding="ascii") as stream:
            stream.write("".join(lines))
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        command = [str(binary), "-V", str(source), "-O", str(output),
                   "--sample-name", "S2", "--remove-unused-alternates",
                   "--output-manifest", str(manifest)]
        samples: list[float] = []
        for iteration in range(7):
            started = time.perf_counter()
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            elapsed = time.perf_counter() - started
            if completed.returncode != 0:
                raise SystemExit(completed.stderr or completed.stdout)
            if iteration > 0:
                samples.append(elapsed)
        ordered = sorted(samples)
        p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-select-variants-file-boundary",
            "tool": "SelectVariants",
            "backend": "Kokkos",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": source.stat().st_size,
            "output_bytes": output.stat().st_size,
            "pl_remap_kernel_calls": metadata["telemetry"]["pl_remap_kernel_calls"],
            "gt_gq_kernel_calls": metadata["telemetry"]["gt_gq_kernel_calls"],
            "genotype_kernel_execution_space": metadata["telemetry"]["genotype_kernel_execution_space"],
            "pl_remap_kernel_seconds": metadata["telemetry"]["pl_remap_kernel_seconds"],
            "allele_field_remap_kernel_calls": metadata["telemetry"]["allele_field_remap_kernel_calls"],
            "allele_field_remap_kernel_seconds": metadata["telemetry"]["allele_field_remap_kernel_seconds"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
