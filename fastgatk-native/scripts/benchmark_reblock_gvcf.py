#!/usr/bin/env python3
"""File-boundary benchmark for the Kokkos ReblockGVCF PL-remap path."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import statistics
import subprocess
import tempfile
import time


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
"""


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY",
        root / "fastgatk-native/build/fastgatk-reblock-gvcf",
    ))
    if not binary.is_file():
        raise SystemExit("native ReblockGVCF binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-benchmark-") as directory:
        work = pathlib.Path(directory)
        source = work / "input.g.vcf.gz"
        lines = [HEADER]
        for index in range(512):
            position = index * 100 + 1
            if index % 16 == 0:
                lines.append(
                    f"chr1\t{position}\t.\tA\tC,G,<NON_REF>\t.\tPASS\tDP=20\t"
                    "GT:DP:AD:PL:GQ\t"
                    "0/1:20:12,8,0,0:0,10,20,30,40,50,99,99,99,99:99\t"
                    "0/1:20:12,8,0,0:0,10,20,30,40,50,99,99,99,99:99\t"
                    "0/1:20:12,8,0,0:0,10,20,30,40,50,99,99,99,99:99\t"
                    "0/1:20:12,8,0,0:0,10,20,30,40,50,99,99,99,99:99\n"
                )
            else:
                lines.append(
                    f"chr1\t{position}\t.\tA\t<NON_REF>\t.\tPASS\tEND={position + 99}\t"
                    "GT:DP:AD:PL:GQ\t"
                    "0/0:10:10,0:0,20,99:20\t0/0:10:10,0:0,20,99:20\t"
                    "0/0:10:10,0:0,20,99:20\t0/0:10:10,0:0,20,99:20\n"
                )
        with gzip.open(source, "wt", encoding="ascii") as stream:
            stream.write("".join(lines))
        samples: list[float] = []
        manifest = work / "output.manifest.json"
        output = work / "output.g.vcf.gz"
        command = [str(binary), "-V", str(source), "-O", str(output),
                   "--output-manifest", str(manifest), "--gvcf-gq-bands", "20",
                   "--gvcf-gq-bands", "100"]
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
            "suite": "fastgatk-reblock-gvcf-file-boundary",
            "tool": "ReblockGVCF",
            "backend": "Kokkos",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": source.stat().st_size,
            "output_bytes": output.stat().st_size,
            "sample_count": metadata["telemetry"]["sample_count"],
            "execution_space": metadata["telemetry"].get(
                "pl_remap_kernel_execution_space",
                metadata["telemetry"].get("non_ref_ad_kernel_execution_space", ""),
            ),
            "pl_remap_kernel_calls": metadata["telemetry"]["pl_remap_kernel_calls"],
            "pl_remap_kernel_execution_space": metadata["telemetry"]["pl_remap_kernel_execution_space"],
            "pl_remap_kernel_seconds": metadata["telemetry"]["pl_remap_kernel_seconds"],
            "allele_field_remap_kernel_calls": metadata["telemetry"]["allele_field_remap_kernel_calls"],
            "allele_field_remap_kernel_seconds": metadata["telemetry"]["allele_field_remap_kernel_seconds"],
            "non_ref_ad_kernel_lifecycle": metadata["telemetry"]["non_ref_ad_kernel_lifecycle"],
            "non_ref_ad_kernel_execution_policy": metadata["telemetry"]["non_ref_ad_kernel_execution_policy"],
            "non_ref_ad_kernel_calls": metadata["telemetry"]["non_ref_ad_kernel_calls"],
            "non_ref_ad_kernel_prepare_seconds": metadata["telemetry"]["non_ref_ad_kernel_prepare_seconds"],
            "non_ref_ad_kernel_execute_seconds": metadata["telemetry"]["non_ref_ad_kernel_execute_seconds"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
