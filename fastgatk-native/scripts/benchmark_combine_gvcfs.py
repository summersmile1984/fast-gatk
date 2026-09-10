#!/usr/bin/env python3
"""File-boundary benchmark for CombineGVCFs Kokkos FORMAT projection."""

from __future__ import annotations

import gzip
import json
import pathlib
import statistics
import subprocess
import tempfile
import time
import argparse


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stream-merge", action="store_true",
                        help="benchmark bounded one-record-per-input k-way merge")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-combine-gvcfs"
    if not binary.is_file():
        raise SystemExit("native CombineGVCFs binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-combine-gvcfs-benchmark-") as directory:
        work = pathlib.Path(directory)
        shard_g = work / "shard-g.g.vcf.gz"
        shard_t = work / "shard-t.g.vcf.gz"
        with gzip.open(shard_g, "wt", encoding="ascii") as first, gzip.open(
            shard_t, "wt", encoding="ascii"
        ) as second:
            first.write(HEADER.format(sample="S1"))
            second.write(HEADER.format(sample="S2"))
            for index in range(128):
                position = index * 10 + 1
                first.write(
                    f"chr1\t{position}\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
                    "GT:DP:AD:PL\t0/1:20:12,8,0:50,0,60,70,80,90\n"
                )
                second.write(
                    f"chr1\t{position}\t.\tA\tT,G,<NON_REF>\t.\tPASS\tDP=21\t"
                    "GT:DP:AD:PL\t1/2:21:2,9,10,0:90,80,70,60,50,0\n"
                )
        output = work / "output.g.vcf.gz"
        manifest = work / "output.manifest.json"
        command = [
            str(binary), "-V", str(shard_g), "-V", str(shard_t),
            "-O", str(output), "--output-manifest", str(manifest),
            "--call-genotypes",
        ]
        if args.stream_merge:
            command.append("--stream-merge")
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
        telemetry = metadata["telemetry"]
        if not args.stream_merge and (telemetry.get("pl_remap_kernel_calls", 0) < 128 or
                                      telemetry.get("genotype_kernel_calls", 0) < 128):
            raise SystemExit("CombineGVCFs Kokkos projection kernels were not exercised")
        ordered = sorted(samples)
        p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]
        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-combine-gvcfs-file-boundary",
            "tool": "CombineGVCFs",
            "backend": "Kokkos",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "samples": 2,
            "records_per_shard": 128,
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": shard_g.stat().st_size + shard_t.stat().st_size,
            "output_bytes": output.stat().st_size,
            "pl_remap_kernel_calls": telemetry.get("pl_remap_kernel_calls", 0),
            "pl_remap_kernel_seconds": telemetry.get("pl_remap_kernel_seconds", 0.0),
            "allele_field_remap_kernel_calls": telemetry.get("allele_field_remap_kernel_calls", 0),
            "allele_field_remap_kernel_seconds": telemetry.get("allele_field_remap_kernel_seconds", 0.0),
            "genotype_kernel_calls": telemetry.get("genotype_kernel_calls", 0),
            "genotype_kernel_seconds": telemetry.get("genotype_kernel_seconds", 0.0),
            "genotype_kernel_execution_space": telemetry.get("genotype_kernel_execution_space", ""),
            "call_genotypes": telemetry.get("call_genotypes", True),
            "stream_merge": args.stream_merge,
            "max_inflight_records": telemetry.get("max_inflight_records", 0),
            "peak_queued_bytes": telemetry.get("peak_queued_bytes", 0),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
