#!/usr/bin/env python3
"""File-boundary benchmark for Kokkos VariantRecalibrator scoring."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_VARIANT_RECALIBRATOR_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-variant-recalibrator")))
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=10_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-gaussians", type=int, default=1)
    parser.add_argument("--max-negative-gaussians", type=int, default=2)
    parser.add_argument("--max-iterations", type=int, default=20)
    parser.add_argument("--k-means-iterations", type=int, default=5)
    parser.add_argument("--maximum-training-variants", type=int, default=2_500_000)
    parser.add_argument("--shrinkage", type=float, default=1.0)
    parser.add_argument("--dirichlet", type=float, default=0.001)
    parser.add_argument("--prior-counts", type=float, default=20.0)
    parser.add_argument("--full-covariance", action="store_true")
    args = parser.parse_args()
    if (args.records < 1 or args.threads < 1 or args.max_gaussians < 1 or
            args.max_negative_gaussians < 1 or args.max_iterations < 1 or
            args.k_means_iterations < 1 or args.maximum_training_variants < 1):
        raise SystemExit("records/threads/gaussian/iteration options must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-bench-") as temporary:
        work = pathlib.Path(temporary)
        input_vcf = work / "input.vcf"
        resource_vcf = work / "resource.vcf"
        output = work / "recal.vcf.gz"
        tranches = work / "recal.tranches"
        manifest = pathlib.Path(str(output) + ".manifest.json")
        with input_vcf.open("w", encoding="utf-8") as data, resource_vcf.open("w", encoding="utf-8") as resource:
            data.write(HEADER)
            resource.write(HEADER)
            for index in range(args.records):
                position = index * 10 + 1
                qd = 20.0 + (index % 100) / 10.0
                mq = 45.0 + (index % 50) / 10.0
                row = f"chr1\t{position}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\n"
                data.write(row)
                if index % 2 == 0:
                    resource.write(row)
        begin = time.perf_counter()
        command = [str(BINARY), "-V", str(input_vcf),
                   "--resource:truth,training=true,truth=true", str(resource_vcf),
                   "-an", "QD", "-an", "MQ", "-O", str(output), "--tranches-file", str(tranches),
                   "--threads", str(args.threads), "--max-gaussians", str(args.max_gaussians)]
        command += ["--max-negative-gaussians", str(args.max_negative_gaussians),
                    "--max-iterations", str(args.max_iterations),
                    "--k-means-iterations", str(args.k_means_iterations),
                    "--maximum-training-variants", str(args.maximum_training_variants),
                    "--shrinkage", str(args.shrinkage), "--dirichlet", str(args.dirichlet),
                    "--prior-counts", str(args.prior_counts)]
        if args.full_covariance:
            command.append("--full-covariance")
        result = subprocess.run(
            command,
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in gzip.open(output, "rt", encoding="utf-8")
                      if line and not line.startswith("#"))
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        print(json.dumps({
            "status": "pass",
            "tool": "VariantRecalibrator",
            "input_records": args.records,
            "output_records": records,
            "wall_seconds": elapsed,
            "records_per_second": args.records / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "kernel": {
                "lifecycle": manifest_payload["kernel_lifecycle"],
                "execution_policy": manifest_payload["kernel_execution_policy"],
                "batches": manifest_payload["kernel_batches"],
                "observations": manifest_payload["kernel_observations"],
                "prepare_seconds": manifest_payload["kernel_prepare_seconds"],
                "execute_seconds": manifest_payload["kernel_execute_seconds"],
            },
            "summary": json.loads(result.stdout),
            "training_rows": manifest_payload["training_records"],
            "background_rows": manifest_payload["background_records"],
            "bad_cutoff_selected": manifest_payload["bad_cutoff_selected"],
            "bad_fallback_selected": manifest_payload["bad_fallback_selected"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
