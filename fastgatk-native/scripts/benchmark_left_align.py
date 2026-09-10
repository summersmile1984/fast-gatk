#!/usr/bin/env python3
"""File-boundary benchmark for the Kokkos-backed LeftAlignAndTrim path."""
from __future__ import annotations

import gzip
import json
import os
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=20000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    index = min(len(ordered) - 1, int(round((len(ordered) - 1) * fraction)))
    return ordered[index]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_LEFT_ALIGN_BINARY",
        str(root / "fastgatk-native/build/fastgatk-left-align-trim"),
    ))
    assert binary.is_file() and os.access(binary, os.X_OK)
    repetitions = max(2, int(os.environ.get("FASTGATK_LEFT_ALIGN_BENCH_REPEATS", "5")))
    with tempfile.TemporaryDirectory(prefix="fastgatk-left-align-bench-") as directory:
        work = Path(directory)
        reference = work / "ref.fa"
        source = work / "input.vcf"
        reference.write_text(">chr1\n" + "A" * 20000 + "\n", encoding="utf-8")
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20000\n", encoding="utf-8")
        reference.with_suffix(".fa.fai").write_text("chr1\t20000\t6\t20000\t20001\n", encoding="utf-8")
        with source.open("w", encoding="utf-8") as handle:
            handle.write(HEADER)
            for index in range(1, 1001):
                position = 10 + index * 10
                # Most rows are biallelic SNPs; every tenth record is a
                # triploid multiallelic split that exercises the remap kernel.
                if index % 10 == 0:
                    handle.write(
                        f"chr1\t{position}\tid{index}\tA\tC,G\t50\tPASS\tDP=20;AC=1,1;AN=3;AF=0.333333,0.333333\t"
                        "GT:AD:PL:GQ\t0/1/2:10,3,4:90,40,70,30,10,80,20,60,50,0:0\n")
                else:
                    handle.write(
                        f"chr1\t{position}\tid{index}\tA\tC\t50\tPASS\tDP=10;AC=1;AN=2;AF=0.5\t"
                        "GT:AD:PL:GQ\t0/1:8,2:20,0,30:10\n")

        def run_once() -> tuple[float, dict[str, object]]:
            output = work / f"out-{time.time_ns()}.vcf.gz"
            started = time.perf_counter()
            result = subprocess.run([
                str(binary), "-V", str(source), "-R", str(reference), "-O", str(output),
                "--split-multi-allelics",
            ], text=True, capture_output=True, check=False)
            elapsed = time.perf_counter() - started
            assert result.returncode == 0, result.stderr
            summary = json.loads(result.stdout.splitlines()[-1])
            manifest = json.loads(Path(f"{output}.manifest.json").read_text(encoding="utf-8"))
            assert output.is_file() and Path(f"{output}.tbi").is_file()
            assert summary["pl_remap_kernel_calls"] == 200
            assert manifest["telemetry"]["max_ploidy"] == 3
            return elapsed, {"summary": summary, "manifest": manifest, "output_bytes": output.stat().st_size}

        warmup_seconds, _ = run_once()
        samples: list[float] = []
        last: dict[str, object] = {}
        for _ in range(repetitions):
            elapsed, last = run_once()
            samples.append(elapsed)
        manifest = last["manifest"]
        telemetry = manifest["telemetry"]
        report = {
            "schema_version": 1,
            "status": "pass",
            "suite": "fastgatk-left-align-file-boundary",
            "tool": "LeftAlignAndTrimVariants",
            "records": 1000,
            "split_records": 1100,
            "repetitions": repetitions,
            "warmup_seconds": warmup_seconds,
            "p50_seconds": statistics.median(samples),
            "p95_seconds": percentile(samples, 0.95),
            "output_bytes": last["output_bytes"],
            "pl_remap_kernel_calls": telemetry["pl_remap_kernel_calls"],
            "allele_field_remap_kernel_calls": telemetry["allele_field_remap_kernel_calls"],
            "execution_space": telemetry["pl_remap_kernel_execution_space"],
            "max_ploidy": telemetry["max_ploidy"],
        }
        print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
