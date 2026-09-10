#!/usr/bin/env python3
"""File-boundary benchmark for the native GetPileupSummaries pipeline."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GET_PILEUP_BINARY",
        str(root / "fastgatk-native/build/fastgatk-get-pileup-summaries"),
    ))
    records = int(os.environ.get("FASTGATK_PILEUP_BENCH_RECORDS", "20000"))
    sites = int(os.environ.get("FASTGATK_PILEUP_BENCH_SITES", "100"))
    iterations = int(os.environ.get("FASTGATK_PILEUP_BENCH_ITERATIONS", "3"))
    if min(records, sites, iterations) < 1:
        raise SystemExit("benchmark sizes must be positive")

    with tempfile.TemporaryDirectory(prefix="fastgatk-pileup-benchmark-") as directory:
        work = Path(directory)
        sam = work / "reads.sam"
        variants = work / "sites.vcf"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000000\n"
            "@RG\tID:rg1\tSM:BENCH\n",
            encoding="utf-8",
        )
        with sam.open("a", encoding="utf-8") as stream:
            for index in range(records):
                start = (index % (sites * 50)) + 1
                stream.write(
                    f"r{index}\t0\tchr1\t{start}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\n"
                )
        with variants.open("w", encoding="utf-8") as stream:
            stream.write(
                "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n"
                "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            )
            stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for index in range(sites):
                stream.write(f"chr1\t{2 + index * 50}\t.\tC\tA\t.\tPASS\tAF=0.1\n")

        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        timings: list[float] = []
        manifest_payload: dict[str, object] | None = None
        for index in range(iterations):
            output = work / f"pileups-{index}.table"
            manifest = work / f"pileups-{index}.manifest.json"
            start = time.perf_counter()
            completed = subprocess.run([
                str(native), "-I", str(sam), "-V", str(variants), "-L", str(variants),
                "-O", str(output), "--batch-records", "1024", "--threads", "4",
                "--output-manifest", str(manifest),
            ], check=True, text=True, capture_output=True, env=env)
            timings.append(time.perf_counter() - start)
            manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
            telemetry = manifest_payload["telemetry"]
            if telemetry["pipeline_lifecycle"] != "decode->compute->encode->sink":
                raise AssertionError(manifest_payload)
            if telemetry["pipeline_decoded_items"] != telemetry["pipeline_encoded_items"]:
                raise AssertionError(manifest_payload)
            if output.stat().st_size == 0 or not completed.stdout.strip():
                raise AssertionError("GetPileupSummaries produced an empty output")

        timings.sort()
        assert manifest_payload is not None
        telemetry = manifest_payload["telemetry"]
        p50 = timings[len(timings) // 2]
        print(json.dumps({
            "status": "pass",
            "tool": "GetPileupSummaries",
            "records": records,
            "sites": sites,
            "iterations": iterations,
            "p50_seconds": p50,
            "p95_seconds": timings[min(len(timings) - 1, int(len(timings) * 0.95))],
            "records_per_second": records / p50 if p50 else None,
            "kernel": {
                "execution_space": telemetry["count_kernel_execution_space"],
                "batches": telemetry["count_kernel_batches"],
                "observations": telemetry["count_kernel_observations"],
                "prepare_seconds": telemetry["count_kernel_prepare_seconds"],
                "execute_seconds": telemetry["count_kernel_seconds"],
            },
            "pipeline": {
                "lifecycle": telemetry["pipeline_lifecycle"],
                "decoded_items": telemetry["pipeline_decoded_items"],
                "computed_items": telemetry["pipeline_computed_items"],
                "encoded_items": telemetry["pipeline_encoded_items"],
                "peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
                "peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
                "peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
