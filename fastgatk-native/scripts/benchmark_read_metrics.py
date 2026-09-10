#!/usr/bin/env python3
"""File-boundary benchmark for native CountReads and FlagStat."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
BINARIES = {
    "CountReads": Path(os.environ.get("FASTGATK_COUNT_READS_BINARY", ROOT / "fastgatk-native/build/fastgatk-count-reads")),
    "FlagStat": Path(os.environ.get("FASTGATK_FLAG_STAT_BINARY", ROOT / "fastgatk-native/build/fastgatk-flag-stat")),
}
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK_JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--input-shards", type=int, default=1,
                        help="repeat the fixture this many times as -I shards")
    parser.add_argument("--include-java", action="store_true",
                        help="also time the pinned GATK Java implementation")
    args = parser.parse_args()
    if args.input_shards < 1:
        parser.error("--input-shards must be positive")
    if not BAM.is_file() or any(not path.is_file() for path in BINARIES.values()):
        print(json.dumps({"status": "skip", "suite": "read-metrics-benchmark"}))
        return 0
    results = []
    input_count = args.input_shards
    input_records = 493 * input_count
    shard_arguments = []
    for _ in range(input_count):
        shard_arguments.extend(["-I", str(BAM)])
    with tempfile.TemporaryDirectory(prefix="fastgatk-read-metrics-bench-") as directory:
        work = Path(directory)
        for tool, binary in BINARIES.items():
            timings: list[float] = []
            for _ in range(max(1, args.iterations)):
                output = work / f"{tool}.txt"
                manifest = work / f"{tool}.manifest.json"
                start = time.perf_counter()
                result = subprocess.run(
                    [str(binary), *shard_arguments, "-O", str(output), "--threads", str(args.threads),
                     "--output-manifest", str(manifest)],
                    text=True, capture_output=True,
                )
                if result.returncode != 0:
                    raise RuntimeError(result.stderr)
                timings.append(time.perf_counter() - start)
            timings.sort()
            p50 = timings[len(timings) // 2]
            native_manifest = json.loads(manifest.read_text(encoding="utf-8"))
            telemetry = native_manifest["telemetry"]
            results.append({"tool": tool, "threads": args.threads, "iterations": len(timings),
                            "input_count": input_count, "input_records": input_records,
                            "p50_seconds": p50, "p95_seconds": timings[min(len(timings) - 1, int(len(timings) * 0.95))],
                            "records_per_second": input_records / p50 if p50 else None,
                            "backend": "native-kokkos",
                            "pipeline_lifecycle": telemetry["pipeline_lifecycle"],
                            "pipeline_decoded_items": telemetry["pipeline_decoded_items"],
                            "pipeline_computed_items": telemetry["pipeline_computed_items"],
                            "pipeline_encoded_items": telemetry["pipeline_encoded_items"],
                            "pipeline_peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
                            "pipeline_peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
                            "pipeline_peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
                            "effective_threads": native_manifest["effective_threads"]})
        if args.include_java and JAVA.is_file() and GATK_JAR.is_file():
            for tool in BINARIES:
                timings = []
                for _ in range(max(1, args.iterations)):
                    start = time.perf_counter()
                    result = subprocess.run(
                        [str(JAVA), "-jar", str(GATK_JAR), tool, *shard_arguments],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                    )
                    if result.returncode != 0:
                        raise RuntimeError(f"GATK {tool} failed with exit code {result.returncode}")
                    timings.append(time.perf_counter() - start)
                timings.sort()
                p50 = timings[len(timings) // 2]
                results.append({"tool": tool, "iterations": len(timings),
                                "input_count": input_count, "input_records": input_records,
                                "p50_seconds": p50,
                                "p95_seconds": timings[min(len(timings) - 1, int(len(timings) * 0.95))],
                                "records_per_second": input_records / p50 if p50 else None,
                                "backend": "gatk-java-4.6.2.0"})
    print(json.dumps({"status": "pass", "suite": "read-metrics-benchmark", "input_records": input_records,
                      "input_count": input_count,
                      "input_bytes": BAM.stat().st_size * input_count, "results": results}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
