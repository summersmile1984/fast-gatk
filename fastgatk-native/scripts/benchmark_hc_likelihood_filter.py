#!/usr/bin/env python3
"""Benchmark the Host-side HC read-disqualification contract."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    cases = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-likelihood-benchmark-") as directory:
        work = Path(directory)
        for label, rate in (("gatk-default", "0.02"), ("permissive", "0.5")):
            timings = []
            manifest_data = {}
            for repeat in range(3):
                output = work / f"{label}-{repeat}.vcf"
                manifest = work / f"{label}-{repeat}.manifest.json"
                start = time.perf_counter()
                subprocess.run([
                    str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
                    "-O", str(output), "--threads", "2", "--min-depth", "1",
                    "--min-alt-support", "1", "--expected-mismatch-rate-for-read-disqualification",
                    rate, "--output-manifest", str(manifest),
                ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                timings.append(time.perf_counter() - start)
                manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            timings.sort()
            telemetry = manifest_data["telemetry"]
            cases.append({"label": label, "expected_error_rate_per_base": float(rate),
                          "warmup_seconds": timings[0], "p50_seconds": timings[1],
                          "p95_seconds": timings[-1], "output_bytes": output.stat().st_size,
                          "pairhmm_reads_disqualified": telemetry["pairhmm_reads_disqualified"],
                          "pairhmm_reads_clipped": telemetry["pairhmm_reads_clipped"],
                          "pairhmm_execution_space": telemetry["pairhmm_execution_space"]})
    assert cases[0]["pairhmm_reads_disqualified"] >= 1
    assert cases[1]["pairhmm_reads_disqualified"] == 0
    print(json.dumps({"schema_version": 1, "status": "pass", "tool": "HaplotypeCaller",
                      "benchmark": "pairhmm-read-disqualification-file-boundary",
                      "cases": cases}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
