#!/usr/bin/env python3
"""Replay GATK's exact PairHMM debug rows through the Kokkos kernel.

This is deliberately narrower than the Mutect2 VCF oracle.  It asks GATK
4.6.2.0 to emit the read/haplotype/quality rows it actually sent to PairHMM,
then replays those rows in the native bucketed kernel.  A pass therefore
isolates arithmetic compatibility; a remaining Mutect2 TLOD delta is an
assembly/request-set difference and cannot be hidden by changing a tolerance
on the VCF oracle.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_PAIRHMM_RESULTS_ORACLE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-pairhmm-results-oracle"),
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK PairHMM oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK PairHMM oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-pairhmm-results-oracle-") as directory:
        work = Path(directory)
        pair_results = work / "gatk.pairhmm.txt"
        gatk_output = work / "gatk.vcf.gz"
        gatk_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(gatk_output), "--pair-hmm-results-file", str(pair_results),
        ])
        if gatk_run.returncode != 0:
            raise AssertionError(f"GATK PairHMM debug run failed: {gatk_run.stderr[-4000:]}")
        if not pair_results.is_file() or pair_results.stat().st_size == 0:
            raise AssertionError("GATK produced no PairHMM debug rows")
        replay = run([str(binary), str(pair_results)])
        if replay.returncode != 0:
            raise AssertionError(f"Kokkos PairHMM replay failed: {replay.stderr[-4000:]}")
        report = None
        for line in reversed(replay.stdout.splitlines()):
            if line.startswith("{") and line.endswith("}"):
                report = json.loads(line)
                break
        if report is None or report.get("status") != "pass":
            raise AssertionError(f"invalid PairHMM replay report: {replay.stdout[-4000:]}")
        report.update({
            "gatk_version": "4.6.2.0",
            "pairhmm_results_rows": sum(
                1 for line in pair_results.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")
            ),
            "assembly_boundary": "debug-request-set-isolated",
        })
        print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
