#!/usr/bin/env python3
"""Pin paired-fragment key construction without MC against Picard.

Picard's MarkDuplicates pairs mapped mates by read name and uses the actual
unclipped five-prime coordinate of each alignment.  A native implementation
that derives the second end only from SAM mate fields splits ordinary pairs
when MC is absent (and can miss terminal clipping).  This fixture deliberately
omits MC and includes terminal soft clips, then compares flags, tags and the
paired metrics row with GATK/Picard 4.6.2.0.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

from verify_mark_duplicates_gatk_oracle import metric_semantics, semantic_records
import oracle_guard


SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
@RG\tID:rg\tSM:S\tLB:lib
A:1:1:100:100\t99\tchr1\t11\t60\t2S4M\t=\t31\t24\tAACGTA\tIIIIII\tRG:Z:rg
B:1:1:101:101\t99\tchr1\t11\t20\t2S4M\t=\t31\t24\tAACGTA\t!!!!!!\tRG:Z:rg
A:1:1:100:100\t147\tchr1\t31\t60\t4M2S\t=\t11\t-24\tTGCAAA\t!!!!!!\tRG:Z:rg
B:1:1:101:101\t147\tchr1\t31\t20\t4M2S\t=\t11\t-24\tTGCAAA\tIIIIII\tRG:Z:rg
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_MARK_DUPLICATES_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-mark-duplicates"
    if not (java.is_file() and gatk.is_file() and native.is_file() and os.access(native, os.X_OK)):
        oracle_guard.oracle_not_verified('verify_mark_duplicates_pair_key_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned GATK/Picard MarkDuplicates oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned GATK/Picard oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mark-duplicates-pair-key-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        source.write_text(SAM, encoding="utf-8")
        java_output = work / "picard.sam"
        native_output = work / "native.sam"
        java_metrics = work / "picard.metrics"
        native_metrics = work / "native.metrics"
        run([str(java), "-Xmx1g", "-jar", str(gatk), "MarkDuplicates",
             "-I", str(source), "-O", str(java_output), "-M", str(java_metrics),
             "--ASSUME_SORTED", "true", "--CREATE_INDEX", "false",
             "--TAGGING_POLICY", "All"])
        native_result = run([str(native), "-I", str(source), "-O", str(native_output),
                             "--metrics-file", str(native_metrics),
                             "--create-output-bam-index=false", "--tagging-policy", "All"])
        if semantic_records(java_output) != semantic_records(native_output):
            raise AssertionError({"picard": semantic_records(java_output),
                                  "native": semantic_records(native_output)})
        if metric_semantics(java_metrics) != metric_semantics(native_metrics):
            raise AssertionError({"picard_metrics": metric_semantics(java_metrics),
                                  "native_metrics": metric_semantics(native_metrics)})
        summary = json.loads(native_result.stdout.strip().splitlines()[-1])
        if summary.get("duplicate_records") != 2:
            raise AssertionError({"summary": summary})

    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "records_compared": 4, "duplicate_records": 2,
                      "no_mc_pair_key_exact": True, "unclipped_coordinate_exact": True},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
