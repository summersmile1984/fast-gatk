#!/usr/bin/env python3
"""Pin Picard's OpticalOnly duplicate-tagging policy against the native path.

The broad MarkDuplicates oracle exercises ``DontTag``/``All`` and removal.
Picard also exposes ``OpticalOnly``: duplicate records caused by an optical
cluster receive ``DT:Z:SQ`` while ordinary library duplicates receive no DT
tag.  This is an output-semantic boundary, not merely a parser alias, and the
fixture contains both kinds of paired duplicate plus unpaired duplicates.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

from verify_mark_duplicates_gatk_oracle import (
    SAM,
    histogram_semantics,
    metric_semantics,
    pg_header,
    run,
    semantic_records,
)
import oracle_guard


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get("FASTGATK_MARK_DUPLICATES_BINARY",
                                os.environ.get("FASTGATK_NATIVE_BUILD",
                                               root / "fastgatk-native/build")))
    if native.is_dir():
        native /= "fastgatk-mark-duplicates"
    if not (java.is_file() and gatk.is_file() and native.is_file() and os.access(native, os.X_OK)):
        oracle_guard.oracle_not_verified('verify_mark_duplicates_tagging_policy_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned GATK/Picard MarkDuplicates oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned GATK/Picard oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mark-duplicates-optical-only-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        source.write_text(SAM, encoding="utf-8")
        java_output = work / "picard.sam"
        native_output = work / "native.sam"
        java_metrics = work / "picard.metrics"
        native_metrics = work / "native.metrics"
        java_result = subprocess.run(
            [str(java), "-Xmx1g", "-jar", str(gatk), "MarkDuplicates",
             "-I", str(source), "-O", str(java_output), "-M", str(java_metrics),
             "--ASSUME_SORTED", "true", "--CREATE_INDEX", "false",
             "--TAGGING_POLICY", "OpticalOnly"],
            text=True, capture_output=True, check=False)
        if java_result.returncode != 0:
            raise AssertionError(java_result.stderr)
        native_result = run(
            [str(native), "-I", str(source), "-O", str(native_output),
             "--metrics-file", str(native_metrics),
             "--create-output-bam-index=false", "--tagging-policy", "OpticalOnly"])
        if semantic_records(java_output) != semantic_records(native_output):
            raise AssertionError({"records": (semantic_records(java_output), semantic_records(native_output))})
        if metric_semantics(java_metrics) != metric_semantics(native_metrics):
            raise AssertionError({"metrics": (metric_semantics(java_metrics), metric_semantics(native_metrics))})
        if histogram_semantics(java_metrics) != histogram_semantics(native_metrics):
            raise AssertionError({"histogram": (histogram_semantics(java_metrics), histogram_semantics(native_metrics))})
        if not any(header.get("ID") == "MarkDuplicates" for header in pg_header(native_output)):
            raise AssertionError("native output is missing the MarkDuplicates @PG header")

        rows = [line.split("\t") for line in native_output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("@")]
        optical_tags = [row for row in rows if "DT:Z:SQ" in row[11:]]
        ordinary_tags = [row for row in rows if "DT:Z:LB" in row[11:]]
        if len(optical_tags) != 2 or ordinary_tags:
            raise AssertionError({"optical_sq_records": len(optical_tags),
                                  "ordinary_lb_records": len(ordinary_tags)})

    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "records_compared": 14, "optical_dt_records": 2,
                      "ordinary_dt_records": 0,
                      "metrics_and_histogram_semantically_identical": True,
                      "tagging_policy": "OpticalOnly"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
