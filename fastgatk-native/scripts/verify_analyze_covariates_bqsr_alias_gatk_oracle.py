#!/usr/bin/env python3
"""Pinned GATK oracle for AnalyzeCovariates' ``-bqsr`` selector alias.

GATK 4.6.2.0 spells the short alias in lower case (``-bqsr``), while the
native prototype historically accepted only its legacy upper-case ``-BQSR``
spelling.  This gate checks the parser boundary without conflating it with
the separate, already-bounded numerical AnalyzeCovariates report oracle:
Java's lower/long spellings must produce the same CSV and native's lower/
upper spellings must produce the same CSV.  Barclay rejects an embedded '='
for this short alias, so that failure is pinned as well.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def invoke(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed ({result.returncode}): {result.stderr[-3000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_ANALYZE_COVARIATES_BINARY",
        root / "fastgatk-native/build/fastgatk-analyze-covariates"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    # This is a checked-in GATK v1.1 report from the pinned BQSR test corpus;
    # using it keeps this CLI-only gate fast and independent of recalibration.
    report = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/BQSR/"
        "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt")
    required = (native, java, jar, report)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_analyze_covariates_bqsr_alias_gatk_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("AnalyzeCovariates -bqsr oracle inputs are required")
        print('{"status":"skip","reason":"GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-analyze-covariates-bqsr-alias-") as directory:
        work = Path(directory)
        java_lower = work / "java.lower.csv"
        java_long = work / "java.long.csv"
        native_lower = work / "native.lower.csv"
        native_upper = work / "native.upper.csv"
        manifest = work / "native.manifest.json"

        lower = invoke([
            str(java), "-jar", str(jar), "AnalyzeCovariates", "-bqsr", str(report),
            "-csv", str(java_lower)])
        require_success(lower, "GATK -bqsr")
        long = invoke([
            str(java), "-jar", str(jar), "AnalyzeCovariates",
            "--bqsr-recal-file", str(report), "-csv", str(java_long)])
        require_success(long, "GATK --bqsr-recal-file")
        if java_lower.read_bytes() != java_long.read_bytes():
            raise AssertionError("GATK -bqsr and --bqsr-recal-file CSV differ")

        native_result = invoke([
            str(native), "-bqsr", str(report), "-csv", str(native_lower),
            "--output-manifest", str(manifest)])
        require_success(native_result, "native -bqsr")
        native_legacy = invoke([
            str(native), "-BQSR", str(report), "-csv", str(native_upper)])
        require_success(native_legacy, "native -BQSR")
        if native_lower.read_bytes() != native_upper.read_bytes():
            raise AssertionError("native -bqsr and -BQSR CSV differ")

        # GATK's short-option parser does not accept -bqsr=...; keeping this
        # explicit prevents a launcher rewrite from silently changing scope.
        embedded = invoke([
            str(native), f"-bqsr={report}", "-csv", str(work / "embedded.csv")])
        if embedded.returncode == 0:
            raise AssertionError("native accepted GATK-invalid embedded -bqsr= form")

        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        if metadata["tool"] != "AnalyzeCovariates":
            raise AssertionError("manifest tool mismatch")
        if metadata["telemetry"]["input_reports"] != 1:
            raise AssertionError("manifest input report count mismatch")
        if not all(item["complete"] for item in metadata["outputs"]):
            raise AssertionError("native alias output is incomplete")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "gatk_lower_alias": True,
            "gatk_lower_long_csv_identical": True,
            "native_lower_upper_csv_identical": True,
            "embedded_short_alias_rejected": True,
            "native_rows": metadata["telemetry"]["rows"],
            "execution_space": metadata["telemetry"]["execution_space"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
