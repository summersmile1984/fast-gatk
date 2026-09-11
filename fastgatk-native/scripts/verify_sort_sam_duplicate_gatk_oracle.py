#!/usr/bin/env python3
"""Pinned Picard/GATK 4.6.2.0 oracle for SortSam duplicate order.

This exercises HTSJDK's SAMRecordDuplicateComparator boundary: library order,
unclipped read/mate coordinates (including MC), pair orientation, mapped-end
precedence, duplicate-score tie breaking, and deterministic read-name order.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


SAM = """@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:200
@SQ\tSN:chr2\tLN:200
@RG\tID:rgA\tSM:S1\tLB:alpha
@RG\tID:rgZ\tSM:S1\tLB:zeta
zeta-pair\t99\tchr1\t40\t60\t10M\t=\t90\t60\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgZ\tMC:Z:10M
alpha-low\t99\tchr1\t12\t60\t1H2S8M\t=\t50\t48\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA\tMC:Z:10M
unknown-fragment\t0\tchr2\t20\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII
alpha-high\t99\tchr1\t12\t60\t1H2S10M\t=\t50\t50\tAAAAAAAAAAAA\tIIIIIIIIIIII\tRG:Z:rgA\tMC:Z:10M
alpha-low\t147\tchr1\t50\t60\t10M\t=\t12\t-48\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA\tMC:Z:1H2S8M
alpha-reverse\t83\tchr1\t70\t60\t8M2S1H\t=\t30\t-48\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA\tMC:Z:10M
alpha-forward\t99\tchr1\t30\t60\t10M\t=\t70\t48\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA\tMC:Z:8M2S1H
zeta-pair\t147\tchr1\t90\t60\t10M\t=\t40\t-60\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgZ\tMC:Z:10M
alpha-high\t147\tchr1\t50\t60\t10M\t=\t12\t-50\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA\tMC:Z:1H2S10M
alpha-unpaired\t0\tchr1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA
alpha-unmapped-mate\t73\tchr1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rgA
alpha-unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t!!!!\tRG:Z:rgA
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        result.append(tuple(fields[:11]) + tuple(sorted(fields[11:])))
    return result


def sort_order(path: Path) -> str | None:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("@HD\t"):
            for field in line.split("\t")[1:]:
                if field.startswith("SO:"):
                    return field[3:]
    return None


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_SORT_SAM_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-sort-sam"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not (native.is_file() and os.access(native, os.X_OK) and java.is_file() and gatk.is_file()):
        oracle_guard.oracle_not_verified('verify_sort_sam_duplicate_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned SortSam duplicate-order oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-sort-sam-duplicate-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        picard = work / "picard.sam"
        actual = work / "native.sam"
        manifest = work / "native.manifest.json"
        source.write_text(SAM, encoding="utf-8")
        picard_command = [str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
                          "-I", str(source), "-O", str(picard), "-SO", "duplicate",
                          "--CREATE_INDEX", "false"]
        native_command = [str(native), "-I", str(source), "-O", str(actual),
                          "-SO", "duplicate", "--MAX_RECORDS_IN_RAM", "3",
                          "--TMP_DIR", str(work / "spill"),
                          "--output-manifest", str(manifest)]
        picard_result = run(picard_command)
        if picard_result.returncode != 0:
            raise AssertionError(picard_result.stderr)
        native_result = run(native_command)
        if native_result.returncode != 0:
            raise AssertionError(native_result.stderr)
        actual_records = records(actual)
        if records(picard) != actual_records:
            raise AssertionError({"picard": records(picard), "native": actual_records})
        assert sort_order(picard) == "duplicate" and sort_order(actual) == "duplicate"
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["sort_order"] == "duplicate"
        assert metadata["compatibility"]["bam_index"] is False
        assert metadata["telemetry"]["spill_runs"] == 4

    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "sort_order": "duplicate", "records": len(actual_records),
                      "spill_runs": 4}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
