#!/usr/bin/env python3
"""Compare native GatherPileupSummaries with bundled GATK output bytes."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def table(rows: list[str]) -> str:
    return (
        "#<METADATA>SAMPLE=TUMOR\n"
        "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        + ("\n".join(rows) + "\n" if rows else "")
    )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_GATHER_PILEUP_BINARY",
        str(root / "fastgatk-native/build/fastgatk-gather-pileup-summaries"),
    ))
    required = (java, gatk, binary)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_gather_pileup_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK GatherPileupSummaries oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-pileup-gatk-oracle-") as directory:
        work = Path(directory)
        dictionary = work / "reference.dict"
        dictionary.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100\n",
            encoding="utf-8",
        )
        first = work / "shard-1.table"
        first.write_text(table([
            "chr1\t20\t10\t2\t0\t0.12345678901234567",
            "chr1\t21\t10\t2\t0\t0.0",
            "chr1\t22\t10\t2\t0\t0.0001",
            "chr1\t23\t10\t2\t0\t0.001",
            "chr2\t4\t8\t4\t0\t0.20",
        ]), encoding="utf-8")
        second = work / "shard-2.table"
        second.write_text(table(["chr1\t3\t12\t1\t0\t0.05"]), encoding="utf-8")
        empty = work / "empty.table"
        empty.write_text(table([]), encoding="utf-8")
        gatk_output = work / "gatk.table"
        native_output = work / "native.table"
        gatk_run = run([
            str(java), "-jar", str(gatk), "GatherPileupSummaries",
            "--I", str(first), "--I", str(empty), "--I", str(second),
            "--sequence-dictionary", str(dictionary), "--O", str(gatk_output),
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = run([
            str(binary), "--I", str(first), "--I", str(empty), "--I", str(second),
            "-SD", str(dictionary), "--O", str(native_output),
        ])
        assert native_run.returncode == 0, native_run.stderr
        expected = gatk_output.read_text(encoding="utf-8")
        actual = native_output.read_text(encoding="utf-8")
        assert actual == expected, {"expected": expected, "actual": actual}

        # Empty shards are removed before GATK validates sample metadata.  An
        # all-empty gather therefore succeeds even with conflicting SAMPLE
        # lines and emits only the table header (no metadata line).
        empty_a = work / "empty-a.table"
        empty_b = work / "empty-b.table"
        empty_a.write_text(table([]).replace("SAMPLE=TUMOR", "SAMPLE=A"), encoding="utf-8")
        empty_b.write_text(table([]).replace("SAMPLE=TUMOR", "SAMPLE=B"), encoding="utf-8")
        gatk_empty = work / "gatk-empty.table"
        native_empty = work / "native-empty.table"
        gatk_empty_run = run([
            str(java), "-jar", str(gatk), "GatherPileupSummaries",
            "-I", str(empty_a), "-I", str(empty_b),
            "--sequence-dictionary", str(dictionary), "-O", str(gatk_empty),
        ])
        assert gatk_empty_run.returncode == 0, gatk_empty_run.stderr
        native_empty_run = run([
            str(binary), "-I", str(empty_a), "-I", str(empty_b),
            "-SD", str(dictionary), "-O", str(native_empty),
        ])
        assert native_empty_run.returncode == 0, native_empty_run.stderr
        assert native_empty.read_text(encoding="utf-8") == gatk_empty.read_text(encoding="utf-8")
        assert native_empty.read_text(encoding="utf-8") == (
            "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        )

        # An empty shard is also ignored when non-empty shards are present;
        # the contributing sample controls the output metadata regardless of
        # input order.
        gatk_mixed = work / "gatk-mixed-empty-first.table"
        native_mixed = work / "native-mixed-empty-first.table"
        gatk_mixed_run = run([
            str(java), "-jar", str(gatk), "GatherPileupSummaries",
            "-I", str(empty_b), "-I", str(first),
            "--sequence-dictionary", str(dictionary), "-O", str(gatk_mixed),
        ])
        assert gatk_mixed_run.returncode == 0, gatk_mixed_run.stderr
        native_mixed_run = run([
            str(binary), "-I", str(empty_b), "-I", str(first),
            "-SD", str(dictionary), "-O", str(native_mixed),
        ])
        assert native_mixed_run.returncode == 0, native_mixed_run.stderr
        assert native_mixed.read_text(encoding="utf-8") == gatk_mixed.read_text(encoding="utf-8")
        assert native_mixed.read_text(encoding="utf-8").startswith("#<METADATA>SAMPLE=TUMOR\n")
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "input_files": 3,
            "empty_inputs": 1,
            "file_order_exact": True,
            "table_bytes_exact": True,
            "all_empty_header_exact": True,
            "empty_sample_conflict_ignored": True,
            "mixed_empty_sample_ignored": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
