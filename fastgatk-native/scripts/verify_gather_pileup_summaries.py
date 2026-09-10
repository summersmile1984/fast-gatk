#!/usr/bin/env python3
"""Verify dictionary-ordered, sample-safe GatherPileupSummaries."""
from __future__ import annotations

import json
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True, check=False)


def table(sample: str, rows: list[str]) -> str:
    return (
        f"#<METADATA>SAMPLE={sample}\n"
        "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        + ("\n".join(rows) + "\n" if rows else "")
    )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-gather-pileup-summaries"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-pileup-summaries-") as directory:
        work = Path(directory)
        dictionary = work / "reference.dict"
        dictionary.write_text("@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
        first = work / "shard-1.table"
        first.write_text(table("TUMOR", [
            "chr1\t20\t10\t2\t0\t0.10",
            "chr2\t4\t8\t4\t0\t0.20",
        ]), encoding="utf-8")
        second = work / "shard-2.table"
        second.write_text(table("TUMOR", ["chr1\t3\t12\t1\t0\t0.05"]), encoding="utf-8")
        empty = work / "empty.table"
        empty.write_text(table("TUMOR", []), encoding="utf-8")
        output = work / "gathered.table"
        manifest = work / "gathered.manifest.json"
        result = run(binary, [
            "-I", str(first), "-I", str(empty), "-I", str(second),
            "-SD", str(dictionary), "-O", str(output), "--output-manifest", str(manifest),
        ])
        assert result.returncode == 0, result.stderr
        lines = output.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "#<METADATA>SAMPLE=TUMOR"
        # GATK sorts shards by each file's first record, then preserves the
        # row order inside each shard (it does not globally sort by dictionary).
        assert lines[2].startswith("chr1\t3\t12\t1")
        assert lines[3].startswith("chr1\t20\t10\t2")
        assert lines[4].startswith("chr2\t4\t8\t4")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["file_first_record_order"] is True
        assert metadata["telemetry"]["empty_inputs"] == 1
        assert metadata["telemetry"]["output_sites"] == 3

        # GatherPileupSummaries is declared with Barclay short names --I and
        # --O.  Keep those exact spellings working alongside the more common
        # -I/-O and --input/--output aliases used by the native dispatcher.
        alias_output = work / "gathered-alias.table"
        alias_result = run(binary, [
            "--I", str(first), "--I", str(empty), "--I", str(second),
            "-SD", str(dictionary), "--O", str(alias_output),
        ])
        assert alias_result.returncode == 0, alias_result.stderr
        assert alias_output.read_text(encoding="utf-8") == output.read_text(encoding="utf-8")

        duplicate = work / "duplicate.table"
        duplicate.write_text(table("TUMOR", ["chr1\t20\t9\t2\t0\t0.05"]), encoding="utf-8")
        gatk_default = run(binary, [
            "-I", str(first), "-I", str(duplicate), "-SD", str(dictionary),
            "-O", str(work / "invalid.table"),
        ])
        assert gatk_default.returncode == 0, gatk_default.stderr
        assert len(gatk_default.stdout.strip().splitlines()) == 1
        assert json.loads((work / "invalid.table.manifest.json").read_text(encoding="utf-8"))["telemetry"]["duplicate_sites"] == 1

        failed = run(binary, [
            "-I", str(first), "-I", str(duplicate), "-SD", str(dictionary),
            "-O", str(work / "strict-invalid.table"), "--reject-overlaps",
        ])
        assert failed.returncode != 0
        assert "overlapping pileup summaries" in failed.stderr

        allowed = run(binary, [
            "-I", str(first), "-I", str(duplicate), "-SD", str(dictionary),
            "-O", str(work / "allowed.table"), "--allow-overlaps",
        ])
        assert allowed.returncode == 0, allowed.stderr
        assert len(allowed.stdout.strip().splitlines()) == 1
        assert len((work / "allowed.table").read_text(encoding="utf-8").splitlines()) == 5

        # Scattered workflows commonly compress pileup summaries between
        # stages.  Exercise mixed plain/.gz inputs and a compressed output,
        # while retaining the same file-first-record ordering contract.
        compressed_first = work / "shard-1.table.gz"
        compressed_second = work / "shard-2.table.gz"
        with gzip.open(compressed_first, "wt", encoding="utf-8") as stream:
            stream.write(first.read_text(encoding="utf-8"))
        with gzip.open(compressed_second, "wt", encoding="utf-8") as stream:
            stream.write(empty.read_text(encoding="utf-8"))
        compressed_output = work / "gathered.table.gz"
        compressed_manifest = work / "gathered-compressed.manifest.json"
        compressed_result = run(binary, [
            "-I", str(compressed_first), "-I", str(second), "-I", str(compressed_second),
            "-SD", str(dictionary), "-O", str(compressed_output),
            "--output-manifest", str(compressed_manifest),
        ])
        assert compressed_result.returncode == 0, compressed_result.stderr
        with gzip.open(compressed_output, "rt", encoding="utf-8") as stream:
            assert stream.read().splitlines() == lines
        compressed_metadata = json.loads(compressed_manifest.read_text(encoding="utf-8"))
        assert compressed_metadata["compatibility"]["compressed_input_count"] == 2
        assert compressed_metadata["compatibility"]["compressed_output"] is True

        out_of_range = work / "out-of-range.table"
        out_of_range.write_text(table("TUMOR", ["chr1\t101\t1\t1\t0\t0.1"]), encoding="utf-8")
        failed_range = run(binary, [
            "-I", str(out_of_range), "-SD", str(dictionary),
            "-O", str(work / "range-invalid.table"),
        ])
        assert failed_range.returncode != 0
        assert "exceeds sequence dictionary length" in failed_range.stderr

        unknown = work / "unknown.table"
        unknown.write_text(table("TUMOR", ["chrX\t1\t1\t1\t0\t0.1"]), encoding="utf-8")
        failed_unknown = run(binary, [
            "-I", str(unknown), "-SD", str(dictionary),
            "-O", str(work / "unknown-invalid.table"),
        ])
        assert failed_unknown.returncode != 0
        assert "absent from sequence dictionary" in failed_unknown.stderr
        allowed_unknown = run(binary, [
            "-I", str(unknown), "-SD", str(dictionary),
            "-O", str(work / "unknown-allowed.table"),
            "--disable-sequence-dictionary-validation",
        ])
        assert allowed_unknown.returncode == 0, allowed_unknown.stderr

    print(json.dumps({"status": "pass", "output_sites": 3, "empty_inputs": 1}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
