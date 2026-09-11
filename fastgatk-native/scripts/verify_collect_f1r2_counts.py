#!/usr/bin/env python3
"""Verify native CollectF1R2Counts against the bundled GATK contract."""

from __future__ import annotations

import collections
import gzip
import json
import os
import shutil
import struct
import subprocess
import sys
import tarfile
import tempfile
import zlib
from pathlib import Path
import oracle_guard


def archive_member(path: Path, suffix: str) -> str:
    with tarfile.open(path, "r:gz") as archive:
        members = [member for member in archive.getmembers() if member.name.endswith(suffix)]
        assert len(members) == 1, (path, suffix, [member.name for member in archive.getmembers()])
        return archive.extractfile(members[0]).read().decode("utf-8")


def histogram_counts(path: Path, suffix: str) -> dict[str, int]:
    lines = archive_member(path, suffix).splitlines()
    header = next(line.split("\t")[1:] for line in lines if line.startswith("depth\t"))
    counts: dict[str, int] = collections.Counter()
    for line in lines:
        fields = line.split("\t")
        if fields and fields[0].isdigit():
            for label, value in zip(header, fields[1:]):
                counts[label] += int(value)
    return dict(counts)


def alt_rows(path: Path) -> list[str]:
    return [line for line in archive_member(path, ".alt_table").splitlines()
            if line and not line.startswith("#") and not line.startswith("context\t")]


def _bgzf_block(payload: bytes) -> bytes:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    block_size = 18 + len(compressed) + 8
    assert block_size <= 65536, block_size
    header = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff" +
              struct.pack("<H", 6) + b"BC" + struct.pack("<H", 2) +
              struct.pack("<H", block_size - 1))
    footer = struct.pack("<II", zlib.crc32(payload) & 0xffffffff, len(payload) & 0xffffffff)
    return header + compressed + footer


def rewrite_bam_with_reordered_references(source: Path, destination: Path) -> None:
    """Create a valid sequential BAM whose reference header order is reversed.

    The fixture has only a small number of references, so a byte-level BAM
    rewrite keeps this regression self-contained without requiring samtools or
    pysam. Record core.tid values are remapped to the new header order.
    """
    raw = gzip.decompress(source.read_bytes())
    assert raw[:4] == b"BAM\1"
    text_length = struct.unpack_from("<i", raw, 4)[0]
    cursor = 8 + text_length
    reference_count = struct.unpack_from("<i", raw, cursor)[0]
    cursor += 4
    references: list[tuple[bytes, int]] = []
    for _ in range(reference_count):
        name_length = struct.unpack_from("<i", raw, cursor)[0]
        cursor += 4
        name = raw[cursor:cursor + name_length]
        cursor += name_length
        length = struct.unpack_from("<i", raw, cursor)[0]
        cursor += 4
        references.append((name, length))
    remap = {old: reference_count - old - 1 for old in range(reference_count)}
    output = bytearray(raw[:8 + text_length])
    output.extend(struct.pack("<i", reference_count))
    for name, length in reversed(references):
        output.extend(struct.pack("<i", len(name)))
        output.extend(name)
        output.extend(struct.pack("<i", length))
    while cursor < len(raw):
        block_size = struct.unpack_from("<i", raw, cursor)[0]
        block_start = cursor + 4
        block_end = block_start + block_size
        block = bytearray(raw[block_start:block_end])
        old_tid = struct.unpack_from("<i", block, 0)[0]
        if old_tid >= 0:
            struct.pack_into("<i", block, 0, remap[old_tid])
        output.extend(struct.pack("<i", block_size))
        output.extend(block)
        cursor = block_end
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    encoded = bytearray()
    cursor = 0
    # Keep each BGZF payload comfortably below the 64 KiB block limit.
    while cursor < len(output):
        chunk_size = min(60000, len(output) - cursor)
        while True:
            try:
                block = _bgzf_block(bytes(output[cursor:cursor + chunk_size]))
            except AssertionError:
                chunk_size //= 2
                continue
            if len(block) <= 65536:
                break
            chunk_size //= 2
        encoded.extend(block)
        cursor += chunk_size
    destination.write_bytes(bytes(encoded) + eof)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_COLLECT_F1R2_COUNTS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-collect-f1r2-counts")))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk"
    java = root / "third_party/jdk17/bin/java"
    assert native.exists() and reference.exists() and bam.exists()
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-f1r2-") as directory:
        work = Path(directory)
        native_output = work / "native.tar.gz"
        native_manifest = work / "native.manifest.json"
        native_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
            "-O", str(native_output), "--output-manifest", str(native_manifest), "--threads", "2",
            "--batch-records", "128",
        ], env=env, text=True, capture_output=True)
        assert native_run.returncode == 0, native_run.stderr[-4000:]
        native_summary = json.loads(native_run.stdout.splitlines()[-1])
        assert native_summary["records_seen"] == 493
        assert native_summary["accepted_loci"] == 683
        assert native_output.exists() and native_output.stat().st_size > 0
        with tarfile.open(native_output, "r:gz") as archive:
            assert sorted(member.name for member in archive.getmembers()) == [
                "NA12878.alt_histogram", "NA12878.alt_table", "NA12878.ref_histogram"]
        native_manifest_data = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert native_manifest_data["compatibility"]["standard_tar_members"] is True
        assert native_manifest_data["compatibility"]["bam_pileup"] is True
        assert native_summary["indexed_inputs"] == 1
        assert native_summary["sequential_inputs"] == 0
        assert native_summary["iterator_intervals"] == 1
        assert native_manifest_data["telemetry"]["indexed_traversal"] is True
        telemetry = native_manifest_data["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] >= 2
        assert telemetry["kernel_observations"] == native_summary["observations"]
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        assert telemetry["batch_records"] == 128
        assert telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink")
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] == telemetry["kernel_batches"]
        assert telemetry["pipeline_decoded_bytes"] > 0
        assert telemetry["pipeline_computed_bytes"] > 0
        assert telemetry["pipeline_encoded_bytes"] > 0
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0

        # Barclay exposes the short-name form as --I/--O (in addition to
        # -I/-O).  This spelling is common in generated GATK workflows and
        # must remain a direct native replacement.  Compare the complete
        # archive payload semantically so tar/gzip metadata does not obscure
        # the CLI compatibility boundary.
        alias_output = work / "alias.tar.gz"
        alias_run = subprocess.run([
            str(native), "-R", str(reference), "--I", str(bam), "-L", "17:69000-70000",
            "--O", str(alias_output), "--threads", "2", "--batch-records", "128",
        ], env=env, text=True, capture_output=True)
        assert alias_run.returncode == 0, alias_run.stderr[-4000:]
        alias_summary = json.loads(alias_run.stdout.splitlines()[-1])
        assert alias_summary["records_seen"] == native_summary["records_seen"]
        assert alias_summary["observations"] == native_summary["observations"]
        for suffix in (".ref_histogram", ".alt_histogram"):
            assert histogram_counts(alias_output, suffix) == histogram_counts(native_output, suffix)
        assert alt_rows(alias_output) == alt_rows(native_output)

        # Removing the sidecar index must remain a supported, explicit
        # sequential fallback rather than silently changing the count contract.
        unindexed_bam = work / "unindexed.bam"
        shutil.copyfile(bam, unindexed_bam)
        unindexed_output = work / "unindexed.tar.gz"
        unindexed_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(unindexed_bam), "-L", "17:69000-70000",
            "-O", str(unindexed_output), "--threads", "1", "--batch-records", "128",
        ], env=env, text=True, capture_output=True)
        assert unindexed_run.returncode == 0, unindexed_run.stderr[-4000:]
        unindexed_summary = json.loads(unindexed_run.stdout.splitlines()[-1])
        assert unindexed_summary["records_seen"] == 493
        assert unindexed_summary["indexed_inputs"] == 0
        assert unindexed_summary["sequential_inputs"] == 1
        assert unindexed_summary["iterator_intervals"] == 0

        # A second input with the same contig names in the opposite header
        # order must merge into the same canonical (sample, contig, position)
        # loci.  The rewritten BAM is intentionally left without a .bai so the
        # test covers indexed + sequential mixed ingestion as well.
        reordered_bam = work / "reordered-header.bam"
        rewrite_bam_with_reordered_references(bam, reordered_bam)
        reordered_output = work / "reordered-merged.tar.gz"
        reordered_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-I", str(reordered_bam),
            "-L", "17:69000-70000", "-O", str(reordered_output), "--threads", "1",
            "--batch-records", "128",
        ], env=env, text=True, capture_output=True)
        assert reordered_run.returncode == 0, reordered_run.stderr[-4000:]
        reordered_summary = json.loads(reordered_run.stdout.splitlines()[-1])
        assert reordered_summary["records_seen"] == native_summary["records_seen"] * 2
        assert reordered_summary["observations"] == native_summary["observations"] * 2
        assert reordered_summary["loci"] == native_summary["loci"]
        assert reordered_summary["indexed_inputs"] == 1
        assert reordered_summary["sequential_inputs"] == 1
        assert reordered_summary["iterator_intervals"] == 1

        oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
        if oracle_guard.oracle_ready('verify_collect_f1r2_counts.py', gatk, java):
            oracle_env = env.copy()
            oracle_env["JAVA_HOME"] = str(java.parent.parent)
            oracle_env["PATH"] = str(java.parent) + os.pathsep + oracle_env.get("PATH", "")
            gatk_output = work / "gatk.tar.gz"
            gatk_run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "CollectF1R2Counts",
                "-R", str(reference), "--I", str(bam), "-L", "17:69000-70000",
                "--O", str(gatk_output),
            ], env=oracle_env, text=True, capture_output=True)
            assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
            for suffix in (".ref_histogram", ".alt_histogram"):
                expected = histogram_counts(gatk_output, suffix)
                actual = histogram_counts(native_output, suffix)
                assert expected == actual, (suffix, [(key, expected[key], actual[key])
                                                     for key in expected if expected[key] != actual[key]][:20])
            assert alt_rows(gatk_output) == alt_rows(native_output)

            # -XL affects the site stream, rather than excluding whole reads.
            # Compare member contents (not tar/gzip bytes) so this covers the
            # observable GATK archive boundary while retaining deterministic
            # native member ordering.
            excluded_native = work / "native-excluded.tar.gz"
            excluded_java = work / "gatk-excluded.tar.gz"
            excluded_manifest = work / "native-excluded.manifest.json"
            excluded_native_run = subprocess.run([
                str(native), "-R", str(reference), "--I", str(bam),
                "-L", "17:69000-70000", "-XL", "17:69500-69600",
                "-O", str(excluded_native), "--output-manifest", str(excluded_manifest),
                "--threads", "2", "--batch-records", "128",
            ], env=env, text=True, capture_output=True)
            assert excluded_native_run.returncode == 0, excluded_native_run.stderr[-4000:]
            excluded_java_run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "CollectF1R2Counts",
                "-R", str(reference), "--I", str(bam), "-L", "17:69000-70000",
                "-XL", "17:69500-69600", "--O", str(excluded_java),
            ], env=oracle_env, text=True, capture_output=True)
            assert excluded_java_run.returncode == 0, excluded_java_run.stderr[-4000:]
            for suffix in (".ref_histogram", ".alt_histogram"):
                assert histogram_counts(excluded_java, suffix) == histogram_counts(excluded_native, suffix)
            assert alt_rows(excluded_java) == alt_rows(excluded_native)
            assert json.loads(excluded_manifest.read_text(encoding="utf-8"))["telemetry"]["excluded_intervals"] == 1
            oracle = {
                "status": "pass", "ref_histogram_labels": 64,
                "alt_histogram_labels": 384, "alt_table_rows": len(alt_rows(native_output)),
                "exclude_interval_archive_boundary": True,
            }
        print(json.dumps({
            "status": "pass", "records_seen": native_summary["records_seen"],
            "observations": native_summary["observations"], "loci": native_summary["loci"],
            "accepted_loci": native_summary["accepted_loci"], "archive_members": 3,
            "gatk_oracle": oracle,
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
