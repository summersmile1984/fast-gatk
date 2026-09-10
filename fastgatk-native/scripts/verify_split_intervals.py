#!/usr/bin/env python3
"""Contract/oracle test for native SplitIntervals scatter output."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "gatk-source/src/test/resources/hg19micro.fasta"
BINARY = Path(os.environ.get("FASTGATK_SPLIT_INTERVALS_BINARY", ROOT / "fastgatk-native/build/fastgatk-split-intervals"))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def body(path: Path) -> list[tuple[str, int, int]]:
    result: list[tuple[str, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        result.append((fields[0], int(fields[1]), int(fields[2])))
    return result


def shard_files(directory: Path) -> list[Path]:
    return sorted(directory.glob("*-scattered.interval_list"))


def run_java(directory: Path, *extra: str) -> None:
    result = subprocess.run(
        [str(JAVA), "-jar", str(JAR), "SplitIntervals", "-R", str(REFERENCE),
         "-O", str(directory), *extra], text=True, capture_output=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)


def main() -> int:
    if not BINARY.is_file():
        print(json.dumps({"status": "skip", "suite": "split-intervals", "reason": f"missing binary: {BINARY}"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-split-intervals-") as directory:
        work = Path(directory)
        native_dir = work / "native"
        java_dir = work / "java"
        manifest = work / "native.manifest.json"
        native = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-O", str(native_dir), "--scatter-count", "3",
             "--output-manifest", str(manifest)], text=True, capture_output=True,
        )
        assert native.returncode == 0, (native.stdout, native.stderr)
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        assert payload["output_shards"] == 3 and payload["bases"] == 32000
        native_files = shard_files(native_dir)
        assert len(native_files) == 3
        run_java(java_dir, "--scatter-count", "3")
        java_files = shard_files(java_dir)
        assert len(java_files) == 3
        assert [body(path) for path in native_files] == [body(path) for path in java_files]
        assert native_files[0].read_text(encoding="utf-8") == java_files[0].read_text(encoding="utf-8")

        # Interval-count mode keeps interval records intact and distributes
        # the remainder deterministically; compare the semantic body to Java.
        count_native = work / "count-native"
        count_java = work / "count-java"
        native_count = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-L", "1:1-100", "-L", "2:1-50",
             "-O", str(count_native), "--scatter-count", "2", "--subdivision-mode", "INTERVAL_COUNT"],
            text=True, capture_output=True,
        )
        assert native_count.returncode == 0, (native_count.stdout, native_count.stderr)
        run_java(count_java, "-L", "1:1-100", "-L", "2:1-50", "--scatter-count", "2",
                 "--subdivision-mode", "INTERVAL_COUNT")
        assert [body(path) for path in shard_files(count_native)] == [body(path) for path in shard_files(count_java)]

        # Exercise every Picard scatter mode on a non-uniform interval list,
        # including scatter-count > interval-count and the overflow policy.
        mode_intervals = work / "modes.interval_list"
        mode_intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:1\tLN:16000\n@SQ\tSN:2\tLN:16000\n"
            "1\t1\t100\t+\ta\n1\t201\t500\t+\tb\n2\t1\t200\t+\tc\n",
            encoding="utf-8",
        )
        for mode in (
            "INTERVAL_SUBDIVISION",
            "INTERVAL_COUNT",
            "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER",
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
        ):
            mode_native = work / f"mode-native-{mode}"
            mode_java = work / f"mode-java-{mode}"
            native_mode = subprocess.run(
                [str(BINARY), "-R", str(REFERENCE), "-L", str(mode_intervals),
                 "-O", str(mode_native), "--scatter-count", "8",
                 "--subdivision-mode", mode],
                text=True, capture_output=True,
            )
            assert native_mode.returncode == 0, (mode, native_mode.stdout, native_mode.stderr)
            run_java(mode_java, "-L", str(mode_intervals), "--scatter-count", "8",
                     "--subdivision-mode", mode)
            assert [body(path) for path in shard_files(mode_native)] == [
                body(path) for path in shard_files(mode_java)
            ], mode

        # The same interval-list and BED semantics must survive HTSlib's
        # compressed text reader; compare the shard bodies with the plain
        # files above so compression cannot change coordinate conversion or
        # deterministic scatter ordering.
        compressed_interval = work / "modes.interval_list.gz"
        with gzip.open(compressed_interval, "wt", encoding="utf-8") as stream:
            stream.write(mode_intervals.read_text(encoding="utf-8"))
        compressed_native = work / "compressed-native"
        compressed_result = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-L", str(compressed_interval),
             "-O", str(compressed_native), "--scatter-count", "8",
             "--subdivision-mode", "INTERVAL_COUNT"],
            text=True, capture_output=True,
        )
        assert compressed_result.returncode == 0, (compressed_result.stdout, compressed_result.stderr)
        plain_count_native = work / "mode-native-INTERVAL_COUNT"
        assert [body(path) for path in shard_files(compressed_native)] == [
            body(path) for path in shard_files(plain_count_native)
        ]

        compressed_bed = work / "select.bed.gz"
        with gzip.open(compressed_bed, "wt", encoding="utf-8") as stream:
            stream.write("track name=targets\n1\t0\t100\n2\t0\t50\n")
        compressed_bed_native = work / "compressed-bed-native"
        compressed_bed_result = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-L", str(compressed_bed), "-XL", "1:1-10",
             "-O", str(compressed_bed_native), "--scatter-count", "2"],
            text=True, capture_output=True,
        )
        assert compressed_bed_result.returncode == 0, (
            compressed_bed_result.stdout, compressed_bed_result.stderr)
        compressed_selected = [item for path in shard_files(compressed_bed_native) for item in body(path)]
        assert compressed_selected == [
            ("1", 11, 80), ("1", 81, 100), ("2", 1, 50)
        ]

        # BED and excluded intervals exercise the file parser and subtraction
        # boundary used by scatter workflows.
        bed = work / "select.bed"
        bed.write_text("1\t0\t100\n2\t0\t50\n", encoding="utf-8")
        bed_native = work / "bed-native"
        bed_result = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-L", str(bed), "-XL", "1:1-10",
             "-O", str(bed_native), "--scatter-count", "2"],
            text=True, capture_output=True,
        )
        assert bed_result.returncode == 0, (bed_result.stdout, bed_result.stderr)
        selected = [item for path in shard_files(bed_native) for item in body(path)]
        assert selected == [("1", 11, 80), ("1", 81, 100), ("2", 1, 50)]

        mixed_native = work / "mixed-native"
        mixed_result = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-O", str(mixed_native), "--scatter-count", "3",
             "--dont-mix-contigs"], text=True, capture_output=True,
        )
        assert mixed_result.returncode == 0, (mixed_result.stdout, mixed_result.stderr)
        assert all(len({row[0] for row in body(path)}) <= 1 for path in shard_files(mixed_native))

        # --min-contig-size is applied to reference-derived intervals only.
        # With the fixture's 16,000-base contigs, a threshold of 16,001
        # leaves a deterministic empty shard rather than silently retaining
        # short contigs.
        filtered_native = work / "filtered-native"
        filtered_manifest = work / "filtered.manifest.json"
        filtered_result = subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-O", str(filtered_native),
             "--scatter-count", "1", "--min-contig-size", "16001",
             "--output-manifest", str(filtered_manifest)],
            text=True, capture_output=True,
        )
        assert filtered_result.returncode == 0, (filtered_result.stdout, filtered_result.stderr)
        filtered_payload = json.loads(filtered_manifest.read_text(encoding="utf-8"))
        assert filtered_payload["intervals"] == 0 and filtered_payload["bases"] == 0
        assert filtered_payload["output_shards"] == 0 and not shard_files(filtered_native)
    print(json.dumps({"status": "pass", "suite": "split-intervals", "scatter_shards": 3,
                      "scatter_modes_checked": 5, "java_oracle": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
