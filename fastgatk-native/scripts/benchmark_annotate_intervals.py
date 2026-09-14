#!/usr/bin/env python3
"""Reproducible file-boundary benchmark for native AnnotateIntervals."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
import time


def write_reference(path: pathlib.Path, fai: pathlib.Path, length: int) -> None:
    width = 80
    sequence = ("ACGT" * ((length + 3) // 4))[:length]
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(">chr1\n")
        for offset in range(0, length, width):
            stream.write(sequence[offset:offset + width] + "\n")
    # FASTA offsets: the first base starts immediately after '>chr1\n'.
    fai.write_text(f"chr1\t{length}\t6\t{width}\t{width + 1}\n", encoding="ascii")
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:{length}\n", encoding="ascii",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", type=int, default=10000)
    parser.add_argument("--interval-size", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--with-tracks", action="store_true",
        help="also benchmark length-weighted mappability/segmental BED tracks",
    )
    parser.add_argument(
        "--include-java", action="store_true",
        help="also time the pinned GATK Java implementation when available",
    )
    args = parser.parse_args()
    if args.intervals < 1 or args.interval_size < 1 or args.threads < 1:
        raise SystemExit("intervals/interval-size/threads must be positive")

    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_ANNOTATE_INTERVALS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-annotate-intervals"),
    ))
    reference_length = args.intervals * args.interval_size
    with tempfile.TemporaryDirectory(prefix="fastgatk-annotate-intervals-benchmark-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fa"
        write_reference(reference, work / "reference.fa.fai", reference_length)
        interval_list = work / "targets.interval_list"
        with interval_list.open("w", encoding="ascii") as stream:
            stream.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:" + str(reference_length) + "\n")
            for index in range(args.intervals):
                start = index * args.interval_size + 1
                stream.write(f"chr1\t{start}\t{start + args.interval_size - 1}\t+\tI{index}\n")
        mappability = work / "mappability.bed"
        segmental = work / "segmental.bed"
        if args.with_tracks:
            # Keep rows non-overlapping and sorted like production tracks.  A
            # row per interval deliberately exercises the per-contig binary
            # search index in the Kokkos kernel instead of the old full scan.
            with mappability.open("w", encoding="ascii") as stream:
                for index in range(args.intervals):
                    start = index * args.interval_size
                    end = min(start + args.interval_size, reference_length)
                    if end - start > 1:
                        stream.write(f"chr1\t{start}\t{end}\tmap{index}\t0.75\n")
            with segmental.open("w", encoding="ascii") as stream:
                for index in range(0, args.intervals, 3):
                    start = index * args.interval_size
                    end = min(start + args.interval_size, reference_length)
                    if end - start > 1:
                        stream.write(f"chr1\t{start}\t{end}\tseg{index}\t0.25\n")
        output = work / "annotated.tsv"
        manifest = work / "annotated.manifest.json"
        command = [str(binary), "-R", str(reference), "-L", str(interval_list),
                   "--interval-merging-rule", "OVERLAPPING_ONLY",
                   "-O", str(output), "--threads", str(args.threads),
                   "--output-manifest", str(manifest)]
        if args.with_tracks:
            command.extend(["--mappability-track", str(mappability),
                            "--segmental-duplication-track", str(segmental)])
        begin = time.perf_counter()
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        elapsed = time.perf_counter() - begin
        summary = json.loads(completed.stdout.splitlines()[-1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        result = {
            "status": "pass",
            "tool": "AnnotateIntervals",
            "intervals": args.intervals,
            "interval_size": args.interval_size,
            "threads": args.threads,
            "wall_seconds": elapsed,
            "intervals_per_second": args.intervals / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "output_records": summary["intervals"],
            "tracks": args.with_tracks,
            "kernel": {
                "execution_space": telemetry["kernel_execution_space"],
                "records": telemetry["kernel_records"],
                "prepare_seconds": telemetry["kernel_prepare_seconds"],
                "execute_seconds": telemetry["kernel_execute_seconds"],
            },
        }
        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if args.include_java and java.is_file() and jar.is_file():
            if args.with_tracks:
                for track in (mappability, segmental):
                    subprocess.run([
                        str(java), "-jar", str(jar), "IndexFeatureFile", "-I", str(track),
                    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            java_output = work / "annotated-java.tsv"
            java_command = [
                str(java), "-Xmx1g", "-jar", str(jar), "AnnotateIntervals",
                "-R", str(reference), "-L", str(interval_list),
                "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(java_output),
            ]
            if args.with_tracks:
                java_command.extend([
                    "--mappability-track", str(mappability),
                    "--segmental-duplication-track", str(segmental),
                ])
            java_started = time.perf_counter()
            subprocess.run(java_command, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            java_wall = time.perf_counter() - java_started
            byte_identical = output.read_bytes() == java_output.read_bytes()
            if not byte_identical:
                raise AssertionError("native/Java AnnotateIntervals benchmark output differs")
            result["java_baseline"] = {
                "implementation": "gatk-java-4.6.2.0",
                "wall_seconds": java_wall,
                "output_bytes": java_output.stat().st_size,
                "byte_identical": byte_identical,
                "speedup_native_over_java": java_wall / elapsed if elapsed else None,
            }
        elif args.include_java:
            result["java_baseline"] = {"status": "skip", "reason": "pinned Java/GATK unavailable"}
        print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
