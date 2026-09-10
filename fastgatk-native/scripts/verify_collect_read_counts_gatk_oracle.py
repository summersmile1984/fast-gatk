#!/usr/bin/env python3
"""Compare the native CollectReadCounts TSV boundary with pinned GATK."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_READ_COUNTS_BINARY", str(BUILD / "fastgatk-collect-read-counts")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def data_rows(path: pathlib.Path) -> list[tuple[str, int, int, int]]:
    rows: list[tuple[str, int, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@") or line.startswith("CONTIG"):
            continue
        fields = line.split("\t")
        rows.append((fields[0], int(fields[1]), int(fields[2]), int(float(fields[3]))))
    return rows


def run_json(command: list[str], env: dict[str, str]) -> dict:
    result = subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    return json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    if not NATIVE.exists():
        raise SystemExit(f"missing native binary: {NATIVE}")
    if not JAVA.exists() or not JAR.exists():
        raise SystemExit("missing pinned GATK 4.6.2.0 runtime")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-read-counts-gatk-") as directory:
        work = pathlib.Path(directory)
        sam = work / "reads.sam"
        intervals = work / "targets.interval_list"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:SAMPLE\n"
            # read starts before the interval but its alignment span overlaps it;
            # this catches the common (and incorrect) read-start-only shortcut.
            "span\t0\tchr1\t5\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n"
            "inside\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "outside\t0\tchr1\t21\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "duplicate\t1024\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "lowmq\t0\tchr1\t10\t20\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
            "chr1\t10\t20\t+\ttarget\n",
            encoding="utf-8",
        )
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        # GATK's interval traversal requires an indexed input.  Materialize a
        # coordinate-sorted BAM with the pinned Picard implementation so the
        # comparison exercises the same indexed HTSlib/htsjdk boundary.
        bam = work / "reads.bam"
        sort_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "SortSam", f"I={sam}", f"O={bam}",
            "SORT_ORDER=coordinate", "CREATE_INDEX=true", "VALIDATION_STRINGENCY=LENIENT",
        ], text=True, capture_output=True, env=env)
        if sort_result.returncode != 0:
            raise AssertionError("Picard SortSam failed:\n" + sort_result.stderr[-4000:])
        bam_index = pathlib.Path(str(bam) + ".bai")
        if not bam_index.exists():
            bam_index = bam.with_suffix(".bai")
        if not bam.exists() or not bam_index.exists():
            raise AssertionError("Picard SortSam did not create indexed BAM")
        fastq = work / "reads.fastq"
        fastq_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "SamToFastq", f"I={bam}", f"FASTQ={fastq}",
        ], text=True, capture_output=True, env=env)
        if fastq_result.returncode != 0:
            raise AssertionError("Picard BAM decode failed:\n" + fastq_result.stderr[-8000:])
        # Verify the materialized BAM really contains mapped records before
        # attributing any disagreement to CollectReadCounts itself.
        count_binary = BUILD / "fastgatk-count-reads"
        count_probe_text = ""
        if count_binary.exists():
            count_probe = subprocess.run([str(count_binary), "-I", str(bam)],
                                         text=True, capture_output=True, env=env)
            if count_probe.returncode != 0:
                raise AssertionError("native BAM probe failed:\n" + count_probe.stderr[-4000:])
            count_probe_text = count_probe.stdout
            if count_probe_text.strip() == "0":
                raise AssertionError("Picard BAM contains no records; FASTQ bytes="
                                     + str(fastq.stat().st_size if fastq.exists() else 0)
                                     + "\n" + sort_result.stderr[-8000:])
        native_output = work / "native.tsv"
        native_manifest = work / "native.manifest.json"
        java_output = work / "java.tsv"
        java_tmp = work / "java-tmp"
        native_tmp = work / "native-tmp"
        java_tmp.mkdir()
        native_tmp.mkdir()
        native = run_json([
            str(NATIVE), "-I", str(bam), "-L", str(intervals), "-O", str(native_output),
            "--format", "TSV", "--sample", "SAMPLE",
            "--output-manifest", str(native_manifest), "--threads", "2",
            "--QUIET", "true", "--tmp-dir", str(native_tmp),
            "--read-validation-stringency", "SILENT",
            "--disable-bam-index-caching", "false",
            "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
            "--verbosity", "WARNING",
        ], env)
        java_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "CollectReadCounts", "-I", str(bam),
            "-L", str(intervals), "-O", str(java_output), "--format", "TSV",
            "--interval-merging-rule", "OVERLAPPING_ONLY", "--QUIET", "true",
            "--tmp-dir", str(java_tmp), "--read-validation-stringency", "SILENT",
            "--disable-bam-index-caching", "false", "--use-jdk-deflater", "false",
            "--use-jdk-inflater", "false", "--verbosity", "WARNING",
        ], text=True, capture_output=True, env=env)
        if java_result.returncode != 0:
            raise AssertionError("GATK CollectReadCounts failed:\n" + java_result.stderr[-4000:])
        native_rows = data_rows(native_output)
        java_rows = data_rows(java_output)
        if native_rows != java_rows:
            raise AssertionError({"native_rows": native_rows, "java_rows": java_rows,
                                  "native_text": native_output.read_text(encoding="utf-8"),
                                  "java_text": java_output.read_text(encoding="utf-8"),
                                  "bam_probe": count_probe_text})
        if native_rows != [("chr1", 10, 20, 1)]:
            raise AssertionError({"unexpected_rows": native_rows})
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        if manifest["compatibility"]["bit_identical_to_gatk"] is not False:
            raise AssertionError(manifest)
        if manifest["telemetry"]["reads_used"] != 1:
            raise AssertionError(manifest)
        if manifest["compatibility"]["quiet"] is not True or manifest["telemetry"]["quiet"] is not True:
            raise AssertionError(manifest)
        if manifest["compatibility"]["read_validation_stringency"] != "SILENT" or \
                manifest["telemetry"]["read_validation_stringency"] != "SILENT":
            raise AssertionError(manifest)
        # GATK's CNV tool validates OVERLAPPING_ONLY even though the generic
        # interval help advertises ALL as its default; preserve that observed
        # release-specific contract and verify overlap union explicitly.
        adjacent = work / "adjacent.interval_list"
        adjacent.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
            "chr1\t1\t5\t+\ta\nchr1\t6\t8\t+\tb\n", encoding="utf-8")
        adjacent_native = work / "adjacent.tsv"
        adjacent_java = work / "adjacent-java.tsv"
        run_json([
            str(NATIVE), "-I", str(bam), "-L", str(adjacent), "-O", str(adjacent_native),
            "--format", "TSV", "--interval-merging-rule", "OVERLAPPING_ONLY",
        ], env)
        adjacent_java_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "CollectReadCounts", "-I", str(bam),
            "-L", str(adjacent), "-O", str(adjacent_java), "--format", "TSV",
            "--interval-merging-rule", "OVERLAPPING_ONLY",
        ], text=True, capture_output=True, env=env)
        if adjacent_java_result.returncode != 0:
            raise AssertionError("GATK adjacency oracle failed:\n" + adjacent_java_result.stderr[-4000:])
        if data_rows(adjacent_native) != data_rows(adjacent_java):
            raise AssertionError({"native_adjacent": data_rows(adjacent_native),
                                  "java_adjacent": data_rows(adjacent_java)})

        # -XL is applied to the output target set, not as a whole-read
        # filter.  The native writer must therefore split the original target
        # at the excluded span and retain read starts in either remaining
        # piece.  This is a pinned Java command boundary, including the
        # single-dash Barclay form.
        excluded_native = work / "excluded.tsv"
        excluded_java = work / "excluded-java.tsv"
        excluded_manifest = work / "excluded.manifest.json"
        excluded_native_summary = run_json([
            str(NATIVE), "-I", str(bam), "-L", "chr1:1-20", "-XL", "chr1:9-12",
            "-O", str(excluded_native), "--format", "TSV", "-imr", "OVERLAPPING_ONLY",
            "--output-manifest", str(excluded_manifest),
        ], env)
        excluded_java_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "CollectReadCounts", "-I", str(bam),
            "-L", "chr1:1-20", "-XL", "chr1:9-12", "-O", str(excluded_java),
            "--format", "TSV", "-imr", "OVERLAPPING_ONLY",
        ], text=True, capture_output=True, env=env)
        if excluded_java_result.returncode != 0:
            raise AssertionError("GATK CollectReadCounts exclusion oracle failed:\n" +
                                 excluded_java_result.stderr[-4000:])
        excluded_native_rows = data_rows(excluded_native)
        excluded_java_rows = data_rows(excluded_java)
        if excluded_native_rows != excluded_java_rows:
            raise AssertionError({"native_excluded": excluded_native_rows,
                                  "java_excluded": excluded_java_rows})
        if excluded_native_rows != [("chr1", 1, 8, 1), ("chr1", 13, 20, 0)]:
            raise AssertionError({"unexpected_excluded_rows": excluded_native_rows})
        excluded_meta = json.loads(excluded_manifest.read_text(encoding="utf-8"))
        if excluded_native_summary["excluded_intervals"] != 1 or \
                excluded_meta["telemetry"]["excluded_intervals"] != 1:
            raise AssertionError(excluded_meta)

        for option in ("-ip", "-ixp"):
            expected_error = "Interval exclusion padding must be set to 0" if option == "-ixp" \
                else "Interval padding must be set to 0"
            native_reject = subprocess.run([
                str(NATIVE), "-I", str(bam), "-L", "chr1:1-20", "-O", str(work / f"native-{option}.tsv"),
                "--format", "TSV", "-imr", "OVERLAPPING_ONLY", option, "1",
            ], text=True, capture_output=True, env=env)
            if native_reject.returncode == 0 or expected_error not in native_reject.stderr:
                raise AssertionError({"option": option, "native": native_reject.stderr})
            java_reject = subprocess.run([
                str(JAVA), "-jar", str(JAR), "CollectReadCounts", "-I", str(bam),
                "-L", "chr1:1-20", "-O", str(work / f"java-{option}.tsv"), "--format", "TSV",
                "-imr", "OVERLAPPING_ONLY", option, "1",
            ], text=True, capture_output=True, env=env)
            if java_reject.returncode == 0 or expected_error not in java_reject.stderr:
                raise AssertionError({"option": option, "java": java_reject.stderr[-4000:]})
        print(json.dumps({
            "status": "pass", "gatk_version": "4.6.2.0", "records_compared": len(native_rows),
            "read_start_semantics": True, "default_filter_semantics": True,
            "adjacent_interval_rule": "OVERLAPPING_ONLY", "exclude_intervals": True,
            "gatk_utility_cli_controls": True,
            "manifest_execution_space": manifest["execution_space"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
