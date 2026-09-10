#!/usr/bin/env python3
"""Pin GATK's multi-tumor sample-role and FORMAT evidence contract.

Mutect2 treats every input sample not named by --normal-sample as a tumor,
even when --tumor-sample names one preferred sample.  This DREAM window has
real evidence in both tumors at 20:10022820 and an empty second-tumor matrix
at 20:10023442.  It pins the complete GATK VCF, including the joint EventMap
candidate set, per-sample FORMAT columns, and aggregate/stream equivalence.
"""

from __future__ import annotations

import gzip
import difflib
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "downloads" / "reference" / "hs37d5.fa.gz"
TUMOR = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "tumor.bam"
NORMAL = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "normal.bam"
SECOND_TUMOR = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "tumor_1.bam"
SECOND_NORMAL = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "normal_1.bam"
TUMOR_SAMPLE = "synthetic.challenge.set1.tumor"
NORMAL_SAMPLE = "synthetic.challenge.set1.normal"
SECOND_TUMOR_SAMPLE = "tumor sample"
WINDOW = "20:10022000-10023500"
FORMAT_FIELDS = ("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB")
EVIDENCE_SITE = ("20", "10022820", "A", "T")
EMPTY_SITE = ("20", "10023442", "G", "A")


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream
                if not line.startswith("##GATKCommandLine=")]


def rows_by_key(content: list[str]) -> dict[tuple[str, str, str, str], list[str]]:
    result: dict[tuple[str, str, str, str], list[str]] = {}
    for line in content:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        result[tuple(fields[index] for index in (0, 1, 3, 4))] = fields
    return result


def evidence_fields(row: list[str], sample_index: int) -> dict[str, str]:
    format_ids = row[8].split(":")
    values = row[9 + sample_index].split(":")
    return {field: values[format_ids.index(field)] for field in FORMAT_FIELDS}


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL, SECOND_TUMOR,
                SECOND_NORMAL)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError(f"missing multi-tumor oracle inputs: {missing}")
        print(json.dumps({"status": "skip", "reason": "missing inputs", "paths": missing}))
        return 0

    common = [
        "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(SECOND_TUMOR),
        "--tumor-sample", TUMOR_SAMPLE, "-L", WINDOW,
        "--add-output-vcf-command-line", "false",
    ]
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-multitumor-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        native_stream_output = work / "native-stream.vcf.gz"
        gatk_run = run([
            str(JAVA), "-jar", str(GATK), "Mutect2", *common,
            "-O", str(gatk_output), "--QUIET", "true",
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
        native_run = run([
            str(NATIVE), *common, "-O", str(native_output),
            "--native-pair-hmm-threads", "1",
        ])
        assert native_run.returncode == 0, native_run.stderr[-4000:]
        native_stream_run = run([
            str(NATIVE), *common, "-O", str(native_stream_output),
            "--native-pair-hmm-threads", "1", "--stream-by-region", "2000",
        ])
        assert native_stream_run.returncode == 0, native_stream_run.stderr[-4000:]

        gatk = lines(gatk_output)
        native = lines(native_output)
        native_stream = lines(native_stream_output)
        if native != gatk:
            diff = "\n".join(difflib.unified_diff(
                gatk, native, fromfile="gatk", tofile="native", lineterm=""))
            raise AssertionError(f"multi-tumor aggregate VCF mismatch:\n{diff}")
        if native_stream != gatk:
            diff = "\n".join(difflib.unified_diff(
                gatk, native_stream, fromfile="gatk", tofile="native-stream", lineterm=""))
            raise AssertionError(f"multi-tumor streaming VCF mismatch:\n{diff}")
        gatk_header = next(line for line in gatk if line.startswith("#CHROM")).split("\t")
        native_header = next(line for line in native if line.startswith("#CHROM")).split("\t")
        assert native_header == gatk_header, (gatk_header, native_header)
        assert native_header[9:] == [TUMOR_SAMPLE, SECOND_TUMOR_SAMPLE]

        gatk_rows = rows_by_key(gatk)
        native_rows = rows_by_key(native)
        for key in (EVIDENCE_SITE, EMPTY_SITE):
            assert key in gatk_rows and key in native_rows, (key, gatk_rows.keys(), native_rows.keys())
            for sample_index in range(2):
                native_evidence = evidence_fields(native_rows[key], sample_index)
                gatk_evidence = evidence_fields(gatk_rows[key], sample_index)
                assert native_evidence.pop("GT") == gatk_evidence.pop("GT"), (
                    key, sample_index, gatk_rows[key], native_rows[key])
                assert native_evidence == gatk_evidence, (key, sample_index,
                                                          gatk_rows[key], native_rows[key])
        # The first site proves the second tumor's column is scored from its
        # own non-zero likelihood matrix, rather than being a copied/empty
        # placeholder.  The second pins GATK's flat-AF empty-matrix behavior.
        assert evidence_fields(native_rows[EVIDENCE_SITE], 1)["AD"] == "19,14"
        assert evidence_fields(native_rows[EMPTY_SITE], 1) == {
            "GT": "0/1", "AD": "0,0", "AF": "0.500", "DP": "0",
            "F1R2": "0,0", "F2R1": "0,0", "FAD": "0,0", "SB": "0,0,0,0",
        }

        # Keep the same two tumor matrices while adding a selected matched
        # normal.  The ActivityProfile must pool both tumor IDs, apply the
        # normal gate only to the normal ID, and preserve every sample-local
        # fragment identity through joint assembly and streamed replay.
        matched_common = [
            "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "-I", str(SECOND_TUMOR), "--tumor-sample", TUMOR_SAMPLE,
            "--normal-sample", NORMAL_SAMPLE, "-L", WINDOW,
            "--add-output-vcf-command-line", "false",
        ]
        gatk_matched = work / "gatk-matched.vcf.gz"
        native_matched = work / "native-matched.vcf.gz"
        native_matched_stream = work / "native-matched-stream.vcf.gz"
        for command, output in (
            ([str(JAVA), "-jar", str(GATK), "Mutect2", *matched_common,
              "-O", str(gatk_matched), "--QUIET", "true"], gatk_matched),
            ([str(NATIVE), *matched_common, "-O", str(native_matched),
              "--native-pair-hmm-threads", "1"], native_matched),
            ([str(NATIVE), *matched_common, "-O", str(native_matched_stream),
              "--native-pair-hmm-threads", "1", "--stream-by-region", "2000"],
             native_matched_stream),
        ):
            completed = run(command)
            assert completed.returncode == 0, completed.stderr[-4000:]
        gatk_matched_lines = lines(gatk_matched)
        native_matched_lines = lines(native_matched)
        native_matched_stream_lines = lines(native_matched_stream)
        assert native_matched_lines == gatk_matched_lines, "matched multi-tumor VCF mismatch"
        assert native_matched_stream_lines == gatk_matched_lines, (
            "matched multi-tumor streaming VCF mismatch")
        matched_header = next(line for line in native_matched_lines
                              if line.startswith("#CHROM")).split("\t")
        assert matched_header[9:] == [NORMAL_SAMPLE, TUMOR_SAMPLE, SECOND_TUMOR_SAMPLE]

        # The same normal SM can arrive through several BAM inputs.  GATK's
        # callRegion() overlap correction groups them by sample and read name,
        # including duplicate same-end records from the two inputs.  This
        # changes the complete MBQ median at the evidence site, so pin both
        # aggregate and bounded streaming output to the Java oracle.
        repeated_normal_common = [
            "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "-I", str(SECOND_TUMOR), "-I", str(SECOND_NORMAL),
            "--tumor-sample", TUMOR_SAMPLE, "--normal-sample", NORMAL_SAMPLE,
            "-L", WINDOW, "--add-output-vcf-command-line", "false",
        ]
        gatk_repeated_normal = work / "gatk-repeated-normal.vcf.gz"
        native_repeated_normal = work / "native-repeated-normal.vcf.gz"
        native_repeated_normal_stream = work / "native-repeated-normal-stream.vcf.gz"
        for command in (
            [str(JAVA), "-jar", str(GATK), "Mutect2", *repeated_normal_common,
             "-O", str(gatk_repeated_normal), "--QUIET", "true"],
            [str(NATIVE), *repeated_normal_common, "-O", str(native_repeated_normal),
             "--native-pair-hmm-threads", "1"],
            [str(NATIVE), *repeated_normal_common, "-O", str(native_repeated_normal_stream),
             "--native-pair-hmm-threads", "1", "--stream-by-region", "2000"],
        ):
            completed = run(command)
            assert completed.returncode == 0, completed.stderr[-4000:]
        gatk_repeated_normal_lines = lines(gatk_repeated_normal)
        native_repeated_normal_lines = lines(native_repeated_normal)
        native_repeated_normal_stream_lines = lines(native_repeated_normal_stream)
        assert native_repeated_normal_lines == gatk_repeated_normal_lines, (
            "repeated-normal multi-input VCF mismatch")
        assert native_repeated_normal_stream_lines == gatk_repeated_normal_lines, (
            "repeated-normal multi-input streaming VCF mismatch")
        repeated_normal_rows = rows_by_key(native_repeated_normal_lines)
        assert "MBQ=29,32" in repeated_normal_rows[EVIDENCE_SITE][7].split(";"), (
            repeated_normal_rows[EVIDENCE_SITE])

    print(json.dumps({
        "status": "pass",
        "fixture": "dream_multitumor_role_and_format",
        "gatk_vcf_exact": True,
        "gatk_sample_header_exact": True,
        "gatk_per_sample_format_exact": True,
        "second_tumor_real_evidence_exact": True,
        "empty_tumor_matrix_exact": True,
        "matched_normal_multitumor_vcf_exact": True,
        "repeated_normal_multi_input_vcf_exact": True,
        "aggregate_stream_vcf_exact": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
