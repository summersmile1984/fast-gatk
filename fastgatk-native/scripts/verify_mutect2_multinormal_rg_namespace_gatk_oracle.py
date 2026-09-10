#!/usr/bin/env python3
"""Strict GATK oracle for multi-normal inputs with colliding @RG IDs.

GATK resolves each read's RG tag against the header of its originating -I
file.  This fixture intentionally combines the DREAM normal and tumor_1 BAMs:
both contain @RG ID ``C09DF.1`` but map it to different SM values.  The second
selected normal has no reads in the target interval, which also pins GATK's
zero-count / flat-AF genotype serialization for an empty matched normal.
"""

from __future__ import annotations

import difflib
import gzip
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
SECOND_NORMAL = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "tumor_1.bam"
TUMOR_SAMPLE = "synthetic.challenge.set1.tumor"
NORMAL_SAMPLE = "synthetic.challenge.set1.normal"
SECOND_NORMAL_SAMPLE = "tumor sample"
WINDOW = "20:10004000-10004400"
TARGET = ("20", "10004241", "G", "A")
EMPTY_NORMAL_FORMAT = "0/0:0,0:0.500:0:0,0:0,0:0,0:0,0,0,0"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def text(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream
                if not line.startswith("##GATKCommandLine=")]


def data_rows(lines: list[str]) -> list[str]:
    return [line for line in lines if line and not line.startswith("#")]


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL, SECOND_NORMAL)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError(f"missing multi-normal oracle inputs: {missing}")
        print(json.dumps({"status": "skip", "reason": "missing inputs", "paths": missing}))
        return 0

    common = [
        "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL), "-I", str(SECOND_NORMAL),
        "--tumor-sample", TUMOR_SAMPLE,
        "--normal-sample", NORMAL_SAMPLE,
        "--normal-sample", SECOND_NORMAL_SAMPLE,
        "-L", WINDOW, "--add-output-vcf-command-line", "false",
    ]
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-multinormal-") as directory:
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
            "--native-pair-hmm-threads", "1", "--stream-by-region", "500",
        ])
        assert native_stream_run.returncode == 0, native_stream_run.stderr[-4000:]

        gatk_text = text(gatk_output)
        native_text = text(native_output)
        native_stream_text = text(native_stream_output)
        if native_text != gatk_text:
            diff = "\n".join(difflib.unified_diff(
                gatk_text, native_text, fromfile="gatk", tofile="native", lineterm=""))
            raise AssertionError(f"multi-normal VCF mismatch:\n{diff}")
        if native_stream_text != gatk_text:
            diff = "\n".join(difflib.unified_diff(
                gatk_text, native_stream_text, fromfile="gatk", tofile="native-stream", lineterm=""))
            raise AssertionError(f"multi-normal streaming VCF mismatch:\n{diff}")
        header = next(line for line in native_text if line.startswith("#CHROM")).split("\t")
        assert header[9:] == [NORMAL_SAMPLE, TUMOR_SAMPLE, SECOND_NORMAL_SAMPLE], header
        rows = data_rows(native_text)
        assert len(rows) == 1, rows
        fields = rows[0].split("\t")
        assert tuple(fields[index] for index in (0, 1, 3, 4)) == TARGET, fields
        assert fields[8] == "GT:AD:AF:DP:F1R2:F2R1:FAD:SB", fields[8]
        assert fields[-1] == EMPTY_NORMAL_FORMAT, fields[-1]

    print(json.dumps({
        "status": "pass",
        "fixture": "dream_multinormal_colliding_rg",
        "vcf_header_exact": True,
        "vcf_row_exact": True,
        "normal_samples": [NORMAL_SAMPLE, SECOND_NORMAL_SAMPLE],
        "rg_namespace_per_input": True,
        "empty_normal_genotype_exact": True,
        "streaming_vcf_exact": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
