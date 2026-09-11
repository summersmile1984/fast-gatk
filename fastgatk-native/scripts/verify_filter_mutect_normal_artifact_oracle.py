#!/usr/bin/env python3
"""Pin the NormalArtifactFilter posterior threshold to GATK 4.6.2.0.

The normal-artifact posterior is an ARTIFACT error probability and must be
compared with FilterMutectCalls' learned/constant error threshold.  A
positive posterior below a 0.1 threshold is not a filter.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[tuple[str, str, str]]:
    result: list[tuple[str, str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = {}
            for token in fields[7].split(";"):
                if "=" in token:
                    key, value = token.split("=", 1)
                    info[key] = value
            result.append((fields[1], fields[6], info.get("NALOD", "")))
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    required = (java, gatk, reference, binary)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_filter_mutect_normal_artifact_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK normal-artifact oracle unavailable")
        print(json.dumps({"status": "skip", "reason": "bundled GATK normal-artifact oracle unavailable"}))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##normal_sample=NORMAL
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=NALOD,Number=A,Type=Float,Description=Negative log10 normal artifact odds>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>
##INFO=<ID=MBQ,Number=A,Type=Integer,Description=Median base quality>
##INFO=<ID=MMQ,Number=R,Type=Integer,Description=Median mapping quality>
##INFO=<ID=MPOS,Number=A,Type=Integer,Description=Median read position>
##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>
##INFO=<ID=NCount,Number=1,Type=Integer,Description=Number of N bases>
##INFO=<ID=ECNT,Number=1,Type=Integer,Description=Potential somatic events in assembly region>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
"""
    # The normal AF is safely above 10% of the tumour AF, and the normal
    # pileup p-value is above 0.001.  NALOD=-3 therefore produces a small but
    # positive normal-artifact posterior (~8.92e-3), which must pass at a 0.1
    # constant error threshold.
    body = (
        "17\t69110\t.\tA\tG\t.\tPASS\t"
        "TLOD=30;NALOD=-3;AF=0.1;MBQ=30;MMQ=40,40;MPOS=25;MFRL=100,100;NCount=0;ECNT=1\t"
        "GT:AD\t0/0:98,2\t0/1:90,10\n"
    )
    with tempfile.TemporaryDirectory(prefix="fastgatk-normal-artifact-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        stats = work / "calls.stats"
        input_path.write_text(header + body, encoding="utf-8")
        stats.write_text("#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n", encoding="utf-8")
        index_run = subprocess.run([
            str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path),
        ], text=True, capture_output=True, check=False)
        assert index_run.returncode == 0, index_run.stderr
        common = ["--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1"]
        gatk_run = subprocess.run([
            str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
            "-V", str(input_path), "-O", str(gatk_output), "--stats", str(stats),
            "--lenient", "true", *common,
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(binary), "-V", str(input_path), "-O", str(native_output),
            "--stats", str(stats), *common,
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        gatk_records = records(gatk_output)
        native_records = records(native_output)
        assert len(gatk_records) == len(native_records) == 1
        assert gatk_records[0][0] == native_records[0][0] == "69110"
        # This is intentionally expected to fail until native compares the
        # posterior with the configured error threshold instead of > 0.
        assert gatk_records[0][1] == "PASS", gatk_records
        assert native_records[0][1] == "PASS", native_records
        expected_posterior = 0.9 / (0.9 + 100.0)
        assert 0.0 < expected_posterior < 0.1
        print(json.dumps({
            "status": "pass",
            "gatk_records": gatk_records,
            "native_records": native_records,
            "expected_normal_artifact_posterior": expected_posterior,
            "normal_artifact_posterior_threshold_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
