#!/usr/bin/env python3
"""Pinned GATK oracle for ReblockGVCF overlap trim/split semantics."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

from verify_reblock_gatk_oracle import read_records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    binary = Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY", str(root / "fastgatk-native/build/fastgatk-reblock-gvcf")
    ))
    required = (java, gatk, reference, binary)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK overlap oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-overlap-oracle-") as directory:
        work = Path(directory)
        source = work / "overlap.g.vcf"
        gatk_output = work / "gatk.g.vcf"
        native_output = work / "native.g.vcf"
        native_manifest = work / "native.manifest.json"
        source.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
17\t69000\t.\tT\t<NON_REF>\t.\t.\tEND=69010\tGT:DP:GQ:MIN_DP:PL\t0/0:10:20:10:0,20,200
17\t69005\t.\tT\tG,<NON_REF>\t50\t.\tDP=10\tGT:DP:AD:PL:GQ\t0/1:10:5,5,0:50,0,80,99,99,99:50
""",
            encoding="utf-8",
        )
        java_run = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(gatk_output),
            "--add-output-vcf-command-line", "false",
        ], text=True, capture_output=True, check=False)
        assert java_run.returncode == 0, java_run.stderr
        native_run = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source),
            "-O", str(native_output), "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr

        java_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert native_records == java_records
        assert len(native_records) == 3
        assert native_records[0]["site"][1] == "69000"
        assert native_records[0]["info"]["END"] == "69004"
        assert native_records[1]["site"][1] == "69005"
        assert native_records[2]["site"][1] == "69006"
        assert native_records[2]["info"]["END"] == "69010"
        metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["overlapping_ref_block_trim_split"] is True
        assert metadata["telemetry"]["overlapping_ref_block_split"] == 1

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "overlapping_ref_block_split_exact": True,
        "output_records": 3,
        "coverage": ["69000-69004", "69005", "69006-69010"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
