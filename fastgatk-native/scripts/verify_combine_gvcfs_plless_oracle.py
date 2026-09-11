#!/usr/bin/env python3
"""Pinned GATK/native oracle for PL-less CombineGVCFs records."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of a reference block>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Genotype likelihoods>
"""


def write(path: Path, sample: str, body: str) -> None:
    with path.open("wt", encoding="utf-8") as stream:
        stream.write(HEADER)
        stream.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        stream.write(body)


def records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n").split("\t") for line in stream
                if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_COMBINE_GVCFS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-combine-gvcfs")))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.is_file() for path in (java, gatk, binary, reference)):
        oracle_guard.oracle_not_verified('verify_combine_gvcfs_plless_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("PL-less CombineGVCFs oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-combine-plless-oracle-") as directory:
        work = Path(directory)
        first = work / "first.g.vcf"
        second = work / "second.g.vcf"
        gatk_output = work / "gatk.g.vcf.gz"
        native_output = work / "native.g.vcf.gz"
        native_manifest = work / "native.manifest.json"
        native_stream_output = work / "native-stream.g.vcf.gz"
        native_stream_manifest = work / "native-stream.manifest.json"
        body = "17\t69000\t.\tA\t<NON_REF>\t.\tPASS\tEND=69002\tGT:GQ:DP:AD\t0/0:60:10:10,0\n"
        write(first, "SAMPLE1", body)
        write(second, "SAMPLE2", body.replace("10:10,0", "12:12,0").replace(":60:", ":55:"))
        for input_path in (first, second):
            run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)],
                "GATK IndexFeatureFile")
        run([str(java), "-jar", str(gatk), "CombineGVCFs", "-R", str(reference),
             "-V", str(first), "-V", str(second), "-O", str(gatk_output)],
            "GATK CombineGVCFs")
        run([str(binary), "-R", str(reference), "-V", str(first), "-V", str(second),
             "-O", str(native_output), "--output-manifest", str(native_manifest)],
            "native CombineGVCFs")
        run([str(binary), "-R", str(reference), "-V", str(first), "-V", str(second),
             "-O", str(native_stream_output), "--stream-merge",
             "--output-manifest", str(native_stream_manifest)],
            "native CombineGVCFs --stream-merge")
        gatk_rows = records(gatk_output)
        native_rows = records(native_output)
        native_stream_rows = records(native_stream_output)
        if gatk_rows != native_rows or gatk_rows != native_stream_rows:
            raise AssertionError({"gatk": gatk_rows, "native": native_rows,
                                  "native_stream": native_stream_rows})
        assert gatk_rows and "PL" not in gatk_rows[0][8].split(":")
        assert "DP=" not in gatk_rows[0][7]
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        stream_manifest = json.loads(native_stream_manifest.read_text(encoding="utf-8"))
        assert manifest["compatibility"]["gatk_format_semantics"] is True
        assert stream_manifest["compatibility"]["stream_merge"] is True
        print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                          "records": len(gatk_rows), "pl_less_rows_exact": True,
                          "stream_rows_exact": True, "sample_count": 2}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
