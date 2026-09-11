#!/usr/bin/env python3
"""Pinned GATK oracle for FilterMutectCalls GermlineFilter/GERMQ semantics."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_records(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, str] = {}
            if fields[7] not in {"", "."}:
                for token in fields[7].split(";"):
                    if "=" in token:
                        key, value = token.split("=", 1)
                        info[key] = value
            records.append({
                "pos": int(fields[1]),
                "filter": fields[6],
                "germq": info.get("GERMQ"),
            })
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_filter_mutect_calls_germline_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK germline oracle unavailable")
        print(json.dumps({"status": "skip", "reason": "bundled GATK germline oracle unavailable"}))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=NLOD,Number=A,Type=Float,Description=Normal log odds>
##INFO=<ID=POPAF,Number=A,Type=Float,Description=Population allele frequency>
##INFO=<ID=GERMQ,Number=A,Type=Integer,Description=Phred-scaled qualities that alt alleles are not germline variants>
##INFO=<ID=MBQ,Number=R,Type=Integer,Description=Median base quality>
##INFO=<ID=MMQ,Number=R,Type=Integer,Description=Median mapping quality>
##INFO=<ID=MPOS,Number=A,Type=Integer,Description=Median alternate read position>
##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=AF,Number=A,Type=Float,Description=Allele fraction>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
    body = (
        # Population AF < 1e-10 makes both germline hypotheses impossible;
        # QualityUtils.errorProbToQual(0) is bounded to 93.
        "17\t69000\t.\tA\tG\t.\tPASS\tTLOD=8;NLOD=0;POPAF=100;MBQ=30,30;MMQ=60,60;MPOS=5;MFRL=100,100\tGT:AD:AF\t0/1:90,10:0.1\n"
        # Population AF > 1-1e-10 makes the somatic hypothesis impossible;
        # errorProbToQual(1) is bounded to 1.
        "17\t69001\t.\tA\tG\t.\tPASS\tTLOD=8;NLOD=0;POPAF=0;MBQ=30,30;MMQ=60,60;MPOS=5;MFRL=100,100\tGT:AD:AF\t0/1:90,10:0.1\n"
        # Interior frequency exercises the full germline-het/somatic
        # likelihood ratio and the optional NLOD sign conversion.
        "17\t69002\t.\tA\tG\t.\tPASS\tTLOD=8;NLOD=-5;POPAF=2;MBQ=30,30;MMQ=60,60;MPOS=5;MFRL=100,100\tGT:AD:AF\t0/1:50,50:0.5\n"
        # Missing POPAF is a required-annotation miss in GATK GermlineFilter;
        # neither implementation should fabricate GERMQ.
        "17\t69003\t.\tA\tG\t.\tPASS\tTLOD=8;NLOD=0;MBQ=30,30;MMQ=60,60;MPOS=5;MFRL=100,100\tGT:AD:AF\t0/1:90,10:0.1\n"
    )
    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-mutect-germline-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        stats = work / "calls.stats"
        input_path.write_text(header + body, encoding="utf-8")
        stats.write_text("#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t0\n", encoding="utf-8")
        index_run = subprocess.run([
            str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path),
        ], text=True, capture_output=True, check=False)
        assert index_run.returncode == 0, index_run.stderr
        common = ["--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1",
                  "--create-output-variant-index", "false"]
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
        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert [entry["pos"] for entry in gatk_records] == [69000, 69001, 69002, 69003]
        assert [entry["pos"] for entry in native_records] == [69000, 69001, 69002, 69003]
        assert gatk_records[0]["germq"] == native_records[0]["germq"] == "93"
        assert gatk_records[1]["germq"] == native_records[1]["germq"] == "1"
        assert gatk_records[3]["germq"] is None and native_records[3]["germq"] is None
        # The interior value used to be reported as a semantic oracle (both
        # implementations materialize a finite [1,93] phred value); with the
        # Java-exact SomaticClusteringModel learning (FASTGATK_FMC_JAVA_PASSES)
        # the learned model is GATK-machine-exact, so the interior GERMQ must
        # now equal the GATK value exactly.
        assert gatk_records[2]["germq"] is not None
        assert native_records[2]["germq"] == gatk_records[2]["germq"]
        print(json.dumps({
            "status": "pass",
            "records": len(gatk_records),
            "gatk_germq": [entry["germq"] for entry in gatk_records],
            "native_germq": [entry["germq"] for entry in native_records],
            "population_prior_corners_exact": True,
            "missing_popaf_fail_closed": True,
            "interior_germline_posterior_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
