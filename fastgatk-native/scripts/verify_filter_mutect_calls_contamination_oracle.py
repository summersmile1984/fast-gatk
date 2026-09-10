#!/usr/bin/env python3
"""Compare the native contamination posterior path with GATK 4.6.2.0."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_records(path: Path) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = {}
            if fields[7] not in {"", "."}:
                for token in fields[7].split(";"):
                    if "=" in token:
                        key, value = token.split("=", 1)
                        info[key] = value
            result.append({"pos": fields[1], "filter": fields[6], "contq": info.get("CONTQ", "")})
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
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK contamination oracle unavailable")
        print(json.dumps({"status": "skip", "reason": "bundled GATK contamination oracle unavailable"}))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=POPAF,Number=A,Type=Float,Description=Population allele frequency>
##INFO=<ID=PSOMATIC,Number=A,Type=Float,Description=Somatic posterior>
##INFO=<ID=PGERMLINE,Number=A,Type=Float,Description=Germline posterior>
##INFO=<ID=PARTIFACT,Number=A,Type=Float,Description=Artifact posterior>
##INFO=<ID=CONTQ,Number=1,Type=Float,Description=Phred-scaled qualities that alt allele are not due to contamination>
##INFO=<ID=MBQ,Number=1,Type=Integer,Description=Median alternate base quality>
##INFO=<ID=MMQ,Number=R,Type=Integer,Description=Median reference/alternate mapping quality>
##INFO=<ID=MPOS,Number=1,Type=Integer,Description=Median alternate read position>
##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
    # The first record is compatible with contamination, the second has a
    # strong alt count at a low population AF and should be filtered.
    body = (
        "17\t69000\t.\tA\tG\t.\tPASS\tTLOD=8;POPAF=2;PSOMATIC=1;PGERMLINE=0;PARTIFACT=0\tGT:AD\t0/1:90,10\n"
        "17\t69001\t.\tA\tG\t.\tPASS\tTLOD=8;POPAF=2;PSOMATIC=1;PGERMLINE=0;PARTIFACT=0\tGT:AD\t0/1:50,50\n"
        # GATK's argument field uses -1 as a parser sentinel, but the
        # M2FiltersArgumentCollection getter resolves it to an effective 30 in
        # normal mode. This record is a direct-replacement guard: native must
        # add map_qual for the ALT MMQ=29 in the Number=R vector.
        "17\t69002\t.\tA\tG\t.\tPASS\tTLOD=8;POPAF=2;PSOMATIC=1;PGERMLINE=0;PARTIFACT=0;MBQ=20;MMQ=60,29;MPOS=5;MFRL=100,100\tGT:AD\t0/1:50,50\n"
    )
    with tempfile.TemporaryDirectory(prefix="fastgatk-contamination-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        gatk_unmatched_table_output = work / "gatk-unmatched-table.vcf.gz"
        native_unmatched_table_output = work / "native-unmatched-table.vcf.gz"
        gatk_sites_only = work / "gatk-sites-only.vcf.gz"
        native_sites_only = work / "native-sites-only.vcf.gz"
        stats = work / "calls.stats"
        # ContaminationRecord is a per-sample lookup.  This deliberately
        # contains a valid estimate for a different sample: GATK must fall
        # back to its configured default (0.0) for TUMOR rather than applying
        # OTHER's 0.9 estimate globally.
        unmatched_table = work / "unmatched-contamination.table"
        unmatched_table.write_text(
            "sample\tcontamination\terror\nOTHER\t0.90\t0.001\n", encoding="utf-8")
        stats.write_text("#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n", encoding="utf-8")
        input_path.write_text(header + body, encoding="utf-8")
        index_run = subprocess.run([
            str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path),
        ], text=True, capture_output=True, check=False)
        assert index_run.returncode == 0, index_run.stderr
        common = ["--contamination-estimate", "0.1", "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1"]
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
        gatk_sites_only_run = subprocess.run([
            str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
            "-V", str(input_path), "-O", str(gatk_sites_only), "--stats", str(stats),
            "--lenient", "true", "--sites-only-vcf-output", "true", *common,
        ], text=True, capture_output=True, check=False)
        assert gatk_sites_only_run.returncode == 0, gatk_sites_only_run.stderr
        native_sites_only_run = subprocess.run([
            str(binary), "-V", str(input_path), "-O", str(native_sites_only),
            "--stats", str(stats), "--sites-only-vcf-output", "true", *common,
        ], text=True, capture_output=True, check=False)
        assert native_sites_only_run.returncode == 0, native_sites_only_run.stderr
        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert len(gatk_records) == len(native_records) == 3
        assert [entry["pos"] for entry in gatk_records] == [entry["pos"] for entry in native_records]
        assert "contamination" in gatk_records[0]["filter"]
        assert "contamination" in native_records[0]["filter"]
        assert "contamination" not in gatk_records[1]["filter"]
        assert "contamination" not in native_records[1]["filter"]
        assert "map_qual" in gatk_records[2]["filter"]
        assert "map_qual" in native_records[2]["filter"]
        assert native_records[0]["contq"] == "1"
        assert native_records[1]["contq"] == "93"

        # The table path must retain GATK's sample-aware fallback.  Before
        # this guard native reduced all rows to a global maximum, so OTHER's
        # 0.90 incorrectly changed TUMOR's FILTER/CONTQ.
        table_common = [
            "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1",
            "--contamination-table", str(unmatched_table),
        ]
        gatk_table_run = subprocess.run([
            str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
            "-V", str(input_path), "-O", str(gatk_unmatched_table_output),
            "--stats", str(stats), "--lenient", "true", *table_common,
        ], text=True, capture_output=True, check=False)
        assert gatk_table_run.returncode == 0, gatk_table_run.stderr
        native_table_run = subprocess.run([
            str(binary), "-V", str(input_path), "-O", str(native_unmatched_table_output),
            "--stats", str(stats), *table_common,
        ], text=True, capture_output=True, check=False)
        assert native_table_run.returncode == 0, native_table_run.stderr
        gatk_table_records = read_records(gatk_unmatched_table_output)
        native_table_records = read_records(native_unmatched_table_output)
        assert len(gatk_table_records) == len(native_table_records) == 3
        assert [entry["pos"] for entry in gatk_table_records] == [
            entry["pos"] for entry in native_table_records]
        assert ["contamination" in entry["filter"] for entry in gatk_table_records] == [
            "contamination" in entry["filter"] for entry in native_table_records]
        assert "contamination" not in gatk_table_records[0]["filter"]
        assert gatk_table_records[0]["contq"] in {"", "93"}
        def assert_sites_only(path: Path) -> None:
            with gzip.open(path, "rt", encoding="utf-8") as stream:
                lines = stream.read().splitlines()
            header_line = next(line for line in lines if line.startswith("#CHROM"))
            assert len(header_line.split("\t")) == 8
            records = [line for line in lines if line and not line.startswith("#")]
            assert len(records) == 3 and all(len(line.split("\t")) == 8 for line in records)
        assert_sites_only(gatk_sites_only)
        assert_sites_only(native_sites_only)
        print(json.dumps({
            "status": "pass",
            "gatk_records": gatk_records,
            "native_records": native_records,
            "contamination_filter_exact": (
                ("contamination" in gatk_records[0]["filter"] and
                 "contamination" in native_records[0]["filter"]) and
                ("contamination" not in gatk_records[1]["filter"] and
                 "contamination" not in native_records[1]["filter"])
            ),
            "gatk_contq_materialized": any(entry["contq"] for entry in gatk_records),
            "native_contq_materialized": any(entry["contq"] for entry in native_records),
            "contamination_table_unmatched_sample_exact": (
                ["contamination" in entry["filter"] for entry in gatk_table_records] ==
                ["contamination" in entry["filter"] for entry in native_table_records]
            ),
            "contamination_table_uses_sample_lookup": "contamination" not in gatk_table_records[0]["filter"],
            "gatk_default_mmq_effective_30": "map_qual" in gatk_records[2]["filter"],
            "native_default_mmq_effective_30": "map_qual" in native_records[2]["filter"],
            "sites_only_vcf_output_shape_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
