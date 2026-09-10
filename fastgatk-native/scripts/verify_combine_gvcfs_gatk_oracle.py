#!/usr/bin/env python3
"""Compare the bounded CombineGVCFs materialization against GATK."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=All filters passed>
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of a reference block>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Genotype likelihoods>
##GVCFBlock=<GVCFBlock name="GVCFBlock1" minGQ=0 maxGQ=10>
##GVCFBlock=<GVCFBlock name="GVCFBlock2" minGQ=11 maxGQ=20>
##GVCFBlock=<GVCFBlock name="GVCFBlock3" minGQ=21 maxGQ=30>
##GVCFBlock=<GVCFBlock name="GVCFBlock4" minGQ=31 maxGQ=40>
##GVCFBlock=<GVCFBlock name="GVCFBlock5" minGQ=41 maxGQ=50>
##GVCFBlock=<GVCFBlock name="GVCFBlock6" minGQ=51 maxGQ=60>
##GVCFBlock=<GVCFBlock name="GVCFBlock7" minGQ=61 maxGQ=70>
##GVCFBlock=<GVCFBlock name="GVCFBlock8" minGQ=71 maxGQ=80>
##GVCFBlock=<GVCFBlock name="GVCFBlock9" minGQ=81 maxGQ=90>
##GVCFBlock=<GVCFBlock name="GVCFBlock10" minGQ=91 maxGQ=99>
"""


def write(path: Path, sample: str, body: str) -> None:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "wt", encoding="utf-8") as handle:
        handle.write(HEADER)
        handle.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        handle.write(body)


def records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                # INFO/FORMAT order is release-specific; compare stable site,
                # ALT, and sample payload fields after the common columns.
                result.append((fields[0], fields[1], fields[3], fields[4], fields[8], *fields[9:]))
    return result


def site_only_records(path: Path) -> tuple[list[str], list[tuple[str, ...]]]:
    """Read the writer-boundary shape without assuming FORMAT columns."""
    headers: list[str] = []
    result: list[tuple[str, ...]] = []
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                headers.append(line.rstrip("\n"))
            elif line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                result.append(tuple(fields))
    return headers, result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_COMBINE_GVCFS_BINARY", str(root / "fastgatk-native/build/fastgatk-combine-gvcfs")
    ))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK CombineGVCFs oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-combine-gvcfs-gatk-oracle-") as directory:
        work = Path(directory)
        first = work / "first.g.vcf"
        second = work / "second.g.vcf"
        write(first, "SAMPLE1",
              "17\t69000\t.\tA\t<NON_REF>\t.\tPASS\tEND=69002\tGT:GQ:DP:AD:PL\t0/0:60:10:10,0:0,99,99\n"
              "17\t69003\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT:GQ:DP:AD:PL\t0/1:50:20:12,8,0:80,0,80,99,99,99\n")
        write(second, "SAMPLE2",
              "17\t69000\t.\tA\t<NON_REF>\t.\tPASS\tEND=69002\tGT:GQ:DP:AD:PL\t0/0:55:12:12,0:0,99,99\n"
              "17\t69003\t.\tA\tT,<NON_REF>\t.\tPASS\t.\tGT:GQ:DP:AD:PL\t1/1:45:18:0,0,18:99,99,0,99,99,99\n")
        for input_path in (first, second):
            index_run = subprocess.run([
                str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path),
            ], text=True, capture_output=True, check=False)
            assert index_run.returncode == 0, index_run.stderr
        gatk_output = work / "gatk.g.vcf.gz"
        native_output = work / "native.g.vcf.gz"
        native_manifest = work / "native.manifest.json"
        gatk_run = subprocess.run([
            str(java), "-jar", str(gatk), "CombineGVCFs", "-R", str(reference),
            "-V", str(first), "-V", str(second), "-O", str(gatk_output),
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(binary), "-V", str(first), "-V", str(second), "-O", str(native_output),
            "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        native_metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert native_metadata["determinism"] == "strict"
        assert native_metadata["compatibility"]["gatk_format_semantics"] is True
        assert native_metadata["compatibility"]["native_genotype_materialization"] is False
        assert native_metadata["compatibility"]["bit_identical_to_gatk"] is True
        gatk_records = records(gatk_output)
        native_records = records(native_output)
        gatk_sites = {(row[0], row[1], row[2]) for row in gatk_records}
        native_sites = {(row[0], row[1], row[2]) for row in native_records}
        assert gatk_sites == native_sites, {"gatk": gatk_records, "native": native_records}
        assert gatk_records == native_records, {
            "gatk": gatk_records, "native": native_records,
        }
        assert Path(f"{gatk_output}.tbi").is_file() and Path(f"{native_output}.tbi").is_file()

        # --call-genotypes is the public GATK spelling for materializing GT;
        # compare this mode separately because GATK preserves an input GQ
        # when present while deriving GT from PL.  The native path must retain
        # that boundary and still emit the same site/sample payload.
        gatk_called = work / "gatk-called.g.vcf.gz"
        native_called = work / "native-called.g.vcf.gz"
        gatk_called_run = subprocess.run([
            str(java), "-jar", str(gatk), "CombineGVCFs", "-R", str(reference),
            "-V", str(first), "-V", str(second), "-O", str(gatk_called),
            "--call-genotypes",
        ], text=True, capture_output=True, check=False)
        assert gatk_called_run.returncode == 0, gatk_called_run.stderr
        native_called_run = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(first), "-V", str(second),
            "-O", str(native_called), "--call-genotypes",
        ], text=True, capture_output=True, check=False)
        assert native_called_run.returncode == 0, native_called_run.stderr
        assert records(gatk_called) == records(native_called), {
            "gatk_called": records(gatk_called),
            "native_called": records(native_called),
        }
        assert Path(f"{gatk_called}.tbi").is_file() and Path(f"{native_called}.tbi").is_file()

        # GATK's sites-only writer removes all sample/FORMAT columns after
        # combining.  Compare the exact 8-column site shape and stable site
        # payload while allowing release-specific metadata headers to differ.
        gatk_sites_only = work / "gatk-sites-only.g.vcf.gz"
        native_sites_only = work / "native-sites-only.g.vcf.gz"
        native_sites_only_manifest = work / "native-sites-only.manifest.json"
        gatk_sites_only_run = subprocess.run([
            str(java), "-jar", str(gatk), "CombineGVCFs", "-R", str(reference),
            "-V", str(first), "-V", str(second), "-O", str(gatk_sites_only),
            "--sites-only-vcf-output", "true",
        ], text=True, capture_output=True, check=False)
        assert gatk_sites_only_run.returncode == 0, gatk_sites_only_run.stderr
        native_sites_only_run = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(first), "-V", str(second),
            "-O", str(native_sites_only), "--sites-only-vcf-output=true",
            "--output-manifest", str(native_sites_only_manifest),
        ], text=True, capture_output=True, check=False)
        assert native_sites_only_run.returncode == 0, native_sites_only_run.stderr
        gatk_site_header, gatk_site_rows = site_only_records(gatk_sites_only)
        native_site_header, native_site_rows = site_only_records(native_sites_only)
        assert gatk_site_header and native_site_header
        assert gatk_site_header[0].split("\t") == native_site_header[0].split("\t") == [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
        ]
        assert gatk_site_rows == native_site_rows, {
            "gatk_sites_only": gatk_site_rows,
            "native_sites_only": native_site_rows,
        }
        native_sites_only_metadata = json.loads(native_sites_only_manifest.read_text(encoding="utf-8"))
        assert native_sites_only_metadata["compatibility"]["sites_only_vcf_output"] is True
        assert native_sites_only_metadata["telemetry"]["sites_only_vcf_output"] is True
        assert Path(f"{gatk_sites_only}.tbi").is_file() and Path(f"{native_sites_only}.tbi").is_file()
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(native_records),
            "site_set_exact": True,
            "sample_union": 2,
            "indexed": True,
            "full_sample_payload_bit_identical": True,
            "call_genotypes_payload_bit_identical": True,
            "sites_only_vcf_output_shape_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
