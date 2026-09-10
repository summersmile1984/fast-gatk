#!/usr/bin/env python3
"""Pinned oracle for GenotypeGVCFs posterior FORMAT cleanup.

GenotypeGVCFs always uses ``PREFER_PLS`` internally.  When an upstream gVCF
contains the optional FORMAT/GP and FORMAT/PG posterior annotations, GATK's
allele-subsetting utility removes those vectors before writing the final
genotyped record (they are indexed by the pre-subset ``<NON_REF>`` allele
list).  This test uses deliberately populated GP/PG values and requires the
native aggregate and stream-by-locus paths to match GATK's rows exactly while
retaining the input header declarations.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return stream.read().splitlines()


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in lines(path)
            if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    java = Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK GP-input oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-gp-oracle-") as directory:
        work = Path(directory)
        source = work / "input.g.vcf"
        gatk_output = work / "gatk.g.vcf.gz"
        gatk_baseline = work / "gatk.baseline.vcf.gz"
        source.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled genotype likelihoods>
##FORMAT=<ID=GP,Number=G,Type=Float,Description=Genotype posterior in Phred Scale>
##FORMAT=<ID=PG,Number=G,Type=Float,Description=Genotype prior in Phred Scale>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t69067\t.\tT\tG,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL:GP:PG:GQ\t0/0:10:10,0,0:0,50,100:0,50,100:0,20,80:50\t0/1:10:5,5,0:50,0,50:90,0,10:80,0,20:0
""", encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(source)],
            "GATK IndexFeatureFile")
        common = ["-R", str(reference), "-V", str(source),
                  "--create-output-variant-index", "false"]
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs", *common,
             "-O", str(gatk_output), "--gp-qual"], "GATK GenotypeGVCFs --gp-qual")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs", *common,
             "-O", str(gatk_baseline)], "GATK GenotypeGVCFs baseline")
        gatk_rows = rows(gatk_output)
        baseline_rows = rows(gatk_baseline)
        assert gatk_rows == baseline_rows, "GATK --gp-qual changed this PREFER_PLS output unexpectedly"
        assert len(gatk_rows) == 1
        gatk_format = gatk_rows[0][8].split(":")
        assert gatk_format == ["GT", "AD", "DP", "GQ", "PL"]
        assert "GP" not in gatk_format and "PG" not in gatk_format

        comparisons: dict[str, int] = {}
        for traversal in ("aggregate", "stream-by-locus"):
            native_output = work / f"native.{traversal}.vcf.gz"
            manifest = work / f"native.{traversal}.manifest.json"
            command = [str(native), "-R", str(reference), "-V", str(source),
                       "-O", str(native_output), "--gp-qual",
                       "--gatk-compatible-annotations",
                       "--create-output-variant-index=false",
                       "--output-manifest", str(manifest)]
            if traversal == "stream-by-locus":
                command.append("--stream-by-locus")
            run(command, f"native GenotypeGVCFs {traversal}")
            native_rows = rows(native_output)
            assert native_rows == gatk_rows, {
                "traversal": traversal, "gatk": gatk_rows, "native": native_rows}
            assert native_rows[0][8].split(":") == gatk_format
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            assert metadata["compatibility"]["gatk_annotation_compatibility"] is True
            assert metadata["compatibility"]["posterior_qual_opt_in"] is True
            comparisons[traversal] = len(native_rows)

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "input_samples": 2,
        "input_posterior_formats": ["GP", "PG"],
        "gatk_gp_qual_gate_exact": True,
        "posterior_formats_removed_after_allele_subset": True,
        "aggregate_rows": comparisons["aggregate"],
        "stream_by_locus_rows": comparisons["stream-by-locus"],
        "java_native_rows_exact": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
