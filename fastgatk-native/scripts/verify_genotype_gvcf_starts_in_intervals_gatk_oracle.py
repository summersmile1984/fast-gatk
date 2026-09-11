#!/usr/bin/env python3
"""Pinned GATK oracle for GenotypeGVCFs writer-side STARTS_IN filtering.

The long deletion starts immediately before ``-L`` but overlaps it.  GATK
must traverse and genotype that record in the default mode, then remove it
only when the deprecated ``--only-output-calls-starting-in-intervals`` writer
mode is enabled.  A SNP whose POS is inside the interval remains.  Both native
aggregate and stream-by-locus paths are compared to the GATK 4.6.2.0 rows.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-5000:]}")


def rows(path: Path) -> list[list[str]]:
    with path.open("rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def write_fixture(work: Path) -> tuple[Path, Path]:
    reference = work / "reference.fasta"
    sequence = "ACGT" * 25
    reference.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
    (work / "reference.fasta.fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    (work / "reference.dict").write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")

    source = work / "input.g.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        "##ALT=<ID=NON_REF,Description=Represents any possible alternative allele>\n"
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=Approximate read depth>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Approximate read depth>\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred likelihoods>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t10\t.\tCGT\tC,<NON_REF>\t.\t.\tDP=30\tGT:AD:DP:GQ:PL\t./.:0,30,0:30:99:200,0,300,200,300,300\n"
        "chr1\t12\t.\tT\tA,<NON_REF>\t.\t.\tDP=25\tGT:AD:DP:GQ:PL\t./.:8,17,0:25:99:180,0,250,180,250,250\n",
        encoding="utf-8")
    return reference, source


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    java = Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not all(path.is_file() for path in (native, java, gatk)):
        oracle_guard.oracle_not_verified('verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GenotypeGVCFs STARTS_IN oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-starts-in-oracle-") as directory:
        work = Path(directory)
        reference, source = write_fixture(work)
        run([str(java), "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile",
             "-I", str(source)], "GATK IndexFeatureFile")
        region = "chr1:11-12"
        common = ["-R", str(reference), "-V", str(source), "-L", region,
                  "--standard-min-confidence-threshold-for-calling", "0",
                  "--create-output-variant-index", "false"]
        java_overlap = work / "java.overlap.vcf"
        java_starts = work / "java.starts.vcf"
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs", *common,
             "-O", str(java_overlap)], "GATK GenotypeGVCFs default interval writer")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs", *common,
             "--only-output-calls-starting-in-intervals", "true",
             "-O", str(java_starts)], "GATK GenotypeGVCFs STARTS_IN writer")

        overlap = rows(java_overlap)
        expected = rows(java_starts)
        assert [int(row[1]) for row in overlap] == [10, 12], overlap
        assert [int(row[1]) for row in expected] == [12], expected

        for traversal in ("aggregate", "stream"):
            output = work / f"native.{traversal}.vcf"
            manifest = work / f"native.{traversal}.manifest.json"
            command = [str(native), *common,
                       "--only-output-calls-starting-in-intervals"]
            if traversal == "aggregate":
                command.append("true")
            command += ["--gatk-compatible-annotations", "-O", str(output),
                       "--output-manifest", str(manifest)]
            if traversal == "stream":
                command.append("--stream-by-locus")
            run(command, f"native GenotypeGVCFs STARTS_IN {traversal}")
            observed = rows(output)
            assert observed == expected, {
                "traversal": traversal, "java": expected, "native": observed}
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            assert metadata["compatibility"]["output_interval_starts_in"] is True
            assert metadata["telemetry"]["only_output_calls_starting_in_intervals"] is True
            assert metadata["telemetry"]["output_interval_skipped"] == 1
            assert metadata["telemetry"]["output_records"] == 1

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "interval": region,
        "default_positions": [10, 12],
        "starts_in_positions": [12],
        "java_native_rows_exact": True,
        "aggregate_and_stream": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
