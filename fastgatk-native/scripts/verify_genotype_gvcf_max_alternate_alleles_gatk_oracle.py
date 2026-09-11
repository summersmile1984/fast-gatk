#!/usr/bin/env python3
"""Pinned GATK/native oracle for GenotypeGVCFs maximum-ALT subsetting."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY", str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")
    ))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK max-alternate-alleles oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    # Seven independent concrete ALTs exercise the GenotypingEngine reduction
    # boundary.  Each sample's best PL genotype names one ALT; the seventh ALT
    # is dropped at --max-alternate-alleles 6, making that sample a no-call.
    alts = ["C", "G", "T", "AA", "AC", "AG", "AT", "<NON_REF>"]
    allele_count = 1 + len(alts)

    def rank(first: int, second: int) -> int:
        if first > second:
            first, second = second, first
        return second * (second + 1) // 2 + first

    width = allele_count * (allele_count + 1) // 2
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-max-alt-oracle-") as directory:
        work = Path(directory)
        input_path = work / "cohort.g.vcf"
        header = ("##fileformat=VCFv4.2\n"
                  "##contig=<ID=17,length=1000000>\n"
                  "##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n"
                  "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
                  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
                  "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Quality>\n"
                  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
                  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(f"S{sample}" for sample in range(1, 8)) + "\n")
        samples: list[str] = []
        for sample in range(1, 8):
            pl = [99] * width
            pl[rank(sample, sample)] = 0
            ad = [0] * allele_count
            ad[sample] = 20
            samples.append("./.:20:" + ",".join(map(str, ad)) + ":" +
                           ",".join(map(str, pl)))
        row = ("17\t69000\t.\tA\t" + ",".join(alts) +
               "\t.\tPASS\tDP=140\tGT:DP:AD:PL\t" + "\t".join(samples) + "\n")
        input_path.write_text(header + row, encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)],
            "GATK IndexFeatureFile")

        gatk_output = work / "gatk.vcf.gz"
        run([str(java), "-jar", str(gatk), "GenotypeGVCFs", "-R", str(reference),
             "-V", str(input_path), "-O", str(gatk_output),
             "--max-alternate-alleles", "6", "--stand-call-conf", "0",
             "--annotate-with-num-discovered-alleles", "true",
             "--create-output-variant-index", "false"], "GATK GenotypeGVCFs")
        gatk_rows = records(gatk_output)
        assert len(gatk_rows) == 1
        gatk_row = gatk_rows[0]
        assert gatk_row[4] == "C,G,T,AA,AC,AG"
        assert gatk_row[8] == "GT:AD:DP:GQ:PL"
        assert "NDA=7" in gatk_row[7].split(";")
        assert gatk_row[-1].startswith("./.:0,0,0,0,0,0,0:20:.:0,0")

        results: dict[str, bool] = {}
        for traversal in ("aggregate", "stream-by-locus"):
            native_output = work / f"native-{traversal}.vcf.gz"
            manifest = work / f"native-{traversal}.manifest.json"
            command = [str(binary), "-R", str(reference), "-V", str(input_path),
                       "-O", str(native_output), "--gatk-compatible-annotations",
                       "--max-alternate-alleles", "6", "--stand-call-conf", "0",
                       "--annotate-with-num-discovered-alleles", "true",
                       "--create-output-variant-index=false", "--output-manifest",
                       str(manifest)]
            if traversal == "stream-by-locus":
                command.append("--stream-by-locus")
            run(command, f"native GenotypeGVCFs {traversal}")
            native_rows = records(native_output)
            if native_rows != gatk_rows:
                raise AssertionError({"traversal": traversal, "gatk": gatk_rows,
                                      "native": native_rows})
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            compatibility = metadata["compatibility"]
            telemetry = metadata["telemetry"]
            assert compatibility["max_alternate_alleles_likelihood_subset"] is True
            assert compatibility["gatk_annotation_compatibility"] is True
            assert telemetry["max_alternate_alleles"] == 6
            assert telemetry["max_alt_pruning_calls"] == 1
            assert telemetry["max_alt_alleles_pruned"] == 1
            assert telemetry["max_alt_score_kernel_calls"] == 1
            assert telemetry["max_alt_score_kernel_execution_space"]
            assert telemetry["annotate_with_num_discovered_alleles"] is True
            assert telemetry["num_discovered_alleles_calls"] == 1
            results[traversal] = True

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "max_alternate_alleles": 6,
            "input_concrete_alt_count": 7,
            "output_concrete_alt_count": 6,
            "java_native_rows_exact": True,
            "aggregate_rows_exact": results["aggregate"],
            "stream_rows_exact": results["stream-by-locus"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
