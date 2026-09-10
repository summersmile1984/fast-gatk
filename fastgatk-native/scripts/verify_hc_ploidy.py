#!/usr/bin/env python3
"""Regression for HaplotypeCaller reference-backed arbitrary sample ploidy."""
from __future__ import annotations

import json
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def parse_records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()
            if line and not line.startswith("#")]


def parse_records_gzip(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-hc-call"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not binary.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing HC binary or bundled fixture")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-ploidy-") as directory:
        work = Path(directory)
        reports: dict[int, dict[str, object]] = {}
        for ploidy in (1, 3, 8):
            output = work / f"ploidy-{ploidy}.vcf"
            manifest = work / f"ploidy-{ploidy}.manifest.json"
            subprocess.run([
                str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
                "-O", str(output), "--threads", "2", "--min-depth", "1",
                "--min-alt-support", "1", "--sample-ploidy", str(ploidy),
                "--standard-min-confidence-threshold-for-calling", "0",
                "--output-manifest", str(manifest),
            ], check=True, stdout=subprocess.DEVNULL)
            records = parse_records(output)
            assert records, f"ploidy {ploidy} emitted no variant records"
            record = next(row for row in records if row[1] == "69067")
            format_keys = record[8].split(":")
            sample_values = dict(zip(format_keys, record[9].split(":")))
            genotype = sample_values["GT"].split("/")
            assert len(genotype) == ploidy, (ploidy, sample_values["GT"])
            assert len(sample_values["PL"].split(",")) == ploidy + 1
            info = dict(item.split("=", 1) for item in record[7].split(";") if "=" in item)
            assert int(info["AN"]) == ploidy
            metadata = json.loads(manifest.read_text())
            assert metadata["telemetry"]["sample_ploidy"] == ploidy
            assert metadata["compatibility"]["arbitrary_ploidy_vcf"] is (ploidy != 2)
            reports[ploidy] = {"records": len(records), "gt": sample_values["GT"],
                              "pl": len(sample_values["PL"].split(","))}

        poly_gvcf = work / "poly.g.vcf.gz"
        subprocess.run([
            str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(poly_gvcf), "--emit-ref-confidence", "GVCF",
            "--sample-ploidy", "3",
            "--min-depth", "1", "--min-alt-support", "1",
        ], check=True, stdout=subprocess.DEVNULL)
        gvcf_records = parse_records_gzip(poly_gvcf)
        assert gvcf_records
        block = next(row for row in gvcf_records if row[4] == "<NON_REF>")
        block_values = dict(zip(block[8].split(":"), block[9].split(":")))
        assert len(block_values["GT"].split("/")) == 3
        assert len(block_values["PL"].split(",")) == 4
        candidate = next(row for row in gvcf_records if row[4] != "<NON_REF>")
        candidate_values = dict(zip(candidate[8].split(":"), candidate[9].split(":")))
        assert len(candidate_values["GT"].split("/")) == 3
        assert len(candidate_values["PL"].split(",")) == 10

    print(json.dumps({"status": "pass", "ploidies": reports,
                      "gvcf_non_diploid": "reference-backed"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
