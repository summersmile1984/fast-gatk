#!/usr/bin/env python3
"""Pin ReadOrientationFilter's contribution to learned joint posteriors."""

from __future__ import annotations

import gzip
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_filter_mutect_orientation_gatk_oracle import prior_table  # noqa: E402


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def read_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                records.append({"pos": fields[1], "filter": fields[6]})
    return records


def orientation_positions(records: list[dict[str, str]]) -> list[str]:
    return [
        record["pos"] for record in records
        if "orientation" in record["filter"].split(";")
    ]


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
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK orientation joint oracle unavailable")
        print(json.dumps({
            "status": "skip",
            "reason": "bundled GATK orientation joint oracle unavailable",
        }))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=F1R2 orientation counts>
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=F2R1 orientation counts>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
    body = (
        "17\t69000\t.\tT\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,10:20,0\n"
        "17\t69001\t.\tT\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,5:20,5\n"
        "17\t69002\t.\tTT\tCG\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,10:20,0\n"
    )

    with tempfile.TemporaryDirectory(prefix="fastgatk-orientation-joint-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        input_path.write_text(header + body, encoding="utf-8")
        stats = work / "stats.table"
        stats.write_text(
            "#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n",
            encoding="utf-8",
        )
        prior = work / "TUMOR.orientation_priors.tar.gz"
        with tarfile.open(prior, "w:gz") as archive:
            payload = prior_table().encode("utf-8")
            member = tarfile.TarInfo("TUMOR.orientation_priors")
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
        index_run = run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)])
        assert index_run.returncode == 0, index_run.stderr

        common = ["--orientation-bias-artifact-priors", str(prior)]
        gatk_output = work / "gatk.vcf.gz"
        gatk_stats = work / "gatk.filtering.table"
        gatk_run = run([
            str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
            "-V", str(input_path), "-O", str(gatk_output), "--stats", str(stats),
            "--filtering-stats", str(gatk_stats), "--lenient", "true", *common,
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr

        native_output = work / "native.vcf.gz"
        native_stats = work / "native.stats.json"
        native_run = run([
            str(binary), "-R", str(reference), "-V", str(input_path),
            "-O", str(native_output), "--stats", str(stats),
            "--filtering-stats", str(native_stats), *common,
        ])
        assert native_run.returncode == 0, native_run.stderr

        baseline_output = work / "baseline.vcf.gz"
        baseline_stats = work / "baseline.stats.json"
        baseline_run = run([
            str(binary), "-R", str(reference), "-V", str(input_path),
            "-O", str(baseline_output), "--stats", str(stats),
            "--filtering-stats", str(baseline_stats),
        ])
        assert baseline_run.returncode == 0, baseline_run.stderr

        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        baseline_records = read_records(baseline_output)
        assert len(gatk_records) == len(native_records) == len(baseline_records) == 3
        expected_orientation = ["69002"]
        assert orientation_positions(gatk_records) == expected_orientation
        assert orientation_positions(native_records) == expected_orientation
        assert orientation_positions(baseline_records) == []

        def threshold_from_table(path: Path) -> float:
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.startswith("#<METADATA>threshold="):
                    return float(line.split("=", 1)[1])
            raise AssertionError(f"missing GATK threshold in {path}")

        gatk_threshold = threshold_from_table(gatk_stats)
        native_json = json.loads(native_stats.read_text(encoding="utf-8"))
        baseline_json = json.loads(baseline_stats.read_text(encoding="utf-8"))
        assert gatk_threshold == 0.0
        assert native_json["orientation_prior_records"] == 3
        assert native_json["orientation_prior_sample_evaluations"] == 3
        assert native_json["joint_error_threshold_observations"] == 3
        assert native_json["joint_error_threshold"] < 1.0e-10
        # Without the orientation contributor the same native callset has no
        # learned joint threshold (all three records pass at threshold 1.0).
        assert baseline_json["joint_error_threshold"] == 1.0

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "gatk_threshold": gatk_threshold,
            "native_threshold": native_json["joint_error_threshold"],
            "baseline_without_orientation_threshold": baseline_json["joint_error_threshold"],
            "orientation_filter_pattern_exact": True,
            "orientation_joint_learning": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
