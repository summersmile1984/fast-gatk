#!/usr/bin/env python3
"""Pin the FilterMutectCalls contamination contribution to joint learning.

The Java implementation adds ContaminationFilter to the NON_SOMATIC error
type before ThresholdCalculator learns an OPTIMAL_F_SCORE threshold.  This
oracle uses a small deterministic AD/POPAF callset to ensure the native
pre-pass does the same.  The complete release-specific somatic model is not
claimed bit-identical here; the oracle instead checks the stable contamination
boundary and that the learned threshold is finite and in the pinned GATK
range rather than the pre-fix no-contamination sentinel (1.0).
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path


def read_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            records.append({"pos": fields[1], "filter": fields[6]})
    return records


def filtering_threshold(path: Path) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"#<METADATA>threshold=([0-9.eE+-]+)$", line)
        if match:
            return float(match.group(1))
    raise AssertionError(f"missing threshold metadata in {path}")


def contamination_positions(records: list[dict[str, str]]) -> list[str]:
    return [
        record["pos"] for record in records
        if "contamination" in record["filter"].split(";")
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
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK contamination joint oracle unavailable")
        print(json.dumps({
            "status": "skip",
            "reason": "bundled GATK contamination joint oracle unavailable",
        }))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=POPAF,Number=A,Type=Float,Description=Population allele frequency>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
    # The low-to-high AD series crosses the contamination/germline posterior
    # boundary while keeping the model inputs and ALT cardinality constant.
    alt_depths = (1, 2, 3, 5, 10, 15, 20, 25, 30, 40)
    body = "".join(
        f"17\t{69000 + index}\t.\tA\tG\t.\tPASS\tTLOD=8;POPAF=2\t"
        f"GT:AD\t0/1:{100 - alt_depth},{alt_depth}\n"
        for index, alt_depth in enumerate(alt_depths)
    )

    with tempfile.TemporaryDirectory(prefix="fastgatk-contamination-joint-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        input_path.write_text(header + body, encoding="utf-8")
        stats = work / "calls.stats"
        stats.write_text(
            "#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n",
            encoding="utf-8",
        )
        index_run = subprocess.run(
            [str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        assert index_run.returncode == 0, index_run.stderr

        gatk_output = work / "gatk.vcf.gz"
        gatk_stats = work / "gatk.filtering.table"
        gatk_run = subprocess.run(
            [
                str(java), "-jar", str(gatk), "FilterMutectCalls",
                "-R", str(reference), "-V", str(input_path), "-O", str(gatk_output),
                "--stats", str(stats), "--filtering-stats", str(gatk_stats),
                "--lenient", "true", "--contamination-estimate", "0.1",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        assert gatk_run.returncode == 0, gatk_run.stderr

        native_output = work / "native.vcf.gz"
        native_stats = work / "native.stats"
        native_run = subprocess.run(
            [
                str(binary), "-R", str(reference), "-V", str(input_path),
                "-O", str(native_output), "--stats", str(stats),
                "--filtering-stats", str(native_stats),
                "--contamination-estimate", "0.1",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        assert native_run.returncode == 0, native_run.stderr

        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert len(gatk_records) == len(native_records) == len(alt_depths)
        expected_contamination = [str(69000 + index) for index in range(7)]
        assert contamination_positions(gatk_records) == expected_contamination
        assert contamination_positions(native_records) == expected_contamination

        gatk_threshold = filtering_threshold(gatk_stats)
        native_json = json.loads(native_stats.read_text(encoding="utf-8"))
        native_threshold = float(native_json["joint_error_threshold"])
        assert abs(gatk_threshold - 0.193) < 1.0e-12
        # Native and Java use different release-specific somatic calibration
        # for the remaining filters, so allow that bounded model difference;
        # the old implementation produced exactly 1.0 by omitting this
        # posterior from the empirical pre-pass.
        assert 0.15 < native_threshold < 0.30
        assert abs(native_threshold - gatk_threshold) < 0.05
        assert native_threshold != 1.0
        # Java-parity path (FASTGATK_FMC_JAVA_PASSES=1): the SomaticClustering
        # Model learning trajectory is now GATK-machine-exact, so pin the
        # learned threshold against the EngineProbe/GATK value with a tight
        # band (R17: 5 EM(false) per split fix).  Guards against any future
        # off-by-N EM or cluster-learning regression.
        if os.environ.get("FASTGATK_FMC_JAVA_PASSES") == "1":
            assert abs(native_threshold - 0.1925917465406033) < 1.0e-8
        assert native_json["joint_error_threshold_learned"] is True
        assert native_json["joint_error_threshold_observations"] == len(alt_depths)
        assert native_json["contamination_posterior_records"] == len(alt_depths)
        assert native_json["contamination_posterior_alleles"] == len(alt_depths)

        print(json.dumps({
            "status": "pass",
            "gatk_threshold": gatk_threshold,
            "native_threshold": native_threshold,
            "threshold_delta": native_threshold - gatk_threshold,
            "gatk_contamination_positions": expected_contamination,
            "native_contamination_positions": expected_contamination,
            "contamination_joint_learning": True,
            "gatk_release": "4.6.2.0",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
