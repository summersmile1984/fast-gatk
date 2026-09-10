#!/usr/bin/env python3
"""Verify HaplotypeCaller BP_RESOLUTION records against pinned GATK."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_records(path: Path) -> tuple[list[str], list[tuple[str, ...]]]:
    samples: list[str] = []
    rows: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            values = dict(zip(fields[8].split(":"), fields[9].split(":")))
            rows.append((fields[0], fields[1], fields[3], fields[4],
                         *(values.get(key, "") for key in ("GT", "AD", "DP", "GQ", "PL"))))
    return samples, rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", str(root / "fastgatk-native/build/fastgatk-hc-call")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, native, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HaplotypeCaller BP_RESOLUTION oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-bp-oracle-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.bp.g.vcf.gz"
        native_output = work / "native.bp.g.vcf.gz"
        region = "17:69000-69020"
        gatk_run = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "BP_RESOLUTION", "-O", str(gatk_output),
            "--create-output-variant-index", "false", "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "BP_RESOLUTION", "-O", str(native_output),
            "--min-depth", "1", "--min-alt-support", "1",
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        gatk_samples, gatk_records = read_records(gatk_output)
        native_samples, native_records = read_records(native_output)
        assert gatk_samples == native_samples == ["NA12878"]
        assert gatk_records == native_records, {
            "gatk_records": len(gatk_records), "native_records": len(native_records),
            "first_gatk": gatk_records[:3], "first_native": native_records[:3],
        }
        positions = [int(row[1]) for row in native_records]
        assert positions == list(range(69000, 69021))
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "mode": "BP_RESOLUTION",
            "records": len(native_records),
            "sample_contract_exact": True,
            "gt_ad_dp_gq_pl_exact": True,
            "positions_contiguous": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
