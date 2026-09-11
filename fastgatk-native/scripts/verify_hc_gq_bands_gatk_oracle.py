#!/usr/bin/env python3
"""Verify HaplotypeCaller custom GVCF GQ partitions against GATK."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_gvcf(path: Path) -> tuple[list[str], list[str]]:
    blocks: list[str] = []
    records: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("##GVCFBlock"):
                blocks.append(line.rstrip("\n"))
            elif line and not line.startswith("#"):
                records.append(line.rstrip("\n"))
    return blocks, records


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
        oracle_guard.oracle_not_verified('verify_hc_gq_bands_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HaplotypeCaller GQ-band oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    bands = ("10", "20", "50")
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-gq-bands-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.g.vcf.gz"
        native_output = work / "native.g.vcf.gz"
        native_streamed_output = work / "native.streamed.g.vcf.gz"
        native_streamed_manifest = work / "native.streamed.manifest.json"
        common = [
            "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
            "--emit-ref-confidence", "GVCF", "-O", str(gatk_output),
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ]
        gatk_run = subprocess.run(
            [str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
             *sum((["-GQB", band] for band in bands), [])],
            text=True, capture_output=True, check=False,
        )
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run(
            [str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
             "-ERC", "GVCF", *sum((["-GQB", band] for band in bands), []),
             "-O", str(native_output), "--create-output-variant-index", "false",
             "--min-depth", "1", "--min-alt-support", "1"],
            text=True, capture_output=True, check=False,
        )
        assert native_run.returncode == 0, native_run.stderr
        streamed_run = subprocess.run(
            [str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
             "-ERC", "GVCF", *sum((["-GQB", band] for band in bands), []),
             "-O", str(native_streamed_output), "--stream-by-region", "8",
             "--create-output-variant-index", "false", "--min-depth", "1",
             "--min-alt-support", "1", "--output-manifest", str(native_streamed_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert streamed_run.returncode == 0, streamed_run.stderr
        invalid_run = subprocess.run(
            [str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
             "-ERC", "GVCF", "-GQB", "20", "-GQB", "10", "-O", str(work / "invalid.g.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert invalid_run.returncode != 0 and "strictly increasing" in invalid_run.stderr
        gatk_blocks, gatk_records = read_gvcf(gatk_output)
        native_blocks, native_records = read_gvcf(native_output)
        streamed_blocks, streamed_records = read_gvcf(native_streamed_output)
        expected = [
            "##GVCFBlock0-10=minGQ=0(inclusive),maxGQ=10(exclusive)",
            "##GVCFBlock10-20=minGQ=10(inclusive),maxGQ=20(exclusive)",
            "##GVCFBlock20-50=minGQ=20(inclusive),maxGQ=50(exclusive)",
            "##GVCFBlock50-100=minGQ=50(inclusive),maxGQ=100(exclusive)",
        ]
        assert gatk_blocks == native_blocks == sorted(expected)
        assert gatk_records == native_records, {
            "gatk_records": gatk_records, "native_records": native_records,
        }
        assert gatk_blocks == streamed_blocks == sorted(expected)
        assert gatk_records == streamed_records, {
            "gatk_records": gatk_records, "streamed_records": streamed_records,
        }
        streamed_telemetry = json.loads(
            native_streamed_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert streamed_telemetry["stream_by_region"] is True
        assert streamed_telemetry["streamed_regions"] == 3
        # GATK's writer-only --floor-blocks mode retains all calling/RCM work
        # but serializes ordinary reference blocks without MIN_DP or PL.  Keep
        # the pinned Java comparison here because the GQ-band fixture already
        # exercises real blocks and avoids turning a parser-only check into a
        # compatibility claim.
        gatk_floor_output = work / "gatk.floor.g.vcf.gz"
        native_floor_output = work / "native.floor.g.vcf.gz"
        native_floor_manifest = work / "native.floor.manifest.json"
        floor_common = [
            "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
            "--emit-ref-confidence", "GVCF", "-O", str(gatk_floor_output),
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ]
        floor_gatk_run = subprocess.run(
            [str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *floor_common,
             *sum((["-GQB", band] for band in bands), []), "--floor-blocks", "true"],
            text=True, capture_output=True, check=False,
        )
        assert floor_gatk_run.returncode == 0, floor_gatk_run.stderr
        floor_native_run = subprocess.run(
            [str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-69020",
             "-ERC", "GVCF", *sum((["-GQB", band] for band in bands), []),
             "--floor-blocks", "true", "-O", str(native_floor_output),
             "--create-output-variant-index", "false", "--min-depth", "1",
             "--min-alt-support", "1", "--output-manifest", str(native_floor_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert floor_native_run.returncode == 0, floor_native_run.stderr
        floor_gatk_blocks, floor_gatk_records = read_gvcf(gatk_floor_output)
        floor_native_blocks, floor_native_records = read_gvcf(native_floor_output)
        assert floor_gatk_blocks == floor_native_blocks == sorted(expected)
        assert floor_gatk_records == floor_native_records, {
            "gatk_records": floor_gatk_records, "native_records": floor_native_records,
        }
        floor_ref_blocks = [record.split("\t") for record in floor_native_records
                            if record.split("\t")[4] == "<NON_REF>"]
        assert floor_ref_blocks
        assert all(row[8] == "GT:DP:GQ" and row[9].count(":") == 2
                   for row in floor_ref_blocks), floor_ref_blocks
        floor_manifest = json.loads(native_floor_manifest.read_text(encoding="utf-8"))
        compatibility = floor_manifest["compatibility"]
        assert compatibility["floor_blocks"] is True
        assert compatibility["gvcf_floor_blocks"] is True
        print(json.dumps({
            "status": "pass", "gatk_version": "4.6.2.0",
            "gq_bands": [10, 20, 50, 100], "gvcf_blocks": len(native_blocks),
            "records": len(native_records), "header_exact": True, "records_exact": True,
            "streamed_records_exact": True, "streamed_regions": 3,
            "floor_blocks_records_exact": True,
            "floor_blocks_reference_format_exact": True,
            "floor_blocks_manifest_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
