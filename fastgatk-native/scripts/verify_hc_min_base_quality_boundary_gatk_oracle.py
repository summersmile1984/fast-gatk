#!/usr/bin/env python3
"""Pinned GATK gVCF oracle for the Q0/Q1 HC terminal-quality boundary.

HaplotypeCaller converts ``--min-base-quality-score`` to the signed
``minTailQuality = score - 1`` before it finalizes a region.  It then builds
FragmentCollection from those finalized reads, so a reverted soft clip may
change both a read's end and its coordinate-sort position before paired-overlap
quality correction.  This fixture exercises the Q0/Q1 boundary and compares
every decompressed gVCF record, including reference-block DP medians.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def normalized_vcf(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith("##GATKCommandLine=")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    if not native.exists() or not reference.exists() or not bam.exists():
        raise SystemExit("missing native build or HC min-base-quality boundary fixture")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    region = "17:69000-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-min-base-quality-boundary-") as directory:
        work = Path(directory)
        for quality in (0, 1):
            gatk_gvcf = work / f"gatk-q{quality}.g.vcf.gz"
            native_gvcf = work / f"native-q{quality}.g.vcf.gz"
            manifest = work / f"native-q{quality}.manifest.json"
            common = [
                "-R", str(reference), "-I", str(bam), "-L", region,
                "--min-base-quality-score", str(quality), "-ERC", "GVCF",
                "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false",
            ]
            subprocess.run([
                java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
                *common, "-O", str(gatk_gvcf), "--native-pair-hmm-threads", "2",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run([
                str(native), *common, "-O", str(native_gvcf), "--threads", "2",
                "--output-manifest", str(manifest),
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            assert normalized_vcf(gatk_gvcf) == normalized_vcf(native_gvcf), (
                f"Q{quality} gVCF differs from GATK at the terminal-quality boundary")
            data = json.loads(manifest.read_text(encoding="utf-8"))
            telemetry = data["telemetry"]
            assert telemetry["min_base_quality"] == quality
            assert telemetry["overlapping_quality_correction_used"] is True
            assert telemetry["overlapping_pairs"] > 0
            assert data["compatibility"]["reference_confidence_kokkos"] is True

        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0", "region": region,
            "min_base_qualities": [0, 1], "gvcf_exact": True,
            "host_finalization_fragment_collection": True,
            "kokkos_reference_confidence": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
