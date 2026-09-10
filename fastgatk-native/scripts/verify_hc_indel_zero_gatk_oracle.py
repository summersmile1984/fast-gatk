#!/usr/bin/env python3
"""Pinned gVCF oracle for HC --indel-size-to-eliminate-in-ref-model 0.

GATK's zero value does not remove reference-confidence evidence.  Its indel
checker enumerates no candidate indel sizes and consequently regards every
ordinary (non insertion/deletion-anchor) pileup element as informative.  This
test locks that public CLI behavior to a full decompressed gVCF comparison.
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
        raise SystemExit("missing native build or HC indel-zero fixture")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    region = "17:69000-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-indel-zero-") as directory:
        work = Path(directory)
        common = [
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--indel-size-to-eliminate-in-ref-model", "0",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ]
        cases = (
            ("default", []),
            ("dont-use-soft-clipped-bases", ["--dont-use-soft-clipped-bases"]),
        )
        for label, extra in cases:
            gatk_gvcf = work / f"gatk.{label}.g.vcf.gz"
            native_gvcf = work / f"native.{label}.g.vcf.gz"
            manifest = work / f"native.{label}.manifest.json"
            subprocess.run([
                java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
                *common, *extra, "-O", str(gatk_gvcf),
                "--native-pair-hmm-threads", "2",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run([
                str(native), *common, *extra, "-O", str(native_gvcf),
                "--threads", "2", "--output-manifest", str(manifest),
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            assert normalized_vcf(gatk_gvcf) == normalized_vcf(native_gvcf), (
                f"--indel-size-to-eliminate-in-ref-model 0 gVCF differs from GATK for {label}")
            data = json.loads(manifest.read_text(encoding="utf-8"))
            assert data["compatibility"]["reference_confidence_kokkos"] is True
            assert data["telemetry"]["indel_size_to_eliminate_in_ref_model"] == 0
            assert data["telemetry"]["gvcf_reference_blocks"] > 4
            assert data["telemetry"]["use_soft_clipped_bases"] is (label == "default")
        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0", "region": region,
            "indel_size_to_eliminate_in_ref_model": 0,
            "cases": [case[0] for case in cases],
            "gvcf_exact": True,
            "host_cigar_informativeness": True,
            "kokkos_reference_confidence": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
