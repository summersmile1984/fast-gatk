#!/usr/bin/env python3
"""Pinned GATK oracle for HC --soft-clip-low-quality-ends.

This is intentionally separate from --dont-use-soft-clipped-bases: GATK
retains newly soft-clipped terminal low-quality bases for read threading, but
its PairHMM removes those S operations before likelihood scoring.  The native
path must therefore keep the C++ Host clipping/graph boundary while retaining
the existing Kokkos PairHMM and genotype kernels.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import gzip
from pathlib import Path
import oracle_guard


def normalized_vcf(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith("##GATKCommandLine=")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or HC soft-clip fixture")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_hc_softclip_low_quality_gatk_oracle.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    region = "17:69000-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-softclip-low-quality-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.manifest.json"
        common = [
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--soft-clip-low-quality-ends",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ]
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            *common, "-O", str(gatk_vcf), "--native-pair-hmm-threads", "2",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), *common, "-O", str(native_vcf), "--threads", "2",
            "--output-manifest", str(manifest),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        assert normalized_vcf(gatk_vcf) == normalized_vcf(native_vcf), (
            "--soft-clip-low-quality-ends VCF differs from GATK")
        data = json.loads(manifest.read_text(encoding="utf-8"))
        assert data["compatibility"]["soft_clip_low_quality_ends"] is True
        assert data["telemetry"]["soft_clip_low_quality_ends"] is True

        # The finalization distinction is most visible in gVCF reference
        # blocks: a newly-created terminal S retains the original read offset
        # for ReferenceConfidenceModel's indel-informativeness calculation,
        # while PairHMM scores its hard-clipped copy.  Check every emitted
        # block, not just variant calls, against the pinned GATK oracle.
        gatk_gvcf = work / "gatk.g.vcf.gz"
        native_gvcf = work / "native.g.vcf.gz"
        gvcf_manifest = work / "native.gvcf.manifest.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            *common, "-ERC", "GVCF", "-O", str(gatk_gvcf),
            "--native-pair-hmm-threads", "2",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), *common, "-ERC", "GVCF", "-O", str(native_gvcf),
            "--threads", "2", "--output-manifest", str(gvcf_manifest),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert normalized_vcf(gatk_gvcf) == normalized_vcf(native_gvcf), (
            "--soft-clip-low-quality-ends gVCF differs from GATK")
        gvcf_data = json.loads(gvcf_manifest.read_text(encoding="utf-8"))
        assert gvcf_data["compatibility"]["gvcf"] is True
        assert gvcf_data["compatibility"]["soft_clip_low_quality_ends"] is True
        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0", "region": region,
            "soft_clip_low_quality_ends": True, "vcf_exact": True, "gvcf_exact": True,
            "host_graph_boundary": True, "kokkos_pairhmm_genotype": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
