#!/usr/bin/env python3
"""Pinned GATK oracle for Mutect2's overloaded-pool downsampling policy.

The issue3845 BAM has 701 post-filter reads beginning at chrM:1.  HC uses
AssemblyRegionWalker's ordinary reservoir, but Mutect2 substitutes
``MutectDownsampler``: an overloaded stride samples only reads with MAPQ > 50
and an optional suspicious-read limit drops its whole stride.  This fixture
exposes each Host-side distinction through the complete somatic VCF at chrM:73.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{completed.stderr}"
        )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    native = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")))
    fixture = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/"
        "tools/haplotypecaller/issue3845_revertSoftClip_bug")
    reference = root / "gatk-source/src/test/resources/Homo_sapiens_assembly38_chrM_only.fasta"
    bam = fixture / "issue3845_bug.bam"
    required = (java, gatk, native, reference, bam)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_mutect2_issue3845_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled Mutect2 issue3845 oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-issue3845-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        gatk_stride_vcf = work / "gatk-stride.vcf"
        native_stride_vcf = work / "native-stride.vcf"
        stride_manifest = work / "native-stride.manifest.json"
        gatk_suspicious_vcf = work / "gatk-suspicious.vcf"
        native_suspicious_vcf = work / "native-suspicious.vcf"
        manifest = work / "native.manifest.json"
        shared = ["-R", str(reference), "-I", str(bam), "--tumor-sample", "FOO",
                  "-L", "chrM", "--native-pair-hmm-threads", "1",
                  "--add-output-vcf-command-line", "false"]
        run([str(java), "-jar", str(gatk), "Mutect2", *shared, "-O", str(gatk_vcf)])
        run([str(native), *shared, "-O", str(native_vcf),
             "--output-manifest", str(manifest)])
        if gatk_vcf.read_bytes() != native_vcf.read_bytes():
            raise AssertionError("native Mutect2 VCF differs from GATK for the issue3845 fixture")
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        if telemetry.get("read_filter_max_reads_per_locus") != 50:
            raise AssertionError("native did not apply GATK's default 50-read start cap")
        if telemetry.get("downsampled_reads") != 651:
            raise AssertionError("native MutectDownsampler selected a different read set")
        # Pools begin at their first observed alignment start, and their
        # capacity is max-reads-per-start * stride.  This 701-read fixture
        # therefore retains 100 reads with stride two.  Beyond the C++ kernel
        # smoke's per-read reservoir sequence, compare the complete GATK VCF
        # to cover the Host's trimmed EventMap/PairHMM boundary as well.
        stride_args = [*shared, "-stride", "2"]
        run([str(java), "-jar", str(gatk), "Mutect2", *stride_args,
             "-O", str(gatk_stride_vcf)])
        run([str(native), *stride_args, "-O", str(native_stride_vcf),
             "--output-manifest", str(stride_manifest)])
        if gatk_stride_vcf.read_bytes() != native_stride_vcf.read_bytes():
            raise AssertionError("native Mutect2 stride VCF differs from GATK for issue3845")
        stride_telemetry = json.loads(stride_manifest.read_text(encoding="utf-8"))["telemetry"]
        if stride_telemetry.get("downsampled_reads") != 601:
            raise AssertionError("native Mutect2 stride did not retain the GATK-sized 100-read pool")
        # A single MAPQ <= 50 read rejects the entire pending stride.  All
        # retained reads in this fixture share the same start, so the emitted
        # data section must be empty rather than merely lower-depth.
        suspicious_args = [*shared,
                           "--max-suspicious-reads-per-alignment-start", "1"]
        run([str(java), "-jar", str(gatk), "Mutect2", *suspicious_args,
             "-O", str(gatk_suspicious_vcf)])
        run([str(native), *suspicious_args, "-O", str(native_suspicious_vcf)])
        if gatk_suspicious_vcf.read_bytes() != native_suspicious_vcf.read_bytes():
            raise AssertionError("native Mutect2 suspicious-read downsampling differs from GATK")
        if any(not line.startswith("#") for line in native_suspicious_vcf.read_text(
                encoding="utf-8").splitlines()):
            raise AssertionError("suspicious-read limit should reject the issue3845 stride")
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "fixture": "Mutect2 issue3845",
            "vcf_byte_identical": True,
            "downsampling": "MutectDownsampler MAPQ>50 reservoir, stride, suspicious limit",
            "downsampled_reads": telemetry["downsampled_reads"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
