#!/usr/bin/env python3
"""Pinned GATK oracle for issue3845's chrM-only reference dictionary.

The upstream fixture deliberately supplies a reference whose dictionary is a
subset of the input BAM dictionary.  GATK accepts that input, applies its
default 50-read alignment-start reservoir, and annotates the call with the
candidate-local retained AlleleLikelihoods evidence.  This protects all three
Host-side contracts together: dictionary subset validation, Java-compatible
downsampling, and post-realignment MQ/DP evidence ownership.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout={completed.stdout[-2000:]}\nstderr={completed.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    fixture = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/"
        "tools/haplotypecaller/issue3845_revertSoftClip_bug")
    reference = root / "gatk-source/src/test/resources/Homo_sapiens_assembly38_chrM_only.fasta"
    bam = fixture / "issue3845_bug.bam"
    if not all(path.is_file() for path in (native, java, gatk, reference, bam)):
        raise SystemExit("missing HC issue3845 oracle assets")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-issue3845-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        gatk_gvcf = work / "gatk.g.vcf"
        native_gvcf = work / "native.g.vcf"
        gatk_bp = work / "gatk.bp.vcf"
        native_bp = work / "native.bp.vcf"
        gatk_cap100_vcf = work / "gatk.cap100.vcf"
        native_cap100_vcf = work / "native.cap100.vcf"
        cap100_manifest = work / "native.cap100.json"
        manifest = work / "native.json"
        gatk_shared = [
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", "chrM",
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ]
        native_shared = [
            str(native), "-R", str(reference), "-I", str(bam), "-L", "chrM",
            "--threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ]
        run([*gatk_shared, "-O", str(gatk_vcf)])
        run([*native_shared, "-O", str(native_vcf), "--output-manifest", str(manifest)])
        if gatk_vcf.read_bytes() != native_vcf.read_bytes():
            raise AssertionError("native VCF differs from GATK for issue3845 fixture")
        run([*gatk_shared, "-ERC", "GVCF", "-O", str(gatk_gvcf)])
        run([*native_shared, "-ERC", "GVCF", "-O", str(native_gvcf)])
        if gatk_gvcf.read_bytes() != native_gvcf.read_bytes():
            raise AssertionError("native GVCF differs from GATK for issue3845 fixture")
        run([*gatk_shared, "-ERC", "BP_RESOLUTION", "-O", str(gatk_bp)])
        run([*native_shared, "-ERC", "BP_RESOLUTION", "-O", str(native_bp)])
        if gatk_bp.read_bytes() != native_bp.read_bytes():
            raise AssertionError("native BP_RESOLUTION output differs from GATK for issue3845 fixture")
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        if telemetry.get("read_filter_max_reads_per_locus") != 50:
            raise AssertionError("native did not apply GATK's default 50-read start cap")
        if telemetry.get("downsampled_reads") != 651:
            raise AssertionError("native reservoir result differs for issue3845 fixture")
        # GenotypeLikelihoods.getAsPLs does not cap ordinary VCF PL values
        # at 999.  With 100 retained reads, this fixture emits PL=1356 and
        # QUAL=1342.06; compare the complete VCF so the C++ Host's PL/QUAL
        # materialization cannot silently reintroduce a legacy 999 cap.
        cap100_args = ["--max-reads-per-alignment-start", "100"]
        run([*gatk_shared, *cap100_args, "-O", str(gatk_cap100_vcf)])
        run([*native_shared, *cap100_args, "-O", str(native_cap100_vcf),
             "--output-manifest", str(cap100_manifest)])
        if gatk_cap100_vcf.read_bytes() != native_cap100_vcf.read_bytes():
            raise AssertionError("native high-depth HC VCF differs from GATK for issue3845 fixture")
        cap100_telemetry = json.loads(cap100_manifest.read_text(encoding="utf-8"))["telemetry"]
        if cap100_telemetry.get("downsampled_reads") != 601:
            raise AssertionError("native 100-read reservoir result differs for issue3845 fixture")
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "fixture": "HaplotypeCaller issue3845",
            "vcf_byte_identical": True,
            "gvcf_byte_identical": True,
            "bp_resolution_byte_identical": True,
            "reference_dictionary": "BAM-dictionary subset accepted",
            "default_max_reads_per_alignment_start": 50,
            "downsampled_reads": 651,
            "high_depth_pl_qual": "100-read cap VCF byte-identical (PL > 999)",
            "scope": "Host dictionary, reservoir, and retained MQ evidence",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
