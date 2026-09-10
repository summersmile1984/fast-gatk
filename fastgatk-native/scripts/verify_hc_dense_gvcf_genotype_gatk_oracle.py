#!/usr/bin/env python3
"""Strict HC gVCF -> GenotypeGVCFs oracle on a dense phased chr20 region.

The ordinary GenotypeGVCFs fixtures exercise cohort math and multi-sample
merging, but their source gVCFs do not contain a compact run of phased calls
interleaved with reference-confidence blocks.  This fixture does: it validates
the complete default C++ Host handoff from EventMap/CIGAR/RCM production to
the Kokkos genotype kernels and back to GATK's phased VCF writer semantics.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


REGION = "20:10019901-10020710"
EXPECTED_POSITIONS = [
    10019967, 10019969, 10020228, 10020229, 10020429, 10020431,
    10020434, 10020438, 10020679, 10020680, 10020681,
]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise AssertionError(f"{label} failed: {result.stderr[-2000:]}")


def decompressed_bytes(path: Path) -> bytes:
    with gzip.open(path, "rb") as stream:
        return stream.read()


def normalized_joint_text(path: Path) -> list[str]:
    # Command provenance, source/date and contig dictionary order are writer
    # metadata rather than GenotypeGVCFs call semantics.  All remaining header
    # lines and every data row must be byte-for-byte equal.
    ignored = (
        "##GATKCommandLine=", "##source=", "##fileDate=", "##contig=",
        "##fastgatk_genotype_gvcfs_status=",
    )
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if not line.startswith(ignored)]


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def format_sample(record: list[str]) -> dict[str, str]:
    return dict(zip(record[8].split(":"), record[9].split(":")))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    hc = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    genotype = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        root / "fastgatk-native/build/fastgatk-genotype-gvcf"))
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "fixtures/chr20/ref20mnp.fasta"
    bam = root / "fixtures/chr20/mnp.bam"
    required = (hc, genotype, java, gatk, reference, bam)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("dense HC -> GenotypeGVCFs oracle inputs are required")
        print(json.dumps({"status": "skipped", "reason": "oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-dense-gvcf-genotype-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf.gz"
        native_gvcf = work / "native.g.vcf.gz"
        native_hc_manifest = work / "native.hc.json"
        gatk_joint = work / "gatk.joint.vcf"
        native_joint = work / "native.joint.vcf"
        renamed_gvcf = work / "sample2.g.vcf.gz"
        cohort_gvcf = work / "cohort.g.vcf.gz"
        gatk_cohort_joint = work / "gatk.cohort.joint.vcf"
        native_cohort_joint = work / "native.cohort.joint.vcf"

        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", REGION,
            "-ERC", "GVCF", "-O", str(gatk_gvcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "true",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK HaplotypeCaller")
        run([
            str(hc), "-R", str(reference), "-I", str(bam), "-L", REGION,
            "-ERC", "GVCF", "-O", str(native_gvcf),
            "--output-manifest", str(native_hc_manifest), "--threads", "1",
            "--create-output-variant-index", "true",
            "--add-output-vcf-command-line", "false",
        ], "native HaplotypeCaller")
        assert decompressed_bytes(gatk_gvcf) == decompressed_bytes(native_gvcf), (
            "dense default HC gVCF differs from GATK")
        hc_manifest = json.loads(native_hc_manifest.read_text(encoding="utf-8"))
        assert hc_manifest["telemetry"]["rcm_haplotype_realignment_used"] is True
        assert hc_manifest["telemetry"]["rcm_realigned_observations"] > 0

        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-L", REGION,
            "-O", str(gatk_joint), "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK GenotypeGVCFs")
        run([
            str(genotype), "-R", str(reference), "-V", str(native_gvcf),
            "-L", REGION, "-O", str(native_joint),
            "--create-output-variant-index", "false",
            "--gatk-compatible-annotations",
        ], "native GenotypeGVCFs")
        assert normalized_joint_text(gatk_joint) == normalized_joint_text(native_joint), (
            "dense GenotypeGVCFs header or records differ from GATK")
        joint_rows = rows(native_joint)
        assert [int(row[1]) for row in joint_rows] == EXPECTED_POSITIONS
        by_position = {int(row[1]): format_sample(row) for row in joint_rows}
        assert by_position[10019967]["GT"] == "0|1"
        assert by_position[10019967]["PGT"] == "0|1"
        assert by_position[10019967]["PID"] == "10019967_C_G"
        assert by_position[10020680]["GT"] == "1|1"
        assert by_position[10020680]["PGT"] == "1|1"
        assert list(format_sample(joint_rows[0])) == [
            "GT", "AD", "DP", "GQ", "PGT", "PID", "PL", "PS"]

        # Exercise Host FORMAT materialization for more than one sample.
        # This is deliberately a GATK-generated combined input: the assertion
        # isolates the native GenotypeGVCFs sample-union boundary from
        # CombineGVCFs, which is outside this core caller oracle.
        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "RenameSampleInVcf",
            "--INPUT", str(gatk_gvcf), "--OUTPUT", str(renamed_gvcf),
            "--NEW_SAMPLE_NAME", "SAMPLE2", "--CREATE_INDEX", "true",
            "--VALIDATION_STRINGENCY", "STRICT",
        ], "GATK RenameSampleInVcf")
        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "CombineGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-V", str(renamed_gvcf),
            "-O", str(cohort_gvcf), "--create-output-variant-index", "true",
            "--seconds-between-progress-updates", "1",
        ], "GATK CombineGVCFs")
        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(cohort_gvcf), "-L", REGION,
            "-O", str(gatk_cohort_joint), "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK cohort GenotypeGVCFs")
        run([
            str(genotype), "-R", str(reference), "-V", str(cohort_gvcf),
            "-L", REGION, "-O", str(native_cohort_joint),
            "--create-output-variant-index", "false",
            "--gatk-compatible-annotations",
        ], "native cohort GenotypeGVCFs")
        assert normalized_joint_text(gatk_cohort_joint) == normalized_joint_text(
            native_cohort_joint), "dense two-sample GenotypeGVCFs differs from GATK"
        cohort_rows = rows(native_cohort_joint)
        assert len(cohort_rows) == len(EXPECTED_POSITIONS)
        cohort_first = format_sample(cohort_rows[0])
        cohort_second = dict(zip(cohort_rows[0][8].split(":"),
                                 cohort_rows[0][10].split(":")))
        assert cohort_first["GT"] == cohort_second["GT"] == "0|1"
        assert cohort_first["PGT"] == cohort_second["PGT"] == "0|1"
        cohort_hom_var = next(row for row in cohort_rows if int(row[1]) == 10020680)
        cohort_hom_var_samples = [
            dict(zip(cohort_hom_var[8].split(":"), sample.split(":")))
            for sample in cohort_hom_var[9:]
        ]
        assert all(sample["GT"] == sample["PGT"] == "1|1"
                   for sample in cohort_hom_var_samples)
        print(json.dumps({
            "status": "pass",
            "region": REGION,
            "gvcf_decompressed_text_exact": True,
            "joint_normalized_full_text_exact": True,
            "joint_records": len(joint_rows),
            "cohort_joint_records": len(cohort_rows),
            "phased_gt_preserved": True,
            "hom_var_pgt_normalized": True,
            "multisample_phasing_exact": True,
            "rcm_haplotype_realignment_default": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
