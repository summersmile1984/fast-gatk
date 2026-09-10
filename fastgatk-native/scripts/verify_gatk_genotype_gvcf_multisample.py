#!/usr/bin/env python3
"""Strict two-sample GenotypeGVCFs oracle.

The ordinary GenotypeGVCFs oracle exercises a single-sample gVCF.  This test
builds a real two-sample cohort with GATK's RenameSampleInVcf and
CombineGVCFs, then requires the native joint output to be text-identical under
the GATK-compatible annotation profile.  It intentionally keeps the fixture
small and deterministic while exercising cohort AF/QUAL, FORMAT merging, and
sample-name/order preservation.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return stream.read().splitlines()


def fields(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in read_lines(path)
            if line and not line.startswith("#")]


def semantic_header(path: Path) -> list[str]:
    ignored = ("##GATKCommandLine=", "##source=", "##fileDate=",
               "##contig=", "##GVCFBlock", "##fastgatk_genotype_gvcfs_status=")
    return sorted({line for line in read_lines(path)
                   if line.startswith("##") and not line.startswith(ignored)})


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or GATK fixtures; run build_native.sh first")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    region = "17:69000-69100"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-multisample-oracle-") as directory:
        work = Path(directory)
        first = work / "sample1.g.vcf.gz"
        second = work / "sample2.g.vcf.gz"
        cohort = work / "cohort.g.vcf.gz"
        gatk_joint = work / "gatk.joint.vcf.gz"
        native_joint = work / "native.joint.vcf.gz"
        native_manifest = work / "native.joint.manifest.json"

        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-O", str(first), "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "true",
            "--seconds-between-progress-updates", "1",
        ], "GATK HaplotypeCaller")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "RenameSampleInVcf",
            "--INPUT", str(first), "--OUTPUT", str(second),
            "--NEW_SAMPLE_NAME", "SAMPLE2", "--CREATE_INDEX", "true",
            "--VALIDATION_STRINGENCY", "STRICT",
        ], "GATK RenameSampleInVcf")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "CombineGVCFs",
            "-R", str(reference), "-V", str(first), "-V", str(second),
            "-O", str(cohort), "--create-output-variant-index", "true",
            "--seconds-between-progress-updates", "1",
        ], "GATK CombineGVCFs")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(cohort), "-L", region,
            "-O", str(gatk_joint), "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK GenotypeGVCFs")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(cohort), "-L", region,
            "-O", str(native_joint), "--create-output-variant-index", "false",
            "--gatk-compatible-annotations", "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True)
        if native_result.returncode != 0:
            raise RuntimeError(f"native GenotypeGVCFs failed:\n{native_result.stderr[-4000:]}")
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))

        gatk_rows = fields(gatk_joint)
        native_rows = fields(native_joint)
        gatk_header = semantic_header(gatk_joint)
        native_header = semantic_header(native_joint)

    assert gatk_rows, "GATK produced no variant rows for the cohort fixture"
    assert native_rows == gatk_rows, "native multi-sample data rows are not text-identical"
    assert gatk_header == native_header, "native multi-sample semantic header differs"
    assert len(gatk_rows[0]) == 11, "expected two sample columns in the joint output"
    assert gatk_rows[0][9].startswith("1/1:") and gatk_rows[0][10].startswith("1/1:")
    assert native_summary["status"] == "contract-compatible"
    assert manifest["telemetry"]["sample_count"] == 2
    assert manifest["telemetry"]["cohort_af_kernel_calls"] > 0
    assert manifest["telemetry"]["cohort_af_converged"] > 0
    assert manifest["compatibility"]["multi_sample_joint_merge"] is True
    assert manifest["compatibility"]["gatk_annotation_compatibility"] is True
    print(json.dumps({
        "status": "pass",
        "gatk_records": len(gatk_rows),
        "native_records": len(native_rows),
        "record_text_exact": len(native_rows),
        "header_semantic_exact": True,
        "sample_count": manifest["telemetry"]["sample_count"],
        "cohort_iterations": manifest["telemetry"]["cohort_af_iterations"],
        "cohort_converged": manifest["telemetry"]["cohort_af_converged"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
