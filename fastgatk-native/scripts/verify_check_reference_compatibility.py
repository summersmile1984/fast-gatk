#!/usr/bin/env python3
"""GATK 4.6.2.0 contract/oracle checks for reference compatibility."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_CHECK_REFERENCE_COMPATIBILITY_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-check-reference-compatibility"),
    ))
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    base = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/reference"
    compare = base / "CompareReferences"
    checks = [
        (base / "CheckReferenceCompatibility/reads_data_source_test1_withmd5s.bam",
         compare / "hg19mini.fasta", "expected.testReferenceCompatibilityBAMWithMD5s_exactmatch.table"),
        (base / "CheckReferenceCompatibility/reads_data_source_test1_withmd5s_missingchr1.bam",
         compare / "hg19mini.fasta", "expected.testReferenceCompatibilityBAMWithMD5s_subset.table"),
        (base / "CheckReferenceCompatibility/reads_data_source_test1_withmd5s.bam",
         compare / "hg19mini_missingchr3.fasta", "expected.testReferenceCompatibilityBAMWithMD5s_notcompatible.table"),
        (base / "CheckReferenceCompatibility/reads_data_source_test1_withoutmd5s.bam",
         compare / "hg19mini.fasta", "expected.testReferenceCompatibilityBAMWithoutMD5s_compatible.table"),
    ]
    with tempfile.TemporaryDirectory(prefix="fastgatk-check-reference-") as directory:
        work = pathlib.Path(directory)
        for index, (query, reference, expected_name) in enumerate(checks):
            output = work / f"compat-{index}.table"
            manifest = work / f"compat-{index}.manifest.json"
            result = subprocess.run([
                str(native), "-I", str(query), "-refcomp", str(reference),
                "-O", str(output), "--output-manifest", str(manifest),
            ], check=True, text=True, capture_output=True)
            expected = (base / "CheckReferenceCompatibility" / expected_name).read_text(encoding="ascii")
            if output.read_text(encoding="ascii") != expected:
                raise AssertionError({"case": index, "expected": expected, "actual": output.read_text(encoding="ascii")})
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            assert metadata["status"] == "contract-compatible"
            assert metadata["query_md5s_present"] is ("WithoutMD5s" not in expected_name)
            assert json.loads(result.stderr.splitlines()[-1])["status"] == "contract-compatible"

        vcf_query = base / "CheckReferenceCompatibility/example_variants_withSequenceDict_withmd5.vcf"
        vcf_output = work / "compat-vcf.table"
        subprocess.run([
            str(native), "-V", str(vcf_query), "-refcomp", str(compare / "hg19mini.fasta"),
            "-O", str(vcf_output),
        ], check=True, text=True, capture_output=True)
        expected_vcf = (base / "CheckReferenceCompatibility/expected.testReferenceCompatibilityVCFWithMD5s.table").read_text(encoding="ascii")
        if vcf_output.read_text(encoding="ascii") != expected_vcf:
            raise AssertionError("VCF compatibility table differs from GATK fixture")

        # HTSlib recognizes VCF content even when a workflow/object-store
        # staging layer gives it a non-standard `.bgz` suffix.  This used to
        # fall through to the BAM parser because dispatch relied on filename
        # extensions; keep the same GATK table contract for the content-
        # detected path.
        vcf_bgz_query = work / "example_variants_withSequenceDict_withmd5.vcf.bgz"
        vcf_bgz_query.write_bytes(vcf_query.read_bytes())
        vcf_bgz_output = work / "compat-vcf-bgz.table"
        vcf_bgz_result = subprocess.run([
            str(native), "-V", str(vcf_bgz_query), "-refcomp", str(compare / "hg19mini.fasta"),
            "-O", str(vcf_bgz_output),
        ], check=True, text=True, capture_output=True)
        expected_vcf_bgz = expected_vcf.replace(
            "#Current Reference: example_variants_withSequenceDict_withmd5.vcf",
            "#Current Reference: example_variants_withSequenceDict_withmd5.vcf.bgz",
        )
        if vcf_bgz_output.read_text(encoding="ascii") != expected_vcf_bgz:
            raise AssertionError("content-detected VCF compatibility table differs from GATK fixture")
        assert '"status":"contract-compatible"' in vcf_bgz_result.stderr

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_oracle = False
        if java.exists() and jar.exists():
            java_output = work / "java.table"
            java_result = subprocess.run([
                str(java), "-jar", str(jar), "CheckReferenceCompatibility", "-I", str(checks[0][0]),
                "-refcomp", str(checks[0][1]), "-O", str(java_output),
            ], check=True, text=True, capture_output=True)
            if java_output.read_text(encoding="ascii") != (base / "CheckReferenceCompatibility/expected.testReferenceCompatibilityBAMWithMD5s_exactmatch.table").read_text(encoding="ascii"):
                raise AssertionError("native/Java compatibility tables differ")
            java_oracle = True
        print(json.dumps({"status": "pass", "bam_cases": len(checks), "vcf_cases": 2,
                          "content_detected_vcf": True, "java_oracle": java_oracle}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
