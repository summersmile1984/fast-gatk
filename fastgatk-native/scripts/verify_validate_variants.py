#!/usr/bin/env python3
"""Verify ValidateVariants reference, allele, chromosome-count and GVCF gates."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>
##INFO=<ID=END,Number=1,Type=Integer,Description=End>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-validate-variants"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-validate-variants-") as directory:
        work = Path(directory)
        reference = work / "ref.fa"
        reference.write_text(">chr1\nACG" + "A" * 97 + "\n", encoding="utf-8")
        (work / "ref.fa.fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")
        (work / "ref.dict").write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
        valid = work / "valid.vcf"
        valid.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n"
            + "chr1\t2\trs2\tC\tT\t50\tPASS\tAC=2;AN=2\tGT\t1/1\n"
            + "chr1\t3\trs3\tG\tGA\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n"
            + "chr1\t4\trs4\tA\tC\t50\tLowQual\tAC=1;AN=2\tGT\t0/1\n",
            encoding="utf-8",
        )
        report = work / "valid.report.tsv"
        manifest = work / "valid.manifest.json"
        result = run(binary, [
            "-V", str(valid), "-R", str(reference), "-L", "chr1:1-3",
            "-O", str(report), "--output-manifest", str(manifest),
        ])
        assert result.returncode == 0, result.stderr
        assert "\tPASS\t4\t3\t0\t0" in report.read_text(encoding="utf-8")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["status"] == "pass"
        assert metadata["telemetry"]["skipped_interval"] == 1
        assert metadata["compatibility"]["ref"] is True
        assert metadata["compatibility"]["field_cardinality"] is True

        # Repeatable -L selectors use GATK's explicit interval-set rule.  The
        # intersection of chr1:1-3 and chr1:3-3 contains only the POS=3
        # record; the manifest and summary expose the selected rule.
        intersection_report = work / "intersection.report.tsv"
        intersection_manifest = work / "intersection.manifest.json"
        intersection_result = run(binary, [
            "-V", str(valid), "-R", str(reference),
            "-L", "chr1:1-3", "-L", "chr1:3-3",
            "--interval-set-rule", "INTERSECTION",
            "-O", str(intersection_report), "--output-manifest", str(intersection_manifest),
        ])
        assert intersection_result.returncode == 0, intersection_result.stderr
        assert "\tPASS\t4\t1\t0\t0" in intersection_report.read_text(encoding="utf-8")
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["skipped_interval"] == 3

        # GATK applies exclusion intervals after inclusion intervals.  Verify
        # that an excluded record is not validated, and that -ip/-ixp expand
        # the respective sets in the same zero-based half-open space.
        excluded_invalid = work / "excluded-invalid.vcf"
        excluded_invalid.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n"
            + "chr1\t2\trs2\tT\tC\t50\tPASS\tAC=2;AN=2\tGT\t1/1\n"
            + "chr1\t3\trs3\tG\tGA\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n",
            encoding="utf-8",
        )
        excluded_manifest = work / "excluded.manifest.json"
        excluded_result = run(binary, [
            "-V", str(excluded_invalid), "-R", str(reference),
            "-L", "chr1:1-3", "-XL", "chr1:2-2",
            "--output-manifest", str(excluded_manifest),
        ])
        assert excluded_result.returncode == 0, excluded_result.stderr
        excluded_metadata = json.loads(excluded_manifest.read_text(encoding="utf-8"))
        assert excluded_metadata["telemetry"]["validated_records"] == 2
        assert excluded_metadata["telemetry"]["skipped_excluded"] == 1
        assert excluded_metadata["compatibility"]["interval_exclusion"] is True

        padded_manifest = work / "padded-exclusion.manifest.json"
        padded_result = run(binary, [
            "-V", str(excluded_invalid), "-R", str(reference),
            "-L", "chr1:2-2", "-ip", "1", "-XL", "chr1:2-2", "-ixp", "1",
            "--output-manifest", str(padded_manifest),
        ])
        assert padded_result.returncode == 0, padded_result.stderr
        padded_metadata = json.loads(padded_manifest.read_text(encoding="utf-8"))
        assert padded_metadata["telemetry"]["validated_records"] == 0
        assert padded_metadata["telemetry"]["skipped_excluded"] == 3
        assert padded_metadata["compatibility"]["interval_padding"] == 1
        assert padded_metadata["compatibility"]["interval_exclusion_padding"] == 1

        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if oracle_guard.oracle_ready('verify_validate_variants.py', java, gatk_jar):
            index_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "IndexFeatureFile",
                "-I", str(excluded_invalid),
            ], text=True, capture_output=True, check=False)
            assert index_result.returncode == 0, index_result.stderr
            gatk_excluded = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "ValidateVariants",
                "-V", str(excluded_invalid), "-R", str(reference),
                "-L", "chr1:1-3", "-XL", "chr1:2-2",
            ], text=True, capture_output=True, check=False)
            assert gatk_excluded.returncode == 0, gatk_excluded.stderr

        cardinality_header = HEADER.replace(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
            "##INFO=<ID=ADINFO,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n",
        )
        cardinality = work / "bad-cardinality.vcf"
        cardinality.write_text(
            cardinality_header
            + "chr1\t1\trs1\tA\tG\t50\tPASS\tAC=1;AN=2;ADINFO=1\tGT:AD:PL\t0/1:1:0,1\n",
            encoding="utf-8",
        )
        cardinality_manifest = work / "cardinality.manifest.json"
        cardinality_result = run(binary, [
            "-V", str(cardinality), "-R", str(reference), "--warn-on-errors",
            "--output-manifest", str(cardinality_manifest),
        ])
        assert cardinality_result.returncode == 0, cardinality_result.stderr
        cardinality_metadata = json.loads(cardinality_manifest.read_text(encoding="utf-8"))
        assert cardinality_metadata["telemetry"]["field_cardinality_errors"] >= 3
        assert cardinality_metadata["compatibility"]["field_cardinality"] is True
        excluded_format = run(binary, [
            "-V", str(cardinality), "-R", str(reference),
            "--validation-type-to-exclude", "FORMAT",
        ])
        assert excluded_format.returncode == 0, excluded_format.stderr

        invalid = work / "invalid.vcf"
        invalid.write_text(
            HEADER + "chr1\t1\trs1\tA\tG,G\t50\tPASS\tAC=4;AN=4\tGT\t0/1\n",
            encoding="utf-8",
        )
        warning_manifest = work / "warning.manifest.json"
        warning = run(binary, [
            "-V", str(invalid), "-R", str(reference), "--warn-on-errors",
            "--output-manifest", str(warning_manifest),
        ])
        assert warning.returncode == 0, warning.stderr
        warning_metadata = json.loads(warning_manifest.read_text(encoding="utf-8"))
        assert warning_metadata["status"] == "pass"
        assert warning_metadata["telemetry"]["warnings"] >= 1

        failure = run(binary, ["-V", str(invalid), "-R", str(reference)])
        assert failure.returncode != 0
        assert "VALIDATION_FAILURE" in failure.stderr

        # GATK's strict ALLELES validation also rejects a concrete ALT that is
        # never present in any called genotype.  Disable CHR_COUNTS here so
        # the regression isolates the allele-usage contract.
        unused_alt = work / "unused-alt.vcf"
        unused_alt.write_text(
            HEADER + "chr1\t1\trs1\tA\tC,G\t50\tPASS\tAC=1,0;AN=2\tGT\t0/1\n",
            encoding="utf-8",
        )
        unused_failure = run(binary, [
            "-V", str(unused_alt), "-R", str(reference),
            "--validation-type-to-exclude", "CHR_COUNTS",
        ])
        assert unused_failure.returncode != 0
        assert "not observed in any called genotype" in unused_failure.stderr
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if oracle_guard.oracle_ready('verify_validate_variants.py', java, gatk_jar):
            oracle_unused = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "ValidateVariants",
                "-V", str(unused_alt), "--validation-type-to-exclude", "REF",
                "--validation-type-to-exclude", "CHR_COUNTS",
                "--validation-type-to-exclude", "IDS",
            ], text=True, capture_output=True, check=False)
            assert oracle_unused.returncode != 0
            assert "not observed at all" in oracle_unused.stderr

        excluded = run(binary, [
            "-V", str(invalid), "--validation-type-to-exclude", "ALL",
        ])
        assert excluded.returncode == 0, excluded.stderr

        mismatched = work / "mismatched-dictionary.vcf"
        mismatched.write_text(
            HEADER.replace("length=100", "length=101")
            + "chr1\t1\trs1\tA\tG\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n",
            encoding="utf-8",
        )
        dictionary_failure = run(binary, [
            "-V", str(mismatched), "-R", str(reference),
        ])
        assert dictionary_failure.returncode != 0
        assert "sequence dictionary" in dictionary_failure.stderr
        dictionary_disabled = run(binary, [
            "-V", str(mismatched), "-R", str(reference),
            "--disable-sequence-dictionary-validation",
        ])
        assert dictionary_disabled.returncode == 0, dictionary_disabled.stderr

        gvcf = work / "input.g.vcf"
        gvcf.write_text(
            HEADER + "chr1\t1\trs1\tA\tG\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n",
            encoding="utf-8",
        )
        gvcf_manifest = work / "gvcf.manifest.json"
        gvcf_result = run(binary, [
            "-V", str(gvcf), "--gvcf", "--warn-on-errors",
            "--output-manifest", str(gvcf_manifest),
        ])
        assert gvcf_result.returncode == 0, gvcf_result.stderr
        gvcf_metadata = json.loads(gvcf_manifest.read_text(encoding="utf-8"))
        assert gvcf_metadata["compatibility"]["gvcf"] is True
        assert gvcf_metadata["telemetry"]["warnings"] >= 1

        gvcf_header = HEADER.replace(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n",
            "##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n",
        )
        complete_gvcf = work / "complete.g.vcf"
        complete_gvcf.write_text(
            gvcf_header
            + "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=50\tGT\t0/0\n"
            + "chr1\t51\t.\tA\t<NON_REF>\t.\tPASS\tEND=100\tGT\t0/0\n",
            encoding="utf-8",
        )
        complete_manifest = work / "complete-gvcf.manifest.json"
        complete_result = run(binary, [
            "-V", str(complete_gvcf), "-R", str(reference), "--validate-gvcf",
            "--output-manifest", str(complete_manifest),
        ])
        assert complete_result.returncode == 0, complete_result.stderr
        complete_metadata = json.loads(complete_manifest.read_text(encoding="utf-8"))
        assert complete_metadata["telemetry"]["gvcf_blocks"] == 2
        assert complete_metadata["telemetry"]["gvcf_uncovered_gaps"] == 0
        # ValidateVariants.calculateValidationTypesToApply() disables strict
        # ALLELES checks for GVCF mode, because <NON_REF> is a
        # reference-confidence allele rather than an ordinary concrete ALT.
        assert complete_metadata["compatibility"]["alleles"] is False

        # Interval selection is record-span based for GVCFs: the first block
        # starts before chr1:25 but covers it through END=50.  A POS-only
        # selector would incorrectly validate zero records here.
        span_manifest = work / "span-gvcf.manifest.json"
        span_result = run(binary, [
            "-V", str(complete_gvcf), "-R", str(reference), "--validate-gvcf",
            "-L", "chr1:25-25", "--output-manifest", str(span_manifest),
        ])
        assert span_result.returncode == 0, span_result.stderr
        span_metadata = json.loads(span_manifest.read_text(encoding="utf-8"))
        assert span_metadata["telemetry"]["validated_records"] == 1
        assert span_metadata["telemetry"]["skipped_interval"] == 1
        assert span_metadata["telemetry"]["gvcf_blocks"] == 1

        # Coverage gaps must be detected even when both neighboring records
        # are concrete variant records (the old adjacent-reference-block-only
        # check silently missed this case).
        concrete_gap = work / "concrete-gap.g.vcf"
        concrete_gap.write_text(
            gvcf_header
            + "chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT\t0/1\n"
            + "chr1\t3\t.\tG\tA,<NON_REF>\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        concrete_gap_manifest = work / "concrete-gap.manifest.json"
        concrete_gap_result = run(binary, [
            "-V", str(concrete_gap), "-R", str(reference), "--validate-gvcf",
            "--warn-on-errors", "--output-manifest", str(concrete_gap_manifest),
        ])
        assert concrete_gap_result.returncode == 0, concrete_gap_result.stderr
        concrete_gap_metadata = json.loads(concrete_gap_manifest.read_text(encoding="utf-8"))
        assert concrete_gap_metadata["telemetry"]["gvcf_uncovered_gaps"] >= 1

        incomplete_gvcf = work / "incomplete.g.vcf"
        incomplete_gvcf.write_text(
            gvcf_header
            + "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=50\tGT\t0/0\n",
            encoding="utf-8",
        )
        incomplete_result = run(binary, [
            "-V", str(incomplete_gvcf), "-R", str(reference), "--validate-gvcf",
        ])
        assert incomplete_result.returncode != 0
        assert "does not reach the end" in incomplete_result.stderr

    print(json.dumps({"status": "pass", "validated_records": 3, "warning_cases": 2}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
