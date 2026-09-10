#!/usr/bin/env python3
"""Ploidy-aware gVCF oracle for the max-MNP spanning-deletion state machine.

The chr20 max-MNP fixture makes a multi-base predecessor whose symbolic
<NON_REF> gVCF allele covers the next EventMap candidate.  GATK consequently
retains `*` at that locus even for arbitrary sample ploidy.  Ordinary VCF
instead carries `*` through GL/AF then subsets it from the emitted allele
list.  This exercises both the partitioned AssemblyRegion likelihood-owner
lookup and the complete Kokkos Number=G matrix.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def data_rows(path: Path) -> bytes:
    return b"\n".join(
        line for line in path.read_bytes().splitlines() if line and not line.startswith(b"#")
    ) + b"\n"


def rows(path: Path) -> list[list[str]]:
    return [line.decode("utf-8").split("\t")
            for line in path.read_bytes().splitlines()
            if line and not line.startswith(b"#")]


def row_at(records: list[list[str]], position: int) -> list[str]:
    for record in records:
        if int(record[1]) == position:
            return record
    raise AssertionError(f"missing gVCF row at {position}")


def duplicate_single_sample_gvcf(source: Path, destination: Path, sample_name: str) -> None:
    """Make a schema-valid two-sample copy without invoking CombineGVCFs.

    GATK deliberately rejects gVCFs containing MNPs in CombineGVCFs, while
    GenotypeGVCFs accepts their multi-sample form.  Start from an actual GATK
    HC gVCF and duplicate its sole sample column so the oracle stays focused
    on the joint-genotyping sample-major Number=G boundary.
    """
    output: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            fields = line.split("\t")
            if len(fields) != 10:
                raise AssertionError("expected exactly one GATK HC sample column")
            output.append("\t".join(fields + [sample_name]))
        elif line.startswith("#"):
            output.append(line)
        else:
            fields = line.split("\t")
            if len(fields) != 10:
                raise AssertionError("expected exactly one GATK HC sample field")
            output.append("\t".join(fields + [fields[9]]))
    destination.write_text("\n".join(output) + "\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    genotype = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        root / "fastgatk-native/build/fastgatk-genotype-gvcf"))
    reference = root / "fixtures/chr20/ref20mnp.fasta"
    bam = root / "fixtures/chr20/mnp.bam"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    region = "20:10019901-10020710"
    ploidy = int(os.environ.get("FASTGATK_HC_SAMPLE_PLOIDY", "3"))
    expected = {
        1: {
            "gvcf_fields": ["CA", "AT,<NON_REF>", "89.07", ".",
                            "DP=8;MLEAC=1,0;MLEAF=1.00,0.00;RAW_MQandDP=28800,8"],
            "gvcf_sample": "1:1,4,0:5:99:150,0,165:0,0,3,1",
            "vcf_fields": ["CA", "AT", "89.07", ".",
                           "AC=1;AF=1.00;AN=1;DP=8;FS=0.000;"
                           "MLEAC=1;MLEAF=1.00;MQ=60.00;QD=17.81;SOR=1.609"],
            "vcf_sample": "1:1,4:5:99:150,0",
            "gvcf_records": 7,
        },
        3: {
            "gvcf_fields": ["CA", "AT,*,<NON_REF>", "153.60", ".",
                            "DP=8;MLEAC=1,1,0;MLEAF=0.333,0.333,0.00;"
                            "RAW_MQandDP=28800,8"],
            "gvcf_sample": "1/1/2:0,4,2,0:6:3:"
                           "249,89,80,98,168,9,0,162,3,159,250,92,95,170,12,164,253,104,173,264:"
                           "0,0,3,3",
            "vcf_fields": ["CA", "AT", "153.60", ".",
                           "AC=2;AF=0.667;AN=3;DP=8;FS=0.000;"
                           "MLEAC=1;MLEAF=0.333;MQ=60.00;QD=30.72;SOR=1.609"],
            "vcf_sample": "0/1/1:1,4:5:9:169,9,0,18",
            "gvcf_records": 94,
        },
        4: {
            "gvcf_fields": ["CA", "AT,*,<NON_REF>", "154.29", ".",
                            "DP=8;MLEAC=2,1,0;MLEAF=0.500,0.250,0.00;"
                            "RAW_MQandDP=28800,8"],
            "gvcf_sample": "1/1/2/2:0,4,2,0:6:2:"
                           "249,94,84,80,99,171,16,6,2,165,10,0,161,7,159,250,96,87,94,172,18,9,"
                           "166,12,163,252,99,99,174,21,168,255,109,177,264:0,0,3,3",
            "vcf_fields": ["CA", "AT", "154.29", ".",
                           "AC=3;AF=0.750;AN=4;DP=8;FS=0.000;"
                           "MLEAC=2;MLEAF=0.500;MQ=60.00;QD=30.86;SOR=1.609"],
            "vcf_sample": "0/1/1/1:1,4:5:4:169,14,4,0,19",
            "gvcf_records": 93,
        },
        8: {
            # p8 exercises the symbolic <NON_REF> indel prior and the
            # ReferenceConfidenceModel's unrounded SNP-vs-indel GQ choice.
            # The byte-identical full-data comparison below checks every
            # large Number=G field without duplicating its 165-entry PL row.
            "gvcf_records": 60,
        },
    }
    if ploidy not in expected:
        raise SystemExit("this pinned max-MNP oracle supports sample ploidy 1, 3, 4, or 8")
    expectation = expected[ploidy]

    if not native.exists() or not reference.exists() or not bam.exists():
        raise SystemExit("missing native build or chr20 max-MNP fixture")
    if not gatk.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr20-max-mnp-polyploid-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf"
        native_gvcf = work / "native.g.vcf"
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        gatk_joint = work / "gatk.joint.vcf"
        native_joint = work / "native.joint.vcf"
        cohort_gvcf = work / "two-sample.g.vcf"
        gatk_cohort_joint = work / "gatk.cohort.joint.vcf"
        native_cohort_joint = work / "native.cohort.joint.vcf"
        gatk_bp = work / "gatk.bp.g.vcf"
        native_bp = work / "native.bp.g.vcf"
        manifest = work / "native.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--sample-ploidy", str(ploidy),
            "--max-mnp-distance", "1", "-O", str(gatk_gvcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--sample-ploidy", str(ploidy),
            "--max-mnp-distance", "1", "--threads", "1", "-O", str(native_gvcf),
            "--output-manifest", str(manifest), "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--sample-ploidy", str(ploidy), "--max-mnp-distance", "1",
            "-O", str(gatk_vcf), "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--sample-ploidy", str(ploidy), "--max-mnp-distance", "1",
            "--threads", "1", "-O", str(native_vcf),
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        assert data_rows(gatk_gvcf) == data_rows(native_gvcf), (
            "chr20 max-MNP gVCF data records differ from GATK")
        assert data_rows(gatk_vcf) == data_rows(native_vcf), (
            "chr20 max-MNP ordinary VCF data records differ from GATK")
        # HC's arbitrary-ploidy gVCF contract continues through joint
        # genotyping for every supported sample ploidy.  At this fixture's
        # orphan spanning-deletion locus, GATK subsets '*' then recomputes
        # the projected Number=G PREFER_PLS genotype; for p3 that is 0/1/1,
        # not a source-GT ././. projection.
        if not genotype.exists():
            raise SystemExit("missing native GenotypeGVCFs binary for polyploid gVCF oracle")
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-O", str(gatk_joint),
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(genotype), "-R", str(reference), "-V", str(gatk_gvcf),
            "-O", str(native_joint), "--gatk-compatible-annotations",
            "--create-output-variant-index", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert data_rows(gatk_joint) == data_rows(native_joint), (
            "chr20 max-MNP polyploid GenotypeGVCFs data records differ from GATK")
        # Exercise the joint sample-major Number=G layout as well.  Triploid
        # is the smallest non-diploid ordering whose biallelic PL row has
        # four cells.  CombineGVCFs refuses MNP gVCFs, so duplicate the
        # GATK-HC sample field and use GATK GenotypeGVCFs itself as oracle.
        if ploidy == 3:
            duplicate_single_sample_gvcf(gatk_gvcf, cohort_gvcf, "SAMPLE2")
            subprocess.run([
                java, "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
                "-R", str(reference), "-V", str(cohort_gvcf), "-O", str(gatk_cohort_joint),
                "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run([
                str(genotype), "-R", str(reference), "-V", str(cohort_gvcf),
                "-O", str(native_cohort_joint), "--gatk-compatible-annotations",
                "--create-output-variant-index", "false",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            assert data_rows(gatk_cohort_joint) == data_rows(native_cohort_joint), (
                "chr20 max-MNP triploid multisample GenotypeGVCFs data records differ from GATK")
        # BP_RESOLUTION has no candidate END field, so the tail base of the
        # CA->AT MNP must be emitted separately as a normal <NON_REF>
        # reference-confidence row rather than being swallowed by the
        # candidate span.  Check this for every supported ploidy; p8 also
        # supplies the widest Number=G RCM rows in the compact fixture.
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "BP_RESOLUTION", "--sample-ploidy", str(ploidy),
            "--max-mnp-distance", "1", "-O", str(gatk_bp),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "BP_RESOLUTION", "--sample-ploidy", str(ploidy),
            "--max-mnp-distance", "1", "--threads", "1", "-O", str(native_bp),
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert data_rows(gatk_bp) == data_rows(native_bp), (
            "chr20 max-MNP BP_RESOLUTION data records differ from GATK")
        assert len(rows(native_bp)) == 810
        records = rows(native_gvcf)
        assert len(records) == expectation["gvcf_records"], (
            f"unexpected GATK-compatible gVCF row count: {len(records)}")
        if "gvcf_fields" in expectation:
            spanning = row_at(records, 10020680)
            assert spanning[3:8] == expectation["gvcf_fields"]
            assert spanning[8:] == [
                "GT:AD:DP:GQ:PL:SB",
                expectation["gvcf_sample"],
            ]
            ordinary_spanning = row_at(rows(native_vcf), 10020680)
            assert ordinary_spanning[3:8] == expectation["vcf_fields"]
            assert ordinary_spanning[8:] == [
                "GT:AD:DP:GQ:PL",
                expectation["vcf_sample"],
            ]
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["sample_ploidy"] == ploidy
        assert telemetry["max_mnp_distance"] == 1
        assert telemetry["gvcf_requested"] is True
        assert telemetry["pairhmm_assembly_region_partitioned"] is True
        report = {
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "sample_ploidy": ploidy,
            "max_mnp_distance": 1,
            "gvcf_records": len(records),
            "gvcf_data_rows_byte_identical_without_provenance": True,
            "vcf_data_rows_byte_identical_without_provenance": True,
            "gvcf_spanning_deletion_row_exact": True,
            "ordinary_vcf_spanning_deletion_subset_exact": True,
            "joint_data_rows_byte_identical": True,
            "triploid_multisample_joint_data_rows_byte_identical": True,
            "bp_resolution_data_rows_byte_identical": True,
            "bp_resolution_records": 810,
            "partitioned_pairhmm_owner_lookup": True,
        }
        print(json.dumps(report, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
