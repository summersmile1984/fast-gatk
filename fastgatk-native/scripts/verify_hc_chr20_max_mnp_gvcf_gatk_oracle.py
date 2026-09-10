#!/usr/bin/env python3
"""Byte oracle for GATK HC gVCF spanning-deletion output after MNP merging.

At chr20:10020680, --max-mnp-distance 1 creates a preceding multi-base
EventMap candidate.  In reference-confidence mode GATK's emitted
<NON_REF> allele records that predecessor as an upstream deletion, so the
covered event retains its symbolic `*` allele.  The complete REF/ALT/*/
<NON_REF> PL matrix, QUAL, AD and phasing fields are all part of the output
contract; ordinary VCF intentionally does not emit this `*` row.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def rows(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def row_at(records: list[list[str]], position: int) -> list[str]:
    for record in records:
        if int(record[1]) == position:
            return record
    raise AssertionError(f"missing gVCF row at {position}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    reference = root / "fixtures/chr20/ref20mnp.fasta"
    bam = root / "fixtures/chr20/mnp.bam"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    region = "20:10019901-10020710"

    if not native.exists() or not reference.exists() or not bam.exists():
        raise SystemExit("missing native build or chr20 max-MNP fixture")
    if not gatk.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr20-max-mnp-gvcf-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf"
        native_gvcf = work / "native.g.vcf"
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        feature = work / "spanning.multiallelic.vcf"
        overlapping_feature = work / "spanning.overlapping-mnp.vcf"
        gatk_forced_gvcf = work / "gatk.alleles.g.vcf"
        native_forced_gvcf = work / "native.alleles.g.vcf"
        gatk_overlapping_gvcf = work / "gatk.overlapping-alleles.g.vcf"
        native_overlapping_gvcf = work / "native.overlapping-alleles.g.vcf"
        gatk_capped_gvcf = work / "gatk.max-alt.g.vcf"
        native_capped_gvcf = work / "native.max-alt.g.vcf"
        manifest = work / "native.json"
        # At the existing spanning-deletion locus, retain the genuinely
        # supported AT through GenotypeGivenAlleles.  The full gVCF event is
        # AT,*,<NON_REF>; at a cap of one the oracle chooses AT and removes
        # `*`, exercising symbolic Number=G subsetting rather than the simpler
        # no-spanning-deletion max-ALT case.
        feature.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=20,length=63025520>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "20\t10020680\t.\tCA\tAT\t.\tPASS\t.\n",
            encoding="utf-8")
        bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
        tabix = root / "third_party/htslib-build/htslib-src/tabix"
        subprocess.run([str(bgzip), "-f", str(feature)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        feature = feature.with_suffix(feature.suffix + ".gz")
        subprocess.run([str(tabix), "-f", "-p", "vcf", str(feature)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # CG is a same-start sibling of AT.  Here AT is already represented
        # by a retained MNP/EventMap haplotype.  GATK tries to inject CG only
        # into its five selected base haplotypes and drops it when each is
        # occupied by the overlapping MNP.  This pins the source boundary:
        # a post-assembly synthetic candidate would incorrectly emit CG.
        overlapping_feature.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=20,length=63025520>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "20\t10020680\t.\tCA\tAT,CG\t.\tPASS\t.\n",
            encoding="utf-8")
        subprocess.run([str(bgzip), "-f", str(overlapping_feature)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        overlapping_feature = overlapping_feature.with_suffix(overlapping_feature.suffix + ".gz")
        subprocess.run([str(tabix), "-f", "-p", "vcf", str(overlapping_feature)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--max-mnp-distance", "1", "-O", str(gatk_gvcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1", "--threads", "1",
            "-O", str(native_gvcf), "--output-manifest", str(manifest),
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--max-mnp-distance", "1", "--alleles", str(overlapping_feature),
            "-O", str(gatk_overlapping_gvcf), "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1",
            "--alleles", str(overlapping_feature), "--threads", "1",
            "-O", str(native_overlapping_gvcf), "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--max-mnp-distance", "1", "--alleles", str(feature),
            "-O", str(gatk_forced_gvcf), "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1",
            "--alleles", str(feature), "--threads", "1", "-O", str(native_forced_gvcf),
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--max-mnp-distance", "1", "--alleles", str(feature),
            "--max-alternate-alleles", "1", "-O", str(gatk_capped_gvcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        capped_native_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1",
            "--alleles", str(feature), "--max-alternate-alleles", "1", "--threads", "1",
            "-O", str(native_capped_gvcf), "--add-output-vcf-command-line", "false",
        ], text=True, capture_output=True)
        if capped_native_run.returncode != 0:
            raise RuntimeError(
                "native max-ALT spanning-deletion invocation failed:\n" +
                capped_native_run.stderr[-5000:])
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--max-mnp-distance", "1", "-O", str(gatk_vcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--max-mnp-distance", "1", "--threads", "1",
            "-O", str(native_vcf), "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        assert gatk_gvcf.read_bytes() == native_gvcf.read_bytes(), (
            "chr20 max-MNP-distance=1 gVCF is not byte-identical with provenance disabled")
        assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
            "chr20 max-MNP-distance=1 VCF is not byte-identical with provenance disabled")
        records = rows(native_gvcf)
        assert len(records) == 48, f"unexpected GATK-compatible gVCF row count: {len(records)}"

        spanning = row_at(records, 10020680)
        assert spanning[3:8] == ["CA", "AT,*,<NON_REF>", "154.25", ".",
                                  "DP=8;ExcessHet=0.0000;MLEAC=1,1,0;"
                                  "MLEAF=0.500,0.500,0.00;RAW_MQandDP=28800,8"]
        assert spanning[8:] == [
            "GT:AD:DP:GQ:PGT:PID:PL:PS:SB",
            "1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:"
            "249,84,99,165,0,159,252,99,168,264:10020679:0,0,3,3",
        ]
        ordinary_spanning = row_at(rows(native_vcf), 10020680)
        assert ordinary_spanning[3:5] == ["CA", "AT"], (
            "ordinary VCF must suppress the otherwise spurious spanning-deletion allele")
        native_forced_spanning = row_at(rows(native_forced_gvcf), 10020680)
        gatk_forced_spanning = row_at(rows(gatk_forced_gvcf), 10020680)
        assert native_forced_spanning == gatk_forced_spanning, (
            "the unsubsetted forced spanning candidate must first agree with GATK: "
            f"native={native_forced_spanning!r} gatk={gatk_forced_spanning!r}")
        native_overlapping_spanning = row_at(rows(native_overlapping_gvcf), 10020680)
        gatk_overlapping_spanning = row_at(rows(gatk_overlapping_gvcf), 10020680)
        assert native_overlapping_spanning == gatk_overlapping_spanning, (
            "an overlapping --alleles sibling must follow GATK's bounded base-haplotype "
            f"injection result: native={native_overlapping_spanning!r} "
            f"gatk={gatk_overlapping_spanning!r}")
        assert native_overlapping_spanning[4] == "AT,*,<NON_REF>", (
            "the non-injectable overlapping CG sibling must not create a GATK-absent ALT")
        gatk_capped_spanning = row_at(rows(gatk_capped_gvcf), 10020680)
        native_capped_spanning = row_at(rows(native_capped_gvcf), 10020680)
        assert native_capped_spanning == gatk_capped_spanning, (
            "GVCF --max-alternate-alleles must subset the original "
            "REF/ALT/*/<NON_REF> Number=G matrix exactly like GATK: "
            f"native={native_capped_spanning!r} gatk={gatk_capped_spanning!r}")
        assert native_capped_spanning[4] == "AT,<NON_REF>", (
            "max-ALT spanning candidate must remove the unselected `*` allele")

        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["max_mnp_distance"] == 1
        assert telemetry["gvcf_requested"] is True
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "max_mnp_distance": 1,
            "gvcf_records": len(records),
            "gvcf_byte_identical_without_provenance": True,
            "vcf_byte_identical_without_provenance": True,
            "spanning_deletion_row_exact": True,
            "overlapping_given_allele_bounded_haplotype_injection_exact": True,
            "max_alternate_alleles_spanning_deletion_row_exact": True,
            "ordinary_vcf_spanning_deletion_suppressed": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
