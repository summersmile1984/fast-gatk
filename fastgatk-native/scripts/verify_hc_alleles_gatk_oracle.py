#!/usr/bin/env python3
"""Pinned GATK oracle for HaplotypeCaller --alleles / GenotypeGivenAlleles.

The feature marks a SNP in a reference-only pileup.  It therefore verifies
that the Host turns the feature into a Kokkos activity seed and injects the
concrete allele into the local haplotype set; parsing the VCF alone would not
produce a call in this fixture.
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def write_fixture(work: Path) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    random.seed(91)
    reference_text = "".join(random.choice("ACGT") for _ in range(1200))
    reference = work / "ref.fa"
    reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
    # The FASTA is deliberately a single sequence line so its .fai geometry
    # stays deterministic without relying on a local samtools executable.
    (work / "ref.fa.fai").write_text(
        f"chr1\t{len(reference_text)}\t6\t{len(reference_text)}\t{len(reference_text) + 1}\n",
        encoding="utf-8")
    position = 600  # 1-based VCF coordinate, central in every read.
    ref_base = reference_text[position - 1]
    alt_base = next(base for base in "ACGT" if base != ref_base)
    sam = work / "input.sam"
    sequence = reference_text[499:699]
    quality = "I" * len(sequence)
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1200",
             "@RG\tID:rg\tSM:S1"]
    for index in range(20):
        lines.append(
            f"ref{index}\t0\tchr1\t500\t60\t200M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:rg")
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def write_feature(path: Path, filter_text: str, site: int = position,
                      alternates: str | None = None,
                      allele_pair: tuple[str, str] | None = None) -> Path:
        site_ref = reference_text[site - 1]
        site_alt = next(base for base in "ACGT" if base != site_ref)
        if alternates is None:
            alternates = site_alt
        if allele_pair is not None:
            site_ref, alternates = allele_pair
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1200>\n"
            "##FILTER=<ID=LowQual,Description=deliberately filtered feature>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"chr1\t{site}\t.\t{site_ref}\t{alternates}\t.\t{filter_text}\t.\n",
            encoding="utf-8")
        root = Path(__file__).resolve().parents[2]
        bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
        tabix = root / "third_party/htslib-build/htslib-src/tabix"
        subprocess.run([str(bgzip), "-f", str(path)], check=True, stdout=subprocess.DEVNULL)
        compressed = path.with_suffix(path.suffix + ".gz")
        subprocess.run([str(tabix), "-f", "-p", "vcf", str(compressed)], check=True,
                       stdout=subprocess.DEVNULL)
        return compressed

    passing = write_feature(work / "alleles.pass.vcf", "PASS")
    filtered = write_feature(work / "alleles.filtered.vcf", "LowQual")
    uncovered = write_feature(work / "alleles.uncovered.vcf", "PASS", 900)
    ref_base = reference_text[position - 1]
    multi_alt = ",".join(base for base in "ACGT" if base != ref_base)
    multiallelic = write_feature(work / "alleles.multiallelic.vcf", "PASS",
                                 alternates=multi_alt)
    indel_ref = reference_text[position - 1:position + 1]
    indel = write_feature(work / "alleles.indel.vcf", "PASS",
                          allele_pair=(indel_ref, indel_ref[0]))
    return reference, sam, work / "input.bam", passing, filtered, uncovered, multiallelic, indel


def run_pair(java: str, gatk: Path, native: Path, reference: Path, bam: Path,
             feature: Path, output_stem: str, work: Path,
             filtered_flag: str | None = None, interval: str = "chr1:500-700",
             extra: list[str] | None = None
             ) -> tuple[list[list[str]], list[list[str]], dict]:
    java_vcf = work / f"gatk.{output_stem}.vcf"
    native_vcf = work / f"native.{output_stem}.vcf"
    manifest = work / f"native.{output_stem}.json"
    common = ["-R", str(reference), "-I", str(bam), "-L", interval,
              "--alleles", str(feature), "--min-pruning", "1",
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false"]
    if filtered_flag is not None:
        common.append(filtered_flag)
    if extra is not None:
        common.extend(extra)
    subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                    "-O", str(java_vcf)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run([str(native), *common, "--threads", "2", "--output-manifest",
                    str(manifest), "-O", str(native_vcf)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return records(java_vcf), records(native_vcf), json.loads(manifest.read_text(encoding="utf-8"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    if not native.exists() or not gatk.exists():
        oracle_guard.oracle_not_verified('verify_hc_alleles_gatk_oracle.py', gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit("missing HC binary or pinned GATK 4.6.2.0 jar")
        print(json.dumps({"status": "skipped", "reason": "GATK jar or native binary absent"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-alleles-oracle-") as directory:
        work = Path(directory)
        reference, sam, bam, passing, filtered, uncovered, multiallelic, indel = write_fixture(work)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
            "-R", str(reference), "-O", str(work / "ref.dict"),
            "--TRUNCATE_NAMES_AT_WHITESPACE", "true",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam), "-O", str(bam),
            "-SO", "coordinate", "--CREATE_INDEX", "true",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        pass_java, pass_native, pass_manifest = run_pair(
            java, gatk, native, reference, bam, passing, "pass", work)
        assert pass_java == pass_native, (pass_java, pass_native, pass_manifest["telemetry"])
        assert len(pass_native) == 1, pass_native
        assert pass_native[0][1] == "600" and pass_native[0][9].startswith("0/0:"), pass_native
        telemetry = pass_manifest["telemetry"]
        assert telemetry["active_loci"] >= 1 and telemetry["assembly_regions"] == 1, telemetry
        compatibility = pass_manifest["compatibility"]
        assert compatibility["alleles_feature"] is True
        assert compatibility["forced_alleles"] == 1
        assert compatibility["force_call_filtered_alleles"] is False

        # HaplotypeCaller ignores filtered given-allele records by default.
        skipped_java, skipped_native, skipped_manifest = run_pair(
            java, gatk, native, reference, bam, filtered, "filtered-default", work)
        assert skipped_java == skipped_native == [], (skipped_java, skipped_native)
        assert skipped_manifest["telemetry"]["active_loci"] == 0

        # GATK exposes --force-call-filtered-alleles; its historical alias is
        # --genotype-filtered-alleles.  Native accepts both spellings and must
        # make them observationally identical to a passing feature record.
        for label, option in (("force-call-filtered-alleles", "--force-call-filtered-alleles"),
                              ("genotype-filtered-alleles", "--genotype-filtered-alleles")):
            java_rows, native_rows, manifest = run_pair(
                java, gatk, native, reference, bam, filtered, label, work, option)
            assert java_rows == native_rows == pass_native, (label, java_rows, native_rows)
            assert manifest["telemetry"]["active_loci"] >= 1
            assert manifest["compatibility"]["force_call_filtered_alleles"] is True

        # GATK's feature-driven activity state also exists at an otherwise
        # empty pileup. Keep the Java row visible while the native sparse-locus
        # projection is being implemented.
        uncovered_java, uncovered_native, uncovered_manifest = run_pair(
            java, gatk, native, reference, bam, uncovered, "uncovered", work,
            interval="chr1:800-1000")
        assert uncovered_java == uncovered_native, (
            uncovered_java, uncovered_native, uncovered_manifest["telemetry"])
        assert uncovered_java == []

        # One feature record may contain several concrete ALT alleles. GATK
        # splits it into EventMap events but emits one multi-allelic site.
        multi_java, multi_native, multi_manifest = run_pair(
            java, gatk, native, reference, bam, multiallelic, "multiallelic", work)
        assert multi_java == multi_native, (multi_java, multi_native, multi_manifest["telemetry"])
        assert len(multi_native) == 1 and multi_native[0][4].count(",") == 2, multi_native

        # Given alleles are EventMap objects, not just SNP labels. Keep an
        # anchored deletion on the same reference-only pileup in this oracle.
        indel_java, indel_native, indel_manifest = run_pair(
            java, gatk, native, reference, bam, indel, "indel", work)
        assert indel_java == indel_native, (indel_java, indel_native, indel_manifest["telemetry"])
        assert len(indel_native) == 1 and len(indel_native[0][3]) == 2, indel_native

        # gVCF must retain the forced concrete ALT beside <NON_REF> and
        # preserve the same reference-confidence block boundaries.
        gvcf_java, gvcf_native, gvcf_manifest = run_pair(
            java, gatk, native, reference, bam, passing, "gvcf", work,
            extra=["-ERC", "GVCF"])
        assert gvcf_java == gvcf_native, (gvcf_java, gvcf_native, gvcf_manifest["telemetry"])
        assert any(row[1] == "600" and row[4].endswith(",<NON_REF>")
                   for row in gvcf_native), gvcf_native

        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0",
            "fixture": "20x reference-only 200M reads; one central forced SNP",
            "java_native_records_exact": True,
            "forced_feature_active": True,
            "filtered_feature_default_ignored": True,
            "force_call_filtered_alleles_and_legacy_alias_gatk_exact": True,
            "empty_pileup_feature_no_call_gatk_exact": True,
            "multi_alt_and_anchored_indel_gatk_exact": True,
            "gvcf_forced_alt_and_non_ref_gatk_exact": True,
            "scope": "HaplotypeCaller GenotypeGivenAlleles active-site and allele injection boundary",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
