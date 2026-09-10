#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for the strict HC min-pruning boundary.

The chr17 fixture contains low-support alternates whose fate changes with the
pruning floor.  This check keeps both the strict min-pruning=3 decision and
the min-pruning=1 low-depth annotation path tied to the release-pinned Java
implementation.  The latter is important because GATK's exact Mann-Whitney
branch emits finite rank-sum values with one REF and one ALT observation.  It
also covers the reference-confidence boundary where an EventMap site can be a
concrete gVCF record despite its final genotype being 0/0.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def records(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in
            path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def key(record: list[str]) -> tuple[str, int, str, str]:
    return record[0], int(record[1]), record[3], record[4]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    region = "17:69000-70000"
    min_pruning = "3"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or strict HC fixture")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-min-pruning-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_vcf),
            "-L", region, "--min-pruning", min_pruning,
            "--num-pruning-samples", "1", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_vcf), "--output-manifest", str(manifest),
            "--threads", "2", "--min-pruning", min_pruning,
            "--num-pruning-samples", "1",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gatk_records = records(gatk_vcf)
        native_records = records(native_vcf)
        assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
            "strict min-pruning VCF is not byte-identical with provenance disabled")
        gatk_keys = {key(record) for record in gatk_records}
        native_keys = {key(record) for record in native_records}
        assert gatk_keys == native_keys, (
            f"strict HC call-set mismatch: gatk-only={sorted(gatk_keys - native_keys)} "
            f"native-only={sorted(native_keys - gatk_keys)}")
        expected = {
            ("17", 69368, "G", "C"),
            ("17", 69631, "C", "T"),
        }
        assert gatk_keys == expected, f"unexpected pinned GATK call set: {sorted(gatk_keys)}"
        assert all(position not in {69067, 69803, 69807}
                   for _, position, _, _ in native_keys)

        manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = manifest_data["telemetry"]
        assert telemetry["graph_min_pruning"] == 3
        assert telemetry["graph_num_pruning_samples"] == 1
        # The oracle assertion is the byte-identical strict-pruning VCF
        # above.  HC now drops non-EventMap pileup candidates before
        # genotyping, so this implementation-detail counter need not be
        # non-zero when the graph itself has pruned the alternate.

        # GATK min-pruning=1 emits two additional 1xREF/1xALT SNPs.  This
        # enters RankSumTest's exact permutation branch rather than dropping
        # the annotations due to their total depth of two.
        min_pruning_one = "1"
        gatk_one_vcf = work / "gatk.min-pruning-one.vcf"
        native_one_vcf = work / "native.min-pruning-one.vcf"
        one_manifest = work / "native.min-pruning-one.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_one_vcf),
            "-L", region, "--min-pruning", min_pruning_one,
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_one_vcf), "--output-manifest", str(one_manifest),
            "--threads", "2", "--min-pruning", min_pruning_one,
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        gatk_one_records = records(gatk_one_vcf)
        native_one_records = records(native_one_vcf)
        assert gatk_one_vcf.read_bytes() == native_one_vcf.read_bytes(), (
            "min-pruning=1 VCF is not byte-identical with provenance disabled")
        expected_one = {
            ("17", 69067, "T", "G"),
            ("17", 69368, "G", "C"),
            ("17", 69631, "C", "T"),
            ("17", 69803, "G", "A"),
            ("17", 69807, "A", "C"),
        }
        assert {key(record) for record in native_one_records} == expected_one
        low_depth_rows = {int(record[1]): record for record in native_one_records
                          if int(record[1]) in {69803, 69807}}
        assert set(low_depth_rows) == {69803, 69807}
        for record in low_depth_rows.values():
            assert "BaseQRankSum=0.674" in record[7]
            assert "MQRankSum=0.000" in record[7]
            assert "ReadPosRankSum=0.674" in record[7]
        one_telemetry = json.loads(one_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert one_telemetry["graph_min_pruning"] == 1

        # In gVCF mode HaplotypeCaller uses its complete REF/ALT/<NON_REF> AF
        # posterior, and a retained EventMap site may remain a concrete row
        # even when the final GT is 0/0.  It blocks a reference-confidence
        # band and carries the complete Number=G PL envelope.  The GATK oracle
        # emits exactly one such row in this region (17:69657); 17:69388 is
        # monomorphic and must remain folded into a reference block.
        gatk_one_gvcf = work / "gatk.min-pruning-one.g.vcf"
        native_one_gvcf = work / "native.min-pruning-one.g.vcf"
        gvcf_manifest = work / "native.min-pruning-one.g.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_one_gvcf),
            "-L", region, "--min-pruning", min_pruning_one,
            "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_one_gvcf), "--output-manifest", str(gvcf_manifest),
            "--threads", "2", "--min-pruning", min_pruning_one,
            "--emit-ref-confidence", "GVCF",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert gatk_one_gvcf.read_bytes() == native_one_gvcf.read_bytes(), (
            "min-pruning=1 gVCF is not byte-identical with provenance disabled")
        gatk_gvcf_records = records(gatk_one_gvcf)
        native_gvcf_records = records(native_one_gvcf)
        assert len(gatk_gvcf_records) == 111
        assert len(native_gvcf_records) == len(gatk_gvcf_records)
        gvcf_by_position = {int(record[1]): record for record in native_gvcf_records}
        assert 69388 not in gvcf_by_position, (
            "monomorphic EventMap site must not become a concrete gVCF row")
        hom_ref_event = gvcf_by_position[69657]
        assert hom_ref_event[3:6] == ["T", "G,<NON_REF>", "0"]
        assert "MLEAC=0,0" in hom_ref_event[7]
        assert "MLEAF=0.00,0.00" in hom_ref_event[7]
        assert hom_ref_event[8] == "GT:AD:DP:GQ:PL:SB"
        assert hom_ref_event[9] == "0/0:19,1,0:20:26:0,26,692,57,695,726:1,18,0,1"
        # The symbolic-allele posterior is observable even at a concrete
        # alternate call: GATK rounds this gVCF QUAL to 60.31, rather than the
        # 60.32 produced by its ordinary biallelic VCF AF calculation.
        assert gvcf_by_position[69067][5] == "60.31"
        gvcf_telemetry = json.loads(gvcf_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert gvcf_telemetry["graph_min_pruning"] == 1

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "min_pruning": 3,
            "num_pruning_samples": 1,
            "calls": len(gatk_records),
            "call_set_exact": True,
            "normal_vcf_byte_identical_without_provenance": True,
            "min_pruning_one_calls": len(gatk_one_records),
            "min_pruning_one_normal_vcf_byte_identical_without_provenance": True,
            "min_pruning_one_low_depth_ranksums_exact": True,
            "min_pruning_one_gvcf_records": len(gatk_gvcf_records),
            "min_pruning_one_gvcf_byte_identical_without_provenance": True,
            "min_pruning_one_gvcf_hom_ref_event_exact": True,
            "candidate_low_support_suppressed": telemetry[
                "candidate_low_support_suppressed"],
            "scope": "normal VCF and gVCF min-pruning=1 plus strict min-pruning=3",
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
