#!/usr/bin/env python3
"""Verify an end-to-end CIGAR-derived insertion through HC VCF output."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-hc-call"
    with tempfile.TemporaryDirectory(prefix="fastgatk-indel-") as directory:
        work = Path(directory)
        reference = work / "ref.fa"
        sam = work / "reads.sam"
        vcf = work / "calls.vcf"
        manifest = work / "calls.vcf.manifest.json"
        gvcf = work / "calls.g.vcf.gz"
        gvcf_manifest = work / "calls.g.vcf.gz.manifest.json"
        full_gvcf = work / "calls.full.g.vcf.gz"
        full_gvcf_manifest = work / "calls.full.g.vcf.gz.manifest.json"
        # Keep the CIGAR insertion in an assembly-valid reference context.
        # A 32 bp all-A reference contains no unique 10/25-mers, so GATK's
        # default ReadThreadingAssembler rejects every k-mer retry and emits
        # no EventMap call.  This non-repetitive sequence is independently
        # checked against GATK: 10M1I10M produces 17:10 A>AC and 10M1D10M
        # produces 17:10 AG>A.
        reference_bases = "ACGTGCACTAGTCAGTACGATCGTACGTTAGC"
        reference.write_text(">17\n" + reference_bases + "\n", encoding="utf-8")
        sequence = reference_bases[:10] + "C" + reference_bases[10:20]
        qualities = "I" * len(sequence)
        lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:17\tLN:32",
                 "@RG\tID:rg1\tSM:SAMPLE", "@RG\tID:rg2\tSM:OTHER"]
        for index in range(4):
            lines.append(f"read{index}\t0\t17\t1\t60\t10M1I10M\t*\t0\t0\t{sequence}\t{qualities}\tRG:Z:rg1")
        # A second sample carries the same CIGAR event.  The selected-sample
        # path must exclude it before Host batching.
        lines.append(f"other\t0\t17\t1\t60\t10M1I10M\t*\t0\t0\t{sequence}\t{qualities}\tRG:Z:rg2")
        sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        result = subprocess.run([
            str(binary), "-I", str(sam), "-R", str(reference), "-L", "17:1-32",
            "-O", str(vcf), "--output-manifest", str(manifest), "--min-depth", "1",
            "--min-alt-support", "2", "--max-candidates", "16", "--sample-name", "SAMPLE",
        ], check=True, text=True, capture_output=True, env=env)
        data = json.loads(result.stdout.strip().splitlines()[-1])
        records = [line.split("\t") for line in vcf.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")]
        assert "##FORMAT=<ID=PL,Number=G" in vcf.read_text(encoding="utf-8")
        assert any(fields[3] == "A" and fields[4] == "AC" for fields in records)
        assert any("AC=" in fields[7] and "AN=2" in fields[7] and "AF=" in fields[7]
                   for fields in records)
        assert all(len(fields[9].split(":")[-1].split(",")) == 3 for fields in records)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["cigar_aware_projection"] is True
        assert metadata["compatibility"]["sample_name_selection"] is True
        assert metadata["compatibility"]["pcr_indel_model"] == "CONSERVATIVE"
        assert metadata["compatibility"]["pcr_error_rate_factor"] == 3.0
        assert metadata["telemetry"]["pairhmm_pcr_indel_model"] == "CONSERVATIVE"
        assert metadata["telemetry"]["pairhmm_pcr_error_rate_factor"] == 3.0
        assert metadata["telemetry"]["pairhmm_pcr_adjusted_positions"] > 0
        assert metadata["sample_name"] == "SAMPLE"
        assert data["reads"] == 4
        assert data["variant_calls"] >= 1
        missing_sample = subprocess.run([
            str(binary), "-I", str(sam), "-R", str(reference), "-L", "17:1-32",
            "-O", str(work / "missing.vcf"), "--sample-name", "MISSING",
        ], text=True, capture_output=True, env=env)
        assert missing_sample.returncode != 0
        assert "requested sample is absent" in missing_sample.stderr
        subprocess.run([
            str(binary), "-I", str(sam), "-R", str(reference), "-L", "17:1-32",
            "-O", str(gvcf), "-ERC", "GVCF", "--output-manifest", str(gvcf_manifest),
            "--min-depth", "1", "--min-alt-support", "2", "--max-candidates", "16",
            "--sample-name", "SAMPLE",
        ], check=True, text=True, capture_output=True, env=env)
        gvcf_text = subprocess.check_output(["gzip", "-dc", str(gvcf)], text=True)
        gvcf_records = [line.split("\t") for line in gvcf_text.splitlines()
                        if line and not line.startswith("#")]
        assert "##FORMAT=<ID=PL,Number=G" in gvcf_text
        assert "##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)" in gvcf_text
        assert "##GVCFBlock99-100=minGQ=99(inclusive),maxGQ=100(exclusive)" in gvcf_text
        for fields in gvcf_records:
            if ",<NON_REF>" in fields[4]:
                candidate_values = dict(zip(fields[8].split(":"), fields[9].split(":")))
                candidate_pl = [int(value) for value in candidate_values["PL"].split(",")]
                assert len(candidate_pl) == 6
                # GATK's candidate-site gVCF PL vector is not subject to the
                # normal VCF writer's 999 phred cap; retain only the
                # non-negative/integer contract here.
                assert all(value >= 0 for value in candidate_pl)
            elif fields[4] == "<NON_REF>":
                info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                        for item in fields[7].split(";") if "=" in item}
                sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
                assert "END" in info and "DP" in sample and "MIN_DP" in sample
                assert int(sample["MIN_DP"]) <= int(sample["DP"])
                block_ad = []
                if "AD" in sample:
                    block_ad = [int(value) for value in sample["AD"].split(",")]
                if block_ad:
                    assert len(block_ad) == 2
                block_pl = [int(value) for value in sample["PL"].split(",")]
                assert len(block_pl) == 3 and block_pl[0] == 0
                assert 0 <= int(sample["GQ"]) <= 99
        # GATK starts the zero-depth reference block immediately after this
        # insertion EventMap record, even though the source CIGAR has aligned
        # bases farther right; retain that gVCF boundary rather than assuming
        # raw BAM coverage is the block authority.
        assert any(fields[0] == "17" and fields[1] == "11" and
                   "DP" in dict(zip(fields[8].split(":"), fields[9].split(":"))) and
                   dict(zip(fields[8].split(":"), fields[9].split(":")))["DP"] == "0"
                   and "END=32" in fields[7]
                   for fields in gvcf_records)
        assert Path(f"{gvcf}.tbi").exists()
        gvcf_metadata = json.loads(gvcf_manifest.read_text(encoding="utf-8"))
        assert gvcf_metadata["telemetry"]["gvcf_interval_reference_blocks"] is True
        assert gvcf_metadata["telemetry"]["reference_block_lookup_indexed"] is True

        # Without -L, reference-backed GVCF mode now emits depth-0 blocks for
        # uncovered contig tails instead of silently dropping them.
        subprocess.run([
            str(binary), "-I", str(sam), "-R", str(reference),
            "-O", str(full_gvcf), "-ERC", "GVCF",
            "--output-manifest", str(full_gvcf_manifest), "--min-depth", "1",
            "--min-alt-support", "2", "--max-candidates", "16", "--sample-name", "SAMPLE",
        ], check=True, text=True, capture_output=True, env=env)
        full_text = subprocess.check_output(["gzip", "-dc", str(full_gvcf)], text=True)
        full_records = [line.split("\t") for line in full_text.splitlines()
                        if line and not line.startswith("#")]
        assert full_records and any(
            "END=32" in fields[7] and
            dict(zip(fields[8].split(":"), fields[9].split(":"))).get("DP") == "0"
            for fields in full_records)
        full_metadata = json.loads(full_gvcf_manifest.read_text(encoding="utf-8"))
        assert full_metadata["telemetry"]["gvcf_interval_reference_blocks"] is False
        assert full_metadata["telemetry"]["reference_block_lookup_indexed"] is True

        # A deletion candidate is left-normalized through the homopolymer and
        # occupies its full REF span in a gVCF record; reference blocks must
        # not overlap that canonical span.
        deletion_sam = work / "deletion.sam"
        deletion_gvcf = work / "deletion.g.vcf.gz"
        deletion_lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:17\tLN:32",
                          "@RG\tID:rg1\tSM:SAMPLE"]
        deletion_sequence = reference_bases[:10] + reference_bases[11:21]
        for index in range(4):
            deletion_lines.append(
                f"del{index}\t0\t17\t1\t60\t10M1D10M\t*\t0\t0\t{deletion_sequence}\t"
                f"{'I' * 20}\tRG:Z:rg1")
        deletion_sam.write_text("\n".join(deletion_lines) + "\n", encoding="utf-8")
        subprocess.run([
            str(binary), "-I", str(deletion_sam), "-R", str(reference), "-L", "17:1-32",
            "-O", str(deletion_gvcf), "-ERC", "GVCF", "--min-depth", "1",
            "--min-alt-support", "2", "--max-candidates", "16",
        ], check=True, text=True, capture_output=True, env=env)
        deletion_text = subprocess.check_output(["gzip", "-dc", str(deletion_gvcf)], text=True)
        deletion_records = [line.split("\t") for line in deletion_text.splitlines()
                            if line and not line.startswith("#")]
        deletion_candidates = [fields for fields in deletion_records
                               if fields[3] == "AG" and fields[4].startswith("A,")]
        assert deletion_candidates
        deletion_fields = deletion_candidates[0]
        info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in deletion_fields[7].split(";") if "=" in item}
        assert deletion_fields[1] == "10"
        # GATK emits a concrete gVCF variant record verbatim: INFO/END belongs to
        # reference blocks only (GVCFBlock.java sets it; GVCFBlockCombiner.submit
        # adds a variant VC unchanged).  This assertion used to require END=11 on
        # this concrete deletion candidate, i.e. it pinned native-specific
        # behaviour that diverged from GATK.  Assert the GATK contract instead.
        assert "END" not in info, f"concrete gVCF variant record must not carry INFO/END: {info}"
        for fields in deletion_records:
            if fields[4] != "<NON_REF>":
                continue
            block_start = int(fields[1])
            block_end = int({item.split("=", 1)[1] for item in fields[7].split(";")
                             if item.startswith("END=")}.pop())
            assert block_end < 10 or block_start > 11

        # Complex multiallelic calling is covered by the strict real-BAM GATK
        # oracle in verify_hc_complex_multiallelic_oracle.py.  The former
        # 32bp, 21M synthetic trace does not construct an EventMap in GATK,
        # so treating it as an expected call would test non-GATK behaviour.

        # Soft-clip semantics are covered by dedicated real-BAM GATK-oracle
        # tests: fastgatk-hc-softclip-gatk-oracle and
        # fastgatk-hc-softclip-contig-start-gatk-oracle.  The former
        # micro-fixture was not a GATK-reproducible calling scenario.

        # All GATK PCR model enum values share the same Host preparation path;
        # exercise a non-default rate factor so the CLI does not regress to the
        # former CONSERVATIVE/NONE-only parser.
        hostile_vcf = work / "hostile.vcf"
        hostile_manifest = work / "hostile.manifest.json"
        subprocess.run([
            str(binary), "-I", str(sam), "-R", str(reference), "-L", "17:1-32",
            "-O", str(hostile_vcf), "--output-manifest", str(hostile_manifest),
            "--min-depth", "1", "--min-alt-support", "2", "--max-candidates", "16",
            "--sample-name", "SAMPLE", "--pcr-indel-model", "HOSTILE",
        ], check=True, text=True, capture_output=True, env=env)
        hostile_metadata = json.loads(hostile_manifest.read_text(encoding="utf-8"))
        assert hostile_metadata["compatibility"]["pcr_indel_model"] == "HOSTILE"
        assert hostile_metadata["compatibility"]["pcr_error_rate_factor"] == 1.0
        assert hostile_metadata["telemetry"]["pairhmm_default_indel_quality"] == 40
        assert hostile_metadata["telemetry"]["pairhmm_pcr_adjusted_positions"] > 0
        print(json.dumps({"status": "pass", "variant_calls": data["variant_calls"],
                          "records": len(records), "gvcf_records": len(gvcf_records),
                          "deletion_records": len(deletion_records),
                          "softclip_coverage":
                              "fastgatk-hc-softclip-gatk-oracle"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
