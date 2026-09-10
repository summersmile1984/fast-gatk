#!/usr/bin/env python3
"""Verify deterministic GVCF block merging and allele compaction."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=TREE_SCORE,Number=1,Type=Float,Description=Tree score>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-reblock-gvcf")))
    java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-gvcf-") as temp:
        work = Path(temp)
        source = work / "input.g.vcf.gz"
        output = work / "output.g.vcf.gz"
        manifest = work / "output.manifest.json"
        body = (
            "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT:DP:AD:PL:GQ\t0/0:10:10,0:0,10,200:10\n"
            "chr1\t6\t.\tC\t<NON_REF>\t.\tPASS\tEND=10\tGT:DP:AD:PL:GQ\t0/0:8:8,0:0,15,200:15\n"
            "chr1\t12\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL:GQ\t0/1:20:12,8,0,0:50,0,80,99,99,99,99,99,99,20:.\n"
            # PL[0], rather than the pre-existing GQ field, controls the
            # --rgq-threshold-to-no-call decision.
            "chr1\t20\t.\tC\tT,<NON_REF>\t.\tPASS\tDP=5\tGT:DP:AD:PL:GQ\t0/1:5:4,1,0:2,0,40,99,99,99:99\n"
        )
        with gzip.open(source, "wt", encoding="utf-8") as stream:
            stream.write(HEADER + body)
        result = subprocess.run(
            [str(binary), "-V", str(source), "-GQB", "20", "-GQB", "100",
             "--drop-low-quals", "--rgq-threshold", "10", "--floor-blocks", "-O", str(output),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["merged_blocks"] == 1
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert len(records) == 3
        assert records[0][7] == "END=10"
        format_names = records[0][8].split(":")
        format_values = records[0][9].split(":")
        assert format_values[format_names.index("GQ")] == "0"
        # --floor-blocks is also a storage optimization in GATK: reference
        # blocks retain GT/DP/GQ but drop PL and MIN_DP.  Variant records keep
        # their likelihood arrays.
        assert "PL" not in format_names
        assert "MIN_DP" not in format_names
        # The multiallelic site is compacted to the called G allele plus NON_REF.
        assert records[1][4] == "G,<NON_REF>"
        assert len(records[1][9].split(":")[3].split(",")) == 6
        # Missing GQ is materialized from the complete PL vector after
        # compaction (the compacted PL has min/second = 0/20).
        assert records[1][9].split(":")[-1] == "20"
        # The low-quality site was converted into a hom-ref reference block.
        assert records[2][4] == "<NON_REF>"
        assert "END=20" in records[2][7]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["gq_bands"] is True
        assert metadata["compatibility"]["reference_block_merge"] is True
        assert metadata["compatibility"]["rgq_threshold_uses_pl0"] is True
        assert metadata["compatibility"]["min_dp_merge"] is True
        assert metadata["telemetry"]["pl_remap_kernel_calls"] == 1
        assert metadata["telemetry"]["allele_field_remap_kernel_calls"] == 1
        assert metadata["telemetry"]["allele_field_remap_kernel_execution_space"] in {"OpenMP", "Serial"}
        assert metadata["telemetry"]["pl_remap_kernel_execution_space"] in {"OpenMP", "Serial"}
        assert metadata["telemetry"]["pl_remap_kernel_seconds"] >= 0.0
        assert all(item["complete"] for item in metadata["outputs"])

        # Plain GVCF output is indexed with the shared Tribble writer and can
        # be queried by GATK directly, preserving the drop-in file contract.
        plain_output = work / "output-plain.g.vcf"
        plain_manifest = work / "output-plain.manifest.json"
        plain_result = subprocess.run([
            str(binary), "-V", str(source), "-GQB", "20", "-GQB", "100",
            "--drop-low-quals", "--rgq-threshold", "10", "--floor-blocks",
            "-O", str(plain_output), "--output-manifest", str(plain_manifest),
        ], text=True, capture_output=True, check=False)
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.exists() and Path(f"{plain_output}.idx").exists()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "output-plain-query.vcf"
        plain_query_result = subprocess.run([
            str(java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
            "-L", "chr1:12-12", "-O", str(plain_query),
        ], text=True, capture_output=True, check=False)
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True

        # The output-index switch follows GATK's optional-boolean contract:
        # an explicit false (whether separated or inline) suppresses the TBI
        # sidecar and records the choice in the manifest, while a bare flag
        # enables it.  Invalid literals must fail instead of being coerced to
        # true, which would make Nextflow output declarations nondeterministic.
        no_index_output = work / "no-index.g.vcf.gz"
        no_index_manifest = work / "no-index.manifest.json"
        no_index_result = subprocess.run(
            [str(binary), "-V", str(source), "--create-output-variant-index", "false",
             "-O", str(no_index_output), "--output-manifest", str(no_index_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert no_index_result.returncode == 0, no_index_result.stderr
        no_index_summary = json.loads(no_index_result.stdout.splitlines()[-1])
        assert no_index_summary["vcf_index"] is False
        assert no_index_output.exists() and not Path(f"{no_index_output}.tbi").exists()
        no_index_metadata = json.loads(no_index_manifest.read_text(encoding="utf-8"))
        assert no_index_metadata["compatibility"]["vcf_index"] is False
        assert not any(item["kind"] == "vcf-index" for item in no_index_metadata["outputs"])
        bare_index_output = work / "bare-index.g.vcf.gz"
        bare_index_result = subprocess.run(
            [str(binary), "-V", str(source), "--create-output-variant-index",
             "-O", str(bare_index_output)], text=True, capture_output=True, check=False,
        )
        assert bare_index_result.returncode == 0, bare_index_result.stderr
        bare_index_summary = json.loads(bare_index_result.stdout.splitlines()[-1])
        assert bare_index_summary["vcf_index"] is True
        assert Path(f"{bare_index_output}.tbi").exists()
        invalid_index_result = subprocess.run(
            [str(binary), "-V", str(source), "--create-output-variant-index=maybe",
             "-O", str(work / "invalid-index.g.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert invalid_index_result.returncode != 0
        assert "invalid boolean" in invalid_index_result.stderr

        # TREE_SCORE below the configured threshold (including a missing
        # value, which GATK treats as zero) converts a variant to a GQ0
        # hom-ref block.  The threshold option must not be a parser-only
        # compatibility shim.
        tree_output = work / "tree-score.g.vcf.gz"
        tree_result = subprocess.run(
            [str(binary), "-V", str(source), "--tree-score-threshold-to-no-call", "0.5",
             "-O", str(tree_output)], text=True, capture_output=True, check=False,
        )
        assert tree_result.returncode == 0, tree_result.stderr
        tree_summary = json.loads(tree_result.stdout.splitlines()[-1])
        assert tree_summary["status"] == "prototype"
        assert tree_summary["converted_tree_score"] >= 1
        tree_records = [line.split("\t") for line in gzip.open(
            tree_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert any(record[1] == "12" and record[4] == "<NON_REF>" for record in tree_records)

        # QUALapprox is derived from the post-compaction PL[0] and is exposed
        # only when requested, matching ReblockGVCF's optional annotation
        # contract.
        approx_output = work / "qual-approx.g.vcf.gz"
        approx_result = subprocess.run(
            [str(binary), "-V", str(source), "--do-qual-score-approximation",
             "-O", str(approx_output)],
            text=True, capture_output=True, check=False,
        )
        assert approx_result.returncode == 0, approx_result.stderr
        approx_records = [line.split("\t") for line in gzip.open(
            approx_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        approx_variant = next(record for record in approx_records if record[1] == "12")
        assert "QUALapprox=50" in approx_variant[7]
        assert "AS_QUALapprox=|50|0" in approx_variant[7]
        assert "VarDP=20" in approx_variant[7]
        assert "AS_VarDP=12|8|0" in approx_variant[7]
        assert "RAW_GT_COUNT=0,1,0" in approx_variant[7]
        approx_header = gzip.open(approx_output, "rt", encoding="utf-8").read()
        assert "##INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer" in approx_header
        approx_manifest = json.loads((approx_output.with_name(approx_output.name + ".manifest.json")).read_text(encoding="utf-8"))
        assert approx_manifest["compatibility"]["raw_genotype_count"] is True
        assert approx_manifest["telemetry"]["raw_genotype_count_annotations"] == 2
        assert approx_manifest["compatibility"]["qual_approx_depth"] is True
        assert approx_manifest["telemetry"]["qual_approx_depth_annotations"] == 2

        # Repeated interval selectors use GATK's union semantics and retain
        # a gVCF block when its END span overlaps the requested interval.
        interval_output = work / "interval.g.vcf.gz"
        interval_manifest = work / "interval.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-V", str(source), "-L", "chr1:12-12",
             "--intervals", "chr1:20-20", "-O", str(interval_output),
             "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert interval_result.returncode == 0, interval_result.stderr
        interval_summary = json.loads(interval_result.stdout.splitlines()[-1])
        assert interval_summary["interval_skipped"] == 2
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["12", "20"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["intervals"] == 2

        # INTERSECTION is a true set operation across repeated selectors;
        # only the locus at 12 survives the overlap of 1..12 and 12..12.
        intersection_output = work / "interval-intersection.g.vcf.gz"
        intersection_manifest = work / "interval-intersection.manifest.json"
        intersection_result = subprocess.run(
            [str(binary), "-V", str(source), "-L", "chr1:1-12",
             "-L", "chr1:12-12", "--interval-set-rule", "INTERSECTION",
             "-O", str(intersection_output), "--output-manifest", str(intersection_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_summary = json.loads(intersection_result.stdout.splitlines()[-1])
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        intersection_records = [line.split("\t") for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in intersection_records] == ["12"]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        invalid_rule = subprocess.run(
            [str(binary), "-V", str(source), "-L", "chr1:1-2",
             "--interval-set-rule=NOT_A_RULE", "-O", str(work / "invalid-rule.g.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert invalid_rule.returncode != 0 and "interval-set-rule" in invalid_rule.stderr
        interval_list = work / "reblock.interval_list"
        interval_list.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
                                 "chr1\t12\t12\t+\tfirst\nchr1\t20\t20\t+\tsecond\n",
                                 encoding="utf-8")
        interval_file_output = work / "interval-file.g.vcf.gz"
        interval_file_manifest = work / "interval-file.manifest.json"
        interval_file_result = subprocess.run(
            [str(binary), "-V", str(source), "-L", str(interval_list),
             "-O", str(interval_file_output), "--output-manifest", str(interval_file_manifest)],
            text=True, capture_output=True, check=False)
        assert interval_file_result.returncode == 0, interval_file_result.stderr
        interval_file_metadata = json.loads(interval_file_manifest.read_text(encoding="utf-8"))
        assert interval_file_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_file_metadata["telemetry"]["interval_list_records"] == 2
        # Reblocking must use the VCF Number=G rank for the sample ploidy,
        # not the diploid triangular shortcut.  With three alleles and
        # triploid GT, the old 20-entry PL is compacted to 10 entries after
        # dropping the unused concrete ALT while retaining <NON_REF>.
        triploid_source = work / "triploid.g.vcf.gz"
        triploid_output = work / "triploid.reblocked.g.vcf.gz"
        with gzip.open(triploid_source, "wt", encoding="utf-8") as stream:
            stream.write(HEADER)
            stream.write(
                "chr1\t30\t.\tA\tC,G,<NON_REF>\t.\tPASS\tDP=17\t"
                "GT:DP:AD:PL:GQ\t0/1/1:17:12,5,0,0:"
                "0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190:99\n"
            )
        triploid_result = subprocess.run(
            [str(binary), "-V", str(triploid_source), "-O", str(triploid_output)],
            text=True, capture_output=True, check=False,
        )
        assert triploid_result.returncode == 0, triploid_result.stderr
        triploid_records = [line.split("\t") for line in gzip.open(
            triploid_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(triploid_records) == 1
        assert triploid_records[0][4] == "C,<NON_REF>"
        triploid_values = triploid_records[0][9].split(":")
        assert triploid_values[0] == "0/1/1"
        assert triploid_values[2] == "12,5,0"
        assert triploid_values[3] == "0,10,20,30,100,110,120,160,170,190"

        # Multi-sample records use sample-major Kokkos remap for Number=R AD
        # and Number=G PL.  Reference-block merging must keep each sample's
        # MIN_DP/GQ band independent rather than collapsing to sample 0.
        multisample_source = work / "multisample.g.vcf.gz"
        multisample_output = work / "multisample.reblocked.g.vcf.gz"
        multisample_header = HEADER.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
        with gzip.open(multisample_source, "wt", encoding="utf-8") as stream:
            stream.write(multisample_header)
            stream.write(
                "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\t"
                "GT:DP:AD:PL:GQ\t0/0:10:10,0:0,10,200:10\t"
                "0/0:8:8,0:0,15,200:15\n"
                "chr1\t6\t.\tC\t<NON_REF>\t.\tPASS\tEND=10\t"
                "GT:DP:AD:PL:GQ\t0/0:8:8,0:0,15,200:15\t"
                "0/0:7:7,0:0,12,200:12\n"
                "chr1\t12\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=20\t"
                "GT:DP:AD:PL:GQ\t0/1:20:12,8,0,0:50,0,80,99,99,99:20\t"
                "0/1:18:10,6,0,0:60,0,90,99,99,99:30\n"
            )
        multisample_result = subprocess.run(
            [str(binary), "-V", str(multisample_source), "--do-qual-approx",
             "-O", str(multisample_output)], text=True, capture_output=True, check=False,
        )
        assert multisample_result.returncode == 0, multisample_result.stderr
        multisample_records = [line.split("\t") for line in gzip.open(
            multisample_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(multisample_records) == 2
        block = multisample_records[0]
        block_names = block[8].split(":")
        block_s1 = block[9].split(":")
        block_s2 = block[10].split(":")
        assert block[7] == "END=10"
        assert block_s1[block_names.index("MIN_DP")] == "8"
        assert block_s2[block_names.index("MIN_DP")] == "7"
        # Existing hom-ref blocks preserve the element-wise minimum PL/GQ
        # rather than being rewritten to synthetic GQ0 calls.  This mirrors
        # GVCFBlock's merge semantics in GATK.
        assert block_s1[block_names.index("GQ")] == "10"
        assert block_s2[block_names.index("GQ")] == "12"
        assert block_s1[block_names.index("DP")] == "9"
        assert block_s2[block_names.index("DP")] == "8"
        variant = multisample_records[1]
        assert variant[4] == "G,<NON_REF>"
        variant_names = variant[8].split(":")
        variant_s1 = variant[9].split(":")
        variant_s2 = variant[10].split(":")
        assert variant_s1[variant_names.index("AD")] == "12,8,0"
        assert variant_s2[variant_names.index("AD")] == "10,6,0"
        assert len(variant_s1[variant_names.index("PL")].split(",")) == 6
        assert len(variant_s2[variant_names.index("PL")].split(",")) == 6
        assert "QUALapprox=110" in variant[7]
        assert "AS_QUALapprox=|110|0" in variant[7]
        assert "RAW_GT_COUNT=0,2,0" in variant[7]
        multisample_manifest = json.loads((multisample_output.with_name(
            multisample_output.name + ".manifest.json")).read_text(encoding="utf-8"))
        assert multisample_manifest["telemetry"]["sample_count"] == 2
        assert multisample_manifest["telemetry"]["raw_genotype_count_annotations"] == 1

        # Low-quality deletion calls are converted to a reference block with
        # the REF trimmed to its leading base while END retains the original
        # deletion span, matching ReblockGVCF's changeCallToHomRefVersusNonRef.
        deletion_source = work / "deletion.g.vcf.gz"
        deletion_output = work / "deletion.reblocked.g.vcf.gz"
        with gzip.open(deletion_source, "wt", encoding="utf-8") as stream:
            stream.write(HEADER)
            stream.write(
                "chr1\t25\t.\tAT\tA,<NON_REF>\t.\tPASS\tEND=26\t"
                "GT:DP:AD:PL:GQ\t0/1:10:6,4,0:1,0,20,99,99,99:99\n"
            )
        deletion_result = subprocess.run(
            [str(binary), "-V", str(deletion_source), "--rgq-threshold", "10",
             "-O", str(deletion_output)],
            text=True, capture_output=True, check=False,
        )
        assert deletion_result.returncode == 0, deletion_result.stderr
        deletion_records = [line.split("\t") for line in gzip.open(
            deletion_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(deletion_records) == 1
        assert deletion_records[0][1] == "25"
        assert deletion_records[0][3] == "A"
        assert deletion_records[0][4] == "<NON_REF>"
        assert "END=26" in deletion_records[0][7]
        # GATK's changeCallToHomRefVersusNonRef() clears AD/noAttributes when
        # a low-quality deletion is converted to a reference block.  Keep the
        # complete scalar/vector contract explicit so stale allele depths do
        # not silently leak into downstream CombineGVCFs/GenotypeGVCFs.
        deletion_names = deletion_records[0][8].split(":")
        deletion_values = deletion_records[0][9].split(":")
        assert "AD" not in deletion_names
        assert deletion_values[deletion_names.index("GT")] == "0/0"
        assert deletion_values[deletion_names.index("PL")] == "0,0,0"
        assert deletion_values[deletion_names.index("GQ")] == "0"
        assert deletion_values[deletion_names.index("DP")] == "10"
        assert deletion_values[deletion_names.index("MIN_DP")] == "10"

        # After an unused concrete ALT is removed, GATK reverse-trims a
        # common suffix from the remaining high-quality alleles.  The
        # symbolic <NON_REF> allele is preserved and Number=G/AD fields keep
        # their remapped ordering.
        trim_source = work / "trim.g.vcf.gz"
        trim_output = work / "trim.reblocked.g.vcf.gz"
        trim_manifest = work / "trim.manifest.json"
        trim_reference = work / "trim.fa"
        trim_reference.write_text(">chr1\n" + ("A" * 100) + "\n", encoding="utf-8")
        # FASTA layout is one 100-base line: header is six bytes, sequence
        # starts at offset six, line bases=100 and line bytes=101.
        trim_reference.with_suffix(".fa.fai").write_text(
            "chr1\t100\t6\t100\t101\n", encoding="utf-8")
        with gzip.open(trim_source, "wt", encoding="utf-8") as stream:
            stream.write(HEADER)
            stream.write(
                "chr1\t50\t.\tAC\tGC,TC,<NON_REF>\t.\tPASS\tDP=15\t"
                "GT:DP:AD:PL:GQ\t0/1:19:10,5,0,4:"
                "0,10,20,30,40,50,60,70,80,90:30\n"
            )
        trim_result = subprocess.run(
            [str(binary), "-V", str(trim_source), "-O", str(trim_output),
             "--output-manifest", str(trim_manifest), "-R", str(trim_reference)],
            text=True, capture_output=True, check=False,
        )
        assert trim_result.returncode == 0, trim_result.stderr
        trim_records = [line.split("\t") for line in gzip.open(
            trim_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(trim_records) == 2
        assert trim_records[0][1] == "50"
        assert trim_records[0][3] == "A"
        assert trim_records[0][4] == "G,<NON_REF>"
        trim_values = trim_records[0][9].split(":")
        assert trim_values[0] == "0/1"
        assert trim_values[1] == "15"
        assert trim_values[2] == "10,5,0"
        assert trim_values[3] == "0,10,20,60,70,90"
        gap_record = trim_records[1]
        assert gap_record[1] == "51"
        assert gap_record[3] == "A"
        assert gap_record[4] == "<NON_REF>"
        assert "END=51" in gap_record[7]
        gap_names = gap_record[8].split(":")
        gap_values = gap_record[9].split(":")
        assert gap_values[gap_names.index("DP")] == "15"
        assert gap_values[gap_names.index("PL")] == "0,30,50"
        assert gap_values[gap_names.index("GQ")] == "30"
        trim_metadata = json.loads(trim_manifest.read_text(encoding="utf-8"))
        assert trim_metadata["compatibility"]["reverse_allele_trimming"] is True
        assert trim_metadata["compatibility"]["deletion_gap_ref_block"] is True
        assert trim_metadata["compatibility"]["non_ref_ad_cleanup"] is True
        assert trim_metadata["telemetry"]["reverse_allele_trimmed_variants"] == 1
        assert trim_metadata["telemetry"]["deletion_gap_ref_blocks"] == 1
        assert trim_metadata["telemetry"]["non_ref_ad_zeroed"] == 1
        assert trim_metadata["telemetry"]["non_ref_ad_kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert trim_metadata["telemetry"]["non_ref_ad_kernel_execution_policy"] == "RangePolicy"
        assert trim_metadata["telemetry"]["non_ref_ad_kernel_calls"] == 1
        assert trim_metadata["telemetry"]["non_ref_ad_kernel_prepare_seconds"] >= 0.0
        assert trim_metadata["telemetry"]["non_ref_ad_kernel_execute_seconds"] >= 0.0

        # GATK removes stale single-sample INFO annotations and supports an
        # explicit keep/remove contract for caller-specific fields.
        annotation_source = work / "annotations.g.vcf.gz"
        annotation_output = work / "annotations.reblocked.g.vcf.gz"
        annotation_manifest = work / "annotations.manifest.json"
        annotation_header = HEADER.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=stale>\n"
            "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=stale>\n"
            "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=stale>\n"
            "##INFO=<ID=GVCFBlock10-20,Number=0,Type=Flag,Description=stale>\n"
            "##INFO=<ID=KEEP_ME,Number=1,Type=Integer,Description=caller>\n").replace(
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
            "##FORMAT=<ID=ZZ,Number=1,Type=String,Description=remove-me>\n")
        with gzip.open(annotation_source, "wt", encoding="utf-8") as stream:
            stream.write(annotation_header)
            stream.write(
                "chr1\t60\t.\tA\tG,<NON_REF>\t.\tPASS\t"
                "DP=10;HaplotypeScore=3.0;InbreedingCoeff=0.1;MLEAC=1;"
                "GVCFBlock10-20;KEEP_ME=7\t"
                "GT:DP:AD:PL:GQ:ZZ\t0/1:10:6,4,0:20,0,40,99,99,99:20:drop\n"
            )
        annotation_result = subprocess.run(
            [str(binary), "-V", str(annotation_source), "-O", str(annotation_output),
             "--annotations-to-keep", "KEEP_ME", "--format-annotations-to-remove", "ZZ",
             "--output-manifest", str(annotation_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert annotation_result.returncode == 0, annotation_result.stderr
        annotation_text = gzip.open(annotation_output, "rt", encoding="utf-8").read()
        annotation_record = next(line.split("\t") for line in annotation_text.splitlines()
                                 if line and not line.startswith("#"))
        assert "KEEP_ME=7" in annotation_record[7]
        assert "HaplotypeScore" not in annotation_record[7]
        assert "InbreedingCoeff" not in annotation_record[7]
        assert "MLEAC" not in annotation_record[7]
        assert "GVCFBlock10-20" not in annotation_record[7]
        assert "ZZ" not in annotation_record[8]
        assert "KEEP_ME" in annotation_text
        assert "HaplotypeScore" not in annotation_text
        assert "GVCFBlock10-20" not in annotation_text
        annotation_metadata = json.loads(annotation_manifest.read_text(encoding="utf-8"))
        assert annotation_metadata["compatibility"]["annotation_cleanup"] is True
        assert annotation_metadata["compatibility"]["format_annotation_removal"] is True
        assert annotation_metadata["telemetry"]["info_annotations_removed"] == 4
        assert annotation_metadata["telemetry"]["format_annotations_removed"] == 1
        # The pre-existing native spelling remains a compatibility alias for
        # callers that adopted it before the GATK formal option was added.
        legacy_annotation_output = work / "annotations-legacy-alias.g.vcf.gz"
        legacy_annotation_result = subprocess.run(
            [str(binary), "-V", str(annotation_source), "-O", str(legacy_annotation_output),
             "--annotations-to-remove", "ZZ"],
            text=True, capture_output=True, check=False,
        )
        assert legacy_annotation_result.returncode == 0, legacy_annotation_result.stderr
        legacy_annotation_text = gzip.open(
            legacy_annotation_output, "rt", encoding="utf-8").read()
        assert "ZZ" not in legacy_annotation_text
        # GATK rejects unknown annotations and Number=A annotations in
        # --annotations-to-keep because allele compaction would make them
        # stale.  Native must fail closed before publishing an output.
        keep_allele_specific = subprocess.run(
            [str(binary), "-V", str(annotation_source), "--annotations-to-keep", "MLEAC",
             "-O", str(work / "invalid-allele-annotation.g.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert keep_allele_specific.returncode != 0
        assert "allele-specific" in keep_allele_specific.stderr
        keep_unknown = subprocess.run(
            [str(binary), "-V", str(annotation_source), "--annotations-to-keep", "MISSING_TAG",
             "-O", str(work / "invalid-unknown-annotation.g.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert keep_unknown.returncode != 0
        assert "INFO annotation not found" in keep_unknown.stderr

        # GATK's updateMQAnnotations() emits the modern RAW_MQandDP tuple for
        # high-quality variants.  A legacy RAW_MQ value is retained and also
        # gets the deprecated MQ_DP depth field for downstream compatibility.
        mq_source = work / "mq-annotations.g.vcf.gz"
        mq_output = work / "mq-annotations.reblocked.g.vcf.gz"
        mq_header = HEADER.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=MQ,Number=1,Type=Float,Description=RMS mapping quality>\n"
            "##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description=legacy raw MQ>\n")
        with gzip.open(mq_source, "wt", encoding="utf-8") as stream:
            stream.write(mq_header)
            stream.write(
                "chr1\t70\t.\tA\tG,<NON_REF>\t.\tPASS\tMQ=50;DP=10\t"
                "GT:DP:AD:PL:GQ\t0/1:10:6,4,0:20,0,40,99,99,99:20\n"
                "chr1\t80\t.\tA\tG,<NON_REF>\t.\tPASS\tRAW_MQ=123.4;DP=10\t"
                "GT:DP:AD:PL:GQ\t0/1:10:6,4,0:20,0,40,99,99,99:20\n"
            )
        mq_result = subprocess.run(
            [str(binary), "-V", str(mq_source), "-O", str(mq_output)],
            text=True, capture_output=True, check=False,
        )
        assert mq_result.returncode == 0, mq_result.stderr
        mq_text = gzip.open(mq_output, "rt", encoding="utf-8").read()
        mq_records = [line.split("\t") for line in mq_text.splitlines()
                      if line and not line.startswith("#")]
        assert len(mq_records) == 2
        assert "RAW_MQandDP=25000,10" in mq_records[0][7]
        assert "RAW_MQandDP=123,10" in mq_records[1][7]
        assert "MQ_DP=10" in mq_records[1][7]
        assert "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer" in mq_text
        assert "##INFO=<ID=MQ_DP,Number=1,Type=Integer" in mq_text

        # Site FILTER handling follows ReblockGVCF's keep/add options.
        filter_source = work / "filters.g.vcf.gz"
        filter_header = HEADER.replace(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "##FILTER=<ID=LowQual,Description=low quality site>\n")
        with gzip.open(filter_source, "wt", encoding="utf-8") as stream:
            stream.write(filter_header)
            stream.write(
                "chr1\t40\t.\tA\tG,<NON_REF>\t.\tLowQual\tDP=20\t"
                "GT:DP:AD:PL:GQ\t0/1:20:12,8,0:50,0,80,99,99,99:50\n"
            )
        filter_default_output = work / "filters-default.g.vcf.gz"
        filter_default = subprocess.run(
            [str(binary), "-V", str(filter_source), "-O", str(filter_default_output)],
            text=True, capture_output=True, check=False,
        )
        assert filter_default.returncode == 0, filter_default.stderr
        filter_default_record = next(line.split("\t") for line in gzip.open(
            filter_default_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#"))
        assert filter_default_record[6] in {"PASS", "."}
        filter_keep_output = work / "filters-keep.g.vcf.gz"
        filter_keep_manifest = work / "filters-keep.manifest.json"
        filter_keep = subprocess.run(
            [str(binary), "-V", str(filter_source), "--keep-site-filters",
             "--add-site-filters-to-genotype", "-O", str(filter_keep_output),
             "--output-manifest", str(filter_keep_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert filter_keep.returncode == 0, filter_keep.stderr
        filter_keep_record = next(line.split("\t") for line in gzip.open(
            filter_keep_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#"))
        assert filter_keep_record[6] == "LowQual"
        filter_keep_names = filter_keep_record[8].split(":")
        filter_keep_values = filter_keep_record[9].split(":")
        assert filter_keep_values[filter_keep_names.index("FT")] == "LowQual"
        filter_keep_metadata = json.loads(filter_keep_manifest.read_text(encoding="utf-8"))
        assert filter_keep_metadata["compatibility"]["keep_site_filters"] is True
        assert filter_keep_metadata["compatibility"]["add_site_filters_to_genotype"] is True
    print(json.dumps({"status": "pass", "output_records": summary["output_records"],
                      "merged_blocks": summary["merged_blocks"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
