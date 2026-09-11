#!/usr/bin/env python3
"""Verify the native VariantsToTable field/row compatibility contract."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##FILTER=<ID=LowQual,Description=Low quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tZETA\tALPHA
"""


def parse_table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-variants-to-table"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-variants-to-table-") as directory:
        work = Path(directory)
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        source = work / "input.vcf"
        source.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG,T\t50\tPASS\tDP=20;AC=2,3;AF=0.200000000,3e-1\tGT:AD:PL:DP\t0/1:8,4,2:50,0,60,80,90,100:14\t0/2:8,1,5:70,80,90,60,70,0:14\n"
            + "chr1\t2\t.\tC\tT\t10\tLowQual\tDP=5;AC=1;AF=1e-1\tGT:AD:PL:DP\t0/0:5,0:0,20,60:5\t./.:.,.:.,.,.:0\n",
            encoding="utf-8",
        )

        table = work / "table.tsv"
        manifest = work / "table.manifest.json"
        result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(table),
            "-F", "CHROM", "-F", "POS", "-F", "ALT", "-F", "MULTI-ALLELIC",
            "-F", "AC", "-F", "QUAL", "-F", "TYPE", "-F", "TRANSITION",
            "-F", "EVENTLENGTH", "-F", "HET", "-F", "HOM-REF", "-F", "HOM-VAR",
            "-F", "NO-CALL", "-F", "VAR", "-F", "NSAMPLES", "-F", "NCALLED",
            "-F", "SAMPLE_NAME", "-GF", "GT", "-GF", "AD", "-GF", "DP",
            "--output-manifest", str(manifest),
        ], text=True, capture_output=True, check=False)
        assert result.returncode == 0, result.stderr
        header, rows = parse_table(table)
        assert header == ["CHROM", "POS", "ALT", "MULTI-ALLELIC", "AC", "QUAL", "TYPE",
                          "TRANSITION", "EVENTLENGTH", "HET", "HOM-REF", "HOM-VAR",
                          "NO-CALL", "VAR", "NSAMPLES", "NCALLED", "SAMPLE_NAME",
                          "ALPHA.GT", "ALPHA.AD", "ALPHA.DP", "ZETA.GT", "ZETA.AD", "ZETA.DP"]
        assert len(rows) == 1  # filtered record is ignored by default
        assert rows[0][2] == "G,T" and rows[0][4] == "2,3"
        assert rows[0][5:17] == ["50.0", "SNP", "-1", "0", "2", "0", "0", "0", "2", "2", "2", "ZETA"]
        assert rows[0][17] == "A/T" and rows[0][20] == "A/G"
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["field_extraction"] is True
        assert metadata["compatibility"]["sample_name_sort"] is True
        assert metadata["compatibility"]["raw_info_tokens"] is True
        assert metadata["telemetry"]["skipped_filtered"] == 1

        # HTSlib stores VCF Float values as IEEE-754 float32, while GATK keeps
        # the original text INFO token.  The native Host reader retains the
        # token so long decimals and exponent spelling are not lost.
        float_table = work / "float.tsv"
        float_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(float_table),
            "-F", "AF", "--show-filtered",
        ], text=True, capture_output=True, check=False)
        assert float_result.returncode == 0, float_result.stderr
        float_oracle = work / "float-oracle.tsv"
        float_gatk_result = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
            "-V", str(source), "-O", str(float_oracle), "-F", "AF", "--show-filtered",
        ], text=True, capture_output=True, check=False)
        assert float_gatk_result.returncode == 0, float_gatk_result.stderr
        assert float_table.read_text(encoding="utf-8") == float_oracle.read_text(encoding="utf-8")

        # Literal missing INFO values and flag attributes are also lexical
        # fields in htsjdk's VariantContext table path.
        edge_source = work / "info-edge.vcf"
        edge_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "##INFO=<ID=FLAG,Number=0,Type=Flag,Description=Flag>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\tG\t50\tPASS\tAF=.;FLAG\n"
            "chr1\t2\t.\tA\tG\t50\tPASS\tFLAG\n"
            "chr1\t3\t.\tA\tG\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        edge_table = work / "info-edge.tsv"
        edge_result = subprocess.run([
            str(binary), "-V", str(edge_source), "-O", str(edge_table),
            "-F", "AF", "-F", "FLAG", "--show-filtered",
        ], text=True, capture_output=True, check=False)
        assert edge_result.returncode == 0, edge_result.stderr
        edge_oracle = work / "info-edge-oracle.tsv"
        edge_gatk_result = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
            "-V", str(edge_source), "-O", str(edge_oracle),
            "-F", "AF", "-F", "FLAG", "--show-filtered",
        ], text=True, capture_output=True, check=False)
        assert edge_gatk_result.returncode == 0, edge_gatk_result.stderr
        assert edge_table.read_text(encoding="utf-8") == edge_oracle.read_text(encoding="utf-8")

        raw_table = work / "raw.tsv"
        raw_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(raw_table),
            "-F", "CHROM", "-F", "FILTER", "--show-filtered",
        ], text=True, capture_output=True, check=False)
        assert raw_result.returncode == 0, raw_result.stderr
        _, raw_rows = parse_table(raw_table)
        assert [row[1] for row in raw_rows] == ["PASS", "LowQual"]

        excluded_table = work / "excluded.tsv"
        excluded_manifest = work / "excluded.manifest.json"
        excluded_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(excluded_table),
            "-F", "POS", "--show-filtered", "-XL", "chr1:1-1",
            "--output-manifest", str(excluded_manifest),
        ], text=True, capture_output=True, check=False)
        assert excluded_result.returncode == 0, excluded_result.stderr
        _, excluded_rows = parse_table(excluded_table)
        assert excluded_rows == [["2"]]
        excluded_metadata = json.loads(excluded_manifest.read_text(encoding="utf-8"))
        assert excluded_metadata["compatibility"]["interval_exclusion"] is True
        assert excluded_metadata["telemetry"]["skipped_excluded"] == 1
        intersection_table = work / "intersection.tsv"
        intersection_manifest = work / "intersection.manifest.json"
        intersection_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(intersection_table), "-F", "POS",
            "--show-filtered", "-L", "chr1:1-2", "-L", "chr1:2-2",
            "--interval-set-rule", "INTERSECTION",
            "--output-manifest", str(intersection_manifest),
        ], text=True, capture_output=True, check=False)
        assert intersection_result.returncode == 0, intersection_result.stderr
        _, intersection_rows = parse_table(intersection_table)
        assert intersection_rows == [["2"]]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"

        # GATK exposes --variant-output-filtering through the shared interval
        # argument collection.  VariantsToTable writes a table directly (it
        # does not install IntervalFilteringVcfWriter), so the pinned Java
        # implementation validates the mode but keeps normal OVERLAPS
        # traversal.  Exercise every accepted mode and retain the selected
        # value in the manifest so an adapter cannot silently drop it.
        for mode in ("STARTS_IN", "ENDS_IN", "OVERLAPS", "CONTAINED", "ANYWHERE"):
            mode_table = work / f"variant-output-{mode}.tsv"
            mode_manifest = work / f"variant-output-{mode}.manifest.json"
            mode_result = subprocess.run([
                str(binary), "-V", str(source), "-O", str(mode_table), "-F", "POS",
                "--show-filtered", "-L", "chr1:1-2",
                "--variant-output-filtering", mode,
                "--output-manifest", str(mode_manifest),
            ], text=True, capture_output=True, check=False)
            assert mode_result.returncode == 0, mode_result.stderr
            _, mode_rows = parse_table(mode_table)
            assert mode_rows == [["1"], ["2"]]
            mode_metadata = json.loads(mode_manifest.read_text(encoding="utf-8"))
            assert mode_metadata["compatibility"]["variant_output_filtering"] == mode

        # Non-ANYWHERE modes require an interval selector in GATK's shared
        # argument validation; keep the same fail-closed behavior natively.
        missing_selector = subprocess.run([
            str(binary), "-V", str(source), "-O", str(work / "missing-selector.tsv"),
            "-F", "POS", "--variant-output-filtering", "CONTAINED",
        ], text=True, capture_output=True, check=False)
        assert missing_selector.returncode != 0
        assert "requires at least one" in missing_selector.stderr

        # Verify the interval-merging spelling/manifest independently of the
        # output set: abutting selectors have the same union, but ALL versus
        # OVERLAPPING_ONLY is a compatibility-visible parser contract.
        merging_manifest = work / "merging.manifest.json"
        merging_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(work / "merging.tsv"),
            "-F", "POS", "--show-filtered", "-L", "chr1:1-1", "-L", "chr1:2-2",
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--output-manifest", str(merging_manifest),
        ], text=True, capture_output=True, check=False)
        assert merging_result.returncode == 0, merging_result.stderr
        _, merging_rows = parse_table(work / "merging.tsv")
        assert merging_rows == [["1"], ["2"]]
        merging_metadata = json.loads(merging_manifest.read_text(encoding="utf-8"))
        assert merging_metadata["compatibility"]["interval_merging_rule"] == "OVERLAPPING_ONLY"
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            # GATK requires a random-access feature index for interval
            # exclusion; build the sidecar only for this oracle invocation.
            index_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "IndexFeatureFile",
                "-I", str(source),
            ], text=True, capture_output=True, check=False)
            assert index_result.returncode == 0, index_result.stderr
            oracle_excluded = work / "oracle-excluded.tsv"
            oracle_excluded_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle_excluded), "-F", "POS",
                "--show-filtered", "-XL", "chr1:1-1",
            ], text=True, capture_output=True, check=False)
            assert oracle_excluded_result.returncode == 0, oracle_excluded_result.stderr
            assert oracle_excluded.read_text(encoding="utf-8") == excluded_table.read_text(encoding="utf-8")

            padded = work / "padded.tsv"
            padded_result = subprocess.run([
                str(binary), "-V", str(source), "-O", str(padded), "-F", "POS",
                "--show-filtered", "-L", "chr1:2-2", "-ip", "1",
            ], text=True, capture_output=True, check=False)
            assert padded_result.returncode == 0, padded_result.stderr
            _, padded_rows = parse_table(padded)
            assert padded_rows == [["1"], ["2"]]
            oracle_padded = work / "oracle-padded.tsv"
            oracle_padded_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle_padded), "-F", "POS",
                "--show-filtered", "-L", "chr1:2-2", "-ip", "1",
            ], text=True, capture_output=True, check=False)
            assert oracle_padded_result.returncode == 0, oracle_padded_result.stderr
            assert oracle_padded.read_text(encoding="utf-8") == padded.read_text(encoding="utf-8")

            # Java's table writer accepts the inherited mode and emits the
            # same overlap traversal.  Compare one representative mode
            # byte-for-byte against the native output above.
            oracle_mode = work / "oracle-variant-output-STARTS_IN.tsv"
            oracle_mode_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle_mode), "-F", "POS",
                "--show-filtered", "-L", "chr1:1-2",
                "--variant-output-filtering", "STARTS_IN",
            ], text=True, capture_output=True, check=False)
            assert oracle_mode_result.returncode == 0, oracle_mode_result.stderr
            assert oracle_mode.read_text(encoding="utf-8") == (work / "variant-output-STARTS_IN.tsv").read_text(encoding="utf-8")

        # Interval selection uses record-span overlap for gVCF blocks.  A
        # reference block with POS=1 and END=5 is selected by chr1:4-4.
        block_source = work / "block.vcf"
        block_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\n",
            encoding="utf-8")
        block_table = work / "block.tsv"
        block_result = subprocess.run([
            str(binary), "-V", str(block_source), "-O", str(block_table),
            "-L", "chr1:4-4", "-F", "CHROM", "-F", "POS", "-F", "END",
        ], text=True, capture_output=True, check=False)
        assert block_result.returncode == 0, block_result.stderr
        _, block_rows = parse_table(block_table)
        assert block_rows == [["chr1", "1", "5"]]

        # Keep the standard-field surface aligned with the pinned GATK
        # release.  TRANSVERSION/INDEL/N_ALLELES are not VariantsToTable
        # standard fields and must therefore render NA, while spanning
        # deletion `*` is classified by VariantContext as a one-base allele.
        standard_source = work / "standard-fields.vcf"
        standard_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\tstar\tA\tG,*\t50\tPASS\t.\n"
            "chr1\t2\tindel\tAT\tA,*\t50\tPASS\t.\n"
            "chr1\t3\tmissing\tA\tG\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        standard_table = work / "standard-fields.tsv"
        standard_args = [
            "-F", "ALT", "-F", "TYPE", "-F", "TRANSITION", "-F", "TRANSVERSION",
            "-F", "INDEL", "-F", "MULTI-ALLELIC", "-F", "N_ALLELES", "-F", "EVENTLENGTH", "-F", "QUAL",
        ]
        standard_result = subprocess.run(
            [str(binary), "-V", str(standard_source), "-O", str(standard_table), *standard_args],
            text=True, capture_output=True, check=False)
        assert standard_result.returncode == 0, standard_result.stderr
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            oracle_standard = work / "oracle-standard-fields.tsv"
            oracle_standard_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(standard_source), "-O", str(oracle_standard), *standard_args,
            ], text=True, capture_output=True, check=False)
            assert oracle_standard_result.returncode == 0, oracle_standard_result.stderr
            assert oracle_standard.read_text(encoding="utf-8") == standard_table.read_text(encoding="utf-8")

        # Namespace prefixes are intentionally not interpreted as bare VCF
        # keys; this is the Java tool's missing-field behavior.
        prefixed = work / "prefixed.tsv"
        prefixed_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(prefixed), "-F", "INFO/DP", "-GF", "FORMAT/GT",
        ], text=True, capture_output=True, check=False)
        assert prefixed_result.returncode == 0, prefixed_result.stderr
        prefixed_header, prefixed_rows = parse_table(prefixed)
        assert prefixed_header == ["INFO/DP", "ALPHA.FORMAT/GT", "ZETA.FORMAT/GT"]
        assert prefixed_rows == [["NA", "NA", "NA"]]

        wildcard_source = work / "wildcard.vcf"
        wildcard_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=ac>\n"
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=an>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\ta\tA\tG,T\t50\tPASS\tAC=2,3;AN=5\n",
            encoding="utf-8",
        )
        wildcard_table = work / "wildcard.tsv"
        wildcard_result = subprocess.run([
            str(binary), "-V", str(wildcard_source), "-O", str(wildcard_table), "-F", "A*",
        ], text=True, capture_output=True, check=False)
        assert wildcard_result.returncode == 0, wildcard_result.stderr
        assert wildcard_table.read_text(encoding="utf-8").splitlines() == ["A*", "5,[2, 3]"]
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            oracle_wildcard = work / "oracle-wildcard.tsv"
            oracle_wildcard_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(wildcard_source), "-O", str(oracle_wildcard), "-F", "A*",
            ], text=True, capture_output=True, check=False)
            assert oracle_wildcard_result.returncode == 0, oracle_wildcard_result.stderr
            assert oracle_wildcard.read_text(encoding="utf-8") == wildcard_table.read_text(encoding="utf-8")

        # The pinned GATK oracle is available in this workspace.  Compare the
        # complete table, including getter formatting and allele-valued GT,
        # rather than only checking a few cells.
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            oracle = work / "oracle.tsv"
            oracle_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle),
                "-F", "CHROM", "-F", "POS", "-F", "ALT", "-F", "MULTI-ALLELIC",
                "-F", "AC", "-F", "QUAL", "-F", "TYPE", "-F", "TRANSITION",
                "-F", "EVENTLENGTH", "-F", "HET", "-F", "HOM-REF", "-F", "HOM-VAR",
                "-F", "NO-CALL", "-F", "VAR", "-F", "NSAMPLES", "-F", "NCALLED",
                "-F", "SAMPLE_NAME", "-GF", "GT", "-GF", "AD", "-GF", "DP",
            ], text=True, capture_output=True, check=False)
            assert oracle_result.returncode == 0, oracle_result.stderr
            assert oracle.read_text(encoding="utf-8") == table.read_text(encoding="utf-8")
            oracle_raw = work / "oracle-raw.tsv"
            oracle_raw_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle_raw),
                "-F", "CHROM", "-F", "FILTER", "--show-filtered",
            ], text=True, capture_output=True, check=False)
            assert oracle_raw_result.returncode == 0, oracle_raw_result.stderr
            assert oracle_raw.read_text(encoding="utf-8") == raw_table.read_text(encoding="utf-8")

        split = work / "split.tsv"
        split_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(split), "-SMA",
            "-F", "CHROM", "-F", "ALT", "-F", "AC", "-GF", "AD",
        ], text=True, capture_output=True, check=False)
        assert split_result.returncode == 0, split_result.stderr
        _, split_rows = parse_table(split)
        assert len(split_rows) == 2
        assert [row[1] for row in split_rows] == ["G", "T"]
        assert [row[2] for row in split_rows] == ["2", "3"]

        split_as = work / "split-as.tsv"
        split_as_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(split_as), "-SMA",
            "-F", "CHROM", "-F", "ALT", "-GF", "PL", "--asGenotypeFieldsToTake", "AD",
        ], text=True, capture_output=True, check=False)
        assert split_as_result.returncode == 0, split_as_result.stderr
        split_as_header, split_as_rows = parse_table(split_as)
        assert len(split_as_rows) == 2
        # Number=G is retained on every ALT row; AD is the special R field
        # rendered as ref,selected-alt by GATK's ASGF path.
        assert split_as_header == ["CHROM", "ALT", "ALPHA.PL", "ALPHA.AD", "ZETA.PL", "ZETA.AD"]
        assert split_as_rows[0][2] == split_as_rows[1][2]
        assert split_as_rows[0][4] == split_as_rows[1][4]
        assert split_as_rows[0][3] == "8,1" and split_as_rows[1][3] == "8,5"
        assert split_as_rows[0][5] == "8,4" and split_as_rows[1][5] == "8,2"
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            oracle_split_as = work / "oracle-split-as.tsv"
            oracle_split_as_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(source), "-O", str(oracle_split_as), "-SMA",
                "-F", "CHROM", "-F", "ALT", "-GF", "PL", "-ASGF", "AD",
            ], text=True, capture_output=True, check=False)
            assert oracle_split_as_result.returncode == 0, oracle_split_as_result.stderr
            assert oracle_split_as.read_text(encoding="utf-8") == split_as.read_text(encoding="utf-8")

        type_source = work / "types.vcf"
        type_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\tsnp\tA\tG\t50\tPASS\t.\n"
            "chr1\t2\tins\tA\tAT\t50\tPASS\t.\n"
            "chr1\t3\tdel\tAT\tA\t50\tPASS\t.\n"
            "chr1\t4\tmnp\tAT\tGC\t50\tPASS\t.\n"
            "chr1\t5\tsymbolic\tA\t<NON_REF>\t50\tPASS\t.\n"
            "chr1\t6\tmixed\tA\tG,AT\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        type_table = work / "types.tsv"
        type_command = [
            str(binary), "-V", str(type_source), "-O", str(type_table),
            "-F", "CHROM", "-F", "POS", "-F", "TYPE", "-F", "EVENTLENGTH", "-F", "TRANSITION",
        ]
        type_result = subprocess.run(type_command, text=True, capture_output=True, check=False)
        assert type_result.returncode == 0, type_result.stderr
        if oracle_guard.oracle_ready('verify_variants_to_table.py', java, gatk_jar):
            oracle_types = work / "oracle-types.tsv"
            oracle_types_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "VariantsToTable",
                "-V", str(type_source), "-O", str(oracle_types),
                "-F", "CHROM", "-F", "POS", "-F", "TYPE", "-F", "EVENTLENGTH", "-F", "TRANSITION",
            ], text=True, capture_output=True, check=False)
            assert oracle_types_result.returncode == 0, oracle_types_result.stderr
            assert oracle_types.read_text(encoding="utf-8") == type_table.read_text(encoding="utf-8")

        molten = work / "molten.tsv"
        molten_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(molten), "-F", "CHROM",
            "-F", "POS", "-GF", "GT", "--moltenize",
        ], text=True, capture_output=True, check=False)
        assert molten_result.returncode == 0, molten_result.stderr
        molten_lines = molten.read_text(encoding="utf-8").splitlines()
        assert molten_lines[0] == "RecordID\tSample\tVariable\tValue"
        assert "1\tsite\tCHROM\tchr1" in molten_lines
        assert any(line.startswith("1\tALPHA\tGT\t") for line in molten_lines)

        missing = work / "missing.tsv"
        missing_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(missing), "-F", "MISSING",
            "--error-if-missing-data",
        ], text=True, capture_output=True, check=False)
        assert missing_result.returncode != 0
        assert "missing VariantsToTable" in missing_result.stderr

    print(json.dumps({"status": "pass", "records": 1, "split_rows": 2,
                      "gatk_oracle": (java.is_file() and gatk_jar.is_file())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
