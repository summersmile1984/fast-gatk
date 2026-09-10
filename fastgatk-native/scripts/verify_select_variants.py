#!/usr/bin/env python3
"""Verify native SelectVariants type/sample/filter/allele-subset contracts."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FILTER=<ID=q10,Description=Low quality>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AS_INT,Number=A,Type=Integer,Description=Per-ALT score>
##INFO=<ID=AS_FLOAT,Number=A,Type=Float,Description=Per-ALT frequency>
##INFO=<ID=ADINFO,Number=R,Type=Integer,Description=Per-allele score>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def normalized_records(path: Path) -> list[dict[str, object]]:
    """Normalize record fields so GATK header/FORMAT ordering is immaterial."""
    sample_names: list[str] = []
    records: list[dict[str, object]] = []
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, object] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, _, value = item.partition("=")
                    if not value:
                        info[key] = True
                    else:
                        try:
                            info[key] = float(value) if "." in value else int(value)
                        except ValueError:
                            info[key] = value
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            samples = {}
            for name, sample in zip(sample_names, fields[9:]):
                samples[name] = dict(zip(format_keys, sample.split(":")))
            records.append({
                "site": tuple(fields[index] for index in (0, 1, 2, 3, 4, 5, 6)),
                "info": info,
                "samples": samples,
            })
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SELECT_VARIANTS_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-select-variants")))
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-select-variants-") as temp:
        work = Path(temp)
        source = work / "input.vcf.gz"
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            # Keep S2's second alternate; S1's first alternate is removed by
            # --remove-unused-alternates after sample subsetting.
            handle.write(
                "chr1\t1\trs1\tA\tG,T\t50\tPASS\tDP=20;AS_INT=11,22;AS_FLOAT=0.1,0.2;ADINFO=100,11,22\tGT:AD:PL:GQ\t0/1:10,10,0:50,0,50,99,99,99:50\t0/2:8,0,8:80,99,80,99,99,0:1\n"
            )
            # This filtered SNP is intentionally dropped.
            handle.write(
                "chr1\t2\tfoo\tC\tT\t10\tq10\tDP=5\tGT:AD:PL\t0/1:2,3:20,0,20\t0/0:5,0:0,20,80\n"
            )
            # A SNP+INDEL record exercises MIXED classification.
            handle.write(
                "chr1\t3\tbar\tA\tG,AT\t60\tPASS\tDP=30\tGT:AD:PL\t0/1:15,15,0:50,0,80,99,99,99\t0/0:30,0,0:0,50,70,99,99,99\n"
            )

        result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(output),
             "--sample-name", "S2", "--select-type-to-include", "SNP",
             "--select-expression", 'vc.getAttribute("DP") >= 10 && QUAL > 40',
             "--exclude-filtered", "--remove-unused-alternates",
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["output_records"] == 1
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert len(records) == 1
        assert "\tFORMAT\tS2" in text
        assert "\tFORMAT\tS1" not in text
        assert records[0][4] == "T"
        assert records[0][9].split(":")[0] == "0/1"
        assert len(records[0][9].split(":")[1].split(",")) == 2
        assert len(records[0][9].split(":")[2].split(",")) == 3
        assert records[0][9].split(":")[-1] == "80"
        assert "AC=1" in records[0][7] and "AN=2" in records[0][7]
        assert "AF=0.5" in records[0][7]
        # GATK intentionally preserves custom Number=A/R INFO payloads when
        # removing unused alternates; only the genotype-dependent annotations
        # (AC/AN/AF) and FORMAT fields are rewritten.
        assert "AS_INT=11,22" in records[0][7]
        assert "AS_FLOAT=0.1,0.2" in records[0][7]
        assert "ADINFO=100,11,22" in records[0][7]

        # The same supported boundary is compared with GATK's real
        # SelectVariants implementation.  Record comparison is semantic:
        # htsjdk may reorder INFO/FORMAT fields and format AF with three
        # decimals, while the allele/PL/GT/AD/GQ values must be identical.
        gatk_java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        source_plain = work / "input.vcf"
        source_plain.write_text(gzip.open(source, "rt", encoding="utf-8").read(), encoding="utf-8")
        gatk_output = work / "gatk-output.vcf"
        assert gatk_java.exists() and gatk_jar.exists()
        gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(source_plain), "-O", str(gatk_output),
             "--sample-name", "S2", "--select-type-to-include", "SNP",
             "--select", 'vc.getAttribute("DP") >= 10 && QUAL > 40',
             "--exclude-filtered", "--remove-unused-alternates"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_result.returncode == 0, gatk_result.stderr
        assert normalized_records(output) == normalized_records(gatk_output), (
            normalized_records(output), normalized_records(gatk_output)
        )

        # Plain-text VCF output uses Tribble LinearIndex v3, just like the
        # Java VariantContextWriter.  Prove that the bundled GATK reader can
        # perform an interval query through the native `.idx` sidecar.
        plain_output = work / "output-plain.vcf"
        plain_manifest = work / "output-plain.manifest.json"
        plain_result = subprocess.run(
            [str(binary), "-V", str(source_plain), "-O", str(plain_output),
             "--sample-name", "S2", "--select-type-to-include", "SNP",
             "--select-expression", 'vc.getAttribute("DP") >= 10 && QUAL > 40',
             "--exclude-filtered", "--remove-unused-alternates",
             "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.is_file() and Path(f"{plain_output}.idx").is_file()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "output-plain-query.vcf"
        plain_query_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(plain_output), "-L", "chr1:1-1", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert len(normalized_records(plain_query)) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True
        assert all(item["complete"] for item in plain_metadata["outputs"])
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["type_filter"] is True
        assert metadata["compatibility"]["sample_subset"] is True
        assert metadata["compatibility"]["remove_unused_alternates"] is True
        assert metadata["compatibility"]["site_annotation_recompute"] is True
        assert metadata["compatibility"]["jexl_expression_subset"] is True
        assert metadata["compatibility"]["jexl_info_vector_indexing"] is True
        assert metadata["compatibility"]["jexl_boolean_comparisons"] is True
        assert metadata["compatibility"]["jexl_regex_operators"] is True
        assert metadata["compatibility"]["jexl_arithmetic"] is True
        assert metadata["telemetry"]["pl_remap_kernel_calls"] == 1
        assert metadata["telemetry"]["allele_field_remap_kernel_calls"] == 1
        expected_execution_space = os.environ.get("FASTGATK_EXPECTED_EXECUTION_SPACE", "OpenMP")
        assert metadata["telemetry"]["allele_field_remap_kernel_execution_space"] == expected_execution_space
        assert metadata["telemetry"]["gt_gq_kernel_calls"] == 1
        assert metadata["telemetry"]["genotype_kernel_execution_space"] == expected_execution_space
        assert all(item["complete"] for item in metadata["outputs"])

        # --exclude-non-variants follows VariantContext.isVariant(): a
        # sample-bearing site with only hom-ref/no-call genotypes is dropped
        # (including reference-confidence <NON_REF> blocks), while a called
        # ALT remains selectable.
        nonvariant_source = work / "nonvariants.vcf"
        nonvariant_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate allele>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t1\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\n"
            "chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\n"
            "chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t./.\n"
            "chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            "chr1\t5\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/0\n"
            "chr1\t6\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        nonvariant_native = work / "nonvariants-native.vcf"
        nonvariant_result = subprocess.run(
            [str(binary), "-V", str(nonvariant_source), "-O", str(nonvariant_native),
             "--exclude-non-variants", "--output-manifest", str(work / "nonvariants.manifest.json")],
            text=True, capture_output=True, check=False,
        )
        assert nonvariant_result.returncode == 0, nonvariant_result.stderr
        nonvariant_gatk = work / "nonvariants-gatk.vcf"
        nonvariant_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(nonvariant_source), "-O", str(nonvariant_gatk),
             "--exclude-non-variants", "true"],
            text=True, capture_output=True, check=False,
        )
        assert nonvariant_gatk_result.returncode == 0, nonvariant_gatk_result.stderr
        assert normalized_records(nonvariant_native) == normalized_records(nonvariant_gatk)
        nonvariant_metadata = json.loads(
            (work / "nonvariants.manifest.json").read_text(encoding="utf-8"))
        assert nonvariant_metadata["compatibility"]["exclude_non_variants"] is True
        assert nonvariant_metadata["telemetry"]["non_variant_skipped"] == 4

        # Absolute REF/ALT length bounds are applied only to indel records;
        # compare the minimum and maximum paths independently with GATK.
        indel_source = work / "indel-sizes.vcf"
        indel_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t1\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1\n"
            "chr1\t2\t.\tATG\tA\t50\tPASS\t.\tGT\t0/1\n"
            "chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        for bound_args, stem, expected_positions in (
            (["--min-indel-size", "2"], "min-indel", ["2"]),
            (["--max-indel-size", "1"], "max-indel", ["1"]),
        ):
            native_bound = work / f"{stem}-native.vcf"
            native_bound_result = subprocess.run(
                [str(binary), "-V", str(indel_source), "-O", str(native_bound),
                 "--select-type-to-include", "INDEL", *bound_args],
                text=True, capture_output=True, check=False,
            )
            assert native_bound_result.returncode == 0, native_bound_result.stderr
            gatk_bound = work / f"{stem}-gatk.vcf"
            gatk_bound_result = subprocess.run(
                [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
                 "-V", str(indel_source), "-O", str(gatk_bound),
                 "--select-type-to-include", "INDEL", *bound_args],
                text=True, capture_output=True, check=False,
            )
            assert gatk_bound_result.returncode == 0, gatk_bound_result.stderr
            assert normalized_records(native_bound) == normalized_records(gatk_bound)
            assert [record["site"][1] for record in normalized_records(native_bound)] == expected_positions

        # Keep-original annotations snapshot the pre-subset values while the
        # regular AC/AN/AF fields reflect the selected sample/ALT space.
        original_source = work / "original-annotations.vcf"
        original_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr1\t1\t.\tA\tG,T\t50\tPASS\tAC=2,1;AN=4;AF=0.5,0.25;DP=20\tGT\t0/1\t0/2\n",
            encoding="utf-8",
        )
        original_native = work / "original-native.vcf"
        original_result = subprocess.run(
            [str(binary), "-V", str(original_source), "-O", str(original_native),
             "--sample-name", "S1", "--keep-original-ac", "true",
             "--keep-original-dp", "true", "--remove-unused-alternates"],
            text=True, capture_output=True, check=False,
        )
        assert original_result.returncode == 0, original_result.stderr
        original_gatk = work / "original-gatk.vcf"
        original_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(original_source), "-O", str(original_gatk),
             "--sample-name", "S1", "--keep-original-ac", "true",
             "--keep-original-dp", "true", "--remove-unused-alternates"],
            text=True, capture_output=True, check=False,
        )
        assert original_gatk_result.returncode == 0, original_gatk_result.stderr
        assert normalized_records(original_native) == normalized_records(original_gatk)
        original_record = normalized_records(original_native)[0]
        assert original_record["info"]["AC_Orig"] == 2
        assert original_record["info"]["AF_Orig"] == 0.5
        assert original_record["info"]["AN_Orig"] == 4
        assert original_record["info"]["DP_Orig"] == 20
        original_metadata = json.loads(
            Path(str(original_native) + ".manifest.json").read_text(encoding="utf-8"))
        assert original_metadata["compatibility"]["keep_original_ac"] is True
        assert original_metadata["compatibility"]["keep_original_dp"] is True

        # Region selection is span-overlap based for gVCF records.  A
        # reference block at POS=1 with END=5 is retained by chr1:4-4.
        block_source = work / "block.vcf.gz"
        with gzip.open(block_source, "wt", encoding="utf-8") as handle:
            handle.write(
                "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\n")
        block_output = work / "block-output.vcf.gz"
        block_result = subprocess.run([
            str(binary), "-V", str(block_source), "-O", str(block_output),
            "-L", "chr1:4-4"], text=True, capture_output=True, check=False)
        assert block_result.returncode == 0, block_result.stderr
        block_records = [line for line in gzip.open(
            block_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(block_records) == 1 and block_records[0].split("\t")[1] == "1"

        # Repeated genotype JEXL is evaluated over all samples and selects the
        # site when any sample satisfies it, matching SelectVariants' default
        # genotype-expression semantics.
        genotype_output = work / "genotype-selected.vcf.gz"
        genotype_manifest = work / "genotype-selected.manifest.json"
        genotype_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(genotype_output),
             "--select-genotype", "GQ >= 50", "--output-manifest", str(genotype_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert genotype_result.returncode == 0, genotype_result.stderr
        genotype_records = [line.split("\t") for line in gzip.open(
            genotype_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(genotype_records) == 1 and genotype_records[0][1] == "1"
        genotype_metadata = json.loads(genotype_manifest.read_text(encoding="utf-8"))
        assert genotype_metadata["compatibility"]["genotype_expression_subset"] is True

        het_output = work / "het-selected.vcf.gz"
        het_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(het_output),
             "--select-genotype", "isHet()"],
            text=True, capture_output=True, check=False,
        )
        assert het_result.returncode == 0, het_result.stderr
        het_records = [line.split("\t") for line in gzip.open(
            het_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in het_records] == ["1", "2", "3"]

        # Multiple site-level -select expressions are OR-combined by GATK;
        # conjunction remains available inside one expression.
        or_output = work / "or-selected.vcf.gz"
        or_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(or_output),
             "--select-expression", "QUAL > 55", "--select-expression", "QUAL < 15"],
            text=True, capture_output=True, check=False,
        )
        assert or_result.returncode == 0, or_result.stderr
        or_records = [line.split("\t") for line in gzip.open(
            or_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in or_records] == ["2", "3"]

        mixed_output = work / "mixed.vcf.gz"
        mixed = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(mixed_output),
             "--select-type-to-include", "MIXED"],
            text=True, capture_output=True, check=False,
        )
        assert mixed.returncode == 0, mixed.stderr
        mixed_text = gzip.open(mixed_output, "rt", encoding="utf-8").read()
        assert sum(1 for line in mixed_text.splitlines() if line and not line.startswith("#")) == 1
        assert "\t3\t" in mixed_text

        # Common GATK VariantContext methods are evaluated directly against
        # the native record, including boolean predicates, INFO presence,
        # coordinates, and string-valued type comparisons.
        def method_records(expression: str, stem: str, input_path: Path = source,
                           genotype_expression: bool = False) -> list[str]:
            method_output = work / f"{stem}.vcf.gz"
            expression_option = "--select-genotype" if genotype_expression else "--select-expression"
            method_result = subprocess.run(
                [str(binary), "-V", str(input_path), "-O", str(method_output),
                 expression_option, expression],
                text=True, capture_output=True, check=False,
            )
            assert method_result.returncode == 0, method_result.stderr
            return [line.split("\t")[1] for line in gzip.open(
                method_output, "rt", encoding="utf-8").read().splitlines()
                if line and not line.startswith("#")]

        # GATK also exposes genotype predicates as numeric convenience
        # fields.  These must be evaluated per sample and then use the
        # SelectVariants any-sample selection rule.
        assert method_records("isHet == 1", "numeric-het", genotype_expression=True) == ["1", "2", "3"]
        assert method_records("isHomRef == 1", "numeric-homref", genotype_expression=True) == ["2", "3"]
        assert method_records("GQ + 1 >= 51", "arithmetic-genotype", genotype_expression=True) == ["1"]
        assert method_records('vc.getGenotype("S2").isHomRef() == true',
                              "explicit-genotype-bool") == ["2", "3"]

        assert method_records("vc.isSNP()", "is-snp") == ["1", "2"]
        assert method_records("vc.isSNP() == true", "is-snp-bool") == ["1", "2"]
        assert method_records("vc.isMixed()", "is-mixed") == ["3"]
        assert method_records("vc.isMultiallelic()", "is-multiallelic") == ["1", "3"]
        assert method_records("vc.isTransition()", "is-transition") == ["2"]
        assert method_records("vc.isFiltered()", "is-filtered") == ["2"]
        assert method_records('vc.getFilters().contains("q10")', "filter-contains") == ["2"]
        assert method_records("vc.isPass()", "is-pass") == ["1", "3"]
        assert method_records('vc.hasAttribute("DP")', "has-dp") == ["1", "2", "3"]
        assert method_records("vc.getStart() >= 2", "start") == ["2", "3"]
        assert method_records('vc.getType() == "SNP"', "type-snp") == ["1", "2"]
        assert method_records("vc.getType() == VariantContext.Type.SNP", "type-enum") == ["1", "2"]
        assert method_records('vc.getContig() == "chr1"', "contig") == ["1", "2", "3"]
        assert method_records("vc.getNAlleles() >= 3", "n-alleles") == ["1", "3"]
        assert method_records("vc.hasAlternateAllele(1)", "alternate-index") == ["1", "3"]
        assert method_records('vc.getID() == "rs1"', "id") == ["1"]
        assert method_records('vc.getID() !~ "rs.*"', "id-regex-not") == ["2", "3"]
        assert method_records('vc.getAttribute("DP").contains("2")', "numeric-contains") == ["1"]
        # Compound numeric JEXL operands use Java precedence and can read
        # scalar fields (INFO vector arithmetic is intentionally left to the
        # explicit indexed predicate path below, matching SelectVariants' JEXL
        # evaluator boundary).
        assert method_records("QUAL / 2 >= 25", "arithmetic-qual") == ["1", "3"]
        assert method_records("(QUAL - 10) / 2 >= 20", "arithmetic-paren") == ["1", "3"]
        # INFO vectors are addressable with the same zero-based JEXL index
        # used by HTSJDK.  Exercise both compact TAG[i] and explicit
        # vc.getAttribute("TAG")[i] spellings, including Number=A and
        # Number=R payloads.
        assert method_records("AS_INT[0] >= 11", "info-array-int") == ["1"]
        assert method_records('vc.getAttribute("AS_INT")[1] > 20',
                              "info-array-int-explicit") == ["1"]
        assert method_records('vc.getAttribute("AS_INT").1 > 20',
                              "info-array-int-dot") == ["1"]
        assert method_records('vc.getAttribute("AS_INT").get(1) > 20',
                              "info-array-int-get") == ["1"]
        assert method_records("AS_FLOAT[1] >= 0.2", "info-array-float") == ["1"]
        assert method_records("ADINFO[1] == 11", "info-array-number-r") == ["1"]
        assert method_records('vc.getType() =~ "SNP"', "type-regex") == ["1", "2"]
        assert method_records('vc.getType() !~ "SNP"', "type-regex-not") == ["3"]
        assert method_records('vc.getAttribute("MISSING") == null', "missing-null") == ["1", "2", "3"]
        assert method_records("vc.isNotFiltered()", "not-filtered") == ["1", "3"]
        assert method_records("vc.getAlleles().size() >= 3", "alleles-size") == ["1", "3"]
        assert method_records("vc.getGenotypes().size() == 2", "genotypes-size") == ["1", "2", "3"]
        assert method_records("vc.getFilters().isEmpty()", "filters-empty") == ["1", "3"]
        assert method_records("vc.getReference().isReference()", "reference-allele") == ["1", "2", "3"]
        assert method_records("vc.getAlternateAllele(0).isNonReference()", "non-reference-allele") == ["1", "2", "3"]
        assert method_records("vc.getAlternateAllele(0).isNonRefAllele()", "non-ref-allele", block_source) == ["1"]
        assert method_records("vc.getAlternateAllele(0).isCalled()", "called-allele", block_source) == ["1"]
        assert method_records("vc.getAlternateAllele(0).isSymbolic() == true",
                              "symbolic-allele-bool", block_source) == ["1"]

        # Explicit boolean comparisons and !~ are checked against the real
        # GATK SelectVariants evaluator as well as the compact assertions
        # above.  This closes the gap where a parser could accept the syntax
        # but silently route it through a numeric/missing field path.
        for expression, stem in (
            ("vc.isSNP() == true", "bool-site-oracle"),
            ('vc.getID() !~ "rs.*"', "regex-not-oracle"),
            ('vc.getType() !~ "SNP"', "type-regex-not-oracle"),
            ("vc.getAlternateAllele(0).isSymbolic() == true", "bool-allele-oracle"),
            ('vc.getGenotype("S2").isHomRef() == true', "explicit-genotype-oracle"),
            ("QUAL / 2 >= 25", "arithmetic-qual-oracle"),
            ("(QUAL - 10) / 2 >= 20", "arithmetic-paren-oracle"),
        ):
            native_path = work / f"{stem}-native.vcf.gz"
            native_result = subprocess.run(
                [str(binary), "-V", str(source), "-O", str(native_path),
                 "--select-expression", expression],
                text=True, capture_output=True, check=False)
            assert native_result.returncode == 0, native_result.stderr
            gatk_path = work / f"{stem}-gatk.vcf"
            gatk_result = subprocess.run(
                [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
                 "-V", str(source_plain), "-O", str(gatk_path),
                 "--select", expression],
                text=True, capture_output=True, check=False)
            assert gatk_result.returncode == 0, gatk_result.stderr
            # Without --remove-unused-alternates GATK leaves pre-existing
            # genotype-derived AC/AN/AF annotations untouched, while the
            # native path deterministically recomputes them.  The oracle for
            # a pure selection expression is therefore the selected site and
            # allele/sample payload, not incidental annotation refresh.
            native_sites = [(row["site"], row["samples"]) for row in normalized_records(native_path)]
            gatk_sites = [(row["site"], row["samples"]) for row in normalized_records(gatk_path)]
            assert native_sites == gatk_sites

        # HTSJDK Allele methods are also valid SelectVariants JEXL operands.
        # Symbolic alleles have length zero in HTSJDK, while primitive REF/ALT
        # lengths retain their sequence length.
        allele_method_source = work / "allele-methods.vcf.gz"
        allele_method_header = HEADER.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "##ALT=<ID=DEL,Description=Deletion>\n")
        with gzip.open(allele_method_source, "wt", encoding="utf-8") as handle:
            handle.write(allele_method_header)
            handle.write("chr1\t1\t.\tA\tG\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
            handle.write("chr1\t2\t.\tA\t<DEL>\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
            handle.write("chr1\t3\t.\tAT\tA\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
            handle.write("chr1\t4\t.\tAT\tGC\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
        allele_method_output = work / "allele-methods-native.vcf.gz"
        allele_method_result = subprocess.run(
            [str(binary), "-V", str(allele_method_source), "-O", str(allele_method_output),
             "--select-expression",
             "vc.getAlternateAllele(0).isSymbolic() || "
             "vc.getAlternateAllele(0).length() > 1 || "
             "vc.getReference().length() > 1"],
            text=True, capture_output=True, check=False,
        )
        assert allele_method_result.returncode == 0, allele_method_result.stderr
        allele_method_records = [line.split("\t") for line in gzip.open(
            allele_method_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in allele_method_records] == ["2", "3", "4"]
        allele_method_plain = work / "allele-methods.vcf"
        allele_method_plain.write_text(
            gzip.open(allele_method_source, "rt", encoding="utf-8").read(), encoding="utf-8")
        allele_method_gatk = work / "allele-methods-gatk.vcf"
        allele_method_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(allele_method_plain), "-O", str(allele_method_gatk),
             "--select",
             "vc.getAlternateAllele(0).isSymbolic() || "
             "vc.getAlternateAllele(0).length() > 1 || "
             "vc.getReference().length() > 1"],
            text=True, capture_output=True, check=False,
        )
        assert allele_method_gatk_result.returncode == 0, allele_method_gatk_result.stderr
        assert normalized_records(allele_method_output) == normalized_records(allele_method_gatk)

        breakend_source = work / "breakend-methods.vcf.gz"
        breakend_header = HEADER.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n")
        with gzip.open(breakend_source, "wt", encoding="utf-8") as handle:
            handle.write(breakend_header)
            handle.write("chr1\t5\t.\tA\tN]chr1:10]\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
            handle.write("chr1\t6\t.\tA\t.N\t50\tPASS\tDP=10;AC=2;AN=4;AF=0.5\tGT\t0/1\t0/1\n")
        breakend_output = work / "breakend-methods-native.vcf.gz"
        breakend_result = subprocess.run(
            [str(binary), "-V", str(breakend_source), "-O", str(breakend_output),
             "--select-expression",
             "vc.getAlternateAllele(0).isBreakpoint() || vc.getAlternateAllele(0).isSingleBreakend()"],
            text=True, capture_output=True, check=False,
        )
        assert breakend_result.returncode == 0, breakend_result.stderr
        breakend_records = [line.split("\t") for line in gzip.open(
            breakend_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in breakend_records] == ["5", "6"]
        breakend_plain = work / "breakend-methods.vcf"
        breakend_plain.write_text(gzip.open(breakend_source, "rt", encoding="utf-8").read(), encoding="utf-8")
        breakend_gatk = work / "breakend-methods-gatk.vcf"
        breakend_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(breakend_plain), "-O", str(breakend_gatk), "--select",
             "vc.getAlternateAllele(0).isBreakpoint() || vc.getAlternateAllele(0).isSingleBreakend()"],
            text=True, capture_output=True, check=False,
        )
        assert breakend_gatk_result.returncode == 0, breakend_gatk_result.stderr
        assert normalized_records(breakend_output) == normalized_records(breakend_gatk)

        genotype_method_source = work / "genotype-methods.vcf.gz"
        genotype_method_header = HEADER
        with gzip.open(genotype_method_source, "wt", encoding="utf-8") as handle:
            handle.write(genotype_method_header)
            handle.write("chr1\t1\t.\tA\tG\t50\tPASS\tDP=20\tGT:AD:PL:GQ\t0/1:5,7:50,0,50:50\t0/0:10,0:0,30,60:20\n")
            handle.write("chr1\t2\t.\tC\tT\t10\tq10\tDP=5\tGT:AD:PL:GQ\t0/0:10,0:0,30,60:10\t0/0:10,0:0,30,60:20\n")
            handle.write("chr1\t3\t.\tG\tA\t60\tPASS\tDP=30\tGT:AD:PL:GQ\t1/1:0,20:80,40,0:60\t0/0:10,0:0,30,60:20\n")
        assert method_records('vc.getGenotype("S1").getAD()[1] >= 7',
                              "genotype-ad", genotype_method_source, True) == ["1", "3"]
        assert method_records('vc.getGenotype("S1").getAD().1 >= 7',
                              "genotype-ad-dot", genotype_method_source, True) == ["1", "3"]
        assert method_records('vc.getGenotype("S1").getPL()[2] < 10',
                              "genotype-pl", genotype_method_source, True) == ["3"]
        assert method_records('vc.getGenotype("S1").getPL().2 < 10',
                              "genotype-pl-dot", genotype_method_source, True) == ["3"]
        assert method_records('vc.getGenotype("S1").isHom()',
                              "genotype-hom", genotype_method_source, True) == ["2", "3"]
        assert method_records('vc.getGenotype("S1").isAvailable()',
                              "genotype-available", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").isCalled()',
                              "genotype-called", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").hasAD()',
                              "genotype-has-ad", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").hasPL()',
                              "genotype-has-pl", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").hasGQ()',
                              "genotype-has-gq", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").getPloidy() == 2',
                              "genotype-ploidy", genotype_method_source, True) == ["1", "2", "3"]
        assert method_records('vc.getGenotype("S1").getPhredScaledQual() >= 50',
                              "genotype-phred", genotype_method_source, True) == ["1", "3"]
        min_dp_source = work / "min-dp-genotype.vcf.gz"
        min_dp_header = genotype_method_header.replace(
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n",
            "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum depth>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n")
        with gzip.open(min_dp_source, "wt", encoding="utf-8") as handle:
            handle.write(min_dp_header)
            handle.write("chr1\t1\t.\tA\tG\t50\tPASS\tDP=5\tGT:MIN_DP:AD:PL:GQ\t0/1:5:5,7:50,0,50:50\t0/0:20:10,0:0,30,60:20\n")
            handle.write("chr1\t2\t.\tC\tT\t50\tPASS\tDP=20\tGT:MIN_DP:AD:PL:GQ\t0/0:20:10,0:0,30,60:20\t0/0:20:10,0:0,30,60:20\n")
        assert method_records("MIN_DP < 10", "genotype-min-dp", min_dp_source, True) == ["1"]
        assert method_records("vc.hasGenotypes()", "has-genotypes", genotype_method_source) == ["1", "2", "3"]
        assert method_records("vc.getCalledChrCount() == 4", "called-chr", genotype_method_source) == ["1", "2", "3"]
        assert method_records("vc.getNoCallCount() == 0", "no-call-count", genotype_method_source) == ["1", "2", "3"]
        assert method_records("vc.getHomRefCount() == 2", "hom-ref-count", genotype_method_source) == ["2"]
        assert method_records("vc.getHetCount() == 1", "het-count", genotype_method_source) == ["1"]
        assert method_records("vc.getHomVarCount() == 1", "hom-var-count", genotype_method_source) == ["3"]
        assert method_records("vc.isPolymorphicInSamples()", "polymorphic", genotype_method_source) == ["1", "3"]
        assert method_records("vc.isMonomorphicInSamples()", "monomorphic", genotype_method_source) == ["2"]

        interval_output = work / "interval-selected.vcf.gz"
        interval_manifest = work / "interval-selected.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(interval_output),
             "-L", "chr1:1-1", "--intervals", "chr1:3-3",
             "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert interval_result.returncode == 0, interval_result.stderr
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["1", "3"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["interval_skipped"] == 1

        intersection_output = work / "interval-intersection.vcf.gz"
        intersection_manifest = work / "interval-intersection.manifest.json"
        intersection_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(intersection_output),
             "-L", "chr1:1-3", "-L", "chr1:3-3",
             "--interval-set-rule=INTERSECTION",
             "--output-manifest", str(intersection_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_summary = json.loads(intersection_result.stdout.splitlines()[-1])
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        intersection_records = [line.split("\t") for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in intersection_records] == ["3"]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"

        interval_list = work / "select.interval_list"
        interval_list.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t1\t+\tfirst\n"
            "chr1\t3\t3\t+\tthird\n", encoding="utf-8")
        interval_file_output = work / "interval-file-selected.vcf.gz"
        interval_file_manifest = work / "interval-file-selected.manifest.json"
        interval_file_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(interval_file_output),
             "-L", str(interval_list), "--output-manifest", str(interval_file_manifest)],
            text=True, capture_output=True, check=False)
        assert interval_file_result.returncode == 0, interval_file_result.stderr
        interval_file_records = [line.split("\t") for line in gzip.open(
            interval_file_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_file_records] == ["1", "3"]
        interval_file_metadata = json.loads(interval_file_manifest.read_text(encoding="utf-8"))
        assert interval_file_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_file_metadata["telemetry"]["interval_list_records"] == 2

        bed = work / "select.bed"
        bed.write_text("chr1\t0\t1\nchr1\t2\t3\n", encoding="utf-8")
        bed_output = work / "bed-selected.vcf.gz"
        bed_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(bed_output), "-L", str(bed)],
            text=True, capture_output=True, check=False)
        assert bed_result.returncode == 0, bed_result.stderr
        bed_records = [line.split("\t") for line in gzip.open(
            bed_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in bed_records] == ["1", "3"]

        # Concordance/discordance use a deterministic contig/POS/REF/ALT key
        # and are evaluated before the ordinary type/JEXL filters.
        comparison = work / "comparison.vcf.gz"
        with gzip.open(comparison, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t1\t.\tA\tG,T\t50\tPASS\tDP=20\tGT:AD:PL:GQ\t0/1:10,10,0:50,0,50,99,99,99:50\t0/2:8,0,8:80,99,80,99,99,0:1\n")
            handle.write("chr1\t3\t.\tA\tG,AT\t60\tPASS\tDP=30\tGT:AD:PL\t0/1:15,15,0:50,0,80,99,99,99\t0/0:30,0,0:0,50,70,99,99,99\n")
        concordant = work / "concordant.vcf.gz"
        concordant_manifest = work / "concordant.manifest.json"
        concordant_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(concordant),
             "--concordance", str(comparison), "--output-manifest", str(concordant_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert concordant_result.returncode == 0, concordant_result.stderr
        concordant_records = [line.split("\t") for line in gzip.open(
            concordant, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in concordant_records] == ["1", "3"]
        concordant_metadata = json.loads(concordant_manifest.read_text(encoding="utf-8"))
        assert concordant_metadata["compatibility"]["concordance_discordance_subset"] is True
        assert concordant_metadata["telemetry"]["comparison_skipped"] == 1

        discordant = work / "discordant.vcf.gz"
        discordant_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(discordant),
             "--discordance", str(comparison)],
            text=True, capture_output=True, check=False,
        )
        assert discordant_result.returncode == 0, discordant_result.stderr
        discordant_records = [line.split("\t") for line in gzip.open(
            discordant, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in discordant_records] == ["2"]

        # Explicit sample-level mode compares unphased GTs for shared sample
        # names, while retaining ordinary site discordance for absent sites.
        genotype_comparison = work / "genotype-comparison.vcf.gz"
        with gzip.open(genotype_comparison, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t1\t.\tA\tG,T\t50\tPASS\tDP=20\tGT:AD:PL:GQ\t0/1:10,10,0:50,0,50,99,99,99:50\t0/2:8,0,8:80,99,80,99,99,0:1\n")
            # Same site key as source chr1:3, but both sample GTs disagree.
            handle.write("chr1\t3\t.\tA\tG,AT\t60\tPASS\tDP=30\tGT:AD:PL\t1/1:15,15,0:50,0,80,99,99,99\t0/1:30,0,0:0,50,70,99,99,99\n")
        genotype_concordant = work / "genotype-concordant.vcf.gz"
        genotype_manifest = work / "genotype-concordant.manifest.json"
        genotype_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(genotype_concordant),
            "--concordance", str(genotype_comparison), "--concordance-genotypes",
            "--output-manifest", str(genotype_manifest),
        ], text=True, capture_output=True, check=False)
        assert genotype_result.returncode == 0, genotype_result.stderr
        genotype_records = [line.split("\t") for line in gzip.open(
            genotype_concordant, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in genotype_records] == ["1"]
        genotype_metadata = json.loads(genotype_manifest.read_text(encoding="utf-8"))
        assert genotype_metadata["compatibility"]["sample_level_concordance"] is True

        genotype_discordant = work / "genotype-discordant.vcf.gz"
        genotype_discordant_result = subprocess.run([
            str(binary), "-V", str(source), "-O", str(genotype_discordant),
            "--discordance", str(genotype_comparison), "--discordance-genotypes",
        ], text=True, capture_output=True, check=False)
        assert genotype_discordant_result.returncode == 0, genotype_discordant_result.stderr
        genotype_discordant_records = [line.split("\t") for line in gzip.open(
            genotype_discordant, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in genotype_discordant_records] == ["2", "3"]

        # PL is Number=G and therefore depends on both allele count and
        # ploidy.  Removing an unused ALT from a triploid record must compact
        # the 10-entry old PL vector to the 4-entry diploid-allele vector in
        # VCF genotype ordering (000, 001, 011, 111), not assume diploidy.
        triploid_source = work / "triploid.vcf.gz"
        triploid_output = work / "triploid-output.vcf.gz"
        with gzip.open(triploid_source, "wt", encoding="utf-8") as handle:
            handle.write("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t10\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:PL:GQ\t0/1/1:30,8,0:0,10,20,30,40,50,60,70,80,90:90
""")
        triploid_result = subprocess.run(
            [str(binary), "-V", str(triploid_source), "-O", str(triploid_output),
             "--remove-unused-alternates"],
            text=True, capture_output=True, check=False,
        )
        assert triploid_result.returncode == 0, triploid_result.stderr
        triploid_records = [line.split("\t") for line in gzip.open(
            triploid_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(triploid_records) == 1
        triploid_sample = triploid_records[0][9].split(":")
        assert triploid_records[0][4] == "C"
        assert triploid_sample[0] == "0/1/1"
        assert triploid_sample[1] == "30,8"
        assert triploid_sample[2] == "0,10,20,30"
        assert triploid_sample[3] == "10"
        assert "AC=2" in triploid_records[0][7] and "AN=3" in triploid_records[0][7]

        # If every sample is hom-ref, removing unused alternates produces a
        # ref-only record, which GATK drops instead of attempting an invalid
        # one-allele Number=G/Number=R remap.
        ref_only_source = work / "ref-only-unused.vcf.gz"
        ref_only_output = work / "ref-only-unused-output.vcf.gz"
        with gzip.open(ref_only_source, "wt", encoding="utf-8") as handle:
            handle.write("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:PL:GQ\t0/0:30,0,0:0,30,60,40,70,80:30
""")
        ref_only_result = subprocess.run(
            [str(binary), "-V", str(ref_only_source), "-O", str(ref_only_output),
             "--remove-unused-alternates"], text=True, capture_output=True, check=False)
        assert ref_only_result.returncode == 0, ref_only_result.stderr
        ref_only_records = [line for line in gzip.open(
            ref_only_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert ref_only_records == []

    print(json.dumps({"status": "pass", "output_records": summary["output_records"],
                      "gatk_oracle_exact": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
