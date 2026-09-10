#!/usr/bin/env python3
"""Verify deterministic simple-expression VariantFiltration semantics."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def normalized_filter_records(path: Path) -> list[tuple[tuple[str, ...], tuple[str, ...], str]]:
    records = []
    with open_vcf(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            filters = tuple(sorted(filter_name for filter_name in fields[6].split(";")
                                   if filter_name and filter_name != "PASS"))
            records.append((tuple(fields[:6]), filters,
                            fields[8] + "\t" + "\t".join(fields[9:])))
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_FILTRATION_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-variant-filtration")))
    assert binary.is_file() and os.access(binary, os.X_OK)
    gatk_java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    assert gatk_java.exists() and gatk_jar.exists()
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-filtration-") as temp:
        work = Path(temp)
        source = work / "input.vcf.gz"
        output = work / "filtered.vcf.gz"
        manifest = work / "filtered.manifest.json"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t1\trs1\tA\tG\t20\tPASS\tDP=5;QD=1.0\tGT\t0/1\n")
            handle.write("chr1\t2\tfoo\tC\tT\t50\tPASS\tDP=20;QD=3.0\tGT\t0/0\n")
            handle.write("chr1\t3\t.\tG\tA\t50\tPASS\tDP=20\tGT\t0/1\n")
        result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(output),
             "--filter-expression", "QUAL < 30", "--filter-name", "LowQual",
             "--filter-expression", "DP < 10", "--filter-name", "LowDP",
             "--filter-expression", 'vc.getAttribute("QD") < 2 || DP < 10', "--filter-name", "LowQD",
             "--filter-expression", 'vc.getGenotype("S1").isHet()', "--filter-name", "Het",
             "--missing-values-evaluate-as-failing", "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["input_records"] == 3 and summary["filtered_records"] == 2
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        assert "##FILTER=<ID=LowQual" in text
        records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert records[0][6] == "Het;LowDP;LowQD;LowQual"
        assert records[0][6] != "PASS"
        assert records[1][6] == "PASS"
        assert "LowQD" in records[2][6] and "Het" in records[2][6]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["jexl_boolean_subset"] is True
        assert metadata["compatibility"]["genotype_predicates"] is True
        assert metadata["compatibility"]["jexl_info_vector_indexing"] is True
        assert metadata["compatibility"]["jexl_boolean_comparisons"] is True
        assert metadata["compatibility"]["jexl_regex_operators"] is True
        assert metadata["compatibility"]["jexl_arithmetic"] is True
        assert metadata["compatibility"]["missing_values_failing"] is True
        assert all(item["complete"] for item in metadata["outputs"])

        # Java VariantFiltration preserves existing genotype FT labels and
        # canonicalizes the combined set lexically before adding a new label.
        # This is important when several filtration stages are chained.
        ft_source = work / "existing-ft.vcf"
        ft_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##FORMAT=<ID=FT,Number=1,Type=String,Description=Filter>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT:DP:FT\t0/1:5:Prior\t0/0:20:Existing2\n",
            encoding="utf-8",
        )
        ft_native = work / "existing-ft-native.vcf.gz"
        ft_result = subprocess.run(
            [str(binary), "-V", str(ft_source), "-O", str(ft_native),
             "--genotype-filter-expression", "DP < 10", "--genotype-filter-name", "New"],
            text=True, capture_output=True, check=False,
        )
        assert ft_result.returncode == 0, ft_result.stderr
        ft_lines = [line.split("\t") for line in gzip.open(
            ft_native, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert ft_lines[0][9].endswith(":New;Prior")
        assert ft_lines[0][10].endswith(":Existing2")
        ft_gatk = work / "existing-ft-gatk.vcf"
        ft_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(ft_source), "-O", str(ft_gatk),
             "--genotype-filter-expression", "DP < 10", "--genotype-filter-name", "New"],
            text=True, capture_output=True, check=False,
        )
        assert ft_gatk_result.returncode == 0, ft_gatk_result.stderr
        assert normalized_filter_records(ft_native) == normalized_filter_records(ft_gatk)

        # GVCF interval selection follows the shared record-span contract: a
        # block beginning at POS=1 and ending at END=5 must be selected by
        # -L chr1:4-4 even though its POS is outside the interval.
        gvcf_source = work / "span-input.g.vcf.gz"
        gvcf_header = HEADER.replace(
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n",
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate allele>\n")
        with gzip.open(gvcf_source, "wt", encoding="utf-8") as handle:
            handle.write(gvcf_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t20\tPASS\tEND=5;DP=5\tGT\t0/0\n")
        gvcf_output = work / "span-output.g.vcf.gz"
        gvcf_result = subprocess.run(
            [str(binary), "-V", str(gvcf_source), "-O", str(gvcf_output),
             "-L", "chr1:4-4", "--filter-expression", "DP < 10", "--filter-name", "LowDP"],
            text=True, capture_output=True, check=False,
        )
        assert gvcf_result.returncode == 0, gvcf_result.stderr
        gvcf_records = [line.split("\t") for line in gzip.open(
            gvcf_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(gvcf_records) == 1 and "LowDP" in gvcf_records[0][6]

        # Compare the supported site/genotype predicate boundary to GATK's
        # real VariantFiltration. htsjdk emits FILTER labels in lexical order,
        # while the native writer preserves rule order, so FILTER is compared
        # as a set; coordinates, QUAL, INFO, FORMAT and genotype values remain
        # exact.
        source_plain = work / "input.vcf"
        source_plain.write_text(gzip.open(source, "rt", encoding="utf-8").read(), encoding="utf-8")
        gatk_output = work / "gatk-filtered.vcf"
        native_oracle_output = work / "native-oracle-filtered.vcf.gz"
        native_oracle_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(native_oracle_output),
             "--filter-expression", "QUAL < 30", "--filter-name", "LowQual",
             "--filter-expression", "DP < 10", "--filter-name", "LowDP",
             "--filter-expression", 'vc.getGenotype("S1").isHet()',
             "--filter-name", "Het"],
            text=True, capture_output=True, check=False,
        )
        assert native_oracle_result.returncode == 0, native_oracle_result.stderr
        gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(source_plain), "-O", str(gatk_output),
             "--filter-expression", "QUAL < 30", "--filter-name", "LowQual",
             "--filter-expression", "DP < 10", "--filter-name", "LowDP",
             "--filter-expression", 'vc.getGenotype("S1").isHet()',
             "--filter-name", "Het"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_result.returncode == 0, gatk_result.stderr
        assert normalized_filter_records(native_oracle_output) == normalized_filter_records(gatk_output), (
            normalized_filter_records(native_oracle_output), normalized_filter_records(gatk_output)
        )

        # The uncompressed writer must publish a Tribble `.idx` that GATK can
        # use for a real interval traversal; `.tbi` is only valid for BGZF.
        plain_output = work / "plain-filtered.vcf"
        plain_manifest = work / "plain-filtered.manifest.json"
        plain_result = subprocess.run(
            [str(binary), "-V", str(source_plain), "-O", str(plain_output),
             "--filter-expression", "QUAL < 30", "--filter-name", "LowQual",
             "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.is_file() and Path(f"{plain_output}.idx").is_file()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "plain-filtered-query.vcf"
        plain_query_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
             "-V", str(plain_output), "-L", "chr1:1-1", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert len(normalized_filter_records(plain_query)) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True
        assert all(item["complete"] for item in plain_metadata["outputs"])

        # VariantContext site methods and string-valued accessors use the
        # same common JEXL subset as SelectVariants.
        def method_filter(expression: str, name: str, input_path: Path = source) -> list[list[str]]:
            method_output = work / f"{name}.vcf.gz"
            method_result = subprocess.run(
                [str(binary), "-V", str(input_path), "-O", str(method_output),
                 "--filter-expression", expression, "--filter-name", name],
                text=True, capture_output=True, check=False,
            )
            assert method_result.returncode == 0, method_result.stderr
            return [line.split("\t") for line in gzip.open(
                method_output, "rt", encoding="utf-8").read().splitlines()
                if line and not line.startswith("#")]

        method_records = method_filter("vc.isSNP()", "IsSNP")
        assert all("IsSNP" in record[6] for record in method_records)
        method_records = method_filter("vc.isTransition()", "IsTransition")
        assert all("IsTransition" in record[6] for record in method_records)
        method_records = method_filter("vc.isPass()", "IsPass")
        assert all("IsPass" in record[6] for record in method_records)
        method_records = method_filter('vc.getType() == "SNP"', "TypeSNP")
        assert all("TypeSNP" in record[6] for record in method_records)
        method_records = method_filter('vc.hasAttribute("QD")', "HasQD")
        assert "HasQD" in method_records[0][6] and "HasQD" in method_records[1][6]
        assert method_records[2][6] == "PASS"
        method_records = method_filter("vc.getStart() >= 2", "StartGE2")
        assert method_records[0][6] == "PASS"
        assert "StartGE2" in method_records[1][6] and "StartGE2" in method_records[2][6]
        method_records = method_filter("vc.isNotFiltered()", "IsNotFiltered")
        assert all("IsNotFiltered" in record[6] for record in method_records)
        method_records = method_filter("vc.hasAlternateAllele(0)", "HasAlt0")
        assert all("HasAlt0" in record[6] for record in method_records)
        method_records = method_filter("vc.getNAlleles() == 2", "NAlleles2")
        assert all("NAlleles2" in record[6] for record in method_records)
        method_records = method_filter("vc.getAlleles().size() == 2", "AllelesSize2")
        assert all("AllelesSize2" in record[6] for record in method_records)
        method_records = method_filter('vc.getContig().startsWith("chr")', "ContigPrefix")
        assert all("ContigPrefix" in record[6] for record in method_records)
        method_records = method_filter('vc.getID().contains("rs")', "IDContains")
        assert "IDContains" in method_records[0][6]
        assert method_records[1][6] == "PASS" and method_records[2][6] == "PASS"
        method_records = method_filter('vc.getReference().getBaseString().matches("[ACG]")', "RefRegex")
        assert all("RefRegex" in record[6] for record in method_records)
        method_records = method_filter('vc.getAttribute("DP").contains("2")', "DPString")
        assert method_records[0][6] == "PASS" and "DPString" in method_records[1][6]
        assert "DPString" in method_records[2][6]
        method_records = method_filter('vc.getAttribute("MISSING") == null', "MissingNull")
        assert all("MissingNull" in record[6] for record in method_records)
        method_records = method_filter("vc.getType() == VariantContext.Type.SNP", "EnumSNP")
        assert all("EnumSNP" in record[6] for record in method_records)
        # GATK/JEXL also permits explicit boolean comparisons and the regex
        # binary operators.  Keep these separate from the method-call tests
        # because they exercise the AST's typed RHS and string reader paths.
        bool_site_records = method_filter("vc.isSNP() == true", "BoolSNP")
        assert all("BoolSNP" in record[6] for record in bool_site_records)
        bool_genotype_records = method_filter(
            'vc.getGenotype("S1").isHet() == true', "BoolHet")
        assert "BoolHet" in bool_genotype_records[0][6]
        assert bool_genotype_records[1][6] == "PASS"
        regex_records = method_filter('vc.getID() =~ "rs.*"', "RegexID")
        assert "RegexID" in regex_records[0][6]
        assert regex_records[1][6] == "PASS" and regex_records[2][6] == "PASS"
        regex_not_records = method_filter('vc.getID() !~ "rs.*"', "RegexNotID")
        assert regex_not_records[0][6] == "PASS"
        assert "RegexNotID" in regex_not_records[1][6] and "RegexNotID" in regex_not_records[2][6]
        # Compare all four forms to the real GATK writer; FILTER order is
        # normalized by normalized_filter_records below.
        for expression, name in (
            ("vc.isSNP() == true", "BoolSNPOracle"),
            ('vc.getGenotype("S1").isHet() == true', "BoolHetOracle"),
            ('vc.getID() =~ "rs.*"', "RegexIDOracle"),
            ('vc.getID() !~ "rs.*"', "RegexNotIDOracle"),
        ):
            native_oracle = work / f"{name}-native.vcf.gz"
            native_result = subprocess.run(
                [str(binary), "-V", str(source), "-O", str(native_oracle),
                 "--filter-expression", expression, "--filter-name", name],
                text=True, capture_output=True, check=False)
            assert native_result.returncode == 0, native_result.stderr
            gatk_oracle = work / f"{name}-gatk.vcf"
            gatk_result = subprocess.run(
                [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
                 "-V", str(source_plain), "-O", str(gatk_oracle),
                 "--filter-expression", expression, "--filter-name", name],
                text=True, capture_output=True, check=False)
            assert gatk_result.returncode == 0, gatk_result.stderr
            assert normalized_filter_records(native_oracle) == normalized_filter_records(gatk_oracle)
        method_records = method_filter("vc.isPolymorphicInSamples()", "Polymorphic")
        assert "Polymorphic" in method_records[0][6] and method_records[1][6] == "PASS" and "Polymorphic" in method_records[2][6]
        method_records = method_filter("vc.isMonomorphicInSamples()", "Monomorphic")
        assert "Monomorphic" in method_records[1][6] and method_records[0][6] == "PASS"
        method_records = method_filter("vc.getReference().isReference()", "ReferenceAllele")
        assert all("ReferenceAllele" in record[6] for record in method_records)
        method_records = method_filter("vc.getAlternateAllele(0).isNonReference()", "NonReferenceAllele")
        assert all("NonReferenceAllele" in record[6] for record in method_records)
        method_records = method_filter("vc.getAlternateAllele(0).isCalled()", "CalledAllele")
        assert all("CalledAllele" in record[6] for record in method_records)

        # INFO vectors use HTSJDK's zero-based JEXL indexing.  Keep the
        # explicit and compact spellings on an oracle fixture so a parser
        # regression cannot silently fall back to element zero.
        info_array_source = work / "info-array.vcf.gz"
        info_array_header = HEADER.replace(
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n",
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
            "##INFO=<ID=VEC_INT,Number=A,Type=Integer,Description=Per-ALT score>\n"
            "##INFO=<ID=VEC_FLOAT,Number=A,Type=Float,Description=Per-ALT frequency>\n")
        with gzip.open(info_array_source, "wt", encoding="utf-8") as handle:
            handle.write(info_array_header)
            handle.write("chr1\t1\trs1\tA\tG,T\t50\tPASS\tVEC_INT=11,22;VEC_FLOAT=0.1,0.2\tGT\t0/1\n")
            handle.write("chr1\t2\t.\tC\tT\t50\tPASS\tVEC_INT=5;VEC_FLOAT=0.05\tGT\t0/1\n")
        info_array_plain = work / "info-array.vcf"
        info_array_plain.write_text(gzip.open(info_array_source, "rt", encoding="utf-8").read(),
                                    encoding="utf-8")
        info_array_native = work / "info-array-native.vcf.gz"
        info_array_result = subprocess.run(
            [str(binary), "-V", str(info_array_source), "-O", str(info_array_native),
             "--filter-expression", 'vc.getAttribute("VEC_FLOAT").get(1) >= 0.2',
             "--filter-name", "ArrayHi"], text=True, capture_output=True, check=False)
        assert info_array_result.returncode == 0, info_array_result.stderr
        info_array_gatk = work / "info-array-gatk.vcf"
        info_array_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(info_array_plain), "-O", str(info_array_gatk),
             "--filter-expression", 'vc.getAttribute("VEC_FLOAT").get(1) >= 0.2',
             "--filter-name", "ArrayHi"], text=True, capture_output=True, check=False)
        assert info_array_gatk_result.returncode == 0, info_array_gatk_result.stderr
        assert normalized_filter_records(info_array_native) == normalized_filter_records(info_array_gatk)
        # HTSJDK's historical dot-index shorthand is still used by existing
        # production command lines (for example getAttribute("AF").1).
        dot_array_native = work / "info-array-dot-native.vcf.gz"
        dot_array_result = subprocess.run(
            [str(binary), "-V", str(info_array_source), "-O", str(dot_array_native),
             "--filter-expression", 'vc.getAttribute("VEC_FLOAT").1 >= 0.2',
             "--filter-name", "ArrayDot"], text=True, capture_output=True, check=False)
        assert dot_array_result.returncode == 0, dot_array_result.stderr
        dot_array_gatk = work / "info-array-dot-gatk.vcf"
        dot_array_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(info_array_plain), "-O", str(dot_array_gatk),
             "--filter-expression", 'vc.getAttribute("VEC_FLOAT").1 >= 0.2',
             "--filter-name", "ArrayDot"], text=True, capture_output=True, check=False)
        assert dot_array_gatk_result.returncode == 0, dot_array_gatk_result.stderr
        assert normalized_filter_records(dot_array_native) == normalized_filter_records(dot_array_gatk)
        compact_array_native = work / "info-array-compact-native.vcf.gz"
        compact_array_result = subprocess.run(
            [str(binary), "-V", str(info_array_source), "-O", str(compact_array_native),
             "--filter-expression", "VEC_INT[0] >= 11", "--filter-name", "ArrayCompact"],
            text=True, capture_output=True, check=False)
        assert compact_array_result.returncode == 0, compact_array_result.stderr
        compact_array_records = [line.split("\t") for line in gzip.open(
            compact_array_native, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "ArrayCompact" in compact_array_records[0][6]
        assert compact_array_records[1][6] == "PASS"

        # JEXL arithmetic keeps normal Java precedence and supports both
        # scalar fields and indexed INFO vectors.  Compare the compound
        # operands to GATK so a parser regression cannot silently treat the
        # expression as a literal field name.
        arithmetic_cases = [
            ('QUAL / 2 >= 25', "ArithmeticQual"),
            ('(QUAL - 10) / 2 >= 20', "ArithmeticParen"),
            ('vc.getAttribute("VEC_FLOAT").1 * 10 >= 2', "ArithmeticVector"),
        ]
        for expression, name in arithmetic_cases:
            native_path = work / f"{name}-native.vcf.gz"
            native_result = subprocess.run(
                [str(binary), "-V", str(info_array_source), "-O", str(native_path),
                 "--filter-expression", expression, "--filter-name", name],
                text=True, capture_output=True, check=False)
            assert native_result.returncode == 0, native_result.stderr
            gatk_path = work / f"{name}-gatk.vcf"
            gatk_result = subprocess.run(
                [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
                 "-V", str(info_array_plain), "-O", str(gatk_path),
                 "--filter-expression", expression, "--filter-name", name],
                text=True, capture_output=True, check=False)
            assert gatk_result.returncode == 0, gatk_result.stderr
            assert normalized_filter_records(native_path) == normalized_filter_records(gatk_path)

        # HTSJDK Allele accessors are a frequent production JEXL boundary:
        # symbolic ALT detection and reference/ALT length comparisons must
        # use the zero-based ALT index, including symbolic and indel records.
        allele_method_source = work / "allele-methods.vcf.gz"
        allele_method_header = HEADER.replace(
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n",
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
            "##ALT=<ID=DEL,Description=Deletion>\n")
        with gzip.open(allele_method_source, "wt", encoding="utf-8") as handle:
            handle.write(allele_method_header)
            handle.write("chr1\t1\t.\tA\tG\t50\tPASS\tDP=10;QD=5\tGT\t0/1\n")
            handle.write("chr1\t2\t.\tA\t<DEL>\t50\tPASS\tDP=10;QD=5\tGT\t0/1\n")
            handle.write("chr1\t3\t.\tAT\tA\t50\tPASS\tDP=10;QD=5\tGT\t0/1\n")
            handle.write("chr1\t4\t.\tAT\tGC\t50\tPASS\tDP=10;QD=5\tGT\t0/1\n")
        allele_method_output = work / "allele-methods-native.vcf.gz"
        allele_method_result = subprocess.run(
            [str(binary), "-V", str(allele_method_source), "-O", str(allele_method_output),
             "--filter-expression", "vc.getAlternateAllele(0).isSymbolic()",
             "--filter-name", "AltSymbolic",
             "--filter-expression", "vc.getAlternateAllele(0).length() > 1",
             "--filter-name", "AltLong",
             "--filter-expression", "vc.getReference().length() > 1",
             "--filter-name", "RefLong"],
            text=True, capture_output=True, check=False,
        )
        assert allele_method_result.returncode == 0, allele_method_result.stderr
        allele_method_records = [line.split("\t") for line in gzip.open(
            allele_method_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert allele_method_records[0][6] == "PASS"
        assert "AltSymbolic" in allele_method_records[1][6] and "AltLong" not in allele_method_records[1][6]
        assert "RefLong" in allele_method_records[2][6]
        assert "AltLong" in allele_method_records[3][6] and "RefLong" in allele_method_records[3][6]
        allele_method_plain = work / "allele-methods.vcf"
        allele_method_plain.write_text(
            gzip.open(allele_method_source, "rt", encoding="utf-8").read(), encoding="utf-8")
        allele_method_gatk = work / "allele-methods-gatk.vcf"
        allele_method_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(allele_method_plain), "-O", str(allele_method_gatk),
             "--filter-expression", "vc.getAlternateAllele(0).isSymbolic()",
             "--filter-name", "AltSymbolic",
             "--filter-expression", "vc.getAlternateAllele(0).length() > 1",
             "--filter-name", "AltLong",
             "--filter-expression", "vc.getReference().length() > 1",
             "--filter-name", "RefLong"],
            text=True, capture_output=True, check=False,
        )
        assert allele_method_gatk_result.returncode == 0, allele_method_gatk_result.stderr
        assert normalized_filter_records(allele_method_output) == normalized_filter_records(allele_method_gatk)

        breakend_source = work / "breakend-methods.vcf.gz"
        with gzip.open(breakend_source, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t5\t.\tA\tN]chr1:10]\t50\tPASS\tDP=10\tGT\t0/1\n")
            handle.write("chr1\t6\t.\tA\t.N\t50\tPASS\tDP=10\tGT\t0/1\n")
        breakend_output = work / "breakend-methods-native.vcf.gz"
        breakend_result = subprocess.run(
            [str(binary), "-V", str(breakend_source), "-O", str(breakend_output),
             "--filter-expression", "vc.getAlternateAllele(0).isBreakpoint()",
             "--filter-name", "Breakpoint",
             "--filter-expression", "vc.getAlternateAllele(0).isSingleBreakend()",
             "--filter-name", "SingleBreakend",
             "--filter-expression", "vc.getAlternateAllele(0).isSymbolic()",
             "--filter-name", "Symbolic"],
            text=True, capture_output=True, check=False,
        )
        assert breakend_result.returncode == 0, breakend_result.stderr
        breakend_records = [line.split("\t") for line in gzip.open(
            breakend_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "Breakpoint" in breakend_records[0][6] and "SingleBreakend" not in breakend_records[0][6]
        assert "SingleBreakend" in breakend_records[1][6] and "Breakpoint" not in breakend_records[1][6]
        breakend_plain = work / "breakend-methods.vcf"
        breakend_plain.write_text(gzip.open(breakend_source, "rt", encoding="utf-8").read(), encoding="utf-8")
        breakend_gatk = work / "breakend-methods-gatk.vcf"
        breakend_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(breakend_plain), "-O", str(breakend_gatk),
             "--filter-expression", "vc.getAlternateAllele(0).isBreakpoint()",
             "--filter-name", "Breakpoint",
             "--filter-expression", "vc.getAlternateAllele(0).isSingleBreakend()",
             "--filter-name", "SingleBreakend",
             "--filter-expression", "vc.getAlternateAllele(0).isSymbolic()",
             "--filter-name", "Symbolic"],
            text=True, capture_output=True, check=False,
        )
        assert breakend_gatk_result.returncode == 0, breakend_gatk_result.stderr
        assert normalized_filter_records(breakend_output) == normalized_filter_records(breakend_gatk)

        # Common allele-specific (Number=A) annotations are lowered to the
        # per-ALT AS_FilterStatus INFO vector.  AS rules do not alter the site
        # FILTER column, which has no unambiguous ALT index.
        allele_source = work / "allele-specific.vcf.gz"
        allele_header = HEADER.replace(
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n",
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
            "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=Allele quality by depth>\n")
        with gzip.open(allele_source, "wt", encoding="utf-8") as handle:
            handle.write(allele_header)
            handle.write("chr1\t1\t.\tA\tC,G\t50\tPASS\tAS_QD=1.0,3.0\tGT\t0/1\n")
            handle.write("chr1\t2\t.\tC\tT,G\t50\tPASS\tAS_QD=3.0,1.0\tGT\t0/1\n")
        allele_output = work / "allele-specific-filtered.vcf.gz"
        allele_manifest = work / "allele-specific-filtered.manifest.json"
        allele_result = subprocess.run(
            [str(binary), "-V", str(allele_source), "-O", str(allele_output),
             "--filter-expression", 'vc.getAttribute("AS_QD") < 2',
             "--filter-name", "LowASQD", "--apply-allele-specific-filters",
             "--output-manifest", str(allele_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert allele_result.returncode == 0, allele_result.stderr
        allele_records = [line.split("\t") for line in gzip.open(
            allele_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert allele_records[0][6] == "PASS" and allele_records[1][6] == "PASS"
        assert "AS_FilterStatus=LowASQD,PASS" in allele_records[0][7]
        assert "AS_FilterStatus=PASS,LowASQD" in allele_records[1][7]
        allele_metadata = json.loads(allele_manifest.read_text(encoding="utf-8"))
        assert allele_metadata["compatibility"]["allele_specific_filters"] is True
        assert allele_metadata["telemetry"]["allele_filtered"] == 2

        unsupported_output = work / "allele-specific-unsupported.vcf.gz"
        unsupported_result = subprocess.run(
            [str(binary), "-V", str(allele_source), "-O", str(unsupported_output),
             "--filter-expression", 'vc.getAttribute("AS_QD") < 2',
             "--filter-name", "LowASQD"],
            text=True, capture_output=True, check=False,
        )
        assert unsupported_result.returncode != 0
        assert "apply-allele-specific-filters" in unsupported_result.stderr

        genotype_method_source = work / "genotype-methods.vcf.gz"
        genotype_method_header = HEADER.replace(
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n")
        with gzip.open(genotype_method_source, "wt", encoding="utf-8") as handle:
            handle.write(genotype_method_header)
            handle.write("chr1\t1\t.\tA\tG\t20\tPASS\tDP=5\tGT:AD:PL:GQ\t0/1:5,7:50,0,50:50\n")
            handle.write("chr1\t2\t.\tC\tT\t50\tPASS\tDP=20\tGT:AD:PL:GQ\t0/0:10,0:0,30,60:10\n")
            handle.write("chr1\t3\t.\tG\tA\t50\tPASS\tDP=20\tGT:AD:PL:GQ\t1/1:0,20:80,40,0:60\n")
        ad_records = method_filter('vc.getGenotype("S1").getAD()[1] >= 7',
                                    "GTAD", genotype_method_source)
        assert "GTAD" in ad_records[0][6] and "GTAD" in ad_records[2][6]
        assert ad_records[1][6] == "PASS"
        pl_records = method_filter('vc.getGenotype("S1").getPL()[2] < 10',
                                   "GTPL", genotype_method_source)
        assert pl_records[2][6].startswith("GTPL")
        dot_ad_records = method_filter('vc.getGenotype("S1").getAD().1 >= 7',
                                       "GTDotAD", genotype_method_source)
        assert "GTDotAD" in dot_ad_records[0][6] and "GTDotAD" in dot_ad_records[2][6]
        dot_pl_records = method_filter('vc.getGenotype("S1").getPL().2 < 10',
                                       "GTDotPL", genotype_method_source)
        assert dot_pl_records[2][6].startswith("GTDotPL")
        hom_records = method_filter('vc.getGenotype("S1").isHom()',
                                    "GTHOM", genotype_method_source)
        assert "GTHOM" in hom_records[1][6] and "GTHOM" in hom_records[2][6]
        assert hom_records[0][6] == "PASS"
        available_records = method_filter('vc.getGenotype("S1").isAvailable()',
                                          "GTAVAILABLE", genotype_method_source)
        assert all("GTAVAILABLE" in record[6] for record in available_records)
        called_records = method_filter('vc.getGenotype("S1").isCalled()',
                                       "GTCALLED", genotype_method_source)
        assert all("GTCALLED" in record[6] for record in called_records)
        for method, name in (("hasAD()", "GTHASAD"), ("hasPL()", "GTHASPL"), ("hasGQ()", "GTHASGQ")):
            records = method_filter(f'vc.getGenotype("S1").{method}', name, genotype_method_source)
            assert all(name in record[6] for record in records)
        ploidy_records = method_filter('vc.getGenotype("S1").getPloidy() == 2',
                                       "GTPloidy", genotype_method_source)
        assert all("GTPloidy" in record[6] for record in ploidy_records)
        phred_records = method_filter('vc.getGenotype("S1").getPhredScaledQual() >= 50',
                                      "GTPhred", genotype_method_source)
        assert "GTPhred" in phred_records[0][6] and phred_records[1][6] == "PASS" and "GTPhred" in phred_records[2][6]
        min_dp_source = work / "min-dp-genotype-filter.vcf.gz"
        min_dp_header = genotype_method_header.replace(
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n",
            "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum depth>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n")
        with gzip.open(min_dp_source, "wt", encoding="utf-8") as handle:
            handle.write(min_dp_header)
            handle.write("chr1\t1\t.\tA\tG\t50\tPASS\tDP=5\tGT:MIN_DP:AD:PL:GQ\t0/1:5:5,7:50,0,50:50\n")
            handle.write("chr1\t2\t.\tC\tT\t50\tPASS\tDP=20\tGT:MIN_DP:AD:PL:GQ\t0/0:20:10,0:0,30,60:20\n")
        min_dp_records = method_filter('vc.getGenotype("S1").getMIN_DP() < 10',
                                        "GTMinDP", min_dp_source)
        assert "GTMinDP" in min_dp_records[0][6] and min_dp_records[1][6] == "PASS"
        method_records = method_filter("vc.hasGenotypes()", "HasGenotypes", genotype_method_source)
        assert all("HasGenotypes" in record[6] for record in method_records)
        filters_empty_records = method_filter("vc.getFilters().isEmpty()", "FiltersEmpty", genotype_method_source)
        assert all("FiltersEmpty" in record[6] for record in filters_empty_records)
        genotype_size_records = method_filter("vc.getGenotypes().size() == 1", "GenotypesSize", genotype_method_source)
        assert all("GenotypesSize" in record[6] for record in genotype_size_records)
        method_records = method_filter("vc.getCalledChrCount() == 2", "CalledChr", genotype_method_source)
        assert all("CalledChr" in record[6] for record in method_records)
        method_records = method_filter("vc.getNoCallCount() == 0", "NoCall", genotype_method_source)
        assert all("NoCall" in record[6] for record in method_records)
        method_records = method_filter("vc.getHomRefCount() == 1", "HomRef", genotype_method_source)
        assert method_records[0][6] == "PASS" and "HomRef" in method_records[1][6] and method_records[2][6] == "PASS"
        method_records = method_filter("vc.getHetCount() == 1", "HetCount", genotype_method_source)
        assert "HetCount" in method_records[0][6]
        method_records = method_filter("vc.getHomVarCount() == 1", "HomVar", genotype_method_source)
        assert "HomVar" in method_records[2][6]

        # GATK's genotype-filter surface also accepts compact per-sample JEXL
        # without an explicit vc.getGenotype("sample") receiver.  The native
        # AST must bind these expressions to each sample, not interpret DP as
        # a site INFO field or silently skip the filter.
        compact_genotype_cases = [
            ("GQ < 20", "CompactGQ", 1),
            ("AD[1] == 0", "CompactAD", 1),
            ("AD.1 == 0", "CompactDotAD", 1),
            ("isHet == 1", "CompactHet", 0),
        ]
        genotype_method_plain = work / "genotype-methods.vcf"
        genotype_method_plain.write_text(
            gzip.open(genotype_method_source, "rt", encoding="utf-8").read(), encoding="utf-8")
        for expression, name, expected_index in compact_genotype_cases:
            compact_output = work / f"{name}.vcf.gz"
            compact_result = subprocess.run(
                [str(binary), "-V", str(genotype_method_source), "-O", str(compact_output),
                 "--genotype-filter-expression", expression,
                 "--genotype-filter-name", name],
                text=True, capture_output=True, check=False,
            )
            assert compact_result.returncode == 0, compact_result.stderr
            compact_records = [line.split("\t") for line in gzip.open(
                compact_output, "rt", encoding="utf-8").read().splitlines()
                if line and not line.startswith("#")]
            assert name in compact_records[expected_index][9]
            assert all(name not in compact_records[index][9]
                       for index in range(len(compact_records)) if index != expected_index)

            compact_gatk = work / f"{name}-gatk.vcf"
            compact_gatk_result = subprocess.run(
                [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
                 "-V", str(genotype_method_plain), "-O", str(compact_gatk),
                 "--genotype-filter-expression", expression,
                 "--genotype-filter-name", name],
                text=True, capture_output=True, check=False,
            )
            assert compact_gatk_result.returncode == 0, compact_gatk_result.stderr
            gatk_compact_records = [line.split("\t") for line in compact_gatk.read_text(
                encoding="utf-8").splitlines() if line and not line.startswith("#")]
            def ft_values(records: list[list[str]]) -> list[str]:
                result = []
                for record in records:
                    format_fields = record[8].split(":")
                    if "FT" not in format_fields:
                        result.append("PASS")
                    else:
                        result.append(record[9].split(":")[format_fields.index("FT")])
                return result
            native_ft = ft_values(compact_records)
            gatk_ft = ft_values(gatk_compact_records)
            assert native_ft == gatk_ft, (expression, native_ft, gatk_ft)

        arithmetic_genotype_output = work / "CompactArithmeticGQ.vcf.gz"
        arithmetic_genotype_result = subprocess.run(
            [str(binary), "-V", str(genotype_method_source), "-O", str(arithmetic_genotype_output),
             "--genotype-filter-expression", "GQ + 1 >= 61",
             "--genotype-filter-name", "CompactArithmeticGQ"],
            text=True, capture_output=True, check=False)
        assert arithmetic_genotype_result.returncode == 0, arithmetic_genotype_result.stderr
        arithmetic_genotype_records = [line.split("\t") for line in gzip.open(
            arithmetic_genotype_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "CompactArithmeticGQ" in arithmetic_genotype_records[2][9]
        assert all("CompactArithmeticGQ" not in arithmetic_genotype_records[index][9]
                   for index in (0, 1))
        arithmetic_genotype_gatk = work / "CompactArithmeticGQ-gatk.vcf"
        arithmetic_genotype_gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantFiltration",
             "-V", str(genotype_method_plain), "-O", str(arithmetic_genotype_gatk),
             "--genotype-filter-expression", "GQ + 1 >= 61",
             "--genotype-filter-name", "CompactArithmeticGQ"],
            text=True, capture_output=True, check=False)
        assert arithmetic_genotype_gatk_result.returncode == 0, arithmetic_genotype_gatk_result.stderr
        assert ft_values(arithmetic_genotype_records) == ft_values([
            line.split("\t") for line in arithmetic_genotype_gatk.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")])

        compact_min_dp_output = work / "CompactMinDP.vcf.gz"
        compact_min_dp_result = subprocess.run(
            [str(binary), "-V", str(min_dp_source), "-O", str(compact_min_dp_output),
             "--genotype-filter-expression", "MIN_DP < 10",
             "--genotype-filter-name", "CompactMinDP"],
            text=True, capture_output=True, check=False)
        assert compact_min_dp_result.returncode == 0, compact_min_dp_result.stderr
        compact_min_dp_records = [line.split("\t") for line in gzip.open(
            compact_min_dp_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "CompactMinDP" in compact_min_dp_records[0][9]
        assert compact_min_dp_records[1][9].endswith("PASS")

        interval_list = work / "filtration.interval_list"
        interval_list.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t1\t+\tfirst\n"
            "chr1\t3\t3\t+\tthird\n", encoding="utf-8")
        interval_output = work / "interval-filtered.vcf.gz"
        interval_manifest = work / "interval-filtered.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(interval_output),
             "-L", str(interval_list), "--filter-expression", "QUAL < 30",
             "--filter-name", "LowQual", "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False)
        assert interval_result.returncode == 0, interval_result.stderr
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["1", "3"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_metadata["telemetry"]["interval_list_records"] == 2

        # Genotype-level predicates are written to FORMAT/FT without
        # overwriting the site FILTER column.
        genotype_output = work / "genotype-filtered.vcf.gz"
        genotype_manifest = work / "genotype-filtered.manifest.json"
        genotype_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(genotype_output),
             "--genotype-filter-expression", 'vc.getGenotype("S1").isHet()',
             "--genotype-filter-name", "GT_HET",
             "--genotype-filter-expression", 'vc.getGenotype("S1").isHet()',
             "--genotype-filter-name", "GT_HET_2", "--output-manifest", str(genotype_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert genotype_result.returncode == 0, genotype_result.stderr
        genotype_text = gzip.open(genotype_output, "rt", encoding="utf-8").read()
        genotype_records = [line.split("\t") for line in genotype_text.splitlines()
                            if line and not line.startswith("#")]
        assert len(genotype_records) == 3
        genotype_format = genotype_records[0][8].split(":")
        assert "FT" in genotype_format
        assert genotype_records[0][9].split(":")[genotype_format.index("FT")] == "GT_HET;GT_HET_2"
        assert genotype_records[1][9].split(":")[genotype_format.index("FT")] == "PASS"
        genotype_metadata = json.loads(genotype_manifest.read_text(encoding="utf-8"))
        assert genotype_metadata["compatibility"]["genotype_ft"] is True

        invert_output = work / "invert.vcf.gz"
        invert_manifest = work / "invert.manifest.json"
        invert_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(invert_output),
             "--filter-expression", "QUAL < 30", "--filter-name", "NotLowQual",
             "--invert-filter-expression", "--output-manifest", str(invert_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert invert_result.returncode == 0, invert_result.stderr
        invert_records = [line.split("\t") for line in gzip.open(
            invert_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert invert_records[0][6] == "PASS"
        assert "NotLowQual" in invert_records[1][6] and "NotLowQual" in invert_records[2][6]
        invert_metadata = json.loads(invert_manifest.read_text(encoding="utf-8"))
        assert invert_metadata["compatibility"]["invert_filter_expression"] is True

        invert_gt_output = work / "invert-genotype.vcf.gz"
        invert_gt_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(invert_gt_output),
             "--genotype-filter-expression", 'vc.getGenotype("S1").isHet()',
             "--genotype-filter-name", "NOT_HET", "--invert-genotype-filter-expression"],
            text=True, capture_output=True, check=False,
        )
        assert invert_gt_result.returncode == 0, invert_gt_result.stderr
        invert_gt_records = [line.split("\t") for line in gzip.open(
            invert_gt_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        invert_gt_format = invert_gt_records[0][8].split(":")
        assert invert_gt_records[0][9].split(":")[invert_gt_format.index("FT")] == "PASS"
        assert invert_gt_records[1][9].split(":")[invert_gt_format.index("FT")] == "NOT_HET"

        mask = work / "mask.vcf.gz"
        with gzip.open(mask, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\n")
        mask_output = work / "masked.vcf.gz"
        mask_manifest = work / "masked.manifest.json"
        mask_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(mask_output),
             "--mask", str(mask), "--mask-name", "InMask",
             "--output-manifest", str(mask_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert mask_result.returncode == 0, mask_result.stderr
        mask_records = [line.split("\t") for line in gzip.open(
            mask_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert mask_records[0][6] == "PASS"
        assert "InMask" in mask_records[1][6]
        assert mask_records[2][6] == "PASS"
        mask_metadata = json.loads(mask_manifest.read_text(encoding="utf-8"))
        assert mask_metadata["compatibility"]["mask_filter"] is True
        assert mask_metadata["telemetry"]["masked_records"] == 1

        # The input side of --mask also follows the gVCF record-span contract:
        # a block at POS=1 with END=5 must overlap a point mask at POS=4.
        gvcf_mask_source = work / "gvcf-mask-input.vcf.gz"
        gvcf_mask_header = HEADER.replace(
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n",
            "##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate allele>\n")
        with gzip.open(gvcf_mask_source, "wt", encoding="utf-8") as handle:
            handle.write(gvcf_mask_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t20\tPASS\tEND=5;DP=20\tGT\t0/0\n")
        gvcf_mask = work / "gvcf-mask.vcf.gz"
        with gzip.open(gvcf_mask, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t4\t.\tT\tC\t.\tPASS\t.\tGT\t0/0\n")
        gvcf_mask_output = work / "gvcf-masked.vcf.gz"
        gvcf_mask_result = subprocess.run(
            [str(binary), "-V", str(gvcf_mask_source), "-O", str(gvcf_mask_output),
             "--mask", str(gvcf_mask), "--mask-name", "InMask"],
            text=True, capture_output=True, check=False,
        )
        assert gvcf_mask_result.returncode == 0, gvcf_mask_result.stderr
        gvcf_mask_records = [line.split("\t") for line in gzip.open(
            gvcf_mask_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(gvcf_mask_records) == 1 and "InMask" in gvcf_mask_records[0][6]

        cluster_source = work / "cluster-input.vcf.gz"
        with gzip.open(cluster_source, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write("chr1\t1\t.\tA\tG\t20\tPASS\tDP=5;QD=1.0\tGT\t0/1\n")
            handle.write("chr1\t2\t.\tC\tT\t50\tPASS\tDP=20;QD=3.0\tGT\t0/0\n")
            handle.write("chr1\t10\t.\tG\tA\t50\tPASS\tDP=20\tGT\t0/1\n")
        cluster_output = work / "clustered.vcf.gz"
        cluster_manifest = work / "clustered.manifest.json"
        cluster_result = subprocess.run(
            [str(binary), "-V", str(cluster_source), "-O", str(cluster_output),
             "--cluster-size", "2", "--cluster-window-size", "1",
             "--output-manifest", str(cluster_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert cluster_result.returncode == 0, cluster_result.stderr
        cluster_records = [line.split("\t") for line in gzip.open(
            cluster_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "ClusteredEvents" in cluster_records[0][6]
        assert "ClusteredEvents" in cluster_records[1][6]
        assert cluster_records[2][6] == "PASS"
        cluster_metadata = json.loads(cluster_manifest.read_text(encoding="utf-8"))
        assert cluster_metadata["compatibility"]["cluster_filter"] is True
        assert cluster_metadata["telemetry"]["clustered_records"] == 2

        interval_output = work / "interval-filtered.vcf.gz"
        interval_manifest = work / "interval-filtered.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-V", str(source), "-O", str(interval_output),
             "-L", "chr1:1-1", "--intervals", "chr1:3-3",
             "--filter-expression", "QUAL < 30", "--filter-name", "OnlyLowQual",
             "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert interval_result.returncode == 0, interval_result.stderr
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["1", "3"]
        assert "OnlyLowQual" in interval_records[0][6]
        assert interval_records[1][6] == "PASS"
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["interval_skipped"] == 1

    print(json.dumps({"status": "pass", "filtered_records": summary["filtered_records"],
                      "gatk_oracle_exact": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
