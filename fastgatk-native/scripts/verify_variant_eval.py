#!/usr/bin/env python3
"""Verify the native VariantEval count/TiTv/comparison compatibility contract."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##contig=<ID=chr2,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""


def metric(report: Path, table: str, key: str) -> str:
    lines = report.read_text(encoding="utf-8").splitlines()
    active = False
    for line in lines:
        if line == f"##table={table}":
            active = True
            continue
        if active and line.startswith("##table="):
            active = False
        if active and line.startswith(f"{key}\t"):
            return line.split("\t", 1)[1]
    raise AssertionError(f"missing {table}.{key}")


def stratum_row(report: Path, table: str, name: str) -> list[str]:
    lines = report.read_text(encoding="utf-8").splitlines()
    active = False
    for line in lines:
        if line == f"##table={table}":
            active = True
            continue
        if active and line.startswith("##table="):
            break
        if active and line.startswith(name + "\t"):
            return line.split("\t")
    raise AssertionError(f"missing {table} stratum {name}")


def gatk_table_rows(report: Path, table: str) -> tuple[list[str], list[list[str]]]:
    lines = report.read_text(encoding="utf-8").splitlines()
    marker = next(index for index, line in enumerate(lines)
                  if line.startswith(f"#:GATKTable:") and f":{table}:" in line)
    header = lines[marker + 1].split()
    rows: list[list[str]] = []
    index = marker + 2
    while index < len(lines) and lines[index] and not lines[index].startswith("#:"):
        rows.append(lines[index].split())
        index += 1
    return header, rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_EVAL_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-variant-eval")))
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-eval-") as directory:
        work = Path(directory)
        eval_vcf = work / "eval.vcf"
        comp_vcf = work / "comp.vcf"
        eval_vcf.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t2\trs2\tC\tA\t50\tPASS\t.\tGT\t1/1\n"
            + "chr1\t3\trs3\tA\tAT\t50\tPASS\t.\tGT\t0/1\n"
            + "chr2\t10\trs4\tG\tT\t50\tLowQual\t.\tGT\t./.\n",
            encoding="utf-8",
        )
        comp_vcf.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t2\trs2\tC\tG\t50\tPASS\t.\tGT\t1/1\n",
            encoding="utf-8",
        )
        report = work / "eval.report"
        manifest = work / "eval.manifest.json"
        result = subprocess.run(
            [
                str(binary), "-eval", str(eval_vcf), "-comp", str(comp_vcf),
                "-L", "chr1:1-3", "-O", str(report),
                "-EV", "CountVariants", "-S", "Contig",
                "--output-manifest", str(manifest),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
        assert metric(report, "CountVariants", "nProcessedLoci") == "3"
        assert metric(report, "CountVariants", "nSNPs") == "2"
        assert metric(report, "CountVariants", "nInsertions") == "1"
        assert metric(report, "VariantSummary", "nFilteredLoci") == "0"
        assert metric(report, "TiTvVariantEvaluator", "nTi") == "1"
        assert metric(report, "TiTvVariantEvaluator", "nTv") == "1"
        assert metric(report, "CompOverlap", "nEvalCompOverlap") == "2"
        assert metric(report, "CompOverlap", "nStrictAlleleMismatches") == "1"
        assert stratum_row(report, "Contig", "chr1")[1] == "3"
        assert stratum_row(report, "VariantType", "SNP")[2] == "2"
        assert stratum_row(report, "VariantType", "INSERTION")[2] == "1"
        assert metric(report, "GenotypeConcordance", "nComparedGenotypes") == "2"
        assert metric(report, "GenotypeConcordance", "nConcordantGenotypes") == "2"
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["count_variants"] is True
        assert metadata["compatibility"]["standard_stratifications"] is True
        assert metadata["compatibility"]["genotype_concordance"] is True
        assert metadata["telemetry"]["interval_list_inputs"] == 0
        assert metadata["fallback"]["unsupported_modules"] is False

        # Repeated -L selectors use GATK's explicit interval-set rule.  The
        # intersection of chr1:1-3 and chr1:3-3 leaves a single evaluation
        # record, and the selected rule is persisted in both report metadata
        # and OutputManifest telemetry.
        intersection_report = work / "intersection.report"
        intersection_manifest = work / "intersection.manifest.json"
        intersection_result = subprocess.run(
            [str(binary), "-eval", str(eval_vcf), "-O", str(intersection_report),
             "-L", "chr1:1-3", "-L", "chr1:3-3",
             "--interval-set-rule", "INTERSECTION", "-no-ev", "-EV", "CountVariants",
             "--output-manifest", str(intersection_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert intersection_result.returncode == 0, intersection_result.stderr
        assert metric(intersection_report, "CountVariants", "nProcessedLoci") == "1"
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"

        gatk_report = work / "eval.gatk.report"
        gatk_manifest = work / "eval.gatk.manifest.json"
        gatk_result = subprocess.run(
            [str(binary), "-eval", str(eval_vcf), "-comp", str(comp_vcf), "-O", str(gatk_report),
             "--gatk-report", "--output-manifest", str(gatk_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert gatk_result.returncode == 0, gatk_result.stderr
        gatk_text = gatk_report.read_text(encoding="utf-8")
        assert gatk_text.startswith("#:GATKReport.v1.1:")
        assert "#:GATKTable:CountVariants:" in gatk_text
        assert "#:GATKTable:VariantSummary:" in gatk_text
        gatk_metadata = json.loads(gatk_manifest.read_text(encoding="utf-8"))
        assert gatk_metadata["compatibility"]["gatk_report"] is True

        # Sample stratification keeps the per-sample genotype state instead of
        # collapsing every record into an ALL_SAMPLES row.  The GATK-shaped
        # CountVariants table must expose both sample rows and the aggregate.
        multi_eval = work / "multi.eval.vcf"
        multi_eval.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
17\t2\t.\tC\tT\t50\tPASS\t.\tGT\t./.\t1/1
""",
            encoding="utf-8",
        )
        multi_report = work / "multi.gatk.report"
        multi_manifest = work / "multi.manifest.json"
        multi_result = subprocess.run(
            [str(binary), "-eval", str(multi_eval), "-O", str(multi_report),
             "-ST", "Sample", "-EV", "CountVariants", "--gatk-report",
             "--output-manifest", str(multi_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert multi_result.returncode == 0, multi_result.stderr
        multi_header, multi_rows = gatk_table_rows(multi_report, "CountVariants")
        assert "Sample" in multi_header
        sample_index = multi_header.index("Sample")
        by_sample = {row[sample_index]: row for row in multi_rows}
        assert {"S1", "S2", "all"}.issubset(by_sample)
        assert by_sample["S1"][multi_header.index("nVariantLoci")] == "1"
        assert by_sample["S1"][multi_header.index("nNoCalls")] == "1"
        assert by_sample["S2"][multi_header.index("nHomRef")] == "1"
        multi_metadata = json.loads(multi_manifest.read_text(encoding="utf-8"))
        assert multi_metadata["compatibility"]["sample_stratification"] is True
        assert multi_metadata["telemetry"]["sample_names"] == 2

        # GVCF blocks use the same span/END interval contract as the other
        # VCF tools; POS=1/END=5 must be visible to 17:4-4.
        span_eval = work / "span.eval.g.vcf"
        span_eval.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Any alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
17\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT\t0/0
""",
            encoding="utf-8",
        )
        span_report = work / "span.report"
        span_result = subprocess.run(
            [str(binary), "-eval", str(span_eval), "-O", str(span_report),
             "-L", "17:4-4", "-no-ev", "-EV", "CountVariants"],
            text=True, capture_output=True, check=False,
        )
        assert span_result.returncode == 0, span_result.stderr
        assert metric(span_report, "CountVariants", "nProcessedLoci") == "1"

        # Compare the per-sample CountVariants counters with the pinned
        # GATK 4.6.2.0 report.  Header provenance/formatting is intentionally
        # not used as the oracle; the row semantics are.
        gatk_java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        assert gatk_java.is_file() and gatk_jar.is_file()
        oracle_report = work / "multi.oracle.gatk.report"
        oracle_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantEval",
             "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
             "-eval", str(multi_eval), "-O", str(oracle_report),
             "-no-ev", "-EV", "CountVariants", "-no-st", "-ST", "Sample"],
            text=True, capture_output=True, check=False,
        )
        assert oracle_result.returncode == 0, oracle_result.stderr
        oracle_header, oracle_rows = gatk_table_rows(oracle_report, "CountVariants")
        oracle_index = oracle_header.index("Sample")
        oracle_by_sample = {row[oracle_index]: row for row in oracle_rows}
        for sample in ("S1", "S2", "all"):
            assert by_sample[sample][multi_header.index("nVariantLoci")] == \
                   oracle_by_sample[sample][oracle_header.index("nVariantLoci")]
            assert by_sample[sample][multi_header.index("nNoCalls")] == \
                   oracle_by_sample[sample][oracle_header.index("nNoCalls")]

        # VariantAFEvaluator is a compact but real evaluator in the Java
        # graph: default AC=0 SNPs are ignored, each called diploid genotype
        # contributes 0/0.5/1.0 AF, and sites-only records use INFO/AF.  The
        # native GATKReport table is compared directly to the pinned oracle.
        native_af_report = work / "variant-af.native.gatk.report"
        native_af_result = subprocess.run(
            [str(binary), "-eval", str(multi_eval), "-O", str(native_af_report),
             "-no-ev", "-no-st", "-EV", "VariantAFEvaluator", "--gatk-report"],
            text=True, capture_output=True, check=False,
        )
        assert native_af_result.returncode == 0, native_af_result.stderr
        native_af_header, native_af_rows = gatk_table_rows(native_af_report, "VariantAFEvaluator")
        assert len(native_af_rows) == 1
        gatk_af_report = work / "variant-af.gatk.report"
        gatk_af_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantEval",
             "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
             "-eval", str(multi_eval), "-O", str(gatk_af_report),
             "-no-ev", "-EV", "VariantAFEvaluator", "-no-st"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_af_result.returncode == 0, gatk_af_result.stderr
        gatk_af_header, gatk_af_rows = gatk_table_rows(gatk_af_report, "VariantAFEvaluator")
        assert native_af_header == gatk_af_header
        assert len(gatk_af_rows) == 1
        for field in ("avgVarAF", "totalCalledSites", "totalHetSites",
                      "totalHomVarSites", "totalHomRefSites"):
            assert native_af_rows[0][native_af_header.index(field)] == \
                   gatk_af_rows[0][gatk_af_header.index(field)], (field, native_af_rows, gatk_af_rows)

        # ThetaVariantEvaluator uses the same polymorphic SNP stream but
        # reports per-site heterozygosity, average pairwise differences and
        # the harmonic-number region estimate.  Keep this comparison against
        # the pinned Java evaluator so the native implementation cannot drift
        # behind GATK's report contract.
        native_theta_report = work / "theta.native.gatk.report"
        native_theta_result = subprocess.run(
            [str(binary), "-eval", str(work / "multi.eval.vcf"), "-O", str(native_theta_report),
             "-no-ev", "-no-st", "-EV", "ThetaVariantEvaluator", "--gatk-report"],
            text=True, capture_output=True, check=False,
        )
        assert native_theta_result.returncode == 0, native_theta_result.stderr
        native_theta_header, native_theta_rows = gatk_table_rows(
            native_theta_report, "ThetaVariantEvaluator")
        assert len(native_theta_rows) == 1
        gatk_theta_report = work / "theta.gatk.report"
        gatk_theta_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantEval",
             "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
             "-eval", str(work / "multi.eval.vcf"), "-O", str(gatk_theta_report),
             "-no-ev", "-EV", "ThetaVariantEvaluator", "-no-st"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_theta_result.returncode == 0, gatk_theta_result.stderr
        gatk_theta_header, gatk_theta_rows = gatk_table_rows(
            gatk_theta_report, "ThetaVariantEvaluator")
        assert native_theta_header == gatk_theta_header
        assert len(gatk_theta_rows) == 1
        for field in ("avgHet", "avgAvgDiffs", "totalHet", "totalAvgDiffs",
                      "thetaRegionNumSites"):
            assert native_theta_rows[0][native_theta_header.index(field)] == \
                   gatk_theta_rows[0][gatk_theta_header.index(field)], \
                   (field, native_theta_rows, gatk_theta_rows)

        # MendelianViolationEvaluator consumes the standard six-column PED
        # relationship file and the default GATK mvq=50 threshold.  Compare
        # every one of its 31 report fields against the Java evaluator.
        mendel_eval = work / "mendel.eval.vcf"
        mendel_ped = work / "mendel.ped"
        mendel_eval.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype Quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tC
17\t1\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/0:99\t0/0:99\t0/1:99
17\t2\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/1:99\t0/1:99\t1/1:99
17\t3\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t1/1:99\t1/1:99\t0/0:99
""",
            encoding="utf-8",
        )
        mendel_ped.write_text(
            "fam F 0 0 1 1\nfam M 0 0 2 1\nfam C F M 1 1\n",
            encoding="utf-8",
        )
        native_mendel_report = work / "mendel.native.gatk.report"
        native_mendel_manifest = work / "mendel.native.manifest.json"
        native_mendel_result = subprocess.run(
            [str(binary), "-eval", str(mendel_eval), "--pedigree", str(mendel_ped),
             "-O", str(native_mendel_report), "-no-ev", "-no-st",
             "-EV", "MendelianViolationEvaluator", "--gatk-report",
             "--output-manifest", str(native_mendel_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert native_mendel_result.returncode == 0, native_mendel_result.stderr
        native_mendel_header, native_mendel_rows = gatk_table_rows(
            native_mendel_report, "MendelianViolationEvaluator")
        assert len(native_mendel_rows) == 1
        gatk_mendel_report = work / "mendel.gatk.report"
        gatk_mendel_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantEval",
             "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
             "-eval", str(mendel_eval), "--pedigree", str(mendel_ped),
             "-O", str(gatk_mendel_report), "-no-ev",
             "-EV", "MendelianViolationEvaluator", "-no-st"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_mendel_result.returncode == 0, gatk_mendel_result.stderr
        gatk_mendel_header, gatk_mendel_rows = gatk_table_rows(
            gatk_mendel_report, "MendelianViolationEvaluator")
        assert native_mendel_header == gatk_mendel_header
        assert len(gatk_mendel_rows) == 1
        assert native_mendel_rows[0] == gatk_mendel_rows[0], \
               (native_mendel_rows, gatk_mendel_rows)
        mendel_metadata = json.loads(native_mendel_manifest.read_text(encoding="utf-8"))
        assert mendel_metadata["compatibility"]["mendelian_violation_evaluator"] is True
        assert mendel_metadata["compatibility"]["pedigree_trios"] == 1
        assert mendel_metadata["telemetry"]["mendelian_n_violations"] == 2

        # Family stratification subsets CountVariants to each pedigree family
        # and exposes the aggregate "all" row.  Rate columns are compared
        # numerically because GATK uses per-column printf formats.
        native_family_report = work / "family.native.gatk.report"
        native_family_manifest = work / "family.native.manifest.json"
        native_family_result = subprocess.run(
            [str(binary), "-eval", str(mendel_eval), "--pedigree", str(mendel_ped),
             "-O", str(native_family_report), "-no-ev", "-EV", "CountVariants",
             "-no-st", "-ST", "Family", "--gatk-report",
             "--output-manifest", str(native_family_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert native_family_result.returncode == 0, native_family_result.stderr
        native_family_header, native_family_rows = gatk_table_rows(
            native_family_report, "CountVariants")
        gatk_family_report = work / "family.gatk.report"
        gatk_family_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "VariantEval",
             "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
             "-eval", str(mendel_eval), "--pedigree", str(mendel_ped),
             "-O", str(gatk_family_report), "-no-ev", "-EV", "CountVariants",
             "-no-st", "-ST", "Family"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_family_result.returncode == 0, gatk_family_result.stderr
        gatk_family_header, gatk_family_rows = gatk_table_rows(
            gatk_family_report, "CountVariants")
        assert native_family_header == gatk_family_header
        assert {row[native_family_header.index("Family")] for row in native_family_rows} == \
               {"all", "fam"}
        native_by_family = {
            row[native_family_header.index("Family")]: row for row in native_family_rows}
        gatk_by_family = {
            row[gatk_family_header.index("Family")]: row for row in gatk_family_rows}
        assert native_by_family.keys() == gatk_by_family.keys()
        float_fields = {"variantRate", "variantRatePerBp", "heterozygosity",
                        "heterozygosityPerBp", "hetHomRatio", "indelRate",
                        "indelRatePerBp", "insertionDeletionRatio"}
        for family in native_by_family:
            for field in native_family_header:
                left = native_by_family[family][native_family_header.index(field)]
                right = gatk_by_family[family][gatk_family_header.index(field)]
                if field in float_fields:
                    assert abs(float(left) - float(right)) < 1e-8, \
                           (family, field, left, right)
                else:
                    assert left == right, (family, field, left, right)
        family_metadata = json.loads(native_family_manifest.read_text(encoding="utf-8"))
        assert family_metadata["compatibility"]["family_stratification"] is True

        # Novelty stratification uses the comparison track to split the same
        # CountVariants evaluator into the stable GATK rows all/known/novel.
        # The Host comparison matcher is deliberately tested with one shared
        # and one novel site so both branches are observable.
        novelty_comp = work / "novelty.comp.vcf"
        novelty_comp.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t1\tknown\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
""",
            encoding="utf-8",
        )
        novelty_report = work / "novelty.gatk.report"
        novelty_manifest = work / "novelty.manifest.json"
        novelty_result = subprocess.run(
            [str(binary), "-eval", str(multi_eval), "-comp", str(novelty_comp),
             "-O", str(novelty_report), "-no-ev", "-EV", "CountVariants",
             "-no-st", "-ST", "Novelty", "--gatk-report",
             "--output-manifest", str(novelty_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert novelty_result.returncode == 0, novelty_result.stderr
        novelty_header, novelty_rows = gatk_table_rows(novelty_report, "CountVariants")
        assert novelty_header[3] == "Novelty"
        novelty_index = novelty_header.index("Novelty")
        novelty_by_name = {row[novelty_index]: row for row in novelty_rows}
        assert set(novelty_by_name) == {"all", "known", "novel"}
        variant_index = novelty_header.index("nVariantLoci")
        assert novelty_by_name["all"][variant_index] == "2"
        assert novelty_by_name["known"][variant_index] == "1"
        assert novelty_by_name["novel"][variant_index] == "1"
        novelty_metadata = json.loads(novelty_manifest.read_text(encoding="utf-8"))
        assert novelty_metadata["compatibility"]["novelty_stratification"] is True
        assert novelty_metadata["telemetry"]["novelty_known_records"] == 1
        assert novelty_metadata["telemetry"]["novelty_novel_records"] == 1

        # Standard GATK IndelSummary and MultiallelicSummary evaluators count
        # sites separately from concrete ALT alleles and distinguish complete
        # versus partial comparison novelty.
        rich_eval = work / "rich.eval.vcf"
        rich_comp = work / "rich.comp.vcf"
        rich_eval.write_text(
            HEADER
            + "chr1\t20\tmulti-snp\tA\tG,C\t50\tPASS\t.\tGT\t1/2\n"
            + "chr1\t21\tmulti-indel\tA\tAT,AAT\t50\tPASS\t.\tGT\t1/2\n",
            encoding="utf-8",
        )
        rich_comp.write_text(
            HEADER
            + "chr1\t20\tmulti-snp\tA\tG\t50\tPASS\t.\tGT\t1/1\n"
            + "chr1\t21\tmulti-indel\tA\tAT\t50\tPASS\t.\tGT\t1/1\n",
            encoding="utf-8",
        )
        rich_report = work / "rich.report"
        rich_manifest = work / "rich.manifest.json"
        rich_result = subprocess.run(
            [str(binary), "-eval", str(rich_eval), "-comp", str(rich_comp),
             "-O", str(rich_report), "-EV", "IndelSummary", "-EV", "MultiallelicSummary",
             "--output-manifest", str(rich_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert rich_result.returncode == 0, rich_result.stderr
        assert metric(rich_report, "IndelSummary", "n_indels") == "2"
        assert metric(rich_report, "IndelSummary", "n_multiallelic_indel_sites") == "1"
        assert metric(rich_report, "IndelSummary", "n_insertions") == "2"
        assert metric(rich_report, "VariantSummary", "nVariantLoci") == "2"
        assert metric(rich_report, "IndelLengthHistogram", "1") == "0.50"
        assert metric(rich_report, "IndelLengthHistogram", "2") == "0.50"
        assert metric(rich_report, "MultiallelicSummary", "nMultiSNPs") == "1"
        assert metric(rich_report, "MultiallelicSummary", "nMultiIndels") == "1"
        assert metric(rich_report, "MultiallelicSummary", "knownSNPsPartial") == "1"
        rich_metadata = json.loads(rich_manifest.read_text(encoding="utf-8"))
        assert rich_metadata["compatibility"]["indel_summary"] is True
        assert rich_metadata["compatibility"]["multiallelic_summary"] is True
        assert rich_metadata["compatibility"]["indel_length_histogram"] is True
        assert rich_metadata["telemetry"]["indel_histogram_total"] == 2

        ft_header = HEADER.replace(
            "#CHROM", "##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>\n#CHROM")
        ft_eval = work / "ft.eval.vcf"
        ft_eval.write_text(
            ft_header
            + "chr1\t30\tft-pass\tA\tG\t50\tPASS\t.\tGT:FT\t0/1:PASS\n"
            + "chr1\t31\tft-nocall\tA\tG\t50\tPASS\t.\tGT:FT\t./.:LowQual\n"
            + "chr1\t32\tft-filtered\tA\tG\t50\tPASS\t.\tGT:FT\t1/1:LowQual\n",
            encoding="utf-8",
        )
        ft_report = work / "ft.report"
        ft_manifest = work / "ft.manifest.json"
        ft_result = subprocess.run(
            [str(binary), "-eval", str(ft_eval), "-O", str(ft_report), "-no-ev",
             "-EV", "GenotypeFilterSummary", "--output-manifest", str(ft_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert ft_result.returncode == 0, ft_result.stderr
        assert metric(ft_report, "GenotypeFilterSummary", "nCalledNotFiltered") == "1"
        assert metric(ft_report, "GenotypeFilterSummary", "nNoCallOrFiltered") == "2"
        ft_metadata = json.loads(ft_manifest.read_text(encoding="utf-8"))
        assert ft_metadata["compatibility"]["genotype_filter_summary"] is True
        assert ft_metadata["telemetry"]["genotype_called_not_filtered"] == 1
        assert ft_metadata["telemetry"]["genotype_no_call_or_filtered"] == 2

        missing_eval = work / "missing.eval.vcf"
        missing_comp = work / "missing.comp.vcf"
        missing_eval.write_text(HEADER + "chr1\t40\tpresent\tA\tG\t50\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
        missing_comp.write_text(
            HEADER
            + "chr1\t40\tpresent\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t41\tmissing\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t42\tfiltered\tC\tT\t50\tLowQual\t.\tGT\t0/1\n"
            + "chr1\t43\tindel\tA\tAT\t50\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        missing_report = work / "missing.report"
        missing_manifest = work / "missing.manifest.json"
        missing_result = subprocess.run(
            [str(binary), "-eval", str(missing_eval), "-comp", str(missing_comp), "-O", str(missing_report),
             "-no-ev", "-EV", "PrintMissingComp", "--output-manifest", str(missing_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert missing_result.returncode == 0, missing_result.stderr
        assert metric(missing_report, "PrintMissingComp", "nMissing") == "1"
        missing_metadata = json.loads(missing_manifest.read_text(encoding="utf-8"))
        assert missing_metadata["compatibility"]["print_missing_comp"] is True
        assert missing_metadata["telemetry"]["missing_comp_snps"] == 1

        listed = subprocess.run([str(binary), "--list"], text=True, capture_output=True, check=False)
        assert listed.returncode == 0 and "CountVariants" in listed.stdout and \
               "VariantSummary" in listed.stdout and "MendelianViolationEvaluator" in listed.stdout

        unsupported_manifest = work / "unsupported.manifest.json"
        unsupported_report = work / "unsupported.report"
        unsupported = subprocess.run(
            [
                str(binary), "-eval", str(eval_vcf), "-O", str(unsupported_report),
                "-EV", "CustomEvaluator", "--output-manifest", str(unsupported_manifest),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        assert unsupported.returncode == 0, unsupported.stderr
        unsupported_metadata = json.loads(unsupported_manifest.read_text(encoding="utf-8"))
        assert unsupported_metadata["fallback"]["unsupported_modules"] is True

    print(json.dumps({"status": "pass", "eval_records": 3, "overlap": 2,
                      "indel_summary": True, "multiallelic_summary": True,
                      "indel_length_histogram": True,
                      "genotype_filter_summary": True,
                      "print_missing_comp": True, "variant_af": True,
                      "theta": True, "mendelian": True, "family": True,
                      "novelty": True,
                      "gatk_report": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
