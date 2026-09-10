#!/usr/bin/env python3
"""Verify native FilterMutectCalls TLOD/germline filter contract."""

from __future__ import annotations

import gzip
import io
import json
import math
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=TLOD,Number=1,Type=Float,Description=Tumor log odds>
##INFO=<ID=AF,Number=1,Type=Float,Description=Allele fraction>
##INFO=<ID=F1R2,Number=1,Type=Integer,Description=Forward support>
##INFO=<ID=R1F2,Number=1,Type=Integer,Description=Reverse support>
##INFO=<ID=MBQ,Number=1,Type=Integer,Description=Median alternate base quality>
##INFO=<ID=MMQ,Number=1,Type=Integer,Description=Median alternate mapping quality>
##INFO=<ID=MPOS,Number=1,Type=Integer,Description=Median alternate read position>
##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>
##INFO=<ID=NCount,Number=1,Type=Integer,Description=Number of unknown bases>
##INFO=<ID=ECNT,Number=1,Type=Integer,Description=Potential somatic events in assembly region>
##INFO=<ID=ECNTH,Number=A,Type=Integer,Description=Somatic events in best supporting haplotype>
##FILTER=<ID=germline,Description=Normal support>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
    body = "chr1\t2\t.\tA\tG\t.\tPASS\tTLOD=4.0;AF=0.5;F1R2=1;R1F2=1\n"
    body += "chr1\t3\t.\tA\tC\t.\tPASS\tTLOD=-1.0;AF=0.05;F1R2=1;R1F2=0\n"
    body += "chr1\t4\t.\tA\tT\t.\tgermline\tTLOD=5.0;AF=0.2;F1R2=2;R1F2=0\n"
    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-mutect-") as directory:
        work = Path(directory)
        input_path = work / "mutect.vcf.gz"
        output = work / "filtered.vcf.gz"
        manifest = work / "filtered.vcf.gz.manifest.json"
        with gzip.open(input_path, "wt", encoding="utf-8") as stream:
            stream.write(header + body)
        command = [str(binary), "-V", str(input_path), "-O", str(output),
                   "--min-tlod", "0", "--contamination-fraction", "0.1",
                   "--min-orientation-balance", "0.25", "--output-manifest", str(manifest)]
        env = os.environ.copy()
        result = json.loads(subprocess.check_output(command, text=True, env=env).splitlines()[-1])
        assert result["tool"] == "FilterMutectCalls" and result["status"] == "prototype"
        # FilterMutectCalls clears the provisional input FILTER column before
        # applying its own model.  The third fixture row has no POPAF/NLOD
        # evidence for GermlineFilter, so its inherited ``germline`` value
        # must not be preserved as a final result.
        assert result["input_records"] == 3 and result["pass_records"] == 2, result
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line for line in text.splitlines() if line and not line.startswith("#")]
        assert "\tPASS\tTLOD=4" in records[0], records
        assert "\tFAIL\tTLOD=-1" in records[1], records
        assert "\tPASS\tTLOD=5" in records[2], records
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["tlod_filter"] is True
        assert metadata["compatibility"]["germline_filter"] is True
        assert metadata["compatibility"]["contamination_filter"] is True
        assert metadata["compatibility"]["contamination_posterior_filter"] is False
        assert metadata["compatibility"]["orientation_filter"] is True
        assert metadata["compatibility"]["kokkos_filter_kernel"] is True
        telemetry = metadata["telemetry"]
        assert telemetry["kernel_lifecycle"] == "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] > 0
        assert telemetry["kernel_observations"] >= telemetry["kernel_batches"]
        assert telemetry["kernel_batch_calls"] >= 0
        assert telemetry["kernel_batch_observations"] >= telemetry["kernel_batch_calls"]
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        assert all(item["complete"] for item in metadata["outputs"])

        # FilterMutectCalls follows the same VariantContextWriter boundary:
        # plain VCF gets a Tribble `.idx`, and the GATK reader must be able to
        # query it without a generated Tabix sidecar.
        plain_input = work / "mutect.vcf"
        plain_input.write_text(gzip.open(input_path, "rt", encoding="utf-8").read(), encoding="utf-8")
        plain_output = work / "filtered.vcf"
        plain_manifest = work / "filtered.vcf.manifest.json"
        plain_result = subprocess.run([
            str(binary), "-V", str(plain_input), "-O", str(plain_output),
            "--min-tlod", "0", "--output-manifest", str(plain_manifest),
        ], text=True, capture_output=True, check=False, env=env)
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.is_file() and Path(f"{plain_output}.idx").is_file()
        assert not Path(f"{plain_output}.tbi").exists()
        gatk_java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        plain_query = work / "filtered-query.vcf"
        query_result = subprocess.run([
            str(gatk_java), "-jar", str(gatk_jar), "SelectVariants",
            "-V", str(plain_output), "-L", "chr1:2-2", "-O", str(plain_query),
        ], text=True, capture_output=True, check=False, env=env)
        assert query_result.returncode == 0, query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True
        assert all(item["complete"] for item in plain_metadata["outputs"])

        # GATK's sites-only writer suppresses FORMAT/sample columns only at the
        # output boundary; filtering above still consumed tumor AD evidence.
        sites_only_output = work / "sites-only.filtered.vcf.gz"
        sites_only_manifest = work / "sites-only.filtered.vcf.gz.manifest.json"
        sites_only_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(sites_only_output),
            "--min-tlod", "0", "--sites-only-vcf-output", "true",
            "--output-manifest", str(sites_only_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert sites_only_result["input_records"] == 3
        sites_only_lines = gzip.open(sites_only_output, "rt", encoding="utf-8").read().splitlines()
        sites_only_header = next(line for line in sites_only_lines if line.startswith("#CHROM"))
        assert len(sites_only_header.split("\t")) == 8
        sites_only_records = [line for line in sites_only_lines if line and not line.startswith("#")]
        assert sites_only_records and all(len(line.split("\t")) == 8 for line in sites_only_records)
        sites_only_metadata = json.loads(sites_only_manifest.read_text(encoding="utf-8"))
        assert sites_only_metadata["compatibility"]["sites_only_vcf_output"] is True
        assert sites_only_metadata["telemetry"]["sites_only_vcf_output"] is True
        assert all(item["complete"] for item in sites_only_metadata["outputs"])

        # M2FiltersArgumentCollection's raw field uses -1 as a sentinel, but
        # getMinMedianMappingQuality() resolves it to an effective default of
        # 30 (20 in microbial mode). Verify the effective runtime defaults.
        defaults_input = work / "default-hard-filters.vcf.gz"
        defaults_output = work / "default-hard-filters.out.vcf.gz"
        defaults_manifest = work / "default-hard-filters.out.vcf.gz.manifest.json"
        defaults_body = (
            "chr1\t10\t.\tA\tG\t.\tPASS\tTLOD=4;AF=0.5;F1R2=5;R1F2=5;MBQ=19;MMQ=30;MPOS=5;MFRL=100,100\n"
            "chr1\t11\t.\tA\tC\t.\tPASS\tTLOD=4;AF=0.5;F1R2=5;R1F2=5;MBQ=20;MMQ=29;MPOS=5;MFRL=100,100\n"
            "chr1\t12\t.\tA\tT\t.\tPASS\tTLOD=4;AF=0.5;F1R2=5;R1F2=5;MBQ=20;MMQ=30;MPOS=0;MFRL=100,100\n"
            "chr1\t13\t.\tA\tG\t.\tPASS\tTLOD=4;AF=0.5;F1R2=5;R1F2=5;MBQ=20;MMQ=30;MPOS=5;MFRL=100,10101\n"
            "chr1\t14\t.\tA\tC\t.\tPASS\tTLOD=4;AF=0.5;F1R2=5;R1F2=5;MBQ=20;MMQ=30;MPOS=5;MFRL=100,100\n"
        )
        with gzip.open(defaults_input, "wt", encoding="utf-8") as stream:
            stream.write(header + defaults_body)
        defaults_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(defaults_input), "-O", str(defaults_output),
            "--output-manifest", str(defaults_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert defaults_result["input_records"] == 5
        assert defaults_result["pass_records"] == 1
        defaults_records = [line.split("\t") for line in gzip.open(
            defaults_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "base_qual" in defaults_records[0][6]
        assert "map_qual" in defaults_records[1][6]
        assert "position" in defaults_records[2][6]
        assert "fragment" in defaults_records[3][6]
        assert defaults_records[4][6] == "PASS"
        defaults_metadata = json.loads(defaults_manifest.read_text(encoding="utf-8"))
        assert defaults_metadata["telemetry"]["max_n_ratio"] == "Infinity"

        # GATK's NRatioFilter also skips a present NCount annotation when the
        # summed alternate AD depth is zero. This must remain PASS under the
        # +Infinity default rather than becoming a synthetic n_ratio failure.
        ncount_input = work / "default-n-ratio.vcf.gz"
        ncount_output = work / "default-n-ratio.out.vcf.gz"
        ncount_manifest = work / "default-n-ratio.out.vcf.gz.manifest.json"
        with gzip.open(ncount_input, "wt", encoding="utf-8") as stream:
            stream.write(header + "chr1\t15\t.\tA\tG\t.\tPASS\tTLOD=4;AF=0.5;NCount=100\n")
        ncount_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(ncount_input), "-O", str(ncount_output),
            "--output-manifest", str(ncount_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert ncount_result["input_records"] == 1 and ncount_result["pass_records"] == 1
        ncount_records = [line.split("\t") for line in gzip.open(
            ncount_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert ncount_records[0][6] == "PASS"

        # M2FiltersArgumentCollection mode semantics: mitochondrial and
        # microbial modes remove the genomic-only clustered/multiallelic/
        # fragment/haplotype filters. Microbial mode also lowers the implicit
        # MMQ default from 30 to 20, so MMQ=25 passes there but not in normal
        # or mitochondrial mode.
        mode_body = (
            "chr1\t16\t.\tA\tC,G\t.\tPASS\t"
            "TLOD=4;AF=0.5;MBQ=20;MMQ=25;MPOS=5;MFRL=100,20000;ECNT=4;ECNTH=3\n"
        )
        mode_input = work / "mode-filters.vcf.gz"
        with gzip.open(mode_input, "wt", encoding="utf-8") as stream:
            stream.write(header + mode_body)

        def run_mode(name: str, *flags: str):
            mode_output = work / f"mode-{name}.out.vcf.gz"
            mode_manifest = work / f"mode-{name}.out.vcf.gz.manifest.json"
            payload = json.loads(subprocess.check_output([
                str(binary), "-V", str(mode_input), "-O", str(mode_output),
                "--output-manifest", str(mode_manifest), *flags,
            ], text=True, env=env).splitlines()[-1])
            mode_records = [line.split("\t") for line in gzip.open(
                mode_output, "rt", encoding="utf-8").read().splitlines()
                if line and not line.startswith("#")]
            mode_metadata = json.loads(mode_manifest.read_text(encoding="utf-8"))
            stats_path = next(
                Path(item["path"])
                for item in mode_metadata["outputs"]
                if item["kind"] == "filter-stats"
            )
            mode_stats = json.loads(stats_path.read_text(encoding="utf-8"))
            return payload, mode_records[0], mode_metadata, mode_stats

        normal_mode, normal_record, normal_metadata, normal_stats = run_mode("normal")
        assert normal_mode["pass_records"] == 0
        assert "AS_FilterStatus=map_qual|SITE" in normal_record[7], normal_record
        assert "multiallelic" in normal_record[6]
        assert "fragment" in normal_record[6]
        assert "clustered_events" in normal_record[6]
        assert normal_metadata["compatibility"]["genomic_hard_filters_disabled"] is False
        assert normal_stats["mitochondria_mode"] is False
        assert normal_stats["microbial_mode"] is False

        microbial_mode, microbial_record, microbial_metadata, microbial_stats = run_mode(
            "microbial", "--microbial-mode")
        assert microbial_mode["pass_records"] == 1 and microbial_record[6] == "PASS"
        assert microbial_metadata["compatibility"]["microbial_mode"] is True
        assert microbial_metadata["compatibility"]["genomic_hard_filters_disabled"] is True
        assert microbial_stats["min_median_mapping_quality"] == 20
        assert microbial_stats["log_snv_prior"] == normal_stats["log_snv_prior"]
        assert microbial_stats["log_indel_prior"] == normal_stats["log_indel_prior"]

        mitochondria_mode, mitochondria_record, mitochondria_metadata, mitochondria_stats = run_mode(
            "mitochondria", "--mitochondria-mode")
        assert mitochondria_mode["pass_records"] == 1, (mitochondria_mode, mitochondria_record)
        assert mitochondria_record[6] == "PASS"
        assert "AS_FilterStatus=map_qual|SITE" in mitochondria_record[7], mitochondria_record
        assert "multiallelic" not in mitochondria_record[6]
        assert "fragment" not in mitochondria_record[6]
        assert "clustered_events" not in mitochondria_record[6]
        assert mitochondria_metadata["compatibility"]["mitochondria_mode"] is True
        assert abs(mitochondria_stats["log_snv_prior"] - (-2.5 * math.log(10.0))) < 1.0e-12
        assert abs(mitochondria_stats["log_indel_prior"] - (-3.75 * math.log(10.0))) < 1.0e-12

        # Boolean arguments retain GATK/Picard's optional-value surface.
        # Explicit true and false spellings must be consumed as values rather
        # than leaking into the dispatcher as unknown positional arguments.
        microbial_inline, microbial_inline_record, microbial_inline_metadata, _ = run_mode(
            "microbial-inline", "--microbial-mode=true")
        assert microbial_inline["pass_records"] == 1
        assert microbial_inline_record[6] == "PASS"
        assert microbial_inline_metadata["compatibility"]["microbial_mode"] is True
        microbial_false, microbial_false_record, microbial_false_metadata, _ = run_mode(
            "microbial-false", "--microbial-mode", "false")
        assert microbial_false["pass_records"] == 0
        assert "multiallelic" in microbial_false_record[6]
        assert microbial_false_metadata["compatibility"]["microbial_mode"] is False

        # GATK-style interval selectors are applied to every FilterMutectCalls
        # learning/filtering pass.  INTERSECTION of chr1:1-3 and chr1:2-2
        # keeps only the first fixture record and reports the skipped records.
        interval_output = work / "interval.filtered.vcf.gz"
        interval_manifest = work / "interval.filtered.vcf.gz.manifest.json"
        interval_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(interval_output),
            "-L", "chr1:1-3", "-L", "chr1:2-2",
            "--interval-set-rule", "INTERSECTION",
            "--output-manifest", str(interval_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert interval_result["input_records"] == 1
        assert interval_result["interval_skipped_records"] == 2
        interval_records = [line for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(interval_records) == 1 and "\tPOS=" not in interval_records[0]
        assert "\t2\t" in interval_records[0]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert interval_metadata["telemetry"]["interval_skipped_records"] == 2

        # A disjoint INTERSECTION is an active empty subset, not an omitted
        # selector: all records must be skipped rather than treated as unfiltered.
        empty_output = work / "empty-interval.filtered.vcf.gz"
        empty_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(empty_output),
            "-L", "chr1:1-1", "-L", "chr1:10-10",
            "--interval-set-rule", "INTERSECTION",
        ], text=True, env=env).splitlines()[-1])
        assert empty_result["input_records"] == 0
        assert empty_result["interval_skipped_records"] == 3

        # GATK's --stats is an input Mutect2 stats table.  The native adapter
        # must not overwrite it; it creates a separate filtering-stats JSON
        # unless --filtering-stats is supplied.
        gatk_stats = work / "calls.stats"
        gatk_stats_text = "#<METADATA>threshold=0.234\nstatistic\tvalue\ncallable\t1000000.0\n"
        gatk_stats.write_text(gatk_stats_text, encoding="utf-8")
        stats_input_output = work / "stats-input.filtered.vcf.gz"
        stats_input_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(stats_input_output),
            "--stats", str(gatk_stats),
        ], text=True, env=env).splitlines()[-1])
        assert stats_input_result["input_records"] == 3
        assert gatk_stats.read_text(encoding="utf-8") == gatk_stats_text
        stats_input_json = Path(f"{stats_input_output}.stats.json")
        assert stats_input_json.exists()
        stats_input_manifest = json.loads(Path(f"{stats_input_output}.manifest.json").read_text(encoding="utf-8"))
        assert stats_input_manifest["compatibility"]["input_stats_preserved"] is True

        explicit_filtering_stats = work / "explicit-filtering.stats.json"
        explicit_stats_output = work / "explicit-stats.filtered.vcf.gz"
        explicit_stats_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(explicit_stats_output),
            "--stats", str(gatk_stats), "--filtering-stats", str(explicit_filtering_stats),
        ], text=True, env=env).splitlines()[-1])
        assert explicit_stats_result["input_records"] == 3
        assert explicit_filtering_stats.exists()
        explicit_stats_data = json.loads(explicit_filtering_stats.read_text(encoding="utf-8"))
        assert explicit_stats_data["input_stats"] == str(gatk_stats)
        assert gatk_stats.read_text(encoding="utf-8") == gatk_stats_text
        assert explicit_stats_data["joint_threshold_strategy"] == "OPTIMAL_F_SCORE"
        assert "filtering_output_stats" in explicit_stats_data
        assert explicit_stats_data["filtering_output_stats"]["filters"]["weak_evidence"]["FDR"] >= 0.0
        explicit_stats_manifest = json.loads(Path(f"{explicit_stats_output}.manifest.json").read_text(encoding="utf-8"))
        assert explicit_stats_manifest["telemetry"]["filtering_output_filter_count"] >= 20

        # Direct GATK replacement accepts a table-shaped filtering-stats path;
        # native JSON remains available as a sidecar for telemetry consumers.
        gatk_filtering_table = work / "explicit-filtering.table"
        table_stats_output = work / "explicit-table.filtered.vcf.gz"
        subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(table_stats_output),
            "--stats", str(gatk_stats), "--filtering-stats", str(gatk_filtering_table),
        ], text=True, env=env)
        table_text = gatk_filtering_table.read_text(encoding="utf-8")
        assert "#<METADATA>threshold=" in table_text
        assert "filter\tFP\tFDR\tFN\tFNR" in table_text
        assert Path(f"{gatk_filtering_table}.json").exists()
        table_manifest = json.loads(Path(f"{table_stats_output}.manifest.json").read_text(encoding="utf-8"))
        assert table_manifest["compatibility"]["filtering_stats_format"] == "gatk-filter-stats-table"

        # GATK's contamination resource is accepted directly.  The native
        # adapter records the source and applies the same deterministic AF
        # floor used by the explicit fraction option.
        resource = work / "contamination.table"
        resource.write_text("sample contamination error\nTUMOR whole_bam 0.10 0.001\n", encoding="utf-8")
        table_output = work / "filtered-table.vcf.gz"
        table_manifest = work / "filtered-table.vcf.gz.manifest.json"
        table_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(input_path), "-O", str(table_output),
            "--min-tlod", "0", "--contamination-table", str(resource),
            "--output-manifest", str(table_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert table_result["pass_records"] == 2, table_result
        table_metadata = json.loads(table_manifest.read_text(encoding="utf-8"))
        assert table_metadata["compatibility"]["contamination_resource"] is True
        assert table_metadata["compatibility"]["contamination_source"] == "contamination-table"
        assert table_metadata["compatibility"]["contamination_sample_estimates"] == 1

        # Common Mutect2 hard filters operate on the tumor AD/F1R2/R1F2
        # evidence, not only on TLOD.  Exercise all four filters together and
        # verify deterministic FILTER ordering and sidecar counters.
        hard_header = header.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>\n"
            "##INFO=<ID=NCount,Number=1,Type=Integer,Description=Unknown base count>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n",
        )
        hard_body = (
            "chr1\t10\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;F1R2=3;R1F2=3;MBQ=10;MMQ=30;MPOS=5;ECNT=1;ECNTH=1\tGT:AD\t0/1:8,2\n"
            "chr1\t11\t.\tA\tC,G\t.\tPASS\tTLOD=5;AF=0.4;F1R2=3;R1F2=3;MBQ=30;MMQ=10;MPOS=5;ECNT=1;ECNTH=1,1\tGT:AD\t0/1:8,1,1\n"
            "chr1\t12\t.\tA\tT\t.\tPASS\tTLOD=5;AF=0.1;F1R2=0;R1F2=3;MBQ=30;MMQ=30;MPOS=2;ECNT=4;ECNTH=1\tGT:AD\t0/1:10,1\n"
        )
        hard_input = work / "hard-filters.vcf.gz"
        hard_output = work / "hard-filters.out.vcf.gz"
        hard_stats = work / "hard-filters.stats.json"
        hard_manifest = work / "hard-filters.manifest.json"
        with gzip.open(hard_input, "wt", encoding="utf-8") as stream:
            stream.write(hard_header + hard_body)
        hard_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(hard_input), "-O", str(hard_output),
            "--stats", str(hard_stats), "--min-allele-fraction", "0.2",
            "--min-reads-per-strand", "2", "--unique-alt-read-count", "3",
            "--max-alt-allele-count", "1", "--min-median-base-quality", "20",
            "--min-median-mapping-quality", "20", "--min-median-read-position", "3",
            "--output-manifest", str(hard_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert hard_result["pass_records"] == 0
        hard_text = gzip.open(hard_output, "rt", encoding="utf-8").read()
        hard_records = [line.split("\t") for line in hard_text.splitlines()
                        if line and not line.startswith("#")]
        assert "duplicate" in hard_records[0][6], hard_records
        assert "base_qual" in hard_records[0][6]
        assert "multiallelic" in hard_records[1][6]
        assert "map_qual" in hard_records[1][7], hard_records
        assert "low_allele_frac" in hard_records[2][6]
        assert "position" in hard_records[2][6]
        assert "clustered_events" in hard_records[2][6]
        hard_stats_data = json.loads(hard_stats.read_text(encoding="utf-8"))
        assert hard_stats_data["low_af_records"] == 1
        assert hard_stats_data["strand_records"] == 1
        assert hard_stats_data["low_alt_read_records"] == 3
        assert hard_stats_data["multiallelic_records"] == 1
        assert hard_stats_data["low_median_base_quality_records"] == 1
        assert hard_stats_data["low_median_mapping_quality_records"] == 1
        assert hard_stats_data["low_median_read_position_records"] == 1
        assert hard_stats_data["clustered_events_records"] == 1
        hard_metadata = json.loads(hard_manifest.read_text(encoding="utf-8"))
        assert hard_metadata["compatibility"]["min_allele_fraction_filter"] is True
        assert hard_metadata["compatibility"]["strand_support_filter"] is True
        assert hard_metadata["compatibility"]["unique_alt_read_filter"] is True
        assert hard_metadata["compatibility"]["multiallelic_filter"] is True
        assert hard_metadata["compatibility"]["median_base_quality_filter"] is True
        assert hard_metadata["compatibility"]["median_mapping_quality_filter"] is True
        assert hard_metadata["compatibility"]["median_read_position_filter"] is True
        assert hard_metadata["compatibility"]["clustered_events_filter"] is True

        # Symbolic <NON_REF> alleles may precede concrete ALTs in gVCF-shaped
        # Mutect records.  FilteringOutputStats must key per-ALT probabilities
        # by the original allele index rather than by the compact concrete-ALT
        # offset; otherwise the weak-evidence probability for G is replaced by
        # the zero placeholder for <NON_REF>.
        symbolic_header = hard_header.replace(
            "##INFO=<ID=TLOD,Number=1,Type=Float,Description=Tumor log odds>",
            "##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>",
        ).replace(
            "##INFO=<ID=AF,Number=1,Type=Float,Description=Allele fraction>",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>",
        )
        symbolic_input = work / "symbolic-filtering.vcf.gz"
        symbolic_output = work / "symbolic-filtering.out.vcf.gz"
        symbolic_stats = work / "symbolic-filtering.stats.json"
        with gzip.open(symbolic_input, "wt", encoding="utf-8") as stream:
            stream.write(symbolic_header +
                         "chr1\t15\t.\tA\t<NON_REF>,G\t.\tPASS\t"
                         "TLOD=-10,-10;AF=0,0.5;F1R2=2,2;R1F2=2,2;MBQ=30;MMQ=30;MPOS=5;"
                         "\tGT:AD\t0/1:10,0,5\n")
        subprocess.check_output([
            str(binary), "-V", str(symbolic_input), "-O", str(symbolic_output),
            "--filtering-stats", str(symbolic_stats),
        ], text=True, env=env)
        symbolic_stats_data = json.loads(symbolic_stats.read_text(encoding="utf-8"))
        weak_stats = symbolic_stats_data["filtering_output_stats"]["filters"]["weak_evidence"]
        assert symbolic_stats_data["weak_evidence_records"] == 1
        assert set(("FP", "FDR", "FN", "FNR")) <= set(weak_stats)

        # ErrorProbabilities groups PARTIFACT and orientation bias under the
        # correlated ARTIFACT type.  Verify the native joint reduction takes
        # max(PARTIFACT, OBP) within that type, then applies the independent
        # 1-prod(1-p) reduction across SEQUENCING/NON_SOMATIC/ARTIFACT.
        joint_header = header.replace(
            "##INFO=<ID=TLOD,Number=1,Type=Float,Description=Tumor log odds>",
            "##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>\n"
            "##INFO=<ID=PSOMATIC,Number=A,Type=Float,Description=Somatic posterior>\n"
            "##INFO=<ID=PGERMLINE,Number=A,Type=Float,Description=Germline posterior>\n"
            "##INFO=<ID=PARTIFACT,Number=A,Type=Float,Description=Artifact posterior>\n"
            "##INFO=<ID=OBP,Number=A,Type=Float,Description=Orientation posterior>",
        ).replace(
            "##INFO=<ID=AF,Number=1,Type=Float,Description=Allele fraction>",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>",
        )
        joint_input = work / "joint-error.vcf.gz"
        joint_output = work / "joint-error.out.vcf.gz"
        joint_stats = work / "joint-error.stats.json"
        with gzip.open(joint_input, "wt", encoding="utf-8") as stream:
            stream.write(joint_header +
                         "chr1\t16\t.\tA\tG\t.\tPASS\t"
                         "TLOD=5;AF=0.5;PSOMATIC=0.8;PGERMLINE=0.5;"
                         "PARTIFACT=0.4;OBP=0.3\n")
        subprocess.check_output([
            str(binary), "-V", str(joint_input), "-O", str(joint_output),
            "--filtering-stats", str(joint_stats),
        ], text=True, env=env)
        joint_stats_data = json.loads(joint_stats.read_text(encoding="utf-8"))
        expected_joint_fn = (1.0 - 0.2) * (1.0 - 0.5) * (1.0 - max(0.4, 0.3))
        assert abs(joint_stats_data["filtering_output_stats"]["FN"] - expected_joint_fn) < 2.0e-7

        # A deterministic hard ARTIFACT filter participates in the same
        # ErrorProbabilities max-within-type reduction.  The calibrated
        # posterior alone is below threshold, but low MBQ must raise the
        # joint artifact probability to one.  GATK reports the responsible
        # BaseQualityFilter, not a synthetic combined-error FILTER ID.
        hard_joint_header = joint_header.replace(
            "##INFO=<ID=MBQ,Number=1,Type=Integer,Description=Median alternate base quality>",
            "##INFO=<ID=MBQ,Number=A,Type=Integer,Description=Median alternate base quality>",
        )
        hard_joint_input = work / "hard-joint-error.vcf.gz"
        hard_joint_output = work / "hard-joint-error.out.vcf.gz"
        hard_joint_stats = work / "hard-joint-error.stats.json"
        with gzip.open(hard_joint_input, "wt", encoding="utf-8") as stream:
            stream.write(hard_joint_header +
                         "chr1\t17\t.\tA\tG\t.\tPASS\t"
                         "TLOD=5;AF=0.5;MBQ=10;PSOMATIC=0.99;PGERMLINE=0.0;"
                         "PARTIFACT=0.01;OBP=0.0\n")
        subprocess.check_output([
            str(binary), "-V", str(hard_joint_input), "-O", str(hard_joint_output),
            "--filtering-stats", str(hard_joint_stats), "--threshold-strategy", "CONSTANT",
            "--initial-threshold", "0.5", "--min-median-base-quality", "20",
        ], text=True, env=env)
        hard_joint_records = [line.split("\t") for line in gzip.open(
            hard_joint_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert "base_qual" in hard_joint_records[0][6]
        hard_joint_stats_data = json.loads(hard_joint_stats.read_text(encoding="utf-8"))
        assert hard_joint_stats_data["joint_error_records"] == 1
        assert hard_joint_stats_data["filtering_output_stats"]["FN"] == 0.0

        extended_body = (
            "chr1\t13\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;MFRL=100,140;NCount=3\tGT:AD\t0/1:8,10\n"
            "chr1\t14\t.\tA\tC\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;MFRL=100,105;NCount=0\tGT:AD\t0/1:8,10\n"
        )
        extended_input = work / "extended-hard-filters.vcf.gz"
        extended_output = work / "extended-hard-filters.out.vcf.gz"
        extended_stats = work / "extended-hard-filters.stats.json"
        with gzip.open(extended_input, "wt", encoding="utf-8") as stream:
            stream.write(hard_header + extended_body)
        extended_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(extended_input), "-O", str(extended_output),
            "--stats", str(extended_stats),
            "--max-median-fragment-length-difference", "10", "--max-n-ratio", "0.2",
        ], text=True, env=env).splitlines()[-1])
        assert extended_result["pass_records"] == 1
        extended_records = [line.split("\t") for line in gzip.open(extended_output, "rt", encoding="utf-8").read().splitlines()
                            if line and not line.startswith("#")]
        assert "fragment" in extended_records[0][6] and "n_ratio" in extended_records[0][6]
        assert extended_records[1][6] == "PASS"
        extended_stats_data = json.loads(extended_stats.read_text(encoding="utf-8"))
        assert extended_stats_data["fragment_length_records"] == 1
        assert extended_stats_data["n_ratio_records"] == 1

        # GATK's allele-specific annotations are common in Mutect2 VCFs.  The
        # native adapter must use the ALT pair from AS_SB_TABLE when scalar
        # F1R2/R1F2 are absent, and prefer AS_UNIQ_ALT_READ_COUNT over total AD
        # for duplicate-alt filtering.
        as_header = hard_header.replace(
            "##INFO=<ID=NCount,Number=1,Type=Integer,Description=Unknown base count>\n",
            "##INFO=<ID=NCount,Number=1,Type=Integer,Description=Unknown base count>\n"
            "##INFO=<ID=AS_SB_TABLE,Number=A,Type=String,Description=Allele-specific strand table>\n"
            "##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=String,Description=Unique ALT read count>\n",
        )
        as_input = work / "allele-specific-filters.vcf.gz"
        as_output = work / "allele-specific-filters.out.vcf.gz"
        as_manifest = work / "allele-specific-filters.manifest.json"
        with gzip.open(as_input, "wt", encoding="utf-8") as stream:
            stream.write(as_header)
            stream.write(
                "chr1\t15\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;AS_SB_TABLE=8,1|1,4;AS_UNIQ_ALT_READ_COUNT=2|4\tGT:AD\t0/1:8,10\n"
            )
        as_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(as_input), "-O", str(as_output),
            "--min-reads-per-strand", "2", "--unique-alt-read-count", "3",
            "--output-manifest", str(as_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert as_result["pass_records"] == 0
        as_records = [line.split("\t") for line in gzip.open(as_output, "rt", encoding="utf-8").read().splitlines()
                      if line and not line.startswith("#")]
        assert "duplicate" in as_records[0][7], as_records
        as_metadata = json.loads(as_manifest.read_text(encoding="utf-8"))
        assert as_metadata["telemetry"]["allele_specific_strand_records"] == 1
        assert as_metadata["telemetry"]["allele_specific_unique_records"] == 1

        # Rust/GATK's engine-free hard filters are allele-specific.  Verify
        # natural-log TLOD semantics (4.0 still passes the hard 5.0 ln-LOD
        # threshold), AS_FilterStatus serialization, PON presence, aggregate
        # tumour+normal NRatio, and the MBQ/MMQ/MPOS reference/ALT layouts.
        semantic_header = hard_header.replace(
            "##INFO=<ID=TLOD,Number=1,Type=Float,Description=Tumor log odds>",
            "##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>",
        ).replace(
            "##INFO=<ID=AF,Number=1,Type=Float,Description=Allele fraction>",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>",
        ).replace(
            "##INFO=<ID=MBQ,Number=1,Type=Integer,Description=Median alternate base quality>",
            "##INFO=<ID=MBQ,Number=R,Type=Integer,Description=Median base quality>",
        ).replace(
            "##INFO=<ID=MMQ,Number=1,Type=Integer,Description=Median alternate mapping quality>",
            "##INFO=<ID=MMQ,Number=R,Type=Integer,Description=Median mapping quality>",
        ).replace(
            "##INFO=<ID=MPOS,Number=1,Type=Integer,Description=Median alternate read position>",
            "##INFO=<ID=MPOS,Number=A,Type=Integer,Description=Median read position>",
        ).replace(
            "##INFO=<ID=NCount,Number=1,Type=Integer,Description=Unknown base count>\n",
            "##INFO=<ID=NCount,Number=1,Type=Integer,Description=Unknown base count>\n"
            "##INFO=<ID=PON,Number=0,Type=Flag,Description=Panel of normals>\n"
            "##INFO=<ID=AS_SB_TABLE,Number=A,Type=String,Description=Allele-specific strand table>\n"
            "##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=String,Description=Unique ALT read count>\n",
        ).replace(
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=AF,Number=A,Type=Float,Description=Allele fractions>\n",
        ).replace("FORMAT\tTUMOR\n", "FORMAT\tNORMAL\tTUMOR\n")
        semantic_input = work / "rust-semantic-filters.vcf.gz"
        semantic_output = work / "rust-semantic-filters.out.vcf.gz"
        semantic_stats = work / "rust-semantic-filters.stats.json"
        with gzip.open(semantic_input, "wt", encoding="utf-8") as stream:
            stream.write(semantic_header)
            stream.write(
                "chr1\t40\t.\tA\tC,G\t.\tPASS\t"
                "TLOD=6,4;AF=0.5,0.5;MBQ=30,10,30;MMQ=40,10,30;MPOS=2,5;"
                "AS_SB_TABLE=5,5|0,7|4,4;AS_UNIQ_ALT_READ_COUNT=2|9;"
                "MFRL=100,105;NCount=3;PON;ECNT=1;ECNTH=1,1\t"
                "GT:AD:AF\t0/0:90,0,4:0.5,0.5\t0/1:8,1,1:0.5,0.5\n"
            )
        semantic_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(semantic_input), "-O", str(semantic_output),
            "--stats", str(semantic_stats), "--min-allele-fraction", "0.2",
            "--min-reads-per-strand", "2", "--unique-alt-read-count", "3",
            "--max-alt-allele-count", "1", "--min-median-base-quality", "20",
            "--min-median-mapping-quality", "20", "--min-median-read-position", "3",
            "--max-n-ratio", "0.3", "--max-median-fragment-length-difference", "10",
        ], text=True, env=env).splitlines()[-1])
        assert semantic_result["pass_records"] == 0
        semantic_record = next(line for line in gzip.open(semantic_output, "rt", encoding="utf-8").read().splitlines()
                               if line and not line.startswith("#"))
        semantic_fields = semantic_record.split("\t")
        assert "multiallelic" in semantic_fields[6]
        assert "panel_of_normals" in semantic_fields[6]
        info_fields = dict(item.split("=", 1) for item in semantic_fields[7].split(";") if "=" in item)
        assert "AS_FilterStatus" in info_fields
        as_status = info_fields["AS_FilterStatus"].split("|")
        assert "base_qual" in as_status[0] and "map_qual" in as_status[0]
        assert "position" in as_status[0] and "strict_strand" in as_status[0]
        assert "duplicate" in as_status[0]
        assert as_status[1] == "SITE"
        semantic_stats_data = json.loads(semantic_stats.read_text(encoding="utf-8"))
        assert semantic_stats_data["strict_strand_records"] == 1
        assert semantic_stats_data["duplicate_records"] == 1
        assert semantic_stats_data["panel_of_normals_records"] == 1
        assert semantic_stats_data["as_filter_records"] == 1

        # Number=A evidence must be reduced across concrete ALTs instead of
        # silently taking ALT #1.  Here ALT #1 is weak/low-AF while ALT #2
        # passes; the site remains PASS and only AS_FilterStatus marks ALT #1.
        mixed_header = hard_header.replace(
            "##INFO=<ID=TLOD,Number=1,Type=Float,Description=Tumor log odds>",
            "##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>",
        ).replace(
            "##INFO=<ID=AF,Number=1,Type=Float,Description=Allele fraction>",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>",
        )
        mixed_input = work / "mixed-multialt.vcf.gz"
        mixed_output = work / "mixed-multialt.out.vcf.gz"
        with gzip.open(mixed_input, "wt", encoding="utf-8") as stream:
            stream.write(mixed_header)
            stream.write(
                "chr1\t18\t.\tA\tC,G\t.\tPASS\t"
                "TLOD=-1,5;AF=0.05,0.5;F1R2=2;R1F2=2;MBQ=30;MMQ=30;MPOS=5;ECNT=1;ECNTH=1,1\t"
                "GT:AD\t0/1:99,1,20\n"
            )
        subprocess.check_output([
            str(binary), "-V", str(mixed_input), "-O", str(mixed_output),
            "--min-tlod", "0", "--min-allele-fraction", "0.2",
            "--contamination-fraction", "0.1",
        ], text=True, env=env)
        mixed_record = next(line.split("\t") for line in gzip.open(
            mixed_output, "rt", encoding="utf-8").read().splitlines()
                            if line and not line.startswith("#"))
        assert mixed_record[6] == "PASS", mixed_record
        mixed_info = dict(item.split("=", 1) for item in mixed_record[7].split(";") if "=" in item)
        mixed_status = mixed_info["AS_FilterStatus"].split("|")
        assert "low_tlod" in mixed_status[0] and "low_allele_frac" in mixed_status[0]
        assert mixed_status[1] == "SITE"

        # NALOD plus matched-normal AD exercises the Rust/GATK
        # NormalArtifactFilter boundary: a normal AF below 10% of tumour AF
        # is clean, otherwise the NALOD posterior (or an extreme normal
        # pileup tail) contributes a positive error probability.
        normal_header = semantic_header.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n",
            "##normal_sample=NORMAL\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n",
        ).replace(
            "##INFO=<ID=AS_SB_TABLE,Number=A,Type=String,Description=Allele-specific strand table>\n",
            "##INFO=<ID=NALOD,Number=A,Type=Float,Description=Negative log10 normal artifact odds>\n"
            "##INFO=<ID=AS_SB_TABLE,Number=A,Type=String,Description=Allele-specific strand table>\n",
        )
        normal_input = work / "normal-artifact-filters.vcf.gz"
        normal_output = work / "normal-artifact-filters.out.vcf.gz"
        normal_stats = work / "normal-artifact-filters.stats.json"
        with gzip.open(normal_input, "wt", encoding="utf-8") as stream:
            stream.write(normal_header)
            stream.write(
                "chr1\t41\t.\tA\tC\t.\tPASS\t"
                "TLOD=3;NALOD=2;AF=0.03;MBQ=30,30;MMQ=40,40;MPOS=25;NCount=0\t"
                "GT:AD:AF\t0/0:99,1:0.01\t0/1:97,3:0.03\n"
                "chr1\t42\t.\tA\tC\t.\tPASS\t"
                "TLOD=30;NALOD=2;AF=0.2;MBQ=30,30;MMQ=40,40;MPOS=25;NCount=0\t"
                "GT:AD:AF\t0/0:99,1:0.01\t0/1:80,20:0.2\n"
                "chr1\t43\t.\tA\tC\t.\tPASS\t"
                "TLOD=30;NALOD=-3;AF=0.5;MBQ=30,30;MMQ=40,40;MPOS=25;NCount=0\t"
                "GT:AD:AF\t0/1:60,40:0.4\t0/1:50,50:0.5\n"
            )
        normal_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(normal_input), "-O", str(normal_output),
            "--tumor-sample", "TUMOR", "--stats", str(normal_stats),
        ], text=True, env=env).splitlines()[-1])
        assert normal_result["pass_records"] == 1
        normal_records = [line.split("\t") for line in gzip.open(normal_output, "rt", encoding="utf-8").read().splitlines()
                          if line and not line.startswith("#")]
        assert "normal_artifact" in normal_records[0][6]
        assert normal_records[1][6] == "PASS"
        assert "normal_artifact" in normal_records[2][6]
        normal_stats_data = json.loads(normal_stats.read_text(encoding="utf-8"))
        assert normal_stats_data["normal_artifact_records"] == 2

        # RPA/RU and an indel exercise the deterministic PolymeraseSlippage
        # boundary: one repeat contraction at the default minimum length is
        # evaluated by the Kokkos regularized-beta/initial-cluster kernel and
        # emits the STRQ annotation plus the site filter.
        slippage_header = normal_header.replace(
            "##INFO=<ID=NALOD,Number=A,Type=Float,Description=Negative log10 normal artifact odds>\n",
            "##INFO=<ID=NALOD,Number=A,Type=Float,Description=Negative log10 normal artifact odds>\n"
            "##INFO=<ID=RPA,Number=R,Type=Integer,Description=Repeats per allele>\n"
            "##INFO=<ID=RU,Number=1,Type=String,Description=Repeat unit>\n",
        )
        slippage_input = work / "slippage.vcf.gz"
        slippage_output = work / "slippage.out.vcf.gz"
        slippage_stats = work / "slippage.stats.json"
        with gzip.open(slippage_input, "wt", encoding="utf-8") as stream:
            stream.write(slippage_header)
            stream.write(
                "chr1\t44\t.\tAA\tA\t.\tPASS\t"
                "TLOD=0.5;NALOD=2;RPA=10,9;RU=A;AF=0.03;MBQ=30,30;MMQ=40,40;MPOS=25;NCount=0\t"
                "GT:AD:AF\t0/0:99,1:0.01\t0/1:97,3:0.03\n"
            )
        slippage_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(slippage_input), "-O", str(slippage_output),
            "--tumor-sample", "TUMOR", "--stats", str(slippage_stats),
        ], text=True, env=env).splitlines()[-1])
        assert slippage_result["pass_records"] == 0
        slippage_record = next(line.split("\t") for line in gzip.open(slippage_output, "rt", encoding="utf-8").read().splitlines()
                               if line and not line.startswith("#"))
        assert "slippage" in slippage_record[6]
        slippage_info = dict(item.split("=", 1) for item in slippage_record[7].split(";") if "=" in item)
        assert "STRQ" in slippage_info
        slippage_stats_data = json.loads(slippage_stats.read_text(encoding="utf-8"))
        assert "weak_evidence" in slippage_record[6]
        assert slippage_stats_data["slippage_records"] == 1

        # Sample order is not a contract in GATK VCFs.  Explicit tumor-sample
        # selection must make AD-based filters use the requested sample rather
        # than silently taking the first column.
        reordered_header = hard_header.replace("FORMAT\tTUMOR\n", "FORMAT\tNORMAL\tTUMOR\n")
        reordered_input = work / "reordered-samples.vcf.gz"
        reordered_output = work / "reordered-samples.out.vcf.gz"
        reordered_manifest = work / "reordered-samples.manifest.json"
        with gzip.open(reordered_input, "wt", encoding="utf-8") as stream:
            stream.write(reordered_header)
            stream.write(
                "chr1\t30\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;F1R2=3;R1F2=3\tGT:AD\t0/1:8,20\t0/1:8,1\n"
            )
        reordered_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(reordered_input), "-O", str(reordered_output),
            "--tumor-sample", "TUMOR", "--unique-alt-read-count", "3",
            "--output-manifest", str(reordered_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert reordered_result["pass_records"] == 0
        reordered_text = gzip.open(reordered_output, "rt", encoding="utf-8").read()
        assert "duplicate" in next(line for line in reordered_text.splitlines() if not line.startswith("#"))
        reordered_metadata = json.loads(reordered_manifest.read_text(encoding="utf-8"))
        assert reordered_metadata["compatibility"]["tumor_sample_selection"] is True

        # Mutect2's native posterior annotations can be consumed directly by
        # deterministic probability filters. Missing posterior values are
        # fail-closed when the corresponding threshold is enabled.
        posterior_header = header.replace(
            "##FILTER=<ID=germline,Description=Normal support>\n",
            "##FILTER=<ID=germline,Description=Normal support>\n"
            "##INFO=<ID=PSOMATIC,Number=1,Type=Float,Description=Somatic posterior>\n"
            "##INFO=<ID=PGERMLINE,Number=1,Type=Float,Description=Germline posterior>\n"
            "##INFO=<ID=PARTIFACT,Number=1,Type=Float,Description=Artifact posterior>\n",
        )
        posterior_body = (
            "chr1\t20\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;PSOMATIC=0.05;PGERMLINE=0.1;PARTIFACT=0.1\n"
            "chr1\t21\t.\tA\tC\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;PSOMATIC=0.9;PGERMLINE=0.8;PARTIFACT=0.9\n"
        )
        posterior_input = work / "posterior.vcf.gz"
        posterior_output = work / "posterior.out.vcf.gz"
        posterior_stats = work / "posterior.stats.json"
        posterior_manifest = work / "posterior.manifest.json"
        with gzip.open(posterior_input, "wt", encoding="utf-8") as stream:
            stream.write(posterior_header + posterior_body)
        posterior_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(posterior_input), "-O", str(posterior_output),
            "--stats", str(posterior_stats), "--min-somatic-probability", "0.5",
            "--max-germline-probability", "0.5", "--max-artifact-probability", "0.5",
            "--output-manifest", str(posterior_manifest),
        ], text=True, env=env).splitlines()[-1])
        # These compatibility-only thresholds do not replace GATK's
        # per-filter ErrorProbabilities decision boundary.  Raw INFO
        # posterior annotations alone therefore remain PASS here.
        assert posterior_result["pass_records"] == 2, posterior_result
        posterior_text = gzip.open(posterior_output, "rt", encoding="utf-8").read()
        posterior_records = [line.split("\t") for line in posterior_text.splitlines()
                             if line and not line.startswith("#")]
        assert all(record[6] == "PASS" for record in posterior_records)
        posterior_stats_data = json.loads(posterior_stats.read_text(encoding="utf-8"))
        assert posterior_stats_data["low_somatic_probability_records"] == 1
        assert posterior_stats_data["high_germline_probability_records"] == 1
        assert posterior_stats_data["high_artifact_probability_records"] == 1
        posterior_metadata = json.loads(posterior_manifest.read_text(encoding="utf-8"))
        assert posterior_metadata["compatibility"]["somatic_probability_filter"] is True
        assert posterior_metadata["compatibility"]["germline_probability_filter"] is True
        assert posterior_metadata["compatibility"]["artifact_probability_filter"] is True

        # When the caller asks for a threshold strategy, reconstruct the
        # ErrorProbabilities product over the exposed posterior types in a
        # streaming pre-pass.  The second record's combined error is above
        # 0.96 while the first is below it.  The joint value is retained for
        # threshold telemetry, while GATK still reports only concrete filter
        # names in the VCF FILTER column.
        joint_output = work / "joint-error.out.vcf.gz"
        joint_stats = work / "joint-error.stats.json"
        joint_manifest = work / "joint-error.manifest.json"
        joint_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(posterior_input), "-O", str(joint_output),
            "--stats", str(joint_stats), "--threshold-strategy", "CONSTANT",
            "--initial-threshold", "0.96", "--output-manifest", str(joint_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert joint_result["input_records"] == 2
        joint_text = gzip.open(joint_output, "rt", encoding="utf-8").read()
        joint_records = [line.split("\t") for line in joint_text.splitlines()
                         if line and not line.startswith("#")]
        assert joint_records[0][6] == "PASS"
        assert joint_records[1][6] == "PASS"
        joint_stats_data = json.loads(joint_stats.read_text(encoding="utf-8"))
        assert joint_stats_data["joint_error_records"] == 1
        assert joint_stats_data["joint_error_threshold_learned"] is True
        assert joint_stats_data["joint_error_threshold_observations"] == 0, joint_stats_data
        joint_metadata = json.loads(joint_manifest.read_text(encoding="utf-8"))
        assert joint_metadata["compatibility"]["joint_error_probability_filter"] is True

        model_header = header.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n",
        )
        model_body = (
            "chr1\t22\t.\tA\tG\t.\tPASS\tTLOD=0.5;AF=0.2;F1R2=2;R1F2=2\tGT:AD\t0/1:8,2\n"
            "chr1\t23\t.\tA\tC\t.\tPASS\tTLOD=5;AF=0.7;F1R2=2;R1F2=2\tGT:AD\t0/1:8,20\n"
        )
        model_input = work / "empirical-model.vcf.gz"
        model_output = work / "empirical-model.out.vcf.gz"
        model_input_stats = work / "empirical-model.input.stats"
        model_stats = work / "empirical-model.stats.json"
        with gzip.open(model_input, "wt", encoding="utf-8") as stream:
            stream.write(model_header + model_body)
        model_input_stats.write_text("statistic\tvalue\ncallable\t1000000\n", encoding="utf-8")
        subprocess.check_output([
            str(binary), "-V", str(model_input), "-O", str(model_output),
            "--stats", str(model_input_stats), "--filtering-stats", str(model_stats),
            "--threshold-strategy", "CONSTANT",
            "--initial-threshold", "0.1", "--log-snv-prior", "-13.815510557964274",
            "--log-indel-prior", "-16.11809565095832", "--log-artifact-prior", "-2.302585092994046",
            "--normal-p-value-threshold", "0.001", "--min-pcr-slippage-size", "8",
            "--pcr-slippage-rate", "0.1",
        ], text=True, env=env)
        model_stats_data = json.loads(model_stats.read_text(encoding="utf-8"))
        assert model_stats_data["empirical_somatic_model_learned"] is True
        assert model_stats_data["empirical_somatic_model_records"] == 0, model_stats_data
        assert model_stats_data["callable_sites"] == 1000000
        assert len(model_stats_data["empirical_variant_priors"]) == 21
        assert abs(model_stats_data["log_snv_prior"] + 13.815510557964274) < 1.0e-9
        assert abs(model_stats_data["log_indel_prior"] + 16.11809565095832) < 1.0e-9
        assert model_stats_data["empirical_cluster_count"] >= 2
        assert len(model_stats_data["empirical_cluster_weights"]) == model_stats_data["empirical_cluster_count"]
        assert len(model_stats_data["empirical_cluster_means"]) == model_stats_data["empirical_cluster_count"]

        # A multimodal AF spectrum exercises GATK's probability-weighted peak
        # initialization and BIC-gated dynamic cluster splitting.
        cluster_values = [5, 10, 15, 20, 25] + [60] * 90
        cluster_body = "".join(
            f"chr1\t{index}\t.\tA\tG\t.\tPASS\tTLOD=8;AF={af / 100:.2f};F1R2=20;R1F2=20\tGT:AD\t0/1:{100 - af},{af}\n"
            for index, af in enumerate(cluster_values, start=1)
        )
        cluster_input = work / "empirical-clusters.vcf.gz"
        cluster_output = work / "empirical-clusters.out.vcf.gz"
        cluster_stats = work / "empirical-clusters.stats.json"
        with gzip.open(cluster_input, "wt", encoding="utf-8") as stream:
            stream.write(model_header + cluster_body)
        subprocess.check_output([
            str(binary), "-V", str(cluster_input), "-O", str(cluster_output),
            "--stats", str(model_input_stats), "--filtering-stats", str(cluster_stats),
            "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1",
        ], text=True, env=env)
        cluster_stats_data = json.loads(cluster_stats.read_text(encoding="utf-8"))
        assert cluster_stats_data["empirical_dynamic_peak_clustering"] is True
        assert cluster_stats_data["empirical_cluster_count"] >= 2
        assert len(cluster_stats_data["empirical_cluster_means"]) >= 2

        # FilteredHaplotypeFilter is a stateful GATK filter: it learns the
        # combined artifact probability for each PGT+PID haplotype in a
        # streaming first pass, then answers from the learned loci with an
        # inclusive distance window.  The distant zero-probability record is
        # retained, while both nearby records see the learned 0.8 maximum.
        haplotype_header = posterior_header.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=AF,Number=A,Type=Float,Description=Genotype allele fraction>\n"
            "##FORMAT=<ID=PGT,Number=1,Type=String,Description=Phased genotype>\n"
            "##FORMAT=<ID=PID,Number=1,Type=String,Description=Phasing identifier>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n",
        )
        haplotype_body = (
            "chr1\t40\t.\tA\tG\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;PARTIFACT=0.8\tGT:AD:AF:PGT:PID\t0/1:5,5:0.5:0|1:100_A_C\n"
            "chr1\t80\t.\tA\tC\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;PARTIFACT=0.2\tGT:AD:AF:PGT:PID\t0/1:5,5:0.5:0|1:100_A_C\n"
            "chr1\t200\t.\tA\tT\t.\tPASS\tTLOD=5;AF=0.5;F1R2=2;R1F2=2;PARTIFACT=0.0\tGT:AD:AF:PGT:PID\t0/1:5,5:0.5:0|1:100_A_C\n"
        )
        haplotype_input = work / "haplotype-input.vcf.gz"
        haplotype_output = work / "haplotype-output.vcf.gz"
        haplotype_stats = work / "haplotype.stats.json"
        haplotype_manifest = work / "haplotype.manifest.json"
        with gzip.open(haplotype_input, "wt", encoding="utf-8") as stream:
            stream.write(haplotype_header + haplotype_body)
        haplotype_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(haplotype_input), "-O", str(haplotype_output),
            "--stats", str(haplotype_stats), "--max-intra-haplotype-distance", "100",
            "--output-manifest", str(haplotype_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert haplotype_result["input_records"] == 3
        haplotype_text = gzip.open(haplotype_output, "rt", encoding="utf-8").read()
        haplotype_records = [line.split("\t") for line in haplotype_text.splitlines()
                             if line and not line.startswith("#")]
        assert "haplotype" in haplotype_records[0][6]
        assert "haplotype" in haplotype_records[1][6]
        assert haplotype_records[2][6] == "PASS"
        haplotype_stats_data = json.loads(haplotype_stats.read_text(encoding="utf-8"))
        assert haplotype_stats_data["haplotype_records"] == 2
        assert haplotype_stats_data["haplotype_learning_records"] == 3
        assert haplotype_stats_data["learned_haplotypes"] == 1
        haplotype_metadata = json.loads(haplotype_manifest.read_text(encoding="utf-8"))
        assert haplotype_metadata["compatibility"]["haplotype_filter"] is True
        assert haplotype_metadata["compatibility"]["haplotype_learning"] is True

        # Standard GATK LearnReadOrientationModel output is consumed directly:
        # use a real three-mer reference context and FORMAT F1R2/F2R1 arrays,
        # then exercise the Kokkos 12-state beta-binomial posterior path.
        reference = work / "reference.fa"
        reference.write_text(">chr1\nAACGTT\n", encoding="utf-8")
        (work / "reference.fa.fai").write_text("chr1\t6\t6\t6\t7\n", encoding="utf-8")
        columns = [
            "context", "rev_comp", "f1r2_a", "f1r2_c", "f1r2_g", "f1r2_t",
            "f2r1_a", "f2r1_c", "f2r1_g", "f2r1_t", "hom_ref", "germline_het",
            "somatic_het", "hom_var", "num_examples", "num_alt_examples",
        ]
        complements = str.maketrans("ACGT", "TGCA")
        prior_rows = []
        for left in "ACGT":
            for middle in "ACGT":
                for right in "ACGT":
                    context = left + middle + right
                    masses = [0.0] * 12
                    masses[1] = 0.90  # F1R2_C artifact
                    masses[5] = 0.01  # F2R1_C artifact
                    masses[8:] = [0.01, 0.01, 0.03, 0.01]
                    normalizer = sum(masses)
                    values = [value / normalizer for value in masses]
                    prior_rows.append("\t".join([
                        context, context.translate(complements)[::-1],
                        *[f"{value:.12g}" for value in values], "10", "10",
                    ]))
        prior_text = "#<METADATA>SAMPLE=TUMOR\n" + "\t".join(columns) + "\n" + "\n".join(prior_rows) + "\n"
        prior_archive = work / "TUMOR.orientation_priors.tar.gz"
        with tarfile.open(prior_archive, "w:gz") as archive:
            member = tarfile.TarInfo("TUMOR.orientation_priors")
            member.size = len(prior_text.encode("utf-8"))
            archive.addfile(member, io.BytesIO(prior_text.encode("utf-8")))
        orientation_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=6>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=F1R2,Number=A,Type=Integer,Description=Forward orientation alternate support>
##FORMAT=<ID=F2R1,Number=A,Type=Integer,Description=Reverse orientation alternate support>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
"""
        orientation_body = (
            "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:10:0\n"
            "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:5:5\n"
            "chr1\t2\t.\tAC\tCG\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:10:0\n"
        )
        orientation_input = work / "orientation-input.vcf.gz"
        with gzip.open(orientation_input, "wt", encoding="utf-8") as stream:
            stream.write(orientation_header + orientation_body)
        orientation_output = work / "orientation-output.vcf.gz"
        orientation_stats = work / "orientation.stats.json"
        orientation_manifest = work / "orientation.manifest.json"
        orientation_result = json.loads(subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(orientation_input),
            "-O", str(orientation_output), "--stats", str(orientation_stats),
            "--orientation-bias-artifact-priors", str(prior_archive),
            "--max-orientation-artifact-probability", "0.5",
            "--output-manifest", str(orientation_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert orientation_result["input_records"] == 3
        assert orientation_result["pass_records"] == 1
        orientation_text = gzip.open(orientation_output, "rt", encoding="utf-8").read()
        orientation_records = [line.split("\t") for line in orientation_text.splitlines()
                               if line and not line.startswith("#")]
        assert "orientation" in orientation_records[0][6]
        assert orientation_records[1][6] == "PASS"
        assert "ROQ=" in orientation_records[0][7] and "ROQ=" in orientation_records[1][7]
        orientation_stats_data = json.loads(orientation_stats.read_text(encoding="utf-8"))
        assert orientation_stats_data["orientation_prior_records"] == 3
        assert orientation_stats_data["orientation_records"] == 2
        orientation_metadata = json.loads(orientation_manifest.read_text(encoding="utf-8"))
        assert orientation_metadata["compatibility"]["orientation_prior_table"] is True
        assert orientation_metadata["compatibility"]["orientation_prior_filter"] is True

        # Standard Mutect2 emits F1R2/F2R1 as Number=R (REF, then each ALT).
        # The orientation posterior must use those vectors for both the ALT
        # count and its total depth; FORMAT/AD is only the weighted-median
        # sample weight.  Keep the orientation vectors identical while making
        # AD deliberately disagree and require identical ROQ values.
        standard_orientation_header = orientation_header.replace(
            "##FORMAT=<ID=F1R2,Number=A,", "##FORMAT=<ID=F1R2,Number=R,").replace(
            "##FORMAT=<ID=F2R1,Number=A,", "##FORMAT=<ID=F2R1,Number=R,")
        standard_orientation_input = work / "orientation-standard-input.vcf.gz"
        with gzip.open(standard_orientation_input, "wt", encoding="utf-8") as stream:
            stream.write(standard_orientation_header)
            stream.write(
                "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1"
                "\t0/1:1,10:90,10:10,0\n"
                "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1"
                "\t0/1:1000,10:90,10:10,0\n"
            )
        standard_orientation_output = work / "orientation-standard-output.vcf.gz"
        standard_orientation_result = json.loads(subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(standard_orientation_input),
            "-O", str(standard_orientation_output),
            "--orientation-bias-artifact-priors", str(prior_archive),
            "--max-orientation-artifact-probability", "0.5",
        ], text=True, env=env).splitlines()[-1])
        assert standard_orientation_result["input_records"] == 2
        standard_orientation_records = [line.split("\t") for line in gzip.open(
            standard_orientation_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        standard_orientation_roq = [
            next((item[4:] for item in record[7].split(";") if item.startswith("ROQ=")), None)
            for record in standard_orientation_records
        ]
        assert standard_orientation_roq[0] is not None and standard_orientation_roq[0] == standard_orientation_roq[1], \
            "Number=R orientation posterior must not depend on FORMAT/AD depth"

        # AD and orientation observations are separate GATK evidence vectors.
        # A retained orientation observation may legitimately make the two
        # depths differ; the native posterior must consume F1R2/F2R1 as its
        # own beta-binomial trials and expose the discrepancy in telemetry.
        disagreement_input = work / "orientation-disagreement-input.vcf.gz"
        disagreement_output = work / "orientation-disagreement-output.vcf.gz"
        disagreement_manifest = work / "orientation-disagreement.manifest.json"
        with gzip.open(disagreement_input, "wt", encoding="utf-8") as stream:
            stream.write(orientation_header)
            stream.write(
                "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,9:10:1\n"
            )
        disagreement_result = json.loads(subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(disagreement_input),
            "-O", str(disagreement_output),
            "--orientation-bias-artifact-priors", str(prior_archive),
            "--max-orientation-artifact-probability", "0.5",
            "--output-manifest", str(disagreement_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert disagreement_result["input_records"] == 1
        disagreement_metadata = json.loads(disagreement_manifest.read_text(encoding="utf-8"))
        assert disagreement_metadata["telemetry"]["orientation_depth_disagreements"] == 1

        # GATK ReadOrientationFilter aggregates all non-normal genotypes with
        # a weighted median, where the weight is the total alternate AD depth.
        # The first record would be filtered by either tumor alone; the second
        # must pass because the high-depth second tumor has balanced evidence.
        tumor1_archive = work / "TUMOR1.orientation_priors.tar.gz"
        tumor2_archive = work / "TUMOR2.orientation_priors.tar.gz"
        for archive_path, sample_name in ((tumor1_archive, "TUMOR1"), (tumor2_archive, "TUMOR2")):
            sample_text = prior_text.replace("SAMPLE=TUMOR", f"SAMPLE={sample_name}", 1)
            with tarfile.open(archive_path, "w:gz") as archive:
                member = tarfile.TarInfo(f"{sample_name}.orientation_priors")
                encoded = sample_text.encode("utf-8")
                member.size = len(encoded)
                archive.addfile(member, io.BytesIO(encoded))
        multi_header = orientation_header.replace(
            "##contig=<ID=chr1,length=6>\n",
            "##contig=<ID=chr1,length=6>\n##normal_sample=NORMAL\n",
        ).replace("FORMAT\tTUMOR\n", "FORMAT\tTUMOR1\tTUMOR2\tNORMAL\n")
        multi_body = (
            "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1"
            "\t0/1:20,100:100:0\t0/1:20,1:1:1\t0/0:100,0:0:0\n"
            "chr1\t2\t.\tA\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1"
            "\t0/1:20,1:1:0\t0/1:20,100:50:50\t0/0:100,0:0:0\n"
        )
        multi_input = work / "orientation-multi-input.vcf.gz"
        with gzip.open(multi_input, "wt", encoding="utf-8") as stream:
            stream.write(multi_header + multi_body)
        multi_output = work / "orientation-multi-output.vcf.gz"
        multi_stats = work / "orientation-multi.stats.json"
        multi_manifest = work / "orientation-multi.manifest.json"
        multi_result = json.loads(subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(multi_input),
            "-O", str(multi_output), "--stats", str(multi_stats),
            "--orientation-bias-artifact-priors", str(tumor1_archive),
            "--orientation-bias-artifact-priors", str(tumor2_archive),
            "--max-orientation-artifact-probability", "0.5",
            "--output-manifest", str(multi_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert multi_result["input_records"] == 2 and multi_result["pass_records"] == 1
        multi_text = gzip.open(multi_output, "rt", encoding="utf-8").read()
        multi_records = [line.split("\t") for line in multi_text.splitlines()
                         if line and not line.startswith("#")]
        assert "orientation" in multi_records[0][6]
        assert multi_records[1][6] == "PASS"
        multi_stats_data = json.loads(multi_stats.read_text(encoding="utf-8"))
        assert multi_stats_data["orientation_prior_records"] == 2
        assert multi_stats_data["orientation_prior_sample_evaluations"] == 4
        multi_metadata = json.loads(multi_manifest.read_text(encoding="utf-8"))
        assert multi_metadata["compatibility"]["orientation_weighted_median"] is True

        # The learned threshold path performs a streaming pre-pass and then
        # reopens the VCF for final writing.  Exercise GATK's OPTIMAL_F_SCORE
        # strategy without retaining the records in native memory.
        learned_output = work / "orientation-learned-output.vcf.gz"
        learned_stats = work / "orientation-learned.stats.json"
        learned_manifest = work / "orientation-learned.manifest.json"
        learned_result = json.loads(subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(multi_input),
            "-O", str(learned_output), "--stats", str(learned_stats),
            "--orientation-bias-artifact-priors", str(tumor1_archive),
            "--orientation-bias-artifact-priors", str(tumor2_archive),
            "--threshold-strategy", "OPTIMAL_F_SCORE", "--initial-threshold", "0.1",
            "--f-score-beta", "1.0", "--output-manifest", str(learned_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert learned_result["input_records"] == 2
        learned_stats_data = json.loads(learned_stats.read_text(encoding="utf-8"))
        assert learned_stats_data["orientation_threshold_strategy"] == "OPTIMAL_F_SCORE"
        assert learned_stats_data["orientation_threshold_learned"] is True
        assert learned_stats_data["orientation_threshold_observations"] == 2
        assert 0.0 <= learned_stats_data["max_orientation_artifact_probability"] <= 1.0
        learned_metadata = json.loads(learned_manifest.read_text(encoding="utf-8"))
        assert learned_metadata["compatibility"]["orientation_threshold_learning"] is True

        fdr_output = work / "orientation-fdr-output.vcf.gz"
        fdr_stats = work / "orientation-fdr.stats.json"
        subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(multi_input),
            "-O", str(fdr_output), "--stats", str(fdr_stats),
            "--orientation-bias-artifact-priors", str(tumor1_archive),
            "--orientation-bias-artifact-priors", str(tumor2_archive),
            "--threshold-strategy", "FALSE_DISCOVERY_RATE",
            "--false-discovery-rate", "0.05",
        ], text=True, env=env)
        fdr_stats_data = json.loads(fdr_stats.read_text(encoding="utf-8"))
        assert fdr_stats_data["orientation_threshold_strategy"] == "FALSE_DISCOVERY_RATE"
        assert fdr_stats_data["orientation_threshold_learned"] is True
        assert fdr_stats_data["orientation_threshold_observations"] == 2

        constant_output = work / "orientation-constant-output.vcf.gz"
        constant_stats = work / "orientation-constant.stats.json"
        subprocess.check_output([
            str(binary), "-R", str(reference), "-V", str(multi_input),
            "-O", str(constant_output), "--stats", str(constant_stats),
            "--orientation-bias-artifact-priors", str(tumor1_archive),
            "--orientation-bias-artifact-priors", str(tumor2_archive),
            "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.25",
        ], text=True, env=env)
        constant_stats_data = json.loads(constant_stats.read_text(encoding="utf-8"))
        assert constant_stats_data["orientation_threshold_strategy"] == "CONSTANT"
        assert constant_stats_data["orientation_threshold_learned"] is False
        assert abs(constant_stats_data["max_orientation_artifact_probability"] - 0.25) < 1.0e-6
        print(json.dumps({"status": "pass", "records": len(records), "indexed": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
