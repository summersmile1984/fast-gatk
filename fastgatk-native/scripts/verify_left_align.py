#!/usr/bin/env python3
"""Verify reference-backed allele trimming and repeat left-shifting."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=12>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-left-align-trim"
    assert binary.is_file() and os.access(binary, os.X_OK)
    java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-left-align-") as temp:
        work = Path(temp)
        reference = work / "ref.fa"
        source = work / "input.vcf.gz"
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        reference.write_text(">chr1\nAAAAACATGCAT\n", encoding="utf-8")
        # GATK requires the standard FASTA companions for the oracle path.
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\n", encoding="utf-8")
        # One-line FASTA has a deterministic FAI offset/width tuple.
        reference.with_suffix(".fa.fai").write_text("chr1\t12\t6\t12\t13\n", encoding="utf-8")
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(HEADER)
            # Homopolymer deletion shifts from position 4 to the leftmost A.
            handle.write("chr1\t4\t.\tAA\tA\t50\tPASS\tDP=10\tGT\t0/1\n")
            # GATK only normalizes primitive indels; same-length MNPs are
            # emitted unchanged by LeftAlignAndTrimVariants.
            handle.write("chr1\t7\t.\tATG\tACG\t60\tPASS\tDP=12\tGT\t0/1\n")
            handle.write("chr1\t10\t.\tA\tC,G\t70\tPASS\tDP=20;AC=1,1;AN=2;AF=0.5,0.5\tGT:AD:PL:GQ\t0/2:8,2,10:80,50,90,60,70,0:1\n")
        result = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference), "-O", str(output),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["input_records"] == 3 and summary["normalized_records"] == 1
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert [int(record[1]) for record in records] == [1, 7, 10]
        assert records[0][3:5] == ["AA", "A"]
        assert records[1][3:5] == ["ATG", "ACG"]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["common_trim"] is True
        assert metadata["compatibility"]["repeat_left_shift"] is True
        assert all(item["complete"] for item in metadata["outputs"])

        # Plain VCF output uses the shared Tribble LinearIndex v3 writer.  A
        # bundled GATK interval traversal is the authoritative read-back.
        plain_source = work / "input-plain.vcf"
        plain_source.write_text(gzip.open(source, "rt", encoding="utf-8").read(), encoding="utf-8")
        plain_output = work / "output-plain.vcf"
        plain_result = subprocess.run(
            [str(binary), "-V", str(plain_source), "-R", str(reference), "-O", str(plain_output)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.is_file() and Path(f"{plain_output}.idx").is_file()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "output-plain-query.vcf"
        plain_query_result = subprocess.run(
            [str(java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
             "-L", "chr1:1-10", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == 3

        split_output = work / "split.vcf.gz"
        split = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference), "-O", str(split_output),
             "--split-multi-allelics"],
            text=True, capture_output=True, check=False,
        )
        assert split.returncode == 0, split.stderr
        split_text = gzip.open(split_output, "rt", encoding="utf-8").read()
        split_records = [line.split("\t") for line in split_text.splitlines() if line and not line.startswith("#")]
        assert len(split_records) == 4
        split_alts = [record[4] for record in split_records if record[1] == "10"]
        assert split_alts == ["C", "G"]
        assert all(record[9].split(":")[0] in {"0/0", "0/1", "0/.", "./."} for record in split_records if record[1] == "10")
        split_info = {record[4]: record[7] for record in split_records if record[1] == "10"}
        assert "AN=2" in split_info["C"] and "AC=0" in split_info["C"]
        assert "AN=2" in split_info["G"] and "AC=1" in split_info["G"]
        split_gq = {record[4]: record[9].split(":")[-1]
                    for record in split_records if record[1] == "10"}
        assert split_gq["C"] == "30" and split_gq["G"] == "60"

        # A multiallelic site containing a symbolic ALT cannot safely project
        # Number=R/Number=G annotations.  HTSJDK's allele-subsetting path
        # therefore emits one record per ALT with only a no-call GT and clears
        # site annotations; native must follow that fail-closed contract.
        symbolic_header = HEADER.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
            "##ALT=<ID=DEL,Description=Deletion>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        symbolic_source = work / "symbolic-multiallelic.vcf.gz"
        with gzip.open(symbolic_source, "wt", encoding="utf-8") as handle:
            handle.write(symbolic_header)
            handle.write(
                "chr1\t6\t.\tC\tG,<DEL>\t70\tPASS\tDP=20;AC=1,1;AN=2;AF=0.5,0.5\t"
                "GT:AD:PL:GQ\t1/2:8,2,10:80,50,90,60,70,0:1\n")
        symbolic_output = work / "symbolic-multiallelic-native.vcf.gz"
        symbolic_manifest = work / "symbolic-multiallelic.manifest.json"
        symbolic_result = subprocess.run(
            [str(binary), "-V", str(symbolic_source), "-R", str(reference),
             "-O", str(symbolic_output), "--split-multi-allelics",
             "--output-manifest", str(symbolic_manifest)],
            text=True, capture_output=True, check=False)
        assert symbolic_result.returncode == 0, symbolic_result.stderr
        symbolic_records = [line.split("\t") for line in gzip.open(
            symbolic_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[4] for record in symbolic_records] == ["G", "<DEL>"]
        assert all(record[7] == "." and record[8] == "GT" and record[9] == "./."
                   for record in symbolic_records)
        symbolic_metadata = json.loads(symbolic_manifest.read_text(encoding="utf-8"))
        assert symbolic_metadata["telemetry"]["symbolic_split_records"] == 2
        assert symbolic_metadata["telemetry"]["symbolic_info_fields_cleared"] > 0
        assert symbolic_metadata["telemetry"]["symbolic_format_fields_cleared"] > 0
        if oracle_guard.oracle_ready('verify_left_align.py', java, gatk_jar):
            symbolic_plain = work / "symbolic-multiallelic.vcf"
            symbolic_plain.write_text(
                gzip.open(symbolic_source, "rt", encoding="utf-8").read(), encoding="utf-8")
            symbolic_gatk = work / "symbolic-multiallelic-gatk.vcf"
            symbolic_gatk_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "LeftAlignAndTrimVariants",
                "-V", str(symbolic_plain), "-R", str(reference), "-O", str(symbolic_gatk),
                "--split-multi-allelics"], text=True, capture_output=True, check=False)
            assert symbolic_gatk_result.returncode == 0, symbolic_gatk_result.stderr
            symbolic_gatk_records = [line.split("\t") for line in symbolic_gatk.read_text(
                encoding="utf-8").splitlines() if line and not line.startswith("#")]
            assert [record[:10] for record in symbolic_records] == [record[:10]
                                                                    for record in symbolic_gatk_records]

        original_output = work / "split-original-ac.vcf.gz"
        original_result = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference), "-O", str(original_output),
             "--split-multi-allelics", "--keep-original-ac"],
            text=True, capture_output=True, check=False,
        )
        assert original_result.returncode == 0, original_result.stderr
        original_records = [line.split("\t") for line in gzip.open(
            original_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#") and line.split("\t")[1] == "10"]
        assert all("AC_Orig=1" in record[7] and "AF_Orig=0.5" in record[7] and
                   "AN_Orig=2" in record[7] for record in original_records)

        # Repeated intervals are unioned before normalization.  The source
        # MNP at position 7 remains unchanged after selection.
        interval_output = work / "interval.vcf.gz"
        interval_manifest = work / "interval.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference),
             "-L", "chr1:7-7", "--intervals", "chr1:10-10",
             "-O", str(interval_output), "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert interval_result.returncode == 0, interval_result.stderr
        interval_summary = json.loads(interval_result.stdout.splitlines()[-1])
        assert interval_summary["interval_skipped"] == 1
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["7", "10"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["intervals"] == 2
        intersection_output = work / "interval-intersection.vcf.gz"
        intersection_manifest = work / "interval-intersection.manifest.json"
        intersection_result = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference),
             "-L", "chr1:1-10", "-L", "chr1:10-10",
             "--interval-set-rule", "INTERSECTION",
             "-O", str(intersection_output), "--output-manifest", str(intersection_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_summary = json.loads(intersection_result.stdout.splitlines()[-1])
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        intersection_records = [line.split("\t") for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in intersection_records] == ["10"]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        interval_list = work / "left-align.interval_list"
        interval_list.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
                                 "chr1\t7\t7\t+\tfirst\nchr1\t10\t10\t+\tsecond\n",
                                 encoding="utf-8")
        interval_file_output = work / "interval-file.vcf.gz"
        interval_file_manifest = work / "interval-file.manifest.json"
        interval_file_result = subprocess.run(
            [str(binary), "-V", str(source), "-R", str(reference), "-L", str(interval_list),
             "-O", str(interval_file_output), "--output-manifest", str(interval_file_manifest)],
            text=True, capture_output=True, check=False)
        assert interval_file_result.returncode == 0, interval_file_result.stderr
        interval_file_metadata = json.loads(interval_file_manifest.read_text(encoding="utf-8"))
        assert interval_file_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_file_metadata["telemetry"]["interval_list_records"] == 2

        # GATK's realignment bounds are part of the compatibility surface:
        # max-indel-length refuses oversized alleles, while max-leading-bases
        # limits the repeat walk without disabling common trimming.
        long_reference = work / "long-ref.fa"
        long_reference.write_text(">chrL\n" + "A" * 1000 + "\n", encoding="utf-8")
        long_reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chrL\tLN:1000\n", encoding="utf-8")
        long_reference.with_suffix(".fa.fai").write_text(
            "chrL\t1000\t6\t1000\t1001\n", encoding="utf-8")
        long_source = work / "long-indel.vcf"
        long_source.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chrL,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chrL\t500\t.\tAAAA\tA\t50\tPASS\t.\n", encoding="utf-8")

        def first_long_record(path: Path) -> list[str]:
            opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
            with opener as handle:
                return [line.split("\t") for line in handle
                        if line and not line.startswith("#")][0]

        for label, option_args, expected_pos in (
            ("oversized", ["--max-indel-length", "2"], "500"),
            ("inclusive", ["--max-indel-length", "3"], "1"),
            ("lookback2", ["--max-leading-bases", "2"], "498"),
            ("lookback100", ["--max-leading-bases", "100"], "400"),
        ):
            bounded_output = work / f"long-{label}-native.vcf.gz"
            bounded_manifest = work / f"long-{label}.manifest.json"
            bounded_result = subprocess.run(
                [str(binary), "-V", str(long_source), "-R", str(long_reference),
                 "-O", str(bounded_output), "--output-manifest", str(bounded_manifest)] + option_args,
                text=True, capture_output=True, check=False)
            assert bounded_result.returncode == 0, bounded_result.stderr
            bounded_record = first_long_record(bounded_output)
            assert bounded_record[1] == expected_pos
            bounded_metadata = json.loads(bounded_manifest.read_text(encoding="utf-8"))
            if label == "oversized":
                assert bounded_metadata["telemetry"]["oversized_indel_records"] == 1
                assert bounded_metadata["telemetry"]["left_shift_bases"] == 0
            else:
                assert bounded_metadata["telemetry"]["oversized_indel_records"] == 0
            if oracle_guard.oracle_ready('verify_left_align.py', java, gatk_jar):
                bounded_gatk = work / f"long-{label}-gatk.vcf"
                bounded_gatk_result = subprocess.run([
                    str(java), "-Xmx1g", "-jar", str(gatk_jar), "LeftAlignAndTrimVariants",
                    "-V", str(long_source), "-R", str(long_reference), "-O", str(bounded_gatk),
                    *option_args], text=True, capture_output=True, check=False)
                assert bounded_gatk_result.returncode == 0, bounded_gatk_result.stderr
                assert bounded_record[:6] == first_long_record(bounded_gatk)[:6]

        # A gVCF reference block is selected by its END span, not only by its
        # POS.  This keeps LeftAlignAndTrimVariants consistent with the other
        # VCF tools when a workflow passes a locus inside a block.
        gvcf_source = work / "span.g.vcf.gz"
        gvcf_header = HEADER.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate allele>\n")
        with gzip.open(gvcf_source, "wt", encoding="utf-8") as handle:
            handle.write(gvcf_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT\t0/0\n")
        span_output = work / "span.g.vcf.gz.out.vcf.gz"
        span_result = subprocess.run(
            [str(binary), "-V", str(gvcf_source), "-R", str(reference),
             "-L", "chr1:4-4", "-O", str(span_output)],
            text=True, capture_output=True, check=False,
        )
        assert span_result.returncode == 0, span_result.stderr
        span_records = [line.split("\t") for line in gzip.open(
            span_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(span_records) == 1 and span_records[0][1] == "1"

        # Arbitrary-ploidy allele splitting must use the shared Kokkos
        # Number=G remap rather than a diploid triangular shortcut.  For a
        # triploid, three-allele source the target PL rows are [000,001,011,111]
        # and [000,002,022,222] for the two split ALT records.
        triploid_source = work / "triploid.vcf"
        with triploid_source.open("w", encoding="utf-8") as handle:
            handle.write(HEADER)
            handle.write(
                "chr1\t10\ttri\tA\tC,G\t70\tPASS\tDP=20;AC=1,1;AN=3;AF=0.3333333,0.3333333\t"
                "GT:AD:PL:GQ\t0/1/2:10,3,4:90,40,70,30,10,80,20,60,50,0:0\n")
        triploid_output = work / "triploid-split.vcf.gz"
        triploid_result = subprocess.run(
            [str(binary), "-V", str(triploid_source), "-R", str(reference),
             "-O", str(triploid_output), "--split-multi-allelics"],
            text=True, capture_output=True, check=False,
        )
        assert triploid_result.returncode == 0, triploid_result.stderr
        triploid_records = [line.split("\t") for line in gzip.open(
            triploid_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(triploid_records) == 2
        by_alt = {record[4]: record for record in triploid_records}
        assert by_alt["C"][9].split(":") == ["0/1/0", "10,3", "60,10,40,0", "10"]
        assert by_alt["G"][9].split(":") == ["0/0/1", "10,4", "90,10,60,0", "10"]
        triploid_manifest = json.loads(
            Path(f"{triploid_output}.manifest.json").read_text(encoding="utf-8"))
        assert triploid_manifest["telemetry"]["max_ploidy"] == 3
        assert triploid_manifest["telemetry"]["pl_remap_kernel_calls"] == 2
        assert triploid_manifest["telemetry"]["allele_field_remap_kernel_calls"] == 2
        if oracle_guard.oracle_ready('verify_left_align.py', java, gatk_jar):
            gatk_triploid = work / "gatk-triploid-split.vcf"
            gatk_triploid_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "LeftAlignAndTrimVariants",
                "-V", str(triploid_source), "-R", str(reference), "-O", str(gatk_triploid),
                "--split-multi-allelics",
            ], text=True, capture_output=True, check=False)
            assert gatk_triploid_result.returncode == 0, gatk_triploid_result.stderr
            gatk_triploid_records = [line.split("\t") for line in gatk_triploid.read_text(
                encoding="utf-8").splitlines()
                if line and not line.startswith("#")]
            assert [record[0:7] for record in triploid_records] == [record[0:7] for record in gatk_triploid_records]
            def format_values(record: list[str]) -> dict[str, str]:
                names = record[8].split(":")
                values = record[9].split(":")
                return dict(zip(names, values))
            for native_record, gatk_record in zip(triploid_records, gatk_triploid_records):
                native_values = format_values(native_record)
                gatk_values = format_values(gatk_record)
                assert native_values["GT"] == gatk_values["GT"]
                assert native_values["AD"] == gatk_values["AD"]
                assert native_values["PL"] == gatk_values["PL"]
                assert native_values["GQ"] == gatk_values["GQ"]
                native_info = dict(item.split("=", 1) for item in native_record[7].split(";") if "=" in item)
                gatk_info = dict(item.split("=", 1) for item in gatk_record[7].split(";") if "=" in item)
                assert native_info["AC"] == gatk_info["AC"]
                assert native_info["AN"] == gatk_info["AN"]
                assert abs(float(native_info["AF"]) - float(gatk_info["AF"])) < 1e-3

    print(json.dumps({"status": "pass", "normalized_records": summary["normalized_records"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
