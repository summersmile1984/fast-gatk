#!/usr/bin/env python3
"""Compare native ReblockGVCF against bundled GATK on a real HC gVCF."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_records(path: Path) -> list[dict[str, object]]:
    """Normalize INFO/FORMAT ordering while retaining all record values."""
    opener = gzip.open if path.suffix == ".gz" else Path.open
    records: list[dict[str, object]] = []
    sample_names: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, str | bool] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    info[key] = value if separator else True
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            samples = {
                name: dict(zip(format_keys, sample.split(":")))
                for name, sample in zip(sample_names, fields[9:])
            }
            records.append({
                "site": tuple(fields[index] for index in range(7)),
                "info": dict(sorted(info.items())),
                "samples": samples,
            })
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY", str(root / "fastgatk-native/build/fastgatk-reblock-gvcf")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "hc.g.vcf.gz"
        gatk_output = work / "gatk.g.vcf.gz"
        native_output = work / "native.g.vcf.gz"
        hc = subprocess.run([
            str(java), "-jar", str(gatk), "HaplotypeCaller", "-R", str(reference),
            "-I", str(bam), "-L", "17:69000-70000", "-ERC", "GVCF", "-O", str(source),
        ], text=True, capture_output=True, check=False)
        assert hc.returncode == 0, hc.stderr
        gatk_run = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(gatk_output),
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source), "-O", str(native_output),
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        native_records = read_records(native_output)
        gatk_records = read_records(gatk_output)
        assert native_records == gatk_records, (native_records, gatk_records)

        # ReblockGVCF's rgq-threshold is a Java double, not an integer.  A
        # fractional value must be accepted and applied without truncation;
        # compare the complete output against the pinned Java implementation.
        fractional_gatk_output = work / "gatk-rgq-fractional.g.vcf.gz"
        fractional_native_output = work / "native-rgq-fractional.g.vcf.gz"
        fractional_args = ["--rgq-threshold", "10.5"]
        fractional_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(fractional_gatk_output), *fractional_args,
        ], text=True, capture_output=True, check=False)
        assert fractional_gatk.returncode == 0, fractional_gatk.stderr
        fractional_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source),
            "-O", str(fractional_native_output), *fractional_args,
        ], text=True, capture_output=True, check=False)
        assert fractional_native.returncode == 0, fractional_native.stderr
        assert read_records(fractional_gatk_output) == read_records(fractional_native_output)

        # QUAL approximation is an exposed GATK option (including the formal
        # --do-qual-score-approximation spelling).  Compare the Java/native
        # annotation vectors as well; Number=A must retain GATK's leading REF
        # placeholder (for example |43|0).
        qual_gatk_output = work / "gatk-qual.g.vcf.gz"
        qual_native_output = work / "native-qual.g.vcf.gz"
        qual_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(qual_gatk_output),
            "--do-qual-score-approximation",
        ], text=True, capture_output=True, check=False)
        assert qual_gatk.returncode == 0, qual_gatk.stderr
        qual_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source),
            "-O", str(qual_native_output), "--do-qual-score-approximation",
        ], text=True, capture_output=True, check=False)
        assert qual_native.returncode == 0, qual_native.stderr
        assert read_records(qual_gatk_output) == read_records(qual_native_output)

        # --floor-blocks changes the storage shape of reference blocks (GQ is
        # floored and PL/MIN_DP are removed) while preserving variant records.
        # Keep this option in the real Java oracle so both execution spaces
        # exercise the same output contract rather than only the default mode.
        floor_gatk_output = work / "gatk-floor.g.vcf.gz"
        floor_native_output = work / "native-floor.g.vcf.gz"
        floor_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(floor_gatk_output), "--floor-blocks",
        ], text=True, capture_output=True, check=False)
        assert floor_gatk.returncode == 0, floor_gatk.stderr
        floor_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source),
            "-O", str(floor_native_output), "--floor-blocks",
        ], text=True, capture_output=True, check=False)
        assert floor_native.returncode == 0, floor_native.stderr
        assert read_records(floor_gatk_output) == read_records(floor_native_output)

        # With --drop-low-quals, GATK removes existing GQ0 reference blocks
        # before combining adjacent bands.  This must leave holes rather than
        # allowing a retained block to bridge across a dropped span.
        drop_gatk_output = work / "gatk-drop.g.vcf.gz"
        drop_native_output = work / "native-drop.g.vcf.gz"
        drop_native_manifest = work / "native-drop.manifest.json"
        low_qual_source = work / "hc-low-qual.g.vcf"
        with gzip.open(source, "rt", encoding="utf-8") as input_stream, \
                low_qual_source.open("wt", encoding="utf-8") as output_stream:
            for line in input_stream:
                if line.startswith("17\t69067\t"):
                    fields = line.rstrip("\n").split("\t")
                    # A concrete ALT with a hom-ref PL minimum and low
                    # confidence is removed by GATK's drop-mode genotyping
                    # threshold.  Keep the record otherwise well-formed.
                    fields[9] = "0/1:0,1,0:1:0:0,0,0,0,0,0:0,0,0,0"
                    line = "\t".join(fields) + "\n"
                output_stream.write(line)
        index_run = subprocess.run([
            str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(low_qual_source),
        ], text=True, capture_output=True, check=False)
        assert index_run.returncode == 0, index_run.stderr
        drop_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), "-O", str(drop_gatk_output), "--drop-low-quals",
        ], text=True, capture_output=True, check=False)
        assert drop_gatk.returncode == 0, drop_gatk.stderr
        drop_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(source),
            "-O", str(drop_native_output), "--drop-low-quals",
            "--output-manifest", str(drop_native_manifest),
        ], text=True, capture_output=True, check=False)
        assert drop_native.returncode == 0, drop_native.stderr
        drop_gatk_records = read_records(drop_gatk_output)
        drop_native_records = read_records(drop_native_output)
        assert drop_gatk_records == drop_native_records
        assert len(drop_native_records) > len(native_records)
        drop_metadata = json.loads(drop_native_manifest.read_text(encoding="utf-8"))
        assert drop_metadata["telemetry"]["dropped_low_quality_blocks"] > 0

        low_gatk_output = work / "gatk-drop-low-variant.g.vcf.gz"
        low_native_output = work / "native-drop-low-variant.g.vcf.gz"
        low_native_manifest = work / "native-drop-low-variant.manifest.json"
        low_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(low_qual_source), "-O", str(low_gatk_output), "--drop-low-quals",
        ], text=True, capture_output=True, check=False)
        assert low_gatk.returncode == 0, low_gatk.stderr
        low_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(low_qual_source),
            "-O", str(low_native_output), "--drop-low-quals",
            "--output-manifest", str(low_native_manifest),
        ], text=True, capture_output=True, check=False)
        assert low_native.returncode == 0, low_native.stderr
        assert read_records(low_gatk_output) == read_records(low_native_output)
        low_metadata = json.loads(low_native_manifest.read_text(encoding="utf-8"))
        assert low_metadata["telemetry"]["dropped_low_quality_variants"] >= 1

        # A concrete site whose PL-derived call is confidently hom-ref is
        # re-genotyped by GATK and emitted as a reference block, retaining the
        # projected PL/GQ for REF versus <NON_REF> rather than the all-zero
        # low-quality encoding.  This catches the distinct high-confidence
        # conversion branch in the native writer.
        high_ref_source = work / "hc-high-ref.g.vcf"
        with gzip.open(source, "rt", encoding="utf-8") as input_stream, \
                high_ref_source.open("wt", encoding="utf-8") as output_stream:
            for line in input_stream:
                if line.startswith("17\t69067\t"):
                    fields = line.rstrip("\n").split("\t")
                    fields[9] = "0/0:0,0,0:1:50:0,50,100,100,100,100:0,0,0,0"
                    line = "\t".join(fields) + "\n"
                output_stream.write(line)
        high_index_run = subprocess.run([
            str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(high_ref_source),
        ], text=True, capture_output=True, check=False)
        assert high_index_run.returncode == 0, high_index_run.stderr
        high_gatk_output = work / "gatk-drop-high-ref.g.vcf.gz"
        high_native_output = work / "native-drop-high-ref.g.vcf.gz"
        high_native_manifest = work / "native-drop-high-ref.manifest.json"
        high_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(high_ref_source), "-O", str(high_gatk_output), "--drop-low-quals",
        ], text=True, capture_output=True, check=False)
        assert high_gatk.returncode == 0, high_gatk.stderr
        high_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(high_ref_source),
            "-O", str(high_native_output), "--drop-low-quals",
            "--output-manifest", str(high_native_manifest),
        ], text=True, capture_output=True, check=False)
        assert high_native.returncode == 0, high_native.stderr
        assert read_records(high_gatk_output) == read_records(high_native_output)
        high_metadata = json.loads(high_native_manifest.read_text(encoding="utf-8"))
        assert high_metadata["telemetry"]["converted_high_confidence_variants"] >= 1

        # ReblockingGVCFBlockCombiner buffers reference blocks because a
        # concrete call may overlap an input block (most importantly for
        # deletion calls, but this SNP-shaped fixture isolates the writer
        # boundary).  A block spanning 69000-69010 must be split around the
        # concrete call at 69005.  The reference bases are pinned to the
        # bundled human_g1k_v37 chr17 fixture (T at both positions).
        overlap_source = work / "overlap-block.g.vcf"
        overlap_gatk_output = work / "gatk-overlap-block.g.vcf"
        overlap_native_output = work / "native-overlap-block.g.vcf"
        overlap_native_manifest = work / "native-overlap-block.manifest.json"
        overlap_source.write_text(
            """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
17\t69000\t.\tT\t<NON_REF>\t.\t.\tEND=69010\tGT:DP:GQ:MIN_DP:PL\t0/0:10:20:10:0,20,200
17\t69005\t.\tT\tG,<NON_REF>\t50\t.\tDP=10\tGT:DP:AD:PL:GQ\t0/1:10:5,5,0:50,0,80,99,99,99:50
""",
            encoding="utf-8",
        )
        overlap_gatk = subprocess.run([
            str(java), "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(overlap_source), "-O", str(overlap_gatk_output),
            "--add-output-vcf-command-line", "false",
        ], text=True, capture_output=True, check=False)
        assert overlap_gatk.returncode == 0, overlap_gatk.stderr
        overlap_native = subprocess.run([
            str(binary), "-R", str(reference), "-V", str(overlap_source),
            "-O", str(overlap_native_output), "--output-manifest",
            str(overlap_native_manifest),
        ], text=True, capture_output=True, check=False)
        assert overlap_native.returncode == 0, overlap_native.stderr
        overlap_gatk_records = read_records(overlap_gatk_output)
        overlap_native_records = read_records(overlap_native_output)
        assert overlap_native_records == overlap_gatk_records
        assert len(overlap_native_records) == 3
        assert overlap_native_records[0]["site"][1] == "69000"
        assert overlap_native_records[0]["info"]["END"] == "69004"
        assert overlap_native_records[1]["site"][1] == "69005"
        assert overlap_native_records[2]["site"][1] == "69006"
        assert overlap_native_records[2]["info"]["END"] == "69010"
        overlap_metadata = json.loads(overlap_native_manifest.read_text(encoding="utf-8"))
        assert overlap_metadata["compatibility"]["overlapping_ref_block_trim_split"] is True
        assert overlap_metadata["telemetry"]["overlapping_ref_block_split"] == 1

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(native_records),
            "normalized_info_order": True,
            "numeric_fields_exact": True,
            "fractional_rgq_threshold_exact": True,
            "qual_approx_exact": True,
            "floor_blocks_exact": True,
            "drop_low_quals_exact": True,
            "drop_low_quals_records": len(drop_native_records),
            "drop_low_qual_variant_exact": True,
            "high_confidence_homref_exact": True,
            "overlapping_ref_block_split_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
