#!/usr/bin/env python3
"""Verify streaming VCF shard gather and header/order contracts."""
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
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
"""


def write_shard(path: Path, extra_header: str, body: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(HEADER + extra_header)
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        handle.write(body)


def data_records(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-gather-vcfs"
    assert binary.is_file() and os.access(binary, os.X_OK)
    gatk_java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-vcfs-") as temp:
        work = Path(temp)
        shard1 = work / "shard-1.vcf.gz"
        shard2 = work / "shard-2.vcf.gz"
        reference = work / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
        (work / "reference.fa.fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")
        (work / "reference.dict").write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
        output = work / "gathered.vcf.gz"
        manifest = work / "gathered.manifest.json"
        write_shard(shard1, "##FILTER=<ID=q10,Description=Low quality>\n",
                    "chr1\t1\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n")
        # The second header proves that later FILTER/INFO definitions are
        # merged before the output header is written.
        write_shard(shard2, "##INFO=<ID=ZZ,Number=1,Type=String,Description=Shard tag>\n",
                    "chr1\t2\t.\tC\tT\t60\tPASS\tDP=12;ZZ=two\tGT\t0/0\n")
        result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2), "-O", str(output),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["input_files"] == 2 and summary["input_records"] == 2
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        assert "##FILTER=<ID=q10" in text
        assert "##INFO=<ID=ZZ" in text
        records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert [int(record[1]) for record in records] == [1, 2]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["tool"] == "GatherVcfs"
        assert metadata["implementation"] == "fastgatk-gather-vcfs"
        assert metadata["execution_space"] == "Host"
        assert metadata["determinism"] == "strict"
        assert metadata["primary_output"] == str(output)
        assert metadata["compatibility"]["stream_gather"] is True
        assert metadata["compatibility"]["header_merge"] is True
        assert all(item["complete"] for item in metadata["outputs"])
        assert metadata["telemetry"]["input_files"] == 2
        assert metadata["telemetry"]["input_records"] == 2
        assert metadata["telemetry"]["output_bytes"] == output.stat().st_size
        assert metadata["telemetry"]["index_bytes"] == Path(f"{output}.tbi").stat().st_size
        assert metadata["telemetry"]["wall_seconds"] >= 0.0

        # Plain VCF output uses the shared Tribble LinearIndex v3 writer and
        # remains directly queryable by the pinned GATK reader.
        plain_output = work / "gathered-plain.vcf"
        plain_manifest = work / "gathered-plain.manifest.json"
        plain_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2), "-O", str(plain_output),
             "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        plain_summary = json.loads(plain_result.stdout.splitlines()[-1])
        assert plain_summary["input_records"] == summary["input_records"]
        assert plain_output.exists() and Path(f"{plain_output}.idx").exists()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "gathered-plain-query.vcf"
        plain_query_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
             "-L", "chr1:2-2", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert len(data_records(plain_query)) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True

        # Picard/GATK's metadata and index controls are part of the command
        # line contract.  Keep both the short/long comment spellings and the
        # uppercase CREATE_INDEX spelling accepted by the native binary.
        commented_output = work / "commented.vcf.gz"
        commented_manifest = work / "commented.manifest.json"
        commented_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-CO", "hello", "--COMMENT", "world", "--CREATE_INDEX", "false",
             "--COMPRESSION_LEVEL", "0", "-O", str(commented_output),
             "--output-manifest", str(commented_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert commented_result.returncode == 0, commented_result.stderr
        commented_text = gzip.open(commented_output, "rt", encoding="utf-8").read()
        assert "##GatherVcfs.comment=hello" in commented_text
        assert "##GatherVcfs.comment=world" in commented_text
        assert not Path(f"{commented_output}.tbi").exists()
        commented_metadata = json.loads(commented_manifest.read_text(encoding="utf-8"))
        assert commented_metadata["compatibility"]["vcf_index"] is False
        assert commented_metadata["compatibility"]["comments"] == 2
        assert commented_metadata["compatibility"]["compression_level"] == 0

        # Real Picard/GATK GatherVcfs oracle for the same two-shard boundary.
        # Compare records exactly; header ordering is intentionally validated
        # semantically above because HTSlib and htsjdk canonicalize metadata in
        # different stable orders.
        assert gatk_java.exists() and gatk_jar.exists()

        shard1_plain = work / "shard-1.vcf"
        shard2_plain = work / "shard-2.vcf"
        # Picard/GATK requires every shard header to define every INFO key;
        # the native contract additionally accepts a later-shard definition,
        # which is tested above.  Give the oracle a complete common header.
        shard1_text = gzip.open(shard1, "rt", encoding="utf-8").read()
        shard1_text = shard1_text.replace(
            "#CHROM\tPOS", "##INFO=<ID=ZZ,Number=1,Type=String,Description=Shard tag>\n#CHROM\tPOS", 1)
        shard1_plain.write_text(shard1_text, encoding="utf-8")
        shard2_plain.write_text(gzip.open(shard2, "rt", encoding="utf-8").read(), encoding="utf-8")
        gatk_output = work / "gatk-gathered.vcf"
        gatk_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "GatherVcfs",
             "-I", str(shard1_plain), "-I", str(shard2_plain), "-O", str(gatk_output)],
            text=True, capture_output=True, check=False,
        )
        assert gatk_result.returncode == 0, gatk_result.stderr
        assert data_records(output) == data_records(gatk_output)

        # GatherVcfs exposes the Picard/GATK uppercase long spelling
        # --REFERENCE_SEQUENCE in addition to -R. The native adapter must
        # accept the exact workflow spelling and preserve gathered rows.
        reference_native_output = work / "reference-native.vcf.gz"
        reference_native = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "--REFERENCE_SEQUENCE", str(reference), "-O", str(reference_native_output)],
            text=True, capture_output=True, check=False,
        )
        assert reference_native.returncode == 0, reference_native.stderr
        reference_gatk_output = work / "reference-gatk.vcf"
        reference_gatk = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "GatherVcfs",
             "-I", str(shard1_plain), "-I", str(shard2_plain),
             "--REFERENCE_SEQUENCE", str(reference), "-O", str(reference_gatk_output)],
            text=True, capture_output=True, check=False,
        )
        assert reference_gatk.returncode == 0, reference_gatk.stderr
        assert data_records(reference_native_output) == data_records(reference_gatk_output)

        # GATK scatter callers may provide shard paths in an arbitrary order.
        # REORDER_INPUT_BY_FIRST_VARIANT must sort by the first record using
        # the VCF sequence dictionary, while the default path remains strict.
        reversed_output = work / "reversed-gathered.vcf.gz"
        reversed_result = subprocess.run(
            [str(binary), "-I", str(shard2), "-I", str(shard1),
             "--reorder-input-by-first-variant", "true", "-O", str(reversed_output)],
            text=True, capture_output=True, check=False,
        )
        assert reversed_result.returncode == 0, reversed_result.stderr
        assert [int(record.split("\t")[1]) for record in data_records(reversed_output)] == [1, 2]
        reversed_summary = json.loads(reversed_result.stdout.splitlines()[-1])
        assert reversed_summary["reorder_input_by_first_variant"] is True
        gatk_reversed = work / "gatk-reversed.vcf"
        gatk_reversed_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "GatherVcfs",
             "-I", str(shard2_plain), "-I", str(shard1_plain),
             "-RI", "true", "-O", str(gatk_reversed)],
            text=True, capture_output=True, check=False,
        )
        assert gatk_reversed_result.returncode == 0, gatk_reversed_result.stderr
        assert data_records(reversed_output) == data_records(gatk_reversed)
        strict_reversed = subprocess.run(
            [str(binary), "-I", str(shard2), "-I", str(shard1),
             "-O", str(work / "strict-reversed.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert strict_reversed.returncode != 0
        assert "coordinate-sorted/disjoint" in strict_reversed.stderr

        # HTSlib's native BCF writer emits BCF 2.2, while the pinned
        # GATK/htsjdk reader accepts the BCF 2.1 form produced by GATK.  A
        # native .bcf output must therefore fail closed instead of publishing
        # a file that a downstream GATK stage cannot read; dispatcher
        # --fallback owns that boundary.
        bcf_output = work / "gathered.bcf"
        bcf_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2), "-O", str(bcf_output)],
            text=True, capture_output=True, check=False,
        )
        assert bcf_result.returncode != 0
        assert "BACKEND_UNAVAILABLE" in bcf_result.stderr
        assert not bcf_output.exists()

        # Sequence dictionaries are part of the GatherVcfs compatibility
        # contract, not just a sample/header check.  A conflicting contig
        # length fails closed unless the explicit GATK compatibility override
        # is supplied.
        dictionary_mismatch = work / "dictionary-mismatch.vcf.gz"
        with gzip.open(dictionary_mismatch, "wt", encoding="utf-8") as handle:
            handle.write(HEADER.replace("length=100", "length=101"))
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            handle.write("chr1\t2\t.\tC\tT\t60\tPASS\tDP=12\tGT\t0/0\n")
        dictionary_failed = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(dictionary_mismatch),
             "-O", str(work / "dictionary-failed.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert dictionary_failed.returncode != 0
        assert "sequence dictionary mismatch" in dictionary_failed.stderr
        dictionary_allowed = work / "dictionary-allowed.vcf.gz"
        dictionary_allowed_manifest = work / "dictionary-allowed.manifest.json"
        dictionary_allowed_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(dictionary_mismatch),
             "--disable-sequence-dictionary-validation", "-O", str(dictionary_allowed),
             "--output-manifest", str(dictionary_allowed_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert dictionary_allowed_result.returncode == 0, dictionary_allowed_result.stderr
        dictionary_allowed_metadata = json.loads(dictionary_allowed_manifest.read_text(encoding="utf-8"))
        assert dictionary_allowed_metadata["compatibility"]["sequence_dictionary_validation"] is False

        interval_output = work / "interval-gathered.vcf.gz"
        interval_manifest = work / "interval-gathered.manifest.json"
        interval_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", "chr1:2-2", "--intervals", "chr1:2-2",
             "-O", str(interval_output), "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert interval_result.returncode == 0, interval_result.stderr
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(interval_records) == 1 and interval_records[0][1] == "2"
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["interval_skipped"] == 1

        # GATK -L accepts Picard interval_list files as well as literal
        # intervals.  The native adapter must preserve 1-based inclusive
        # coordinates and ignore SAM-style header lines.
        interval_list = work / "targets.interval_list"
        interval_list.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t2\t2\t+\ttarget\n",
            encoding="utf-8",
        )
        list_output = work / "interval-list-gathered.vcf.gz"
        list_manifest = work / "interval-list-gathered.manifest.json"
        list_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", str(interval_list), "-O", str(list_output),
             "--output-manifest", str(list_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert list_result.returncode == 0, list_result.stderr
        list_records = [line.split("\t") for line in gzip.open(
            list_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(list_records) == 1 and list_records[0][1] == "2"
        list_metadata = json.loads(list_manifest.read_text(encoding="utf-8"))
        assert list_metadata["compatibility"]["interval_list_inputs"] == 1
        assert list_metadata["compatibility"]["interval_list_records"] == 1

        # Repeated selectors support the explicit GATK interval-set rule.
        # The intersection of chr1:1-2 and chr1:2-2 keeps only shard-2's
        # POS=2 record and is surfaced in summary and manifest telemetry.
        intersection_output = work / "intersection-gathered.vcf.gz"
        intersection_manifest = work / "intersection-gathered.manifest.json"
        intersection_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", "chr1:1-2", "-L", "chr1:2-2",
             "--interval-set-rule", "INTERSECTION",
             "-O", str(intersection_output), "--output-manifest", str(intersection_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_records = [line.split("\t") for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(intersection_records) == 1 and intersection_records[0][1] == "2"
        intersection_summary = json.loads(intersection_result.stdout.splitlines()[-1])
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_skipped"] == 1

        opaque_interval_list = work / "targets"
        opaque_interval_list.write_text("chr1\t2\t2\t+\ttarget\n", encoding="utf-8")
        opaque_output = work / "opaque-interval-gathered.vcf.gz"
        opaque_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", str(opaque_interval_list), "-O", str(opaque_output)],
            text=True, capture_output=True, check=False,
        )
        assert opaque_result.returncode == 0, opaque_result.stderr
        opaque_records = [line.split("\t") for line in gzip.open(
            opaque_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(opaque_records) == 1 and opaque_records[0][1] == "2"

        # GATK .intervals/.list files may use one literal per line instead of
        # a Picard interval-list row.  The native parser accepts both forms.
        literal_intervals = work / "targets.intervals"
        literal_intervals.write_text("# targets\nchr1:2-2\n", encoding="utf-8")
        literal_output = work / "literal-interval-gathered.vcf.gz"
        literal_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", str(literal_intervals), "-O", str(literal_output)],
            text=True, capture_output=True, check=False,
        )
        assert literal_result.returncode == 0, literal_result.stderr
        literal_records = [line.split("\t") for line in gzip.open(
            literal_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(literal_records) == 1 and literal_records[0][1] == "2"

        # BED is the one common interval format with zero-based half-open
        # coordinates; keep it explicit by extension instead of guessing from
        # a three-column text file.
        bed = work / "targets.bed.gz"
        with gzip.open(bed, "wt", encoding="utf-8") as handle:
            handle.write("track name=targets\n")
            handle.write("browser position chr1:2-2\n")
            handle.write("chr1\t1\t2\n")
        bed_output = work / "bed-gathered.vcf.gz"
        bed_result = subprocess.run(
            [str(binary), "-I", str(shard1), "-I", str(shard2),
             "-L", str(bed), "-O", str(bed_output)],
            text=True, capture_output=True, check=False,
        )
        assert bed_result.returncode == 0, bed_result.stderr
        bed_records = [line.split("\t") for line in gzip.open(
            bed_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(bed_records) == 1 and bed_records[0][1] == "2"

        # Coordinate validation is interval-aware: a gVCF END span cannot
        # overlap the next shard even when its POS is monotonically increasing.
        block_header = "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
        block1 = work / "block-1.vcf.gz"
        block2 = work / "block-2.vcf.gz"
        write_shard(block1, block_header,
                    "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT\t0/0\n")
        write_shard(block2, block_header,
                    "chr1\t4\t.\tA\tG\t40\tPASS\t.\tGT\t0/1\n")
        overlap = subprocess.run(
            [str(binary), "-I", str(block1), "-I", str(block2),
             "-O", str(work / "overlap.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert overlap.returncode != 0
        assert "overlapping intervals" in overlap.stderr
        allowed = subprocess.run(
            [str(binary), "-I", str(block1), "-I", str(block2),
             "--allow-overlaps", "-O", str(work / "overlap-allowed.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert allowed.returncode == 0, allowed.stderr

        # Interval selection is overlap-based, including gVCF reference
        # blocks.  A block beginning at POS=1 with END=5 must be retained by
        # a locus selector at chr1:4 even though its POS is outside the point.
        block_subset = work / "block-subset.vcf.gz"
        block_subset_manifest = work / "block-subset.manifest.json"
        block_subset_result = subprocess.run(
            [str(binary), "-I", str(block1), "-L", "chr1:4-4",
             "-O", str(block_subset), "--output-manifest", str(block_subset_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert block_subset_result.returncode == 0, block_subset_result.stderr
        block_subset_records = [line for line in gzip.open(
            block_subset, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(block_subset_records) == 1 and block_subset_records[0].split("\t")[1] == "1"
        block_subset_metadata = json.loads(block_subset_manifest.read_text(encoding="utf-8"))
        assert block_subset_metadata["telemetry"]["interval_skipped"] == 0

        # A backward shard is rejected by the default disjoint-coordinate
        # contract; --allow-overlaps is an explicit opt-in.
        bad = work / "backward.vcf.gz"
        write_shard(bad, "", "chr1\t1\t.\tA\tC\t40\tPASS\tDP=5\tGT\t0/1\n")
        rejected = subprocess.run(
            [str(binary), "-I", str(shard2), "-I", str(bad), "-O", str(work / "bad.vcf.gz")],
            text=True, capture_output=True, check=False,
        )
        assert rejected.returncode != 0
        assert "coordinate-sorted/disjoint" in rejected.stderr

    print(json.dumps({"status": "pass", "records": summary["input_records"], "indexed": True,
                      "gatk_oracle_exact": True, "reference_sequence_alias": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
