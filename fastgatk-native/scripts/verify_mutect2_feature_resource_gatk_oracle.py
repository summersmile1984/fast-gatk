#!/usr/bin/env python3
"""GATK oracle for Mutect2 germline-resource/PON activity semantics."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import pysam

from verify_mutect2_zero_lod_indel_activity_gatk_oracle import REFERENCE
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def activity_rows(path: Path) -> list[tuple[str, int, int, str, str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("Chromosome\t"):
            continue
        fields = line.split("\t")
        assert len(fields) == 5, fields
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]))
    return rows


def active_at(rows: list[tuple[str, int, int, str, str]], position: int) -> bool:
    for _contig, start, end, label, value in rows:
        if label.startswith("size=") and start <= position < end:
            return float(value) > 0.0
    raise AssertionError((position, rows))


def write_indexed_vcf(path: Path, body: str, *, with_af: bool) -> Path:
    header = "##fileformat=VCFv4.2\n"
    if with_af:
        header += "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    plain = path.with_suffix("")
    plain.write_text(header + body, encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def write_fixture(work: Path) -> tuple[Path, Path, str, str]:
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + REFERENCE + "\n", encoding="utf-8")
    pysam.faidx(str(reference))
    (work / "reference.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:600\n", encoding="utf-8")

    # Insertion is anchored at zero-based reference coordinate 249. A common
    # resource indel at that locus must veto ActivityProfile before graph
    # assembly, whereas --genotype-germline-sites deliberately disables it.
    start = 150
    anchor = 249
    insertion = "G"
    ref_allele = REFERENCE[anchor]
    alt_allele = ref_allele + insertion
    bam = work / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(REFERENCE)}],
        "RG": [{"ID": "rg1", "SM": "sample", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as output:
        for index in range(50):
            read = pysam.AlignedSegment()
            read.query_name = f"ref{index}"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "249M"
            read.query_sequence = REFERENCE[start:start + 249]
            read.query_qualities = pysam.qualitystring_to_array("I" * 249)
            read.set_tag("RG", "rg1")
            output.write(read)
        for index in range(30):
            read = pysam.AlignedSegment()
            read.query_name = f"ins{index}"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "100M1I149M"
            read.query_sequence = (
                REFERENCE[start:anchor + 1] + insertion + REFERENCE[anchor + 1:start + 249]
            )
            read.query_qualities = pysam.qualitystring_to_array("I" * 250)
            read.set_tag("RG", "rg1")
            output.write(read)
    pysam.index(str(bam))
    return reference, bam, ref_allele, alt_allele


def call(command_prefix: list[str], common: list[str], output: Path,
         extra: list[str]) -> tuple[list[str], list[tuple[str, int, int, str, str]]]:
    activity = output.with_suffix(".igv")
    result = run([*command_prefix, *common, *extra, "-O", str(output),
                  "--assembly-region-out", str(activity)])
    assert result.returncode == 0, result.stderr
    return records(output), activity_rows(activity)


def call_without_activity(command_prefix: list[str], common: list[str], output: Path,
                          extra: list[str]) -> list[str]:
    result = run([*command_prefix, *common, *extra, "-O", str(output)])
    assert result.returncode == 0, result.stderr
    return records(output)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    if not all(path.is_file() for path in (java, gatk, native)):
        oracle_guard.oracle_not_verified('verify_mutect2_feature_resource_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 feature-resource oracle is required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-feature-resource-") as directory:
        work = Path(directory)
        reference, bam, ref, alt = write_fixture(work)
        germline = write_indexed_vcf(
            work / "germline.vcf.gz", f"chr1\t250\t.\t{ref}\t{alt}\t.\tPASS\tAF=0.25\n",
            with_af=True)
        pon = write_indexed_vcf(
            work / "pon.vcf.gz", f"chr1\t250\t.\t{ref}\t{alt}\t.\tPASS\t.\n",
            with_af=False)
        common = [
            "-R", str(reference), "-I", str(bam), "--tumor-sample", "sample",
            "-L", "chr1:1-600", "--initial-tumor-lod", "0", "--tumor-lod-to-emit", "0",
            "--min-pruning", "1", "--create-output-variant-index", "false",
        ]
        gatk_prefix = [str(java), "-jar", str(gatk), "Mutect2"]
        native_prefix = [str(native)]

        def paired(label: str, gatk_extra: list[str], native_extra: list[str]) \
                -> tuple[list[str], list[tuple[str, int, int, str, str]]]:
            expected, expected_activity = call(
                gatk_prefix, common, work / f"gatk-{label}.vcf", gatk_extra)
            observed, observed_activity = call(
                native_prefix, common, work / f"native-{label}.vcf", native_extra)
            assert expected == observed, {"label": label, "gatk": expected, "native": observed}
            assert active_at(expected_activity, 249) == active_at(observed_activity, 249), {
                "label": label, "gatk_activity": expected_activity,
                "native_activity": observed_activity,
            }
            return expected, expected_activity

        baseline, baseline_activity = paired("baseline", [], [])
        assert len(baseline) == 1, baseline
        assert active_at(baseline_activity, 249)
        # A common indel halts activity before assembly by default.
        germline_default, germline_default_activity = paired(
            "germline-default", ["--germline-resource", str(germline)],
            ["--germline-resource", str(germline)])
        assert germline_default == []
        assert not active_at(germline_default_activity, 249)
        # This flag bypasses the early germline/normal activity veto and still
        # preserves the AF-derived POPAF writer annotation.
        germline_genotyped, germline_genotyped_activity = paired(
            "germline-genotyped", ["--germline-resource", str(germline),
                                    "--genotype-germline-sites", "true"],
            ["--germline-resource", str(germline), "--genotype-germline-sites=true"])
        assert len(germline_genotyped) == 1 and "POPAF=0.602" in germline_genotyped[0]
        assert active_at(germline_genotyped_activity, 249)
        # PON has the same pre-assembly default veto. With explicit PON-site
        # genotyping it is emitted and marked by SomaticGenotypingEngine.
        pon_default, pon_default_activity = paired(
            "pon-default", ["-pon", str(pon)], ["-pon", str(pon)])
        assert pon_default == []
        assert not active_at(pon_default_activity, 249)
        pon_genotyped, pon_genotyped_activity = paired(
            "pon-genotyped", ["-pon", str(pon), "--genotype-pon-sites", "true"],
            ["-pon", str(pon), "--genotype-pon-sites=true"])
        assert len(pon_genotyped) == 1 and ";PON;" in pon_genotyped[0]
        assert active_at(pon_genotyped_activity, 249)
        # --force-active schedules the zero-activity resource region without
        # changing its writer semantics. In particular, a default PON site is
        # emitted and annotated rather than being silently dropped by a late
        # output-stage proxy.
        pon_forced, _ = paired(
            "pon-force-active", ["-pon", str(pon), "--force-active", "true"],
            ["-pon", str(pon), "--force-active=true"])
        assert len(pon_forced) == 1 and ";PON;" in pon_forced[0]
        # The streaming Host pipeline uses the same resolved masks and passes
        # resource tables through its tile writer. GATK has no corresponding
        # streaming switch, so its aggregate records above are the oracle.
        germline_streamed = call_without_activity(
            native_prefix, common, work / "native-germline-streamed.vcf",
            ["--germline-resource", str(germline), "--stream-by-region", "600"])
        assert germline_streamed == germline_default, {
            "gatk": germline_default, "native_stream": germline_streamed,
        }
        pon_genotyped_streamed = call_without_activity(
            native_prefix, common, work / "native-pon-genotyped-streamed.vcf",
            ["-pon", str(pon), "--genotype-pon-sites=true", "--stream-by-region", "600"])
        assert pon_genotyped_streamed == pon_genotyped, {
            "gatk": pon_genotyped, "native_stream": pon_genotyped_streamed,
        }
        print(json.dumps({
            "status": "pass",
            "baseline_records": len(baseline),
            "germline_default_records": len(germline_default),
            "germline_genotyped_records": len(germline_genotyped),
            "pon_default_records": len(pon_default),
            "pon_genotyped_records": len(pon_genotyped),
            "pon_force_active_records": len(pon_forced),
            "streamed_germline_default_records": len(germline_streamed),
            "streamed_pon_genotyped_records": len(pon_genotyped_streamed),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
