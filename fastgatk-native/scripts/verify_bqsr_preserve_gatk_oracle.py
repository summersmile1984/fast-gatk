#!/usr/bin/env python3
"""Pin ApplyBQSR's preserve-qscores-less-than boundary to GATK 4.6.2.0.

The native report/apply contract already covers ordinary qualities and
quantizers.  This focused oracle changes the input QUAL field to a deterministic
mixture of Q0..Q40 values, applies the same Java-produced recalibration report
with a non-default preserve threshold, and compares Java/native records after
the same pinned HTSJDK reader.  It also checks the semantic promise directly:
qualities below the threshold remain byte-for-byte unchanged.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str], *, capture: bool = False,
        check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=capture)


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    native = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY",
                                build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_bqsr_preserve_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled BQSR preserve oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK BQSR preserve oracle unavailable"}))
        return 0

    threshold = 10
    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-preserve-oracle-") as directory:
        work = Path(directory)
        source = work / "source.sam"
        modified = work / "modified.sam"
        known_sites = work / "known-sites.vcf"
        report = work / "recal.tsv"
        gatk_bam = work / "gatk.bam"
        native_bam = work / "native.bam"
        gatk_sam = work / "gatk.sam"
        native_sam = work / "native.sam"

        run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(bam),
             "--create-output-bam-index", "false", "-O", str(source)])
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])
        # Keep the GATK report independent of the changed input QUAL values.
        # This isolates ApplyBQSR's preservation branch from BaseRecalibrator's
        # training pass while retaining the exact pinned read/header corpus.
        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(bam), "--known-sites", str(known_sites), "-O", str(report)])

        qualities = (0, 2, 5, 9, 10, 11, 20, 30, 40)
        changed_lines: list[str] = []
        for line in source.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("@"):
                changed_lines.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                raise AssertionError(f"malformed source SAM record: {line[:120]}")
            sequence = fields[9]
            fields[10] = "".join(chr(33 + qualities[index % len(qualities)])
                                  for index in range(len(sequence)))
            changed_lines.append("\t".join(fields))
        modified.write_text("\n".join(changed_lines) + "\n", encoding="utf-8")

        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(modified), "--bqsr-recal-file", str(report),
             "--preserve-qscores-less-than", str(threshold),
             "--create-output-bam-index", "false", "-O", str(gatk_bam)])
        native_run = run([str(native), "-I", str(modified), "-R", str(reference),
                          "--bqsr-recal-file", str(report),
                          "--preserve-qscores-less-than", str(threshold),
                          "--create-output-bam-index=false", "-O", str(native_bam)],
                         capture=True)
        native_summary = json.loads(native_run.stdout.strip().splitlines()[-1])

        for input_bam, output_sam in ((gatk_bam, gatk_sam), (native_bam, native_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(input_bam),
                 "--create-output-bam-index", "false", "-O", str(output_sam)])

        java_records = records(gatk_sam)
        native_records = records(native_sam)
        if java_records != native_records:
            raise AssertionError({
                "record_count": (len(java_records), len(native_records)),
                "first_difference": next((index for index, (left, right) in
                                           enumerate(zip(java_records, native_records))
                                           if left != right), None),
            })
        input_records = records(modified)
        if len(input_records) != len(java_records):
            raise AssertionError({"input_records": len(input_records),
                                  "output_records": len(java_records)})

        preserved = 0
        boundary = 0
        for source_record, output_record in zip(input_records, java_records):
            if len(source_record) < 11 or len(output_record) < 11:
                raise AssertionError("record lost QUAL field")
            for raw, recalibrated in zip(source_record[10], output_record[10]):
                quality = ord(raw) - 33
                if quality < threshold:
                    preserved += 1
                    if raw != recalibrated:
                        raise AssertionError({"preserve_threshold": threshold,
                                              "input_quality": quality,
                                              "input_char": raw,
                                              "output_char": recalibrated})
                elif quality == threshold:
                    boundary += 1
        if preserved == 0 or boundary == 0:
            raise AssertionError({"preserved_bases": preserved,
                                  "threshold_boundary_bases": boundary})

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(java_records),
            "preserve_threshold": threshold,
            "preserved_bases": preserved,
            "threshold_boundary_bases": boundary,
            "java_native_records_exact": True,
            "native_records": native_summary["records"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
