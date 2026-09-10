#!/usr/bin/env python3
"""Verify repeatable BAM/CRAM interval selectors across the main read tools."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def run_json(command: list[str]) -> dict:
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    completed = subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = root / "fastgatk-native/build"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    with tempfile.TemporaryDirectory(prefix="fastgatk-bam-intervals-") as directory:
        work = Path(directory)
        interval_list = work / "regions.interval_list"
        interval_list.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:17\tLN:1000000\n"
            "17\t69000\t69020\t+\tleft\n"
            "17\t69620\t69640\t+\tright\n",
            encoding="utf-8",
        )
        bed = work / "regions.bed.gz"
        with gzip.open(bed, "wt", encoding="utf-8") as handle:
            # UCSC BED headers are valid input and must not be interpreted as
            # malformed interval records by the shared Host parser.
            handle.write("track name=fastgatk\n")
            handle.write("browser position 17:69000-69020\n")
            handle.write("17\t68999\t69020\n17\t69619\t69640\n")

        hc_output = work / "hc.vcf.gz"
        hc_manifest = work / "hc.manifest.json"
        run_json([
            str(build / "fastgatk-hc-call"), "-I", str(bam), "-R", str(reference),
            "-L", str(interval_list), "-L", "17:69900-69910", "-O", str(hc_output),
            "--min-depth", "1", "--min-alt-support", "1", "--batch-records", "16",
            "--output-manifest", str(hc_manifest),
        ])
        hc_meta = json.loads(hc_manifest.read_text(encoding="utf-8"))
        hc_telemetry = hc_meta["telemetry"]
        assert hc_telemetry["interval_list_inputs"] == 1
        assert hc_telemetry["interval_list_records"] == 2
        assert hc_telemetry["interval_subset"] is True
        assert hc_telemetry["reads"] < 493
        for line in gzip.open(hc_output, "rt", encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            assert fields[0] == "17"
            assert 69000 <= int(fields[1]) <= 69910
            assert (69000 <= int(fields[1]) <= 69020 or
                    69620 <= int(fields[1]) <= 69640 or
                    69900 <= int(fields[1]) <= 69910)

        mutect_output = work / "mutect.vcf.gz"
        mutect_manifest = work / "mutect.manifest.json"
        run_json([
            str(build / "fastgatk-mutect2"), "-I", str(bam), "-R", str(reference),
            "-L", str(interval_list), "-L", "17:69900-69910", "-O", str(mutect_output),
            "--normal-input", str(bam), "--min-depth", "1",
            "--min-alt-support", "1", "--max-haplotype-paths", "8",
            "--max-haplotype-depth", "64", "--output-manifest", str(mutect_manifest),
        ])
        mutect_meta = json.loads(mutect_manifest.read_text(encoding="utf-8"))
        mutect_telemetry = mutect_meta["telemetry"]
        assert mutect_telemetry["interval_list_inputs"] == 1
        assert mutect_telemetry["interval_list_records"] == 2
        assert mutect_telemetry["tumor_reads"] < 493
        assert Path(f"{mutect_output}.tbi").exists()

        bqsr_report = work / "bqsr.tsv"
        bqsr_manifest = work / "bqsr.manifest.json"
        bqsr_result = run_json([
            str(build / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", str(bed), "-L", "17:69900-69910", "--batch-records", "3",
            "-O", str(bqsr_report), "--output-manifest", str(bqsr_manifest),
        ])
        bqsr_meta = json.loads(bqsr_manifest.read_text(encoding="utf-8"))
        bqsr_telemetry = bqsr_meta["telemetry"]
        assert bqsr_telemetry["interval_file_inputs"] == 1
        assert bqsr_telemetry["interval_file_records"] == 2
        assert bqsr_result["records"] < 493

        print(json.dumps({
            "status": "pass",
            "interval_list_records": 2,
            "hc_reads": hc_telemetry["reads"],
            "mutect_reads": mutect_telemetry["tumor_reads"],
            "bqsr_records": bqsr_result["records"],
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
