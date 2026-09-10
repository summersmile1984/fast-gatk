#!/usr/bin/env python3
"""Pinned GATK oracle for combined ModelSegments allele-fraction segmentation.

GATK's MultisampleMultidimensionalKernelSegmenter uses the first allelic site
in each copy-ratio interval and its oriented ALT fraction.  It does not average
all sites in the interval and does not fold ALT fraction to minor-allele
fraction.  The two sites in every interval below have opposite orientation, so
their average ALT fraction is always 0.5 and their MAF is always 0.1; only the
first oriented site carries the chr1:1600 changepoint.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_MODEL_SEGMENTS_BINARY", str(BUILD / "fastgatk-model-segments")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(
            "command failed: " + " ".join(command) +
            "\nstdout:\n" + result.stdout[-2000:] +
            "\nstderr:\n" + result.stderr[-4000:])
    return result


def rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing first-ALT-fraction oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-first-aaf-") as directory:
        work = pathlib.Path(directory)
        copy_ratios = work / "copy-ratios.tsv"
        allelic_counts = work / "allelic-counts.tsv"
        copy_lines = [
            "@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:100000",
            "@RG\tID:GATKCopyNumber\tSM:SAMPLE",
            "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO",
        ]
        allele_lines = [
            "@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:100000",
            "@RG\tID:GATKAllelicCounts\tSM:SAMPLE",
            "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE",
        ]
        for index in range(32):
            start = index * 100 + 1
            copy_lines.append(f"chr1\t{start}\t{start + 99}\t0.0")
            ref, alt = ((90, 10) if index < 16 else (10, 90))
            allele_lines.append(f"chr1\t{start + 10}\t{ref}\t{alt}\tA\tC")
            allele_lines.append(f"chr1\t{start + 20}\t{alt}\t{ref}\tA\tC")
        copy_ratios.write_text("\n".join(copy_lines) + "\n", encoding="ascii")
        allelic_counts.write_text("\n".join(allele_lines) + "\n", encoding="ascii")

        java_dir = work / "java"
        native_dir = work / "native"
        java_dir.mkdir()
        native_dir.mkdir()
        manifest_path = native_dir / "manifest.json"
        common = [
            "--denoised-copy-ratios", str(copy_ratios),
            "--allelic-counts", str(allelic_counts),
            "--output-prefix", "sample",
            "--genotyping-homozygous-log-ratio-threshold", "1000",
            "--maximum-number-of-segments-per-chromosome", "2",
        ]
        run([
            str(JAVA), "-jar", str(JAR), "ModelSegments",
            *common, "-O", str(java_dir),
            "--number-of-samples-copy-ratio", "2",
            "--number-of-burn-in-samples-copy-ratio", "1",
            "--number-of-samples-allele-fraction", "2",
            "--number-of-burn-in-samples-allele-fraction", "1",
            "--maximum-number-of-smoothing-iterations", "0",
            "--number-of-smoothing-iterations-per-fit", "0",
        ])
        run([
            str(NATIVE), *common, "-O", str(native_dir),
            "--output-manifest", str(manifest_path),
        ])

        java_rows = rows(java_dir / "sample.modelFinal.seg")
        native_rows = rows(native_dir / "sample.modelFinal.seg")
        expected_partition = [
            ["chr1", "1", "1600", "16"],
            ["chr1", "1601", "3200", "16"],
        ]
        if [row[:4] for row in java_rows] != expected_partition:
            raise AssertionError({"java": java_rows, "expected": expected_partition})
        if [row[:4] for row in native_rows] != expected_partition:
            raise AssertionError({"native": native_rows, "expected": expected_partition})
        # The bounded native modeller aggregates counts per copy-ratio interval
        # (so its NUM_POINTS_ALLELE_FRACTION remains an explicit non-parity
        # surface), but its folded MAF must stay separate from segmentation.
        for row in native_rows:
            for value in row[8:11]:
                if abs(float(value) - 0.1) > 1.0e-12:
                    raise AssertionError({"native_modeling_maf": native_rows})
        manifest = json.loads(manifest_path.read_text(encoding="ascii"))
        if manifest.get("segmentation_allele_fraction_semantics") != "first-oriented-alt-per-copy-ratio-interval":
            raise AssertionError(manifest)
        if manifest.get("heterozygous_allelic_loci") != 64:
            raise AssertionError(manifest)

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "copy_ratio_intervals": 32,
            "allelic_sites": 64,
            "average_alt_fraction_per_interval": 0.5,
            "minor_allele_fraction_per_site": 0.1,
            "expected_partition": expected_partition,
            "java_native_partition_exact": True,
            "first_oriented_alt_fraction_semantics": True,
            "modeling_maf_kept_separate": True,
            "full_per_locus_modeling_posterior_exact": False,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
