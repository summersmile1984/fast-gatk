#!/usr/bin/env python3
"""Pinned GATK oracle for the allele-fraction model initialization.

GATK's ``AlleleFractionInitializer`` integrates the alt-minor/ref-minor
responsibility at f=1/2 and adds one flat-prior pseudocount.  A plain
min(ref, alt)/total average is not equivalent for low-depth asymmetric loci.
This fixture uses (alt, ref)=(2,1) and (1,2), for which the integrated
responsibility is 5/16 and the exact initial segment fraction is 29/64.
The native value is emitted in the manifest while the pinned Java run locks
the same modeling partition and output schema.
"""

from __future__ import annotations

import json
import math
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
        raise AssertionError("command failed: " + " ".join(command) +
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
            raise SystemExit(f"missing allele-fraction initialization oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-af-init-") as directory:
        work = pathlib.Path(directory)
        copy_ratios = work / "copy-ratios.tsv"
        copy_ratios.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
            "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
            "chr1\t1\t10\t0.0\nchr1\t11\t20\t0.0\n",
            encoding="ascii")
        segments = work / "segments.interval_list"
        segments.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "chr1\t1\t20\t+\tseg1\n", encoding="ascii")
        allelic_counts = work / "allelic-counts.tsv"
        allelic_counts.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "@RG\tID:GATKAllelicCounts\tSM:SAMPLE\n"
            "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n"
            "chr1\t5\t1\t2\tA\tC\n"
            "chr1\t15\t2\t1\tA\tC\n", encoding="ascii")

        common = [
            "--segments", str(segments),
            "--allelic-counts", str(allelic_counts),
            "--genotyping-homozygous-log-ratio-threshold", "0",
            "--number-of-samples-allele-fraction", "8",
            "--number-of-burn-in-samples-allele-fraction", "2",
            "--maximum-number-of-smoothing-iterations", "0",
            "--number-of-smoothing-iterations-per-fit", "0",
        ]
        native_dir = work / "native"
        java_dir = work / "java"
        native_dir.mkdir()
        java_dir.mkdir()
        native_manifest = native_dir / "manifest.json"
        run([
            str(NATIVE), "--denoised-copy-ratios", str(copy_ratios),
            "--output-prefix", "sample", "-O", str(native_dir),
            "--output-manifest", str(native_manifest), *common])
        java_common = [
            "--segments", str(segments),
            "--allelic-counts", str(allelic_counts),
            "--genotyping-homozygous-log-ratio-threshold", "0",
            "--number-of-samples-allele-fraction", "8",
            "--number-of-burn-in-samples-allele-fraction", "2",
            "--maximum-number-of-smoothing-iterations", "0",
            "--number-of-smoothing-iterations-per-fit", "0",
        ]
        run([
            str(JAVA), "-jar", str(JAR), "ModelSegments",
            "--denoised-copy-ratios", str(copy_ratios),
            "--output-prefix", "sample", "-O", str(java_dir), *java_common])

        manifest = json.loads(native_manifest.read_text(encoding="ascii"))
        initial = manifest["allele_fraction_initial_segment_means"]
        expected = 29.0 / 64.0
        if len(initial) != 1 or not math.isclose(initial[0], expected, rel_tol=0.0, abs_tol=1.0e-12):
            raise AssertionError({"initial": initial, "expected": expected, "manifest": manifest})
        # The old arithmetic MAF would be 1/3; make the regression detectable
        # even if the exact beta implementation is later refactored.
        if math.isclose(initial[0], 1.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-6):
            raise AssertionError("legacy min(ref,alt)/total initialization was restored")

        native_rows = rows(native_dir / "sample.modelFinal.seg")
        java_rows = rows(java_dir / "sample.modelFinal.seg")
        if len(native_rows) != 1 or len(java_rows) != 1:
            raise AssertionError({"native_rows": native_rows, "java_rows": java_rows})
        if [row[:5] for row in native_rows] != [row[:5] for row in java_rows]:
            raise AssertionError({"native_partition": native_rows, "java_partition": java_rows})
        for row in (native_rows[0], java_rows[0]):
            if len(row) < 11 or not all(math.isfinite(float(row[index])) for index in (7, 8, 9, 10)):
                raise AssertionError({"modeled_row": row})

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "alt_ref_fixture": [[2, 1], [1, 2]],
            "integrated_minor_responsibility": 5.0 / 16.0,
            "expected_initial_minor_fraction": expected,
            "native_initialization_exact": True,
            "java_partition_and_schema_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
