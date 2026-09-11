#!/usr/bin/env python3
"""Pinned GATK oracle for GenotypeGVCFs -XL/--exclude-intervals.

GATK constructs the include interval set first and then subtracts excluded
 loci.  The native path previously accepted -L but rejected -XL, which made a
 direct replacement fail before traversal.  This fixture contains a concrete
 variant inside the excluded span and compares the Java output with both the
 aggregate and bounded stream-by-locus native paths.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def rows(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-5000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    java = Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genotype_gvcf_exclude_intervals_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GenotypeGVCFs exclusion oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    region = "17:69000-70000"
    excluded = "17:69300-69400"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-exclude-oracle-") as directory:
        work = Path(directory)
        source = work / "input.g.vcf.gz"
        java_output = work / "java.vcf.gz"
        native_aggregate = work / "native.aggregate.vcf.gz"
        native_stream = work / "native.stream.vcf.gz"
        aggregate_manifest = work / "native.aggregate.manifest.json"
        stream_manifest = work / "native.stream.manifest.json"
        run([str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
             "-R", str(reference), "-I", str(bam), "-L", region,
             "-O", str(source), "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
             "--create-output-variant-index", "true",
             "--seconds-between-progress-updates", "1"], "GATK HaplotypeCaller")
        common = ["-R", str(reference), "-V", str(source), "-L", region,
                  "-XL", excluded, "--create-output-variant-index", "false"]
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
             *common, "--seconds-between-progress-updates", "1",
             "-O", str(java_output)], "GATK GenotypeGVCFs -XL")
        run([str(native), *common, "-O", str(native_aggregate),
             "--gatk-compatible-annotations", "--output-manifest", str(aggregate_manifest)],
            "native GenotypeGVCFs -XL aggregate")
        run([str(native), *common, "-O", str(native_stream),
             "--gatk-compatible-annotations", "--stream-by-locus",
             "--output-manifest", str(stream_manifest)],
            "native GenotypeGVCFs -XL stream-by-locus")

        expected = rows(java_output)
        aggregate = rows(native_aggregate)
        stream = rows(native_stream)
        assert expected, "fixture must produce at least one non-excluded variant"
        assert all(not (row[0] == "17" and 69300 <= int(row[1]) <= 69400)
                   for row in expected), expected
        assert expected == aggregate, {"java": expected, "native": aggregate}
        assert expected == stream, {"java": expected, "native_stream": stream}
        for manifest_path in (aggregate_manifest, stream_manifest):
            metadata = json.loads(manifest_path.read_text(encoding="utf-8"))
            assert metadata["compatibility"]["interval_exclusion"] is True
            assert metadata["telemetry"]["exclude_intervals"] == 1

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "include_interval": region,
        "exclude_interval": excluded,
        "java_native_rows_exact": True,
        "aggregate_rows": len(expected),
        "stream_by_locus_rows": len(stream),
        "excluded_variant_removed": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
