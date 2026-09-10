#!/usr/bin/env python3
"""GATK/native oracle for multi-sample reference-confidence block splitting."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def body(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\r\n") for line in stream
                if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = root / "third_party/jdk17/bin/java"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.exists() for path in (native, gatk, java, bam, reference)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("multi-sample reference-confidence oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    region = "17:69000-69100"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-multisample-rc-oracle-") as directory:
        work = Path(directory)
        first = work / "sample1.g.vcf.gz"
        second = work / "sample2.g.vcf.gz"
        gatk_cohort = work / "gatk-cohort.g.vcf.gz"
        gatk_output = work / "gatk-all-sites.vcf.gz"
        native_output = work / "native-all-sites.vcf.gz"
        native_stream = work / "native-all-sites-stream.vcf.gz"
        native_manifest = work / "native-all-sites.manifest.json"
        native_stream_manifest = work / "native-all-sites-stream.manifest.json"
        for output in (first, second):
            run([str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                 "-R", str(reference), "-I", str(bam), "-L", region,
                 "-O", str(output), "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
                 "--create-output-variant-index", "true",
                 "--seconds-between-progress-updates", "1"], "GATK HaplotypeCaller")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "RenameSampleInVcf",
             "--INPUT", str(first), "--OUTPUT", str(second),
             "--NEW_SAMPLE_NAME", "SAMPLE2", "--CREATE_INDEX", "true",
             "--VALIDATION_STRINGENCY", "STRICT"], "GATK RenameSampleInVcf")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "CombineGVCFs",
             "-R", str(reference), "-V", str(first), "-V", str(second),
             "-O", str(gatk_cohort), "--create-output-variant-index", "true",
             "--seconds-between-progress-updates", "1"], "GATK CombineGVCFs")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
             "-R", str(reference), "-V", str(gatk_cohort), "-L", region,
             "--include-non-variant-sites", "-O", str(gatk_output),
             "--create-output-variant-index", "false",
             "--seconds-between-progress-updates", "1"], "GATK GenotypeGVCFs")
        for output, manifest, extra in (
                (native_output, native_manifest, []),
                (native_stream, native_stream_manifest, ["--stream-by-locus"])):
            run([str(native), "-R", str(reference), "-V", str(gatk_cohort), "-L", region,
                 "--include-non-variant-sites", "--gatk-compatible-annotations",
                 "--create-output-variant-index=false", "-O", str(output),
                 "--output-manifest", str(manifest), *extra],
                f"native GenotypeGVCFs {'stream' if extra else 'aggregate'}")
        gatk_rows = body(gatk_output)
        native_rows = body(native_output)
        stream_rows = body(native_stream)
        assert gatk_rows and len(gatk_rows) == 101, len(gatk_rows)
        if native_rows != gatk_rows:
            mismatch = next(((index, expected, actual)
                             for index, (expected, actual) in enumerate(zip(gatk_rows, native_rows))
                             if expected != actual), None)
            raise AssertionError({"kind": "aggregate", "gatk_count": len(gatk_rows),
                                  "native_count": len(native_rows), "first_mismatch": mismatch})
        if stream_rows != gatk_rows:
            mismatch = next(((index, expected, actual)
                             for index, (expected, actual) in enumerate(zip(gatk_rows, stream_rows))
                             if expected != actual), None)
            raise AssertionError({"kind": "stream", "gatk_count": len(gatk_rows),
                                  "native_count": len(stream_rows), "first_mismatch": mismatch})
        aggregate_manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        stream_manifest = json.loads(native_stream_manifest.read_text(encoding="utf-8"))
        for manifest in (aggregate_manifest, stream_manifest):
            assert manifest["compatibility"]["include_non_variant_sites"] is True
            assert manifest["compatibility"]["cross_sample_reference_confidence"] is True
            assert manifest["telemetry"]["sample_count"] == 2
        assert stream_manifest["compatibility"]["stream_by_locus"] is True
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(gatk_rows),
            "sample_count": 2,
            "aggregate_rows_exact": True,
            "stream_rows_exact": True,
            "cross_sample_reference_confidence": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
