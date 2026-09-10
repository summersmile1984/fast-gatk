#!/usr/bin/env python3
"""Contract checks for the native IndexFeatureFile HTSlib adapter."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("FASTGATK_INDEX_FEATURE_FILE_BINARY", ROOT / "fastgatk-native/build/fastgatk-index-feature-file"))
BGZF_VCF = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/4featuresHG38Header.unindexed.vcf.gz"
BGZF_BCF = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/test_variants_for_index.BCF22compressed.bcf.blockgz.gz"
BGZF_BED = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/test_bed_for_index.bed.gz"
PLAIN_GVCF = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/test_variants_for_index.g.vcf"
PLAIN_VCF = ROOT / "gatk-source/src/test/resources/oneSNP.vcf"
PLAIN_BED = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/test_bed_for_index.bed"
EMPTY_VCF = ROOT / "gatk-source/src/test/resources/empty.vcf"
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK_JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def main() -> int:
    if not BINARY.is_file():
        print(json.dumps({"status": "skip", "suite": "index-feature-file", "reason": f"missing binary: {BINARY}"}))
        return 0
    assert BGZF_VCF.is_file(), BGZF_VCF
    with tempfile.TemporaryDirectory(prefix="fastgatk-index-feature-file-") as directory:
        work = Path(directory)
        index = work / "sample.g.vcf.gz.tbi"
        manifest = work / "sample.manifest.json"
        result = subprocess.run(
            [str(BINARY), "-I", str(BGZF_VCF), "-O", str(index), "--threads", "2",
             "--output-manifest", str(manifest)], text=True, capture_output=True,
        )
        assert result.returncode == 0, (result.stdout, result.stderr)
        assert index.is_file() and index.stat().st_size > 0
        payload = json.loads(manifest.read_text())
        assert payload["tool"] == "IndexFeatureFile"
        assert payload["status"] == "contract-compatible"
        assert payload["format"] == "VCF"
        assert payload["compression"] == "BGZF"
        assert payload["index_type"] == "TABIX"
        assert payload["telemetry"]["kernel_execution_policy"] == "RangePolicy"

        default_input = work / "default.vcf.gz"
        default_input.write_bytes(BGZF_VCF.read_bytes())
        default_result = subprocess.run(
            [str(BINARY), "-I", str(default_input), "--threads", "1"],
            text=True, capture_output=True,
        )
        assert default_result.returncode == 0, (default_result.stdout, default_result.stderr)
        assert Path(str(default_input) + ".tbi").is_file()

        # Java IndexFeatureFile chooses a linear Tribble index for an
        # uncompressed path ending in .g.vcf.  Compare the complete binary
        # output, then make a real GATK interval query through the native .idx
        # to prove that HTSJDK can consume it.
        plain_input = work / "sample.g.vcf"
        plain_input.write_bytes(PLAIN_GVCF.read_bytes())
        plain_manifest = work / "plain.manifest.json"
        plain_result = subprocess.run(
            [str(BINARY), "-I", str(plain_input), "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True,
        )
        assert plain_result.returncode == 0, (plain_result.stdout, plain_result.stderr)
        plain_index = Path(str(plain_input) + ".idx")
        assert plain_index.is_file() and plain_index.stat().st_size > 0
        plain_payload = json.loads(plain_manifest.read_text())
        assert plain_payload["compression"] == "NONE"
        assert plain_payload["index_type"] == "TRIBBLE_LINEAR_GVCF"
        assert plain_payload["linear_records"] == 1
        assert plain_payload["linear_contigs"] == 1
        assert plain_payload["linear_blocks"] == 1

        if JAVA.is_file() and GATK_JAR.is_file():
            java_index = work / "java.idx"
            java_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "IndexFeatureFile",
                 "-I", str(plain_input), "-O", str(java_index)],
                text=True, capture_output=True,
            )
            assert java_result.returncode == 0, java_result.stderr
            assert plain_index.read_bytes() == java_index.read_bytes()
            query_output = work / "plain-query.vcf"
            query_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "SelectVariants",
                 "-V", str(plain_input), "-L", "1:150-160", "-O", str(query_output),
                 "--create-output-variant-index", "false"],
                text=True, capture_output=True,
            )
            assert query_result.returncode == 0, query_result.stderr
            assert any(line.startswith("1\t100\t") for line in query_output.read_text().splitlines())

        # Ordinary uncompressed VCFs now follow HTSJDK DynamicIndexCreator's
        # FOR_SEEK_TIME choice: sparse files use an adaptive LinearIndex, while
        # dense files use the 75-record interval-tree candidate.  The sparse
        # fixture is kept as a query contract because HTSJDK's adaptive merge
        # retains an implementation-specific empty leading block.
        plain_vcf_input = work / "sample.vcf"
        plain_vcf_input.write_bytes(PLAIN_VCF.read_bytes())
        plain_vcf_manifest = work / "plain-vcf.manifest.json"
        plain_vcf_result = subprocess.run(
            [str(BINARY), "-I", str(plain_vcf_input), "--output-manifest", str(plain_vcf_manifest)],
            text=True, capture_output=True,
        )
        assert plain_vcf_result.returncode == 0, (plain_vcf_result.stdout, plain_vcf_result.stderr)
        plain_vcf_index = Path(str(plain_vcf_input) + ".idx")
        assert plain_vcf_index.is_file() and plain_vcf_index.stat().st_size > 0
        plain_vcf_payload = json.loads(plain_vcf_manifest.read_text())
        assert plain_vcf_payload["compression"] == "NONE"
        assert plain_vcf_payload["index_type"] == "TRIBBLE_LINEAR"
        assert plain_vcf_payload["linear_records"] == 1
        if JAVA.is_file() and GATK_JAR.is_file():
            java_vcf_index = work / "java-plain-vcf.idx"
            java_vcf_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "IndexFeatureFile",
                 "-I", str(plain_vcf_input), "-O", str(java_vcf_index)],
                text=True, capture_output=True,
            )
            assert java_vcf_result.returncode == 0, java_vcf_result.stderr
            assert plain_vcf_index.read_bytes() == java_vcf_index.read_bytes()
            query_output = work / "plain-vcf-query.vcf"
            query_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "SelectVariants",
                 "-V", str(plain_vcf_input), "-L", "1:110201699-110201699",
                 "-O", str(query_output), "--create-output-variant-index", "false"],
                text=True, capture_output=True,
            )
            assert query_result.returncode == 0, query_result.stderr
            assert any(line.startswith("1\t110201699\t") for line in query_output.read_text().splitlines())

        # Force the interval-tree candidate with a dense, low-coordinate VCF.
        # This exercises the same 75-record grouping and proves that HTSJDK
        # can consume the native type-2 index, rather than only accepting the
        # linear sparse path.
        dense_input = work / "dense.vcf"
        source_lines = PLAIN_VCF.read_text().splitlines()
        header = [line for line in source_lines if line.startswith("#")]
        record = next(line for line in source_lines if not line.startswith("#")).split("\t")
        dense_records = []
        for position in range(1, 301):
            fields = list(record)
            fields[0] = "1"
            fields[1] = str(position)
            dense_records.append("\t".join(fields))
        dense_input.write_text("\n".join(header + dense_records) + "\n")
        dense_manifest = work / "dense.manifest.json"
        dense_result = subprocess.run(
            [str(BINARY), "-I", str(dense_input), "--output-manifest", str(dense_manifest)],
            text=True, capture_output=True,
        )
        assert dense_result.returncode == 0, (dense_result.stdout, dense_result.stderr)
        dense_index = Path(str(dense_input) + ".idx")
        assert dense_index.is_file() and dense_index.stat().st_size > 0
        dense_payload = json.loads(dense_manifest.read_text())
        assert dense_payload["index_type"] == "TRIBBLE_INTERVAL_TREE"
        assert dense_payload["linear_records"] == 300
        assert dense_payload["linear_blocks"] == 4
        if JAVA.is_file() and GATK_JAR.is_file():
            java_dense_index = work / "java-dense.idx"
            java_dense_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "IndexFeatureFile",
                 "-I", str(dense_input), "-O", str(java_dense_index)],
                text=True, capture_output=True,
            )
            assert java_dense_result.returncode == 0, java_dense_result.stderr
            assert dense_index.read_bytes() == java_dense_index.read_bytes()
            dense_query = work / "dense-query.vcf"
            dense_query_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "SelectVariants",
                 "-V", str(dense_input), "-L", "1:250-250", "-O", str(dense_query),
                 "--create-output-variant-index", "false"],
                text=True, capture_output=True,
            )
            assert dense_query_result.returncode == 0, dense_query_result.stderr
            assert any(line.startswith("1\t250\t") for line in dense_query.read_text().splitlines())

        # Uncompressed BED is also a Tribble feature stream.  HTSJDK uses the
        # same DynamicIndexCreator/LinearIndex v3 contract (with BEDCodec's
        # zero-based, half-open coordinates converted to 1-based inclusive
        # Feature intervals); native output must retain the .idx default and
        # its source offsets/properties.
        plain_bed_input = work / "sample.bed"
        plain_bed_input.write_bytes(PLAIN_BED.read_bytes())
        plain_bed_result = subprocess.run(
            [str(BINARY), "-I", str(plain_bed_input)], text=True, capture_output=True,
        )
        assert plain_bed_result.returncode == 0, (plain_bed_result.stdout, plain_bed_result.stderr)
        plain_bed_index = Path(str(plain_bed_input) + ".idx")
        assert plain_bed_index.is_file() and plain_bed_index.stat().st_size > 0
        plain_bed_summary = json.loads(plain_bed_result.stderr.splitlines()[-1])
        assert plain_bed_summary["index_type"] == "TRIBBLE_LINEAR"
        assert plain_bed_summary["linear_records"] == 6
        if JAVA.is_file() and GATK_JAR.is_file():
            java_bed_index = work / "java-bed.idx"
            java_bed_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "IndexFeatureFile",
                 "-I", str(plain_bed_input), "-O", str(java_bed_index)],
                text=True, capture_output=True,
            )
            assert java_bed_result.returncode == 0, java_bed_result.stderr
            assert plain_bed_index.read_bytes() == java_bed_index.read_bytes()

        dense_bed_input = work / "dense.bed"
        dense_bed_input.write_text("".join(
            f"1\t{position - 1}\t{position}\n" for position in range(1, 301)))
        dense_bed_result = subprocess.run(
            [str(BINARY), "-I", str(dense_bed_input)], text=True, capture_output=True,
        )
        assert dense_bed_result.returncode == 0, (dense_bed_result.stdout, dense_bed_result.stderr)
        dense_bed_index = Path(str(dense_bed_input) + ".idx")
        assert dense_bed_index.is_file() and dense_bed_index.stat().st_size > 0
        dense_bed_summary = json.loads(dense_bed_result.stderr.splitlines()[-1])
        assert dense_bed_summary["index_type"] == "TRIBBLE_INTERVAL_TREE"
        if JAVA.is_file() and GATK_JAR.is_file():
            java_dense_bed_index = work / "java-dense-bed.idx"
            java_dense_bed_result = subprocess.run(
                [str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "IndexFeatureFile",
                 "-I", str(dense_bed_input), "-O", str(java_dense_bed_index)],
                text=True, capture_output=True,
            )
            assert java_dense_bed_result.returncode == 0, java_dense_bed_result.stderr
            assert dense_bed_index.read_bytes() == java_dense_bed_index.read_bytes()

        empty_input = work / "empty.vcf"
        empty_input.write_bytes(EMPTY_VCF.read_bytes())
        empty_result = subprocess.run(
            [str(BINARY), "-I", str(empty_input)], text=True, capture_output=True,
        )
        assert empty_result.returncode == 0, (empty_result.stdout, empty_result.stderr)
        empty_index = Path(str(empty_input) + ".idx")
        assert empty_index.is_file() and empty_index.stat().st_size > 0

        bcf_index = work / "sample.bcf.csi"
        bcf_result = subprocess.run(
            [str(BINARY), "-I", str(BGZF_BCF), "-O", str(bcf_index), "--threads", "2"],
            text=True, capture_output=True,
        )
        assert bcf_result.returncode == 0, (bcf_result.stdout, bcf_result.stderr)
        assert bcf_index.is_file() and bcf_index.stat().st_size > 0

        bed_index = work / "sample.bed.gz.tbi"
        bed_result = subprocess.run(
            [str(BINARY), "-I", str(BGZF_BED), "-O", str(bed_index), "--threads", "2"],
            text=True, capture_output=True,
        )
        assert bed_result.returncode == 0, (bed_result.stdout, bed_result.stderr)
        assert bed_index.is_file() and bed_index.stat().st_size > 0

    print(json.dumps({"status": "pass", "suite": "index-feature-file", "binary": str(BINARY),
                      "plain_gvcf_linear_idx": True, "plain_vcf_dynamic_idx": True,
                      "dense_vcf_interval_idx": True, "plain_bed_dynamic_idx": True,
                      "dense_bed_interval_idx": True,
                      "empty_vcf_dynamic_idx": True,
                      "java_binary_oracle": JAVA.is_file() and GATK_JAR.is_file()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
