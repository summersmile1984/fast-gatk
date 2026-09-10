#!/usr/bin/env python3
"""Pinned GATK oracle for dense reference-only GenotypeGVCFs output.

GATK's ``--include-non-variant-sites`` path is not just a record-count
switch: for monomorphic reference sites it recovers MIN_DP into DP, moves GQ
to RGQ, and drops PL/MIN_DP.  This oracle compares the complete non-header
rows against GATK 4.6.2.0 for both the aggregate and bounded stream paths.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
HC = Path(os.environ.get(
    "FASTGATK_HC_BINARY",
    Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
    / "fastgatk-hc-call",
))
GENOTYPE = Path(os.environ.get(
    "FASTGATK_GENOTYPE_BINARY",
    Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
    / "fastgatk-genotype-gvcf",
))
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
BAM = ROOT / "gatk-source" / "src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source" / "src/test/resources/human_g1k_v37.chr17_1Mb.fasta"


def run(command: list[str], *, expected: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if expected and result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}"
        )
    return result


def body(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle if line and not line.startswith("#")]


def main() -> int:
    for required in (HC, GENOTYPE, JAVA, GATK_JAR, BAM, REFERENCE):
        if not required.exists():
            raise SystemExit(f"missing oracle asset: {required}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-gvcf-include-oracle-") as directory:
        work = Path(directory)
        gvcf = work / "input.g.vcf.gz"
        gatk_output = work / "gatk.vcf.gz"
        run([
            str(HC), "-I", str(BAM), "-R", str(REFERENCE), "-L", "17:69000-69100",
            "-O", str(gvcf), "-ERC", "GVCF", "--min-depth", "1", "--min-alt-support", "1",
        ])
        run([
            str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "GenotypeGVCFs",
            "-R", str(REFERENCE), "-V", str(gvcf), "--include-non-variant-sites",
            "-O", str(gatk_output),
        ])
        gatk_rows = body(gatk_output)
        assert gatk_rows and len(gatk_rows) == 101
        assert all("<NON_REF>" not in row.split("\t")[4] for row in gatk_rows)
        # Every dense reference row in this fixture has the Java cleanup shape.
        assert all(row.split("\t")[8] == "GT:DP:RGQ" for row in gatk_rows if row.split("\t")[4] == ".")

        def native_output(name: str, extra: list[str]) -> tuple[Path, dict]:
            output = work / f"native-{name}.vcf"
            manifest = work / f"native-{name}.manifest.json"
            summary = json.loads("\n".join(run([
                str(GENOTYPE), "-R", str(REFERENCE), "-V", str(gvcf), "-O", str(output),
                "--include-non-variant-sites", "--gatk-compatible-annotations",
                "--output-manifest", str(manifest), *extra,
            ]).stdout.splitlines()[-1:]))
            return output, {"summary": summary, "manifest": json.loads(manifest.read_text(encoding="utf-8"))}

        native, metadata = native_output("aggregate", [])
        assert body(native) == gatk_rows
        native_stream, stream_metadata = native_output("stream", ["--stream-by-locus"])
        assert body(native_stream) == gatk_rows
        for payload in (metadata, stream_metadata):
            assert payload["manifest"]["compatibility"]["gatk_annotation_compatibility"] is True
            assert payload["manifest"]["compatibility"]["include_non_variant_sites"] is True
        assert stream_metadata["summary"]["stream_by_locus"] is True

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "records": len(gatk_rows),
            "aggregate_rows_exact": True,
            "stream_rows_exact": True,
            "reference_cleanup_shape_exact": True,
            "min_dp_to_dp_rgq_cleanup_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
