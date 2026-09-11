#!/usr/bin/env python3
"""Strict Java/native oracle for GenotypeGVCFs' legacy AF-calculator mode.

HaplotypeCaller gVCFs already contain reducible MLEAC/MLEAF values for the
concrete ALT and ``<NON_REF>``.  GenotypeGVCFs with
``--use-new-qual-calculator false`` does not run the new EM calculator, but it
still projects those Number=A fields when materializing the final ALT list.
This guard keeps the native legacy path from leaking the symbolic allele's
second value into the output (for example ``MLEAC=1,0``).
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return stream.read().splitlines()


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in read_lines(path)
            if line and not line.startswith("#")]


def semantic_header(path: Path) -> list[str]:
    ignored = ("##GATKCommandLine=", "##source=", "##fileDate=",
               "##contig=", "##GVCFBlock", "##fastgatk_genotype_gvcfs_status=")
    return sorted({line for line in read_lines(path)
                   if line.startswith("##") and not line.startswith(ignored)})


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or GATK fixtures; run build_native.sh first")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_gatk_genotype_gvcf_legacy_qual.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    region = "17:69000-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-legacy-qual-oracle-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf.gz"
        gatk_output = work / "gatk.legacy.vcf.gz"
        native_output = work / "native.legacy.vcf.gz"
        native_manifest = work / "native.legacy.manifest.json"
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-O", str(gatk_gvcf), "-ERC", "GVCF",
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "true",
            "--seconds-between-progress-updates", "1",
        ], "GATK HaplotypeCaller")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-L", region,
            "-O", str(gatk_output), "--use-new-qual-calculator", "false",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK GenotypeGVCFs legacy AF calculator")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(gatk_gvcf),
            "-L", region, "-O", str(native_output),
            "--use-new-qual-calculator=false",
            "--gatk-compatible-annotations", "--create-output-variant-index=false",
            "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True)
        if native_result.returncode != 0:
            raise RuntimeError(f"native GenotypeGVCFs legacy path failed:\n{native_result.stderr[-4000:]}")
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        gatk_rows = records(gatk_output)
        native_rows = records(native_output)
        gatk_header = semantic_header(gatk_output)
        native_header = semantic_header(native_output)

    assert gatk_rows, "GATK produced no rows for the legacy AF-calculator fixture"
    assert native_rows == gatk_rows, (
        f"legacy Java/native rows differ:\ngatk={gatk_rows!r}\nnative={native_rows!r}")
    assert native_header == gatk_header, "legacy Java/native semantic headers differ"
    assert native_summary["status"] == "contract-compatible"
    assert manifest["telemetry"]["use_new_qual_calculator"] is False
    assert manifest["telemetry"]["cohort_af_kernel_calls"] == 0
    assert all("MLEAC=" in row[7] and "MLEAF=" in row[7] for row in native_rows)
    assert all(",0" not in token and ",0.00" not in token
               for row in native_rows for token in row[7].split(";")
               if token.startswith(("MLEAC=", "MLEAF=")))
    print(json.dumps({
        "status": "pass",
        "gatk_records": len(gatk_rows),
        "native_records": len(native_rows),
        "record_text_exact": len(native_rows),
        "header_semantic_exact": True,
        "legacy_mle_projection_exact": True,
        "cohort_af_kernel_calls": manifest["telemetry"]["cohort_af_kernel_calls"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
