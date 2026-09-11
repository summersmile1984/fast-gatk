#!/usr/bin/env python3
"""Strict GATK oracle for multiallelic InbreedingCoeff.

This deterministic 20-sample fixture has two concrete ALTs and complete
diploid PL vectors.  It requires exact GATK-compatible row text, including the
multiallelic GenotypeUtils projection and the four-decimal InbreedingCoeff.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return stream.read().splitlines()


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in lines(path)
            if line and not line.startswith("#")]


def header_ids(path: Path) -> set[str]:
    result: set[str] = set()
    for line in lines(path):
        if not line.startswith("##"):
            continue
        for prefix in ("##INFO=<ID=", "##FORMAT=<ID=", "##FILTER=<ID=", "##ALT=<ID="):
            if line.startswith(prefix):
                kind = prefix[2:].split("=", 1)[0]
                result.add(kind + ":" + line[len(prefix):].split(",", 1)[0].split(">", 1)[0])
                break
    return result


def info(record: list[str]) -> dict[str, str]:
    return {item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in record[7].split(";") if "=" in item}


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    combine = Path(os.environ.get(
        "FASTGATK_COMBINE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-combine-gvcfs")))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not combine.exists() or not reference.exists():
        raise SystemExit("missing native build or GATK fixtures; run build_native.sh first")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_gatk_genotype_gvcf_inbreeding.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    sample_values = [
        (31, "3,2,13,13", "38,60,0,38,39,7,69,41,58,72"),
        (21, "12,1,7,1", "0,40,70,62,38,2,71,42,49,75"),
        (22, "6,11,3,2", "87,4,53,61,0,69,42,71,54,70"),
        (28, "9,7,5,7", "72,74,48,12,0,69,62,84,55,64"),
        (32, "13,5,10,4", "0,17,68,66,56,81,63,53,73,39"),
        (33, "15,14,2,2", "56,37,20,0,83,70,71,85,87,55"),
        (37, "12,11,0,14", "78,79,0,39,10,81,79,54,76,71"),
        (36, "15,2,5,14", "50,74,0,66,13,48,84,53,43,82"),
        (15, "4,2,5,4", "57,52,43,0,8,90,70,52,80,61"),
        (26, "11,10,4,1", "44,0,35,66,18,72,46,51,53,35"),
        (17, "3,6,2,6", "41,84,90,0,13,70,60,60,60,60"),
        (23, "3,2,6,12", "41,42,20,0,38,41,35,71,44,69"),
        (18, "9,2,4,3", "64,0,57,73,16,65,42,42,89,66"),
        (22, "3,9,2,8", "68,82,5,65,88,0,45,68,36,48"),
        (32, "6,7,12,7", "75,20,57,84,0,69,69,84,67,56"),
        (38, "14,11,11,2", "47,0,57,81,36,12,85,52,65,51"),
        (24, "3,3,12,6", "35,0,65,47,12,48,65,74,74,88"),
        (13, "5,4,0,4", "6,62,85,0,56,40,86,81,60,64"),
        (7, "4,0,0,3", "44,86,76,44,0,18,73,65,77,57"),
        (38, "7,10,8,13", "48,43,17,90,0,87,90,48,36,51"),
    ]
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-inbreeding-oracle-") as directory:
        work = Path(directory)
        source = work / "source.g.vcf"
        input_gvcf = work / "input.g.vcf.gz"
        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf.gz"
        manifest_path = work / "native.manifest.json"
        samples = "\t".join(f"S{index:02d}" for index in range(len(sample_values)))
        calls = "\t".join(
            f"./.:{depth}:{ad}:{pl}" for depth, ad, pl in sample_values)
        source.write_text(
            """##fileformat=VCFv4.2
##reference=human_g1k_v37.chr17_1Mb.fasta
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"""
            + samples + "\n17\t69067\t.\tT\tG,A,<NON_REF>\t.\tPASS\tDP=200\tGT:DP:AD:PL\t"
            + calls + "\n", encoding="utf-8")
        run([str(combine), "-V", str(source), "-O", str(input_gvcf)],
            "native CombineGVCFs fixture serialization")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(input_gvcf), "-O", str(gatk_output),
            "--create-output-variant-index", "false", "--seconds-between-progress-updates", "1",
        ], "GATK GenotypeGVCFs")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(input_gvcf),
            "-O", str(native_output), "--create-output-variant-index", "false",
            "--gatk-compatible-annotations", "--output-manifest", str(manifest_path),
        ], text=True, capture_output=True)
        if native_result.returncode != 0:
            raise RuntimeError(f"native GenotypeGVCFs failed:\n{native_result.stderr[-4000:]}")
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        gatk_rows = rows(gatk_output)
        native_rows = rows(native_output)
        gatk_ids = header_ids(gatk_output)
        native_ids = header_ids(native_output)

    assert len(gatk_rows) == 1 and native_rows == gatk_rows, (
        f"row mismatch: gatk={gatk_rows!r} native={native_rows!r}")
    required = {
        "INFO:AC", "INFO:AF", "INFO:AN", "INFO:DP", "INFO:ExcessHet",
        "INFO:InbreedingCoeff", "INFO:MLEAC", "INFO:MLEAF", "INFO:QD",
        "FORMAT:GT", "FORMAT:AD", "FORMAT:DP", "FORMAT:GQ", "FORMAT:PL",
    }
    assert required <= gatk_ids and required <= native_ids
    row = native_rows[0]
    annotations = info(row)
    assert row[3:5] == ["T", "G,A"]
    assert annotations["InbreedingCoeff"] == "0.0022"
    assert annotations["ExcessHet"] == "1.5298"
    assert len(row) == 29 and len(row[9:]) == 20
    assert native_summary["status"] == "contract-compatible"
    assert manifest["telemetry"]["sample_count"] == 20
    assert manifest["telemetry"]["cohort_af_kernel_calls"] > 0
    print(json.dumps({
        "status": "pass",
        "records": len(native_rows),
        "sample_count": len(row) - 9,
        "record_text_exact": True,
        "inbreeding_coeff_exact": annotations["InbreedingCoeff"],
        "multiallelic_excess_het_exact": annotations["ExcessHet"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
