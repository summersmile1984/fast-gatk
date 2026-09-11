#!/usr/bin/env python3
"""Strict GenotypeGVCFs oracle for a concrete multi-allelic cohort site.

The fixture has two concrete ALTs plus ``<NON_REF>`` and complete diploid
Number=G PL vectors (10 entries for four alleles).  It exercises allele-union
projection, cohort EM/QUAL, MLEAC/MLEAF, and the multi-allelic ExcessHet
genotype-count reducer rather than only the common biallelic path.
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


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in lines(path)
            if line and not line.startswith("#")]


def header(path: Path) -> list[str]:
    ignored = ("##GATKCommandLine=", "##source=", "##fileDate=",
               "##contig=", "##GVCFBlock", "##fastgatk_genotype_gvcfs_status=")
    return sorted({line for line in lines(path)
                   if line.startswith("##") and not line.startswith(ignored)})


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
        oracle_guard.oracle_not_verified('verify_gatk_genotype_gvcf_multiallelic.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-multiallelic-oracle-") as directory:
        work = Path(directory)
        source = work / "source.g.vcf"
        input_gvcf = work / "input.g.vcf.gz"
        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf.gz"
        native_manifest = work / "native.manifest.json"
        source.write_text(
            """##fileformat=VCFv4.2
##reference=human_g1k_v37.chr17_1Mb.fasta
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=\"Represents any possible alternate allele\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
17\t69067\t.\tT\tG,A,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t0/1:10:5,5,0,0:50,0,50,99,99,99,99,99,99,99,99\t1/2:10:0,4,6,0:90,80,70,60,50,0,99,99,99,99
""", encoding="utf-8")
        # Native CombineGVCFs supplies BGZF+tabix serialization for this
        # deterministic fixture; the oracle under test is GenotypeGVCFs.
        run([str(combine), "-V", str(source), "-O", str(input_gvcf)],
            "native CombineGVCFs fixture serialization")
        run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(input_gvcf),
            "-O", str(gatk_output), "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], "GATK GenotypeGVCFs")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(input_gvcf),
            "-O", str(native_output), "--create-output-variant-index", "false",
            "--gatk-compatible-annotations", "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True)
        if native_result.returncode != 0:
            raise RuntimeError(f"native GenotypeGVCFs failed:\n{native_result.stderr[-4000:]}")
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        gatk_rows = records(gatk_output)
        native_rows = records(native_output)
        gatk_header_ids = header_ids(gatk_output)
        native_header_ids = header_ids(native_output)

    assert len(gatk_rows) == 1 and native_rows == gatk_rows, (
        f"row mismatch: gatk={gatk_rows!r} native={native_rows!r}")
    required_header_ids = {
        "INFO:AC", "INFO:AF", "INFO:AN", "INFO:DP", "INFO:ExcessHet",
        "INFO:MLEAC", "INFO:MLEAF", "INFO:QD", "FORMAT:GT", "FORMAT:AD",
        "FORMAT:DP", "FORMAT:GQ", "FORMAT:PL",
    }
    assert required_header_ids <= gatk_header_ids
    assert required_header_ids <= native_header_ids
    row = native_rows[0]
    assert row[3:5] == ["T", "G,A"]
    assert row[8].split(":") == ["GT", "AD", "DP", "GQ", "PL"]
    assert row[9].split(":")[0] == "0/1"
    assert row[10].split(":")[0] == "2/2"
    assert len(row[9].split(":")[-1].split(",")) == 6
    assert len(row[10].split(":")[-1].split(",")) == 6
    annotations = info(row)
    assert annotations["AC"] == "1,2" and annotations["AN"] == "4"
    assert annotations["AF"] == "0.250,0.500"
    assert annotations["MLEAC"] == "1,2" and annotations["MLEAF"] == "0.250,0.500"
    assert annotations["ExcessHet"] == "0.0000"
    assert native_summary["status"] == "contract-compatible"
    assert manifest["telemetry"]["sample_count"] == 2
    assert manifest["telemetry"]["cohort_af_kernel_calls"] > 0
    assert manifest["telemetry"]["cohort_af_converged"] > 0
    print(json.dumps({
        "status": "pass",
        "records": len(native_rows),
        "record_text_exact": True,
        "header_required_ids_exact": True,
        "concrete_alt_count": 2,
        "number_g_pl_entries": 6,
        "excess_het_exact": True,
        "cohort_iterations": manifest["telemetry"]["cohort_af_iterations"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
