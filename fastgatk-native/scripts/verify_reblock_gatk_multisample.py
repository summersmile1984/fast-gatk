#!/usr/bin/env python3
"""Compare native ReblockGVCF with GATK on a deterministic two-sample GVCF.

The existing ReblockGVCF contract test exercises multi-sample data locally.  This
oracle adds the missing Java comparison for sample-major GQ/DP/AD/PL handling,
while keeping the fixture small enough for every CI backend.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=All filters passed>
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of the variant reference block>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Approximate read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depth>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def read_records(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    samples: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, str | bool] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    info[key] = value if separator else True
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            sample_values = {
                name: dict(zip(format_keys, value.split(":")))
                for name, value in zip(samples, fields[9:])
            }
            records.append({
                "site": tuple(fields[index] for index in range(7)),
                "info": dict(sorted(info.items())),
                "samples": sample_values,
            })
    return records


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    gatk = Path(os.environ.get(
        "FASTGATK_GATK_JAR",
        str(root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"),
    ))
    native = Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY",
        str(root / "fastgatk-native/build/fastgatk-reblock-gvcf"),
    ))
    reference = Path(os.environ.get(
        "FASTGATK_REBLOCK_REFERENCE",
        str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
    ))
    required = (Path(java), gatk, native, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK multisample oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-gatk-multisample-") as directory:
        work = Path(directory)
        # GATK requires an index for feature traversal.  Keep this synthetic
        # fixture plain VCF so IndexFeatureFile can create the portable .idx
        # sidecar without depending on an external tabix executable.
        source = work / "multisample.g.vcf"
        gatk_output = work / "gatk.reblocked.g.vcf.gz"
        native_output = work / "native.reblocked.g.vcf.gz"
        native_manifest = work / "native.reblocked.manifest.json"
        body = (
            # Two adjacent blocks with independent per-sample GQ/DP values.
            "17\t69001\t.\tA\t<NON_REF>\t.\tPASS\tEND=69005\t"
            "GT:DP:AD:PL:GQ\t0/0:10:10,0:0,10,200:10\t"
            "0/0:8:8,0:0,15,200:15\n"
            "17\t69006\t.\tC\t<NON_REF>\t.\tPASS\tEND=69010\t"
            "GT:DP:AD:PL:GQ\t0/0:8:8,0:0,15,200:15\t"
            "0/0:7:7,0:0,12,200:12\n"
            # A concrete multi-ALT record exercises sample-major Number=R/G.
            "17\t69012\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=38\t"
            "GT:DP:AD:PL:GQ\t0/1:20:12,8,0,0:50,0,80,99,99,99,99,99,99:20\t"
            "0/1:18:10,6,0,0:60,0,90,99,99,99,99,99,99:30\n"
        )
        with source.open("wt", encoding="ascii") as stream:
            stream.write(HEADER)
            stream.write(body)

        run([
            java, "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile", "-I", str(source),
        ], "GATK IndexFeatureFile multisample")

        common = ["-GQB", "20", "-GQB", "100"]
        run([
            java, "-Xmx1g", "-jar", str(gatk), "ReblockGVCF", "-R", str(reference),
            "-V", str(source), *common, "--create-output-variant-index", "false",
            "-O", str(gatk_output),
        ], "GATK ReblockGVCF multisample")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(source), *common,
            "--create-output-variant-index", "false", "-O", str(native_output),
            "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True, check=False)
        if native_result.returncode != 0:
            raise RuntimeError(f"native ReblockGVCF multisample failed:\n{native_result.stderr[-4000:]}")
        native_records = read_records(native_output)
        gatk_records = read_records(gatk_output)
        assert native_records == gatk_records, (
            "native/GATK multisample records differ", native_records, gatk_records,
        )
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))

    assert len(native_records) == 2
    assert manifest["compatibility"]["multi_sample"] is True
    assert manifest["telemetry"]["sample_count"] == 2
    assert manifest["telemetry"]["pl_remap_kernel_calls"] > 0
    assert manifest["telemetry"]["allele_field_remap_kernel_calls"] > 0
    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "records": len(native_records),
        "sample_count": manifest["telemetry"]["sample_count"],
        "sample_major_fields_exact": True,
        "record_fields_exact": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
