#!/usr/bin/env python3
"""Compare ReblockGVCF on two non-overlapping shards of one sample.

GATK ReblockGVCF deliberately accepts multiple input files only when they are
non-overlapping shards from the same sample; it rejects a multi-sample VCF.
This oracle covers the supported multi-input boundary instead of claiming
multi-sample support that the Java tool does not provide.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=All filters passed>
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of the variant reference block>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Approximate read depth>
##INFO=<ID=RAW_MQandDP,Number=2,Type=Float,Description=Raw MQ and DP values>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depth>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def read_records(path: Path) -> list[dict[str, object]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    records: list[dict[str, object]] = []
    sample_names: list[str] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, str | bool] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    info[key] = value if separator else True
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            samples = {
                name: dict(zip(format_keys, value.split(":")))
                for name, value in zip(sample_names, fields[9:])
            }
            records.append({
                "site": tuple(fields[index] for index in range(7)),
                "info": dict(sorted(info.items())),
                "samples": samples,
            })
    return records


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
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
    if not all(path.is_file() for path in (java, gatk, native, reference)):
        oracle_guard.oracle_not_verified('verify_reblock_gatk_shards.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK shard oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-gatk-shards-") as directory:
        work = Path(directory)
        shard1 = work / "shard1.g.vcf"
        shard2 = work / "shard2.g.vcf"
        multisample = work / "multisample.g.vcf"
        gatk_output = work / "gatk.reblocked.g.vcf.gz"
        native_output = work / "native.reblocked.g.vcf.gz"
        native_manifest = work / "native.reblocked.manifest.json"
        shard1_body = (
            "17\t69001\t.\tA\t<NON_REF>\t.\t.\tEND=69005\t"
            "GT:DP:PL:GQ\t0/0:10:0,10,200:10\n"
            # A low-quality deletion is converted by both tools using PL[0]
            # (not the stale GQ field), while preserving the original END span.
            "17\t69006\t.\tCT\tC,<NON_REF>\t.\tPASS\tDP=20;RAW_MQandDP=72000,20\t"
            "GT:DP:AD:PL:GQ\t0/1:20:12,8,0:1,0,20,99,99,99:99\n"
        )
        shard2_body = (
            "17\t69008\t.\tT\t<NON_REF>\t.\t.\tEND=69010\t"
            "GT:DP:PL:GQ\t0/0:8:0,15,200:15\n"
            "17\t69011\t.\tA\tT,G,<NON_REF>\t.\tPASS\tDP=18;RAW_MQandDP=64800,18\t"
            "GT:DP:AD:PL:GQ\t0/1:18:10,6,0,0:60,0,90,99,99,99,99,99,99,99:30\n"
        )
        shard1.write_text(HEADER + shard1_body, encoding="ascii")
        shard2.write_text(HEADER + shard2_body, encoding="ascii")
        multisample_header = HEADER.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
        )
        multisample.write_text(
            multisample_header
            + "17\t69001\t.\tA\t<NON_REF>\t.\t.\tEND=69005\t"
            "GT:DP:PL:GQ\t0/0:10:0,10,200:10\t0/0:8:0,15,200:15\n",
            encoding="ascii",
        )
        run([java, "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile", "-I", str(shard1)],
            "GATK IndexFeatureFile shard1")
        run([java, "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile", "-I", str(shard2)],
            "GATK IndexFeatureFile shard2")
        run([java, "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile", "-I", str(multisample)],
            "GATK IndexFeatureFile multisample")
        common = ["-R", str(reference), "-GQB", "20", "-GQB", "100",
                  "--rgq-threshold-to-no-call", "10",
                  "--create-output-variant-index", "false"]
        run([
            java, "-Xmx1g", "-jar", str(gatk), "ReblockGVCF", *common,
            "-V", str(shard1), "-V", str(shard2), "-O", str(gatk_output),
        ], "GATK ReblockGVCF shards")
        native_result = subprocess.run([
            str(native), *common, "-V", str(shard1), "-V", str(shard2),
            "-O", str(native_output), "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True, check=False)
        if native_result.returncode != 0:
            raise RuntimeError(f"native ReblockGVCF shards failed:\n{native_result.stderr[-4000:]}")
        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert native_records == gatk_records, (
            "native/GATK shard records differ", native_records, gatk_records,
        )
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))

        # This is an intentional negative compatibility check: Java GATK
        # rejects a multi-sample ReblockGVCF input during traversal startup.
        rejected = subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "ReblockGVCF", *common,
            "-V", str(multisample), "-O", str(work / "gatk.multisample.g.vcf.gz"),
        ], text=True, capture_output=True, check=False)
        rejection_text = (rejected.stdout + rejected.stderr).lower()
        assert rejected.returncode != 0
        assert "found samples" in rejection_text or "single-sample" in rejection_text

    assert manifest["telemetry"]["input_files"] == 2
    assert manifest["telemetry"]["sample_count"] == 1
    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "records": len(native_records),
        "input_shards": 2,
        "sample_count": 1,
        "record_fields_exact": True,
        "gatk_multisample_rejected": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
