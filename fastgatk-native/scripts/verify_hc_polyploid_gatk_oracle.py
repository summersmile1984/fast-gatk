#!/usr/bin/env python3
"""Triploid HC/GVCF oracle for the shared reference-confidence path.

This is deliberately a focused oracle rather than a whole-file comparison:
GATK and the native writer have different provenance/header metadata and may
choose different block boundaries when a low-depth site changes its GQ band.
The records at the first reference block and the pinned candidate site are
stable, and protect the ploidy-specific Java likelihood/indel merge plus the
bounded MLEAC calculation.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n").split("\t") for line in stream
                if line and not line.startswith("#")]


def by_pos(rows: list[list[str]], pos: int) -> list[str]:
    return next(row for row in rows if int(row[1]) == pos)


def sample(row: list[str]) -> dict[str, str]:
    return dict(zip(row[8].split(":"), row[9].split(":")))


def info(row: list[str]) -> dict[str, str]:
    return {item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in row[7].split(";") if "=" in item}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or triploid HC fixture")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_hc_polyploid_gatk_oracle.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-polyploid-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.g.vcf.gz"
        native_vcf = work / "native.g.vcf.gz"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", "17:69000-69100",
            "-O", str(gatk_vcf), "-ERC", "GVCF", "--sample-ploidy", "3",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(native_vcf), "-ERC", "GVCF", "--sample-ploidy", "3",
            "--min-depth", "1", "--min-alt-support", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gatk_rows = records(gatk_vcf)
        native_rows = records(native_vcf)
        gatk_block = by_pos(gatk_rows, 69000)
        native_block = by_pos(native_rows, 69000)
        # The Java reference-confidence record is text-identical at this
        # sentinel (including the triploid indel-model PL vector).
        assert "\t".join(gatk_block) == "\t".join(native_block), (
            gatk_block, native_block)
        block_sample = sample(native_block)
        assert block_sample["GT"] == "0/0/0"
        assert block_sample["PL"] == "0,5,14,135"

        gatk_candidate = by_pos(gatk_rows, 69067)
        native_candidate = by_pos(native_rows, 69067)
        gatk_sample = sample(gatk_candidate)
        native_sample = sample(native_candidate)
        for field in ("GT", "AD", "DP", "GQ", "PL"):
            assert native_sample.get(field) == gatk_sample.get(field), (
                field, gatk_sample, native_sample)
        gatk_info = info(gatk_candidate)
        native_info = info(native_candidate)
        assert "ExcessHet" not in native_info
        # MLEAC/MLEAF are the bounded EM estimate over the full concrete +
        # <NON_REF> genotype matrix, not counts from the selected GT.
        assert native_info.get("MLEAC") == gatk_info.get("MLEAC")
        assert native_info.get("MLEAF") == gatk_info.get("MLEAF")

        print(json.dumps({
            "status": "pass",
            "sample_ploidy": 3,
            "reference_block_69000_exact": True,
            "reference_block_pl": block_sample["PL"],
            "candidate_69067_sample_exact": True,
            "candidate_69067_mleac": native_info.get("MLEAC"),
            "candidate_69067_mleaf": native_info.get("MLEAF"),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
