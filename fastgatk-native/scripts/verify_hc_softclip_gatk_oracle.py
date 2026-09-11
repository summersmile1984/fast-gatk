#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for HC soft-clip assembly evidence.

The chr17 fixture contains a non-terminal mismatch near clipped reads.  GATK's
``--dont-use-soft-clipped-bases`` keeps that aligned SNP while removing the
terminal-clip assembly evidence.  This verifier protects the read-threading
call boundary and all stable genotype/annotation fields that are independent
of the still-evolving Coverage INFO/DP population.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in
            path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def key(record: list[str]) -> tuple[str, int, str, str]:
    return record[0], int(record[1]), record[3], record[4]


def info(record: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in record[7].split(";"):
        if "=" in item:
            name, value = item.split("=", 1)
            values[name] = value
    return values


def sample(record: list[str]) -> dict[str, str]:
    return dict(zip(record[8].split(":"), record[9].split(":")))


def stable_header(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("#") and not line.startswith("##GATKCommandLine=")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    region = "17:69000-70000"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or soft-clip HC fixture")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_hc_softclip_gatk_oracle.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-softclip-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_vcf),
            "-L", region, "--dont-use-soft-clipped-bases",
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_vcf), "--output-manifest", str(manifest),
            "--threads", "2", "--dont-use-soft-clipped-bases",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gatk_records = records(gatk_vcf)
        native_records = records(native_vcf)
        gatk_keys = {key(record) for record in gatk_records}
        native_keys = {key(record) for record in native_records}
        assert gatk_keys == native_keys, (
            f"soft-clip HC call-set mismatch: gatk-only={sorted(gatk_keys - native_keys)} "
            f"native-only={sorted(native_keys - gatk_keys)}")
        expected = {
            ("17", 69067, "T", "G"),
            ("17", 69298, "A", "T"),
            ("17", 69368, "G", "C"),
            ("17", 69631, "C", "T"),
        }
        assert gatk_keys == expected, f"unexpected pinned GATK call set: {sorted(gatk_keys)}"
        assert stable_header(gatk_vcf) == stable_header(native_vcf)

        gatk_by_key = {key(record): record for record in gatk_records}
        native_by_key = {key(record): record for record in native_records}
        info_dp = {}
        for variant in sorted(gatk_keys):
            gatk = gatk_by_key[variant]
            native_record = native_by_key[variant]
            assert gatk[:7] == native_record[:7], f"stable VCF columns differ at {variant}"
            gatk_sample = sample(gatk)
            native_sample = sample(native_record)
            for field in ("GT", "AD", "DP", "GQ", "PL"):
                assert gatk_sample.get(field) == native_sample.get(field), (
                    f"FORMAT/{field} differs at {variant}: "
                    f"gatk={gatk_sample.get(field)} native={native_sample.get(field)}")
            gatk_info = info(gatk)
            native_info = info(native_record)
            for field, value in gatk_info.items():
                if field == "DP":
                    info_dp[str(variant[1])] = {
                        "gatk": value, "native": native_info.get(field)}
                    continue
                assert native_info.get(field) == value or (
                    field.endswith("RankSum") and
                    math.isclose(float(native_info.get(field, "nan")), float(value), abs_tol=0.0005)
                ), f"INFO/{field} differs at {variant}: gatk={value} native={native_info.get(field)}"

        manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = manifest_data["telemetry"]
        assert telemetry["use_soft_clipped_bases"] is False
        assert telemetry["softclip_candidates"] == 0
        assert telemetry["candidate_softclip_suppressed"] == 0

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "dont_use_soft_clipped_bases": True,
            "calls": len(gatk_records),
            "call_set_exact": True,
            "stable_genotype_fields_exact": True,
            "stable_annotations_exact_except_info_dp": True,
            "info_dp_comparison": info_dp,
            "scope": "read-threading/call fields; Coverage INFO/DP population remains outside",
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
