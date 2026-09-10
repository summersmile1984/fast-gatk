#!/usr/bin/env python3
"""Pinned GATK oracle for HC --min-base-quality-score assembly semantics.

GATK wires this option into both the pileup/calling threshold and
ReadThreadingAssembler's per-base graph segmentation.  In the canonical chr17
fixture, Q20 removes the 17:69067 T>G graph branch: its second supporting
fragment is split at an interior Q14 base, so the remaining branch has only
one observation and --min-pruning=2 removes it.  This is deliberately an
assembly-boundary check; PairHMM and genotyping remain in the Kokkos kernels.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def records(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in
            path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def key(record: list[str]) -> tuple[str, int, str, str]:
    return record[0], int(record[1]), record[3], record[4]


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
        raise SystemExit("missing native build or HC min-base-quality fixture")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    expected = {
        ("17", 69368, "G", "C"),
        ("17", 69631, "C", "T"),
    }
    common = [
        "-R", str(reference), "-I", str(bam), "-L", region,
        "--min-base-quality-score", "20",
        "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
    ]
    cases = (("default", []), ("overlap-correction-disabled", [
        "--do-not-correct-overlapping-quality"]))
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-min-base-quality-") as directory:
        work = Path(directory)
        for label, extra in cases:
            gatk_vcf = work / f"gatk.{label}.vcf"
            native_vcf = work / f"native.{label}.vcf"
            manifest = work / f"native.{label}.json"
            subprocess.run([
                java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
                *common, *extra, "-O", str(gatk_vcf),
                "--native-pair-hmm-threads", "2",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run([
                str(native), *common, *extra, "-O", str(native_vcf),
                "--threads", "2", "--output-manifest", str(manifest),
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
                f"Q20 HC VCF differs from GATK for {label}")
            gatk_keys = {key(record) for record in records(gatk_vcf)}
            native_keys = {key(record) for record in records(native_vcf)}
            assert gatk_keys == expected, (
                f"unexpected GATK Q20 call set for {label}: {sorted(gatk_keys)}")
            assert native_keys == expected, (
                f"unexpected native Q20 call set for {label}: {sorted(native_keys)}")
            assert all(position != 69067 for _, position, _, _ in native_keys)
            telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
            assert telemetry["min_base_quality"] == 20
            assert telemetry["graph_used"] is True
            assert telemetry["pairhmm_used"] is True

        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0", "region": region,
            "min_base_quality": 20, "cases": [case[0] for case in cases],
            "vcf_exact": True, "suppressed_position": 69067,
            "host_graph_segmentation": True, "kokkos_pairhmm_genotype": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
