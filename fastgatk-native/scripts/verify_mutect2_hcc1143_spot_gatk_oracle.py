#!/usr/bin/env python3
"""Strict real-data Mutect2 oracles at two HCC1143 chr20 sites.

These are deliberately narrow, fully pinned real tumor/normal regions rather
than synthetic fixtures. They lock complete emitted VCF rows, including
GATK sample ordering, EventMap ECNT/ECNTH values, and physical phasing fields,
against GATK 4.6.2.0. The broader chr20 integration fixture
remains responsible for the long-running traversal/resource check.
"""

from __future__ import annotations

import gzip
import difflib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "real" / "cnv_somatic" / "human_g1k_v37.chr-20.truncated.fasta"
TUMOR = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_tumor.bam"
NORMAL = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_normal.bam"
TUMOR_SAMPLE = "HCC1143"
NORMAL_SAMPLE = "HCC1143 BL"
ORACLES = (
    ("singleton", "20:67000-69000", (("20", "68037", "G", "T"),)),
    ("phased_event_cluster", "20:157000-162000", (
        ("20", "159395", "C", "A"),
        ("20", "159409", "G", "A"),
        ("20", "159443", "T", "A"),
        ("20", "159459", "G", "A"),
        ("20", "161202", "C", "G"),
    )),
)
PAIRHMM_HAPLOTYPE_ORACLE = "phased_event_cluster"


def run(
    command: list[str], env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False, env=env)


def vcf_header_and_rows(path: Path) -> tuple[str, list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]
    header = next(line for line in lines if line.startswith("#CHROM"))
    return header, [line for line in lines if line and not line.startswith("#")]


def vcf_schema_header(path: Path) -> list[str]:
    """Return VCF metadata plus #CHROM, excluding execution provenance only."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line.startswith("#") and not line.startswith("##GATKCommandLine=")]


def gatk_pairhmm_haplotypes(path: Path) -> set[str]:
    """Return the haplotype column emitted by GATK's PairHMM debug writer."""
    haplotypes: set[str] = set()
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            haplotypes.add(line.split(maxsplit=1)[0])
    return haplotypes


def native_shared_pairhmm_haplotypes(stderr: str) -> set[str]:
    """Return paths from the final shared AssemblyResult PairHMM passes.

    Native performs a preliminary tumor-only assembly to construct the joint
    EventMap. Those provisional likelihood calls are intentionally not part
    of Mutect2Engine's final shared-haplotype likelihood matrix. The debug
    boundary labels the latter with ``external_graph_event_map=1``.
    """
    haplotypes: set[str] = set()
    shared_graph_pass = False
    request_prefix = "[FASTGATK_PAIRHMM_REQUEST] "
    for line in stderr.splitlines():
        if line.startswith("[FASTGATK_PAIRHMM_REQUESTS_BEGIN]"):
            shared_graph_pass = "external_graph_event_map=1" in line
        elif line.startswith("[FASTGATK_PAIRHMM_REQUESTS_END]"):
            shared_graph_pass = False
        elif shared_graph_pass and line.startswith(request_prefix):
            haplotypes.add(line[len(request_prefix):].split("\t", maxsplit=1)[0])
    return haplotypes


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError(f"missing HCC1143 spot oracle inputs: {missing}")
        print(json.dumps({"status": "skip", "reason": "HCC1143 spot inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-hcc1143-spot-") as directory:
        work = Path(directory)
        for name, interval, expected_keys in ORACLES:
            gatk_vcf = work / f"gatk-{name}.vcf.gz"
            native_vcf = work / f"native-{name}.vcf.gz"
            gatk_command = [
                str(JAVA), "-jar", str(GATK), "Mutect2",
                "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
                "-tumor", TUMOR_SAMPLE, "-normal", NORMAL_SAMPLE, "-L", interval,
                "-O", str(gatk_vcf),
            ]
            gatk_pairhmm = work / f"gatk-{name}.pairhmm.txt"
            if name == PAIRHMM_HAPLOTYPE_ORACLE:
                gatk_command.extend(["--pair-hmm-results-file", str(gatk_pairhmm)])
            gatk_run = run(gatk_command)
            assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
            native_command = [
                str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
                "--tumor-sample", TUMOR_SAMPLE, "--normal-sample", NORMAL_SAMPLE,
                "-L", interval, "-O", str(native_vcf), "--native-pair-hmm-threads", "1",
            ]
            native_env = None
            if name == PAIRHMM_HAPLOTYPE_ORACLE:
                native_env = os.environ | {"FASTGATK_DEBUG_PAIRHMM_REQUESTS": "1"}
            native_run = run(native_command, native_env)
            assert native_run.returncode == 0, native_run.stderr[-4000:]

            gatk_header, gatk_rows = vcf_header_and_rows(gatk_vcf)
            native_header, native_rows = vcf_header_and_rows(native_vcf)
            assert gatk_header == native_header, (name, gatk_header, native_header)
            assert vcf_schema_header(gatk_vcf) == vcf_schema_header(native_vcf), (
                name, "VCF schema header mismatch")
            assert len(gatk_rows) == len(native_rows) == len(expected_keys), (
                name, gatk_rows, native_rows)
            if gatk_rows != native_rows:
                diff = "\n".join(difflib.unified_diff(
                    gatk_rows, native_rows, fromfile="gatk", tofile="native", lineterm=""))
                raise AssertionError(f"{name} VCF row mismatch:\n{diff}")
            observed_keys = tuple(
                tuple(row.split("\t")[index] for index in (0, 1, 3, 4))
                for row in native_rows
            )
            assert observed_keys == expected_keys, (name, observed_keys)
            if name == PAIRHMM_HAPLOTYPE_ORACLE:
                gatk_haplotypes = gatk_pairhmm_haplotypes(gatk_pairhmm)
                native_haplotypes = native_shared_pairhmm_haplotypes(native_run.stderr)
                assert gatk_haplotypes == native_haplotypes, {
                    "missing": sorted(gatk_haplotypes - native_haplotypes),
                    "extra": sorted(native_haplotypes - gatk_haplotypes),
                }

    print(json.dumps({
        "status": "pass",
        "fixture": "hcc1143_tumor_normal_real_data",
        "oracles": [name for name, _, _ in ORACLES],
        "record_count": sum(len(keys) for _, _, keys in ORACLES),
        "vcf_header_exact": True,
        "vcf_schema_header_exact_excluding_execution_provenance": True,
        "vcf_row_exact": True,
        "pairhmm_haplotype_set_exact": True,
        "pairhmm_haplotype_oracle": PAIRHMM_HAPLOTYPE_ORACLE,
        "sample_order": [TUMOR_SAMPLE, NORMAL_SAMPLE],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
