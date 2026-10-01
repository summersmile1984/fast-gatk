#!/usr/bin/env python3
"""P2 A1: HaplotypeCaller ReadThreadingAssembler graph oracle vs GATK 4.6.2.0.

Compares, on the pinned chr17 69k-70k fixture:
  * VCF data rows (POS/REF/ALT and GT/AD/DP/GQ/PL/QUAL)
  * GATK --debug-assembly "Using kmer size of N" vs native graph_kmer_sizes_used
  * CLI request list graph_kmer_sizes == [10, 25]
  * native --graph-output DOT haplotypes and GATK --graph-output vertex labels
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def variant_records(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def sample(record: list[str]) -> dict[str, str]:
    return dict(zip(record[8].split(":"), record[9].split(":")))


def gatk_used_kmers(log: str) -> list[int]:
    found = {int(value) for value in re.findall(
        r"Using kmer size of (\d+) in read threading assembler", log)}
    return sorted(found)


def gatk_dot_sequences(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8")
    sequences = set()
    for match in re.finditer(r'\[label="([ACGT]+)', text):
        sequences.add(match.group(1))
    for match in re.finditer(r'_seq_([ACGT]+)', text):
        sequences.add(match.group(1))
    return sequences


def native_dot_sequences(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8")
    return set(re.findall(r'\[label="([ACGT]+)"', text))


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout[-2000:]}\nstderr={result.stderr[-4000:]}"
        )
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.exists() for path in (native, java, gatk, bam, reference)):
        oracle_guard.oracle_not_verified("verify_hc_graph_gatk_oracle.py", java, gatk)
        raise SystemExit("missing pinned HC graph oracle assets")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-graph-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        gatk_dot = work / "gatk.dot"
        native_vcf = work / "native.vcf"
        native_dot = work / "native.dot"
        native_manifest = work / "native.json"
        gatk_run = run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
            "--debug-assembly", "true",
            "--graph-output", str(gatk_dot),
            "-O", str(gatk_vcf),
        ])
        run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
            "--threads", "1", "--add-output-vcf-command-line", "false",
            "--graph-output", str(native_dot),
            "-O", str(native_vcf), "--output-manifest", str(native_manifest),
        ])

        gatk_rows = variant_records(gatk_vcf)
        native_rows = variant_records(native_vcf)
        gatk_keys = [(row[0], int(row[1]), row[3], row[4]) for row in gatk_rows]
        native_keys = [(row[0], int(row[1]), row[3], row[4]) for row in native_rows]
        if gatk_keys != native_keys:
            raise AssertionError(f"graph-oracle VCF keys differ: {gatk_keys} vs {native_keys}")
        for gatk_row, native_row in zip(gatk_rows, native_rows):
            gatk_sample = sample(gatk_row)
            native_sample = sample(native_row)
            for field in ("GT", "AD", "DP", "GQ", "PL"):
                if gatk_sample.get(field) != native_sample.get(field):
                    raise AssertionError(
                        f"{gatk_row[0]}:{gatk_row[1]} {field} "
                        f"{gatk_sample.get(field)} != {native_sample.get(field)}")
            if gatk_row[5] != native_row[5]:
                raise AssertionError(
                    f"{gatk_row[0]}:{gatk_row[1]} QUAL {gatk_row[5]} != {native_row[5]}")

        telemetry = json.loads(native_manifest.read_text(encoding="utf-8"))["telemetry"]
        if telemetry.get("graph_kmer_sizes") != [10, 25]:
            raise AssertionError(f"CLI k-mer list {telemetry.get('graph_kmer_sizes')}")
        if not telemetry.get("graph_used"):
            raise AssertionError("native graph_used is false")
        if int(telemetry.get("graph_haplotype_paths") or 0) <= 0:
            raise AssertionError("native produced no graph haplotypes")
        used = [int(value) for value in telemetry.get("graph_kmer_sizes_used") or []]
        expected_used = gatk_used_kmers(gatk_run.stderr + gatk_run.stdout)
        if used != expected_used:
            raise AssertionError(
                f"accepted k-mer sizes {used} != GATK Using-kmer set {expected_used}")
        if int(telemetry.get("graph_kmer_size_selected") or 0) != max(expected_used):
            raise AssertionError(
                f"selected k {telemetry.get('graph_kmer_size_selected')} "
                f"!= max GATK used {max(expected_used)}")

        gatk_seqs = gatk_dot_sequences(gatk_dot)
        native_seqs = native_dot_sequences(native_dot)
        if not native_dot.is_file() or native_dot.stat().st_size == 0:
            raise AssertionError("native --graph-output was not written")
        if not gatk_seqs:
            raise AssertionError("GATK --graph-output contained no ACGT vertex labels")
        # SeqGraph vertices are fragments of assembled haplotypes.
        missing = [seq for seq in gatk_seqs if not any(seq in hap for hap in native_seqs)]
        if missing:
            raise AssertionError(
                f"GATK SeqGraph sequences absent from native haplotypes: {missing[:3]}")

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": "17:69000-70000",
            "vcf_keys": [f"{c}:{p}:{r}>{a}" for c, p, r, a in gatk_keys],
            "gatk_kmer_sizes_used": expected_used,
            "native_kmer_sizes_used": used,
            "native_kmer_size_selected": telemetry.get("graph_kmer_size_selected"),
            "native_graph_haplotype_paths": telemetry.get("graph_haplotype_paths"),
            "gatk_graph_sequences": len(gatk_seqs),
            "native_graph_sequences": len(native_seqs),
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
