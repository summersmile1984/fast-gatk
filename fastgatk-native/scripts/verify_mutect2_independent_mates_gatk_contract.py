#!/usr/bin/env python3
"""Verify Mutect2 --independent-mates changes only fragment identity.

GATK's SomaticGenotypingEngine calls ``groupEvidence`` with either
``GATKRead::getName`` (default) or ``read -> read`` (--independent-mates).
Use the pinned NA12878 Mutect2 integration corpus so the assertion covers the
actual assembled haplotype and Kokkos PairHMM path, rather than a pileup-only
synthetic shortcut.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def invoke(binary: Path, reference: Path, reads: Path, output: Path, stats: Path,
           independent: bool) -> dict[str, object]:
    command = [
        str(binary), "-R", str(reference), "-I", str(reads),
        "--tumor-sample", "NA12878", "-L", "17:69000-70000", "-O", str(output),
        "--stats", str(stats), "--min-depth", "1", "--min-alt-support", "1",
        "--kmer-size", "7", "--min-kmer-count", "1",
        "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64",
        "--max-haplotype-combination-alleles", "5", "--error-correct-reads",
        "--kmer-length-for-read-error-correction", "5",
        "--min-observations-for-kmer-to-be-solid", "2",
        "--error-correction-log-odds", "3.5",
        "--create-output-variant-index", "false",
    ]
    if independent:
        command.append("--independent-mates")
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"native Mutect2 failed:\n{result.stderr[-5000:]}")
    return json.loads(stats.read_text(encoding="utf-8"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", root / "fastgatk-native/build/fastgatk-mutect2"))
    reads = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not binary.is_file() or not reads.is_file() or not reference.is_file():
        raise SystemExit("Mutect2 binary and pinned NA12878 corpus are required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-independent-mates-") as directory:
        work = Path(directory)
        grouped = invoke(binary, reference, reads, work / "grouped.vcf",
                         work / "grouped.stats.json", False)
        independent = invoke(binary, reference, reads, work / "independent.vcf",
                             work / "independent.stats.json", True)
        assert grouped["independent_mates"] is False, grouped
        assert independent["independent_mates"] is True, independent
        assert grouped["somatic_evidence_groups"] > 0, grouped
        assert independent["somatic_evidence_groups"] > grouped["somatic_evidence_groups"], (
            grouped, independent)
        assert grouped["pairhmm_pairs"] > 0, grouped
        assert grouped["pairhmm_pairs"] == independent["pairhmm_pairs"], (
            grouped, independent)
        print(json.dumps({
            "status": "pass",
            "default_fragment_groups": grouped["somatic_evidence_groups"],
            "independent_read_groups": independent["somatic_evidence_groups"],
            "pairhmm_pairs": independent["pairhmm_pairs"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
