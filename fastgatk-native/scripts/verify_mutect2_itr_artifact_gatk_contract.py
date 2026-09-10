#!/usr/bin/env python3
"""Lock Mutect2's default ITR palindrome clipping to GATK source semantics.

This is the CLI-shaped counterpart of GATK's
``PalindromeArtifactClipReadTransformerUnitTest``.  The synthetic proper pair
has a 10S terminal artifact whose first five aligned bases complete the
minimum five-base palindrome.  Two otherwise identical controls exercise the
0.9 matching-fraction gate: a shifted adaptor boundary and an ordinary
soft-clipped read.  The native sidecar is deliberately used as the Host-boundary
observation point; graph and PairHMM execution remain ordinary Kokkos runs.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_fixture(reference_path: Path, sam_path: Path) -> None:
    # Keep all coordinates one-based here, as in the GATK transformer.  Its
    # forward-read adaptor boundary is readStart + TLEN = 201, and it queries
    # reference [186, 200] against the first 15 read bases in reverse order.
    reference = list("A" * 512)
    itr_reference_start = 186
    itr_reference = "ACGTTAGCGTACGCA"  # 15 bp, intentionally non-periodic
    reference[itr_reference_start - 1:itr_reference_start - 1 + len(itr_reference)] = itr_reference
    artifact_prefix = reverse_complement(itr_reference)
    read_start = 101
    soft_clip = 10
    aligned_bases = 30
    # GATK's read index 10..14 is both aligned and part of the palindrome.
    # Give it the corresponding reference sequence so the CIGAR remains a
    # valid 10S30M alignment while all 15 comparisons match.
    reference[read_start - 1:read_start - 1 + 5] = artifact_prefix[soft_clip:]
    reference_sequence = "".join(reference)
    reference_path.write_text(">1\n" + reference_sequence + "\n", encoding="ascii")

    aligned_suffix = reference_sequence[read_start - 1 + 5:read_start - 1 + aligned_bases]
    artifact_sequence = artifact_prefix + aligned_suffix
    ordinary_sequence = "A" * soft_clip + reference_sequence[read_start - 1:read_start - 1 + aligned_bases]
    # The shifted boundary compares [187, 201] rather than [186, 200].  It
    # must stay below GATK's >=0.9 (14/15 is insufficient in this fixture).
    shifted_expected = reverse_complement(reference_sequence[186:201])
    assert sum(a == b for a, b in zip(artifact_sequence[:15], shifted_expected)) < 14

    def record(name: str, sequence: str, tlen: int) -> str:
        mate_start = read_start + tlen - len(sequence)
        return (
            f"{name}\t99\t1\t{read_start}\t60\t10S30M\t=\t{mate_start}\t{tlen}\t"
            f"{sequence}\t{'I' * len(sequence)}\tRG:Z:RG"
        )

    sam_path.write_text("\n".join((
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:1\tLN:512",
        "@RG\tID:RG\tSM:TUMOR\tPL:ILLUMINA\tLB:LIB\tPU:UNIT",
        record("itr-artifact", artifact_sequence, 100),
        record("shifted-adaptor", artifact_sequence, 101),
        record("ordinary-softclip", ordinary_sequence, 100),
    )) + "\n", encoding="ascii")


def sort_indexed(binary: Path, sam: Path, bam: Path) -> None:
    sorter = binary.with_name("fastgatk-sort-sam")
    if not sorter.is_file():
        raise AssertionError(f"native SortSam sibling is required: {sorter}")
    result = run([
        str(sorter), "-I", str(sam), "-O", str(bam), "--sort-order", "coordinate",
        "--CREATE_INDEX", "true",
    ])
    if result.returncode != 0:
        raise RuntimeError(f"native SortSam failed:\n{result.stderr[-5000:]}")
    if not bam.is_file() or not Path(f"{bam}.bai").is_file():
        raise AssertionError("SortSam did not produce an indexed BAM")


def invoke(binary: Path, reference: Path, reads: Path, output: Path, stats: Path,
           ignore_itr: bool, streamed: bool = False) -> dict[str, object]:
    command = [
        str(binary), "-R", str(reference), "-I", str(reads), "--tumor-sample", "TUMOR",
        "-L", "1:95-140", "-O", str(output), "--stats", str(stats),
        "--min-depth", "1", "--min-alt-support", "1",
        # GATK filters the original 40 bp read, then clips it to 30 bp.  A
        # native implementation that clips before filtering would reject the
        # artifact here and report one filtered read instead of zero.
        "--min-read-length", "35",
        "--create-output-variant-index", "false",
    ]
    if streamed:
        command += ["--stream-by-region", "100"]
    if ignore_itr:
        command.append("--ignore-itr-artifacts")
    result = run(command)
    if result.returncode != 0:
        raise RuntimeError(f"native Mutect2 failed:\n{result.stderr[-5000:]}")
    return json.loads(stats.read_text(encoding="utf-8"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    if not binary.is_file():
        raise SystemExit(f"Mutect2 binary is required: {binary}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-itr-artifact-") as directory:
        work = Path(directory)
        reference, sam = work / "reference.fa", work / "reads.sam"
        write_fixture(reference, sam)
        default_stats = invoke(binary, reference, sam, work / "default.vcf",
                               work / "default.stats.json", False)
        ignored_stats = invoke(binary, reference, sam, work / "ignored.vcf",
                               work / "ignored.stats.json", True)
        assert default_stats["itr_artifact_clipping_enabled"] is True, default_stats
        assert default_stats["itr_artifact_clipped_reads"] == 1, default_stats
        assert default_stats["read_filter_filtered_reads"] == 0, default_stats
        assert ignored_stats["itr_artifact_clipping_enabled"] is False, ignored_stats
        assert ignored_stats["itr_artifact_clipped_reads"] == 0, ignored_stats
        bam = work / "reads.bam"
        sort_indexed(binary, sam, bam)
        stream_default = invoke(binary, reference, bam, work / "stream-default.vcf",
                                work / "stream-default.stats.json", False, True)
        stream_ignored = invoke(binary, reference, bam, work / "stream-ignored.vcf",
                                work / "stream-ignored.stats.json", True, True)
        assert stream_default["stream_by_region"] is True, stream_default
        assert stream_default["itr_artifact_clipped_reads"] == 1, stream_default
        assert stream_ignored["stream_by_region"] is True, stream_ignored
        assert stream_ignored["itr_artifact_clipped_reads"] == 0, stream_ignored
        print(json.dumps({
            "status": "pass",
            "gatk_min_palindrome_size": 5,
            "gatk_matching_fraction": 0.9,
            "default_itr_artifact_clipped_reads": default_stats["itr_artifact_clipped_reads"],
            "ignore_itr_artifacts_clipped_reads": ignored_stats["itr_artifact_clipped_reads"],
            "stream_default_itr_artifact_clipped_reads": stream_default["itr_artifact_clipped_reads"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
