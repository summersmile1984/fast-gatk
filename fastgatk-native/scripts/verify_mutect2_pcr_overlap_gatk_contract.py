#!/usr/bin/env python3
"""Pin Mutect2's PCR overlap-quality arguments to GATK's Host contract.

``Mutect2Engine.callRegion`` calls
``cleanOverlappingReadPairs(..., pcrSnvQual / 2, pcrIndelQual / 2)``.  This
small proper pair contains 20 agreeing overlap bases with Q10 base qualities
and Q10 BI/BD tags: the release defaults (Q40 -> half Q20) leave them alone,
while Q8/Q14 caps respectively change the bases and both indel arrays.  The
sidecar counters observe the Host transformation immediately before the
ordinary graph and Kokkos PairHMM stages.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def write_fixture(reference: Path, reads: Path) -> None:
    sequence = "A" * 256
    reference.write_text(">1\n" + sequence + "\n", encoding="ascii")
    bases = "A" * 60
    # Base qualities and both indel-quality tags are Q10.  The pair spans
    # 1..60 and 41..100, so exactly 20 aligned positions overlap.
    qualities = "+" * 60
    reads.write_text("\n".join((
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:1\tLN:256",
        "@RG\tID:rg\tSM:TUMOR\tPL:ILLUMINA",
        f"pair\t99\t1\t1\t60\t60M\t=\t41\t100\t{bases}\t{qualities}"
        f"\tRG:Z:rg\tBI:Z:{qualities}\tBD:Z:{qualities}",
        f"pair\t147\t1\t41\t60\t60M\t=\t1\t-100\t{bases}\t{qualities}"
        f"\tRG:Z:rg\tBI:Z:{qualities}\tBD:Z:{qualities}",
    )) + "\n", encoding="ascii")


def invoke(binary: Path, reference: Path, reads: Path, output: Path, stats: Path,
           *extra: str) -> dict[str, object]:
    result = subprocess.run([
        str(binary), "-R", str(reference), "-I", str(reads),
        "--tumor-sample", "TUMOR", "-L", "1:1-120", "-O", str(output),
        "--stats", str(stats), "--min-depth", "100", "--min-alt-support", "1",
        "--create-output-variant-index", "false", *extra,
    ], text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"native Mutect2 failed:\n{result.stderr[-5000:]}")
    return json.loads(stats.read_text(encoding="utf-8"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", root / "fastgatk-native/build/fastgatk-mutect2"))
    if not binary.is_file():
        raise SystemExit(f"Mutect2 binary is required: {binary}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-pcr-overlap-") as directory:
        work = Path(directory)
        reference, reads = work / "reference.fa", work / "reads.sam"
        write_fixture(reference, reads)
        defaults = invoke(binary, reference, reads, work / "default.vcf",
                          work / "default.stats.json")
        configured = invoke(binary, reference, reads, work / "configured.vcf",
                            work / "configured.stats.json",
                            "--pcr-snv-qual", "8", "--pcr-indel-qual", "14")
        for stats in (defaults, configured):
            assert stats["overlapping_quality_correction_used"] is True, stats
            assert stats["overlapping_quality_correction_metadata_available"] is True, stats
            assert stats["overlapping_pairs"] == 1, stats
            assert stats["overlapping_bases"] == 20, stats
            assert stats["overlapping_conflicting_bases"] == 0, stats
        # Q10 does not change under the default half-Q20 base/indel caps.
        assert defaults["overlapping_pcr_snv_quality"] == 40, defaults
        assert defaults["overlapping_pcr_indel_quality"] == 40, defaults
        assert defaults["overlapping_quality_caps"] == 0, defaults
        assert defaults["overlapping_indel_quality_caps"] == 0, defaults
        # GATK passes 8 / 2 and 14 / 2, so every one of the 20 agreeing
        # bases in both reads is capped to Q4, while BI and BD each change
        # at every overlapping position in both reads (20 * 2 * 2 = 80).
        assert configured["overlapping_pcr_snv_quality"] == 8, configured
        assert configured["overlapping_pcr_indel_quality"] == 14, configured
        assert configured["overlapping_quality_caps"] == 40, configured
        assert configured["overlapping_indel_quality_caps"] == 80, configured
        print(json.dumps({
            "status": "pass",
            "gatk_pcr_snv_half_quality": 4,
            "gatk_pcr_indel_half_quality": 7,
            "base_quality_caps": configured["overlapping_quality_caps"],
            "indel_quality_caps": configured["overlapping_indel_quality_caps"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
