#!/usr/bin/env python3
"""Pinned GATK oracle for soft-clip reversion at a contig start.

GATK's issue-3466 fixture contains a paired 28S222M read starting at
zero-based reference position 17.  ``ReadClipper.revertSoftClippedBases``
must hard-clip its eleven out-of-contig prefix before PairHMM qualification;
otherwise the remaining match segment is shifted and the MT:152 call loses
one supporting read.  This is deliberately a Host CIGAR-projection contract.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout={completed.stdout[-2000:]}\nstderr={completed.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    fixture = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/"
        "tools/haplotypecaller/issue3466_gatk_cigar_error")
    reference = fixture / "GRCh37_MTonly.fa"
    bam = fixture / "culprit.bam"
    if not all(path.is_file() for path in (native, java, gatk, reference, bam)):
        oracle_guard.oracle_not_verified('verify_hc_softclip_contig_start_gatk_oracle.py', java, gatk)
        raise SystemExit("missing HC soft-clip contig-start oracle assets")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-softclip-contig-start-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.json"
        run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_vcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ])
        run([
            str(native), "-R", str(reference), "-I", str(bam), "-O", str(native_vcf),
            "--threads", "1", "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--output-manifest", str(manifest),
        ])
        if gatk_vcf.read_bytes() != native_vcf.read_bytes():
            raise AssertionError("native VCF differs from GATK for issue3466 fixture")
        compatibility = json.loads(manifest.read_text(encoding="utf-8"))["compatibility"]
        if compatibility.get("pairhmm_reads_disqualified") != 0:
            raise AssertionError("contig-start soft-clipped read was incorrectly disqualified")
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "fixture": "HaplotypeCaller issue3466",
            "vcf_byte_identical": True,
            "pairhmm_reads_disqualified": 0,
            "scope": "Host soft-clip reversion and hard-clip boundary",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
