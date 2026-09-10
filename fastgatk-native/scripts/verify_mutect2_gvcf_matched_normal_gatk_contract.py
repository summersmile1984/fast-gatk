#!/usr/bin/env python3
"""Pin GATK's unsupported Mutect2 gVCF matched-normal combination.

GATK 4.6.2 constructs tumor and normal likelihood matrices, then delegates
reference confidence to a model that only accepts one sample.  It therefore
rejects both ``Mutect2 -ERC GVCF`` and ``BP_RESOLUTION`` with a matched normal
instead of producing a mixed-sample reference-confidence file.  The native
C++ Host must fail before launching its Kokkos PairHMM path as well; accepting
this command would be less compatible than GATK even if the result appeared
plausible.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / \
    "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "real" / "cnv_somatic" / "human_g1k_v37.chr-20.truncated.fasta"
TUMOR = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_tumor.bam"
NORMAL = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_normal.bam"
TUMOR_SAMPLE = "HCC1143"
NORMAL_SAMPLE = "HCC1143 BL"
INTERVAL = "20:67000-69000"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError({"missing": missing})
        print(json.dumps({"status": "skip", "reason": "oracle inputs unavailable"}))
        return 0

    modes = ("GVCF", "BP_RESOLUTION")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-gvcf-normal-") as directory:
        work = Path(directory)
        for mode in modes:
            gatk_output = work / f"gatk-{mode}.g.vcf.gz"
            native_output = work / f"native-{mode}.g.vcf.gz"
            gatk = run([
                str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
                "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
                "-tumor", TUMOR_SAMPLE, "-normal", NORMAL_SAMPLE, "-L", INTERVAL,
                "-ERC", mode, "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false", "-O", str(gatk_output),
            ])
            assert gatk.returncode != 0, (mode, gatk.stderr[-4000:])
            assert "readLikelihoods must contain exactly one sample" in gatk.stderr, (
                mode, gatk.stderr[-4000:])

            native = run([
                str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
                "--tumor-sample", TUMOR_SAMPLE, "--normal-sample", NORMAL_SAMPLE,
                "-L", INTERVAL, "-ERC", mode, "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false", "-O", str(native_output),
            ])
            assert native.returncode != 0, (mode, native.stdout)
            assert "reference confidence requires exactly one sample" in native.stderr, (
                mode, native.stderr)
            assert not native_output.exists(), (mode, native_output)

    print(json.dumps({
        "status": "pass",
        "modes": list(modes),
        "oracle": "GATK-4.6.2.0",
        "rejection": "exactly_one_sample_reference_confidence",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
