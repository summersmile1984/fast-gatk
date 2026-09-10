#!/usr/bin/env python3
"""Pin GATK's single-sample requirement for Mutect2 reference confidence."""

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
REFERENCE = ROOT / "testdata" / "downloads" / "reference" / "hs37d5.fa.gz"
TUMOR = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "tumor.bam"
SECOND_TUMOR = ROOT / "testdata" / "real" / "dream_synthetic" / "chr20" / "tumor_1.bam"
TUMOR_SAMPLE = "synthetic.challenge.set1.tumor"
INTERVAL = "20:10022000-10023500"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, SECOND_TUMOR)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError({"missing": missing})
        print(json.dumps({"status": "skip", "reason": "oracle inputs unavailable"}))
        return 0

    modes = ("GVCF", "BP_RESOLUTION")
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-gvcf-multitumor-") as directory:
        work = Path(directory)
        for mode in modes:
            gatk_output = work / f"gatk-{mode}.g.vcf.gz"
            native_output = work / f"native-{mode}.g.vcf.gz"
            gatk = run([
                str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
                "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(SECOND_TUMOR),
                "-tumor", TUMOR_SAMPLE, "-L", INTERVAL, "-ERC", mode,
                "--create-output-variant-index", "false", "--add-output-vcf-command-line", "false",
                "-O", str(gatk_output),
            ])
            assert gatk.returncode != 0, (mode, gatk.stderr[-4000:])
            assert "readLikelihoods must contain exactly one sample" in gatk.stderr, (
                mode, gatk.stderr[-4000:])

            native = run([
                str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(SECOND_TUMOR),
                "--tumor-sample", TUMOR_SAMPLE, "-L", INTERVAL, "-ERC", mode,
                "--create-output-variant-index", "false", "--add-output-vcf-command-line", "false",
                "-O", str(native_output),
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
