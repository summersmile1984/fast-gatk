#!/usr/bin/env python3
"""Pin Mutect2's GATK AssemblyRegion input halo at an explicit -L edge.

The VCF core is chrM:8860-8890, but GATK's AssemblyRegionWalker fetches the
100-base assembly padding on both sides before constructing its EventMap.
This test deliberately audits that Host read-input contract separately from the
still-in-progress high-depth mitochondrial graph oracle.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "real" / "mitomode" / "mito_shifted_8000.fasta"
BAM = ROOT / "testdata" / "real" / "mitomode" / "mito.bam"
CORE = "chrM:8860-8890"


def data_positions(path: Path) -> list[int]:
    return [
        int(line.split("\t", 2)[1])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, BAM, Path(f"{BAM}.bai"))
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise AssertionError({"missing": missing})
        print(json.dumps({"status": "skip", "reason": "mitochondrial GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-mito-halo-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.vcf.gz"
        gatk_run = subprocess.run([
            str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(BAM), "-L", CORE,
            "--mitochondria-mode", "--max-reads-per-alignment-start", "75",
            "--native-pair-hmm-threads", "1", "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false", "-O", str(gatk_output),
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        processed = re.search(r"(\d+) reads processed", gatk_run.stderr)
        assert processed is not None, gatk_run.stderr
        gatk_reads = int(processed.group(1))
        # This count comes from the GATK command above, but pin it as a
        # release-fixture invariant as well: it makes an accidental crop back
        # to the 3,037 core-only reads immediately diagnosable.
        assert gatk_reads == 6506, gatk_reads

        native_output = work / "native.vcf"
        native_run = subprocess.run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM), "-L", CORE,
            "--mitochondria-mode", "--max-reads-per-alignment-start", "75",
            "--threads", "1", "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false", "-O", str(native_output),
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        native_summary = json.loads(native_run.stdout.splitlines()[-1])
        assert native_summary["tumor_reads"] == gatk_reads, native_summary
        positions = data_positions(native_output)
        assert all(8860 <= position <= 8890 for position in positions), positions

    print(json.dumps({
        "status": "pass",
        "fixture": "mitochondrial_interval_halo",
        "core": CORE,
        "gatk_reads": gatk_reads,
        "native_reads": native_summary["tumor_reads"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
