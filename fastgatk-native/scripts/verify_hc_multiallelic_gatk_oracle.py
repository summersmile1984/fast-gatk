#!/usr/bin/env python3
"""Probe/strict oracle for a real HC multi-ALT assembly site."""
from __future__ import annotations

import json
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get("FASTGATK_HC_BINARY",
                               str(root / "fastgatk-native/build/fastgatk-hc-call")))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    if not native.exists() or not jar.exists() or not java.exists():
        raise SystemExit("missing native HC or GATK/JDK oracle")
    header = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n"
              "@RG\tID:1\tSM:S1\tPL:ILLUMINA\n")
    reads: list[str] = []
    reference_sequence = ("ACGT" * 10)[:40]
    reference_base = reference_sequence[10]
    for index, alternate in enumerate(("A", "T")):
        for copy in range(20):
            sequence = reference_sequence[:10] + alternate + reference_sequence[11:21]
            reads.append(
                f"{alternate.lower()}{index}{copy}\t0\tchr1\t1\t60\t21M\t*\t0\t0\t"
                f"{sequence}\t{'I' * 21}\tRG:Z:1\tNM:i:1\tMD:Z:10G10\n")
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-multiallelic-oracle-") as directory:
        work = Path(directory)
        sam = work / "input.sam"
        bam = work / "input.bam"
        reference = work / "reference.fa"
        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf"
        manifest = work / "native.manifest.json"
        sam.write_text(header + "".join(reads), encoding="utf-8")
        reference_text = ">chr1\n" + reference_sequence + "\n"
        reference.write_text(reference_text, encoding="utf-8")
        # GATK requires both the FASTA index and sequence dictionary even for
        # this tiny synthetic reference.  The sequence occupies one complete
        # line: FASTA offset=6, bases/line=40, bytes/line=41.
        reference.with_suffix(reference.suffix + ".fai").write_text(
            "chr1\t40\t6\t40\t41\n", encoding="utf-8")
        md5 = hashlib.md5(reference_sequence.encode("ascii")).hexdigest()
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\tM5:" + md5 + "\n",
            encoding="utf-8")
        result = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "SortSam",
            "-I", str(sam), "-O", str(bam), "-SO", "coordinate",
            "--CREATE_INDEX", "true", "--VALIDATION_STRINGENCY", "STRICT",
        ], text=True, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(f"GATK SortSam failed:\n{result.stderr[-4000:]}")
        result = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", "chr1:1-40",
            "-O", str(gatk_output), "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(f"GATK HaplotypeCaller failed:\n{result.stderr[-4000:]}")
        result = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", "chr1:1-40",
            "-O", str(native_output), "--threads", "2", "--min-depth", "1",
            "--min-alt-support", "1", "--output-manifest", str(manifest),
        ], text=True, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(f"native HaplotypeCaller failed:\n{result.stderr[-4000:]}")
        print(json.dumps({
            "gatk_rows": rows(gatk_output),
            "native_rows": rows(native_output),
            "native_telemetry": json.loads(manifest.read_text(encoding="utf-8"))["telemetry"],
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
