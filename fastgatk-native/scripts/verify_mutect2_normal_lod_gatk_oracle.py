#!/usr/bin/env python3
"""Pinned GATK oracle for the Mutect2 ``--normal-lod`` emission boundary.

The release-pinned SomaticGenotypingEngine suppresses a candidate when the
matched normal ALT has insufficient normal log-odds.  This verifier builds a
small deterministic tumor/normal BAM pair with a heterozygous normal ALT and
checks the stable boundary (a very low threshold emits, zero suppresses) for
both Java GATK 4.6.2.0 and the native replacement.  It deliberately compares
site presence rather than release-specific posterior text.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def fasta_contig(path: Path, name: str) -> str:
    current = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            continue
        if current == name:
            sequence.append(line.strip())
    result = "".join(sequence)
    if not result:
        raise AssertionError(f"reference contig not found: {name}")
    return result


def vcf_sites(path: Path) -> list[tuple[str, int, str, str]]:
    result: list[tuple[str, int, str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            result.append((fields[0], int(fields[1]), fields[3], fields[4]))
    return result


def write_sam(path: Path, reference: str, sample: str, alternate: str, alternate_reads: int) -> None:
    position = 69005
    start = 68960
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:17\tLN:{len(reference)}",
        f"@RG\tID:RG\tSM:{sample}\tPL:ILLUMINA\tLB:LIB\tPU:UNIT",
    ]
    for read_index in range(20):
        read_start = start + read_index
        sequence = list(reference[read_start - 1:read_start - 1 + 50])
        if read_index < alternate_reads:
            sequence[position - read_start] = alternate
        lines.append(
            f"read{read_index}\t0\t17\t{read_start}\t60\t50M\t*\t0\t0\t"
            f"{''.join(sequence)}\t{'I' * 50}\tRG:Z:RG"
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def sort_indexed(java: Path, gatk: Path, sam: Path, bam: Path) -> None:
    result = run([
        str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
        "-I", str(sam), "-O", str(bam), "-SO", "coordinate",
        "--CREATE_INDEX", "true", "--VALIDATION_STRINGENCY", "STRICT",
    ])
    if result.returncode != 0:
        raise RuntimeError(f"GATK SortSam failed:\n{result.stderr[-4000:]}")
    # Picard writes the index as <output basename>.bai rather than .bam.bai.
    if not bam.exists() or not bam.with_suffix(".bai").exists():
        raise AssertionError(f"SortSam did not create indexed BAM: {bam}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_mutect2_normal_lod_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 normal-lod oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    reference_sequence = fasta_contig(reference, "17")
    # A normal with 3/20 ALT reads is enough for Java's normal evidence to
    # suppress the site at the default/zero threshold but not at -100.
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-normal-lod-") as directory:
        work = Path(directory)
        tumor_sam, normal_sam = work / "tumor.sam", work / "normal.sam"
        tumor_bam, normal_bam = work / "tumor.bam", work / "normal.bam"
        write_sam(tumor_sam, reference_sequence, "TUMOR", "C", 10)
        write_sam(normal_sam, reference_sequence, "NORMAL", "C", 3)
        sort_indexed(java, gatk, tumor_sam, tumor_bam)
        sort_indexed(java, gatk, normal_sam, normal_bam)

        thresholds = (-100.0, 0.0)
        java_sites: dict[float, list[tuple[str, int, str, str]]] = {}
        native_sites: dict[float, list[tuple[str, int, str, str]]] = {}
        native_manifests: dict[float, dict[str, object]] = {}
        for threshold in thresholds:
            label = "neg100" if threshold < 0 else "zero"
            java_output = work / f"gatk-{label}.vcf.gz"
            native_output = work / f"native-{label}.vcf.gz"
            manifest = work / f"native-{label}.manifest.json"
            java_run = run([
                str(java), "-Xmx1g", "-jar", str(gatk), "Mutect2",
                "-R", str(reference), "-I", str(tumor_bam), "-I", str(normal_bam),
                "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL",
                "-L", "17:69000-69020", "--normal-lod", str(threshold),
                "-O", str(java_output), "--create-output-variant-index", "false",
            ])
            if java_run.returncode != 0:
                raise RuntimeError(f"GATK Mutect2 normal-lod={threshold} failed:\n{java_run.stderr[-5000:]}")
            native_run = run([
                str(binary), "-R", str(reference), "-I", str(tumor_bam),
                "--normal-input", str(normal_bam), "--tumor-sample", "TUMOR",
                "--normal-sample", "NORMAL", "-L", "17:69000-69020",
                "--normal-lod", str(threshold), "-O", str(native_output),
                "--create-output-variant-index", "false", "--min-depth", "1",
                "--min-alt-support", "1", "--output-manifest", str(manifest),
            ])
            if native_run.returncode != 0:
                raise RuntimeError(f"native Mutect2 normal-lod={threshold} failed:\n{native_run.stderr[-5000:]}")
            java_sites[threshold] = vcf_sites(java_output)
            native_sites[threshold] = vcf_sites(native_output)
            native_manifests[threshold] = json.loads(manifest.read_text(encoding="utf-8"))

        assert java_sites[-100.0], java_sites
        assert java_sites[0.0] == [], java_sites
        assert native_sites[-100.0], native_sites
        assert native_sites[0.0] == [], native_sites
        assert native_sites[-100.0] == java_sites[-100.0], {
            "java": java_sites[-100.0], "native": native_sites[-100.0]
        }
        assert native_sites[0.0] == java_sites[0.0]
        for threshold, manifest in native_manifests.items():
            compatibility = manifest["compatibility"]
            assert compatibility["normal_lod"] == threshold, compatibility
            assert manifest["telemetry"]["normal_lod"] == threshold
        print(json.dumps({
            "status": "pass",
            "gatk_neg100_sites": java_sites[-100.0],
            "gatk_zero_sites": java_sites[0.0],
            "native_neg100_sites": native_sites[-100.0],
            "native_zero_sites": native_sites[0.0],
            "normal_lod_boundary_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
