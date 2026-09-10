#!/usr/bin/env python3
"""Pinned GATK oracle for Mutect2 --force-active.

GATK keeps ActivityProfile segment boundaries unchanged when force-active is
enabled, but marks every segment active so the AssemblyRegion walker processes
the complete profile.  This verifier checks that bounded contract for both
the IGV profile and the VCF writer, rather than claiming somatic posterior
bit identity.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def read_igv(path: Path) -> list[tuple[str, int, int, str, str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("Chromosome\t"):
            continue
        fields = line.split("\t")
        assert len(fields) == 5, fields
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]))
    return rows


def read_vcf_records(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def vcf_site_keys(records: list[str]) -> list[tuple[str, str, str, str]]:
    return [tuple(record.split("\t", 5)[:5]) for record in records]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 force-active oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-force-active-") as directory:
        work = Path(directory)
        java_igv: dict[bool, Path] = {}
        java_vcf: dict[bool, Path] = {}
        native_igv: dict[bool, Path] = {}
        native_vcf: dict[bool, Path] = {}
        native_manifest: dict[bool, Path] = {}
        for force in (False, True):
            suffix = "true" if force else "false"
            java_vcf[force] = work / f"gatk-{suffix}.vcf.gz"
            java_igv[force] = work / f"gatk-{suffix}.igv"
            native_vcf[force] = work / f"native-{suffix}.vcf.gz"
            native_igv[force] = work / f"native-{suffix}.igv"
            native_manifest[force] = work / f"native-{suffix}.manifest.json"
            java_run = run([
                str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
                "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
                "-O", str(java_vcf[force]), "--assembly-region-out", str(java_igv[force]),
                "--force-active", str(force).lower(),
                "--create-output-variant-index", "false",
            ])
            assert java_run.returncode == 0, java_run.stderr
            native_run = run([
                str(binary), "-R", str(reference), "-I", str(bam),
                "--tumor-sample", "NA12878", "-L", "17:69000-70000",
                "-O", str(native_vcf[force]), "--assembly-region-out", str(native_igv[force]),
                "--force-active=" + str(force).lower(),
                "--create-output-variant-index", "false",
                "--min-depth", "1", "--min-alt-support", "1",
                "--output-manifest", str(native_manifest[force]),
            ])
            assert native_run.returncode == 0, native_run.stderr

        gatk_false = read_igv(java_igv[False])
        gatk_true = read_igv(java_igv[True])
        native_false = read_igv(native_igv[False])
        native_true = read_igv(native_igv[True])
        assert gatk_false, "GATK force-active fixture produced no profile rows"
        assert gatk_false == native_false, {
            "expected": gatk_false,
            "actual": native_false,
        }
        assert gatk_true == native_true, {
            "expected": gatk_true,
            "actual": native_true,
        }
        # The option changes only the effective active state.  Coordinates,
        # row kind and size must remain stable across both modes.
        assert [(row[:4]) for row in gatk_false] == [(row[:4]) for row in gatk_true]
        assert [(row[:4]) for row in native_false] == [(row[:4]) for row in native_true]
        profile_rows = [row for row in gatk_false if row[3].startswith("size=")]
        assert any(row[4] == "-1.00000" for row in profile_rows), profile_rows
        assert all(row[4] == "1.00000" for row in gatk_true
                   if row[3].startswith("size=")), gatk_true
        assert all(row[4] == "1.00000" for row in native_true
                   if row[3].startswith("size=")), native_true

        # GATK's force-active switch does not alter this fixture's VCF.  The
        # native somatic posterior is intentionally not claimed bit-identical;
        # require only that the candidate site set remains stable while the
        # option changes AssemblyRegion ownership.
        java_records_false = read_vcf_records(java_vcf[False])
        java_records_true = read_vcf_records(java_vcf[True])
        native_records_false = read_vcf_records(native_vcf[False])
        native_records_true = read_vcf_records(native_vcf[True])
        assert java_records_false == java_records_true
        assert vcf_site_keys(native_records_false) == vcf_site_keys(native_records_true)
        for force in (False, True):
            manifest = json.loads(native_manifest[force].read_text(encoding="utf-8"))
            assert manifest["compatibility"]["force_active"] is force

        # A selected interval with no reads remains a legal
        # AssemblyRegionWalker traversal.  The source iterator contributes
        # empty pileups at every locus, cuts the constant inactive profile at
        # maxAssemblyRegionSize, and only then applies --force-active.  This
        # exercises the native sparse Host representation without giving
        # Kokkos one artificial evidence row per zero-coverage coordinate.
        zero_interval = "17:800000-800700"
        zero_gatk_rows: dict[bool, list[tuple[str, int, int, str, str]]] = {}
        zero_native_rows: dict[bool, list[tuple[str, int, int, str, str]]] = {}
        for force in (False, True):
            suffix = "true" if force else "false"
            gatk_vcf = work / f"zero-gatk-{suffix}.vcf.gz"
            gatk_igv = work / f"zero-gatk-{suffix}.igv"
            native_vcf = work / f"zero-native-{suffix}.vcf.gz"
            native_igv = work / f"zero-native-{suffix}.igv"
            native_manifest_path = work / f"zero-native-{suffix}.manifest.json"
            java_run = run([
                str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
                "-I", str(bam), "--tumor-sample", "NA12878", "-L", zero_interval,
                "-O", str(gatk_vcf), "--assembly-region-out", str(gatk_igv),
                "--force-active", str(force).lower(),
                "--create-output-variant-index", "false",
            ])
            assert java_run.returncode == 0, java_run.stderr
            native_run = run([
                str(binary), "-R", str(reference), "-I", str(bam),
                "--tumor-sample", "NA12878", "-L", zero_interval,
                "-O", str(native_vcf), "--assembly-region-out", str(native_igv),
                "--force-active=" + str(force).lower(),
                "--create-output-variant-index", "false",
                "--min-depth", "1", "--min-alt-support", "1",
                "--output-manifest", str(native_manifest_path),
            ])
            assert native_run.returncode == 0, native_run.stderr
            zero_gatk_rows[force] = read_igv(gatk_igv)
            zero_native_rows[force] = read_igv(native_igv)
            assert zero_gatk_rows[force] == zero_native_rows[force], {
                "force_active": force,
                "expected": zero_gatk_rows[force],
                "actual": zero_native_rows[force],
            }
            assert read_vcf_records(gatk_vcf) == []
            assert read_vcf_records(native_vcf) == []
            zero_manifest = json.loads(native_manifest_path.read_text(encoding="utf-8"))
            assert zero_manifest["compatibility"]["force_active"] is force
        zero_profile_false = [row for row in zero_gatk_rows[False]
                              if row[3].startswith("size=")]
        zero_profile_true = [row for row in zero_gatk_rows[True]
                             if row[3].startswith("size=")]
        assert [(row[1], row[2], row[3]) for row in zero_profile_false] == [
            (799999, 800299, "size=300"),
            (800299, 800599, "size=300"),
            (800599, 800700, "size=101"),
        ], zero_profile_false
        assert all(row[4] == "-1.00000" for row in zero_profile_false)
        assert [(row[:4]) for row in zero_profile_false] == [
            row[:4] for row in zero_profile_true]
        assert all(row[4] == "1.00000" for row in zero_profile_true)

        # With no -L, AssemblyRegionWalker traverses the whole reference
        # dictionary even if the standard MAPQ filter removes every read.
        # This verifies that Mutect2's Host traversal domain is not inferred
        # only from compact activity rows.
        full_gatk_vcf = work / "full-zero-gatk.vcf.gz"
        full_gatk_igv = work / "full-zero-gatk.igv"
        full_native_vcf = work / "full-zero-native.vcf.gz"
        full_native_igv = work / "full-zero-native.igv"
        full_native_manifest = work / "full-zero-native.manifest.json"
        java_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878",
            "-O", str(full_gatk_vcf), "--assembly-region-out", str(full_gatk_igv),
            "--minimum-mapping-quality", "255",
            "--create-output-variant-index", "false",
        ])
        assert java_run.returncode == 0, java_run.stderr
        native_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-O", str(full_native_vcf),
            "--assembly-region-out", str(full_native_igv),
            "--minimum-mapping-quality", "255",
            "--create-output-variant-index", "false",
            "--min-depth", "1", "--min-alt-support", "1",
            "--output-manifest", str(full_native_manifest),
        ])
        assert native_run.returncode == 0, native_run.stderr
        full_gatk_rows = read_igv(full_gatk_igv)
        full_native_rows = read_igv(full_native_igv)
        assert full_gatk_rows == full_native_rows, {
            "expected": full_gatk_rows,
            "actual": full_native_rows,
        }
        full_profile_rows = [row for row in full_gatk_rows if row[3].startswith("size=")]
        assert len(full_profile_rows) == 3334, len(full_profile_rows)
        assert full_profile_rows[0] == ("17", 0, 300, "size=300", "-1.00000")
        assert full_profile_rows[-1] == ("17", 999900, 1000000, "size=100", "-1.00000")
        assert read_vcf_records(full_gatk_vcf) == []
        assert read_vcf_records(full_native_vcf) == []
        print(json.dumps({
            "status": "pass",
            "gatk_profile_rows": len(gatk_false),
            "native_profile_rows": len(native_false),
            "profile_boundaries_exact": True,
            "force_active_profile_state_exact": True,
            "zero_coverage_profile_boundaries_exact": True,
            "full_contig_zero_coverage_profile_exact": True,
            "gatk_vcf_unchanged_by_force_active": True,
            "native_vcf_site_set_stable": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
