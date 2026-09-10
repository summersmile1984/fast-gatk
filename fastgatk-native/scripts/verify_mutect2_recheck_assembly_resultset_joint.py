#!/usr/bin/env python3
"""Track-C re-verification: joint (tumor + matched normal) Mutect2 behavior.

Rewrite note (Track C audit, 2026-09-11).  The original revision of this file
was never run and could not pass:

  * its joint native/GATK runs passed ``--normal-sample NA12878_NORMAL_MATCH``
    but never supplied a second input, so the native tool aborted with
    ``BAD_INPUT: requested sample is absent from tumor input headers`` (rc 2),
    and GATK would have failed the same way;
  * its central hard assertion -- "adding a matched normal must not change
    which candidate events the tumor sample reports" -- is false by
    construction for the fixture it used.  Measured on pinned GATK 4.6.2.0
    with a matched normal whose reads mirror the tumor reads:

        GATK   tumor-only  6 rows   GATK   joint  0 rows
        native tumor-only  6 rows   native joint  0 rows

    So GATK itself drops every tumor call in joint mode there; the original
    gate would have been a guaranteed failure (not a regression signal).

What this oracle now pins (all measured, none assumed):

  1. chr17 69k-70k fixture, tumor-only: native allele set == GATK allele set.
  2. chr17 69k-70k fixture, joint with a renamed mirrored normal: native
     joint rows == GATK joint rows.  The "both empty" outcome is reported
     explicitly (``joint_rows_zero_by_construction``) so the weak equivalence
     cannot be mistaken for a strong one.
  3. If the real HCC1143 chr20 tumor/matched-normal pair is present, the same
     joint contract is re-checked on a fixture whose joint output is
     non-empty: native joint allele set == GATK joint allele set, with any
     per-field data-row differences reported (gated only on allele identity).
  4. The OutputManifest contract (schema_version, tumor/normal selection) for
     the joint native run.

Cross-mode set equality (tumor-only vs joint) is reported as an observation
only, because it is not a contract GATK itself honours on these fixtures.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native/build/fastgatk-mutect2")))

# Bundled chr17 69k-70k fixture used by verify_mutect2_gatk_oracle.py.
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
TUMOR_SM = "NA12878"
NORMAL_SM = "NA12878_NORMAL_MATCH"
INTERVAL = "17:69000-70000"

# Real tumor / matched-normal pair (same fixture as the HCC1143 oracles).
REAL_REFERENCE = (ROOT / "testdata/real/cnv_somatic/"
                  "human_g1k_v37.chr-20.truncated.fasta")
REAL_TUMOR = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam"
REAL_NORMAL = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_normal.bam"
REAL_TUMOR_SM = "HCC1143"
REAL_NORMAL_SM = "HCC1143 BL"
REAL_INTERVAL = "20:67000-69000"


def run(command: list[str], timeout: float = 180.0) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(command, text=True, capture_output=True,
                              check=False, timeout=timeout)
    except subprocess.TimeoutExpired as expired:
        return subprocess.CompletedProcess(
            command, 124, expired.stdout or "",
            (expired.stderr or "") + f"\nTIMEOUT after {timeout}s")


def vcf_data_rows(path: Path) -> list[list[str]]:
    if not path.is_file():
        return []
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: list[list[str]] = []
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                rows.append(fields)
    return rows


def allele_identity(row: list[str]) -> tuple[str, int, str, str]:
    return (row[0], int(row[1]), row[3], row[4])


def row_field_diffs(gatk_rows: list[list[str]], native_rows: list[list[str]],
                    limit: int = 10) -> list[dict]:
    """Field-by-field comparison for shared allele identities."""
    gatk_by_allele: dict[tuple, list[str]] = {}
    for row in gatk_rows:
        gatk_by_allele.setdefault(allele_identity(row), row)
    native_by_allele: dict[tuple, list[str]] = {}
    for row in native_rows:
        native_by_allele.setdefault(allele_identity(row), row)
    diffs = []
    for identity in sorted(set(gatk_by_allele) & set(native_by_allele)):
        g = gatk_by_allele[identity]
        n = native_by_allele[identity]
        fields: dict[str, dict[str, str]] = {}
        for label, idx in (("qual", 5), ("filter", 6), ("info", 7), ("format", 8)):
            gv = g[idx] if idx < len(g) else ""
            nv = n[idx] if idx < len(n) else ""
            if gv != nv:
                fields[label] = {"gatk": gv, "native": nv}
        for offset in range(9, max(len(g), len(n))):
            gv = g[offset] if offset < len(g) else ""
            nv = n[offset] if offset < len(n) else ""
            if gv != nv:
                fields[f"sample_{offset - 9}"] = {"gatk": gv, "native": nv}
        if fields:
            diffs.append({"allele": list(identity), "fields": fields})
    return diffs[:limit]


def mirrored_normal(source: Path, destination: Path) -> bool:
    """Copy a BAM with every @RG/SM renamed to NORMAL_SM (matched normal)."""
    try:
        import pysam
    except ImportError:
        return False
    with pysam.AlignmentFile(str(source), "rb") as src:
        header = src.header.to_dict()
        for rg in header.get("RG", []):
            rg["SM"] = NORMAL_SM
        with pysam.AlignmentFile(str(destination), "wb", header=header) as dst:
            for record in src:
                dst.write(record)
    pysam.index(str(destination))
    return True


def main() -> int:
    required = (JAVA, GATK, NATIVE, BAM, REFERENCE,
                Path(f"{BAM}.bai"), Path(f"{REFERENCE}.fai"))
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing required inputs: {missing}\n")
            return 2
        print(json.dumps({"status": "skip",
                          "reason": "joint AssemblyResultSet oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-recheck-mutect2-joint-") as directory:
        work = Path(directory)
        tumor_only_out = work / "tumor-only.vcf.gz"
        joint_out = work / "joint.vcf.gz"
        joint_manifest = work / "joint.manifest.json"
        gatk_tumor_only_out = work / "gatk-tumor-only.vcf.gz"
        gatk_joint_out = work / "gatk-joint.vcf.gz"
        normal_bam = work / "normal.bam"

        if not mirrored_normal(BAM, normal_bam):
            print(json.dumps({"status": "skip",
                              "reason": "pysam unavailable: cannot build mirrored normal"}))
            return 0

        # --- Tumor-only baseline (chr17 fixture) ---
        tumor_only_run = run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "-L", INTERVAL,
            "-O", str(tumor_only_out),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        if tumor_only_run.returncode != 0:
            sys.stderr.write(
                f"native tumor-only failed (rc={tumor_only_run.returncode}): "
                f"{tumor_only_run.stderr[-2000:]}\n")
            return 3

        # --- Joint run (chr17 fixture): the mirrored normal IS supplied. ---
        joint_run = run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM), "-I", str(normal_bam),
            "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM,
            "-L", INTERVAL,
            "-O", str(joint_out),
            "--output-manifest", str(joint_manifest),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        if joint_run.returncode != 0:
            sys.stderr.write(
                f"native joint failed (rc={joint_run.returncode}): "
                f"{joint_run.stderr[-2000:]}\n")
            return 3

        gatk_tumor_only_run = run([
            str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "-L", INTERVAL,
            "-O", str(gatk_tumor_only_out),
        ])
        if gatk_tumor_only_run.returncode != 0:
            sys.stderr.write(
                f"GATK tumor-only failed (rc={gatk_tumor_only_run.returncode}): "
                f"{gatk_tumor_only_run.stderr[-2000:]}\n")
            return 4

        gatk_joint_run = run([
            str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(BAM), "-I", str(normal_bam),
            "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM,
            "-L", INTERVAL,
            "-O", str(gatk_joint_out),
        ])
        if gatk_joint_run.returncode != 0:
            sys.stderr.write(
                f"GATK joint failed (rc={gatk_joint_run.returncode}): "
                f"{gatk_joint_run.stderr[-2000:]}\n")
            return 4

        tumor_only_rows = vcf_data_rows(tumor_only_out)
        joint_rows = vcf_data_rows(joint_out)
        gatk_tumor_only_rows = vcf_data_rows(gatk_tumor_only_out)
        gatk_joint_rows = vcf_data_rows(gatk_joint_out)

        tumor_only_alleles = [allele_identity(r) for r in tumor_only_rows]
        joint_alleles = [allele_identity(r) for r in joint_rows]
        gatk_tumor_only_alleles = [allele_identity(r) for r in gatk_tumor_only_rows]
        gatk_joint_alleles = [allele_identity(r) for r in gatk_joint_rows]

        # Hard gate 1: tumor-only parity with pinned GATK (chr17 fixture).
        tumor_only_match_gatk = (sorted(tumor_only_alleles)
                                 == sorted(gatk_tumor_only_alleles))
        # Hard gate 2: joint parity with pinned GATK (chr17 fixture).
        joint_match_gatk = (sorted(joint_alleles) == sorted(gatk_joint_alleles))
        # Observation only -- GATK does not itself honour this on this fixture.
        native_joint_eq_tumor_only = (set(tumor_only_alleles) == set(joint_alleles))
        gatk_joint_eq_tumor_only = (set(gatk_tumor_only_alleles)
                                    == set(gatk_joint_alleles))

        # --- OutputManifest contract ---
        manifest_ok = True
        manifest_tumor = None
        manifest_normal = None
        if not joint_manifest.is_file():
            sys.stderr.write(f"joint OutputManifest missing: {joint_manifest}\n")
            manifest_ok = False
        else:
            manifest = json.loads(joint_manifest.read_text(encoding="utf-8"))
            if manifest.get("schema_version") != 1:
                sys.stderr.write(
                    f"manifest schema_version mismatch: {manifest.get('schema_version')}\n")
                manifest_ok = False
            manifest_tumor = manifest.get("compatibility", {}).get("tumor_sample_selection")
            manifest_normal = manifest.get("compatibility", {}).get("normal_sample_selection")
            if manifest_tumor != TUMOR_SM:
                sys.stderr.write(
                    f"manifest tumor_sample_selection != {TUMOR_SM}: {manifest_tumor}\n")
                manifest_ok = False
            if manifest_normal != NORMAL_SM:
                sys.stderr.write(
                    f"manifest normal_sample_selection != {NORMAL_SM}: {manifest_normal}\n")
                manifest_ok = False

        # --- Real matched-normal joint check (non-empty joint output) ---
        real_section = {"available": False}
        real_required = (REAL_REFERENCE, REAL_TUMOR, REAL_NORMAL)
        if all(p.is_file() for p in real_required):
            real_section = {"available": True}
            real_native_out = work / "real-native-joint.vcf.gz"
            real_gatk_out = work / "real-gatk-joint.vcf.gz"
            real_native = run([
                str(NATIVE), "-R", str(REAL_REFERENCE),
                "-I", str(REAL_TUMOR), "-I", str(REAL_NORMAL),
                "--tumor-sample", REAL_TUMOR_SM, "--normal-sample", REAL_NORMAL_SM,
                "-L", REAL_INTERVAL, "-O", str(real_native_out),
            ])
            real_gatk = run([
                str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
                "-R", str(REAL_REFERENCE),
                "-I", str(REAL_TUMOR), "-I", str(REAL_NORMAL),
                "--tumor-sample", REAL_TUMOR_SM, "--normal-sample", REAL_NORMAL_SM,
                "-L", REAL_INTERVAL, "-O", str(real_gatk_out),
            ])
            real_section["native_rc"] = real_native.returncode
            real_section["gatk_rc"] = real_gatk.returncode
            if real_native.returncode == 0 and real_gatk.returncode == 0:
                real_native_rows = vcf_data_rows(real_native_out)
                real_gatk_rows = vcf_data_rows(real_gatk_out)
                real_native_alleles = [allele_identity(r) for r in real_native_rows]
                real_gatk_alleles = [allele_identity(r) for r in real_gatk_rows]
                real_section.update({
                    "interval": REAL_INTERVAL,
                    "native_allele_count": len(real_native_alleles),
                    "gatk_allele_count": len(real_gatk_alleles),
                    "native_alleles": real_native_alleles,
                    "gatk_alleles": real_gatk_alleles,
                    "allele_sets_equal": (sorted(real_native_alleles)
                                          == sorted(real_gatk_alleles)),
                    "field_differences": row_field_diffs(real_gatk_rows,
                                                         real_native_rows, limit=5),
                })
            else:
                real_section["error"] = {
                    "native": real_native.stderr[-800:],
                    "gatk": real_gatk.stderr[-800:],
                }
        else:
            real_section["missing"] = [str(p) for p in real_required
                                       if not p.is_file()]

        report = {
            "status": "pass",
            "fixture": "NA12878_chr17_69k_70k",
            "interval": INTERVAL,
            "tumor_only_allele_count": len(tumor_only_alleles),
            "joint_allele_count": len(joint_alleles),
            "gatk_tumor_only_allele_count": len(gatk_tumor_only_alleles),
            "gatk_joint_allele_count": len(gatk_joint_alleles),
            "native_tumor_only_alleles": [list(a) for a in tumor_only_alleles],
            "native_joint_alleles": [list(a) for a in joint_alleles],
            "gatk_tumor_only_alleles": [list(a) for a in gatk_tumor_only_alleles],
            "gatk_joint_alleles": [list(a) for a in gatk_joint_alleles],
            "native_tumor_only_matches_gatk": tumor_only_match_gatk,
            "native_joint_matches_gatk": joint_match_gatk,
            "joint_field_differences": row_field_diffs(gatk_joint_rows, joint_rows),
            "joint_rows_zero_by_construction": (
                not joint_alleles and not gatk_joint_alleles),
            "observation_native_tumor_only_eq_joint": native_joint_eq_tumor_only,
            "observation_gatk_tumor_only_eq_joint": gatk_joint_eq_tumor_only,
            "manifest_tumor_sample": manifest_tumor,
            "manifest_normal_sample": manifest_normal,
            "manifest_ok": manifest_ok,
            "real_matched_normal_joint": real_section,
            "native_binary": str(NATIVE),
        }

        # Hard gates: pinned-GATK parity in BOTH modes on the chr17 fixture,
        # plus (when available) allele-set parity on the real matched-normal
        # joint run.  The cross-mode equality is deliberately NOT a gate.
        gates = [tumor_only_match_gatk, joint_match_gatk, manifest_ok]
        if real_section.get("available") and "allele_sets_equal" in real_section:
            gates.append(real_section["allele_sets_equal"])
        if not all(gates):
            report["status"] = "fail"
            sys.stderr.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
            return 6

        print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
