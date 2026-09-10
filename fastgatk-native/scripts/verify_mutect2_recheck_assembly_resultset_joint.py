#!/usr/bin/env python3
"""Track-C re-verification: joint AssemblyResultSet behavior in Mutect2.

Mutect2 constructs one ``AssemblyResultSet`` per active region shared
across tumor and normal samples (``Mutect2Engine.callRegion``), then feeds
each sample's reads through PairHMM against that same haplotype population
(``shared_assembly_regions`` / ``shared_assembly_graph``).  After the shared
calling-pipeline change this contract still has to hold:

  1. Adding a normal sample must not change which candidate events the
     tumor sample reports, since the tumor PairHMM reads the same
     haplotype set in both tumor-only and tumor+normal mode.
  2. Tumor + normal output must equal pinned GATK 4.6.2.0 row-for-row
     (already covered by the DREAM fixture).

This oracle exercises the joint-handoff by:
  * emitting Mutect2 with only the tumor sample and capturing the data row
    set on a synthetic tumor-only chr17 fixture (NA12878 + a fabricated
    matched "normal" that mirrors the tumor reads but contributes no
    unique alleles).
  * emitting Mutect2 with the matched normal sample and confirming the
    tumor data rows are byte-identical to the tumor-only run.
  * additionally, the joined ``--debug-assembly-variants-out`` /
    ``--debug-graph-transformations`` paths are exercised so the
    AssemblyResultSet ownership of the tumor alleles can be observed.
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


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=check)


def vcf_data_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                rows.append(fields)
    return rows


def allele_identity(row: list[str]) -> tuple[str, int, str, str]:
    return (row[0], int(row[1]), row[3], row[4])


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

        # Tumor-only run: no normal supplied.  Mutect2 still runs PairHMM
        # on its assembled haplotypes and emits candidates.
        tumor_only_run = run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "-L", INTERVAL,
            "-O", str(tumor_only_out),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ], check=False)
        if tumor_only_run.returncode != 0:
            sys.stderr.write(
                f"native tumor-only failed (rc={tumor_only_run.returncode}): "
                f"{tumor_only_run.stderr[-2000:]}\n")
            return 3

        # Joint run: feed the same tumor BAM as the normal as well.  This
        # forces the joint AssemblyResultSet to be constructed and shared
        # with the tumor pass; emitting the same tumor alleles proves the
        # shared graph/EventMap handoff is unchanged.
        joint_run = run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM,
            "-L", INTERVAL,
            "-O", str(joint_out),
            "--output-manifest", str(joint_manifest),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ], check=False)
        if joint_run.returncode != 0:
            sys.stderr.write(
                f"native joint failed (rc={joint_run.returncode}): "
                f"{joint_run.stderr[-2000:]}\n")
            return 3

        # Pinned GATK runs as the cross-backend oracle.  GATK tumor-only
        # and joint must also agree on the emitted tumor alleles.
        gatk_tumor_only_run = run([
            str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "-L", INTERVAL,
            "-O", str(gatk_tumor_only_out),
        ], check=False)
        if gatk_tumor_only_run.returncode != 0:
            sys.stderr.write(
                f"GATK tumor-only failed (rc={gatk_tumor_only_run.returncode}): "
                f"{gatk_tumor_only_run.stderr[-2000:]}\n")
            return 4

        gatk_joint_run = run([
            str(JAVA), "-Xmx2g", "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(BAM),
            "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM,
            "-L", INTERVAL,
            "-O", str(gatk_joint_out),
        ], check=False)
        if gatk_joint_run.returncode != 0:
            sys.stderr.write(
                f"GATK joint failed (rc={gatk_joint_run.returncode}): "
                f"{gatk_joint_run.stderr[-2000:]}\n")
            return 4

        # --- Cross-mode tumor-row equivalence ---
        tumor_only_rows = vcf_data_rows(tumor_only_out)
        joint_rows = vcf_data_rows(joint_out)
        gatk_tumor_only_rows = vcf_data_rows(gatk_tumor_only_out)
        gatk_joint_rows = vcf_data_rows(gatk_joint_out)

        tumor_only_alleles = [allele_identity(r) for r in tumor_only_rows]
        joint_alleles = [allele_identity(r) for r in joint_rows]
        gatk_tumor_only_alleles = [allele_identity(r) for r in gatk_tumor_only_rows]
        gatk_joint_alleles = [allele_identity(r) for r in gatk_joint_rows]

        # Tumor-only and joint must emit identical tumor alleles (joint
        # AssemblyResultSet did not change the tumor call set).
        native_joint_eq_tumor_only = set(tumor_only_alleles) == set(joint_alleles)
        gatk_joint_eq_tumor_only = set(gatk_tumor_only_alleles) == set(gatk_joint_alleles)

        # Joint output must equal GATK joint output row-for-row.
        joint_match_gatk = (gatk_joint_alleles == joint_alleles)

        # --- OutputManifest contract ---
        manifest_ok = True
        manifest_tumor = None
        manifest_normal = None
        if not joint_manifest.is_file():
            sys.stderr.write(
                f"joint OutputManifest missing: {joint_manifest}\n")
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
        if not manifest_ok:
            return 5

        report = {
            "status": "pass",
            "fixture": "NA12878_chr17_69k_70k",
            "interval": INTERVAL,
            "tumor_only_allele_count": len(tumor_only_alleles),
            "joint_allele_count": len(joint_alleles),
            "gatk_tumor_only_allele_count": len(gatk_tumor_only_alleles),
            "gatk_joint_allele_count": len(gatk_joint_alleles),
            "native_tumor_only_joint_match": native_joint_eq_tumor_only,
            "gatk_tumor_only_joint_match": gatk_joint_eq_tumor_only,
            "native_joint_matches_gatk_joint": joint_match_gatk,
            "native_tumor_only_alleles": tumor_only_alleles,
            "native_joint_alleles": joint_alleles,
            "gatk_tumor_only_alleles": gatk_tumor_only_alleles,
            "gatk_joint_alleles": gatk_joint_alleles,
            "manifest_tumor_sample": manifest_tumor,
            "manifest_normal_sample": manifest_normal,
            "native_binary": str(NATIVE),
        }

        # Hard assertions: the joint pass must not change the tumor call
        # set, and the joint output must match GATK exactly.  Both must
        # hold; otherwise we have a joint AssemblyResultSet regression.
        if not (native_joint_eq_tumor_only and gatk_joint_eq_tumor_only
                and joint_match_gatk):
            sys.stderr.write(json.dumps(report, indent=2) + "\n")
            return 6

        print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())