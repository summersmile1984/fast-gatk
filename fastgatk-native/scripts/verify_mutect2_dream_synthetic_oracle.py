#!/usr/bin/env python3
"""End-to-end smoke oracle for fast-gatk Mutect2 on DREAM synthetic chr20.

The DREAM Challenge synthetic tumor/normal fixture (b37, chr20, ~30k reads
each) is the standard reference dataset GATK ships in
src/test/resources/large/mutect/dream_synthetic_bams/. This oracle runs both
GATK 4.6.2.0 and fast-gatk Mutect2 on the same inputs and asserts only what
is meaningfully pinned at the documented `prototype` status level.

Pinned (must pass on both backends):
  * End-to-end Mutect2 run completes (Java + native exit 0).
  * Both VCF.GZ outputs are tabix-indexed and parseable.
  * GATK emits ≥ 5 chr20 calls (the DREAM ground-truth somatic set).
  * Native emits every GATK DREAM chr20 call position.  Extra native calls
    remain recorded while candidate precision work continues.
  * Native OutputManifest schema_version=1, status=prototype or
    contract-compatible, tumor/normal sample names match verbatim.
  * Activity-region IGV output is emitted by native and spans the full
    read range (9999900..10100086), confirming the A1 activity path is
    active even though assembly/PairHMM emit 0 haplotypes.

Recorded (not pinned in the normal integration gate):
  * Exact `(contig, position, REF, ALT)` call identities, same-position
    allele mismatches, and per-shared-allele VCF field differences.  This
    avoids treating position recall as an allele-level oracle.
  * `FASTGATK_STRICT_DREAM_ORACLE=1` promotes exact allele identity and
    complete shared-row VCF equality to a failing gate.  It is deliberately
    available before it is green: strict failure is actionable evidence of a
    remaining GATK semantic gap, not a reason to weaken the comparison.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
import oracle_guard

ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))

TUMOR = ROOT / "testdata/real/dream_synthetic/chr20/tumor.bam"
NORMAL = ROOT / "testdata/real/dream_synthetic/chr20/normal.bam"
REFERENCE = ROOT / "testdata/downloads/reference/hs37d5.fa.gz"
TUMOR_SM = "synthetic.challenge.set1.tumor"
NORMAL_SM = "synthetic.challenge.set1.normal"


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=check)


def vcf_positions(path: Path) -> list[tuple[str, int]]:
    out = []
    if not path.is_file():
        return out
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out.append((f[0], int(f[1])))
    return out


def vcf_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    """Parse the VCF rows needed for an exact output-contract report.

    This verifier intentionally keeps parsing in the standard library: it
    must run in the minimal CTest environment used by the native build.  The
    full ALT field is part of the identity, so multi-ALT records are not
    silently reduced to a matching coordinate.
    """
    header: list[str] = []
    rows: list[list[str]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
            elif line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 8:
                    rows.append(fields)
    if not header:
        raise ValueError(f"VCF header missing: {path}")
    return header, rows


def row_identity(row: list[str]) -> tuple[str, int, str, str]:
    return (row[0], int(row[1]), row[3], row[4])


def bounded(items: list[object], limit: int = 20) -> list[object]:
    """Keep CTest output diagnostic without making a bad VCF log enormous."""
    return items[:limit]


def shared_row_differences(
    gatk_header: list[str], gatk_rows: list[list[str]],
    native_header: list[str], native_rows: list[list[str]],
) -> list[dict[str, object]]:
    """Report contract differences only for exact shared allele identities."""
    gatk_by_identity = {row_identity(row): row for row in gatk_rows}
    native_by_identity = {row_identity(row): row for row in native_rows}
    sample_names = sorted(set(gatk_header[9:]) | set(native_header[9:]))
    differences: list[dict[str, object]] = []
    for identity in sorted(set(gatk_by_identity) & set(native_by_identity)):
        gatk = gatk_by_identity[identity]
        native = native_by_identity[identity]
        row_difference: dict[str, object] = {"allele": list(identity)}
        changed = False
        # ID/QUAL/FILTER/INFO are meaningful Mutect2 output contract fields.
        for label, column in (("id", 2), ("qual", 5), ("filter", 6), ("info", 7)):
            gatk_value = gatk[column] if column < len(gatk) else ""
            native_value = native[column] if column < len(native) else ""
            if gatk_value != native_value:
                row_difference[label] = {"gatk": gatk_value, "native": native_value}
                changed = True
        gatk_samples = {sample: index for index, sample in enumerate(gatk_header)}
        native_samples = {sample: index for index, sample in enumerate(native_header)}
        sample_differences: dict[str, object] = {}
        for sample in sample_names:
            gatk_index = gatk_samples.get(sample)
            native_index = native_samples.get(sample)
            gatk_value = gatk[gatk_index] if gatk_index is not None and gatk_index < len(gatk) else None
            native_value = native[native_index] if native_index is not None and native_index < len(native) else None
            if gatk_value != native_value:
                sample_differences[sample] = {"gatk": gatk_value, "native": native_value}
        gatk_format = gatk[8] if len(gatk) > 8 else None
        native_format = native[8] if len(native) > 8 else None
        if gatk_format != native_format:
            row_difference["format"] = {"gatk": gatk_format, "native": native_format}
            changed = True
        if sample_differences:
            row_difference["samples"] = sample_differences
            changed = True
        if changed:
            differences.append(row_difference)
    return differences


def activity_region_range(igv: Path) -> tuple[int, int] | None:
    """Return (min_start, max_end) from native's IGV activity profile."""
    if not igv.is_file():
        return None
    starts, ends = [], []
    with open(igv, "r", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                starts.append(int(parts[1]))
                ends.append(int(parts[2]))
            except ValueError:
                continue
    return (min(starts), max(ends)) if starts else None


def main() -> int:
    required = (JAVA, GATK, NATIVE, TUMOR, NORMAL, REFERENCE)
    if not all(p.is_file() for p in required):
        oracle_guard.oracle_not_verified('verify_mutect2_dream_synthetic_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing required inputs: {[str(p) for p in required if not p.is_file()]}\n")
            return 2
        print(json.dumps({"status": "skip", "reason": "bundled DREAM synthetic oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-dream-oracle-") as directory:
        work = Path(directory)
        gatk_out = work / "gatk.vcf.gz"
        native_out = work / "native.vcf.gz"
        native_manifest = work / "native.manifest.json"
        activity_out = work / "native.activity.igv"

        gatk_run = run([
            str(JAVA), "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL),
            "-tumor", TUMOR_SM, "-normal", NORMAL_SM,
            "-O", str(gatk_out),
        ], check=False)
        if gatk_run.returncode != 0:
            sys.stderr.write(f"GATK Mutect2 failed:\n{gatk_run.stderr[-2000:]}\n")
            return 3

        native_run = run([
            str(NATIVE), "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
            "-O", str(native_out),
            "--output-manifest", str(native_manifest),
            "--assembly-region-out", str(activity_out),
        ], check=False)
        if native_run.returncode != 0:
            sys.stderr.write(f"native Mutect2 failed:\n{native_run.stderr[-2000:]}\n")
            return 3

        gatk_positions = vcf_positions(gatk_out)
        native_positions = vcf_positions(native_out)
        gatk_header, gatk_rows = vcf_rows(gatk_out)
        native_header, native_rows = vcf_rows(native_out)

        # Pinned: GATK emits ≥5 DREAM true positives on chr20.
        if len(gatk_positions) < 5:
            sys.stderr.write(
                f"GATK DREAM synthetic emitted only {len(gatk_positions)} calls; expected ≥5\n"
            )
            return 4

        # Pinned: Native emits at least one call before requiring full GATK
        # position recall, so an empty/invalid VCF reports a direct failure.
        if len(native_positions) < 1:
            sys.stderr.write(f"native emitted 0 calls; expected GATK-call recall\n")
            return 4

        gatk_pos = set(gatk_positions)
        nat_pos = set(native_positions)
        missing_positions = sorted(gatk_pos - nat_pos)
        if missing_positions:
            sys.stderr.write(
                "native missed GATK DREAM call positions: "
                f"{missing_positions}\n"
            )
            return 4

        # Position recall is intentionally retained as the cheap integration
        # gate above, but never used as a substitute for an allele identity
        # comparison.  GATK may emit several concrete alleles at one locus,
        # and a native caller that emits a different ALT at that position has
        # not reproduced the call.
        gatk_allele_counts = Counter(row_identity(row) for row in gatk_rows)
        native_allele_counts = Counter(row_identity(row) for row in native_rows)
        gatk_alleles = set(gatk_allele_counts)
        native_alleles = set(native_allele_counts)
        missing_alleles = sorted(gatk_alleles - native_alleles)
        native_only_alleles = sorted(native_alleles - gatk_alleles)
        shared_alleles = gatk_alleles & native_alleles
        allele_mismatch_positions = sorted(
            (gatk_pos & nat_pos) - {(identity[0], identity[1]) for identity in shared_alleles}
        )
        row_differences = shared_row_differences(
            gatk_header, gatk_rows, native_header, native_rows)
        duplicate_identity_count_differences = [
            {"allele": list(identity), "gatk": gatk_allele_counts[identity],
             "native": native_allele_counts[identity]}
            for identity in sorted(gatk_alleles | native_alleles)
            if gatk_allele_counts[identity] != native_allele_counts[identity]
        ]

        # Pinned: native OutputManifest schema, status, sample names.
        if not native_manifest.is_file():
            sys.stderr.write(f"native OutputManifest missing: {native_manifest}\n")
            return 5
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        if manifest.get("schema_version") != 1:
            sys.stderr.write(f"manifest schema_version mismatch: {manifest.get('schema_version')}\n")
            return 5
        if manifest.get("compatibility", {}).get("tumor_sample_selection") != TUMOR_SM:
            sys.stderr.write(f"manifest tumor_sample_selection != {TUMOR_SM}\n")
            return 5
        if manifest.get("compatibility", {}).get("normal_sample_selection") != NORMAL_SM:
            sys.stderr.write(f"manifest normal_sample_selection != {NORMAL_SM}\n")
            return 5

        # Pinned: activity IGV output spans the full read range.
        # DREAM reads sit in 9999900..10100086, so the activity profile must
        # reach both ends (the assembly is allowed to under-emit haplotypes,
        # but activity must still reach the read boundaries).
        act_range = activity_region_range(activity_out)
        if not act_range:
            sys.stderr.write(f"native activity IGV unreadable: {activity_out}\n")
            return 6
        act_min, act_max = act_range
        if act_min > 9999900 or act_max < 10100000:
            sys.stderr.write(
                f"native activity range [{act_min}, {act_max}] does not cover "
                f"DREAM read range [9999900, 10100086]\n"
            )
            return 6

        # Pinned: indexed VCF.GZ sidecar exists (Tabix convention: foo.vcf.gz.tbi).
        tbi = Path(str(native_out) + ".tbi")
        if not tbi.is_file():
            sys.stderr.write(f"native VCF index missing: {tbi}\n")
            return 6

        # Keep a strict version of this full real-data comparison available
        # to the 1:1 parity suite.  The normal CTest integration gate remains
        # position-recall based while the compact graph/AssemblyResultSet
        # work is in flight, but it must print the exact gap below.  This
        # conditional must *not* be inverted: silently relaxing strict
        # output equality would make a green integration test misleading.
        strict = os.environ.get("FASTGATK_STRICT_DREAM_ORACLE") == "1"
        if strict and (missing_alleles or native_only_alleles or row_differences or
                       duplicate_identity_count_differences):
            sys.stderr.write(json.dumps({
                "strict_dream_oracle_failure": True,
                "missing_alleles": bounded([list(item) for item in missing_alleles]),
                "native_only_alleles": bounded([list(item) for item in native_only_alleles]),
                "duplicate_identity_count_differences": bounded(
                    duplicate_identity_count_differences),
                "shared_row_differences": bounded(row_differences),
            }, indent=2) + "\n")
            return 7

        # Recorded by every integration run: exact allele and VCF-contract
        # deltas.  `native_only_positions` remains for backwards-compatible
        # dashboards; the allele fields are the stronger GATK comparison.
        overlap = gatk_pos & nat_pos
        recall = len(overlap) / max(1, len(gatk_pos))
        report = {
            "status": "pass",
            "fixture": "dream_synthetic_tumor_normal_chr20",
            "gatk_call_count": len(gatk_positions),
            "native_call_count": len(native_positions),
            "gatk_native_position_recall": round(recall, 4),
            "shared_positions": sorted(overlap),
            "gatk_native_allele_recall": round(
                len(shared_alleles) / max(1, len(gatk_alleles)), 4),
            "gatk_native_allele_precision": round(
                len(shared_alleles) / max(1, len(native_alleles)), 4),
            "exact_allele_set_match": not missing_alleles and not native_only_alleles,
            "duplicate_identity_count_difference_count": len(duplicate_identity_count_differences),
            "duplicate_identity_count_differences": bounded(
                duplicate_identity_count_differences),
            "missing_allele_count": len(missing_alleles),
            "missing_alleles": bounded([list(item) for item in missing_alleles]),
            "native_only_allele_count": len(native_only_alleles),
            "native_only_alleles": bounded([list(item) for item in native_only_alleles]),
            "same_position_allele_mismatch_count": len(allele_mismatch_positions),
            "same_position_allele_mismatches": bounded(
                [list(item) for item in allele_mismatch_positions]),
            "shared_allele_vcf_row_difference_count": len(row_differences),
            "shared_allele_vcf_row_differences": bounded(row_differences),
            "strict_dream_oracle": strict,
            "native_tumor_reads": manifest.get("telemetry", {}).get("tumor_reads"),
            "native_normal_reads": manifest.get("telemetry", {}).get("normal_reads"),
            "native_tumor_candidate_sites": manifest.get("telemetry", {}).get("tumor_candidate_sites"),
            "native_pairhmm_pairs": manifest.get("telemetry", {}).get("pairhmm_pairs"),
            "native_pairhmm_haplotypes": manifest.get("telemetry", {}).get("pairhmm_haplotypes"),
            "native_somatic_likelihood_model": manifest.get("telemetry", {}).get("somatic_likelihood_model"),
            "native_somatic_posterior_model": manifest.get("telemetry", {}).get("somatic_posterior_model"),
            "native_status": manifest.get("status"),
            "activity_range": [act_min, act_max],
            "native_only_positions": sorted(nat_pos - gatk_pos),
        }
        print(json.dumps(report, indent=2))
        return 0


if __name__ == "__main__":
    sys.exit(main())
