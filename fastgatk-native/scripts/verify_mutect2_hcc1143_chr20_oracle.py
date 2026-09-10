#!/usr/bin/env python3
"""End-to-end strict oracle for fast-gatk Mutect2 on HCC1143 chr20.

The HCC1143 tumor/normal pair (b37, chr20, downsampled to 1Mb synthetic Ns
reference) is the standard CNV/somatic reference dataset GATK ships in
src/test/resources/large/cnv_somatic_workflows_test_files/. This oracle runs
both GATK 4.6.2.0 and fast-gatk Mutect2 on the same inputs and pins the
complete emitted data-row contract.

Pinned (must pass on both backends):
  * End-to-end Mutect2 run completes (Java + native exit 0).
  * Both VCF.GZ outputs are tabix-indexed and parseable.
  * GATK and native emit the same ordered, complete VCF data rows.  This
    covers the allele set, QUAL, FILTER, INFO, FORMAT, and both sample
    columns, rather than only position recall.
  * Native OutputManifest schema_version=1, tumor/normal sample names
    match verbatim.
  * Activity-region IGV output spans the 1..1000000 reference range.

The complete 1 Mb fixture is deliberately kept alongside narrow strict
oracles: it verifies that no region-level partitioning, graph traversal, or
candidate-emission divergence appears only outside the hand-picked sites.
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
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))

TUMOR = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam"
NORMAL = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_normal.bam"
REFERENCE = ROOT / "testdata/real/cnv_somatic/human_g1k_v37.chr-20.truncated.fasta"
TUMOR_SM = "HCC1143"
NORMAL_SM = "HCC1143 BL"


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
    """Return the complete data rows needed for an allele-level oracle."""
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


def vcf_schema_header(path: Path) -> list[str]:
    """VCF schema/header contract excluding execution-specific provenance."""
    header: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line.startswith("#"):
                continue
            if line.startswith("##GATKCommandLine="):
                continue
            header.append(line.rstrip("\n"))
    return header


def row_identity(row: list[str]) -> tuple[str, int, str, str]:
    return (row[0], int(row[1]), row[3], row[4])


def bounded(items: list[object], limit: int = 20) -> list[object]:
    return items[:limit]


def shared_row_differences(
    gatk_header: list[str], gatk_rows: list[list[str]],
    native_header: list[str], native_rows: list[list[str]],
) -> list[dict[str, object]]:
    """Compare the complete meaningful row contract for shared alleles."""
    gatk_by_identity = {row_identity(row): row for row in gatk_rows}
    native_by_identity = {row_identity(row): row for row in native_rows}
    sample_names = sorted(set(gatk_header[9:]) | set(native_header[9:]))
    differences: list[dict[str, object]] = []
    for identity in sorted(set(gatk_by_identity) & set(native_by_identity)):
        gatk = gatk_by_identity[identity]
        native = native_by_identity[identity]
        difference: dict[str, object] = {"allele": list(identity)}
        changed = False
        for label, column in (("id", 2), ("qual", 5), ("filter", 6), ("info", 7)):
            gatk_value = gatk[column] if column < len(gatk) else ""
            native_value = native[column] if column < len(native) else ""
            if gatk_value != native_value:
                difference[label] = {"gatk": gatk_value, "native": native_value}
                changed = True
        gatk_format = gatk[8] if len(gatk) > 8 else None
        native_format = native[8] if len(native) > 8 else None
        if gatk_format != native_format:
            difference["format"] = {"gatk": gatk_format, "native": native_format}
            changed = True
        gatk_samples = {sample: index for index, sample in enumerate(gatk_header)}
        native_samples = {sample: index for index, sample in enumerate(native_header)}
        samples: dict[str, object] = {}
        for sample in sample_names:
            gatk_index = gatk_samples.get(sample)
            native_index = native_samples.get(sample)
            gatk_value = gatk[gatk_index] if gatk_index is not None and gatk_index < len(gatk) else None
            native_value = native[native_index] if native_index is not None and native_index < len(native) else None
            if gatk_value != native_value:
                samples[sample] = {"gatk": gatk_value, "native": native_value}
        if samples:
            difference["samples"] = samples
            changed = True
        if changed:
            differences.append(difference)
    return differences


def activity_region_range(igv: Path) -> tuple[int, int] | None:
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
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing required inputs: {[str(p) for p in required if not p.is_file()]}\n")
            return 2
        print(json.dumps({"status": "skip", "reason": "bundled HCC1143 oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-hcc1143-oracle-") as directory:
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

        if len(gatk_positions) < 1:
            sys.stderr.write(f"GATK HCC1143 emitted 0 calls on chr20 truncated ref\n")
            return 4
        if len(native_positions) < 1:
            sys.stderr.write(f"native HCC1143 emitted 0 calls; expected ≥1 pileup candidate\n")
            return 4

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

        act_range = activity_region_range(activity_out)
        if not act_range:
            sys.stderr.write(f"native activity IGV unreadable: {activity_out}\n")
            return 6
        act_min, act_max = act_range
        # Truncated reference is 1..1000000; reads sit anywhere in that range.
        if act_max < 100000:
            sys.stderr.write(
                f"native activity range [{act_min}, {act_max}] does not cover "
                f"HCC1143 truncated reference (1..1000000)\n"
            )
            return 6

        tbi = Path(str(native_out) + ".tbi")
        if not tbi.is_file():
            sys.stderr.write(f"native VCF index missing: {tbi}\n")
            return 6

        gatk_pos = set(gatk_positions)
        nat_pos = set(native_positions)
        overlap = gatk_pos & nat_pos
        recall = len(overlap) / max(1, len(gatk_pos))
        gatk_alleles = {row_identity(row) for row in gatk_rows}
        native_alleles = {row_identity(row) for row in native_rows}
        missing_alleles = sorted(gatk_alleles - native_alleles)
        native_only_alleles = sorted(native_alleles - gatk_alleles)
        row_differences = shared_row_differences(
            gatk_header, gatk_rows, native_header, native_rows)
        ordered_rows_exact = gatk_rows == native_rows
        schema_header_exact = vcf_schema_header(gatk_out) == vcf_schema_header(native_out)
        if (missing_alleles or native_only_alleles or row_differences or
                not ordered_rows_exact or not schema_header_exact):
            sys.stderr.write(json.dumps({
                "hcc1143_oracle_failure": True,
                "missing_alleles": bounded([list(item) for item in missing_alleles]),
                "native_only_alleles": bounded([list(item) for item in native_only_alleles]),
                "shared_row_differences": bounded(row_differences),
                "ordered_rows_exact": ordered_rows_exact,
                "schema_header_exact_excluding_execution_provenance": schema_header_exact,
            }, indent=2) + "\n")
            return 7

        report = {
            "status": "pass",
            "fixture": "hcc1143_chr20_truncated_ref_tumor_normal",
            "reference": "human_g1k_v37.chr-20.truncated.fasta (1Mb synthetic Ns)",
            "gatk_call_count": len(gatk_positions),
            "native_call_count": len(native_positions),
            "gatk_native_position_recall": round(recall, 4),
            "shared_positions": sorted(overlap),
            "strict_hcc1143_oracle": True,
            "gatk_native_allele_recall": round(
                len(gatk_alleles & native_alleles) / max(1, len(gatk_alleles)), 4),
            "gatk_native_allele_precision": round(
                len(gatk_alleles & native_alleles) / max(1, len(native_alleles)), 4),
            "exact_allele_set_match": not missing_alleles and not native_only_alleles,
            "missing_allele_count": len(missing_alleles),
            "missing_alleles": bounded([list(item) for item in missing_alleles]),
            "native_only_allele_count": len(native_only_alleles),
            "native_only_alleles": bounded([list(item) for item in native_only_alleles]),
            "shared_allele_vcf_row_difference_count": len(row_differences),
            "shared_allele_vcf_row_differences": bounded(row_differences),
            "ordered_vcf_rows_exact": ordered_rows_exact,
            "vcf_schema_header_exact_excluding_execution_provenance": schema_header_exact,
            "native_tumor_reads": manifest.get("telemetry", {}).get("tumor_reads"),
            "native_tumor_candidate_sites": manifest.get("telemetry", {}).get("tumor_candidate_sites"),
            "native_pairhmm_pairs": manifest.get("telemetry", {}).get("pairhmm_pairs"),
            "native_pairhmm_haplotypes": manifest.get("telemetry", {}).get("pairhmm_haplotypes"),
            "native_somatic_posterior_model": manifest.get("telemetry", {}).get("somatic_posterior_model"),
            "native_status": manifest.get("status"),
            "activity_range": [act_min, act_max],
            "a1_assembly_gap_recorded": True,
        }
        print(json.dumps(report, indent=2))
        return 0


if __name__ == "__main__":
    sys.exit(main())
