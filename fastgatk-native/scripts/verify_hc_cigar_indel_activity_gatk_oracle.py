#!/usr/bin/env python3
"""Pinned GATK oracle for CIGAR-only indel activity ownership.

The reads in this fixture have matching aligned bases but half carry a 2-base
deletion in their CIGAR.  GATK seeds an AssemblyRegion from that explicit
event; a base-mismatch-only activity pass leaves the candidate unassigned and
PairHMM can discard every read.  This check keeps that ownership boundary
observable without claiming complete HC graph equivalence.
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def write_fixture(work: Path) -> tuple[Path, Path, Path]:
    random.seed(4)
    reference = "".join(random.choice("ACGT") for _ in range(250))
    reference += "CA" * 20
    reference += "".join(random.choice("ACGT") for _ in range(950))
    assert len(reference) == 1240
    ref = work / "ref.fa"
    # One-line FASTA makes the .fai offset deterministic.
    ref.write_text(">chr1\n" + reference + "\n", encoding="utf-8")
    (work / "ref.fa.fai").write_text(
        f"chr1\t{len(reference)}\t6\t{len(reference)}\t{len(reference) + 1}\n",
        encoding="utf-8",
    )
    sam = work / "input.sam"
    # POS is 1-based.  The deletion begins after 60 reference bases and is
    # represented by the CIGAR event rather than an aligned-base mismatch.
    ref_read = reference[199:319]
    del_read = reference[199:259] + reference[261:319]
    assert len(ref_read) == 120 and len(del_read) == 118
    qual_ref = "I" * len(ref_read)
    qual_del = "I" * len(del_read)
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1240",
             "@RG\tID:rg\tSM:S1"]
    for i in range(8):
        lines.append(
            f"alt{i}\t0\tchr1\t200\t60\t60M2D58M\t*\t0\t0\t"
            f"{del_read}\t{qual_del}\tRG:Z:rg")
    for i in range(8):
        lines.append(
            f"ref{i}\t0\tchr1\t200\t60\t120M\t*\t0\t0\t"
            f"{ref_read}\t{qual_ref}\tRG:Z:rg")
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return ref, sam, work / "input.bam"


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def info_map(field: str) -> dict[str, str]:
    return {item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in field.split(";") if "=" in item}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    if not native.exists() or not gatk.exists():
        oracle_guard.oracle_not_verified('verify_hc_cigar_indel_activity_gatk_oracle.py', gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit("missing HC binary or pinned GATK 4.6.2.0 jar")
        print(json.dumps({"status": "skipped", "reason": "GATK jar or native binary absent"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-cigar-indel-oracle-") as directory:
        work = Path(directory)
        ref, sam, bam = write_fixture(work)
        dictionary = work / "ref.dict"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
            "-R", str(ref), "-O", str(dictionary), "--TRUNCATE_NAMES_AT_WHITESPACE", "true",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam), "-O", str(bam),
            "-SO", "coordinate", "--CREATE_INDEX", "true",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        java_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.json"
        common = ["-R", str(ref), "-I", str(bam), "-L", "chr1:200-400",
                  "--min-pruning", "1", "--create-output-variant-index", "false",
                  "--add-output-vcf-command-line", "false"]
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                        "-O", str(java_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(native), *common, "--threads", "2", "--output-manifest", str(manifest),
                        "-O", str(native_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        java_rows = records(java_vcf)
        native_rows = records(native_vcf)
        assert len(java_rows) == 1 and len(native_rows) == 1, (java_rows, native_rows)
        j, n = java_rows[0], native_rows[0]
        # Native and Java use different optional INFO annotation sets here;
        # lock the assembly-owned record, likelihood/genotype payload, and
        # stable core annotations shared by both writers.
        assert j[:7] == n[:7], (j[:7], n[:7])
        assert j[8:] == n[8:], (j[8:], n[8:])
        ji, ni = info_map(j[7]), info_map(n[7])
        for key in ("AC", "AF", "AN", "DP", "MLEAC", "MLEAF", "QD"):
            assert ji.get(key) == ni.get(key), (key, ji, ni)
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        # Ref-vs-Any evaluates the CIGAR-adjacent pileup states as ALT evidence.
        # The exact number of active bases is profile geometry, whereas this
        # ownership regression only requires the indel seed to activate at
        # least one locus and stay assigned to the assembly path below.
        assert telemetry["active_loci"] >= 1
        assert telemetry["assembly_regions"] == 1
        assert telemetry["pairhmm_assembly_region_groups"] == 1
        assert telemetry["pairhmm_unassigned_candidates"] == 0
        assert telemetry["assembly_unassigned_candidates"] == 0
        assert telemetry["pairhmm_reads_disqualified"] == 0

        # `callRegion()` is never entered when the activity threshold rejects
        # this CIGAR-only event.  In that case Java must emit no concrete VCF
        # row, and Native must not let its raw CIGAR candidate/compact graph
        # fallback create a contig-wide PairHMM transaction.  This pairs the
        # positive indel-activity contract above with its source ownership
        # boundary rather than merely testing that the same deletion can be
        # represented when a region exists.
        inactive_common = [*common, "--active-probability-threshold", "1.0"]
        inactive_java_vcf = work / "gatk.inactive.vcf"
        inactive_native_vcf = work / "native.inactive.vcf"
        inactive_manifest = work / "native.inactive.json"
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                        *inactive_common, "-O", str(inactive_java_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(native), *inactive_common, "--threads", "2",
                        "--output-manifest", str(inactive_manifest),
                        "-O", str(inactive_native_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert records(inactive_java_vcf) == []
        assert records(inactive_native_vcf) == []
        inactive_telemetry = json.loads(
            inactive_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert inactive_telemetry["assembly_regions"] == 0
        assert inactive_telemetry["candidate_sites"] == 0
        assert inactive_telemetry["variant_calls"] == 0
        assert inactive_telemetry["pairhmm_assembly_region_groups"] == 0
        assert inactive_telemetry["pairhmm_unassigned_candidates"] == 0

        # `--min-base-quality-score` / `-mbq` belongs to the normal HC
        # argument surface (not Mutect2's pileup-only threshold).  All reads
        # in this fixture are Q40, so Q60 must remove the CIGAR deletion
        # before ActivityProfile/assembly in both implementations.  Exercise
        # both GATK spellings and make native telemetry prove that neither was
        # silently accepted and ignored.
        low_quality_common = [*common, "--min-base-quality-score", "60"]
        long_java_vcf = work / "gatk.min-base-quality-score.vcf"
        long_native_vcf = work / "native.min-base-quality-score.vcf"
        long_manifest = work / "native.min-base-quality-score.json"
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                        *low_quality_common, "-O", str(long_java_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(native), *low_quality_common, "--threads", "2",
                        "--output-manifest", str(long_manifest),
                        "-O", str(long_native_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        short_java_vcf = work / "gatk.mbq.vcf"
        short_native_vcf = work / "native.mbq.vcf"
        short_manifest = work / "native.mbq.json"
        short_quality_common = [*common, "-mbq", "60"]
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                        *short_quality_common, "-O", str(short_java_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(native), *short_quality_common, "--threads", "2",
                        "--output-manifest", str(short_manifest),
                        "-O", str(short_native_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert records(long_java_vcf) == records(short_java_vcf) == []
        assert records(long_native_vcf) == records(short_native_vcf) == []
        assert json.loads(long_manifest.read_text(encoding="utf-8"))["telemetry"][
            "min_base_quality"] == 60
        assert json.loads(short_manifest.read_text(encoding="utf-8"))["telemetry"][
            "min_base_quality"] == 60

        # ReferenceConfidenceModel assigns CIGAR deletion pileup elements the
        # configurable --reference-model-deletion-quality.  Prune the eight
        # deletion reads out of the assembly graph so this is a pure gVCF RCM
        # boundary, then compare a low and high quality setting independently
        # against GATK.  The Host prepares the deletion observation while the
        # same Kokkos Ref-vs-Any/PL kernel owns the numeric reduction.
        rcm_common = ["-R", str(ref), "-I", str(bam), "-L", "chr1:200-400",
                      "-ERC", "BP_RESOLUTION", "--min-pruning", "9",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false"]
        rcm_rows: dict[int, list[list[str]]] = {}
        for deletion_quality in (10, 60):
            java_rcm = work / f"gatk.rcm-del-{deletion_quality}.vcf"
            native_rcm = work / f"native.rcm-del-{deletion_quality}.vcf"
            rcm_manifest = work / f"native.rcm-del-{deletion_quality}.json"
            rcm_args = [*rcm_common, "--reference-model-deletion-quality",
                        str(deletion_quality)]
            subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                            *rcm_args, "-O", str(java_rcm)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run([str(native), *rcm_args, "--threads", "2",
                            "--output-manifest", str(rcm_manifest),
                            "-O", str(native_rcm)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            java_rcm_rows = records(java_rcm)
            native_rcm_rows = records(native_rcm)
            assert java_rcm_rows == native_rcm_rows, (
                deletion_quality, java_rcm_rows, native_rcm_rows)
            assert json.loads(rcm_manifest.read_text(encoding="utf-8"))["telemetry"][
                "reference_model_deletion_quality"] == deletion_quality
            rcm_rows[deletion_quality] = native_rcm_rows
        assert rcm_rows[10] != rcm_rows[60]
        print(json.dumps({
            "status": "pass", "release": "GATK 4.6.2.0",
            "fixture": "16 reads; 8x 60M2D58M + 8x 120M",
            "call": {"contig": n[0], "pos": int(n[1]), "ref": n[3], "alt": n[4]},
            "java_native_core_record_exact": True,
            "java_native_genotype_payload_exact": True,
            "activity_indel_seed_active_loci": telemetry["active_loci"],
            "assembly_regions": telemetry["assembly_regions"],
            "pairhmm_unassigned_candidates": telemetry["pairhmm_unassigned_candidates"],
            "inactive_region_has_no_call": True,
            "min_base_quality_score_and_mbq_gatk_exact": True,
            "reference_model_deletion_quality_gatk_exact": True,
            "scope": "CIGAR indel activity and AssemblyRegion/PairHMM ownership",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
