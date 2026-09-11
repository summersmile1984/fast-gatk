#!/usr/bin/env python3
"""GATK oracle for DREAM's low-BQ compound-haplotype Mutect2 event.

GATK 4.6.2.0 calls 20:10004241 G>A even though its two tumor ALT bases have
BQ=2: ReadThreadingAssembler retains their upstream-insertion branch and
SomaticGenotypingEngine subsequently accepts the paired normal replay.  This
small window is a fast regression guard for that assembly/normal-haplotype
contract; the full-chr20 DREAM oracle remains the broader end-to-end check.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import oracle_guard


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native/build/fastgatk-mutect2")))
REFERENCE = ROOT / "testdata/downloads/reference/hs37d5.fa.gz"
TUMOR = ROOT / "testdata/real/dream_synthetic/chr20/tumor.bam"
NORMAL = ROOT / "testdata/real/dream_synthetic/chr20/normal.bam"
TUMOR_SM = "synthetic.challenge.set1.tumor"
NORMAL_SM = "synthetic.challenge.set1.normal"
# Keep the low-BQ compound event, then extend upstream/downstream enough to
# exercise matched-normal haplotype evidence. In particular GATK rejects the
# normal-supported candidates at 20:10006819 and 20:10009871; a normal pass
# that discards its local paths spuriously emits both.
WINDOW = "20:10000000-10010000"
TARGET = ("20", "10004241", "G", "A")
# A nearby independent two-SNV region.  GATK's Coverage annotation retains
# only finite AlleleLikelihoods evidence here, which catches an off-by-one
# read included by a coordinate-only INFO/DP count.
COVERAGE_WINDOW = "20:10016000-10018000"
COVERAGE_TARGET = ("20", "10017264", "T", "G")
# This active region has an 11-read A>T branch beside a one-read C>T error.
# GATK ChainPruner excludes the latter before PairHMM; retaining the native
# recombined graph path adds one false ALT observation and changes every
# tumor annotation at the otherwise correct A>T call.
PRUNING_WINDOW = "20:10022000-10023500"
PRUNING_TARGET = ("20", "10022820", "A", "T")
# This narrow interval isolates the retained 11-read A>T graph branch.  GATK
# gives PairHMM exactly its REF/ALT assembled haplotypes; comparing this
# handoff directly prevents an accidental return to synthetic candidate
# combinations that can leave the VCF target superficially unchanged.
PAIRHMM_WINDOW = "20:10022790-10022850"
PAIRHMM_NATIVE_PREFIX = "[FASTGATK_PAIRHMM_REQUEST] "


def vcf(path: Path) -> tuple[list[str], list[list[str]]]:
    header: list[str] = []
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                header = fields
            elif not line.startswith("#"):
                rows.append(fields)
    assert header, path
    return header, rows


def assert_reference_backed_native_path(path: Path) -> None:
    """Reject a superficially matching call that used a native-only fallback.

    Every invocation in this oracle has an indexed reference and reaches an
    active AssemblyRegion.  GATK therefore has no pileup-only or unassigned
    candidate path at these windows.  The telemetry check complements the VCF
    comparison: it prevents an implementation from restoring pre-PairHMM
    evidence after a missing graph/PairHMM owner and still matching the final
    row by coincidence.
    """
    stats_path = Path(str(path) + ".stats.json")
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    expected = {
        "pairhmm_skip_reason": "executed",
        "somatic_pileup_fallback_candidates": 0,
        "pairhmm_unassigned_candidates": 0,
        "assembly_unassigned_candidates": 0,
    }
    for key, value in expected.items():
        assert stats.get(key) == value, (stats_path, key, stats.get(key), value)


def has_target(rows: list[list[str]]) -> bool:
    return any(tuple(row[index] for index in (0, 1, 3, 4)) == TARGET for row in rows)


def target_row(rows: list[list[str]]) -> list[str]:
    return next(row for row in rows if tuple(row[index] for index in (0, 1, 3, 4)) == TARGET)


def call_identities(rows: list[list[str]]) -> set[tuple[str, str, str, str]]:
    return {tuple(row[index] for index in (0, 1, 3, 4)) for row in rows}


def info_value(row: list[str], key: str) -> float:
    return float(info_fields(row)[key].split(",", 1)[0])


def info_text(row: list[str], key: str) -> str:
    return info_fields(row)[key]


def info_fields(row: list[str]) -> dict[str, str]:
    return dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)


def header_definition(path: Path, kind: str, key: str) -> str:
    prefix = f"##{kind}=<ID={key},"
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(prefix):
                return line.rstrip("\n")
    raise AssertionError((path, kind, key))


def gatk_pairhmm_haplotypes(path: Path) -> set[str]:
    return {
        line.split()[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    }


def native_pairhmm_haplotypes(stderr: str) -> set[str]:
    return {
        line[len(PAIRHMM_NATIVE_PREFIX):].split("\t", 1)[0]
        for line in stderr.splitlines()
        if line.startswith(PAIRHMM_NATIVE_PREFIX)
    }


def sample_format(row: list[str], header: list[str], sample: str) -> dict[str, str]:
    sample_index = header.index(sample)
    keys = row[8].split(":")
    values = row[sample_index].split(":")
    assert len(keys) == len(values), (sample, keys, values)
    return dict(zip(keys, values, strict=True))


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_mutect2_dream_low_bq_gatk_oracle.py', JAVA, GATK)
        message = f"missing required inputs: {[str(path) for path in required if not path.is_file()]}"
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise RuntimeError(message)
        print('{"status":"skip","reason":"' + message.replace('"', "'") + '"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-dream-low-bq-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        coverage_gatk_output = work / "coverage-gatk.vcf.gz"
        coverage_native_output = work / "coverage-native.vcf.gz"
        pruning_gatk_output = work / "pruning-gatk.vcf.gz"
        pruning_native_output = work / "pruning-native.vcf.gz"
        pairhmm_gatk_output = work / "pairhmm-gatk.vcf.gz"
        pairhmm_native_output = work / "pairhmm-native.vcf.gz"
        pairhmm_gatk_results = work / "pairhmm-gatk.txt"
        gatk_command = [
            str(JAVA), "-jar", str(GATK), "Mutect2", "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL), "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM, "-L", WINDOW, "-O", str(gatk_output),
        ]
        native_command = [
            str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
            "-L", WINDOW, "-O", str(native_output),
        ]
        gatk_run = subprocess.run(gatk_command, text=True, capture_output=True)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run(native_command, text=True, capture_output=True)
        assert native_run.returncode == 0, native_run.stderr
        coverage_gatk_run = subprocess.run([
            str(JAVA), "-jar", str(GATK), "Mutect2", "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL), "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM, "-L", COVERAGE_WINDOW,
            "-O", str(coverage_gatk_output),
        ], text=True, capture_output=True)
        assert coverage_gatk_run.returncode == 0, coverage_gatk_run.stderr
        coverage_native_run = subprocess.run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
            "-L", COVERAGE_WINDOW, "-O", str(coverage_native_output),
        ], text=True, capture_output=True)
        assert coverage_native_run.returncode == 0, coverage_native_run.stderr
        pruning_gatk_run = subprocess.run([
            str(JAVA), "-jar", str(GATK), "Mutect2", "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL), "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM, "-L", PRUNING_WINDOW,
            "-O", str(pruning_gatk_output),
        ], text=True, capture_output=True)
        assert pruning_gatk_run.returncode == 0, pruning_gatk_run.stderr
        pruning_native_run = subprocess.run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
            "-L", PRUNING_WINDOW, "-O", str(pruning_native_output),
        ], text=True, capture_output=True)
        assert pruning_native_run.returncode == 0, pruning_native_run.stderr
        pairhmm_gatk_run = subprocess.run([
            str(JAVA), "-jar", str(GATK), "Mutect2", "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL), "--tumor-sample", TUMOR_SM,
            "--normal-sample", NORMAL_SM, "-L", PAIRHMM_WINDOW,
            "-O", str(pairhmm_gatk_output), "--pair-hmm-results-file",
            str(pairhmm_gatk_results),
        ], text=True, capture_output=True)
        assert pairhmm_gatk_run.returncode == 0, pairhmm_gatk_run.stderr
        pairhmm_native_run = subprocess.run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
            "-L", PAIRHMM_WINDOW, "-O", str(pairhmm_native_output),
        ], text=True, capture_output=True,
            env=os.environ | {"FASTGATK_DEBUG_PAIRHMM_REQUESTS": "1"})
        assert pairhmm_native_run.returncode == 0, pairhmm_native_run.stderr
        for native_path in (
            native_output,
            coverage_native_output,
            pruning_native_output,
            pairhmm_native_output,
        ):
            assert_reference_backed_native_path(native_path)
        gatk_header, gatk_rows = vcf(gatk_output)
        native_header, native_rows = vcf(native_output)
        for kind, key in (
            ("INFO", "AS_SB_TABLE"), ("INFO", "DP"), ("INFO", "ECNT"),
            ("INFO", "ECNTH"), ("INFO", "MBQ"), ("INFO", "MFRL"),
            ("INFO", "MMQ"), ("INFO", "MPOS"), ("INFO", "NALOD"),
            ("INFO", "NLOD"), ("INFO", "POPAF"), ("INFO", "TLOD"),
            ("FORMAT", "AD"), ("FORMAT", "AF"), ("FORMAT", "DP"),
            ("FORMAT", "F1R2"), ("FORMAT", "F2R1"), ("FORMAT", "FAD"),
            ("FORMAT", "GT"), ("FORMAT", "PGT"), ("FORMAT", "PID"),
            ("FORMAT", "PS"), ("FORMAT", "SB"),
        ):
            assert header_definition(native_output, kind, key) == \
                header_definition(gatk_output, kind, key), (kind, key)
        assert native_header[9:] == gatk_header[9:], (
            native_header[9:], gatk_header[9:])
        assert call_identities(native_rows) == call_identities(gatk_rows), (
            "called alleles", sorted(call_identities(native_rows)),
            sorted(call_identities(gatk_rows)))
        assert has_target(gatk_rows), gatk_rows
        assert has_target(native_rows), native_rows
        gatk_target = target_row(gatk_rows)
        native_target = target_row(native_rows)
        # This target is the original low-BQ compound-haplotype regression.
        # Its full row includes the local EventMap counts, tumor/normal
        # annotations and all serialized likelihood scores; checking scalar
        # tolerances below would allow a changed field to hide behind the
        # otherwise correct allele identity.
        assert native_target == gatk_target, (
            "low-BQ target VCF row", native_target, gatk_target)
        assert set(info_fields(native_target)) == set(info_fields(gatk_target)), (
            "INFO keys", sorted(info_fields(native_target)), sorted(info_fields(gatk_target)))
        assert native_target[6] == gatk_target[6], (
            "FILTER", native_target[6], gatk_target[6])
        for key in ("TLOD", "NLOD", "NALOD"):
            assert math.isclose(info_value(native_target, key), info_value(gatk_target, key),
                                rel_tol=0.0, abs_tol=0.01), (
                key, info_value(native_target, key), info_value(gatk_target, key))
        for key in ("NALOD", "NLOD"):
            assert info_text(native_target, key) == info_text(gatk_target, key), (
                key, info_text(native_target, key), info_text(gatk_target, key))
        assert info_value(native_target, "DP") == info_value(gatk_target, "DP"), (
            "DP", info_value(native_target, "DP"), info_value(gatk_target, "DP"))
        assert info_text(native_target, "AS_SB_TABLE") == info_text(gatk_target, "AS_SB_TABLE"), (
            "AS_SB_TABLE", info_text(native_target, "AS_SB_TABLE"),
            info_text(gatk_target, "AS_SB_TABLE"))
        for key in ("MBQ", "MFRL", "MMQ", "MPOS"):
            assert info_text(native_target, key) == info_text(gatk_target, key), (
                key, info_text(native_target, key), info_text(gatk_target, key))
        # addGenotypes() initially computes fractional fragment effective
        # counts, but the default Mutect2 annotation engine then overwrites
        # AD/DP from read BestAlleles and writes FAD from fragment
        # BestAlleles. Verify the final per-sample VCF contract rather than
        # only the intermediate somatic likelihoods.
        for sample in (NORMAL_SM, TUMOR_SM):
            gatk_format = sample_format(gatk_target, gatk_header, sample)
            native_format = sample_format(native_target, native_header, sample)
            for key in ("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"):
                assert native_format[key] == gatk_format[key], (
                    sample, key, native_format[key], gatk_format[key])
        _, coverage_gatk_rows = vcf(coverage_gatk_output)
        _, coverage_native_rows = vcf(coverage_native_output)
        assert call_identities(coverage_native_rows) == call_identities(coverage_gatk_rows), (
            "coverage-window called alleles", sorted(call_identities(coverage_native_rows)),
            sorted(call_identities(coverage_gatk_rows)))
        coverage_gatk_target = next(
            row for row in coverage_gatk_rows
            if tuple(row[index] for index in (0, 1, 3, 4)) == COVERAGE_TARGET)
        coverage_native_target = next(
            row for row in coverage_native_rows
            if tuple(row[index] for index in (0, 1, 3, 4)) == COVERAGE_TARGET)
        assert info_text(coverage_native_target, "DP") == info_text(coverage_gatk_target, "DP"), (
            "coverage-window INFO/DP", info_text(coverage_native_target, "DP"),
            info_text(coverage_gatk_target, "DP"))
        _, pruning_gatk_rows = vcf(pruning_gatk_output)
        _, pruning_native_rows = vcf(pruning_native_output)
        assert call_identities(pruning_native_rows) == call_identities(pruning_gatk_rows), (
            "pruning-window called alleles", sorted(call_identities(pruning_native_rows)),
            sorted(call_identities(pruning_gatk_rows)))
        pruning_gatk_target = next(
            row for row in pruning_gatk_rows
            if tuple(row[index] for index in (0, 1, 3, 4)) == PRUNING_TARGET)
        pruning_native_target = next(
            row for row in pruning_native_rows
            if tuple(row[index] for index in (0, 1, 3, 4)) == PRUNING_TARGET)
        assert pruning_native_target == pruning_gatk_target, (
            "pruning-window target row", pruning_native_target, pruning_gatk_target)
        pairhmm_gatk_haplotypes = gatk_pairhmm_haplotypes(pairhmm_gatk_results)
        pairhmm_native_haplotypes = native_pairhmm_haplotypes(pairhmm_native_run.stderr)
        assert pairhmm_native_haplotypes == pairhmm_gatk_haplotypes, (
            "pairhmm-window haplotypes", sorted(pairhmm_native_haplotypes),
            sorted(pairhmm_gatk_haplotypes))
        print('{"status":"pass","fixture":"dream_low_bq_compound_haplotype",'
              '"target":"20:10004241:G>A","gatk_target":true,"native_target":true,'
              '"called_alleles_match":true,"coverage_window_dp_match":true,'
              '"pruning_window_vcf_row_match":true,'
              '"pairhmm_haplotype_set_exact":true,'
              '"target_vcf_row_exact":true,'
              '"filter_and_format_contract":true}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
