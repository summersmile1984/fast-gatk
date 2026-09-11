#!/usr/bin/env python3
"""Pin HaplotypeCaller behavior on the second real-data fixture (chr20 mnp.bam).

The real CEUTrio NA12878-derived reads (fixtures/chr20/mnp.bam, b37 chr20
10019901-10021411) carry genuine variants whose allele identity is pinned
against the bundled GATK 4.6.2.0 execution for the isolated windows.
Isolated single-window runs pin exact (POS, REF, ALT) plus byte determinism.
A contiguous single -L covering those windows must emit every isolated-window
pinned allele with the same REF/ALT (no dropped interior sites, no allele
identity flips), then match GATK's complete ordered VCF data rows and header
schema (except invocation-specific command provenance).  That dense contiguous
cluster also pins the complete,
sequence-unique PairHMM haplotype population emitted by GATK: graph ownership
must not be replaced by a Cartesian product of raw pileup candidates.

The VCF row fields remain an independent annotation/genotyping oracle.  The
same dense span is also run in default GVCF mode: its many non-variant RCM
blocks sit between phased candidate islands, which catches a Host
read-to-haplotype projection or CIGAR-coordinate error that a normal VCF
cannot observe.  The test specifically protects the assembled-path handoff to
the Kokkos PairHMM kernel while keeping graph/event/RCM ownership on the C++
Host.
"""
from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import sys
import oracle_guard

ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_HC_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-hc-call")))
REFERENCE = ROOT / "fixtures" / "chr20" / "ref20mnp.fasta"
BAM = ROOT / "fixtures" / "chr20" / "mnp.bam"
BAM_INDEX = ROOT / "fixtures" / "chr20" / "mnp.bam.bai"
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"

# Interval A contains 10019967/10019969; interval B contains the 10020228
# SNP and the downstream four-SNP haplotype.
# Positions are the
# 1-based genomic coordinates of the reads.  Each interval is run as its own
# single -L window, then the union is re-run as one contiguous -L.
# The dense cluster is still a window-geometry regression surface.  The
# isolated pins lock GATK-compatible native alleles, and the contiguous
# assertion requires those same alleles rather than byte-identical records.
EXPECTED_RECORDS_BY_INTERVAL = {
    ("20:10019901-10019970",): {
        "10019967": ("C", "G"),
        "10019969": ("T", "G"),
    },
    ("20:10020200-10020500",): {
        "10020228": ("G", "A"),
        "10020429": ("C", "G"),
        "10020431": ("T", "G"),
        "10020434": ("A", "G"),
        "10020438": ("T", "G"),
    },
    ("20:10020390-10020480",): {
        "10020429": ("C", "G"),
        "10020431": ("T", "G"),
        "10020434": ("A", "G"),
        "10020438": ("T", "G"),
    },
    ("20:10020650-10020710",): {
        "10020679": ("A", "T"),
        "10020680": ("C", "A"),
        "10020681": ("A", "T"),
    },
}

CONTIGUOUS_INTERVAL = ("20:10019901-10020710",)
# Full bundled-GATK 4.6.2.0 oracle for the contiguous window.  This is kept
# separately from the isolated pins because a graph/EventMap ownership bug
# can add a high-confidence raw pileup call without changing any isolated
# target allele.  In particular, the former native fallback emitted a
# non-EventMap 20:10019947 G>A row here; Java emits only this 11-allele set.
EXPECTED_CONTIGUOUS_RECORDS = {
    "10019967": ("C", "G"),
    "10019969": ("T", "G"),
    "10020228": ("G", "A"),
    "10020229": ("T", "G"),
    "10020429": ("C", "G"),
    "10020431": ("T", "G"),
    "10020434": ("A", "G"),
    "10020438": ("T", "G"),
    "10020679": ("A", "T"),
    "10020680": ("C", "A"),
    "10020681": ("A", "T"),
}

# Proper subspans of the full contiguous interval that fully contain two or
# more isolated pin windows.  The ≥8-island split fires on the 810bp full
# window; these shorter spans are the remaining A-line hole.
CONTIGUOUS_SUBSPANS = (
    "20:10019901-10020250",  # isolated windows 1+2; 350bp, 20228 near -L end
    "20:10020200-10020480",  # isolated windows 2+3; 281bp, below the 310 band
    "20:10019901-10020500",  # isolated windows 1+2+3
    "20:10020200-10020710",  # isolated windows 2+3+4
    "20:10020390-10020710",  # isolated windows 3+4 (321bp, <8 islands)
)



def records(path: pathlib.Path) -> dict[str, tuple[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    out: dict[str, tuple[str, str]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        out[fields[1]] = (fields[3], fields[4])
    return out


def vcf_lines(path: pathlib.Path) -> tuple[list[str], list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        lines = [line.rstrip("\n") for line in stream]
    header = [line for line in lines if line.startswith("#")]
    rows = [line for line in lines if line and not line.startswith("#")]
    # The command line includes the executable path and run environment, so
    # it is intentionally not a semantic HC compatibility assertion.
    header = [line for line in header if not line.startswith("##GATKCommandLine=")]
    return header, rows


def gvcf_candidate_positions(path: pathlib.Path) -> list[int]:
    """Return actual-ALT rows, excluding reference-confidence-only blocks."""
    positions: list[int] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[4] != "<NON_REF>":
            positions.append(int(fields[1]))
    return positions


def gvcf_row(path: pathlib.Path, position: int) -> list[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if int(fields[1]) == position:
            return fields
    raise AssertionError(f"missing gVCF row at {position}")


def run(output: pathlib.Path, manifest: pathlib.Path, intervals: tuple[str, ...],
        *, debug_pairhmm: bool = False) -> str:
    command = [
        str(BINARY), "-R", str(REFERENCE), "-I", str(BAM),
        "-O", str(output), "--output-manifest", str(manifest), "--threads", "4",
        "--add-output-vcf-command-line", "false",
    ]
    for interval in intervals:
        command += ["-L", interval]
    environment = os.environ.copy()
    if debug_pairhmm:
        environment["FASTGATK_DEBUG_PAIRHMM_REQUESTS"] = "1"
    result = subprocess.run(command, text=True, capture_output=True, check=False,
                            env=environment)
    if result.returncode != 0:
        raise AssertionError(f"native HC failed: {result.stderr[-1200:]}")
    return result.stderr


def gatk_pairhmm_haplotypes(path: pathlib.Path) -> set[str]:
    return {
        line.split(maxsplit=1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    }


def native_pairhmm_haplotypes(stderr: str) -> set[str]:
    prefix = "[FASTGATK_PAIRHMM_REQUEST] "
    return {
        line[len(prefix):].split("\t", maxsplit=1)[0]
        for line in stderr.splitlines()
        if line.startswith(prefix)
    }


def main() -> int:
    if not (BINARY.is_file() and REFERENCE.is_file() and BAM.is_file()
            and BAM_INDEX.is_file() and JAVA.is_file() and GATK.is_file()):
        oracle_guard.oracle_not_verified('verify_hc_chr20_real_contract.py', JAVA, GATK)
        raise SystemExit("fastgatk-hc-chr20-real-contract inputs are required")
    import tempfile
    with tempfile.TemporaryDirectory(prefix="hc-chr20-real-") as temporary:
        work = pathlib.Path(temporary)
        total_calls = 0
        deterministic = True
        for ordinal, (intervals, expected_records) in enumerate(EXPECTED_RECORDS_BY_INTERVAL.items()):
            first = work / f"run{ordinal}.vcf.gz"
            first_json = work / f"run{ordinal}.json"
            run(first, first_json, intervals)
            called = records(first)
            assert called == expected_records, {"native": called,
                                                 "expected": expected_records,
                                                 "intervals": intervals}
            second = work / f"run{ordinal}b.vcf.gz"
            second_json = work / f"run{ordinal}b.json"
            run(second, second_json, intervals)
            deterministic = deterministic and (first.read_bytes() == second.read_bytes())
            payload = json.loads(first_json.read_text(encoding="utf-8"))
            total_calls += payload["telemetry"]["variant_calls"]
            assert payload["reads"] >= 10, payload["reads"]
        assert deterministic, "run is not deterministic"
        assert total_calls == sum(len(v) for v in EXPECTED_RECORDS_BY_INTERVAL.values()), total_calls
        isolated_alleles: dict[str, tuple[str, str]] = {}
        for expected_records in EXPECTED_RECORDS_BY_INTERVAL.values():
            isolated_alleles.update(expected_records)
        contiguous = work / "contiguous.vcf.gz"
        contiguous_json = work / "contiguous.json"
        contiguous_stderr = run(contiguous, contiguous_json, CONTIGUOUS_INTERVAL,
                                debug_pairhmm=True)
        contiguous_called = records(contiguous)
        assert contiguous_called == EXPECTED_CONTIGUOUS_RECORDS, {
            "native": contiguous_called,
            "gatk": EXPECTED_CONTIGUOUS_RECORDS,
            "intervals": CONTIGUOUS_INTERVAL,
        }
        missing = {pos: allele for pos, allele in isolated_alleles.items()
                   if pos not in contiguous_called}
        flipped = {pos: (isolated_alleles[pos], contiguous_called[pos])
                   for pos in isolated_alleles
                   if pos in contiguous_called and contiguous_called[pos] != isolated_alleles[pos]}
        assert not missing, {"missing": missing, "contiguous": contiguous_called}
        assert not flipped, {"flipped": flipped, "contiguous": contiguous_called}

        gatk_contiguous = work / "gatk-contiguous.vcf.gz"
        gatk_pairhmm = work / "gatk-contiguous.pairhmm.txt"
        gatk_command = [
            str(JAVA), "-Xmx1g", "-jar", str(GATK), "HaplotypeCaller",
            "-R", str(REFERENCE), "-I", str(BAM),
            "-L", CONTIGUOUS_INTERVAL[0], "-O", str(gatk_contiguous),
            "--pair-hmm-results-file", str(gatk_pairhmm),
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--native-pair-hmm-threads", "1",
        ]
        gatk_result = subprocess.run(gatk_command, text=True, capture_output=True,
                                     check=False)
        if gatk_result.returncode != 0:
            raise AssertionError(f"GATK HC failed: {gatk_result.stderr[-1200:]}")
        gatk_haplotypes = gatk_pairhmm_haplotypes(gatk_pairhmm)
        native_haplotypes = native_pairhmm_haplotypes(contiguous_stderr)
        assert native_haplotypes == gatk_haplotypes, {
            "missing_native": sorted(gatk_haplotypes - native_haplotypes),
            "native_only": sorted(native_haplotypes - gatk_haplotypes),
        }
        native_header, native_rows = vcf_lines(contiguous)
        gatk_header, gatk_rows = vcf_lines(gatk_contiguous)
        assert native_header == gatk_header, {
            "native_header": native_header,
            "gatk_header": gatk_header,
        }
        assert native_rows == gatk_rows, {
            "native_rows": native_rows,
            "gatk_rows": gatk_rows,
        }

        # The normal VCF cannot see how reads are projected through their
        # best haplotype before ReferenceConfidenceModel computes GQ/PL and
        # block boundaries.  Run the exact same candidate-rich interval in
        # the production default gVCF mode and compare the complete text
        # stream.  Command-line provenance is disabled on both sides, so this
        # is a byte-level GATK 4.6.2.0 oracle rather than a selected-field
        # comparison.
        gatk_gvcf = work / "gatk-contiguous.g.vcf"
        native_gvcf = work / "native-contiguous.g.vcf"
        native_gvcf_json = work / "native-contiguous.g.json"
        gatk_gvcf_command = [
            str(JAVA), "-Xmx1g", "-jar", str(GATK), "HaplotypeCaller",
            "-R", str(REFERENCE), "-I", str(BAM),
            "-L", CONTIGUOUS_INTERVAL[0], "-ERC", "GVCF",
            "-O", str(gatk_gvcf), "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ]
        gatk_gvcf_result = subprocess.run(
            gatk_gvcf_command, text=True, capture_output=True, check=False)
        if gatk_gvcf_result.returncode != 0:
            raise AssertionError(
                f"GATK HC gVCF failed: {gatk_gvcf_result.stderr[-1200:]}")
        native_gvcf_result = subprocess.run([
            str(BINARY), "-R", str(REFERENCE), "-I", str(BAM),
            "-L", CONTIGUOUS_INTERVAL[0], "-ERC", "GVCF",
            "-O", str(native_gvcf), "--output-manifest", str(native_gvcf_json),
            "--threads", "1", "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], text=True, capture_output=True, check=False)
        if native_gvcf_result.returncode != 0:
            raise AssertionError(
                f"native HC gVCF failed: {native_gvcf_result.stderr[-1200:]}")
        assert gatk_gvcf.read_bytes() == native_gvcf.read_bytes(), (
            "dense default gVCF is not byte-identical to GATK")
        assert gvcf_candidate_positions(native_gvcf) == [
            10019967, 10019969, 10020228, 10020229, 10020429, 10020431,
            10020434, 10020438, 10020679, 10020680, 10020681,
        ]
        # This is an RCM-only row immediately after a candidate cluster.  Its
        # nonzero GQ/PL is sensitive to the realigned CIGAR cursor; it used to
        # expose the projection off-by-one while normal candidate VCF rows
        # still matched.
        dense_rcm_row = gvcf_row(native_gvcf, 10020230)
        dense_rcm_sample = dict(zip(
            dense_rcm_row[8].split(":"), dense_rcm_row[9].split(":")))
        assert dense_rcm_row[4] == "<NON_REF>"
        assert dense_rcm_sample == {
            "GT": "0/0", "DP": "12", "GQ": "36", "MIN_DP": "12",
            "PL": "0,36,494",
        }
        dense_gvcf_manifest = json.loads(native_gvcf_json.read_text(encoding="utf-8"))
        assert dense_gvcf_manifest["telemetry"]["rcm_haplotype_realignment_used"] is True
        assert dense_gvcf_manifest["telemetry"]["rcm_realigned_observations"] > 0
        assert dense_gvcf_manifest["compatibility"]["gvcf_candidate_standard_fields"] is True

        contiguous_b = work / "contiguous_b.vcf.gz"
        contiguous_b_json = work / "contiguous_b.json"
        run(contiguous_b, contiguous_b_json, CONTIGUOUS_INTERVAL)
        assert contiguous.read_bytes() == contiguous_b.read_bytes(), (
            "contiguous run is not deterministic")

        def interval_bounds(token: str) -> tuple[int, int]:
            start_s, end_s = token.split(":", 1)[1].split("-", 1)
            return int(start_s), int(end_s)

        def covered_isolated_alleles(span: str) -> dict[str, tuple[str, str]]:
            span_start, span_end = interval_bounds(span)
            covered: dict[str, tuple[str, str]] = {}
            for expected_records in EXPECTED_RECORDS_BY_INTERVAL.values():
                for pos, allele in expected_records.items():
                    position = int(pos)
                    if span_start <= position <= span_end:
                        covered[pos] = allele
            return covered

        subspan_calls: dict[str, int] = {}
        for ordinal, span in enumerate(CONTIGUOUS_SUBSPANS):
            expected = covered_isolated_alleles(span)
            assert len(expected) >= 2, {"span": span, "covered": expected}
            first = work / f"subspan{ordinal}.vcf.gz"
            first_json = work / f"subspan{ordinal}.json"
            run(first, first_json, (span,))
            called = records(first)
            missing = {pos: allele for pos, allele in expected.items()
                       if pos not in called}
            flipped = {pos: (expected[pos], called[pos])
                       for pos in expected
                       if pos in called and called[pos] != expected[pos]}
            assert not missing, {"span": span, "missing": missing, "called": called}
            assert not flipped, {"span": span, "flipped": flipped, "called": called}
            second = work / f"subspan{ordinal}b.vcf.gz"
            second_json = work / f"subspan{ordinal}b.json"
            run(second, second_json, (span,))
            assert first.read_bytes() == second.read_bytes(), (
                f"subspan {span} is not deterministic")
            subspan_calls[span] = len(called)

        print(json.dumps({
            "fixture": "chr20-mnp-real",
            "release": "reads-verified",
            "windows": len(EXPECTED_RECORDS_BY_INTERVAL),
            "records_exact": True,
            "deterministic_bytes": True,
            "contiguous_interval": CONTIGUOUS_INTERVAL[0],
            "contiguous_contains_isolated_alleles": True,
            "contiguous_gatk_allele_set_exact": True,
            "contiguous_pairhmm_haplotype_set_exact": True,
            "contiguous_vcf_header_exact_excluding_command_line": True,
            "contiguous_vcf_rows_exact": True,
            "dense_default_gvcf_byte_identical": True,
            "dense_default_gvcf_candidate_positions_exact": True,
            "dense_default_gvcf_rcm_realign_used": True,
            "dense_default_gvcf_rcm_gq_pl_exact": True,
            "contiguous_pairhmm_haplotype_count": len(native_haplotypes),
            "contiguous_deterministic_bytes": True,
            "contiguous_variant_calls": len(contiguous_called),
            "contiguous_subspans": list(CONTIGUOUS_SUBSPANS),
            "contiguous_subspan_variant_calls": subspan_calls,
            "variant_calls": total_calls,
            "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
