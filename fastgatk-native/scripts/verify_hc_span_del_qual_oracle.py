#!/usr/bin/env python3
"""Strict acceptance oracle for the spanning-deletion QUAL/QD divergence.

Observation under test (reported by a prior probe whose fixture no longer
existed, see ``fastgatk-native/evidence/2026-09-11-wave0/round-forced-alleles-emission.md``
section 7.3): **without** ``--alleles``, a pileup containing

* 30 reference-only reads,
* 12 reads carrying a real ``chr1:698 C>A`` SNP,
* 6 reads carrying an assembled 1 bp deletion inside the ``CCC`` run at
  ``chr1:698-700`` (emitted as the left-aligned ``chr1:697 TC>T``),

makes both callers emit a *biallelic* ``chr1 698 . C A`` row whose ALT, every
INFO field except ``QUAL``/``QD``, and every FORMAT field are identical, but
with ``QUAL`` = 282.04 (QD 5.88) for pinned GATK 4.6.2.0 and 291.07 (QD 6.06)
for native.

Hypothesis (SPECULATION until measured -- this oracle exists to test it, not to
assume it): native's ordinary diploid path forces a concrete-only allele set
(``hc_call.cpp`` ``include_spanning_deletion``), so its allele-frequency
calculation excludes the symbolic ``*`` allele that GATK keeps in
``GenotypingEngine.calculateGenotypes``' ``AFresult`` (``GenotypingEngine.java``
lines 134-165), changing the no-variant posterior and therefore ``QUAL``.

Measured outcome (2026-09-11, this repo): the divergence **does** reproduce, and
the hypothesis above is **partially refuted**.  Native does not drop ``*`` before
the AF calculation -- its ``derive_pairhmm_spanning_deletion_pl`` path builds a
three-allele ``[REF, concrete ALT, *]`` genotype-likelihood row and feeds it to
``calculate_allele_frequency_kokkos``.  The residual error was the ``*`` allele's
Dirichlet pseudocount: the native code used ``indel_heterozygosity``
unconditionally, while GATK's ``AlleleFrequencyCalculator`` derives each allele's
pseudocount from its length (``a.length() == refLength ? snpPseudocount :
indelPseudocount``, ``AlleleFrequencyCalculator.java:175-176``) and htsjdk's
``Allele.SPAN_DEL.length()`` is 1, so ``*`` takes the *SNP* prior whenever the
record's REF allele is 1 bp.  Fixing that (``calling_pipeline.cpp``, the
``spanning_pseudocount`` initializer) changes native's QUAL at ``chr1:698`` from
291.07 to GATK's 282.04.

Which rows are compared
-----------------------
Both callers run on the *same* reference/BAM with identical arguments
(``--min-pruning 1``, no ``--alleles``); the native side additionally gets
``--threads``.  Every data row is compared field by field (REF/ALT/QUAL/FILTER,
each INFO sub-field, FORMAT keys and each FORMAT value); ``#`` provenance lines
and the header are ignored by construction.  ``QUAL``/``QD`` are additionally
printed explicitly for every row of every case.

Cases are chosen so the result is attributable rather than anecdotal:

* ``reported-frameshift-span-del`` -- faithful reconstruction of the reported
  fixture; reproduces the reported QUAL pair (282.04 / 291.07) exactly.
* ``no-span-del-control`` -- the same pileup with the 6 deletion reads removed;
  proves the discrepancy needs the deletion (parity expected).
* ``clean-cigar-del-*`` -- the deletion re-expressed as a real ``138M1D161M``
  CIGAR instead of sequence removal, at several deletion-read depths: does the
  divergence survive a non-pathological deletion encoding?
* ``reported-*`` robustness probes: different ``-L`` windows,
  ``--max-mnp-distance`` 0/1/2, higher reference depth.

``OUT_OF_SCOPE_CASES`` holds probes whose rows differ for a *separately
tracked* defect -- at ``--sample-ploidy 3`` the reported pileup makes GATK emit
``ALT=*,A`` where native emits ``ALT=A``, i.e. the ordinary-VCF symbolic
spanning-deletion ALT emission gap.  Those probes are still run and printed
(QUAL/QD included) but their mismatches are reported as ``out_of_scope_notes``
rather than violations, so this oracle neither gates nor hides that gap.

Exit status: 0 when native's data rows are byte-identical to GATK's in every
gated case (``CASES``), non-zero otherwise (so it can be registered in CTest).
With ``--expect-divergence`` the run becomes a pure diagnostic and always exits
0.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path
import oracle_guard

REFERENCE_LENGTH = 1500
HOMOPOLYMER_START = 600          # 1-based, 10 x 'A' -> 600..609
TANDEM_START = 900               # 1-based, 8 x 'CAG' -> 900..923
SNP_POSITION = 698               # 1-based; reference base is 'C' (seed 11)
DELETION_RUN = (698, 699, 700)   # the 'CCC' run the 1 bp deletion slides in
DEFAULT_INTERVAL = "chr1:500-780"

# Read layouts.  ``ref``/``snp``/``del`` describe groups of reads:
#   ref: reference-matching reads, ``count`` of them from ``start``, step ``step``
#   snp: same, but the base at ``position`` is replaced by ``alt``
#   del: ``mode='frameshift'`` removes ``drop_base`` from the sequence and keeps
#        ``length - 1`` M (the encoding used by
#        verify_hc_alleles_overlap_gate_oracle.py, which the reported probe
#        reused); ``mode='cigar'`` emits a real ``<n>M1D<m>M`` CIGAR that deletes
#        ``drop_base``.
CASES: tuple[dict, ...] = (
    {
        "case": "reported-frameshift-span-del",
        "why": "faithful reconstruction of the reported fixture (30 ref reads, "
               "12 SNP reads, 6 deletion reads encoded as sequence removal)",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": (),
    },
    {
        "case": "no-span-del-control",
        "why": "same pileup without the 6 deletion reads: the discrepancy must "
               "require the spanning deletion, otherwise the fixture is not "
               "isolating it",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": None},
        "extra": (),
    },
    {
        "case": "clean-cigar-del-6",
        "why": "same pileup, deletion re-expressed as a real 138M1D161M CIGAR",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "cigar", "drop_base": 698}},
        "extra": (),
    },
    {
        "case": "clean-cigar-del-18",
        "why": "deeper real-CIGAR deletion so the deletion is actually called",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 18, "start": 560, "step": 2, "length": 300,
                           "mode": "cigar", "drop_base": 698}},
        "extra": (),
    },
    {
        "case": "frameshift-del-18",
        "why": "deeper frameshift-encoded deletion (depth robustness)",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 18, "start": 560, "step": 2, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": (),
    },
    {
        "case": "reported-narrow-window",
        "why": "reported pileup restricted to -L chr1:600-760 (window robustness)",
        "interval": "chr1:600-760",
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": (),
    },
    {
        "case": "reported-max-mnp-distance-0",
        "why": "reported pileup with --max-mnp-distance 0",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": ("--max-mnp-distance", "0"),
    },
    {
        "case": "reported-max-mnp-distance-2",
        "why": "reported pileup with --max-mnp-distance 2",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": ("--max-mnp-distance", "2"),
    },
    {
        "case": "reported-deep-ref-60",
        "why": "reported pileup with 60 reference-only reads (depth robustness)",
        "interval": DEFAULT_INTERVAL,
        "layout": {"ref": {"count": 60, "start": 440, "step": 2, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": (),
    },
)

# Probes whose rows differ for a *different*, separately tracked reason and are
# therefore reported but not gated here.  At ploidy 3 the reported pileup makes
# GATK emit `ALT=*,A` while native emits `ALT=A`: that is the ordinary-VCF
# symbolic spanning-deletion ALT emission gap named in
# ``fastgatk-native/evidence/2026-09-11-wave0/round-forced-alleles-emission.md``
# section 7 / recommendation item 3, which this oracle is explicitly not
# allowed to attempt.  Their QUAL/QD values are still printed so the QUAL
# parity (or lack of it) stays visible on every run.
OUT_OF_SCOPE_CASES: tuple[dict, ...] = (
    {
        "case": "reported-ploidy-3",
        "why": "reported pileup genotyped at --sample-ploidy 3; measured: QUAL/QD "
               "already agree (358.64 / 7.47 on both sides), only the ALT set "
               "differs (GATK `*,A` vs native `A`) -- the separately tracked "
               "symbolic-ALT emission gap, out of scope for this oracle",
        "interval": DEFAULT_INTERVAL,
        "ploidy": 3,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
        "extra": (),
    },
)

ALL_CASES: tuple[dict, ...] = CASES + OUT_OF_SCOPE_CASES


def build_reference() -> str:
    random.seed(11)
    bases = list("".join(random.choice("ACGT") for _ in range(REFERENCE_LENGTH)))
    bases[HOMOPOLYMER_START - 1:HOMOPOLYMER_START - 1 + 10] = list("A" * 10)
    bases[TANDEM_START - 1:TANDEM_START - 1 + 24] = list("CAG" * 8)
    return "".join(bases)


def build_reads(layout: dict, reference_text: str) -> list[tuple[int, str, str, str]]:
    """(start, cigar, sequence, tag) for every read of a case, sorted by start."""
    placed: list[tuple[int, str, str, str]] = []
    group = layout["ref"]
    for offset in range(group["count"]):
        start = group["start"] + offset * group["step"]
        placed.append((start, f"{group['length']}M",
                       reference_text[start - 1:start - 1 + group["length"]], "r"))

    group = layout.get("snp")
    if group:
        for offset in range(group["count"]):
            start = group["start"] + offset * group["step"]
            sequence = list(reference_text[start - 1:start - 1 + group["length"]])
            sequence[SNP_POSITION - start] = group.get("alt", "A")
            placed.append((start, f"{group['length']}M", "".join(sequence), "s"))

    group = layout.get("del")
    if group:
        length, drop = group["length"], group["drop_base"]
        for offset in range(group["count"]):
            start = group["start"] + offset * group["step"]
            if group["mode"] == "frameshift":
                # Sequence-level removal with a (length - 1)M CIGAR, exactly the
                # encoding used by the reported probe's fixture.
                sequence = (reference_text[start - 1:drop - 1]
                            + reference_text[drop:start - 1 + length])
                cigar = f"{length - 1}M"
            else:
                # A real deletion: <lead>M1D<trail>M.
                lead = drop - start
                sequence = (reference_text[start - 1:drop - 1]
                            + reference_text[drop:start - 1 + length])
                cigar = f"{lead}M1D{length - lead - 1}M"
            placed.append((start, cigar, sequence, "d"))
    return placed


def write_sam(work: Path, name: str, layout: dict, reference_text: str) -> Path:
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             f"@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}", "@RG\tID:rg\tSM:S1"]
    reads = build_reads(layout, reference_text)
    for index, (start, cigar, sequence, tag) in enumerate(
            sorted(reads, key=lambda item: item[0])):
        lines.append(f"{tag}{index}\t0\tchr1\t{start}\t60\t{cigar}\t*\t0\t0\t"
                     f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg")
    path = work / f"{name}.sam"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def records(path: Path) -> list[list[str]]:
    """Data rows only: header/provenance lines are ignored by construction."""
    if not path.is_file():
        return []
    return [line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def row_field_differences(gatk_row: list[str], native_row: list[str]) -> list[str]:
    """Field-level diff of one VCF row (column 5 is QUAL, column 6 is FILTER)."""
    differences: list[str] = []
    for index, label in ((3, "REF"), (4, "ALT"), (5, "QUAL"), (6, "FILTER")):
        if index < len(gatk_row) and index < len(native_row) \
                and gatk_row[index] != native_row[index]:
            differences.append(f"{label} {gatk_row[index]!r} -> {native_row[index]!r}")
    gatk_info = {} if len(gatk_row) < 8 else dict(
        field.split("=", 1) if "=" in field else (field, "")
        for field in gatk_row[7].split(";"))
    native_info = {} if len(native_row) < 8 else dict(
        field.split("=", 1) if "=" in field else (field, "")
        for field in native_row[7].split(";"))
    for key in sorted(set(gatk_info) | set(native_info)):
        if gatk_info.get(key) != native_info.get(key):
            differences.append(
                f"INFO {key} {gatk_info.get(key)!r} -> {native_info.get(key)!r}")
    if len(gatk_row) > 9 and len(native_row) > 9:
        if gatk_row[8] != native_row[8]:
            differences.append(f"FORMAT keys {gatk_row[8]!r} -> {native_row[8]!r}")
        keys = gatk_row[8].split(":")
        gatk_values = gatk_row[9].split(":")
        native_values = native_row[9].split(":")
        for index, key in enumerate(keys):
            gatk_value = gatk_values[index] if index < len(gatk_values) else None
            native_value = native_values[index] if index < len(native_values) else None
            if gatk_value != native_value:
                differences.append(f"FORMAT {key} {gatk_value!r} -> {native_value!r}")
    return differences


def qual_qd(row: list[str]) -> dict:
    """POS/REF/ALT/QUAL plus the QD INFO sub-field of one data row."""
    info = {} if len(row) < 8 else dict(
        field.split("=", 1) if "=" in field else (field, "")
        for field in row[7].split(";"))
    return {"POS": row[1] if len(row) > 1 else None,
            "REF": row[3] if len(row) > 3 else None,
            "ALT": row[4] if len(row) > 4 else None,
            "QUAL": row[5] if len(row) > 5 else None,
            "QD": info.get("QD")}


def divergent_positions(gatk_rows: list[list[str]],
                        native_rows: list[list[str]]) -> list[str]:
    positions: list[str] = []
    for index in range(max(len(gatk_rows), len(native_rows))):
        gatk_row = gatk_rows[index] if index < len(gatk_rows) else None
        native_row = native_rows[index] if index < len(native_rows) else None
        if gatk_row is None:
            positions.append(f"extra-in-native@{index}:POS={native_row[1]}")
        elif native_row is None:
            positions.append(f"missing-in-native@{index}:POS={gatk_row[1]}")
        elif gatk_row != native_row:
            positions.append(native_row[1] if native_row[1] == gatk_row[1]
                             else f"{gatk_row[1]}/{native_row[1]}")
    return positions


def compare(case: str, gatk_rows: list[list[str]], native_rows: list[list[str]]) -> dict:
    first_row_diff: list[str] = []
    for index in range(min(len(gatk_rows), len(native_rows))):
        if gatk_rows[index] != native_rows[index]:
            first_row_diff = row_field_differences(gatk_rows[index], native_rows[index])
            break
    return {
        "case": case,
        "gatk_rows": len(gatk_rows),
        "native_rows": len(native_rows),
        "data_rows_byte_identical": gatk_rows == native_rows,
        "divergent_positions": divergent_positions(gatk_rows, native_rows),
        "first_divergent_row_field_diff": first_row_diff,
        "gatk_qual_qd": [qual_qd(row) for row in gatk_rows],
        "native_qual_qd": [qual_qd(row) for row in native_rows],
        "gatk_rows_all": gatk_rows,
        "native_rows_all": native_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the HaplotypeCaller spanning-deletion "
                    "QUAL/QD divergence (no --alleles) against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_HC_BINARY"),
        help="native HaplotypeCaller binary (default: $FASTGATK_HC_BINARY or "
             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help=f"run only the named case(s); default is all of "
                             f"{[case['case'] for case in ALL_CASES]}")
    parser.add_argument("--threads", type=int, default=2,
                        help="--threads value for the native run (default 2)")
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(arguments.native) if arguments.native else (
        root / "fastgatk-native/build/fastgatk-hc-call")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))

    assets = [native, Path(java), gatk]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified('verify_hc_span_del_qual_oracle.py', Path(java), gatk)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}))
        return 0

    selected = [case for case in ALL_CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    gated = {case["case"] for case in CASES}

    strict_mode = not arguments.expect_divergence
    results: list[dict] = []
    run_errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-span-del-oracle-") as directory:
        work = Path(directory)
        reference_text = build_reference()
        reference = work / "ref.fa"
        reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
        # One FASTA sequence line keeps the .fai geometry deterministic without
        # depending on a local samtools executable.
        (work / "ref.fa.fai").write_text(
            f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
            encoding="utf-8")
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)

        for case in selected:
            sam = write_sam(work, case["case"], case["layout"], reference_text)
            bam = work / f"{case['case']}.bam"
            subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam",
                            "-I", str(sam), "-O", str(bam), "-SO", "coordinate",
                            "--CREATE_INDEX", "true"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=300)
            common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
                      "--min-pruning", "1",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false"]
            if case.get("ploidy"):
                common += ["--sample-ploidy", str(case["ploidy"])]
            common += list(case["extra"])
            gatk_vcf = work / f"gatk.{case['case']}.vcf"
            native_vcf = work / f"native.{case['case']}.vcf"
            gatk_run = subprocess.run(
                [java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                 "-O", str(gatk_vcf)],
                capture_output=True, text=True, timeout=900)
            native_run = subprocess.run(
                [str(native), *common, "--threads", str(arguments.threads),
                 "-O", str(native_vcf)],
                capture_output=True, text=True, timeout=900)
            if gatk_run.returncode != 0 or native_run.returncode != 0:
                run_errors.append(
                    f"{case['case']}: java_exit={gatk_run.returncode} "
                    f"native_exit={native_run.returncode} "
                    f"java_stderr={gatk_run.stderr[-300:]!r} "
                    f"native_stderr={native_run.stderr[-300:]!r}")
                continue
            comparison = compare(case["case"], records(gatk_vcf), records(native_vcf))
            comparison["why"] = case["why"]
            comparison["interval"] = case["interval"]
            comparison["extra_args"] = list(case["extra"])
            comparison["ploidy"] = case.get("ploidy", 2)
            comparison["gated"] = case["case"] in gated
            results.append(comparison)

    violations: list[str] = list(run_errors)
    out_of_scope_notes: list[str] = []
    for result in results:
        if result["data_rows_byte_identical"]:
            continue
        note = (f"{result['case']}: data rows differ at "
                f"{result['divergent_positions'][:10]} "
                f"(gatk_rows={result['gatk_rows']} native_rows={result['native_rows']})")
        if result["gated"]:
            violations.append(note)
        else:
            out_of_scope_notes.append(note)

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# pileup with a chr1:698 C>A SNP plus one-bp-deletion reads in the "
          "chr1:698-700 CCC run, without --alleles")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated cases (mismatch => exit 1): {[case['case'] for case in CASES]}")
    print(f"# out-of-scope probes (reported, not gated): "
          f"{[case['case'] for case in OUT_OF_SCOPE_CASES]}")
    for result in results:
        print(f"[{result['case']}] interval={result['interval']} "
              f"ploidy={result['ploidy']} extra={result['extra_args']} "
              f"gated={result['gated']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_rows={result['gatk_rows']} native_rows={result['native_rows']} "
              f"data_rows_byte_identical={result['data_rows_byte_identical']} "
              f"divergent_positions={result['divergent_positions'][:10]}")
        for row in result["gatk_rows_all"]:
            print(f"    GATK   {chr(9).join(row)}")
        for row in result["native_rows_all"]:
            print(f"    NATIVE {chr(9).join(row)}")
        for side in ("gatk", "native"):
            for entry in result[f"{side}_qual_qd"]:
                print(f"    {side.upper():6s} QUAL/QD POS={entry['POS']} "
                      f"REF={entry['REF']} ALT={entry['ALT']} "
                      f"QUAL={entry['QUAL']} QD={entry['QD']}")
        if result["first_divergent_row_field_diff"]:
            print("    first divergent row field diff (GATK -> NATIVE): "
                  + "; ".join(result["first_divergent_row_field_diff"]))
    if run_errors:
        print(f"# {len(run_errors)} run error(s):")
        for error in run_errors:
            print(f"#   - {error}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")
    if out_of_scope_notes:
        print(f"# {len(out_of_scope_notes)} out-of-scope note(s) (not counted as "
              f"violations):")
        for note in out_of_scope_notes:
            print(f"#   - {note}")

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller spanning-deletion QUAL/QD parity (no --alleles)",
        "binary": str(native),
        "acceptance_criterion": (
            "native data rows (header/provenance ignored) must be byte-identical "
            "to pinned GATK 4.6.2.0 for every gated case, QUAL and QD included"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": [case["case"] for case in CASES],
        "out_of_scope_cases": [case["case"] for case in OUT_OF_SCOPE_CASES],
        "violations": violations,
        "out_of_scope_notes": out_of_scope_notes,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
