#!/usr/bin/env python3
"""Strict acceptance oracle for the ``--alleles`` forced-allele *emission bypass*
divergence (symbolic ``*`` ALT + ``LowQual`` instead of a dropped record).

Defect under test
-----------------
When a forced ``--alleles`` record is retained in the output although the allele
frequency calculation finds it *implausible*, GATK 4.6.2.0 does **not** drop the
record: ``GenotypingEngine.calculateGenotypes`` guards its emit-threshold null
return with ``&& forcedAlleles.isEmpty()``
(``gatk-source/.../genotyper/GenotypingEngine.java:167-170``) and merely attaches
the ``LowQual`` filter when the phred confidence stays below
``--standard-confidence-for-calling`` (``GenotypingEngine.java:184-186``,
``passesCallThreshold`` at ``GenotypingEngine.java:430-432``).  The record that
survives keeps every allele that
``GenotypingEngine.calculateOutputAlleleSubset`` marked ``toOutput`` -- including
the symbolic spanning-deletion ``*`` allele, which is kept because it *is*
plausible (``GenotypingEngine.java:311-322``) -- so the emitted row is
multi-allelic (``*,A``).

Smallest counterexample (the probe fixture reported in
``fastgatk-native/evidence/2026-09-11-wave0/round-trackb-fix-report.md`` §5.2)
---------------------------------------------------------------------------
A pileup in which six reads carry an assembled 1 bp deletion (left-aligned to
``chr1:697 TC>T``, i.e. deleting ``chr1:698``) **plus** an ``--alleles`` VCF with
two PASS records::

    chr1 698 . C  A . PASS .        <- the forced SNP, on the base the deletion removes
    chr1 700 . CA C . PASS .        <- the D1 forced deletion (equivalent, not equal)

Pinned GATK 4.6.2.0 emits **three** rows::

    chr1 697 . TC T   98.60 .          0/1:30,6:36:99:106,0,961
    chr1 698 . C  *,A 0     LowQual    0/1:30,6,0:36:99:106,0,961,207,1055,1478
    chr1 700 . CA C   144.77 .         0/0:36,0:36:99:0,108,1481

native emits only the ``697`` and ``700`` rows: the forced SNP's own row is
discarded by native's confidence gate (``qual == 0 <
--standard-confidence-for-calling``) instead of being emitted with
``FILTER=LowQual`` and the ``*`` spanning-deletion ALT.

Attribution controls
--------------------
* ``--drop-alleles`` (the ``--alleles`` flag removed) makes **both** sides emit a
  single identical row, proving the ``698`` row exists *only* because of
  ``--alleles`` and is not a discovery/pileup difference.
* the second case keeps the same pileup but omits the forced SNP record, so the
  ``*,A`` row disappears on **both** sides -- proving the forced SNP record (not
  the forced deletion, not the assembled deletion alone) is the trigger.

What this oracle asserts
------------------------
For every case and every mode it runs the pinned GATK 4.6.2.0 jar and the native
binary with identical arguments and requires the native **data rows to be
byte-identical** to GATK's, field by field (header and ``#``/provenance lines are
ignored by construction).  It exits non-zero on any mismatch so it can be
registered in CTest.

``--expect-divergence`` inverts the run into a pure diagnostic report and exits
0; that is how the *current, unfixed* binary can be observed without turning a
test suite red.  The default is the strict assertion.

Out of scope (deliberately not gated here): the D3 ``--alleles`` exit-code
divergence on a record that minimises to an empty ALT (GATK aborts with exit 3,
native exits 0).
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path

REFERENCE_LENGTH = 1500
HOMOPOLYMER_START = 600          # 1-based, 10 x 'A' -> 600..609
TANDEM_START = 900               # 1-based, 8 x 'CAG' -> 900..923
ASSEMBLED_DELETION_POS = 700     # 1-based reference base removed by the 6 alt reads
FORCED_SNP_POS = 698             # 1-based base the (left-aligned) assembled deletion removes
FORCED_DELETION_POS = 700        # 1-based anchor of the D1 forced deletion

# ``feature`` selects the exact ``--alleles`` records of a case (see
# ``case_vcf_rows``); ``mode`` runs each case with the feature injected
# (``alleles``) and with the ``--alleles`` flag removed (``drop-alleles``), which
# is the attribution control: the divergent row must exist only in ``alleles``.
CASES: tuple[dict, ...] = (
    {
        "case": "forced-snp-inside-supported-deletion",
        "interval": "chr1:500-780",
        "feature": "snp+deletion",
        "why": "the trigger: a forced SNP (698 C>A) on the base removed by an "
               "assembled deletion that has real read support (697 TC>T, 6 reads); "
               "GATK keeps the implausible forced allele and emits '*,A' LowQual",
    },
    {
        "case": "forced-deletion-only",
        "interval": "chr1:500-780",
        "feature": "deletion",
        "why": "same pileup without the forced SNP record: the '*,A' row must "
               "disappear on BOTH sides, proving the forced SNP is the trigger",
    },
)

MODES: tuple[tuple[str, bool], ...] = (
    ("alleles", True),
    ("drop-alleles", False),
)


def records(path: Path) -> list[list[str]]:
    """Data rows only: header/provenance lines are ignored by construction."""
    return [line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def build_reference() -> str:
    random.seed(11)
    bases = list("".join(random.choice("ACGT") for _ in range(REFERENCE_LENGTH)))
    bases[HOMOPOLYMER_START - 1:HOMOPOLYMER_START - 1 + 10] = list("A" * 10)
    bases[TANDEM_START - 1:TANDEM_START - 1 + 24] = list("CAG" * 8)
    return "".join(bases)


def write_inputs(work: Path) -> tuple[Path, Path, str]:
    reference_text = build_reference()
    reference = work / "ref.fa"
    reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
    # One FASTA sequence line keeps the .fai geometry deterministic without
    # depending on a local samtools executable.
    (work / "ref.fa.fai").write_text(
        f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
        encoding="utf-8")

    header = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}",
              "@RG\tID:rg\tSM:S1"]
    placed: list[tuple[int, str]] = []
    index = 0

    def add_reads(first_start: int, count: int, length: int, step: int) -> None:
        nonlocal index
        for offset in range(count):
            start = first_start + offset * step
            sequence = reference_text[start - 1:start - 1 + length]
            placed.append((start, f"r{index}\t0\tchr1\t{start}\t60\t{length}M\t*\t0\t0\t"
                                  f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg"))
            index += 1

    # Reference-only coverage over the homopolymer window and over the tandem
    # repeat (the latter only matters to the sibling overlap oracle).
    add_reads(440, 30, 300, 3)
    add_reads(800, 30, 300, 3)

    # Six reads carrying an assembled 1 bp deletion at ASSEMBLED_DELETION_POS.
    # The deletion is left-aligned into the local 'C' run when it is emitted, so
    # the emitted row is chr1:697 TC>T and FORCED_SNP_POS (698) is exactly the
    # base that row removes.
    for offset in range(6):
        start = 560 + offset * 4
        sequence = reference_text[start - 1:start - 1 + 300]
        ref_offset = ASSEMBLED_DELETION_POS - start
        sequence = sequence[:ref_offset] + sequence[ref_offset + 1:]
        placed.append((start, f"d{offset}\t0\tchr1\t{start}\t60\t299M\t*\t0\t0\t{sequence}\t"
                              f"{'I' * len(sequence)}\tRG:Z:rg"))
    lines = header + [line for _, line in sorted(placed, key=lambda item: item[0])]
    sam = work / "input.sam"
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return reference, sam, reference_text


def case_vcf_rows(case: dict, reference_text: str) -> list[tuple[str, int, str, str]]:
    """The exact feature VCF records of a case, in ascending VCF order."""
    rows: list[tuple[str, int, str, str]] = []
    if case["feature"] == "snp+deletion":
        rows.append(("chr1", FORCED_SNP_POS, reference_text[FORCED_SNP_POS - 1], "A"))
    deletion_ref = reference_text[FORCED_DELETION_POS - 1:FORCED_DELETION_POS + 1]
    rows.append(("chr1", FORCED_DELETION_POS, deletion_ref, deletion_ref[0]))
    return rows


def write_feature(work: Path, name: str, mode: str,
                  rows: list[tuple[str, int, str, str]], root: Path) -> Path:
    path = work / f"{name}.{mode}.vcf"
    body = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={REFERENCE_LENGTH}>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, position, ref, alt in rows:
        body.append(f"{chrom}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.")
    path.write_text("\n".join(body) + "\n", encoding="utf-8")
    bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
    tabix = root / "third_party/htslib-build/htslib-src/tabix"
    subprocess.run([str(bgzip), "-f", str(path)], check=True, stdout=subprocess.DEVNULL)
    compressed = path.with_suffix(path.suffix + ".gz")
    subprocess.run([str(tabix), "-f", "-p", "vcf", str(compressed)], check=True,
                   stdout=subprocess.DEVNULL)
    return compressed


def row_field_differences(gatk_row: list[str], native_row: list[str]) -> list[str]:
    """Field-level diff of one VCF row (INFO sub-fields and FORMAT sub-fields)."""
    differences: list[str] = []
    # VCF columns: 3 REF, 4 ALT, 5 QUAL, 6 FILTER.
    for index in (3, 4, 5, 6):
        if index < len(gatk_row) and index < len(native_row) and gatk_row[index] != native_row[index]:
            label = {3: "REF", 4: "ALT", 5: "QUAL", 6: "FILTER"}[index]
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


def compare(case: str, mode: str, gatk_rows: list[list[str]],
            native_rows: list[list[str]]) -> dict:
    positions = divergent_positions(gatk_rows, native_rows)
    first_row_diff: list[str] = []
    for index in range(min(len(gatk_rows), len(native_rows))):
        if gatk_rows[index] != native_rows[index]:
            first_row_diff = row_field_differences(gatk_rows[index], native_rows[index])
            break
    return {
        "case": case,
        "mode": mode,
        "gatk_rows": len(gatk_rows),
        "native_rows": len(native_rows),
        "data_rows_byte_identical": gatk_rows == native_rows,
        "divergent_positions": positions,
        "first_divergent_row_field_diff": first_row_diff,
        "gatk_rows_all": gatk_rows,
        "native_rows_all": native_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the --alleles forced-allele emission bypass "
                    "(LowQual + symbolic '*' ALT) against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_HC_BINARY"),
        help="native HaplotypeCaller binary "
             "(default: $FASTGATK_HC_BINARY or fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument(
        "--case", action="append", default=None,
        help=f"run only the named case(s); default is all of "
             f"{[case['case'] for case in CASES]}")
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
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")

    strict_mode = not arguments.expect_divergence
    results: list[dict] = []
    run_errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-forced-alleles-emission-oracle-") as directory:
        work = Path(directory)
        reference, sam, reference_text = write_inputs(work)
        bam = work / "input.bam"
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam),
                        "-O", str(bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)

        for case in selected:
            feature_rows = case_vcf_rows(case, reference_text)
            for mode, with_alleles in MODES:
                feature = write_feature(work, case["case"], mode, feature_rows, root)
                allele_args = ["--alleles", str(feature)] if with_alleles else []
                common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
                          *allele_args, "--min-pruning", "1",
                          "--create-output-variant-index", "false",
                          "--add-output-vcf-command-line", "false"]
                java_vcf = work / f"gatk.{case['case']}.{mode}.vcf"
                native_vcf = work / f"native.{case['case']}.{mode}.vcf"
                java_run = subprocess.run(
                    [java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                     "-O", str(java_vcf)],
                    capture_output=True, text=True, timeout=300)
                native_run = subprocess.run(
                    [str(native), *common, "--threads", str(arguments.threads),
                     "-O", str(native_vcf)],
                    capture_output=True, text=True, timeout=300)
                label = f"{case['case']} [{mode}]"
                if java_run.returncode != 0 or native_run.returncode != 0:
                    run_errors.append(
                        f"{label}: java_exit={java_run.returncode} "
                        f"native_exit={native_run.returncode} "
                        f"java_stderr={java_run.stderr[-300:]!r} "
                        f"native_stderr={native_run.stderr[-300:]!r}")
                    continue
                comparison = compare(case["case"], mode, records(java_vcf),
                                     records(native_vcf))
                comparison["why"] = case["why"]
                comparison["feature_rows"] = [
                    f"{chrom}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t."
                    for chrom, position, ref, alt in feature_rows]
                results.append(comparison)

    violations: list[str] = list(run_errors)
    for result in results:
        if not result["data_rows_byte_identical"]:
            violations.append(
                f"{result['case']} [{result['mode']}]: data rows differ at "
                f"{result['divergent_positions'][:10]} "
                f"(gatk_rows={result['gatk_rows']} native_rows={result['native_rows']})")
        if result["mode"] == "drop-alleles":
            injected = next((other for other in results
                             if other["case"] == result["case"] and other["mode"] == "alleles"),
                            None)
            # Attribution control: without --alleles both sides must still be
            # byte-identical (they are gated above), and GATK must emit strictly
            # fewer rows than with --alleles, which proves every divergent row is
            # forced-allele driven rather than a pileup/discovery difference.
            if injected is not None and result["gatk_rows"] >= injected["gatk_rows"]:
                violations.append(
                    f"{result['case']}: --drop-alleles GATK emitted "
                    f"{result['gatk_rows']} row(s), not fewer than the "
                    f"{injected['gatk_rows']} row(s) with --alleles; fixture no "
                    f"longer isolates the forced-allele emission path")
    # The trigger control: dropping the forced SNP record must remove every
    # symbolic-'*' forced row on both sides.  If GATK still emits one, the case
    # no longer isolates the trigger and the oracle is measuring something else.
    with_snp = next((result for result in results
                     if result["case"] == "forced-snp-inside-supported-deletion"
                     and result["mode"] == "alleles"), None)
    without_snp = next((result for result in results
                        if result["case"] == "forced-deletion-only"
                        and result["mode"] == "alleles"), None)
    if with_snp is not None and without_snp is not None:
        symbolic_in_with = [row for row in with_snp["gatk_rows_all"]
                            if len(row) > 4 and "*" in row[4].split(",")]
        symbolic_in_without = [row for row in without_snp["gatk_rows_all"]
                               if len(row) > 4 and "*" in row[4].split(",")]
        if not symbolic_in_with:
            violations.append(
                "forced-snp-inside-supported-deletion [alleles]: pinned GATK did "
                "not emit any '*,...' row, so this fixture no longer reproduces "
                "the forced-allele emission bypass")
        if symbolic_in_without:
            violations.append(
                "forced-deletion-only [alleles]: pinned GATK emitted a '*,...' row "
                "without the forced SNP record, so the forced SNP is not the trigger")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native "
          f"({native}); forced --alleles record that stays implausible must be "
          f"emitted (LowQual, symbolic '*' ALT) rather than dropped")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    for result in results:
        print(f"[{result['case']} / {result['mode']}] "
              f"gatk_rows={result['gatk_rows']} native_rows={result['native_rows']} "
              f"data_rows_byte_identical={result['data_rows_byte_identical']} "
              f"divergent_positions={result['divergent_positions'][:10]}")
        print(f"    feature VCF rows: {'; '.join(result['feature_rows'])}")
        for row in result["gatk_rows_all"]:
            print(f"    GATK   {chr(9).join(row)}")
        for row in result["native_rows_all"]:
            print(f"    NATIVE {chr(9).join(row)}")
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

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "--alleles forced-allele emission bypass (LowQual + symbolic '*' ALT)",
        "binary": str(native),
        "acceptance_criterion": (
            "native data rows (header/provenance ignored) must be byte-identical "
            "to pinned GATK 4.6.2.0 for every case, in the --alleles run and in its "
            "--drop-alleles attribution control"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "modes": [mode for mode, _ in MODES],
        "violations": violations,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
