#!/usr/bin/env python3
"""Strict oracle for the gVCF reference-block PL/GQ divergence under forced alleles.

Observation under test
----------------------
``verify_hc_alleles_deep_limits.py`` reports ``gvcf-max-alt-alleles-1`` as
``diverged`` with equal row counts (42 vs 42).  Its own summary only names the
first differing row; the field-level differences are **three reference-block
rows** plus the already-tracked concrete-indel ``END`` INFO (see
``verify_hc_gvcf_indel_end_gatk_oracle.py``)::

    row  6  GATK   chr1 621 . A <NON_REF> . . END=639  GT:DP:GQ:MIN_DP:PL  0/0:36:99:36:0,108,1612
            NATIVE chr1 621 . A <NON_REF> . . END=639  GT:DP:GQ:MIN_DP:PL  0/0:36:90:36:0,90,1350
    row  8  GATK   chr1 641 ... 0,108,1612      NATIVE ... 0,90,1350
    row 10  GATK   chr1 661 ... 0,108,1612      NATIVE ... 0,90,1350

i.e. both ``GQ`` and ``PL`` differ on every reference block that sits between
the forced SNPs; each side's ``PL`` is an internally consistent vector
(``0,108,1612`` = the 36-read Ref-vs-Any SNP model; ``0,90,1350`` = the n=30
indel cache row), so this is a *model-selection* difference, not a rounding one.
The ordinary (non-gVCF) ``--max-alternate-alleles 1`` case passes; only the gVCF
reference-confidence blocks diverge.  The same divergence appears **without**
``--max-alternate-alleles`` (case ``forced-snps-gvcf`` below), so the forced-alt
budget is not what triggers it.

GATK's rule for the affected quantity
------------------------------------
``ReferenceConfidenceModel.doIndelRefConfCalc`` keeps the *least* confident of
the SNP Ref-vs-Any model and an independent indel model keyed on the number of
indel-informative reads in the pileup:

    // gatk-source/.../tools/walkers/haplotypecaller/ReferenceConfidenceModel.java:301-317
    final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(
            homRefCalc.getGenotypeLikelihoodsCappedByHomRefLikelihood());
    final int nIndelInformativeReads = calcNReadsWithNoPlausibleIndelsReads(
            pileup, refOffset, ref, indelInformativeDepthIndelSize);
    final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy, nIndelInformativeReads);
    final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);
    homRefCalc.finalPhredScaledGenotypeLikelihoods = leastConfidenceGLs.getAsPLs();

``getGLwithWorstGQ(gl1, gl2)`` returns ``gl1`` (the indel model) only when its
hom-ref log10 GQ is the *larger* (i.e. the less confident) of the two
(``ReferenceConfidenceModel.java:331-342``).  At depth 36 the two models are
within 0.001 of each other in that comparison, so the emitted vector is decided
entirely by ``nIndelInformativeReads``: ``30`` selects the indel row
``0,90,1350`` (GATK's own value in this fixture *without* forced alleles and in
the no-deletion-read control), ``36`` selects the SNP row ``0,108,1612``.

Measured localization (this round, see ``.diag/round-deeplimits-gvcf.md``)
------------------------------------------------------------------------
The divergence is confined to ``readHasNoPlausibleIdealsOfSize`` /
``traverseEndOfReadForIndelMismatches`` (Java) versus
``aligned_view_is_indel_informative`` / ``read_is_indel_informative`` (native,
``calling_pipeline.cpp:2403-2496``) as they are reached from the *event-trimmed*
reference-confidence segment (``append_segment``, ``calling_pipeline.cpp:13543-13561``):

* ``--indel-size-to-eliminate-in-ref-model 0`` makes the plausible-indel search
  vacuous in both implementations and the divergence disappears exactly (case
  ``forced-snps-gvcf-indel-size-0``);
* removing the six 1 bp-frameshift deletion reads removes it too
  (case ``forced-snps-gvcf-nodeletion``);
* re-encoding those reads as a real one-base CIGAR deletion removes it as well;
* in the same fixture *without* ``--alleles`` both callers agree
  (``0,90,1350`` everywhere), so the trigger is the forced allele's effect on
  which reference window / pileup the reference-confidence model sees, not the
  alleles themselves.

What this oracle gates
----------------------
Every data row of every gated case must be byte-identical to pinned GATK
4.6.2.0; no canonicalization is applied.  The two controls are gated alongside
the two divergent cases so a "fix" that simply forces one of the two PL vectors
cannot pass.  This gate also covers the concrete-indel ``END`` INFO, which is
tracked and gated separately by ``verify_hc_gvcf_indel_end_gatk_oracle.py``;
runs made before that fix landed report it as an additional violation.

Exit status: 0 when every gated case's data rows are byte-identical, non-zero
otherwise (so it can be registered in CTest once the divergence is fixed).
``--expect-divergence`` turns the run into a pure diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_hc_alleles_deep_boundary import (  # noqa: E402
    ASSEMBLED_DELETION_POS, HOMOPOLYMER_INTERVAL, REFERENCE_LENGTH, build_reference,
    records, write_feature,
)

SNP_POSITIONS = (620, 640, 660)

CASES: tuple[dict, ...] = (
    {
        "case": "forced-snps-gvcf",
        "why": "three forced SNPs (chr1:620/640/660) under -ERC GVCF: the "
               "reported shape without the alt budget.  GATK's reference blocks "
               "between the forced SNPs carry the SNP Ref-vs-Any vector "
               "0,108,1612 / GQ 99, native carries the n=30 indel vector "
               "0,90,1350 / GQ 90",
        "alleles": True, "deletion_reads": True, "extra": ("-ERC", "GVCF"),
    },
    {
        "case": "forced-snps-gvcf-max-alt-1",
        "why": "the exact verify_hc_alleles_deep_limits.py case "
               "gvcf-max-alt-alleles-1: same fixture plus "
               "--max-alternate-alleles 1",
        "alleles": True, "deletion_reads": True,
        "extra": ("-ERC", "GVCF", "--max-alternate-alleles", "1"),
    },
    {
        "case": "forced-snps-gvcf-indel-size-0",
        "why": "CONTROL: same command line plus "
               "--indel-size-to-eliminate-in-ref-model 0, which makes the "
               "plausible-indel search vacuous in both implementations; the "
               "reference blocks must then be byte-identical (they are, and both "
               "sides move to the SNP vector)",
        "alleles": True, "deletion_reads": True,
        "extra": ("-ERC", "GVCF", "--indel-size-to-eliminate-in-ref-model", "0"),
    },
    {
        "case": "forced-snps-gvcf-nodeletion",
        "why": "CONTROL: same command line on the same reference with the six "
               "1 bp-frameshift deletion reads removed; the reference blocks must "
               "be byte-identical (they are), which localizes the divergence to "
               "those reads' indel-informativeness classification",
        "alleles": True, "deletion_reads": False, "extra": ("-ERC", "GVCF"),
    },
)


def build_inputs(work: Path, reference_text: str, with_deletion_reads: bool):
    """Reference + SAM for the shared deep-boundary fixture.

    Mirrors ``verify_hc_alleles_deep_boundary.write_inputs`` exactly (same
    reference text, same 30+30 reference reads, same six 1 bp-frameshift
    deletion reads at chr1:560+4k) but lets the deletion-read group be dropped
    for the control case.
    """
    reference = work / "ref.fa"
    reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
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

    add_reads(440, 30, 300, 3)
    add_reads(800, 30, 300, 3)
    if with_deletion_reads:
        for offset in range(6):
            start = 560 + offset * 4
            sequence = reference_text[start - 1:start - 1 + 300]
            ref_offset = ASSEMBLED_DELETION_POS - start
            sequence = sequence[:ref_offset] + sequence[ref_offset + 1:]
            placed.append((start, f"d{offset}\t0\tchr1\t{start}\t60\t299M\t*\t0\t0\t{sequence}\t"
                                  f"{'I' * len(sequence)}\tRG:Z:rg"))
    lines = header + [line for _, line in sorted(placed, key=lambda item: item[0])]
    sam = work / f"input.{'del' if with_deletion_reads else 'nodel'}.sam"
    sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return reference, sam


def info_fields(row: list[str]) -> dict[str, str]:
    if len(row) < 8:
        return {}
    return dict(field.split("=", 1) if "=" in field else (field, "")
                for field in row[7].split(";"))


def format_fields(row: list[str]) -> dict[str, str]:
    if len(row) < 10:
        return {}
    keys = row[8].split(":")
    values = row[9].split(":")
    return {key: values[index] if index < len(values) else None
            for index, key in enumerate(keys)}


def row_field_differences(gatk_row: list[str], native_row: list[str]) -> list[str]:
    differences: list[str] = []
    for index, label in ((3, "REF"), (4, "ALT"), (5, "QUAL"), (6, "FILTER")):
        if index < len(gatk_row) and index < len(native_row) \
                and gatk_row[index] != native_row[index]:
            differences.append(f"{label} {gatk_row[index]!r} -> {native_row[index]!r}")
    gatk_info, native_info = info_fields(gatk_row), info_fields(native_row)
    for key in sorted(set(gatk_info) | set(native_info)):
        if gatk_info.get(key) != native_info.get(key):
            differences.append(
                f"INFO {key} {gatk_info.get(key)!r} -> {native_info.get(key)!r}")
    if len(gatk_row) > 9 and len(native_row) > 9:
        if gatk_row[8] != native_row[8]:
            differences.append(f"FORMAT keys {gatk_row[8]!r} -> {native_row[8]!r}")
        else:
            gatk_format, native_format = format_fields(gatk_row), format_fields(native_row)
            for key in gatk_format:
                if gatk_format.get(key) != native_format.get(key):
                    differences.append(
                        f"FORMAT {key} {gatk_format.get(key)!r} -> "
                        f"{native_format.get(key)!r}")
    return differences


def reference_block_pl(row: list[str]) -> str | None:
    """The ``PL`` text of a ``<NON_REF>``-only reference block row, else None."""
    if len(row) > 4 and row[4] == "<NON_REF>":
        return format_fields(row).get("PL")
    return None


def compare(gatk_rows: list[list[str]], native_rows: list[list[str]]) -> dict:
    violations: list[str] = []
    first: list[str] = []
    if len(gatk_rows) != len(native_rows):
        violations.append(f"ROW_COUNT gatk={len(gatk_rows)} native={len(native_rows)}")
    for index in range(max(len(gatk_rows), len(native_rows))):
        gatk_row = gatk_rows[index] if index < len(gatk_rows) else None
        native_row = native_rows[index] if index < len(native_rows) else None
        if gatk_row is None or native_row is None:
            violations.append(f"row {index}: only one side emitted")
            continue
        if gatk_row == native_row:
            continue
        detail = row_field_differences(gatk_row, native_row)
        if not first:
            first = detail
        violations.append(
            f"row {index} POS={gatk_row[1]}: " + "; ".join(detail or ["<column diff>"]))
    return {"data_rows_byte_identical": not violations,
            "violations": violations,
            "first_divergent_field_diff": first}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for HaplotypeCaller gVCF reference-block PL/GQ "
                    "parity under --alleles against pinned GATK 4.6.2.0.")
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-gvcf-refblock-oracle-") as directory:
        work = Path(directory)
        reference_text = build_reference()
        rows = [("chr1", position, reference_text[position - 1],
                 next(base for base in "ACGT" if base != reference_text[position - 1]))
                for position in SNP_POSITIONS]
        feature = write_feature(work, "forced-snps", rows, root)
        references: dict[bool, Path] = {}
        # The reference FASTA is shared by both read sets, so it is written and
        # indexed once (rewriting it invalidates the sequence dictionary).
        reference = work / "ref.fa"
        reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
        (work / "ref.fa.fai").write_text(
            f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
            encoding="utf-8")
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=300)
        for with_deletion_reads in (True, False):
            _, sam = build_inputs(work, reference_text, with_deletion_reads)
            bam = work / f"input.{'del' if with_deletion_reads else 'nodel'}.bam"
            subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam),
                            "-O", str(bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=300)
            references[with_deletion_reads] = bam

        for case in selected:
            bam = references[case["deletion_reads"]]
            allele_args = ["--alleles", str(feature)] if case["alleles"] else []
            common = ["-R", str(reference), "-I", str(bam), "-L", HOMOPOLYMER_INTERVAL,
                      *allele_args, "--min-pruning", "1",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false", *case["extra"]]
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
            result: dict = {
                "case": case["case"], "why": case["why"],
                "forced_alt_alleles": rows, "deletion_reads": case["deletion_reads"],
                "extra_args": list(case["extra"]),
                "gatk_exit": gatk_run.returncode, "native_exit": native_run.returncode,
                "gatk_stderr_tail": gatk_run.stderr[-400:] if gatk_run.returncode else "",
                "native_stderr_tail": (native_run.stderr[-400:]
                                       if native_run.returncode else ""),
            }
            if gatk_run.returncode != 0 or native_run.returncode != 0:
                result["comparison"] = {
                    "data_rows_byte_identical": False,
                    "violations": [f"run failed: gatk_exit={gatk_run.returncode} "
                                   f"native_exit={native_run.returncode}"],
                    "first_divergent_field_diff": []}
                result["gatk_rows"] = []
                result["native_rows"] = []
                results.append(result)
                continue
            gatk_rows = records(gatk_vcf)
            native_rows = records(native_vcf)
            result["comparison"] = compare(gatk_rows, native_rows)
            result["gatk_rows"] = gatk_rows
            result["native_rows"] = native_rows
            result["gatk_rows_total"] = len(gatk_rows)
            result["native_rows_total"] = len(native_rows)
            result["gatk_reference_block_pl"] = [reference_block_pl(row)
                                                 for row in gatk_rows
                                                 if reference_block_pl(row) is not None]
            result["native_reference_block_pl"] = [reference_block_pl(row)
                                                   for row in native_rows
                                                   if reference_block_pl(row) is not None]
            # Fixture-validity guard: the case must still exercise both a
            # concrete candidate row and `<NON_REF>` reference blocks.
            result["gatk_concrete_candidates"] = [
                row[1] for row in gatk_rows
                if len(row) > 4 and "<NON_REF>" in row[4] and row[4] != "<NON_REF>"]
            results.append(result)

    violations: list[str] = []
    for result in results:
        comparison = result["comparison"]
        if not comparison["data_rows_byte_identical"]:
            violations.append(f"{result['case']}: "
                              + "; ".join(comparison["violations"][:6]))
        if result["gatk_exit"] == 0 and (
                not result.get("gatk_concrete_candidates")
                or not result.get("gatk_reference_block_pl")):
            violations.append(
                f"{result['case']}: fixture no longer produces both a concrete "
                f"candidate and <NON_REF> reference blocks; GATK rows were "
                f"{result.get('gatk_rows_total')}")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: ReferenceConfidenceModel.doIndelRefConfCalc takes "
          "the least-confident of the SNP Ref-vs-Any model and the indel model "
          "built from nInformativeReads (ReferenceConfidenceModel.java:301-317, "
          "331-342), so the block PL/GQ is decided by nInformativeReads")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    for result in results:
        comparison = result["comparison"]
        print(f"[{result['case']}] deletion_reads={result['deletion_reads']} "
              f"extra={result['extra_args']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']} "
              f"data_rows_byte_identical={comparison['data_rows_byte_identical']} "
              f"gatk_rows={result.get('gatk_rows_total')} "
              f"native_rows={result.get('native_rows_total')}")
        if result.get("gatk_concrete_candidates"):
            print(f"    GATK concrete candidate POS={result['gatk_concrete_candidates']}")
        if result.get("gatk_reference_block_pl"):
            print(f"    GATK reference-block PL list: "
                  f"{result['gatk_reference_block_pl']}")
            print(f"    NATIVE reference-block PL list: "
                  f"{result['native_reference_block_pl']}")
        for row in result.get("gatk_rows", []):
            print(f"    GATK   {chr(9).join(row)}")
        for row in result.get("native_rows", []):
            print(f"    NATIVE {chr(9).join(row)}")
        for violation in comparison["violations"][:10]:
            print(f"    VIOLATION: {violation}")
        for tail in ("gatk_stderr_tail", "native_stderr_tail"):
            if result.get(tail):
                print(f"    {tail}: {result[tail]!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload_results: list[dict] = []
    for result in results:
        payload_results.append({key: value for key, value in result.items()
                                if key not in ("gatk_rows", "native_rows")})
    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller gVCF reference-block PL/GQ parity under --alleles",
        "binary": str(native),
        "acceptance_criterion": (
            "every data row (header/provenance ignored) must be byte-identical "
            "to pinned GATK 4.6.2.0 for every gated case, reference-block GQ and "
            "PL included; no canonicalization is applied"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": [case["case"] for case in CASES],
        "violations": violations,
        "results": payload_results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
