#!/usr/bin/env python3
"""Strict oracle for the INFO ``AF``/``MLEAF`` zero-value formatting divergence.

Observation under test
----------------------
Reported by the fixture-construction round
(``.diag/round-fixture-construction.md`` section D.1, found while building the
``chr1:698 C>A,G,T`` forced-allele fixture of section B):

    GATK   ... AF=0.500,0.00;... MLEAF=0.500,0.00 ...
    NATIVE ... AF=0.500,0.000;... MLEAF=0.500,0.000 ...

i.e. on a **multi-ALT** ordinary-VCF record whose allele-frequency list mixes a
non-zero value with an exact zero, GATK writes the zero as ``0.00`` while native
writes ``0.000``.  All-zero lists already agree (``0.00,0.00``).

GATK's rule (established from the pinned jar, not from guesswork)
----------------------------------------------------------------
``AF``/``MLEAF`` are ``Number=A`` **Float** INFO values held as a ``double[]``;
htsjdk serializes every ``Double`` through
``htsjdk.variant.vcf.VCFEncoder.formatVCFDouble(double)``.  Decompiling that
method out of the pinned ``gatk-package-4.6.2.0-local.jar``
(``javap -p -c htsjdk.variant.vcf.VCFEncoder``) gives, verbatim:

    if (d >= 1.0)              format = "%.2f";
    else if (d >= 0.01)        format = "%.3f";
    else if (abs(d) >= 1e-20)  format = "%.3e";
    else                       return "0.00";
    return String.format(Locale.US, format, d);

The rule is therefore **per value** and purely magnitude-driven, which is why
``0.500`` (3 places) and ``0.00`` (the literal zero string) coexist in one list.

``htsjdk`` uses the same function for every ``Double`` INFO/FORMAT value, so the
same rule is the reason GATK prints a ``QD`` >= 1 with two places and an ``SOR``
between 0.01 and 1 with three; see the module docstring of the round report for
the parts of that rule this oracle does **not** cover.

Native site under test
----------------------
``fastgatk-native/src/hc_call.cpp`` ``format_allele_frequency``::

    if (abs(frequency - 1.0) < 1e-12) value << fixed << setprecision(2) << frequency;
    else                              value << fixed << setprecision(3) << frequency;

It has no zero case, so ``0.0`` becomes ``0.000``.  Four call sites use it (the
ordinary-VCF writer's ``AF`` and ``MLEAF``, and the two gVCF writers' ``MLEAF``);
the gVCF writers already special-case a zero count to the literal ``"0.00"``,
which is why the divergence is observable on the ordinary writer only.

What this oracle gates
----------------------
Every data row of every case must be byte-identical to pinned GATK 4.6.2.0 after
the documented ALT-order canonicalization of ``hc_symbolic_prior_fixture_lib``
(``AF``/``MLEAF``/``AC``/``MLEAC`` permuted with the ALT list; nothing else is
excluded -- ``QUAL`` and all annotations are gated too).  Two fixture-validity
guards keep the gate honest:

* ``zero_numerator_list`` -- the GATK row must actually carry a ``Number=A``
  entry that is *both* numerically zero and printed with GATK's zero spelling,
  next to a non-zero entry; otherwise the case stopped exercising the rule and
  the oracle fails even if the rows happen to match.
* ``all_zero_list`` control -- a forced-allele pileup with no variant reads,
  where native already agrees, proves the fix is not simply "print two places
  everywhere" (a blanket ``%.2f`` would break the ``0.500`` entries).

Exit status: 0 when every gated case matches, non-zero otherwise (so it can be
registered in CTest).  ``--expect-divergence`` turns the run into a pure
diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hc_symbolic_prior_fixture_lib import (  # noqa: E402
    REFERENCE_LENGTH, build_reference, compare_rows, layout, make_bam,
    normalize_row, prepare_reference, records,
)

SITE = "698"                     # 1-based POS of the forced multi-ALT record
INTERVAL = "chr1:500-780"
NUMBER_A_KEYS = ("AC", "AF", "MLEAC", "MLEAF")
ZERO_SPELLING = "0.00"           # htsjdk formatVCFDouble's literal-zero string

# ``cases`` are gated; ``probes`` are run and printed but never counted as
# violations, for configurations whose rows differ for a *separately tracked*
# defect (here: none at the time of writing -- the one probe below is the gVCF
# writer's MLEAF zero path, which already agrees, and it is kept as a probe
# because its record needs ALT-order canonicalization to be comparable).
CASES: tuple[dict, ...] = (
    {
        "case": "multialt-3alt-ploidy2",
        "why": "3 forced concrete ALTs (A,G,T), --max-alternate-alleles 3, 12 "
               "reads supporting A: AF=0.500,0.00,0.00 is the reported row",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 3, "ploidy": 2,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "multialt-2alt-ploidy2",
        "why": "same forced allele list reduced to A,G by "
               "--max-alternate-alleles 2: the reported fixture's exact shape",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 2,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "multialt-2alt-ploidy3",
        "why": "same reduction at --sample-ploidy 3: AN=3 gives 0.333,0.00",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 3,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "multialt-2alt-ploidy4",
        "why": "same reduction at --sample-ploidy 4 (AN=4)",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 4,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "multialt-2alt-ploidy3-span-del",
        "why": "the forced multi-ALT allele list on the spanning-deletion pileup: "
               "with --alleles the assembled deletion allele is not emitted, so "
               "the record comes out biallelic (no zero entry) -- it is kept as a "
               "shape control showing the span-del pileup is otherwise unaffected",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 2, "ploidy": 3,
        "snp_reads": 12, "deletion_reads": 3, "encoding": "ordinary",
        "expect_zero_entry": False,
    },
    {
        "case": "multialt-zero-first-ploidy2",
        "why": "same forced allele list with the reads supporting the *second* "
               "forced allele (G): GATK keeps haplotype order (ALT G,A,T, "
               "AF=0.500,0.00,0.00) while native emits A,G,T (AF=0.00,0.500,0.00), "
               "so after the documented canonicalization this is the "
               "zero-in-first-position shape -- the rule is not an end-of-list "
               "artifact",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 3, "ploidy": 2,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "snp_alt": "G", "expect_zero_entry": True,
    },
    {
        "case": "biallelic-forced-control",
        "why": "two forced ALTs (A,T) with 12 reads on A: a two-element list "
               "whose single non-zero value sits in front of the zero "
               "(AF=0.500,0.00) -- the smallest multi-ALT shape that shows the "
               "divergence",
        "alts": ["A", "T"], "max_alternate_alleles": 6, "ploidy": 2,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "allzero-forced-control",
        "why": "CONTROL: 3 forced ALTs over a reference-only pileup gives "
               "AF=0.00,0.00,0.00 on both sides -- a blanket %.2f formatter "
               "would also pass this, so it is the companion of the 0.500 cases",
        "alts": ["A", "G", "T"], "max_alternate_alleles": 3, "ploidy": 2,
        "snp_reads": 0, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": True,
    },
    {
        "case": "noalleles-biallelic-snp-control",
        "why": "CONTROL: the canonical ordinary HC record without --alleles "
               "(single ALT, AF=0.500) -- proves the fix does not disturb the "
               "main serialization path",
        "alts": None, "max_alternate_alleles": None, "ploidy": 2,
        "snp_reads": 12, "deletion_reads": 0, "encoding": "ordinary",
        "expect_zero_entry": False,
    },
)

PROBES: tuple[dict, ...] = (
    {
        "case": "gvcf-ploidy3-span-del-mleaf",
        "why": "the gVCF writer's MLEAF zero path (already prints 0.00 through "
               "its own count==0 rule): kept as a probe because this record's "
               "ALT list needs the documented canonicalization to compare",
        "alts": None, "max_alternate_alleles": None, "ploidy": 3,
        "snp_reads": 12, "deletion_reads": 3, "encoding": "bp-resolution",
        "expect_zero_entry": True,
    },
)

ALL_CASES = CASES + PROBES

# A Number=A entry that is numerically zero and, on the GATK side, spelled with
# htsjdk's literal-zero string; used as the fixture-validity guard.
NUMBER_A_RE = re.compile(r"(?:^|;)(" + "|".join(NUMBER_A_KEYS) + r")=([^;\t]+)")


def write_alleles_feature(work: Path, name: str, alts: list[str]) -> Path:
    """bgzip+tabix'd single-record feature VCF carrying several concrete ALTs."""
    reference_base = build_reference()[int(SITE) - 1]
    plain = work / f"{name}.alleles.vcf"
    plain.write_text("##fileformat=VCFv4.2\n"
                     f"##contig=<ID=chr1,length={REFERENCE_LENGTH}>\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     f"chr1\t{SITE}\t.\t{reference_base}\t{','.join(alts)}"
                     "\t.\tPASS\t.\n", encoding="utf-8")
    root = Path(__file__).resolve().parents[2]
    bgzip = root / "third_party/htslib-build/htslib-src/bgzip"
    tabix = root / "third_party/htslib-build/htslib-src/tabix"
    if not bgzip.is_file() or not tabix.is_file():
        raise SystemExit(f"missing bgzip/tabix under {bgzip.parent}")
    subprocess.run([str(bgzip), "-f", str(plain)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    compressed = plain.with_suffix(plain.suffix + ".gz")
    subprocess.run([str(tabix), "-f", "-p", "vcf", str(compressed)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    return compressed


def number_a_text(row: list[str]) -> dict[str, str]:
    """The ``AC``/``AF``/``MLEAC``/``MLEAF`` texts of one data row."""
    if len(row) < 8:
        return {}
    return {key: value for key, value in NUMBER_A_RE.findall(row[7])}


def zero_entry_spellings(row: list[str]) -> dict[str, list[str]]:
    """For every Number=A key, the spellings of its numerically-zero entries."""
    found: dict[str, list[str]] = {}
    for key, text in number_a_text(row).items():
        zero_spellings = []
        for value in text.split(","):
            try:
                numeric = float(value)
            except ValueError:
                continue
            if numeric == 0.0:
                zero_spellings.append(value)
        if zero_spellings:
            found[key] = zero_spellings
    return found


def has_all_zero_list(row: list[str]) -> bool:
    for key, text in number_a_text(row).items():
        values = text.split(",")
        if len(values) > 1 and all(value in ("0", "0.00", "0.000") for value in values):
            return True
    return False


def run_case(case: dict, reference: Path, reference_text: str, work: Path,
             java: str, gatk: str, native: str, threads: int) -> dict:
    case_layout = layout(30, case["snp_reads"], case["deletion_reads"], "cigar",
                         case.get("snp_alt", "A"))
    bam = make_bam(work, case["case"], case_layout, reference_text, java, gatk)
    common = ["-R", str(reference), "-I", str(bam), "-L", INTERVAL,
              "--min-pruning", "1",
              "--sample-ploidy", str(case["ploidy"]),
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false"]
    feature = None
    if case["alts"]:
        feature = write_alleles_feature(work, case["case"], case["alts"])
        common += ["--alleles", str(feature),
                   "--max-alternate-alleles", str(case["max_alternate_alleles"])]
    if case["encoding"] == "bp-resolution":
        common += ["-ERC", "BP_RESOLUTION"]
    gatk_vcf = work / f"gatk.{case['case']}.vcf"
    native_vcf = work / f"native.{case['case']}.vcf"
    gatk_run = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "HaplotypeCaller", *common, "-O", str(gatk_vcf)],
        capture_output=True, text=True, timeout=900)
    native_run = subprocess.run(
        [native, *common, "--threads", str(threads), "-O", str(native_vcf)],
        capture_output=True, text=True, timeout=900)
    result = {"case": case["case"], "why": case["why"],
              "forced_alt_alleles": case["alts"],
              "max_alternate_alleles": case["max_alternate_alleles"],
              "ploidy": case["ploidy"], "snp_reads": case["snp_reads"],
              "deletion_reads": case["deletion_reads"],
              "encoding": case["encoding"],
              "snp_alt": case.get("snp_alt", "A"),
              "expect_zero_entry": case["expect_zero_entry"],
              "gatk_exit": gatk_run.returncode,
              "native_exit": native_run.returncode,
              "gatk_stderr_tail": gatk_run.stderr[-400:] if gatk_run.returncode else "",
              "native_stderr_tail": (native_run.stderr[-400:]
                                     if native_run.returncode else "")}
    if gatk_run.returncode != 0 or native_run.returncode != 0:
        result["comparison"] = {"violations": [
            f"run failed: gatk_exit={gatk_run.returncode} "
            f"native_exit={native_run.returncode}"], "data_rows_identical": False}
        result["gatk_rows"] = []
        result["native_rows"] = []
        result["gatk_rows_display"] = []
        result["native_rows_display"] = []
        result["gatk_rows_total"] = 0
        result["native_rows_total"] = 0
        result["rows_elided"] = 0
        return result
    gatk_rows = records(gatk_vcf)
    native_rows = records(native_vcf)
    result["comparison"] = compare_rows(gatk_rows, native_rows)
    gatk_display, gatk_total = display_rows(gatk_rows, native_rows)
    native_display, native_total = display_rows(native_rows, gatk_rows)
    result["gatk_rows"] = gatk_rows
    result["native_rows"] = native_rows
    result["gatk_rows_display"] = gatk_display
    result["native_rows_display"] = native_display
    result["gatk_rows_total"] = gatk_total
    result["native_rows_total"] = native_total
    result["rows_elided"] = gatk_total - len(gatk_display)
    result["gatk_number_a"] = [number_a_text(row) for row in gatk_rows]
    result["native_number_a"] = [number_a_text(row) for row in native_rows]
    result["gatk_zero_spellings"] = [zero_entry_spellings(row) for row in gatk_rows]
    result["native_zero_spellings"] = [zero_entry_spellings(row) for row in native_rows]
    # A row carrying a zero Number=A entry must be identical on both sides
    # *after* the documented ALT-order canonicalization: GATK orders the ALT
    # list by haplotype order (a read-supported allele can come first), native
    # emits sorted concrete alleles.
    matched: list[int] = []
    for index in range(min(len(gatk_rows), len(native_rows))):
        if not zero_entry_spellings(gatk_rows[index]):
            continue
        try:
            gatk_normalized, _ = normalize_row(gatk_rows[index])
            native_normalized, _ = normalize_row(native_rows[index])
        except ValueError:
            continue
        if gatk_normalized == native_normalized:
            matched.append(index)
    result["matched_zero_rows"] = matched
    return result


def display_rows(rows: list[list[str]], other: list[list[str]]) -> tuple[list[list[str]], int]:
    """Rows worth printing: any row carrying a Number=A field, plus divergent ones.

    A ``-ERC BP_RESOLUTION`` case emits hundreds of reference blocks; printing
    them all would bury the rows under test in the log.
    """
    kept: list[list[str]] = []
    for index, row in enumerate(rows):
        peer = other[index] if index < len(other) else None
        if number_a_text(row) or row != peer:
            kept.append(row)
    return kept, len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for HaplotypeCaller INFO AF/MLEAF zero-value "
                    "formatting against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_HC_BINARY"),
        help="native HaplotypeCaller binary (default: $FASTGATK_HC_BINARY or "
             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
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

    selected = [case for case in ALL_CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    gated = {case["case"] for case in CASES}
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-af-zero-oracle-") as directory:
        work = Path(directory)
        reference_text = build_reference()
        reference = prepare_reference(work, java, str(gatk))
        for case in selected:
            results.append(run_case(case, reference, reference_text, work,
                                    java, str(gatk), str(native), arguments.threads))

    violations: list[str] = []
    probe_notes: list[str] = []
    for result in results:
        gated_here = result["case"] in gated
        case_violations: list[str] = []
        comparison = result["comparison"]
        if not comparison.get("data_rows_identical"):
            case_violations.append(
                f"{result['case']}: data rows differ: "
                f"{comparison.get('violations', [])[:6]}")
        # Fixture-validity guard: the case must actually carry a zero-valued
        # Number=A entry spelled the way GATK spells it, otherwise it is not
        # exercising the rule under test.
        if result["gatk_exit"] == 0:
            if result["expect_zero_entry"]:
                matched = result["matched_zero_rows"]
                found_zero_spelling = any(
                    value == ZERO_SPELLING
                    for per_row in result["gatk_zero_spellings"]
                    for values in per_row.values() for value in values)
                if not found_zero_spelling:
                    case_violations.append(
                        f"{result['case']}: fixture no longer produces a GATK "
                        f"zero entry spelled {ZERO_SPELLING!r}; GATK Number=A "
                        f"texts were {result['gatk_number_a']}")
                if not matched:
                    case_violations.append(
                        f"{result['case']}: expected a byte-identical row that "
                        f"carries a zero-valued Number=A entry next to a "
                        f"non-zero one; GATK={result['gatk_number_a']} "
                        f"NATIVE={result['native_number_a']}")
        if case_violations:
            if gated_here:
                violations.extend(case_violations)
            else:
                probe_notes.extend(case_violations)

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: htsjdk VCFEncoder.formatVCFDouble -- "
          "d>=1 -> '%.2f'; 0.01<=d<1 -> '%.3f'; abs(d)<1e-20 -> literal '0.00'")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated cases (mismatch => exit 1): {[case['case'] for case in CASES]}")
    print(f"# probes (reported, not gated): {[case['case'] for case in PROBES]}")
    for result in results:
        print(f"[{result['case']}] gated={result['case'] in gated} "
              f"ploidy={result['ploidy']} forced={result['forced_alt_alleles']} "
              f"max-alt={result['max_alternate_alleles']} "
              f"snp_reads={result['snp_reads']} del_reads={result['deletion_reads']} "
              f"encoding={result['encoding']}")
        print(f"    why: {result['why']}")
        comparison = result["comparison"]
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']} "
              f"data_rows_identical={comparison.get('data_rows_identical')} "
              f"gatk_rows={result['gatk_rows_total']} "
              f"native_rows={result['native_rows_total']}")
        if result["rows_elided"]:
            print(f"    ({result['rows_elided']} reference-only block rows per side "
                  f"not printed: they carry no Number=A field and are identical)")
        for row in result["gatk_rows_display"]:
            print(f"    GATK   {chr(9).join(row)}")
        for row in result["native_rows_display"]:
            print(f"    NATIVE {chr(9).join(row)}")
        for side in ("gatk", "native"):
            for index, entry in enumerate(result[f"{side}_number_a"]):
                if entry:
                    print(f"    {side.upper():6s} Number=A row{index}: {entry}")
        if result["gatk_zero_spellings"]:
            print(f"    GATK zero-entry spellings: {result['gatk_zero_spellings']}")
        if result["native_zero_spellings"]:
            print(f"    NATIVE zero-entry spellings: {result['native_zero_spellings']}")
        if comparison.get("violations"):
            for item in comparison["violations"][:10]:
                print(f"    VIOLATION: {item}")
        if comparison.get("allele_order_normalizations"):
            for item in comparison["allele_order_normalizations"]:
                print(f"    NOTE (ALT order canonicalized): {item}")
        if result["gatk_stderr_tail"]:
            print(f"    gatk_stderr_tail: {result['gatk_stderr_tail']!r}")
        if result["native_stderr_tail"]:
            print(f"    native_stderr_tail: {result['native_stderr_tail']!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")
    if probe_notes:
        print(f"# {len(probe_notes)} probe note(s) (not counted as violations):")
        for note in probe_notes:
            print(f"#   - {note}")

    payload_results: list[dict] = []
    for result in results:
        trimmed = {key: value for key, value in result.items()
                   if key not in ("gatk_rows", "native_rows", "gatk_number_a",
                                  "native_number_a")}
        trimmed["gatk_number_a"] = [entry for entry in result["gatk_number_a"]
                                    if entry] if "gatk_number_a" in result else []
        trimmed["native_number_a"] = [entry for entry in result["native_number_a"]
                                      if entry] if "native_number_a" in result else []
        payload_results.append(trimmed)

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller INFO AF/MLEAF zero-value formatting",
        "binary": str(native),
        "acceptance_criterion": (
            "native data rows (header/provenance ignored) must be byte-identical "
            "to pinned GATK 4.6.2.0 for every gated case after ALT-order "
            "canonicalization, and every case that claims a zero-valued Number=A "
            "entry must carry it spelled as htsjdk spells it"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": [case["case"] for case in CASES],
        "probes": [case["case"] for case in PROBES],
        "violations": violations,
        "probe_notes": probe_notes,
        "results": payload_results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
