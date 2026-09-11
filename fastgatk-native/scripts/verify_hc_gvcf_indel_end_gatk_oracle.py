#!/usr/bin/env python3
"""Strict oracle for the gVCF native-only ``END`` INFO on concrete-indel records.

Observation under test
----------------------
``verify_hc_alleles_deep_limits.py`` reports ``gvcf-homopolymer-insertion`` as
``diverged`` with equal row counts (37 vs 37).  The field-level difference is a
single INFO key on the one concrete candidate row::

    GATK   chr1 697 . TC T,<NON_REF> 98.60 . BaseQRankSum=0.000;DP=36;...
    NATIVE chr1 697 . TC T,<NON_REF> 98.60 . END=698;BaseQRankSum=0.000;DP=36;...

Every other field of that row (QUAL, all annotations, ``GT:AD:DP:GQ:PL:SB`` and
every value) and every other row (including all ``<NON_REF>`` block rows and
their own ``END=``) is byte-identical.  The same INFO key appears **without**
``--alleles`` (the assembled ``chr1:697 TC>T`` deletion), so this is not an
``--alleles`` divergence; it is tracked as ``O1`` in
``.diag/track-b-alleles-findings.md``.

GATK's rule (from the pinned 4.6.2.0 sources, not from guesswork)
----------------------------------------------------------------
A HaplotypeCaller gVCF run is ordinary VariantContexts handed to
``GVCFWriter``/``GVCFBlockCombiner``.  A *variant* record is emitted verbatim
with whatever INFO its caller produced:

    // gatk-source/.../utils/variant/writers/GVCFBlockCombiner.java:203-209
    } else {
        // g is variant, so flush the bands and emit vc
        emitCurrentBlock();
        nextAvailableStart = vc.getEnd();
        contigOfNextAvailableStart = vc.getContig();
        toOutput.add(vc);
    }

The ``END`` attribute therefore exists only where GATK itself puts it, which is
the *hom-ref block* path:

    // gatk-source/.../utils/variant/writers/GVCFBlock.java:53
    vcb.attribute(VCFConstants.END_KEY, getEnd());

(plus ``ReblockingGVCFBlockCombiner``/``CombineGVCFs``/``ReblockGVCF``, none of
which the HaplotypeCaller path uses).  A grep of the whole HC genotyping path
(``HaplotypeCallerGenotypingEngine``/``ReferenceConfidenceModel``) finds no
``END`` attribute assignment at all, and
``ReferenceConfidenceModel.calculateRefConfidence`` only ever *reads*
``vc.getAttributeAsInt(END_KEY, vc.getStart())`` when feeding a hom-ref site to
the block combiner (``GVCFBlockCombiner.java:131``).  Consequently a concrete
indel call carries **no** ``END`` in GATK; ``VariantContext.getEnd()`` (which
does span the deletion) is used internally by the combiner but is not
serialized.

Native site under test
----------------------
``fastgatk-native/src/hc_call.cpp`` writes the gVCF candidate rows itself and
unconditionally adds ``END`` to any *concrete-indel* candidate::

    if (!bp_resolution && has_indel) {
        out << "END=" << reference_end;
        info_started = true;
    }

(one copy in the gVCF candidate path, one in the frequency-driven candidate
path).  ``has_indel`` is true whenever any emitted candidate's REF length
differs from its ALT length, so a forced homopolymer insertion, a forced
deletion, or an assembled deletion all pick it up.

What this oracle gates
----------------------
Every data row of every gated case must be byte-identical to pinned GATK
4.6.2.0 -- no canonicalization of any kind is applied.  The cases cover the
reported fixture, the same fixture without ``--alleles`` (the assembled
deletion), an assembled *insertion* candidate, and a ``-ERC BP_RESOLUTION``
control, so the gate cannot be satisfied by removing ``END`` from only the
``--alleles`` code path or only one of the two writer sites.

Exit status: 0 when every gated case's data rows are byte-identical to the
pinned jar, non-zero otherwise (so it can be registered in CTest).
``--expect-divergence`` turns the run into a pure diagnostic that always exits
0.
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
    HOMOPOLYMER_INTERVAL, HOMOPOLYMER_START, TANDEM_INTERVAL,
    records, write_feature, write_inputs,
)
import oracle_guard

# ``cases`` are gated: any difference in any data row is a violation.
CASES: tuple[dict, ...] = (
    {
        "case": "gvcf-homopolymer-insertion",
        "why": "the reported deep-limits case: one forced homopolymer insertion "
               "A>AA at chr1:600 emitted beside <NON_REF>; the concrete row is "
               "the left-aligned chr1:697 TC>T deletion",
        "interval": HOMOPOLYMER_INTERVAL,
        "rows": "forced-insertion",
        "extra": ("-ERC", "GVCF"),
    },
    {
        "case": "gvcf-assembled-deletion",
        "why": "same fixture and interval with no --alleles at all: the six "
               "one-bp-deletion reads still produce the concrete chr1:697 "
               "record, so this is the --alleles-free shape of the same INFO "
               "difference (tracked as O1)",
        "interval": HOMOPOLYMER_INTERVAL,
        "rows": None,
        "extra": ("-ERC", "GVCF"),
    },
    {
        "case": "gvcf-tandem-insertion",
        "why": "a forced insertion inside the tandem-repeat window, so the "
               "concrete candidate is an insertion rather than a deletion",
        "interval": TANDEM_INTERVAL,
        "rows": "tandem-insertion",
        "extra": ("-ERC", "GVCF"),
    },
)

# Run and printed, but never counted as violations: their rows differ for a
# *separately tracked* defect that this oracle is explicitly not allowed to
# touch.
PROBES: tuple[dict, ...] = (
    {
        "case": "bp-resolution-control",
        "why": "the same forced homopolymer insertion at -ERC BP_RESOLUTION: the "
               "concrete candidate carries no END on either side there, but this "
               "fixture additionally exposes an *unrelated*, pre-existing "
               "divergence on the reference row for the deleted base "
               "(chr1:698: GATK FORMAT 0/0:30,6:36:0:0,0,1135 vs native "
               "0/0:30,0:30:90:0,90,1343).  Kept as a probe so this gate stays a "
               "strict oracle for the END rule; that AD/DP/GQ/PL difference is "
               "reported here, neither fixed nor hidden.",
        "interval": HOMOPOLYMER_INTERVAL,
        "rows": "forced-insertion",
        "extra": ("-ERC", "BP_RESOLUTION"),
    },
)

ALL_CASES: tuple[dict, ...] = CASES + PROBES


def feature_rows(kind: str | None, reference_text: str) -> list[tuple] | None:
    if kind is None:
        return None
    if kind == "forced-insertion":
        return [("chr1", HOMOPOLYMER_START, "A", "AA")]
    if kind == "tandem-insertion":
        # A one-unit CAG insertion anchored on the first tandem-repeat base.
        return [("chr1", 900, reference_text[899], reference_text[899:902] + "CAG")]
    raise SystemExit(f"unknown fixture kind: {kind}")


def info_fields(row: list[str]) -> dict[str, str]:
    if len(row) < 8:
        return {}
    return dict(field.split("=", 1) if "=" in field else (field, "")
                for field in row[7].split(";"))


def row_field_differences(gatk_row: list[str], native_row: list[str]) -> list[str]:
    """Field-level diff of one VCF data row (0=CHROM, 3=REF, 4=ALT, 5=QUAL...)."""
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
            gatk_values = gatk_row[9].split(":")
            native_values = native_row[9].split(":")
            for index, key in enumerate(gatk_row[8].split(":")):
                gatk_value = gatk_values[index] if index < len(gatk_values) else None
                native_value = native_values[index] if index < len(native_values) else None
                if gatk_value != native_value:
                    differences.append(f"FORMAT {key} {gatk_value!r} -> {native_value!r}")
    return differences


def compare(gatk_rows: list[list[str]], native_rows: list[list[str]]) -> dict:
    violations: list[str] = []
    first: list[str] = []
    if len(gatk_rows) != len(native_rows):
        violations.append(f"ROW_COUNT gatk={len(gatk_rows)} native={len(native_rows)}")
    for index in range(max(len(gatk_rows), len(native_rows))):
        gatk_row = gatk_rows[index] if index < len(gatk_rows) else None
        native_row = native_rows[index] if index < len(native_rows) else None
        if gatk_row is None or native_row is None:
            violations.append(f"row {index}: only one side emitted (POS="
                              f"{(gatk_row or native_row)[1]})")
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
        description="Strict oracle for the native-only gVCF END INFO on "
                    "concrete-indel records against pinned GATK 4.6.2.0.")
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
        oracle_guard.oracle_not_verified('verify_hc_gvcf_indel_end_gatk_oracle.py', Path(java), gatk)
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-gvcf-end-oracle-") as directory:
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
            rows = feature_rows(case["rows"], reference_text)
            allele_args: list[str] = []
            if rows is not None:
                feature = write_feature(work, case["case"], rows, root)
                allele_args = ["--alleles", str(feature)]
            common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
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
                "case": case["case"], "why": case["why"], "gated": case["case"] in gated,
                "interval": case["interval"], "extra_args": list(case["extra"]),
                "forced_alleles": rows is not None,
                "gatk_exit": gatk_run.returncode,
                "native_exit": native_run.returncode,
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
            # Fixture-validity guard: a byte-identical comparison of two empty
            # files, or of a run whose concrete candidate vanished, would not be
            # exercising this rule.
            concrete = [row for row in gatk_rows
                        if len(row) > 4 and "<NON_REF>" in row[4] and row[4] != "<NON_REF>"]
            result["gatk_concrete_candidate_rows"] = [row[1] for row in concrete]
            result["gatk_concrete_candidate_has_end"] = [
                "END" in info_fields(row) for row in concrete]
            results.append(result)

    violations: list[str] = []
    probe_notes: list[str] = []
    for result in results:
        comparison = result["comparison"]
        case_violations: list[str] = []
        if not comparison["data_rows_byte_identical"]:
            case_violations.append("; ".join(comparison["violations"][:6]))
        if result["gatk_exit"] == 0 and not result.get("gatk_concrete_candidate_rows"):
            case_violations.append(
                f"fixture no longer produces a concrete (non-<NON_REF>) gVCF "
                f"candidate row; GATK rows were {result.get('gatk_rows_total')}")
        if not case_violations:
            continue
        note = f"{result['case']}: " + "; ".join(case_violations)
        if result["gated"]:
            violations.append(note)
        else:
            probe_notes.append(note)

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: a gVCF *variant* record is emitted verbatim by "
          "GVCFBlockCombiner.submit (GVCFBlockCombiner.java:203-209) and only "
          "GVCFBlock.toVariantContext (GVCFBlock.java:53) sets the END attribute, "
          "so concrete indel calls carry no END")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    for result in results:
        comparison = result["comparison"]
        print(f"[{result['case']}] interval={result['interval']} "
              f"extra={result['extra_args']} forced_alleles={result['forced_alleles']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']} "
              f"data_rows_byte_identical={comparison['data_rows_byte_identical']} "
              f"gatk_rows={result.get('gatk_rows_total')} "
              f"native_rows={result.get('native_rows_total')}")
        if result.get("gatk_concrete_candidate_rows"):
            print(f"    GATK concrete candidate POS="
                  f"{result['gatk_concrete_candidate_rows']} "
                  f"carrying END={result['gatk_concrete_candidate_has_end']}")
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
    if probe_notes:
        print(f"# {len(probe_notes)} out-of-scope probe note(s) (not counted as "
              f"violations):")
        for note in probe_notes:
            print(f"#   - {note}")

    payload_results: list[dict] = []
    for result in results:
        payload_results.append({key: value for key, value in result.items()
                                if key not in ("gatk_rows", "native_rows")})
    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller gVCF concrete-indel END INFO parity",
        "binary": str(native),
        "acceptance_criterion": (
            "every data row (header/provenance ignored) must be byte-identical "
            "to pinned GATK 4.6.2.0 for every gated case; no canonicalization "
            "is applied"),
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
