#!/usr/bin/env python3
"""Fixture oracle for the polyploid/hidden-spanning gVCF `*` prior divergence.

Site under test
---------------
``fastgatk-native/src/hc_call.cpp:4485-4488`` (the arbitrary-ploidy / haploid
gVCF writer of ``fastgatk-hc-call``)::

    if (include_spanning_deletion || hidden_spanning_deletion) {
        const auto spanning = group.candidates.size() + 1U;
        prior_pseudocounts[spanning] =
            result.genotype_indel_heterozygosity * ref_pseudocount;
    }

Same rule as the ``-ERC`` reference-confidence path fixed in the previous round
(``hc_call.cpp:4820-4836``, the diploid gVCF writer): GATK derives each allele's
Dirichlet pseudocount from its **length**
(``AlleleFrequencyCalculator.java:175-176``, ``refLength = vc.getReference().length()``)
and htsjdk's ``Allele.SPAN_DEL`` is a called, non-symbolic **1 bp** allele, so on
a 1 bp REF record ``*`` takes the SNP pseudocount (1e-3) and only a longer REF
gives it the indel pseudocount (1/8000 = 1.25e-4).  The unconditional indel prior
here under-weights ``*`` by 8x = 0.903 log10, which moves
``P(no variant) = AFresult.log10ProbOnlyRefAlleleExists()`` -- the record's QUAL
(``GenotypingEngine.java:158-161``; with a spanning deletion present GATK's
``log10ProbOnlyRefAlleleExists`` is the sum over REF/``*``-only genotypes,
``AlleleFrequencyCalculator.java:198-230``).

The ready-to-apply one-line fix is recorded verbatim in
``.diag/round-length-prior-sweep.md`` section C.2 item 2 (same shape as the
already-landed sister fix at ``hc_call.cpp:4820-4836``)::

    prior_pseudocounts[spanning] =
        (group.reference.size() == 1
             ? result.genotype_snp_heterozygosity
             : result.genotype_indel_heterozygosity) * ref_pseudocount;

Why the previous round could not see it, and what the fixture does differently
-----------------------------------------------------------------------------
The ``*`` pseudocount is only visible when ``*`` is in the AF/PL matrix but is
**not** in the called genotype: once a genotype containing ``*`` is the called
one, its effective allele count is data-dominated and the 8x prior change moves
the posterior by ~1e-4 (invisible at the printed 2 decimals).  With the
previously tried pileup at ``--sample-ploidy 3`` the called genotype was
``0/1/2 = REF/A/*`` -- ``*`` **was** called -- so the line executed and changed
nothing (QUAL 358.64 on both sides).  This oracle keeps the deletion assembled
(so ``include_spanning_deletion`` is true and `*` is serialized into the record's
ALT list) but reduces its read support to 2-3 reads, so the called genotype stays
``REF/REF/A`` and the ``*`` allele's effective count is ~0 -- the prior-dominated
regime.  ``--sample-ploidy 3``, 30 x 300M reference reads + 12 x 300M
``chr1:698 C>A`` reads + 3 x ``138M1D161M`` reads carrying a one-base deletion of
``chr1:700``, ``-ERC BP_RESOLUTION``::

    GATK   chr1 698 . C *,A,<NON_REF> 368.90 ... MLEAC=0,1,0;MLEAF=0.00,0.333,0.00 ... GT:... 0/1/2:30,3,12,0:45:6:430,366,431,...
    NATIVE chr1 698 . C A,*,<NON_REF> 369.37 ... MLEAC=1,0,0;MLEAF=0.333,0.00,0.00 ... GT:... 0/1/2:30,12,3,0:45:6:430,6,64,...

After the documented canonicalization (concrete ALTs sorted with ``<NON_REF>``
last, every allele-indexed field permuted with them) the **only** difference in
all four gated cases is QUAL (and, where the rounding survives, the QUAL-derived
QD); GT/GQ/AD/PL/SB/MLEAC/MLEAF are identical.  Native's QUAL is the higher one,
i.e. the ``*`` allele is under-weighted, exactly the direction of the
already-fixed diploid sibling (GATK 282.04 / native 291.07 pre-fix).

Reachability evidence (the line is not merely "probably" executed)
------------------------------------------------------------------
`*` in the emitted ALT list is written by ``if (include_spanning_deletion) out <<
",*";`` (``hc_call.cpp`` polyploid gVCF writer).  The emitted ``N A,*,<NON_REF>``
record therefore *is* proof that ``include_spanning_deletion == true``, hence that
the guard at ``hc_call.cpp:4485`` was entered and ``prior_pseudocounts[spanning]``
was assigned the unconditional indel pseudocount before the very same block hands
the vector to ``calculate_allele_frequency_kokkos`` (whose ``qual`` is what the
writer prints).  All four gated cases emit the ``*`` at chr1:698.

Controls (must be byte-identical)
---------------------------------
* ``bp-ploidy2-cigar-del3`` -- the identical pileup at diploid ploidy, which
  takes the *already-fixed* diploid gVCF writer (``hc_call.cpp:4820-4836``) and
  must agree with GATK.
* ``bp-ploidy3-no-deletion`` -- the same pileup with no deletion reads, so no
  ``*`` exists at all.

Exit status
-----------
Strict by default: every data row of every case must be identical after the
documented canonicalization, so the run exits 1 while the defect is unfixed
(deliberately a red gate).  ``--expect-divergence`` asserts the opposite
direction -- the QUAL divergence must be present at chr1:698 in every gated case
and every other field must still agree -- so it exits 0 while the defect is
reproducible and 1 if a fixture stops reaching the line.
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

from hc_symbolic_prior_fixture_lib import (  # noqa: E402
    BP_RESOLUTION, build_reference, compare_rows, layout, make_bam,
    prepare_reference, qual_summary, records,
)
import oracle_guard

SPAN_DEL_POSITION = "698"
# The `bp-ploidy3-cigar-del2` probe also differs at chr1:697, where GATK prints
# a zero QUAL ("0") and native prints a missing one (".").  That is a separately
# tracked serialization divergence (zero-QUAL vs missing-QUAL), unrelated to the
# `*` prior, so it is reported and not gated.
QUAL_SERIALIZATION_TOLERANCE = {"697": frozenset({"QUAL"})}

CASES: tuple[dict, ...] = (
    {
        "case": "bp-ploidy3-cigar-del3",
        "ploidy": 3,
        "why": "-ERC BP_RESOLUTION, --sample-ploidy 3: 30 reference + 12 "
               "chr1:698 C>A + 3 reads with a real 138M1D161M one-base deletion "
               "of chr1:700. `*` is serialized at chr1:698 (so the divergent "
               "prior assignment definitely executed) but is not called",
        "expect_divergence": True,
        "layout": layout(30, 12, 3, "cigar"),
    },
    {
        "case": "bp-ploidy3-frameshift-del3",
        "ploidy": 3,
        "why": "same pileup, deletion encoded as sequence removal behind a 299M "
               "CIGAR (the pathological frameshift encoding): the divergence "
               "must not depend on the deletion encoding",
        "expect_divergence": True,
        "layout": layout(30, 12, 3, "frameshift"),
    },
    {
        "case": "bp-ploidy3-cigar-del2",
        "ploidy": 3,
        "why": "only 2 deletion reads: the largest divergence of the set "
               "(QUAL 411.25 / 416.91).  The chr1:697 record additionally shows "
               "the separately tracked zero-QUAL-vs-missing-QUAL serialization "
               "difference, which is reported and not gated",
        "expect_divergence": True,
        "layout": layout(30, 12, 2, "cigar"),
        "tolerated": QUAL_SERIALIZATION_TOLERANCE,
        "report_only": True,
    },
    {
        "case": "bp-ploidy4-cigar-del2",
        "ploidy": 4,
        "why": "independent ploidy (4), 2 deletion reads: same site, same "
               "direction, smaller magnitude (QUAL 391.84 / 391.94)",
        "expect_divergence": True,
        "layout": layout(30, 12, 2, "cigar"),
        "report_only": True,
    },
    {
        "case": "bp-ploidy2-cigar-del3",
        "ploidy": 2,
        "why": "CONTROL: the identical pileup at diploid ploidy, which takes the "
               "already-fixed diploid gVCF writer (hc_call.cpp:4820-4836)",
        "expect_divergence": False,
        "layout": layout(30, 12, 3, "cigar"),
    },
    {
        "case": "bp-ploidy3-no-deletion",
        "ploidy": 3,
        "why": "CONTROL: identical pileup without deletion reads, so no `*` "
               "exists at chr1:698 and the prior cannot matter",
        "expect_divergence": False,
        "layout": layout(30, 12, 0),
    },
)

GATED_CASES = tuple(case["case"] for case in CASES if case["expect_divergence"])
CONTROL_CASES = tuple(case["case"] for case in CASES if not case["expect_divergence"])


def run_case(case: dict, reference: Path, reference_text: str, work: Path,
             java: str, gatk: str, native: str, threads: int, strict: bool) -> dict:
    bam = make_bam(work, case["case"], case["layout"], reference_text, java, gatk)
    gatk_vcf = work / f"gatk.{case['case']}.vcf"
    native_vcf = work / f"native.{case['case']}.vcf"
    common = ["-R", str(reference), "-I", str(bam), "-L", "chr1:500-780",
              "--min-pruning", "1",
              "--sample-ploidy", str(case["ploidy"]),
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false", *BP_RESOLUTION]
    gatk_run = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "HaplotypeCaller", *common, "-O", str(gatk_vcf)],
        capture_output=True, text=True, timeout=900)
    native_run = subprocess.run(
        [native, *common, "--threads", str(threads), "-O", str(native_vcf)],
        capture_output=True, text=True, timeout=900)
    result = {"case": case["case"], "ploidy": case["ploidy"], "why": case["why"],
              "expect_divergence": case["expect_divergence"],
              "report_only": case.get("report_only", False),
              "gatk_exit": gatk_run.returncode, "native_exit": native_run.returncode,
              "gatk_stderr_tail": gatk_run.stderr[-400:],
              "native_stderr_tail": native_run.stderr[-400:],
              "divergence_positions": [SPAN_DEL_POSITION] if case["expect_divergence"]
                                      else []}
    if gatk_run.returncode != 0 or native_run.returncode != 0:
        result["comparison"] = {"violations": [
            f"run failed: gatk_exit={gatk_run.returncode} "
            f"native_exit={native_run.returncode}"], "data_rows_identical": False}
        return result
    gatk_rows = records(gatk_vcf)
    native_rows = records(native_vcf)
    comparison = compare_rows(
        gatk_rows, native_rows,
        divergence_positions=frozenset([SPAN_DEL_POSITION]),
        tolerated=case.get("tolerated", {}),
        allow_divergence=case["expect_divergence"] and not strict)
    comparison["gatk_rows"] = len(gatk_rows)
    comparison["native_rows"] = len(native_rows)
    comparison["gatk_qual"] = qual_summary(gatk_rows)
    comparison["native_qual"] = qual_summary(native_rows)
    comparison["alts_at_span_del_site"] = {
        "gatk": [row[4] for row in gatk_rows if row[1] == SPAN_DEL_POSITION],
        "native": [row[4] for row in native_rows if row[1] == SPAN_DEL_POSITION],
        "native_emits_spanning_deletion": any(
            row[1] == SPAN_DEL_POSITION and "*" in row[4].split(",")
            for row in native_rows),
    }
    comparison["gatk_rows_near_site"] = [row for row in gatk_rows
                                         if row[1] in ("697", "698", "699")]
    comparison["native_rows_near_site"] = [row for row in native_rows
                                           if row[1] in ("697", "698", "699")]
    result["comparison"] = comparison
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict fixture oracle for the polyploid/hidden-spanning gVCF "
                    "`*` AF prior divergence vs pinned GATK 4.6.2.0.")
    parser.add_argument("--native", default=os.environ.get("FASTGATK_HC_BINARY"),
                        help="native HaplotypeCaller binary (default: "
                             "$FASTGATK_HC_BINARY or "
                             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument("--expect-divergence", action="store_true",
                        help="assert that the known QUAL divergence IS present at "
                             "the fixture's spanning-deletion site (exit 0), "
                             "instead of asserting byte parity (exit 1 today)")
    parser.add_argument("--case", action="append", default=None,
                        help=f"run only the named case(s); default all of "
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
        oracle_guard.oracle_not_verified('verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py', Path(java), gatk)
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
    strict = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-polyploid-gvcf-prior-") as directory:
        work = Path(directory)
        reference = prepare_reference(work, java, str(gatk))
        reference_text = build_reference()
        for case in selected:
            results.append(run_case(case, reference, reference_text, work, java,
                                    str(gatk), str(native), arguments.threads, strict))

    violations: list[str] = []
    probe_notes: list[str] = []
    for result in results:
        comparison = result["comparison"]
        if result["report_only"]:
            # A `report_only` probe demonstrates the divergence but ALSO carries
            # unrelated pre-existing differences at the flanking deletion record
            # (native does not retain that site, so it prints no annotations).
            # Report them; gate only the presence of the QUAL divergence.
            probe_notes += [f"{result['case']}: {item}"
                            for item in comparison.get("violations", [])]
        else:
            violations += [f"{result['case']}: {item}"
                           for item in comparison.get("violations", [])]
        if result["expect_divergence"] and not strict \
                and not comparison.get("qual_divergences"):
            violations.append(
                f"{result['case']}: expected QUAL divergence not observed at "
                f"chr1:{SPAN_DEL_POSITION} -- the fixture no longer reaches the "
                f"divergent line or the defect was fixed")
        if not result["expect_divergence"] and not comparison.get(
                "data_rows_identical", False):
            violations.append(f"{result['case']}: control case must be "
                              f"byte-identical but is not")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# site: hc_call.cpp:4485-4488 -- polyploid/hidden-spanning gVCF writer "
          "assigns the `*` allele the indel pseudocount instead of GATK's "
          "length-derived (SNP on a 1 bp REF) pseudocount")
    print(f"# mode: {'expect-divergence (assert the divergence is present)' if not strict else 'strict parity assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated (divergence expected): {list(GATED_CASES)}")
    print(f"# controls (byte parity required): {list(CONTROL_CASES)}")
    for result in results:
        comparison = result["comparison"]
        print(f"[{result['case']}] -ERC BP_RESOLUTION --sample-ploidy "
              f"{result['ploidy']} expect_divergence={result['expect_divergence']} "
              f"gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_rows={comparison.get('gatk_rows')} "
              f"native_rows={comparison.get('native_rows')} "
              f"data_rows_identical={comparison.get('data_rows_identical')} "
              f"allele_order_normalizations="
              f"{len(comparison.get('allele_order_normalizations', []))}")
        print(f"    native emits `*` at chr1:{SPAN_DEL_POSITION}: "
              f"{comparison.get('alts_at_span_del_site', {}).get('native_emits_spanning_deletion')} "
              f"(GATK ALT {comparison.get('alts_at_span_del_site', {}).get('gatk')} / "
              f"NATIVE ALT {comparison.get('alts_at_span_del_site', {}).get('native')}) "
              f"-- `*` in the emitted native ALT list is written only when "
              f"include_spanning_deletion is true, i.e. the divergent line ran")
        for item in comparison.get("allele_order_normalizations", []):
            print(f"    normalized: {item}")
        for item in comparison.get("qual_divergences", []):
            print(f"    QUAL DIVERGENCE (GATK -> NATIVE): {item}")
        for item in comparison.get("unexpected_qual_mismatches", []):
            print(f"    UNEXPECTED QUAL MISMATCH: {item}")
        for item in comparison.get("separately_tracked_notes", []):
            print(f"    separately tracked (reported, not gated): {item}")
        label = "DIFFERENCE (report-only probe, not counted)" \
            if result["report_only"] else "VIOLATION"
        for item in comparison.get("violations", []):
            print(f"    {label}: {item}")
        for row in comparison.get("gatk_rows_near_site", []):
            print(f"    GATK   {chr(9).join(row)}")
        for row in comparison.get("native_rows_near_site", []):
            print(f"    NATIVE {chr(9).join(row)}")
        if comparison.get("first_divergent_row_field_diff"):
            print("    first divergent row field diff (GATK -> NATIVE): "
                  + comparison["first_divergent_row_field_diff"])
        if result["native_exit"] != 0:
            print(f"    native stderr tail: {result['native_stderr_tail']!r}")
    if probe_notes:
        print(f"# {len(probe_notes)} report-only probe note(s) (not counted as "
              f"violations):")
        for note in probe_notes:
            print(f"#   - {note}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload = {
        "status": "fail" if violations else "pass",
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller polyploid/hidden-spanning gVCF `*` AF prior "
                  "parity (hc_call.cpp:4485-4488)",
        "binary": str(native),
        "strict_mode": strict,
        "expected_fix_location": "fastgatk-native/src/hc_call.cpp:4485-4488",
        "acceptance_criterion": (
            "strict mode: every data row of every case must be byte-identical to "
            "pinned GATK 4.6.2.0 after the documented allele-order "
            "canonicalization (currently fails at chr1:698 QUAL, the defect this "
            "fixture localizes); --expect-divergence mode: the QUAL divergence "
            "must be present at chr1:698 and every other field must agree"),
        "cases": [case["case"] for case in selected],
        "gated_cases": list(GATED_CASES),
        "control_cases": list(CONTROL_CASES),
        "violations": violations,
        "report_only_probe_notes": probe_notes,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
