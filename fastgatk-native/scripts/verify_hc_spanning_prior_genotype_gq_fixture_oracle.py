#!/usr/bin/env python3
"""Fixture oracle for the `*` genotype-prior divergence in the gVCF GT/GQ path.

Site under test
---------------
``fastgatk-native/src/hc_call.cpp:682-687`` (``genotype_priors_for_group``)::

    if (include_spanning_deletion) {
        const auto spanning_heterozygosity = result.genotype_indel_heterozygosity;
        ...
        allele_priors.emplace_back(log10_spanning_het, 2.0 * log10_spanning_het);
    }

GATK rule (``gatk-source/.../utils/genotyper/GenotypePriorCalculator.java:141-153``)::

    } else if (allele.isCalled() && !allele.isSymbolic()) {
        return allele.length() == referenceLength ? AlleleType.SNP : AlleleType.INDEL;
    } else if (allele.equals(Allele.SV_SIMPLE_INS) || allele.equals(Allele.SV_SIMPLE_DEL)) {
        throw new IllegalArgumentException("cannot handle symbolic indels: " + allele);
    } else { return AlleleType.OTHER; }

with ``hetValues[SNP] = log10(snpHet) - log10(3)``, ``hetValues[INDEL] = log10(indelHet)``,
``hetValues[OTHER] = max(...)`` (``GenotypePriorCalculator.java:62-72, 118-121``).  The
spanning-deletion allele ``*`` is *called and non-symbolic* with htsjdk length 1, so on
a 1 bp REF record GATK classifies it as **SNP** and gives it
``log10(1e-3) - log10(3) = -3.4771``; native gives it ``log10(1/8000) = -3.9031``.
The difference is ``0.426`` in log10 = ``4.26`` in PL/GQ phred units, and it moves the
posterior margin that decides ``GQ`` (and, in a near tie, ``GT``).

The ready-to-apply one-line fix is recorded verbatim in
``.diag/round-length-prior-sweep.md`` section C.2 item 3 (the ``hc_call.cpp`` half;
``calling_pipeline.cpp:405-424`` ``prior_allele_type`` is the sibling half).

Why the previous round could not see it, and what the fixture does differently
------------------------------------------------------------------------------
Two conditions must hold together, and the previous round's fixtures had neither:

1. ``*`` must be in the genotype-prior enumeration, i.e.
   ``genotype_priors_for_group(..., include_spanning_deletion=true)`` -- reached by a
   record that serializes ``*``;
2. the priors must actually be **used**, i.e.
   ``result.genotype_priors_used = options.use_genotype_priors &&
   options.use_posterior_genotype_assignment`` (``calling_pipeline.cpp:14046``) --
   with HaplotypeCaller's default ``USE_PLS_TO_ASSIGN`` the callers take
   ``derive_genotype_gt_gq_kokkos`` and the priors are never read, which is exactly
   why "GT/GQ byte-identical in all six gated fixtures" was observed before.

This oracle therefore adds ``--genotype-assignment-method USE_POSTERIOR_PROBABILITIES``
to the same span-del pileup, at ``--sample-ploidy 3`` and ``-ERC BP_RESOLUTION``
(30 reference + 12 ``chr1:698 C>A`` + N deleting reads).  Measured:

    GATK   chr1 698 . C *,A,<NON_REF> 368.90 ... GT:AD:DP:GP:GQ:PG:PL:SB  0/0/2:30,3,12,0:45:...:29:...:430,366,431,...
    NATIVE chr1 698 . C A,*,<NON_REF> 369.37 ... GT:AD:DP:GQ:PL:SB        0/0/1:30,12,3,0:45:33:430,6,64,...

after canonicalization the called genotype is the same (``REF/REF/A``) and every field
except ``GQ`` agrees: GATK 29 vs native 33 -- ``delta = 4``, the predicted 4.26.

Controls (byte-identical, and they are)
---------------------------------------
* ``priors-no-spanning-deletion`` -- the same pileup without deleting reads: no ``*``
  allele exists, so the priors are concrete-only and ``GQ`` agrees exactly (84/84).
* ``priors-disabled-span-del`` -- the same span-del pileup **without**
  ``--genotype-assignment-method``: priors are not used, ``GQ`` agrees exactly (6/6).
  The pair of controls brackets the divergence to "spanning deletion present AND
  posterior assignment enabled", which is precisely the guarded line.

Separately tracked (reported, not gated)
----------------------------------------
With ``USE_POSTERIOR_PROBABILITIES`` GATK adds the ``GP`` and ``PG`` FORMAT fields
(Number=G posteriors); native emits neither.  That is a distinct FORMAT-emission gap
(gVCF-only; the ordinary-VCF writer has it too), reported per record and excluded from
this gate so the ``GQ`` comparison stays attributable.

Exit status
-----------
Strict by default: every gated case's rows must be identical except ``GQ`` and the
reported ``GP``/``PG`` gap, so the run exits 1 while the defect is unfixed (a red
gate).  ``--expect-divergence`` asserts that the ``GQ`` divergence IS present on every
gated case and that the two controls remain byte-identical, so it exits 0 while the
defect is reproducible and 1 if a fixture stops reaching the line.
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
    prepare_reference, records,
)
import oracle_guard

SPAN_DEL_POSITION = "698"
POSTERIOR_ASSIGNMENT = ("--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES")
# GATK emits GP/PG (Number=G posteriors) only under USE_POSTERIOR_PROBABILITIES;
# native emits neither.  Reported, never gated (see the module docstring).
POSTERIOR_FIELD_TOLERANCE = {"*": frozenset({"GP", "PG", "keys"})}
# The same record also carries the *separately tracked* gVCF `*` AF-pseudocount
# defect in its QUAL (and the QUAL-derived QD): GATK 368.90 / native 369.37 at
# chr1:698.  That one is gated by
# fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py;
# here it is reported and excluded so the GQ comparison stays attributable.
QUAL_TOLERANCE = {SPAN_DEL_POSITION: frozenset({"QUAL", "QD"})}
# `priors-p3-cigar-del2` additionally shows the pre-existing zero-QUAL vs
# missing-QUAL serialization difference at chr1:697 (native does not retain that
# site), reported and not gated.
FLANK_TOLERANCE = {"697": frozenset({"QUAL", "MLEUAC", "MLEAC", "MLEAF", "DP",
                                     "BaseQRankSum", "MQRankSum", "ReadPosRankSum",
                                     "RAW_MQandDP", "SB"})}

CASES: tuple[dict, ...] = (
    {
        "case": "priors-p3-cigar-del3",
        "ploidy": 3, "deletion_reads": 3, "posterior": True,
        "expect_divergence": True,
        "why": "30 reference + 12 chr1:698 C>A + 3 x 138M1D161M deletion reads, "
               "`-ERC BP_RESOLUTION --sample-ploidy 3 "
               "--genotype-assignment-method USE_POSTERIOR_PROBABILITIES`: `*` is "
               "serialized and the priors are used, so GQ carries the 4.26-unit "
               "prior difference (GATK 29 / native 33)",
    },
    {
        "case": "priors-p3-cigar-del2",
        "ploidy": 3, "deletion_reads": 2, "posterior": True,
        "expect_divergence": True,
        "why": "2 deletion reads: same GQ divergence (61/65) with a different "
               "allele balance",
    },
    {
        "case": "priors-p3-cigar-del6",
        "ploidy": 3, "deletion_reads": 6, "posterior": True,
        "expect_divergence": True,
        "why": "6 deletion reads: `*` is now inside the CALLED genotype on both "
               "sides (0/1/2), so the same root cause moves GQ the other way "
               "(GATK 67 / native 63) -- the divergence is not one-signed",
    },
    {
        "case": "priors-p4-cigar-del2",
        "ploidy": 4, "deletion_reads": 2, "posterior": True,
        "expect_divergence": True,
        "why": "independent ploidy 4: same site, same 4-unit GQ divergence "
               "(GATK 27 / native 31)",
    },
    {
        "case": "priors-no-spanning-deletion",
        "ploidy": 3, "deletion_reads": 0, "posterior": True,
        "expect_divergence": False,
        "why": "CONTROL: no deletion reads, so no `*` allele exists and the "
               "genotype priors are concrete-only -- GQ must agree exactly (84/84)",
    },
    {
        "case": "priors-disabled-span-del",
        "ploidy": 3, "deletion_reads": 3, "posterior": False,
        "expect_divergence": False,
        "why": "CONTROL: the span-del pileup without "
               "--genotype-assignment-method USE_POSTERIOR_PROBABILITIES: the "
               "priors are never read, so GQ must agree exactly (6/6)",
    },
)

GATED_CASES = tuple(case["case"] for case in CASES if case["expect_divergence"])
CONTROL_CASES = tuple(case["case"] for case in CASES if not case["expect_divergence"])


def genotype_quality(row: list[str]) -> str:
    if len(row) < 10:
        return "-"
    keys = row[8].split(":")
    values = row[9].split(":")
    return values[keys.index("GQ")] if "GQ" in keys else "-"


def genotype_text(row: list[str]) -> str:
    if len(row) < 10:
        return "-"
    keys = row[8].split(":")
    values = row[9].split(":")
    return values[keys.index("GT")] if "GT" in keys else "-"


def run_case(case: dict, reference: Path, reference_text: str, work: Path,
             java: str, gatk: str, native: str, threads: int, strict: bool) -> dict:
    bam = make_bam(work, case["case"],
                   layout(30, 12, case["deletion_reads"], "cigar"),
                   reference_text, java, gatk)
    gatk_vcf = work / f"gatk.{case['case']}.vcf"
    native_vcf = work / f"native.{case['case']}.vcf"
    common = ["-R", str(reference), "-I", str(bam), "-L", "chr1:500-780",
              "--min-pruning", "1",
              "--sample-ploidy", str(case["ploidy"]),
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false", *BP_RESOLUTION]
    if case["posterior"]:
        common += list(POSTERIOR_ASSIGNMENT)
    gatk_run = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "HaplotypeCaller", *common, "-O", str(gatk_vcf)],
        capture_output=True, text=True, timeout=900)
    native_run = subprocess.run(
        [native, *common, "--threads", str(threads), "-O", str(native_vcf)],
        capture_output=True, text=True, timeout=900)
    result = {"case": case["case"], "why": case["why"], "ploidy": case["ploidy"],
              "deletion_reads": case["deletion_reads"],
              "posterior_assignment": case["posterior"],
              "expect_divergence": case["expect_divergence"],
              "gatk_exit": gatk_run.returncode, "native_exit": native_run.returncode,
              "gatk_stderr_tail": gatk_run.stderr[-400:],
              "native_stderr_tail": native_run.stderr[-400:]}
    if gatk_run.returncode != 0 or native_run.returncode != 0:
        result["comparison"] = {"violations": [
            f"run failed: gatk_exit={gatk_run.returncode} "
            f"native_exit={native_run.returncode}"], "data_rows_identical": False}
        return result
    gatk_rows = records(gatk_vcf)
    native_rows = records(native_vcf)
    tolerated: dict[str, frozenset] = {}
    if case["posterior"]:
        tolerated.update(POSTERIOR_FIELD_TOLERANCE)
    tolerated.update(QUAL_TOLERANCE)
    if case["deletion_reads"] and case["deletion_reads"] < 3:
        tolerated.update(FLANK_TOLERANCE)
    comparison = compare_rows(
        gatk_rows, native_rows,
        divergence_positions=frozenset([SPAN_DEL_POSITION]),
        tolerated=tolerated,
        allow_divergence=case["expect_divergence"] and not strict,
        divergence_fields=frozenset({"GQ"}))
    comparison["gatk_rows"] = len(gatk_rows)
    comparison["native_rows"] = len(native_rows)
    comparison["gatk_rows_near_site"] = [row for row in gatk_rows
                                         if row[1] in ("697", "698")]
    comparison["native_rows_near_site"] = [row for row in native_rows
                                           if row[1] in ("697", "698")]
    comparison["gatk_gq_at_site"] = [genotype_quality(row) for row in gatk_rows
                                     if row[1] == SPAN_DEL_POSITION]
    comparison["native_gq_at_site"] = [genotype_quality(row) for row in native_rows
                                       if row[1] == SPAN_DEL_POSITION]
    comparison["gatk_gt_at_site"] = [genotype_text(row) for row in gatk_rows
                                     if row[1] == SPAN_DEL_POSITION]
    comparison["native_gt_at_site"] = [genotype_text(row) for row in native_rows
                                       if row[1] == SPAN_DEL_POSITION]
    result["comparison"] = comparison
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fixture oracle for the symbolic spanning-deletion `*` "
                    "genotype-prior divergence in the gVCF GT/GQ path "
                    "(hc_call.cpp:682-694) vs pinned GATK 4.6.2.0.")
    parser.add_argument("--native", default=os.environ.get("FASTGATK_HC_BINARY"),
                        help="native HaplotypeCaller binary (default: "
                             "$FASTGATK_HC_BINARY or "
                             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument("--expect-divergence", action="store_true",
                        help="assert that the known GQ divergence IS present on the "
                             "gated cases and that the controls stay byte-identical "
                             "(exit 0), instead of asserting parity (exit 1 today)")
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
        oracle_guard.oracle_not_verified('verify_hc_spanning_prior_genotype_gq_fixture_oracle.py', Path(java), gatk)
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-genotype-prior-") as directory:
        work = Path(directory)
        reference = prepare_reference(work, java, str(gatk))
        reference_text = build_reference()
        for case in selected:
            results.append(run_case(case, reference, reference_text, work, java,
                                    str(gatk), str(native), arguments.threads, strict))

    violations: list[str] = []
    for result in results:
        comparison = result["comparison"]
        violations += [f"{result['case']}: {item}"
                       for item in comparison.get("violations", [])]
        if result["expect_divergence"] and not strict \
                and not comparison.get("qual_divergences"):
            violations.append(
                f"{result['case']}: expected GQ divergence not observed at "
                f"chr1:{SPAN_DEL_POSITION} -- the fixture no longer reaches the "
                f"divergent line or the defect was fixed")
        if not result["expect_divergence"] and not comparison.get(
                "data_rows_identical", False):
            violations.append(f"{result['case']}: control case must be "
                              f"byte-identical but is not")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# site: hc_call.cpp:682-694 -- genotype_priors_for_group assigns the "
          "spanning-deletion `*` allele log10(indelHet) instead of GATK's "
          "GenotypePriorCalculator SNP value log10(snpHet) - log10(3) "
          "(delta = 0.426 log10 = 4.26 phred)")
    print(f"# mode: {'expect-divergence (assert the divergence is present)' if not strict else 'strict parity assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated (GQ divergence expected): {list(GATED_CASES)}")
    print(f"# controls (byte parity required): {list(CONTROL_CASES)}")
    for result in results:
        comparison = result["comparison"]
        print(f"[{result['case']}] -ERC BP_RESOLUTION --sample-ploidy "
              f"{result['ploidy']} deletion_reads={result['deletion_reads']} "
              f"posterior_assignment={result['posterior_assignment']} "
              f"expect_divergence={result['expect_divergence']} "
              f"gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    why: {result['why']}")
        print(f"    data_rows_identical={comparison.get('data_rows_identical')} "
              f"gatk_rows={comparison.get('gatk_rows')} "
              f"native_rows={comparison.get('native_rows')} "
              f"GQ GATK={comparison.get('gatk_gq_at_site')} "
              f"NATIVE={comparison.get('native_gq_at_site')} "
              f"GT GATK={comparison.get('gatk_gt_at_site')} "
              f"NATIVE={comparison.get('native_gt_at_site')}")
        for item in comparison.get("qual_divergences", []):
            print(f"    GQ DIVERGENCE (GATK -> NATIVE): {item}")
        for item in comparison.get("unexpected_qual_mismatches", []):
            print(f"    UNEXPECTED GQ MISMATCH: {item}")
        for item in comparison.get("separately_tracked_notes", []):
            print(f"    separately tracked (reported, not gated): {item}")
        for item in comparison.get("violations", []):
            print(f"    VIOLATION: {item}")
        for row in comparison.get("gatk_rows_near_site", []):
            print(f"    GATK   {chr(9).join(row)}")
        for row in comparison.get("native_rows_near_site", []):
            print(f"    NATIVE {chr(9).join(row)}")
        if result["native_exit"] != 0:
            print(f"    native stderr tail: {result['native_stderr_tail']!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload = {
        "status": "fail" if violations else "pass",
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller gVCF spanning-deletion `*` genotype-prior/GQ "
                  "parity (hc_call.cpp:682-694)",
        "binary": str(native),
        "strict_mode": strict,
        "expected_fix_location": "fastgatk-native/src/hc_call.cpp:682-687 "
                                 "(sibling: calling_pipeline.cpp:405-424 "
                                 "prior_allele_type)",
        "acceptance_criterion": (
            "with `*` serialized and --genotype-assignment-method "
            "USE_POSTERIOR_PROBABILITIES the GQ of the chr1:698 record must equal "
            "GATK's (it differs by 4 today); with `*` absent, or with the priors "
            "unused, every field of every row must be identical"),
        "cases": [case["case"] for case in selected],
        "gated_cases": list(GATED_CASES),
        "control_cases": list(CONTROL_CASES),
        "violations": violations,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
