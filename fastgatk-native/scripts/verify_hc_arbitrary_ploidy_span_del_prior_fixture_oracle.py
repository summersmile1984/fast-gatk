#!/usr/bin/env python3
"""Fixture oracle for the arbitrary-ploidy ordinary-VCF `*` prior divergence.

Site under test
---------------
``fastgatk-native/src/hc_call.cpp:3358-3360`` (the arbitrary-ploidy
ordinary-VCF writer of ``fastgatk-hc-call``)::

    if (has_spanning_deletion)
        prior_pseudocounts.back() = result.genotype_indel_heterozygosity *
            ref_pseudocount;                        // unconditional indel prior

GATK rule (``gatk-source/.../afcalc/AlleleFrequencyCalculator.java:175-176``)::

    final double[] priorPseudocounts = alleles.stream()
        .mapToDouble(a -> a.isReference() ? refPseudocount
                        : (a.length() == refLength ? snpPseudocount : indelPseudocount))
        .toArray();

``refLength = vc.getReference().length()`` and htsjdk's ``Allele.SPAN_DEL`` is a
*called, non-symbolic, 1 bp* allele, so on a 1 bp REF record the spanning
deletion must take the **SNP** pseudocount.  ``hc_call.cpp:3358-3360`` passes the
indel pseudocount (1/8000 = 1.25e-4 instead of 1e-3, an 8x ratio = 0.903 log10),
which changes the Dirichlet posterior for every genotype containing ``*`` and
therefore ``P(no variant) = AFresult.log10ProbOnlyRefAlleleExists()`` -- exactly
the record's QUAL (``GenotypingEngine.java:158-161``).

The ready-to-apply one-line fix is recorded verbatim in
``.diag/round-length-prior-sweep.md`` section C.2 item 1::

    if (has_spanning_deletion)
        prior_pseudocounts.back() =
            (candidate_reference(calls.front()->candidate).size() == 1
                 ? result.genotype_snp_heterozygosity
                 : result.genotype_indel_heterozygosity) * ref_pseudocount;

Why the previous round could not see it, and what the fixture does differently
-----------------------------------------------------------------------------
The ``*`` pseudocount is only visible when ``*`` is in the AF matrix but is **not
in the called genotype**: once a genotype containing ``*`` is the called one, its
effective allele count is data-dominated and an 8x prior change moves the
posterior by ~1e-4 (invisible at the printed 2 decimals).  The previous round's
probe (``--sample-ploidy 3`` on a 6-deletion-read pileup) had
``GT = 0/1/2 = REF/A/*``, i.e. ``*`` **was** called, so the line executed and
changed nothing observable (QUAL 358.64 on both sides).

This oracle therefore uses a pileup in which the deletion is still assembled
(so ``has_spanning_deletion`` is true and `*` enters the AF matrix) but is
supported by only 2-3 of the reads, so the called genotype stays ``REF/A/A`` and
the ``*`` allele's effective count is ~0 -- the prior-dominated regime:

``--sample-ploidy 3``, 30 x 300M reference reads + 12 x 300M ``chr1:698 C>A``
reads + 3 x ``138M1D161M`` reads carrying a one-base deletion of ``chr1:700``::

    GATK   chr1 698 . C A 368.90 ... QD=8.20 ... GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244
    NATIVE chr1 698 . C A 369.37 ... QD=8.21 ... GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244

QUAL (and the QUAL-derived QD) is the **only** field difference; GT/AD/PL/GQ/INFO
are byte-identical.  Native's QUAL is the higher one, i.e. the ``*`` allele is
under-weighted, matching the direction of the already-fixed gVCF site
(``hc_call.cpp:4820-4836``: GATK 282.04 / native 291.07 pre-fix).

Controls (must be byte-identical, and they are)
-----------------------------------------------
* ``ord-ploidy3-no-deletion`` -- the same pileup without the deletion reads, so
  no spanning deletion exists at all (QUAL 421.02 on both sides).  If the
  divergence were caused by anything other than the ``*`` allele this control
  would fail too.
* ``ord-ploidy2-cigar-del3`` -- identical pileup at the diploid ploidy, where
  ``has_spanning_deletion`` is hard-coded false in this writer
  (``hc_call.cpp:3224``) and the diploid path (already fixed in ``4d66251``) is
  taken: byte-identical.
* ``ord-ploidy1-snp30-del6`` -- haploid control (QUAL 423.04 on both sides).

Exit status
-----------
Strict by default: every data row of every case must be identical after the
documented allele-order canonicalization, so the run exits 1 while the defect is
unfixed (this is a red gate, deliberately).  ``--expect-divergence`` asserts the
opposite direction -- the QUAL divergence must be present at the fixture's
spanning-deletion site and every other field must still agree -- so it exits 0
while the defect is reproducible and 1 if the fixture stops reaching the line.
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
    ORDINARY_VCF, compare_rows, layout, make_bam, prepare_reference, qual_summary,
    records,
)

SPAN_DEL_POSITION = "698"

CASES: tuple[dict, ...] = (
    {
        "case": "ord-ploidy3-cigar-del3",
        "ploidy": 3,
        "why": "ordinary VCF, --sample-ploidy 3: 30 reference + 12 chr1:698 C>A "
               "+ 3 reads with a real 138M1D161M one-base deletion of chr1:700. "
               "The deletion is assembled (so the spanning-deletion allele "
               "enters the AF matrix) but not called, so the prior-dominated "
               "regime is reached at chr1:698",
        "expect_divergence": True,
        "layout": layout(30, 12, 3, "cigar"),
    },
    {
        "case": "ord-ploidy4-cigar-del2",
        "ploidy": 4,
        "why": "ordinary VCF, --sample-ploidy 4, 2 deletion reads: independent "
               "ploidy, smaller divergence (QUAL 391.84 / 391.94) but the same "
               "direction and the same site",
        "expect_divergence": True,
        "layout": layout(30, 12, 2, "cigar"),
    },
    {
        "case": "ord-ploidy3-no-deletion",
        "ploidy": 3,
        "why": "CONTROL: identical pileup with no deletion reads at all, so no "
               "spanning deletion exists and the prior cannot matter",
        "expect_divergence": False,
        "layout": layout(30, 12, 0),
    },
    {
        "case": "ord-ploidy2-cigar-del3",
        "ploidy": 2,
        "why": "CONTROL: identical pileup at diploid ploidy, where this writer "
               "sets has_spanning_deletion=false (hc_call.cpp:3224) and the "
               "already-fixed diploid AF path is used",
        "expect_divergence": False,
        "layout": layout(30, 12, 3, "cigar"),
    },
    {
        "case": "ord-ploidy1-snp30-del6",
        "ploidy": 1,
        "why": "CONTROL: haploid pileup (20 reference + 30 C>A + 6 deletion "
               "reads) where the concrete ALT is called outright",
        "expect_divergence": False,
        # The haploid control needs a deeper reference group so the concrete
        # ALT is called; build it explicitly instead of via layout().
        "layout": {"ref": {"count": 20, "start": 440, "step": 6, "length": 300},
                   "snp": {"count": 30, "start": 440, "step": 6, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "cigar", "drop_base": 700}},
    },
)

GATED_CASES = tuple(case["case"] for case in CASES if case["expect_divergence"])
CONTROL_CASES = tuple(case["case"] for case in CASES if not case["expect_divergence"])


def run_case(case: dict, reference: Path, reference_text: str, work: Path,
             java: str, gatk: str, native: str, threads: int,
             strict: bool) -> dict:
    bam = make_bam(work, case["case"], case["layout"], reference_text, java, gatk)
    gatk_vcf = work / f"gatk.{case['case']}.vcf"
    native_vcf = work / f"native.{case['case']}.vcf"
    common = ["-R", str(reference), "-I", str(bam), "-L", "chr1:500-780",
              "--min-pruning", "1",
              "--sample-ploidy", str(case["ploidy"]),
              "--create-output-variant-index", "false",
              "--add-output-vcf-command-line", "false", *ORDINARY_VCF]
    gatk_run = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "HaplotypeCaller", *common, "-O", str(gatk_vcf)],
        capture_output=True, text=True, timeout=900)
    native_run = subprocess.run(
        [native, *common, "--threads", str(threads), "-O", str(native_vcf)],
        capture_output=True, text=True, timeout=900)
    result = {"case": case["case"], "ploidy": case["ploidy"], "why": case["why"],
              "expect_divergence": case["expect_divergence"],
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
        allow_divergence=case["expect_divergence"] and not strict)
    comparison["gatk_rows"] = len(gatk_rows)
    comparison["native_rows"] = len(native_rows)
    comparison["gatk_qual"] = qual_summary(gatk_rows)
    comparison["native_qual"] = qual_summary(native_rows)
    comparison["gatk_rows_at_site"] = [row for row in gatk_rows
                                       if row[1] == SPAN_DEL_POSITION]
    comparison["native_rows_at_site"] = [row for row in native_rows
                                         if row[1] == SPAN_DEL_POSITION]
    result["comparison"] = comparison
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict fixture oracle for the arbitrary-ploidy ordinary-VCF "
                    "spanning-deletion `*` AF prior divergence vs pinned GATK 4.6.2.0.")
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-arbitrary-ploidy-prior-") as directory:
        work = Path(directory)
        reference = prepare_reference(work, java, str(gatk))
        from hc_symbolic_prior_fixture_lib import build_reference
        reference_text = build_reference()
        for case in selected:
            results.append(run_case(case, reference, reference_text, work, java,
                                    str(gatk), str(native), arguments.threads, strict))

    violations: list[str] = []
    for result in results:
        comparison = result["comparison"]
        for item in comparison.get("violations", []):
            violations.append(f"{result['case']}: {item}")
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
    print("# site: hc_call.cpp:3358-3360 -- arbitrary-ploidy ordinary-VCF writer "
          "assigns the `*` allele the indel pseudocount instead of GATK's "
          "length-derived (SNP on a 1 bp REF) pseudocount")
    print(f"# mode: {'expect-divergence (assert the divergence is present)' if not strict else 'strict parity assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated (divergence expected): {list(GATED_CASES)}")
    print(f"# controls (byte parity required): {list(CONTROL_CASES)}")
    for result in results:
        comparison = result["comparison"]
        print(f"[{result['case']}] --sample-ploidy {result['ploidy']} "
              f"expect_divergence={result['expect_divergence']} "
              f"gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_rows={comparison.get('gatk_rows')} "
              f"native_rows={comparison.get('native_rows')} "
              f"data_rows_identical={comparison.get('data_rows_identical')} "
              f"allele_order_normalizations="
              f"{len(comparison.get('allele_order_normalizations', []))}")
        for item in comparison.get("allele_order_normalizations", []):
            print(f"    normalized: {item}")
        for item in comparison.get("qual_divergences", []):
            print(f"    QUAL DIVERGENCE (GATK -> NATIVE): {item}")
        for item in comparison.get("unexpected_qual_mismatches", []):
            print(f"    UNEXPECTED QUAL MISMATCH: {item}")
        for item in comparison.get("separately_tracked_notes", []):
            print(f"    separately tracked (reported, not gated): {item}")
        for item in comparison.get("violations", []):
            print(f"    VIOLATION: {item}")
        for row in comparison.get("gatk_rows_at_site", []):
            print(f"    GATK   {chr(9).join(row)}")
        for row in comparison.get("native_rows_at_site", []):
            print(f"    NATIVE {chr(9).join(row)}")
        if comparison.get("first_divergent_row_field_diff"):
            print("    first divergent row field diff (GATK -> NATIVE): "
                  + comparison["first_divergent_row_field_diff"])
        if result["native_exit"] != 0:
            print(f"    native stderr tail: {result['native_stderr_tail']!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload = {
        "status": "fail" if violations else "pass",
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller arbitrary-ploidy ordinary-VCF spanning-deletion "
                  "`*` AF prior parity (hc_call.cpp:3358-3360)",
        "binary": str(native),
        "strict_mode": strict,
        "expected_fix_location": "fastgatk-native/src/hc_call.cpp:3358-3360",
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
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
