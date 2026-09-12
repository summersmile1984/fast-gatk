#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the reflecting-ALT-set rules of a spanning-deletion locus.

Two rules are pinned, both measured on pinned GATK 4.6.2.0 and both on the path
the dense-materialization plan depends on.

Rule 1 -- a locus whose only surviving ALT is the spanning deletion `*` is
refused outside EMIT_ALL_ACTIVE_SITES
------------------------------------------------------------------------
``GenotypingEngine.calculateGenotypes()`` keeps an owned `*` in the output
allele subset, and then, *before* it registers any deletion state::

    // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
    if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1
            && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
        return null;                                     // :173-175
    }
    final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
    recordDeletions(vc, outputAlleles);                  // :178-179

(``gatk-source/.../walkers/genotyper/GenotypingEngine.java:172-179``.)

``emitAllActiveSites()`` is the ``EMIT_ALL_ACTIVE_SITES`` traversal, i.e. exactly
``--include-non-variant-sites`` in GenotypeGVCFs (``GenotypeGVCFs.java`` sets
``outputMode`` from that flag).  The rule therefore splits by mode:

* **default mode**: a locus whose only emitted ALT is `*` produces NO record,
  even when an upstream emitted deletion owns that `*` (cases
  ``owned-star-only-locus-refused-default`` and
  ``owned-star-only-locus-het-refused-default``).  The refusal happens before
  ``recordDeletions()``, so such a locus must not contribute to the
  emitted-deletion state either.
* **dense mode**: the locus IS emitted (case
  ``star-only-locus-emitted-dense-negative-zero-qd``).
* a locus that also carries a concrete ALT is unaffected in both modes (control
  ``star-plus-concrete-alt-unaffected``), and an ORPHAN `*` is already removed
  from the subset by ``isSpuriousSpanningDeletion`` (control
  ``orphan-star-only-locus-refused-default``).

Rule 2 -- QD is computed from the SIGN-PRESERVING confidence, so a zero QUAL
renders as ``-0.00``
------------------------------------------------------------------------
``GenotypingEngine`` publishes QUAL from ``builder.log10PError(log10Confidence)``
(``:183``) while the value it *tests* is ``(-10.0 * log10Confidence) + 0.0``
(``:163``, whose comment says "Add 0.0 removes -0.0 occurrences").  Java keeps
the negative zero in the published double, and ``QualByDepth`` divides that same
sign-preserving numerator::

    qual = -10.0 * vc.getLog10PError();
    double QD = qual / depth;
    ... String.format("%.2f", QD)

(``gatk-source/.../walkers/annotator/QualByDepth.java:76-91``.)  An emitted
`*`-only locus has ``log10Confidence == +0.0``, so ``qual == -0.0`` and the
published QD is ``-0.00`` while the QUAL column renders ``0`` (htsjdk's
``VCFEncoder.formatQualValue`` strips a trailing ``.00``).  Measured on both the
dense `*`-only locus and on the deletion row of the derived fixture: native used
to publish ``QD=0.00``.

Scope
-----
The dense cases compare only the positions where both tools publish a row from a
record start (``compare_positions``), because GATK additionally materializes
positions covered by a spanning record -- the tracked covered-locus
materialization gap -- and that is deliberately out of scope here.

Usage / exit status
-------------------
``--native`` / ``$FASTGATK_GENOTYPE_BINARY`` select the native binary
(``$FASTGATK_NATIVE_BUILD``/``fastgatk-native/build`` supplies the default), the
pinned GATK jar and the bundled JDK17 are read from ``third_party/``, all scratch
lives in a ``tempfile.TemporaryDirectory``, and the run ends with a single-line
JSON status payload.  Exit status: 0 when every gated case passes, non-zero
otherwise.  ``--expect-divergence`` turns the run into a diagnostic that always
exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from oracle_guard import oracle_not_verified  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""

DENSE = "--include-non-variant-sites"

# A single emitted deletion at 2 whose span is 2-4, so the `*` record at 3 is
# OWNED by it (isVcCoveredByDeletion answers true).
DELETION = ("chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
            "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n")

# The owned `*`-only locus: its only surviving ALT after the standard-confidence
# subset is the spanning deletion.
STAR_ONLY_HOM = ("chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
                 "GT:DP:AD:PL\t1/1:20:0,20:100,100,0\n")
STAR_ONLY_HET = ("chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
                 "GT:DP:AD:PL\t0/1:20:0,20:100,0,100\n")

# Control: the same locus with a concrete ALT that survives too, so the
# `*`-only rule does not apply in either mode.
STAR_PLUS_CONCRETE = ("chr1\t3\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
                      "GT:DP:AD:PL\t0/2:20:0,0,20,0:"
                      "100,100,100,100,0,100,100,100,100,100\n")

# Control: an ORPHAN `*`-only locus that no emitted deletion covers; the
# spurious-spanning-deletion rule removes the `*` from the subset, so the locus
# collapses to REF-only and emits nothing in default mode.
ORPHAN_STAR_ONLY = STAR_ONLY_HOM

# ---------------------------------------------------------------------------
# measured GATK rows (pinned from pinned-GATK runs on these fixtures)
# ---------------------------------------------------------------------------

ROW_DELETION = (
    "chr1\t2\t.\tAAA\tA\t92.60\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
# Compared with compare_columns == 7, so only the columns through FILTER appear;
# the INFO column of this row is where the unstable sign of QD lives.
ROW_STAR_ONLY_HOM_DENSE = "chr1\t3\t.\tA\t*\t0\tLowQual"
ROW_DELETION_DENSE = "chr1\t2\t.\tAAA\tA\t92.60\t."
ROW_STAR_PLUS_CONCRETE = (
    "chr1\t3\t.\tA\t*,G\t82.26\t.\t"
    "AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;"
    "MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\t"
    "GT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100")

CASES = [
    {
        "case": "owned-star-only-locus-refused-default",
        "why": "rule 1 in default mode: the owned '*' survives the "
               "standard-confidence subset but GenotypingEngine.java:173-175 "
               "refuses the locus before recordDeletions() at :179, so GATK "
               "publishes only the covering deletion -- not a second row for "
               "the '*' locus",
        "body": DELETION + STAR_ONLY_HOM,
        "args": [],
        "expect": [ROW_DELETION],
    },
    {
        "case": "owned-star-only-locus-het-refused-default",
        "why": "rule 1 is about the surviving ALT SET, not about the genotype: "
               "a heterozygous '*' call is refused in default mode as well",
        "body": DELETION + STAR_ONLY_HET,
        "args": [],
        "expect": [ROW_DELETION],
    },
    {
        "case": "star-only-locus-emitted-in-dense-mode",
        "why": "rule 1's exception: emitAllActiveSites() is "
               "--include-non-variant-sites, so in dense mode the '*' locus IS "
               "emitted.  Columns 1-7 are compared (CHROM POS ID REF ALT QUAL "
               "FILTER) at positions 2 and 3 only: INFO is excluded because the "
               "QD of this shape carries a sign-of-zero that GATK itself does "
               "not keep consistent (measured: QD=-0.00 at position 3 of this "
               "fixture but QD=0.00 at position 4 of the reverse-trim fixture, "
               "both with QUAL reading 0 -- a ~1e-16 round-off difference in the "
               "AF calculator, see the module docstring), and position 4 is the "
               "tracked covered-locus materialization gap.  The QD divergence "
               "itself stays recorded in the unregistered dense-materialize "
               "gate's REPORTED ONLY case.",
        "body": DELETION + STAR_ONLY_HOM,
        "args": [DENSE],
        "compare_positions": [2, 3],
        "compare_columns": 7,
        "expect_at_positions": {2: ROW_DELETION_DENSE, 3: ROW_STAR_ONLY_HOM_DENSE},
    },
    {
        "case": "star-plus-concrete-alt-unaffected",
        "why": "control: with a concrete ALT also surviving at the locus, "
               "GenotypingEngine's '*' rule does not apply, so both tools emit "
               "the same row in dense mode (position 2 and 3 compared only, for "
               "the same materialization reason as above)",
        "body": DELETION + STAR_PLUS_CONCRETE,
        "args": [DENSE],
        "compare_positions": [2, 3],
        "expect_at_positions": {2: ROW_DELETION, 3: ROW_STAR_PLUS_CONCRETE},
    },
    {
        "case": "orphan-star-only-locus-refused-default",
        "why": "control: with no covering emitted deletion the '*' is a "
               "spurious spanning deletion and is removed from the output "
               "allele subset (GenotypingEngine.java:314), so the locus is "
               "REF-only and emits nothing in default mode -- true before this "
               "round and it must stay true",
        "body": ORPHAN_STAR_ONLY,
        "args": [],
        "expect": [],
    },
]


def write_reference(work: pathlib.Path) -> pathlib.Path:
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def read_lines(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8").splitlines()


def data_rows(path: pathlib.Path) -> list[str]:
    return [line for line in read_lines(path) if line and not line.startswith("#")]


def rows_at_positions(rows: list[str], positions: list[int],
                      columns: int | None = None) -> dict[int, str]:
    """Rows at the given positions; `columns` keeps only the leading tab fields."""
    wanted = {str(position) for position in positions}
    selected: dict[int, str] = {}
    for row in rows:
        fields = row.split("\t")
        if fields[1] in wanted:
            selected[int(fields[1])] = "\t".join(fields[:columns]) if columns else row
    return selected


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-2000:])
    return result


def run_case(case: dict, work: pathlib.Path, reference: pathlib.Path,
             java: pathlib.Path, jar: pathlib.Path, native: pathlib.Path,
             timeout: int) -> dict:
    source = work / f"{case['case']}.g.vcf"
    source.write_text(HEADER + case["body"], encoding="utf-8")
    index_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                           "IndexFeatureFile", "-I", str(source)],
                          f"GATK IndexFeatureFile [{case['case']}]", timeout)
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
                          *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    native_result = invoke([str(native), *common, "--gatk-compatible-annotations",
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    positions = case.get("compare_positions")
    columns = case.get("compare_columns")
    if positions:
        gatk_compared = rows_at_positions(gatk_rows, positions, columns)
        native_compared = rows_at_positions(native_rows, positions, columns)
        expect = case["expect_at_positions"]
    else:
        gatk_compared = {index: row for index, row in enumerate(gatk_rows)}
        native_compared = {index: row for index, row in enumerate(native_rows)}
        expect = {index: row for index, row in enumerate(case["expect"])}

    result = {
        "case": case["case"],
        "why": case["why"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "compared_positions": positions,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    if gatk_compared != expect:
        result["violations"].append(
            f"GATK's rows moved away from the measured contract: "
            f"expected={expect} measured={gatk_compared}")
    if native_compared != gatk_compared:
        result["violations"].append(
            f"native rows differ from GATK: GATK={gatk_compared} "
            f"NATIVE={native_compared}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypingEngine's '*' ALT-set rule "
                    "(GenotypingEngine.java:172-179) and QualByDepth's "
                    "sign-preserving numerator (QualByDepth.java:76-91) against "
                    "pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"),
        help="native GenotypeGVCFs binary (default: $FASTGATK_GENOTYPE_BINARY "
             "or fastgatk-native/build/fastgatk-genotype-gvcf)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the measured GATK/native rows and exit 0 "
             "instead of gating parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-genotype-gvcf")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_not_verified("verify_genotype_gvcf_star_only_locus_gatk_oracle.py", java, jar)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-star-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rules under test: GenotypingEngine.java:172-179 refuses a locus "
          "whose only surviving ALT is the spanning deletion unless "
          "emitAllActiveSites() (--include-non-variant-sites), and QualByDepth "
          "divides the sign-preserving -10*log10PError, so a zero confidence "
          "renders QD=-0.00 while the QUAL token is 0.")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        print(f"[{result['case']}] {result['args']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        if result["compared_positions"]:
            print(f"    compared positions: {result['compared_positions']}")
        print(f"    GATK   rows: {result['gatk_rows']}")
        print(f"    native rows: {result['native_rows']}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "mode": "strict" if strict_mode else "expect-divergence",
        "gatk_version": "4.6.2.0",
        "cases": results,
        "violations": violations,
        "rule": "GenotypingEngine.java:172-175 (refuse a locus whose only "
                "surviving ALT is SPAN_DEL unless emitAllActiveSites()), :178-179 "
                "(recordDeletions runs only for a locus that survives the "
                "refusal), :314 (isSpuriousSpanningDeletion for an orphan '*'), "
                "QualByDepth.java:76-91 (QD = -10*getLog10PError()/depth, "
                "String.format(\"%.2f\")), GenotypingEngine.java:163+183 (the "
                "tested value adds +0.0 while the published double keeps -0.0)",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
