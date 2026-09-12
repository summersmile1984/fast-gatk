#!/usr/bin/env python3
"""Strict pinned-GATK oracle for GenotypeGVCFs' own ``FILTER=LowQual``.

The rule under test
-------------------
GATK 4.6.2.0 ``GenotypingEngine.calculateGenotypes()`` recomputes a phred-scaled
site confidence and applies its own filter when that confidence fails the call
threshold::

    final double log10Confidence =
                !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                        ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;
    final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
    ...
    final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(),
            vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles);
    builder.log10PError(log10Confidence);
    if ( ! passesCallThreshold(phredScaledConfidence) ) {
        builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
    }

(lines 158-163 and 181-186 of
``gatk-source/.../walkers/genotyper/GenotypingEngine.java``), with

    protected final boolean passesCallThreshold(final double conf) {
        return conf >= configuration.genotypeArgs.standardConfidenceForCalling;
    }

(``:430-432``) and ``LOW_QUAL_FILTER_NAME = "LowQual"``
(``GATKVCFConstants.java:179``).  The header declaration is unconditional
(``GenotypeGVCFsEngine.java:416``).

Three properties of that rule are pinned here, because each of them is
separately falsifiable and a plausible implementation can get each of them
wrong:

1. **The value tested is the recomputed confidence, not the source record's
   FILTER and not the printed QUAL token.**  ``builder.log10PError(log10Confidence)``
   publishes QUAL from the SAME double, but htsjdk renders QUAL with two
   decimals and a trailing ``.00`` stripped, so the printed token can be at or
   above the cutoff while the filter is applied.  The ``near-cutoff-*`` case
   pair brackets a real measurement of that: on the weak-locus fixture GATK
   emits ``QUAL=23.14`` and applies ``LowQual`` at
   ``--standard-min-confidence-threshold-for-calling 23.1358`` but not at
   ``23.135``, so the tested confidence lies in ``[23.135, 23.1358)`` -- below
   the token it prints.  A repair that tests the rounded QUAL fails
   ``near-cutoff-just-below``.
2. **The filter is site-level and is decided after the output-allele subset.**
   ``calculateOutputAlleleSubset()`` runs at ``:155`` and the surviving ALT set
   only feeds ``builder`` at ``:181``; ``builder.filter()`` at ``:185`` has no
   per-allele argument.  GenotypeGVCFs never calls
   ``AlleleFilterUtils.addAlleleAndSiteFilters()`` (whose per-allele
   intersection at ``AlleleFilterUtils.java:115-120`` belongs to
   VariantFiltration/Mutect2), and ``--invalidate-previous-filters`` is not a
   GenotypeGVCFs option at all -- measured: GATK exits 1 with
   "is not a recognized option" (case
   ``invalidate-previous-filters-unsupported``).
3. **Rows the engine did not regenotype carry no filter.**  A reference block
   materialized by ``--include-non-variant-sites`` is a passthrough row with a
   missing QUAL (case ``reference-block-row-not-filtered``), and a locus whose
   confidence is infinite never fails the threshold (case
   ``infinite-confidence-never-filtered``).

Reach, measured
---------------
Without ``--include-non-variant-sites`` the filter is unreachable in practice:
``passesEmitThreshold()`` (``:425-427``) requires ``passesCallThreshold()``
before a best-guess-reference locus is emitted at all, so every measured
non-dense fixture emits either a passing call or no record at all
(``default-mode-emits-no-low-confidence-row``).  The divergence this gate pins
is therefore confined to the dense path, which is why the fix is confined to
the GATK-compatibility FILTER-writing path of that tool.

Usage / exit status
-------------------
``--native`` / ``$FASTGATK_GENOTYPE_BINARY`` select the native binary
(``$FASTGATK_NATIVE_BUILD``/``fastgatk-native/build`` supplies the default), the
pinned GATK jar and the bundled JDK17 are read from ``third_party/``, all
scratch lives in a ``tempfile.TemporaryDirectory``, and the run ends with a
single-line JSON status payload.  Exit status: 0 when every gated case passes,
non-zero otherwise.  ``--expect-divergence`` turns the run into a diagnostic
that always exits 0.
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

# --------------------------------------------------------------------------
# Fixtures.  Every body is written under a 100 bp all-``A`` chr1 so the REF
# allele of each record is valid.  The header is the same one the sibling
# spanning-deletion oracle uses, so the two gates describe the same input
# surface; HEADER_GP adds the posterior annotation that
# --use-posteriors-to-calculate-qual consumes.
# --------------------------------------------------------------------------
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

HEADER_GP = HEADER.replace(
    "##FORMAT=<ID=PL",
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior">\n'
    "##FORMAT=<ID=PL")

# A locus whose best guess is the reference: the ALT is pruned by the
# standard-confidence subset, so the dense materialization is the REF-only
# no-call row GATK measures at 23.14.
WEAK_LOCUS_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:19,1:0,0,40\n"
)

# The same locus with FORMAT/GP present, so --use-posteriors-to-calculate-qual
# has an input to consume.  Measured: GATK still writes 23.14 / LowQual.
WEAK_LOCUS_GP_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL:GP\t0/1:20:19,1:0,0,40:0.999,0.001,0.000\n"
)

# A confidently reference locus: the same monomorphic materialization with a
# finite QUAL well above the default cutoff, so no filter is applied.
STRONG_REF_MONOMORPHIC_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20:0,100,100\n"
)

# Two records: the upstream deletion is genotyped normally, the downstream
# locus is materialized REF-only with an INFINITE confidence
# (MathUtils.log10OneMinusPow10(0.0) -> log10PError = -Infinity).  QA'ing the
# token as a string would call it low quality; the rule is the double test.
STAR_ONLY_COVERED_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# A pure reference-confidence block: --include-non-variant-sites materializes
# it without ever calling the genotyping engine, so its QUAL stays missing.
REFERENCE_BLOCK_RECORD = (
    "chr1\t2\t.\tA\t<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:GQ:PL\t0/0:20:99:0,100,100\n"
)

# Two ALT alleles with very different support, both pruned at the cutoff: the
# site-level decision must not be per-allele.
SPLIT_ALT_RECORD = (
    "chr1\t2\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,19,1:0,0,0,100,100,100,100,100,100\n"
)

# A surviving concrete ALT whose site confidence is below the cutoff.  At
# --standard-min-confidence-threshold-for-calling 100 the only ALT is pruned
# too, and GATK publishes QUAL 0 with LowQual.
G_PLAUSIBLE_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tLowQual\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# --------------------------------------------------------------------------
# Measured GATK 4.6.2.0 rows (byte for byte).  The gate asserts GATK against
# these literals AND native against GATK, so a change on either side is caught.
# --------------------------------------------------------------------------
GATK_WEAK_LOWQUAL_ROW = ("chr1\t2\t.\tA\t.\t23.14\tLowQual\t"
                         "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_WEAK_PASS_ROW = ("chr1\t2\t.\tA\t.\t23.14\t.\t"
                      "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_STRONG_REF_ROW = ("chr1\t2\t.\tA\t.\t127.78\t.\t"
                       "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_STAR_UPSTREAM_ROW = (
    "chr1\t2\t.\tAA\tA\t92.60\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
GATK_STAR_DENSE_ROW = ("chr1\t3\t.\tA\t.\tInfinity\t.\t"
                       "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_STAR_UPSTREAM_PRUNED_ROW = ("chr1\t2\t.\tAA\t.\t0\tLowQual\t"
                                 "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_REFERENCE_BLOCK_ROW = ("chr1\t2\t.\tA\t.\t.\t.\tDP=20\tGT:DP:RGQ\t0/0:20:99")
GATK_SPLIT_ROW = ("chr1\t2\t.\tA\t.\t25.59\tLowQual\t"
                  "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_G_PRUNED_ROW = ("chr1\t2\t.\tA\t.\t0\tLowQual\t"
                     "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")

DENSE = "--include-non-variant-sites"
THRESHOLD = "--standard-min-confidence-threshold-for-calling"

CASES: list[dict] = [
    {
        "case": "below-threshold-dense-monomorphic",
        "why": "the reported divergence: the recomputed confidence is 23.14, "
               "below the default cutoff of 30, so GenotypingEngine.java:184-186 "
               "applies LowQual to the REF-only record the dense mode "
               "materializes",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE],
        "expect": [GATK_WEAK_LOWQUAL_ROW],
    },
    {
        "case": "above-threshold-dense-monomorphic",
        "why": "control on the same branch: the identical materialization with a "
               "confidence of 127.78 passes the cutoff and stays unfiltered, so "
               "the filter is not a blanket rule for monomorphic rows",
        "body": STRONG_REF_MONOMORPHIC_RECORD,
        "args": [DENSE],
        "expect": [GATK_STRONG_REF_ROW],
    },
    {
        "case": "near-cutoff-just-above",
        "why": "the lower half of the raw-double bracket: at a cutoff of 23.135 "
               "the weak locus is NOT filtered, so the tested confidence is "
               ">= 23.135",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE, THRESHOLD, "23.135"],
        "expect": [GATK_WEAK_PASS_ROW],
    },
    {
        "case": "near-cutoff-just-below",
        "why": "the upper half of the bracket: at 23.1358 the same locus IS "
               "filtered, so the tested confidence is < 23.1358 even though the "
               "printed QUAL token is 23.14 -- the rule is not a test on the "
               "rendered QUAL",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE, THRESHOLD, "23.1358"],
        "expect": [GATK_WEAK_LOWQUAL_ROW],
    },
    {
        "case": "threshold-lowered-below-confidence",
        "why": "the cutoff also decides the output-allele subset, so this pins "
               "the whole branch at a lower threshold: the locus becomes a "
               "filtered-in call and is not LowQual",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE, THRESHOLD, "10"],
        "expect": [GATK_WEAK_PASS_ROW],
    },
    {
        "case": "threshold-raised-above-confidence",
        "why": "the same fixture at a cutoff of 100 keeps the site LowQual, so "
               "the filter follows the option rather than a hard-coded 30",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE, THRESHOLD, "100"],
        "expect": [GATK_WEAK_LOWQUAL_ROW],
    },
    {
        "case": "surviving-alt-below-threshold",
        "why": "per-allele interaction: at a cutoff of 100 every ALT is pruned "
               "as well, and the site filter is still applied on top of that "
               "pruning (the two decisions are independent)",
        "body": G_PLAUSIBLE_RECORD,
        "args": [DENSE, THRESHOLD, "100"],
        "expect": [GATK_G_PRUNED_ROW],
    },
    {
        "case": "two-alt-locus-below-threshold",
        "why": "a second, independent below-threshold fixture (confidence 25.59) "
               "so the gate does not rest on one calibrated PL row",
        "body": SPLIT_ALT_RECORD,
        "args": [DENSE],
        "expect": [GATK_SPLIT_ROW],
    },
    {
        "case": "infinite-confidence-never-filtered",
        "why": "QUAL=Infinity is a real positive infinity, not a sentinel: "
               "log10ProbVariantPresent() is -Infinity so the confidence is "
               "+Infinity, which passes ANY finite cutoff.  The unrelated "
               "upstream locus in the same run IS filtered once the cutoff is "
               "raised above its confidence, so the case also shows the filter "
               "is per-record",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": [DENSE, THRESHOLD, "1000000000"],
        "expect": [GATK_STAR_UPSTREAM_PRUNED_ROW, GATK_STAR_DENSE_ROW],
    },
    {
        "case": "reference-block-row-not-filtered",
        "why": "a passthrough reference-confidence block is materialized without "
               "running the genotyping engine at all, so it has no recomputed "
               "confidence and no filter; its QUAL is missing, not 0",
        "body": REFERENCE_BLOCK_RECORD,
        "args": [DENSE],
        "expect": [GATK_REFERENCE_BLOCK_ROW],
    },
    {
        "case": "posterior-qual-option-keeps-call-confidence",
        "why": "GenotypingEngine applies the filter at :184-186 BEFORE the "
               "optional posterior QUAL update at :192-199, so the decision is "
               "always made on the pre-posterior confidence; measured on a "
               "GP-bearing fixture with --use-posteriors-to-calculate-qual the "
               "row stays 23.14 / LowQual",
        "body": WEAK_LOCUS_GP_RECORD,
        "header": HEADER_GP,
        "args": [DENSE, "--use-posteriors-to-calculate-qual"],
        "expect": [GATK_WEAK_LOWQUAL_ROW],
    },
    {
        "case": "default-mode-emits-no-low-confidence-row",
        "why": "reach: without --include-non-variant-sites, "
               "passesEmitThreshold() (:425-427) requires the call threshold "
               "before a reference locus is emitted, so a below-threshold locus "
               "produces no record at all and the FILTER column never arises",
        "body": WEAK_LOCUS_RECORD,
        "args": [],
        "expect": [],
    },
    {
        "case": "default-mode-strong-variant-unfiltered",
        "why": "the same reach statement on an emitted call: a confidently "
               "genotyped ALT is emitted with FILTER='.' and no LowQual",
        "body": G_PLAUSIBLE_RECORD,
        "args": [],
        "expect": [
            "chr1\t2\t.\tA\tG\t82.26\t.\t"
            "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
            "GT:AD:DP:PL\t0/1:0,20:20:0,0,0"],
    },
    {
        "case": "invalidate-previous-filters-unsupported",
        "why": "the per-allele filter intersection of AlleleFilterUtils.java "
               "(and the --invalidate-previous-filters flag that drives it) is "
               "NOT part of GenotypeGVCFs: the option belongs to "
               "VariantFiltration/FilterVariantTranches, and pinned GATK "
               "rejects it here with exit 1 and 'is not a recognized option'.  "
               "Native must reject it too, so the two tools stay comparable on "
               "the arguments they share",
        "body": WEAK_LOCUS_RECORD,
        "args": [DENSE, "--invalidate-previous-filters"],
        "expect": [],
        "expect_gatk_failure": True,
        "expect_native_failure": True,
    },
]

LOW_QUAL_FILTER_LINE = '##FILTER=<ID=LowQual,Description="Low quality">'


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
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


def filter_header_lines(path: pathlib.Path) -> list[str]:
    return [line for line in read_lines(path) if line.startswith("##FILTER=")]


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
    source.write_text(case.get("header", HEADER) + case["body"], encoding="utf-8")
    # --include-non-variant-sites switches GATK to a group-by-locus traversal
    # that requires an index; indexing both sides' identical input keeps the two
    # invocations comparable.
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
    # Native's default profile is a non-GATK diagnostic output, so the
    # GATK-compatibility writer is the only byte-comparable surface; the
    # semantic arguments are identical on both sides.
    native_result = invoke([str(native), *common, "--gatk-compatible-annotations",
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    gatk_filters = filter_header_lines(gatk_out) if gatk_out.exists() else []
    native_filters = filter_header_lines(native_out) if native_out.exists() else []

    result = {
        "case": case["case"],
        "why": case["why"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "gatk_qual_filter": [{"qual": row.split("\t")[5], "filter": row.split("\t")[6]}
                             for row in gatk_rows],
        "native_qual_filter": [{"qual": row.split("\t")[5], "filter": row.split("\t")[6]}
                               for row in native_rows],
        "gatk_filter_header_lines": gatk_filters,
        "native_filter_header_lines": native_filters,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    if case.get("expect_gatk_failure"):
        if gatk_result.returncode == 0:
            result["violations"].append(
                "GATK accepted the invocation; the measured contract is that it "
                "rejects it (exit 1, 'is not a recognized option')")
    elif gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if case.get("expect_native_failure"):
        if native_result.returncode == 0:
            result["violations"].append(
                "native accepted the invocation but GATK rejects it")
    elif native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    # GATK's own behaviour is pinned from the measurement, not re-derived.
    if not case.get("expect_gatk_failure") and gatk_rows != case["expect"]:
        result["violations"].append(
            f"GATK's rows moved away from the measured contract: "
            f"expected={case['expect']} measured={gatk_rows}")
    if not case.get("expect_native_failure") and native_rows != gatk_rows:
        result["violations"].append(
            f"native rows differ from GATK: "
            f"GATK={gatk_rows} NATIVE={native_rows}")
    # The FILTER value is only well formed when its declaration is in the
    # output header, so the group is compared as well.
    if not case.get("expect_gatk_failure") and gatk_filters != native_filters:
        result["violations"].append(
            f"##FILTER header lines are not byte-identical: "
            f"GATK={gatk_filters} NATIVE={native_filters}")
    if not case.get("expect_gatk_failure") and LOW_QUAL_FILTER_LINE not in gatk_filters:
        result["violations"].append(
            f"GATK's own output header no longer declares {LOW_QUAL_FILTER_LINE}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypeGVCFs' own FILTER=LowQual rule "
                    "(GenotypingEngine.java:158-163, :181-186, :430-432) against "
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
        oracle_not_verified("verify_genotype_gvcf_lowqual_gatk_oracle.py", java, jar)
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
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-genotype-lowqual-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: GenotypingEngine.calculateGenotypes() recomputes "
          "the site confidence (:158-163) and applies its own LowQual when "
          "!passesCallThreshold(phredScaledConfidence) (:184-186, :430-432) on the "
          "rebuilt REF/ALT-subset record (:181); the filter is site-level and the "
          "published QUAL is the same double, rounded to two decimals by htsjdk.")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        print(f"[{result['case']}] {result['args']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    GATK   QUAL/FILTER: {result['gatk_qual_filter']}")
        print(f"    native QUAL/FILTER: {result['native_qual_filter']}")
        if result["gatk_rows"] != result["native_rows"]:
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
        "rule": "GenotypingEngine.java:158-163 (log10Confidence / "
                "phredScaledConfidence), :181-186 (VariantContextBuilder + "
                "builder.filter(LOW_QUAL_FILTER_NAME)), :430-432 "
                "(passesCallThreshold: conf >= standardConfidenceForCalling), "
                "GATKVCFConstants.java:179 (LOW_QUAL_FILTER_NAME = LowQual), "
                "GATKVCFHeaderLines.java:89 + GenotypeGVCFsEngine.java:416 "
                "(unconditional header declaration)",
        "measured_bracket": "the weak-locus confidence is in [23.135, 23.1358) "
                            "while htsjdk prints QUAL 23.14, so the rule is a "
                            "test on the raw double and not on the QUAL token",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
