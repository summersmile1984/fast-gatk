#!/usr/bin/env python3
"""Strict pinned-GATK oracle for GenotypeGVCFs' merged-depth precondition.

The rule under test
-------------------
``GenotypeGVCFsEngine.regenotypeVC()`` enters its regenotyping block -- the one
that contains the reverse trim at ``:167`` and every site annotation -- only for
a record that is a VARIANT with a positive merged depth::

    if ( originalVC.isVariant() && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 0 ) {
        ... calculateGenotypes / finalizeAnnotations / reverseTrimAlleles ...
    } else {
        result = originalVC;                     // <- passthrough (:174)
    }
    ...
    if (result.isPolymorphicInSamples() && result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 0) {
        ... annotate and return ...
    } else if (includeNonVariants) {
        ... hom-ref cleanup and return ...
    } else {
        return null;                             // <- no record at all (:198)
    }

(``gatk-source/.../walkers/GenotypeGVCFsEngine.java:157-199``.)

The depth it reads is the depth of the MERGED record, and
``ReferenceConfidenceVariantContextMerger.calculateVCDepth()`` (``:352-360``)
takes ``INFO/DP`` when the input record carries the key, and otherwise sums the
per-sample best depth (``MIN_DP`` when present, else ``DP``) and publishes the
result through ``mergeAttributes()`` (``:382``, only when ``depth > 0``).  So an
input that explicitly declares ``INFO/DP=0`` while its genotypes carry real
depth makes the merged depth zero, the condition false, and -- without
``--include-non-variant-sites`` -- the locus emit nothing at all, even though
the genotypes hold a confident call.

The property pinned here is therefore a *reach* property of the whole
regenotyping block, not of the trim: the trimmed alleles of a locus that GATK
never regenotypes must not be published either.

Scope
-----
Only the default (non dense) arm is gated: there GATK emits no record, which is
a bounded, verifiable expectation.  In dense mode GATK instead passes the merged
record through unchanged (``GT:AD`` hom-ref no-call with ``INFO/DP=0``) and
materializes the covered positions as ``ALT='.'`` rows; that arm is entangled
with the known unfixed covered-locus materialization gap and is deliberately
NOT gated here (measured rows are quoted in the case `why` fields for the next
round).

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

GT_FIELDS = "GT:DP:AD:PL"
HET_STRONG = "0/1:20:0,20,0:100,0,100,100,100,100"
DENSE = "--include-non-variant-sites"


def variant_record(info: str) -> str:
    return (f"chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\t{info}\t"
            f"{GT_FIELDS}\t{HET_STRONG}\n")


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

# Merged depth 0 while the genotypes carry depth 20 and a confident het call:
# GATK's regenotypeVC() takes the passthrough arm and the locus emits nothing.
DP_ZERO_VARIANT = variant_record("DP=0")

# Boundary control: the same locus with INFO/DP=1 clears the > 0 test and is
# regenotyped (and therefore reverse-trimmed).
DP_ONE_VARIANT = variant_record("DP=1")

# Same shape as the DP=0 case but homozygous-reference likelihoods, so the check
# that emits nothing cannot be satisfied merely by a pruned ALT.
DP_ZERO_HOM_REF = (
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=0\t"
    "GT:DP:AD:PL\t0/0:20:20,0,0:0,100,100,100,100,100\n")

# Control for calculateVCDepth()'s else branch: no INFO/DP key and zero
# FORMAT/DP, so the merged depth is 0 and the locus must emit nothing as well.
NO_INFO_DP_ZERO_FORMAT_DEPTH = (
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\t.\t"
    "GT:DP:AD:PL\t0/1:0:0,0,0:100,0,100,100,100,100\n")

# INFO/DP wins over the per-sample MIN_DP in calculateVCDepth() (:353-355), so a
# record whose MIN_DP is 20 but whose INFO/DP is 0 still emits nothing.
HEADER_MIN_DP = HEADER.replace(
    "##FORMAT=<ID=DP",
    '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP">\n'
    "##FORMAT=<ID=DP")
INFO_DP_ZERO_WITH_MIN_DP = (
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=0\t"
    "GT:DP:MIN_DP:AD:PL\t0/1:20:20:0,20,0:100,0,100,100,100,100\n")

# ---------------------------------------------------------------------------
# measured GATK rows (pinned from pinned-GATK runs on these fixtures)
# ---------------------------------------------------------------------------

# Note QD is 4.63, not 92.64: GATK's QD uses the FORMAT AD/DP-derived depth, so
# an INFO/DP of 1 only affects the published INFO/DP column.
ROW_DP_ONE = (
    "chr1\t2\t.\tAAA\tAAC\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=1;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")

CASES = [
    {
        "case": "info-dp-zero-variant-emits-nothing",
        "why": "GenotypeGVCFsEngine.java:160 requires the merged INFO/DP to be "
               "positive before the record is regenotyped, and :181 requires it "
               "again before anything is annotated or emitted; with INFO/DP=0 "
               "and no --include-non-variant-sites the traversal returns null at "
               ":198, so this locus must produce no record at all even though "
               "the genotypes hold a confident 0/1 call.  Measured with dense "
               "mode on the same input, GATK instead passes the merged record "
               "through as 'chr1 2 . AAAA AACA . . DP=0' with FORMAT 'GT:AD' "
               "'./.:0,20' and materializes positions 3-5 as ALT='.' rows -- "
               "that arm is out of scope here (covered-locus materialization).",
        "body": DP_ZERO_VARIANT,
        "args": [],
        "expect": [],
    },
    {
        "case": "info-dp-zero-hom-ref-emits-nothing",
        "why": "the same reach property must not depend on the ALT surviving: "
               "the passthrough arm is chosen before any allele subsetting",
        "body": DP_ZERO_HOM_REF,
        "args": [],
        "expect": [],
    },
    {
        "case": "absent-info-dp-and-zero-format-depth-emits-nothing",
        "why": "calculateVCDepth()'s else branch "
               "(ReferenceConfidenceVariantContextMerger.java:352-360): without "
               "INFO/DP the merged depth is the sum of the per-sample best "
               "depth, here 0, so the same precondition rejects the locus",
        "body": NO_INFO_DP_ZERO_FORMAT_DEPTH,
        "args": [],
        "expect": [],
    },
    {
        "case": "info-dp-zero-wins-over-min-dp-emits-nothing",
        "why": "calculateVCDepth() returns INFO/DP whenever the key is present "
               "(:353-355) and never falls back to the genotypes' MIN_DP, so a "
               "record with INFO/DP=0 and MIN_DP=20 is still not regenotyped",
        "body": INFO_DP_ZERO_WITH_MIN_DP,
        "header": HEADER_MIN_DP,
        "args": [],
        "expect": [],
    },
    {
        "case": "info-dp-one-boundary-is-regenotyped",
        "why": "boundary control for the strict '> 0' test: INFO/DP=1 must still "
               "be regenotyped, trimmed (AAAA/AACA -> AAA/AAC) and emitted, with "
               "INFO/DP published as the merged depth 1",
        "body": DP_ONE_VARIANT,
        "args": [],
        "expect": [ROW_DP_ONE],
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
    result = {
        "case": case["case"],
        "why": case["why"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    if gatk_rows != case["expect"]:
        result["violations"].append(
            f"GATK's rows moved away from the measured contract: "
            f"expected={case['expect']} measured={gatk_rows}")
    if native_rows != gatk_rows:
        result["violations"].append(
            f"native rows differ from GATK: GATK={gatk_rows} NATIVE={native_rows}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypeGVCFs' merged INFO/DP > 0 "
                    "precondition (GenotypeGVCFsEngine.java:157-199) against "
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
        oracle_not_verified("verify_genotype_gvcf_depth_gate_gatk_oracle.py", java, jar)
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-depth-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: GenotypeGVCFsEngine.regenotypeVC() regenotypes "
          "only a variant record whose MERGED INFO/DP is positive (:160), and "
          ":181 applies the same depth test before annotating or emitting; "
          "ReferenceConfidenceVariantContextMerger.calculateVCDepth() (:352-360) "
          "takes INFO/DP when present and otherwise sums the per-sample best "
          "depth (MIN_DP else DP).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        print(f"[{result['case']}] {result['args']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
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
        "rule": "GenotypeGVCFsEngine.java:157-174 (isVariant() && INFO/DP > 0 "
                "guards the regenotyping block; else result = originalVC), "
                ":181-199 (the same depth test gates annotation and emission; "
                "otherwise null unless includeNonVariants), "
                "ReferenceConfidenceVariantContextMerger.java:352-360 "
                "(calculateVCDepth) and :382 (DP published only when depth > 0)",
        "scope": "default (non dense) arm only; the dense passthrough shape is "
                 "tracked with the covered-locus materialization gap",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
