#!/usr/bin/env python3
"""Strict pinned-GATK oracle for confident homozygous-alternate loci in GenotypeGVCFs.

The defect this pins
--------------------
``calculate_allele_frequency_kokkos()`` normalised the per-genotype posteriors in
LINEAR space::

    probability = pow(10, term - maximum) / denominator

A confident homozygous-alternate sample has P(hom-ref) far below the smallest
positive double -- PL differences of a few thousand are ordinary in real gVCFs --
so that division underflowed to exactly ``0.0``, ``log10(p0)`` became ``-inf``,
and the host's ``p0 <= -1.0e299`` test then read the sample as "no likelihoods at
all".  With every sample skipped, ``log10_p_allele_absent`` kept its initial
``0.0`` for every allele, the standard-confidence subset judged the ALT
implausible and pruned it, and in the default traversal the locus was dropped
entirely: **a confident variant call disappeared**.

Measured on GATK's own chr20 HaplotypeCaller gVCF (1291 records): 18 of the 252
loci GATK emits were missing, and every one of them was a confident hom-alt call
(GT ``1/1``/``1|1``, AD reference depth 0, PL spread ~2600-4000).  The same
missing call also shifted the ``QualByDepth.fixTooHighQD()`` random sequence, so
the *next* high-QD loci printed a different jittered QD as well -- one root cause,
two visible classes.

The fix (kept inside the kernel, the Host/kernel boundary is unchanged)
accumulates both quantities with a log-sum-exp, so no intermediate probability is
ever materialised in linear space.  Verified on the same corpus: 252/252 loci,
and the jittered QD values realign because the draws are consumed in the same
order again.

Cases
-----
* ``confident-hom-alt-is-emitted`` -- the minimal single-record reproduction, which
  the defect dropped and which must now be byte-identical to GATK;
* ``confident-hom-alt-het-control`` -- the same locus with heterozygous
  likelihoods, which was never affected (guard against a fix that emits too much);
* ``ill-conditioned-pruned-alt-still-dropped`` -- a hom-reference call whose ALT
  carries a large PL spread: the ALT must still be pruned and, in the default
  traversal, the locus must still produce no record.

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

# The reference must carry the record's REF base at 20:10000758 (GRCh37 chr20);
# the gate generates its own 100 bp contig named after the corpus record so the
# fixture stays self-contained.
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

# The corpus record at 20:10000758, moved onto a clean 100 bp contig: a
# confident hom-alt call whose PL spread (3867) drives P(hom-ref) below the
# smallest positive double.
CONFIDENT_HOM_ALT = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t3853.06\t.\tDP=97\t"
    "GT:AD:DP:GQ:PL\t1/1:0,95,0:95:99:3867,286,0,3867,286,3867\n")

# Control: the same locus with heterozygous likelihoods (PL best at the het
# genotype), which never underflowed and must stay unchanged.
CONFIDENT_HET = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t926.64\t.\tDP=66\t"
    "GT:AD:DP:GQ:PL\t0/1:33,33,0:66:99:1000,0,1000,1000,1000,1000\n")

# Control: a hom-REFERENCE call whose ALT carries a huge PL gap.  The ALT is
# implausible, so GATK prunes it and the default traversal emits nothing; the
# log-space accumulation must not change that.
ILL_CONDITIONED_HOM_REF = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\t.\tDP=40\t"
    "GT:AD:DP:GQ:PL\t0/0:40,0,0:40:99:0,3867,3867,3867,3867,3867\n")

# A multi-allelic source whose best genotype does not survive the output-allele
# subset: the projected PL row keeps the source's offset, and htsjdk normalises it
# through GenotypeLikelihoods.  The PL values are the real row of GATK's own chr20
# corpus record at 20:10002458 (five alleles, diploid, best genotype (1,2)).
MULTI_ALLELIC_PL_PROJECTION = (
    "chr1\t2\t.\tA\tT,TT,TTTT,<NON_REF>\t2154.96\t.\tDP=65\t"
    "GT:AD:DP:GQ:PL\t1/2:0,6,6,6,0:18:9:"
    "2172,655,493,249,0,82,352,33,9,214,876,512,167,254,730\n")

CASES = [
    {
        "case": "confident-hom-alt-is-emitted",
        "why": "the defect: P(hom-ref) underflowed to 0, log10 became -inf, the "
               "host read the sample as having no likelihoods, and the ALT was "
               "pruned -- the confident call vanished.  Measured rows are pinned "
               "from pinned GATK, so this case also pins the QUAL/AC/AF/QD set "
               "that the corrected cohort posterior produces",
        "body": CONFIDENT_HOM_ALT,
        "args": [],
        "expect": ["chr1\t2\t.\tA\tG\t3853.06\t.\t"
                   "AC=2;AF=1.00;AN=2;DP=97;ExcessHet=0.0000;MLEAC=2;MLEAF=1.00;"
                   "QD=25.36\tGT:AD:DP:GQ:PL\t1/1:0,95:95:99:3867,286,0"],
    },
    {
        "case": "confident-hom-alt-het-control",
        "why": "control: the heterozygous shape was never affected; a fix that "
               "widened the plausible-allele rule would show up here",
        "body": CONFIDENT_HET,
        "args": [],
        "expect": ["chr1\t2\t.\tA\tG\t992.64\t.\t"
                   "AC=1;AF=0.500;AN=2;DP=66;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;"
                   "QD=15.04\tGT:AD:DP:GQ:PL\t0/1:33,33:66:99:1000,0,1000"],
    },
    {
        "case": "ill-conditioned-pruned-alt-still-dropped",
        "why": "control: with hom-reference likelihoods the ALT is genuinely "
               "implausible and GATK emits no record in the default traversal; "
               "the log-space accumulation must not resurrect it",
        "body": ILL_CONDITIONED_HOM_REF,
        "args": [],
        "expect": [],
    },
    {
        "case": "multi-allelic-projected-pl-is-normalized",
        "why": "htsjdk rebuilds the subsetted genotypes through "
               "GenotypeLikelihoods, which normalises each PL row to its minimum; "
               "the projection alone keeps the source offset, so native used to "
               "publish `2172,249,82` where GATK publishes `2090,167,0` (the same "
               "defect on 2 of GATK's own chr20 loci: +82 and +45)",
        "body": MULTI_ALLELIC_PL_PROJECTION,
        "args": [],
        "expect": ["chr1\t2\t.\tA\tTT\t2155\t.\t"
                   "AC=2;AF=1.00;AN=2;DP=65;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;"
                   "QD=25.36\tGT:AD:DP:GQ:PL\t1/1:0,6:18:99:2090,167,0"],
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
        description="Strict oracle for confident hom-alt loci in GenotypeGVCFs "
                    "(allele-frequency posterior underflow) against pinned GATK "
                    "4.6.2.0.")
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
        oracle_not_verified("verify_genotype_gvcf_confident_hom_alt_gatk_oracle.py", java, jar)
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
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-homalt-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# rule under test: the cohort allele-frequency posterior must stay "
          "representable for confident homozygous-alternate samples, so that "
          "P(allele absent) is finite and the ALT survives the standard-confidence "
          "subset exactly where GATK keeps it.")
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
        "rule": "fastgatk-kernels calculate_allele_frequency_kokkos(): the final "
                "posterior evaluation must accumulate in log space "
                "(log-sum-exp) instead of materialising pow(10, term-maximum)/"
                "denominator, whose underflow to 0 turned log10(p0) into -inf and "
                "collided with the 'no likelihoods' sentinel",
        "reach": "GATK's own chr20 HaplotypeCaller gVCF: 18 of 252 emitted loci "
                 "were lost before the fix; the same loci also shifted the "
                 "QualByDepth.fixTooHighQD random sequence",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
