#!/usr/bin/env python3
"""Strict gate: native GenotypeGVCFs must never SILENTLY DISCARD input.

The defect this gate pins
-------------------------
On a malformed VCF *record* native GenotypeGVCFs exits 0 and writes **zero data
rows** -- including a well-formed record that follows the malformed one.  The
measured shape is a record whose FORMAT column declares four keys while the
sample column carries five values::

    chr1  2  .  A  G,<NON_REF>  .  PASS  DP=20  GT:DP:AD:PL  0/1:20:0,20,0:100,0,100,100,100,100:99
    chr1  5  .  A  G,<NON_REF>  .  PASS  DP=20  GT:DP:AD:PL  0/1:20:0,20,0:100,0,100,100,100,100

Pinned GATK 4.6.2.0 exits 3 and writes no row; native exits 0 and writes no
row *either* -- the record at chr1:5 is lost without a word on stderr beyond
htslib's ``[E::vcf_parse_format_max3]`` log line.

The mechanism (htslib contract, read from ``third_party/htslib-build/htslib-src``)
-------------------------------------------------------------------------------
``bcf_read()`` returns

* ``0``  -- one record was read;
* ``-1`` -- end of input (``vcf.c:3930-3936``, ``hts_getline`` returning EOF);
* ``-2`` -- the record did not parse (``vcf_parse`` returns ``-2`` from its
  ``err:`` label, ``vcf.c:3067``-ish; ``vcf_parse1`` is a macro for
  ``vcf_parse``, ``htslib/vcf.h:274``).

``genotype_gvcf_tool.cpp`` drove its linear traversals with

    while (bcf_read(input, header, record) == 0) process_record(record);

which cannot tell ``-1`` from ``-2``: the first malformed record ends the
traversal as if the file had ended.  The indexed traversal in the *same file*
already distinguishes them (``if (status < -1) throw`` at
``genotype_gvcf_tool.cpp:6505``) -- which is exactly the contract asserted here.

The contract asserted
---------------------
For every malformed fixture in this gate, native must satisfy ONE of:

  (a) exit non-zero, print a diagnostic to stderr, and write NO data row
      (fail loudly -- preferred, and what the fix does); or
  (b) exit 0 and still emit every well-formed record of the input.

What must never happen -- and is what this gate fails on -- is exit 0 with a
well-formed record missing.  ``control-valid-input`` pins the other side: the
same fixture with the malformed record removed must still match pinned GATK row
for row, so a fix that over-rejects valid input fails here too.

Run:  python3 fastgatk-native/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py
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

import oracle_guard  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""

FORMAT = "GT:DP:AD:PL"
GOOD_SAMPLE = "0/1:20:0,20,0:100,0,100,100,100,100"

# The well-formed control row, pinned from the measured GATK 4.6.2.0 run
# (identical on both sides; also the row the malformed fixtures must not lose).
GOOD_ROW = ("chr1\t5\t.\tA\tG\t92.64\t.\t"
            "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;"
            "QD=4.63\tGT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")


def record(pos: str, sample: str, fmt: str = FORMAT,
           alt: str = "G,<NON_REF>") -> str:
    """One data line; ``sample=None`` omits the sample column entirely."""
    columns = ["chr1", pos, ".", "A", alt, ".", "PASS", "DP=20", fmt]
    if sample is not None:
        columns.append(sample)
    return "\t".join(columns) + "\n"


GOOD_RECORD = record("5", GOOD_SAMPLE)

# The same well-formed row one locus earlier (chr1:2), as GATK writes it when
# the malformed record comes *after* it.
GOOD_ROW_AT_2 = GOOD_ROW.replace("\t5\t", "\t2\t")


def bad_first(sample: str, fmt: str = FORMAT) -> str:
    """The malformed record at chr1:2, followed by the well-formed chr1:5 row."""
    return record("2", sample, fmt=fmt) + GOOD_RECORD


CASES = [
    # --- gated: native currently exits 0 and silently loses chr1:5 ----------
    {
        "case": "too-many-format-values",
        "why": "FORMAT declares 4 keys but the sample carries 5 values; htslib "
               "rejects the record (vcf_parse_format_max3) and the traversal "
               "must not treat that as end of input",
        "body": bad_first(GOOD_SAMPLE + ":99"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
    {
        "case": "too-many-format-values-alone",
        "why": "the same malformed record as the only record in the file",
        "body": record("2", GOOD_SAMPLE + ":99"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
    {
        "case": "too-many-format-values-bad-last",
        "why": "the malformed record is the LAST record: GATK has already "
               "written the chr1:2 row when it throws, native reads before it "
               "writes and must not publish a partial file",
        "body": record("2", GOOD_SAMPLE) + record("5", GOOD_SAMPLE + ":99"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 3,
        "gatk_expect": [GOOD_ROW_AT_2],
        "gatk_partial_rows_note": "GATK writes the chr1:2 row and then exits 3; "
                                  "native reads before it writes and must not "
                                  "publish a partial file",
    },
    {
        "case": "empty-sample-column",
        "why": "the sample column is present but empty (columns 0 vs 1, "
               "vcf_parse_format_check7)",
        "body": bad_first(""),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 0,
        "gatk_expect": [GOOD_ROW],
    },
    {
        "case": "dp-non-numeric",
        "why": "FORMAT/DP is not an integer (vcf_parse_format_fill5)",
        "body": bad_first("0/1:abc:0,20,0:100,0,100,100,100,100"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
    {
        "case": "ad-non-numeric",
        "why": "FORMAT/AD holds a non-numeric value "
               "(vcf_parse_format_fill5)",
        "body": bad_first("0/1:20:0,x,0:100,0,100,100,100,100"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 0,
        "gatk_expect": [
            "chr1\t2\t.\tA\tG\t92.64\t.\tAC=1;AF=0.500;AN=2;DP=20;"
            "ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
            "GT:AD:DP:GQ:PL\t0/1:20,0:20:99:100,0,100",
            GOOD_ROW,
        ],
    },
    {
        "case": "pl-non-numeric",
        "why": "FORMAT/PL holds non-integer values "
               "(vcf_parse_format_fill5)",
        "body": bad_first("0/1:20:0,20,0:x,y,z"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 0,
        "gatk_expect": [GOOD_ROW],
    },
    {
        "case": "gt-non-numeric",
        "why": "FORMAT/GT is not an allele index pair "
               "(vcf_parse_format_fill5: could not read GT data)",
        "body": bad_first("x/y:20:0,20,0:100,0,100,100,100,100"),
        "gated": True,
        "fail_loud": True,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
    # --- the control: a fix that rejects valid input fails here --------------
    {
        "case": "control-valid-input",
        "why": "the same fixture with the malformed record removed: native must "
               "still match pinned GATK row for row",
        "body": GOOD_RECORD,
        "gated": True,
        "fail_loud": False,
        "gatk_expected_exit": 0,
        "gatk_expect": [GOOD_ROW],
        "native_expect": [GOOD_ROW],
    },
    # --- reported only: divergence of a different kind, pinned not gated -----
    {
        "case": "missing-sample-column",
        "why": "REPORTED ONLY -- the FORMAT column has no sample columns at "
               "all.  htslib logs vcf_parse_format_empty1 but its caller turns "
               "that -1 into success (vcf.c: vcf_parse_format), so the record "
               "parses with zero samples; native emits a record with no FORMAT "
               "at all where GATK exits 3.  No well-formed record is lost, so "
               "this is not the silent-drop contract -- it is an adjacent "
               "divergence recorded here",
        "body": bad_first(None),
        "gated": False,
        "fail_loud": False,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
        "native_expect": ["chr1\t2\t.\tA\tG\t.\t.\tDP=20", GOOD_ROW],
    },
    {
        "case": "too-few-format-values",
        "why": "REPORTED ONLY -- FORMAT declares 4 keys but the sample carries "
               "3; htslib accepts the record and BOTH tools drop just that site "
               "(GATK warns 'insufficient data ... Site will be skipped')",
        "body": bad_first("0/1:20:0,20,0"),
        "gated": False,
        "fail_loud": False,
        "gatk_expected_exit": 0,
        "gatk_expect": [GOOD_ROW],
        "native_expect": [GOOD_ROW],
    },
    {
        "case": "pl-width-short",
        "why": "REPORTED ONLY -- PL carries 2 values where Number=G wants 6; "
               "GATK exits 3, native pads the row to PL=100,0,. (a wrong-shaped "
               "PL vector, not a dropped record)",
        "body": bad_first("0/1:20:0,20,0:100,0"),
        "gated": False,
        "fail_loud": False,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
        "native_expect": [
            "chr1\t2\t.\tA\tG\t92.64\t.\tAC=1;AF=0.500;AN=2;DP=20;"
            "ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
            "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,.",
            GOOD_ROW,
        ],
    },
    {
        "case": "gt-allele-index-out-of-range",
        "why": "REPORTED ONLY -- GT names allele index 5 while the record has "
               "alleles 0..2; GATK exits 3, native genotypes it as 0/1",
        "body": bad_first("0/5:20:0,20,0:100,0,100,100,100,100"),
        "gated": False,
        "fail_loud": False,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
    {
        "case": "extra-format-key-with-value",
        "why": "REPORTED ONLY -- FORMAT names an undeclared key 'ZZ' "
               "(htslib warns vcf_parse_format_dict2); GATK exits 3, native "
               "emits both rows",
        "body": bad_first("0/1:20:0,20,0:100,0,100,100,100,100:7",
                          fmt="GT:DP:AD:PL:ZZ"),
        "gated": False,
        "fail_loud": False,
        "gatk_expected_exit": 3,
        "gatk_expect": [],
    },
]


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-2000:])
    return result


def data_rows(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


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
    common = ["-R", str(reference), "-V", str(source)]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                          "GenotypeGVCFs", *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    native_result = invoke([str(native), *common,
                            "--gatk-compatible-annotations", "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out)
    native_rows = data_rows(native_out)
    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "fail_loud": case["fail_loud"],
        "input": case["body"],
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_expected_exit": case["gatk_expected_exit"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "native_output_exists": native_out.exists(),
        "native_stderr_tail": native_result.stderr[-600:],
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(
            f"GATK IndexFeatureFile exited {index_result.returncode}")
    # Pinned GATK truth: if GATK moves, the gate must say so rather than
    # silently re-pinning itself.
    if gatk_result.returncode != case["gatk_expected_exit"]:
        result["violations"].append(
            f"GATK moved away from the measured truth: exit={gatk_result.returncode} "
            f"expected={case['gatk_expected_exit']}")
    if gatk_rows != case["gatk_expect"]:
        result["violations"].append(
            f"GATK moved away from the measured truth: rows={gatk_rows} "
            f"expected={case['gatk_expect']}")
    if not case["gated"]:
        expected = case.get("native_expect")
        if expected is not None and native_rows != expected:
            result["violations"].append(
                f"native rows moved: {native_rows} expected={expected}")
        return result

    if case["fail_loud"]:
        # The fix's contract: a malformed record is reported, not swallowed.
        if native_result.returncode == 0:
            result["violations"].append(
                "native exited 0 on a malformed input record: the malformed "
                "record was swallowed instead of reported")
        if native_rows:
            result["violations"].append(
                f"native wrote {len(native_rows)} data row(s) for a failed run "
                f"(no partial output is allowed): {native_rows}")
        if not any("fastgatk-genotype-gvcf:" in line
                   for line in native_result.stderr.splitlines()):
            result["violations"].append(
                "native printed no diagnostic for the malformed record")
    else:
        expected = case.get("native_expect")
        if expected is not None and native_rows != expected:
            result["violations"].append(
                f"native rows moved: {native_rows} expected={expected}")

    # The invariant that holds for EVERY gated case: a well-formed record must
    # never disappear behind a successful exit.
    if (native_result.returncode == 0 and GOOD_ROW not in native_rows
            and case["body"].endswith(GOOD_RECORD)):
        result["violations"].append(
            "native exited 0 but silently dropped the well-formed record at "
            f"chr1:5: rows={native_rows}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict gate: native GenotypeGVCFs must fail loudly on a "
                    "malformed VCF record instead of silently truncating its "
                    "input.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"),
        help="native GenotypeGVCFs binary (default: $FASTGATK_GENOTYPE_BINARY "
             "or fastgatk-native/build/fastgatk-genotype-gvcf)")
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
        oracle_guard.oracle_not_verified(
            "verify_genotype_gvcf_malformed_gatk_oracle.py", java, jar)
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

    results: list[dict] = []
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-genotype-malformed-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}"
                              for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# contract: on a malformed record a gated case must either fail "
          "loudly (non-zero exit, diagnostic, no data row) or still emit every "
          "well-formed record; exit 0 with a lost record is a violation.")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    gatk_exit={result['gatk_exit']} "
              f"(pinned {result['gatk_expected_exit']}) "
              f"native_exit={result['native_exit']} "
              f"native_output_exists={result['native_output_exists']}")
        print(f"    GATK   rows ({len(result['gatk_rows'])}):")
        for row in result["gatk_rows"]:
            print(f"        {row}")
        print(f"    NATIVE rows ({len(result['native_rows'])}):")
        for row in result["native_rows"]:
            print(f"        {row}")
        for line in [line for line in result["native_stderr_tail"].splitlines()
                     if line.startswith(("[E::", "[W::", "fastgatk-genotype-gvcf:"))]:
            print(f"    native stderr: {line}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "gatk_version": "4.6.2.0",
        "contract": "on a malformed VCF record native must fail loudly "
                    "(non-zero exit, diagnostic, no data row) or still emit "
                    "every well-formed record; exit 0 with a lost record is a "
                    "violation",
        "cases": results,
        "violations": violations,
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
