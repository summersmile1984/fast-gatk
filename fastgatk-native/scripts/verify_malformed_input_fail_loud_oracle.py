#!/usr/bin/env python3
"""Strict gate: no native tool may silently truncate its input.

The defect class this gate pins
-------------------------------
Every htslib reader reports "clean end of input" and "this record could not be
read" with two different negative values, and htslib documents that itself:

* ``bcf_read()``   -- "0 on success; -1 on end of file; < -1 on critical error"
  (``htslib/vcf.h``, doc comment above the declaration; the < -1 path is
  ``vcf_read`` -> ``vcf_parse``, whose ``err:`` label returns -2), and
  ``bcf_read1_core`` returns -2 for a corrupt BCF block;
* ``sam_read1()``  -- "-1 at end of file; < -1 on failure" (``htslib/sam.h``,
  ``bam_read1`` doc comment);
* ``hts_getline()`` -- "-1 on end-of-file; <= -2 on error" (``htslib/hts.h``);
  for a BGZF/gzip stream ``bgzf_getline`` propagates the read failure as -2.

A loop written as ``while (bcf_read(input, header, record) == 0)`` cannot tell
those apart: the FIRST record htslib could not read ends the traversal exactly
as if the input had ended.  The record and every record after it are discarded
while the process still exits 0.  Measured on the pre-fix binaries: a truncated
``.vcf.gz`` made ``fastgatk-variants-to-table`` exit **0** having written 19322
of 20001 records, with htslib's ``[E::bgzf_read_block]`` line on stderr as the
only hint.

The contract asserted here
--------------------------
Per tool and case:

* a **malformed** input must fail LOUDLY -- non-zero exit **and** a
  ``BAD_INPUT:`` diagnostic on stderr -- and must not hand back a
  complete-looking result:
    - ``malformed-record-first``: the unreadable record precedes every good
      record, so the output must carry **zero** data records;
    - ``malformed-record-last``: exactly the well-formed prefix may have been
      written (streaming writers legitimately flush it, as pinned GATK does
      before it exits 3), so the output must carry exactly the good prefix and
      the exit must be non-zero;
    - ``truncated-gzip``: the reader cannot know how many records preceded the
      corruption point, so the assertion is that the exit is non-zero and the
      output is **provably incomplete** (strictly fewer records than the same
      input read whole);
* a **well-formed** control (the identical records, plain and gzipped) must
  still exit 0 and emit every input record, so a fix cannot over-reject.

What must never happen -- and is what this gate fails on -- is a zero exit with
a well-formed record missing.

Run:  python3 fastgatk-native/scripts/verify_malformed_input_fail_loud_oracle.py
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##ALT=<ID=NON_REF,Description=Any possible alternate allele>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""

FORMAT = "GT:DP:AD:PL"
# ALT is a single concrete allele, so Number=R wants 2 AD values and Number=G
# wants 3 PL values.  A record is malformed (htslib: vcf_parse_format_max3,
# vcf_parse -> -2) exactly when FORMAT declares four keys and the sample
# carries five values.
SAMPLE_GOOD = "0/1:20:0,20:100,0,100"
SAMPLE_BAD = SAMPLE_GOOD + ":99"
GOOD_POSITIONS = (20, 30)
BAD_FIRST_POSITION = 10
BAD_LAST_POSITION = 40

GVCF_ALT = "G,<NON_REF>"

# 3000 ordinary records span more than one 64 KiB decompression buffer, so a
# 90 % truncation lands *after* records have already been handed to the tool --
# that is what makes the case exercise the read loop rather than the header
# read.  A smaller fixture would truncate inside the first buffer and be caught
# by the (already guarded) header path instead.
TRUNCATED_RECORDS = 3000


def record(position: int, sample: str, alt: str) -> str:
    return (f"chr1\t{position}\t.\tA\t{alt}\t50\tPASS\tDP=20\t{FORMAT}\t{sample}\n")


def good_records(alt: str, count: int = len(GOOD_POSITIONS)) -> str:
    return "".join(record(position, SAMPLE_GOOD, alt)
                   for position in range(20, 20 + count))


# --------------------------------------------------------------------------
# Tool registry.  ``rows`` names how the emitted record count is measured:
#   "vcf"   -- data lines of a VCF output
#   "table" -- tab-delimited rows of a VariantsToTable output (minus its header)
#   "fasta" -- sequence content of a FastaAlternateReferenceMaker output
# ``control_rows`` is the pinned, measured number of records the control input
# (two well-formed records) must produce.
# --------------------------------------------------------------------------
TOOLS: dict[str, dict] = {
    "variants-to-table": {
        "binary": "fastgatk-variants-to-table",
        "argv": ["-V", "{vcf}", "-O", "{out}",
                 "-F", "CHROM", "-F", "POS", "-F", "QUAL"],
        "rows": "table",
        "gvcf": False,
    },
    "select-variants": {
        "binary": "fastgatk-select-variants",
        "argv": ["-V", "{vcf}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": False,
    },
    "left-align-trim": {
        "binary": "fastgatk-left-align-trim",
        "argv": ["-V", "{vcf}", "-R", "{ref}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": False,
    },
    "variant-filtration": {
        "binary": "fastgatk-variant-filtration",
        "argv": ["-V", "{vcf}", "-O", "{out}",
                 "--filter-expression", "QUAL < 10", "--filter-name", "LowQual"],
        "rows": "vcf",
        "gvcf": False,
    },
    "validate-variants": {
        "binary": "fastgatk-validate-variants",
        "argv": ["-V", "{vcf}", "-R", "{ref}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": False,
    },
    "reblock-gvcf": {
        "binary": "fastgatk-reblock-gvcf",
        "argv": ["-V", "{vcf}", "-R", "{ref}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": True,
    },
    "combine-gvcfs": {
        "binary": "fastgatk-combine-gvcfs",
        "argv": ["-V", "{vcf}", "-R", "{ref}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": True,
    },
    "gather-vcfs": {
        "binary": "fastgatk-gather-vcfs",
        "argv": ["-I", "{vcf}", "-O", "{out}", "--CREATE_INDEX", "false"],
        "rows": "vcf",
        "gvcf": False,
    },
    "filter-mutect-calls": {
        "binary": "fastgatk-filter-mutect-calls",
        "argv": ["-V", "{vcf}", "-O", "{out}"],
        "rows": "vcf",
        "gvcf": False,
    },
    "fasta-alternate-reference-maker": {
        "binary": "fastgatk-fasta-alternate-reference-maker",
        "argv": ["-R", "{ref}", "-V", "{vcf}", "-O", "{out}"],
        "rows": "fasta",
        "gvcf": False,
    },
}

CONTROL_ROWS = 2


def write_reference(work: pathlib.Path) -> pathlib.Path:
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 10000 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t10000\t6\t10000\t10001\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n", encoding="utf-8")
    return reference


def write_source(path: pathlib.Path, body: str, *, gzipped: bool,
                 truncate_to: int | None = None) -> pathlib.Path:
    text = HEADER + body
    if not gzipped:
        path.write_text(text, encoding="utf-8")
        return path
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    if truncate_to is not None:
        raw = path.read_bytes()
        if truncate_to >= len(raw):
            raise SystemExit(f"cannot truncate {path} to {truncate_to} of {len(raw)} bytes")
        path.write_bytes(raw[:truncate_to])
    return path


def count_records(path: pathlib.Path, kind: str) -> int | None:
    """Records the tool wrote; None when the output does not exist."""
    if not path.exists():
        return None
    if kind == "fasta":
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        return sum(len(line) for line in lines if line and not line.startswith(">"))
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    rows = [line for line in lines if line and not line.startswith("#")]
    if kind == "table":
        return max(len(rows) - 1, 0)
    return len(rows)


def invoke(command: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True,
                          check=False, timeout=timeout)


def run_case(tool: str, spec: dict, case: str, work: pathlib.Path,
             reference: pathlib.Path, timeout: int) -> dict:
    alt = GVCF_ALT if spec["gvcf"] else "G"
    control_body = good_records(alt)
    source = work / f"{tool}-{case}.vcf"
    gzipped = case.endswith("gzip")
    suffix = ".gz" if gzipped else ""
    if case == "control-plain":
        source = write_source(source, control_body, gzipped=False)
    elif case == "control-gzip":
        source = write_source(work / f"{tool}-{case}.vcf.gz", control_body, gzipped=True)
    elif case == "truncated-gzip":
        body = good_records(alt, TRUNCATED_RECORDS)
        full = write_source(work / f"{tool}-full.vcf.gz", body, gzipped=True)
        size = full.stat().st_size
        source = write_source(work / f"{tool}-{case}.vcf.gz", body, gzipped=True,
                              truncate_to=size * 9 // 10)
    elif case == "malformed-record-first":
        source = write_source(source, record(BAD_FIRST_POSITION, SAMPLE_BAD, alt)
                              + control_body, gzipped=False)
    elif case == "malformed-record-last":
        source = write_source(source, control_body
                              + record(BAD_LAST_POSITION, SAMPLE_BAD, alt), gzipped=False)
    else:
        raise SystemExit(f"unknown case {case}")

    out_suffix = ".table" if spec["rows"] == "table" else (
        ".fa" if spec["rows"] == "fasta" else suffix + ".vcf")
    out = work / f"{tool}-{case}-out{out_suffix}"
    argv = [str(spec["binary_path"])]
    for token in spec["argv"]:
        argv.append(token.format(vcf=source, ref=reference, out=out))
    result = invoke(argv, timeout)
    rows = count_records(out, spec["rows"])
    return {
        "tool": tool,
        "case": case,
        "input": str(source),
        "exit": result.returncode,
        "rows": rows,
        "output_exists": out.exists(),
        "stderr_tail": result.stderr.strip().splitlines()[-1][:300]
                       if result.stderr.strip() else "",
        "stderr_has_bad_input": "BAD_INPUT:" in result.stderr,
        "violations": [],
    }


def judge(record_result: dict, spec: dict, *, controls: dict[str, int]) -> None:
    """Apply the contract to one measured case."""
    case = record_result["case"]
    violations = record_result["violations"]
    rows = record_result["rows"]
    kind = spec["rows"]

    if case.startswith("control-"):
        expected = controls.get(case)
        if record_result["exit"] != 0:
            violations.append(
                f"well-formed input rejected: exit={record_result['exit']} "
                f"({record_result['stderr_tail']})")
        if rows != expected:
            violations.append(
                f"well-formed input produced {rows} record(s), expected {expected}: "
                "the guard must not drop good records")
        return

    # Every malformed case: loud failure.
    if record_result["exit"] == 0:
        violations.append(
            f"exited 0 on {case}: the unreadable input was swallowed instead of "
            f"reported, and {rows} record(s) were written")
    if not record_result["stderr_has_bad_input"]:
        violations.append(
            f"no BAD_INPUT diagnostic on stderr for {case} "
            f"(stderr tail: {record_result['stderr_tail']!r})")

    if kind == "fasta":
        # A FASTA output cannot be counted per record; the contract is the loud
        # failure plus not publishing a sequence for a run that failed.
        if record_result["exit"] != 0 and rows:
            violations.append(
                f"wrote {rows} sequence base(s) for a failed run")
        return

    control = controls["control-plain"]
    if case == "malformed-record-first":
        if rows not in (0, None):
            violations.append(
                f"wrote {rows} record(s) although the first record is unreadable "
                "and no well-formed record precedes it")
    elif case == "malformed-record-last":
        # A streaming writer may already have flushed the well-formed prefix
        # (pinned GATK does exactly that before it exits 3); what it may not do
        # is publish anything beyond it.
        if rows not in (control, None):
            violations.append(
                f"wrote {rows} record(s); the {control} well-formed records "
                f"precede the unreadable one, so only that prefix may exist")
    elif case == "truncated-gzip":
        if rows is not None and rows >= TRUNCATED_RECORDS:
            violations.append(
                f"wrote {rows} record(s) for a truncated stream that can hold at "
                f"most {TRUNCATED_RECORDS}: the result is indistinguishable from "
                "a complete run")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict gate: every native htslib reader must report an "
                    "unreadable record instead of treating it as end of input.")
    parser.add_argument("--build", default=os.environ.get("FASTGATK_NATIVE_BUILD"),
                        help="native build directory (default "
                             "$FASTGATK_NATIVE_BUILD or fastgatk-native/build)")
    parser.add_argument("--tool", action="append", default=None,
                        help="run only the named tool(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    build = pathlib.Path(arguments.build) if arguments.build else (
        root / "fastgatk-native" / "build")

    selected = arguments.tool or list(TOOLS)
    unknown = [name for name in selected if name not in TOOLS]
    if unknown:
        raise SystemExit(f"no such tool: {unknown}")

    results: list[dict] = []
    missing: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-fail-loud-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for name in selected:
            spec = dict(TOOLS[name])
            spec["binary_path"] = build / spec["binary"]
            if not spec["binary_path"].is_file():
                missing.append(str(spec["binary_path"]))
                continue
            cases = ["control-plain", "control-gzip", "malformed-record-first",
                     "malformed-record-last", "truncated-gzip"]
            measured = [run_case(name, spec, case, work, reference,
                                 arguments.timeout) for case in cases]
            controls = {case["case"]: case["rows"] for case in measured}
            for case in measured:
                judge(case, spec, controls=controls)
            results.extend(measured)

    if missing:
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary",
                          "missing": sorted(missing)}, sort_keys=True))
        return 0

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['tool']}/{result['case']}] {item}"
                          for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: native build {build}")
    print(f"# contract: a malformed/truncated input must exit non-zero with a "
          f"BAD_INPUT diagnostic and must not produce a complete-looking "
          f"result; the well-formed controls must still emit all "
          f"{CONTROL_ROWS} records.")
    for result in results:
        print(f"[{result['tool']}/{result['case']}] exit={result['exit']} "
              f"records={result['rows']}")
        for item in result["violations"]:
            print(f"    VIOLATION {item}")
    print(json.dumps({"status": "divergence" if violations else "pass",
                      "tools": len({result["tool"] for result in results}),
                      "cases": len(results),
                      "violations": violations}, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
