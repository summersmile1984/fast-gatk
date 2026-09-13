#!/usr/bin/env python3
"""Reproduce the dense full-region GATK/native residual map for GenotypeGVCFs.

Every round of the dense work in this session re-derived the same comparison by
hand (run GATK, run native, index both outputs by POS, diff the row lists and
classify by ALT shape).  This script makes that one command, so the numbers in
``fastgatk-native/evidence/2026-09-13-round60-dense/full-region-residual-map.md``
can be re-verified after any change.

What it does
------------
Runs pinned GATK 4.6.2.0 and the native binary on the SAME input, the SAME
reference and the SAME ``-L`` window (dense mode by default), then reports:

* row and position counts for both tools;
* positions only GATK emits, only native emits, and positions that differ;
* the differing positions grouped by (GATK ALT shape, native ALT shape, the
  leading columns that differ) -- the same classification the evidence file uses;
* for the ``--positions`` given on the command line, a per-column dump so a
  single locus can be inspected without re-running anything.

Usage
-----
    # default corpus window (GATK's own chr20 HaplotypeCaller gVCF, 100 kb)
    python3 fastgatk-native/scripts/measure_dense_residual.py

    # a narrower window, with per-column detail for two loci
    python3 fastgatk-native/scripts/measure_dense_residual.py \
        -L 20:10008900-10009000 --positions 10008964

    # compare two already-produced outputs instead of running the tools
    python3 fastgatk-native/scripts/measure_dense_residual.py \
        --gatk-vcf gatk.vcf --native-vcf native.vcf

Exit status is 0 whenever the comparison ran (this is a measurement tool, not a
gate); pass ``--require-identical`` to make it exit 1 when any position differs.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from oracle_guard import oracle_not_verified  # noqa: E402

DEFAULT_INPUT = ("gatk-source/src/test/resources/org/broadinstitute/hellbender/"
                 "tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf")
DEFAULT_REFERENCE = "testdata/chr20/reference/GRCh37.chr20.fa"
DEFAULT_WINDOW = "20:10000000-10099999"


def read_rows(path: pathlib.Path) -> dict[int, list[str]]:
    rows: dict[int, list[str]] = collections.defaultdict(list)
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows[int(fields[1])].append(line)
    return rows


def shape(row: str) -> str:
    alt = row.split("\t")[4]
    if alt == "*":
        return "star"
    if alt == ".":
        return "refonly"
    return "multiallelic" if "," in alt else "concrete"


def classify(diff: list[int], gatk: dict[int, list[str]],
             native: dict[int, list[str]]) -> list[tuple[tuple, int]]:
    groups: collections.Counter[tuple] = collections.Counter()
    for position in diff:
        gatk_row = gatk[position][0].split("\t")
        native_row = native[position][0].split("\t")
        columns = tuple(index for index in range(min(len(gatk_row), len(native_row)))
                        if gatk_row[index] != native_row[index])
        groups[(shape(gatk[position][0]), shape(native[position][0]), columns)] += 1
    return groups.most_common()


def run(command: list[str], label: str, timeout: int) -> int:
    result = subprocess.run(command, text=True, capture_output=True, check=False,
                            timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-2000:])
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-L", "--intervals", default=DEFAULT_WINDOW,
                        help=f"interval to compare (default {DEFAULT_WINDOW})")
    parser.add_argument("-V", "--input", default=DEFAULT_INPUT, help="input gVCF")
    parser.add_argument("-R", "--reference", default=DEFAULT_REFERENCE)
    parser.add_argument("--dense", action="store_true", default=True,
                        help="pass --include-non-variant-sites (default)")
    parser.add_argument("--no-dense", dest="dense", action="store_false")
    parser.add_argument("--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"))
    parser.add_argument("--gatk-vcf", default=None,
                        help="use this existing GATK output instead of running GATK")
    parser.add_argument("--native-vcf", default=None,
                        help="use this existing native output instead of running native")
    parser.add_argument("--positions", type=int, action="append", default=[],
                        help="dump every column for these positions")
    parser.add_argument("--require-identical", action="store_true",
                        help="exit 1 if any position differs (not a gate by default)")
    parser.add_argument("--timeout", type=int, default=2400)
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-genotype-gvcf")

    work = pathlib.Path(tempfile.mkdtemp(prefix="fastgatk-dense-residual-"))
    gatk_vcf = pathlib.Path(arguments.gatk_vcf) if arguments.gatk_vcf else work / "gatk.vcf"
    native_vcf = (pathlib.Path(arguments.native_vcf) if arguments.native_vcf
                  else work / "native.vcf")

    if arguments.gatk_vcf is None:
        if not (java.is_file() and jar.is_file()):
            oracle_not_verified("measure_dense_residual.py", java, jar)
            if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
                raise SystemExit("missing pinned GATK assets")
            print(json.dumps({"status": "skipped", "reason": "no pinned GATK"},
                             sort_keys=True))
            return 0
    if arguments.native_vcf is None and not native.is_file():
        raise SystemExit(f"missing native binary: {native}")

    dense = ["--include-non-variant-sites"] if arguments.dense else []
    if arguments.gatk_vcf is None:
        run([str(java), "-Xmx4g", "-jar", str(jar), "GenotypeGVCFs",
             "-R", str(root / arguments.reference), "-V", str(root / arguments.input),
             "-L", arguments.intervals, *dense, "-O", str(gatk_vcf),
             "--create-output-variant-index", "false"], "GATK GenotypeGVCFs", arguments.timeout)
    if arguments.native_vcf is None:
        run([str(native), "-R", str(root / arguments.reference),
             "-V", str(root / arguments.input), "-L", arguments.intervals, *dense,
             "--gatk-compatible-annotations", "-O", str(native_vcf)],
            "native GenotypeGVCFs", arguments.timeout)

    gatk = read_rows(gatk_vcf)
    native_rows = read_rows(native_vcf)
    gatk_positions = set(gatk)
    native_positions = set(native_rows)
    only_gatk = sorted(gatk_positions - native_positions)
    only_native = sorted(native_positions - gatk_positions)
    diff = sorted(position for position in gatk_positions & native_positions
                  if gatk[position] != native_rows[position])

    print(f"# intervals={arguments.intervals} dense={arguments.dense} "
          f"input={arguments.input}")
    print(f"# GATK rows={sum(len(v) for v in gatk.values())} positions={len(gatk)}")
    print(f"# native rows={sum(len(v) for v in native_rows.values())} "
          f"positions={len(native_rows)}")
    print(f"# only GATK={len(only_gatk)} only native={len(only_native)} differ={len(diff)}")
    if only_gatk:
        print(f"#   only-GATK positions: {only_gatk[:12]}")
    if only_native:
        print(f"#   only-native positions: {only_native[:12]}")
    print(f"{'GATK shape':12} {'native shape':12} {'differing columns':24} count")
    for (gatk_shape, native_shape, columns), count in classify(diff, gatk, native_rows):
        print(f"{gatk_shape:12} {native_shape:12} {str(columns):24} {count}")

    for position in arguments.positions:
        print(f"--- position {position}")
        for label, rows in (("GATK", gatk), ("native", native_rows)):
            if position not in rows:
                print(f"  {label}: <absent>")
                continue
            for row in rows[position]:
                fields = row.split("\t")
                print(f"  {label}: POS={fields[1]} REF={fields[3]} ALT={fields[4]} "
                      f"QUAL={fields[5]} FILTER={fields[6]}")
                print(f"         INFO={fields[7]}")
                print(f"         FORMAT={fields[8]} SAMPLE={fields[9] if len(fields) > 9 else ''}")

    payload = {
        "intervals": arguments.intervals,
        "dense": arguments.dense,
        "gatk_positions": len(gatk),
        "native_positions": len(native_rows),
        "only_gatk": len(only_gatk),
        "only_native": len(only_native),
        "differ": len(diff),
        "classes": [
            {"gatk_shape": gatk_shape, "native_shape": native_shape,
             "columns": list(columns), "count": count}
            for (gatk_shape, native_shape, columns), count in classify(diff, gatk, native_rows)
        ],
    }
    print(json.dumps(payload, sort_keys=True))
    if arguments.require_identical and (diff or only_gatk or only_native):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
