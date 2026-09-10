#!/usr/bin/env python3
"""Offline Mutect2 TLOD-gap diagnostic (A4/A1 evidence chain, no CTest wiring).

Runs native Mutect2 on a pinned window with the FASTGATK_DEBUG_TLOD fragment
matrix and request-row dumps, then prints a compact report:

  * per-site TLOD reproduced by the Python transcription of the Java
    SomaticLikelihoodsEngine from the engine's own fragment matrix
    (formula-parity, see verify_mutect2_tlod_formula.py);
  * when a GATK --pair-hmm-results-file is supplied, the per-fragment
    allele-value divergence between native and GATK (fragment membership
    and read/haplotype row differences that the formula check cannot see).

This consolidates the R15-R22 discrimination chain (fragment-matrix dump,
row-level request dump, Java-formula replay) into one reusable command so a
future assembly/haplotype change can be re-audited without rediscovering the
pipeline.  Requires only the native binary and the chr17 fixture (plus an
optional GATK debug file).

Usage:
  diagnose_mutect2_tlod_gap.py [--gatk-pairhmm gatk.pairhmm.txt] [--window 17:69000-70000]
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_mutect2_tlod_formula as formula  # noqa: E402  (run/parse_dump/read_vcf/tlod_biallelic)

ROW_PATTERN = re.compile(
    r"\[FASTGATK_DEBUG_TLOD:req\] cand=(\d+) pos=(-?\d+) ref=(\S+) alt=(\S+)"
    r" read=(\d+) src=(\d+) start=(-?\d+) end=(-?\d+) group=(-?\d+)"
    r" fragonly=(\d+) hap=(\d+) canon=(\d+) state=(\d+) mask=(\d+)"
    r" allele=(\d+) rb=(\S+) rq=(\S+) iq=(\S+) dq=(\S+) gc=(\S+)"
    r" hb=(\S+) value=(-?[\d.eE+]+)")


def decode_quality_rows(text: str, site: int):
    """Parse the request-row dump into per-site native rows (hex->bytes)."""
    rows = []
    for line in text.splitlines():
        match = ROW_PATTERN.search(line)
        if not match:
            continue
        fields = match.groups()
        if int(fields[1]) != site:
            continue
        rows.append({
            "read": int(fields[4]), "group": int(fields[8]),
            "fragonly": int(fields[9]), "state": int(fields[12]),
            "rb": fields[15], "rq": bytes.fromhex(fields[16]),
            "hb": fields[20], "value": float(fields[21]),
        })
    return rows


def row_divergence(native_rows, gatk_path):
    """Per-group allele-level |delta| statistics vs GATK pairhmm rows."""
    gatk_index = collections.defaultdict(list)
    with open(gatk_path, encoding="utf-8") as handle:
        for line in handle:
            columns = line.split()
            if len(columns) != 7:
                continue
            rq = bytes(ord(char) - 33 for char in columns[2])
            gatk_index[(columns[1], rq, columns[0])].append(float(columns[6]))
    groups = collections.defaultdict(list)
    for row in native_rows:
        if row["fragonly"]:
            continue
        gatk_values = gatk_index.get((row["rb"], row["rq"], row["hb"]))
        gatk_value = gatk_values[0] if gatk_values else row["value"]
        groups[row["group"]].append(abs(row["value"] - gatk_value))
    divergent = {group: (len(deltas), round(max(deltas), 6))
                 for group, deltas in groups.items()
                 if max(deltas) > 1e-3}
    return divergent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gatk-pairhmm", type=Path, default=None,
                        help="optional GATK --pair-hmm-results-file")
    parser.add_argument("--window", default="17:69000-70000")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY",
        str(root / "fastgatk-native/build/fastgatk-mutect2")))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, bam, reference, Path(f"{bam}.bai"),
                Path(f"{reference}.fai"))
    if not all(path.is_file() for path in required):
        print("missing required fixture/binary:", [str(p) for p in required
                                                   if not p.is_file()])
        return 1

    with tempfile.TemporaryDirectory(prefix="fastgatk-tlod-diagnose-") as directory:
        work = Path(directory)
        vcf = work / "native.vcf.gz"
        env = dict(os.environ)
        env["FASTGATK_DEBUG_TLOD"] = "1"
        result = formula.run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", args.window,
            "-O", str(vcf), "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ], env)
        if result.returncode != 0:
            print("native run failed:", result.stderr[-2000:])
            return 1

        matrix_path = work / "matrix.txt"
        matrix_path.write_text(result.stderr, encoding="utf-8")
        candidates = formula.parse_dump(matrix_path)
        records = formula.read_vcf(vcf)
        alt_count = {}
        for (position, ref, _alt) in candidates:
            alt_count[(position, ref)] = alt_count.get((position, ref), 0) + 1

        print("=== TLOD formula parity (native fragment matrix vs Java "
              "transcription) ===")
        table = []
        for (position, ref, alt), (refs, alts) in sorted(candidates.items()):
            if alt_count.get((position, ref), 0) > 1:
                continue
            key = (position + 1, ref, alt)
            if key not in records:
                continue
            observed = records[key]
            expected = formula.tlod_biallelic(refs, alts)
            table.append((key[0], f"{ref}>{alt}", expected, observed,
                          abs(expected - observed)))
        for position, variant, expected, observed, delta in table:
            print(f"  {position} {variant}: formula={expected:.6f} "
                  f"native={observed:.6f} delta={delta:.2e}")
        if table:
            print(f"  max formula delta: "
                  f"{max(row[4] for row in table):.2e}")
        if args.gatk_pairhmm is not None:
            print("=== row-level divergence vs GATK (per fragment group, "
                  "|delta|>1e-3) ===")
            for pos in sorted({key[0] for key in candidates}):
                native_rows = decode_quality_rows(result.stderr, pos)
                if not native_rows:
                    continue
                divergent = row_divergence(native_rows, args.gatk_pairhmm)
                print(f"  site {pos + 1}: {len(native_rows)} rows, "
                      f"divergent groups {divergent or '{}'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
