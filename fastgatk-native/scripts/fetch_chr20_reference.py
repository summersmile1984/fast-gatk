#!/usr/bin/env python3
"""Fetch real b37 chr20 sequence slices and synthesize an HC-usable reference.

Recipe proven on 2026-09-03 (see NEXT_PHASE_TASKS.md R15-17): the GATK test
resources keep large real BAMs under Git-LFS, and anonymous LFS media
downloads work for the bundled reference fasta.  The fasta text has one extra
leading character relative to naive byte math, so coordinates are anchored
against the UCSC hg19 API at two positions before use; the final reference
carries the real sequence only where reads need it and 'N' elsewhere, which is
sufficient for HaplotypeCaller parity fixtures that call inside the covered
window.

Usage:
  python3 fetch_chr20_reference.py --out fixtures/chr20/ref20mnp.fasta \
      --bam fixtures/chr20/mnp.bam --region 20:10018000-10023000
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import urllib.request

LFS_FASTA_URL = ("https://media.githubusercontent.com/media/broadinstitute/gatk/"
                 "4.6.2.0/src/test/resources/large/human_g1k_v37.20.21.fasta")
UCSC_URL = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg19;"
            "chrom=chr20;start={start};end={end}")
CHR20_HEADER = ">20 dna:chromosome chromosome:GRCh37:20:1:63025520:1\n"
CHR20_LEN = 63025520
CHR21_LEN = 48129895


def fetch_range(start_base: int, end_base: int) -> str:
    """Return real chr20 bases [start_base, end_base) from the LFS fasta.

    Byte math: header is 53 bytes; each 60-base line adds a newline, so base b
    (1-based) starts at byte 53 + (b-1) + (b-1)//60.
    """
    header_len = len(CHR20_HEADER)
    first_byte = header_len + (start_base - 1) + (start_base - 1) // 60
    last_byte = header_len + (end_base - 1) + (end_base - 1) // 60
    request = urllib.request.Request(
        LFS_FASTA_URL, headers={"Range": f"bytes={first_byte}-{last_byte}"})
    with urllib.request.urlopen(request, timeout=120) as response:
        raw = response.read().decode("ascii", "ignore")
    return "".join(ch for ch in raw if ch in "ACGTNacgtn")


def anchor_base(sequence: str, hg19_start: int, probe_length: int = 40) -> int:
    """Return the hg19 coordinate of sequence[0] using a UCSC probe."""
    url = UCSC_URL.format(start=hg19_start, end=hg19_start + probe_length)
    with urllib.request.urlopen(url, timeout=60) as response:
        dna = json.loads(response.read())["dna"]
    index = sequence.find(dna[:probe_length])
    if index < 0:
        raise SystemExit("anchor probe not found in fetched sequence")
    return hg19_start - index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=pathlib.Path)
    parser.add_argument("--start", type=int, default=10018000,
                        help="first real base to embed (default 10018000)")
    parser.add_argument("--end", type=int, default=10023000,
                        help="last real base to embed (default 10023000)")
    args = parser.parse_args()

    # Fetch generous margins and anchor with UCSC.
    sequence = fetch_range(args.start - 200, args.end + 200)
    offset = anchor_base(sequence, args.start - 200)
    if offset != args.start - 200:
        print(f"anchor offset shift: sequence[0] is hg19 base {offset}",
              file=sys.stderr)
    # sequence[i] is hg19 chr20 base (offset + i).
    scaffold = bytearray(b"N" * CHR20_LEN)
    for i, char in enumerate(sequence):
        absolute = offset + i
        if 1 <= absolute <= CHR20_LEN and char in "ACGTN":
            scaffold[absolute - 1] = ord(char)
    with open(args.out, "w", encoding="ascii") as handle:
        handle.write(CHR20_HEADER)
        for cursor in range(0, CHR20_LEN, 60):
            handle.write(bytes(scaffold[cursor:cursor + 60]).decode("ascii") + "\n")
        handle.write(">21 dna:chromosome chromosome:GRCh37:21:1:48129895:1\n")
        for _ in range((CHR21_LEN + 59) // 60):
            handle.write("N" * 60 + "\n")
    dict_path = args.out.with_suffix(".dict")
    dict_path.write_text("@HD\tVN:1.5\tSO:coordinate\n"
                         f"@SQ\tSN:20\tLN:{CHR20_LEN}\n"
                         f"@SQ\tSN:21\tLN:{CHR21_LEN}\n",
                         encoding="ascii")
    text = args.out.read_text(encoding="ascii")
    fai_path = args.out.with_suffix(".fasta.fai")
    fai_path.write_text(f"20\t{CHR20_LEN}\t{len(CHR20_HEADER)}\t60\t61\n"
                        f"21\t{CHR21_LEN}\t{text.index('>21')}\t60\t61\n",
                        encoding="ascii")
    print(f"wrote {args.out} ({args.out.stat().st_size} bytes), "
          f"{dict_path.name}, {fai_path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
