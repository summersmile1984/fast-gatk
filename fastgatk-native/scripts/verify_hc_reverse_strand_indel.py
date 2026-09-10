#!/usr/bin/env python3
"""Oracle for alignment-oriented reverse-strand CIGAR insertions."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_CALLING_INDEL_SMOKE_BINARY",
        root / "fastgatk-native/build/fastgatk-calling-indel-smoke",
    ))
    if not binary.exists():
        raise SystemExit(f"missing calling smoke binary: {binary}")

    # BAM SEQ is aligned to the reference coordinate system.  The reverse
    # flag changes strand annotations, but the CIGAR 2I slice remains AG and
    # must emit the anchored allele T>TAG at zero-based position 5.
    # Use a non-repeating flank so the branch is retained by the ordinary
    # ReadThreadingAssembler source-to-sink path.  An all-repeat reference
    # cannot promote a raw CIGAR I into an HC EventMap allele by itself.
    reference = "TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTA"
    aligned_read = "TACATTAGTGCTTC"
    cigar_insertion_begin = 6
    cigar_insertion_length = 2
    expected_inserted = aligned_read[cigar_insertion_begin:cigar_insertion_begin + cigar_insertion_length]
    expected_position = 5
    expected_alt = reference[expected_position] + expected_inserted
    assert expected_inserted == "AG"
    assert expected_alt == "TAG"

    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    completed = subprocess.run(
        [str(binary)], check=True, text=True, capture_output=True, env=env
    )
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    assert payload["status"] == "pass"
    assert payload["reverse_cigar_insertion_position"] == expected_position
    assert payload["reverse_cigar_insertion_alt"] == expected_alt

    print(json.dumps({
        "status": "pass",
        "position": expected_position,
        "aligned_inserted": expected_inserted,
        "alt": expected_alt,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
