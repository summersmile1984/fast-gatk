#!/usr/bin/env python3
"""Compare raw IEEE-754 likelihoods emitted by native and the GATK oracle."""
import math
import struct
import sys

def read(path):
    return [int(line.strip(), 16) for line in open(path, encoding="ascii") if line.strip()]

def as_double(bits):
    return struct.unpack("=d", struct.pack("=Q", bits))[0]

def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_pairhmm.py native.hex gatk.hex")
    a, b = read(sys.argv[1]), read(sys.argv[2])
    if len(a) != len(b): raise SystemExit(f"length mismatch: {len(a)} != {len(b)}")
    diffs = [abs(as_double(x) - as_double(y)) for x, y in zip(a, b)]
    bit_different = sum(x != y for x, y in zip(a, b))
    print({"pairs": len(a), "bit_different": bit_different,
           "max_abs": max(diffs, default=0.0),
           "all_within_1e-12": all(x <= 1e-12 for x in diffs)})

if __name__ == "__main__": main()
