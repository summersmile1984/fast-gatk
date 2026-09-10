#!/usr/bin/env python3
"""Compare Kokkos fragment-first aggregation with pinned GATK 4.6.2.0."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> dict:
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr[-4000:]}")
    return json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_FRAGMENT_AGGREGATION_ORACLE_BINARY",
        root / "fastgatk-native/build/fastgatk-fragment-aggregation-gatk-oracle"))
    java = root / "third_party/jdk17/bin/java"
    javac = root / "third_party/jdk17/bin/javac"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    source = root / "pairhmm-demo/tests/FragmentAggregationGatkOracle.java"
    for path in (binary, java, javac, jar, source):
        if not path.exists():
            raise SystemExit(f"missing fragment aggregation oracle asset: {path}")
    native = run([str(binary)])
    with tempfile.TemporaryDirectory(prefix="fastgatk-fragment-oracle-") as temporary:
        subprocess.run([str(javac), "-cp", str(jar), "-d", temporary, str(source)],
                       check=True, text=True, capture_output=True)
        oracle = run([str(java), "-cp", f"{temporary}:{jar}",
                      "FragmentAggregationGatkOracle"])
    if native.get("bits", [])[:6] != oracle.get("fragment_bits"):
        raise AssertionError(f"fragment sums differ: native={native.get('bits')} java={oracle}")
    if native.get("allele_bits") != oracle.get("allele_bits"):
        raise AssertionError(f"fragment-first marginalization differs: native={native} java={oracle}")
    # The first cell is the regression boundary: -1 + -Infinity must remain
    # negative infinity, not the old finite -1 value.
    if native["bits"][0] != "fff0000000000000":
        raise AssertionError("negative-infinity fragment likelihood was treated as missing")
    print(json.dumps({
        "status": "pass",
        "oracle": "GATK 4.6.2.0 AlleleLikelihoods.groupEvidence().marginalize()",
        "execution_space": native["execution_space"],
        "fragment_cells": 6,
        "bit_different": 0,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
