#!/usr/bin/env python3
"""Verify the pinned GATK 4.6.2.0 diploid normal-LOD equation.

The fixture executable emits fixed candidate-major PairHMM likelihood rows.
This verifier independently applies SomaticGenotypingEngine.diploidAltLogOdds:
sum(hom-ref) - sum(log10((10**ref + 10**alt) / 2)).  It intentionally checks
both finite ALT likelihoods and GATK's -Infinity zero-probability boundary so
the variational TLOD evidence cannot accidentally be reused as NLOD.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path


def mix_half(reference: float, alternate: float) -> float:
    if alternate == -math.inf:
        return reference + math.log10(0.5)
    high = max(reference, alternate)
    return high + math.log10(
        0.5 * 10.0 ** (reference - high) +
        0.5 * 10.0 ** (alternate - high))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SOMATIC_ORACLE_BINARY",
        root / "fastgatk-native/build/fastgatk-kernels/fastgatk-somatic-likelihood-oracle"))
    if not binary.is_file():
        raise SystemExit(f"missing somatic oracle executable: {binary}")
    run = subprocess.run([str(binary)], check=False, text=True, capture_output=True)
    if run.returncode != 0:
        raise RuntimeError(f"somatic normal-LOD oracle failed: {run.stderr[-4000:]}")
    payload = json.loads(run.stdout.splitlines()[-1])

    finite_expected = 4.0 * (-1.0 - mix_half(-1.0, -0.01))
    multiallelic_expected = [
        sum(-1.0 - mix_half(-1.0, alt) for alt in (-0.01, -0.01, -3.0, -3.0)),
        sum(-1.0 - mix_half(-1.0, alt) for alt in (-3.0, -3.0, -0.01, -0.01)),
    ]
    one_sided_expected = 0.0 - mix_half(0.0, -math.inf)
    for label, actual, expected in (
        ("finite normal LOD", payload["biallelic_nlod"], finite_expected),
        ("one-sided normal LOD", payload["one_sided_nlod"], one_sided_expected),
    ):
        if not math.isfinite(actual) or abs(actual - expected) > 1.0e-12:
            raise AssertionError(f"{label}: actual={actual!r} expected={expected!r}")
    for index, (actual, expected) in enumerate(zip(
            payload["multiallelic_nlod"], multiallelic_expected)):
        if not math.isfinite(actual) or abs(actual - expected) > 1.0e-12:
            raise AssertionError(
                f"multiallelic normal LOD[{index}]: actual={actual!r} expected={expected!r}")
    if abs(payload["biallelic_nlod"] + payload["biallelic_tlod"]) < 1.0e-6:
        raise AssertionError("NLOD unexpectedly aliases variational TLOD evidence")
    print(json.dumps({
        "status": "pass",
        "execution_space": payload["execution_space"],
        "finite_nlod": payload["biallelic_nlod"],
        "one_sided_nlod": payload["one_sided_nlod"],
        "reference": "GATK-4.6.2.0-SomaticGenotypingEngine-diploidAltLogOdds",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
