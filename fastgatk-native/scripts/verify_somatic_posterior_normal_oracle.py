#!/usr/bin/env python3
"""Pin the Mutect2 tumor/normal evidence-count boundary.

GATK 4.6.2.0's ``SomaticGenotypingEngine`` receives independent tumor and
normal ``AlleleLikelihoods`` matrices.  In particular, the normal's retained
fragment count is not required to equal the tumor's count.  This oracle uses
the Kokkos API fixture emitted by ``fastgatk-somatic-likelihood-oracle`` and
checks that a two-fragment normal contributes exactly once to the posterior
penalty beside four tumor fragments.  The expected values are evaluated from
the release's log10 mixture equations, independently of the C++ kernel.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path


def mix_log10(reference: float, alternate: float, fraction: float) -> float:
    high = max(reference, alternate)
    return high + math.log10(
        (1.0 - fraction) * 10.0 ** (reference - high) +
        fraction * 10.0 ** (alternate - high))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SOMATIC_ORACLE_BINARY",
        root / "fastgatk-native/build/fastgatk-kernels/fastgatk-somatic-likelihood-oracle"))
    if not binary.is_file():
        raise SystemExit(f"missing somatic oracle executable: {binary}")
    run = subprocess.run([str(binary)], check=False, text=True, capture_output=True)
    if run.returncode != 0:
        raise RuntimeError(f"somatic posterior oracle failed: {run.stderr[-4000:]}")
    payload = json.loads(run.stdout.splitlines()[-1])

    # Four tumor rows have best evidence -0.04.  The two normal rows have
    # reference evidence -2.0, hom-alt evidence -0.02, and heterozygous
    # evidence 2*mix(-1,-0.01,0.5).  GATK's normal penalty takes the larger
    # alternate model against the normal reference model.
    tumor_best = 4.0 * -0.01
    normal_reference = 2.0 * -1.0
    normal_germline = 2.0 * mix_log10(-1.0, -0.01, 0.5)
    normal_hom_alt = 2.0 * mix_log10(-1.0, -0.01, 1.0)
    normal_penalty = max(normal_hom_alt, normal_germline) - normal_reference
    expected_with_normal = tumor_best - max(0.0, normal_penalty)
    expected_without_normal = tumor_best
    for label, actual, expected in (
        ("unequal normal somatic evidence", payload["unequal_normal_somatic_evidence"], expected_with_normal),
        ("no-normal somatic evidence", payload["no_normal_somatic_evidence"], expected_without_normal),
    ):
        if not math.isfinite(actual) or abs(actual - expected) > 1.0e-12:
            raise AssertionError(f"{label}: actual={actual!r} expected={expected!r}")
    if payload["unequal_normal_informative"] != 4:
        raise AssertionError(payload)
    if not payload.get("execution_space"):
        raise AssertionError(payload)
    print(json.dumps({
        "status": "pass",
        "execution_space": payload["execution_space"],
        "tumor_evidence_units": 4,
        "normal_evidence_units": 2,
        "normal_count_independent": True,
        "somatic_evidence_with_normal": payload["unequal_normal_somatic_evidence"],
        "somatic_evidence_without_normal": payload["no_normal_somatic_evidence"],
        "reference": "GATK-4.6.2.0-SomaticGenotypingEngine-independent-normal-evidence",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
