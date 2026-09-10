#!/usr/bin/env python3
"""Verify the Kokkos somatic evidence kernel against GATK's equations.

The real Mutect2 Java/native oracle is intentionally a separate integration
test: its PairHMM/assembly inputs can differ before SomaticLikelihoodsEngine
sees a matrix.  This test isolates that engine.  The reference below mirrors
the release-pinned GATK 4.6.2.0 implementation:

* initialize the Dirichlet posterior with ones;
* add one flat pseudocount per allele after each responsibility update;
* stop when relative L1 distance is below 0.001;
* compute log evidence with GATK's 1e-10 responsibility and 1e-8 entropy
  cutoffs.

The executable under test is a Kokkos API consumer and emits only the small
fixed matrices.  Keeping the reference here independent of the C++ kernel
prevents a C++ assertion and implementation from sharing the same mistake.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path


LOG10_TO_NATURAL = math.log(10.0)
NEGLIGIBLE_RESPONSIBILITY = 1.0e-10
ENTROPY_RESPONSIBILITY = 1.0e-8


def digamma_positive(value: float) -> float:
    """Commons Math Gamma.digamma branch used by the pinned Java release."""
    if not value > 0.0:
        raise ValueError("digamma domain")
    if value <= 1.0e-5:
        return -0.5772156649015329 - 1.0 / value
    x = value
    result = 0.0
    while x < 49.0:
        result -= 1.0 / x
        x += 1.0
    inverse = 1.0 / x
    inverse2 = inverse * inverse
    return (result + math.log(x) - 0.5 * inverse - inverse2 *
            (1.0 / 12.0 - inverse2 / 252.0))


def posterior(matrix: list[list[float]], prior: list[float] | None = None) -> tuple[list[float], int]:
    alleles = len(matrix)
    reads = len(matrix[0])
    if prior is None:
        prior = [1.0] * alleles
    alpha = [1.0] * alleles
    for _ in range(1000):
        total = sum(alpha)
        log_weights = [digamma_positive(x) - digamma_positive(total)
                       for x in alpha]
        counts = [0.0] * alleles
        for read in range(reads):
            values = [log_weights[allele] + matrix[allele][read] * LOG10_TO_NATURAL
                      for allele in range(alleles)]
            high = max(values)
            max_index = values.index(high)  # Java maxElementIndex: first tie.
            denominator = 1.0
            for allele, value in enumerate(values):
                if allele != max_index:
                    denominator += math.exp(value - high)
            log_sum = high + (math.log(denominator) if denominator != 1.0 else 0.0)
            for allele, value in enumerate(values):
                counts[allele] += math.exp(value - log_sum)
        new_alpha = [count + prior[index] for index, count in enumerate(counts)]
        distance = sum(abs(new_alpha[i] - alpha[i]) for i in range(alleles))
        alpha = new_alpha
        if distance / sum(alpha) < 0.001:
            return alpha, _ + 1
    raise AssertionError("GATK variational update did not converge")


def log_evidence(matrix: list[list[float]], prior: list[float] | None = None) -> tuple[float, list[float]]:
    if prior is None:
        prior = [1.0] * len(matrix)
    alpha, iterations = posterior(matrix, prior)
    log_weights = [digamma_positive(x) - digamma_positive(sum(alpha))
                   for x in alpha]
    contribution = 0.0
    for read in range(len(matrix[0])):
        values = [log_weights[allele] + matrix[allele][read] * LOG10_TO_NATURAL
                  for allele in range(len(matrix))]
        high = max(values)
        max_index = values.index(high)
        denominator = 1.0
        for allele, value in enumerate(values):
            if allele != max_index:
                denominator += math.exp(value - high)
        log_sum = high + (math.log(denominator) if denominator != 1.0 else 0.0)
        for allele, value in enumerate(values):
            responsibility = math.exp(value - log_sum)
            if responsibility >= NEGLIGIBLE_RESPONSIBILITY:
                contribution += responsibility * matrix[allele][read] * LOG10_TO_NATURAL
            if responsibility >= ENTROPY_RESPONSIBILITY:
                contribution -= responsibility * math.log(responsibility)

    prior_normalization = math.lgamma(sum(prior)) - sum(math.lgamma(x) for x in prior)
    posterior_normalization = math.lgamma(sum(alpha)) - sum(math.lgamma(x) for x in alpha)
    return ((prior_normalization - posterior_normalization + contribution) / LOG10_TO_NATURAL,
            alpha)


def expected_values() -> dict[str, float | list[float]]:
    biallelic = [[-1.0] * 4, [-0.01] * 4]
    biallelic_evidence, biallelic_alpha = log_evidence(biallelic)
    biallelic_reference = sum(biallelic[0])
    multiallelic = [
        [-1.0] * 4,
        [-0.01, -0.01, -3.0, -3.0],
        [-3.0, -3.0, -0.01, -0.01],
    ]
    all_evidence, _ = log_evidence(multiallelic)
    without_alt1, _ = log_evidence([multiallelic[0], multiallelic[2]])
    without_alt2, _ = log_evidence([multiallelic[0], multiallelic[1]])
    # The ALT row is zero-probability in the one-sided case.  It is represented
    # by -Infinity in Java and contributes no likelihood term; the posterior is
    # [2,1], hence AF=1/3 and evidence=-log10(2).  Compute this exact boundary
    # directly to avoid introducing a NaN into the reference implementation.
    one_sided_tlod = -math.log10(2.0)
    min_af_prior = [1.0, 1.0 - math.log(2.0) / math.log(0.1)]
    min_af_evidence, min_af_alpha = log_evidence(biallelic, min_af_prior)
    return {
        "biallelic_tlod": biallelic_evidence - biallelic_reference,
        "biallelic_af": biallelic_alpha[1] / sum(biallelic_alpha),
        "multiallelic_tlod": [all_evidence - without_alt1,
                               all_evidence - without_alt2],
        "one_sided_tlod": one_sided_tlod,
        "one_sided_af": 1.0 / 3.0,
        "min_af_tlod": min_af_evidence - biallelic_reference,
        "min_af_af": min_af_alpha[1] / sum(min_af_alpha),
    }


def assert_close(actual: float, expected: float, label: str, tolerance: float = 2.0e-10) -> None:
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(f"{label}: actual={actual!r} expected={expected!r}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SOMATIC_ORACLE_BINARY",
        root / "fastgatk-native/build/fastgatk-somatic-likelihood-oracle"))
    if not binary.is_file():
        raise SystemExit(f"missing somatic oracle executable: {binary}")
    run = subprocess.run([str(binary)], check=False, text=True,
                         capture_output=True)
    if run.returncode != 0:
        raise RuntimeError(f"somatic oracle failed: {run.stderr[-4000:]}")
    payload = json.loads(run.stdout.splitlines()[-1])
    if payload.get("status") != "pass" or not payload.get("execution_space"):
        raise AssertionError(f"unexpected somatic oracle payload: {payload}")
    expected = expected_values()
    assert_close(payload["biallelic_tlod"], expected["biallelic_tlod"], "biallelic TLOD")
    assert_close(payload["biallelic_af"], expected["biallelic_af"], "biallelic AF")
    for index, (actual, reference) in enumerate(zip(
            payload["multiallelic_tlod"], expected["multiallelic_tlod"])):
        assert_close(actual, reference, f"multiallelic TLOD[{index}]")
    assert_close(payload["one_sided_tlod"], expected["one_sided_tlod"],
                 "one-sided TLOD")
    assert_close(payload["one_sided_af"], expected["one_sided_af"],
                 "one-sided AF")
    assert_close(payload["min_af_tlod"], expected["min_af_tlod"],
                 "minimum-allele-fraction TLOD")
    assert_close(payload["min_af_af"], expected["min_af_af"],
                 "minimum-allele-fraction AF")
    print(json.dumps({
        "status": "pass",
        "execution_space": payload["execution_space"],
        "cases": ["biallelic", "multiallelic", "one-sided-zero-probability", "minimum-allele-fraction"],
        "max_abs_error": 2.0e-10,
        "reference": "GATK-4.6.2.0-SomaticLikelihoodsEngine-equations",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
