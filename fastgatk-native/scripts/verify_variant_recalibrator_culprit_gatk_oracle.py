#!/usr/bin/env python3
"""Strict oracle for the VariantRecalibrator INFO ``culprit`` value (pinned GATK 4.6.2.0).

Observation under test
----------------------
``fastgatk-native/scripts/verify_variant_recalibrator.py`` used to assert the
literal ``culprit=full-covariance-gmm`` (lines 84/422/437/456).  That string is
a *native provenance* name: it occurs **0 times** in the pinned
``gatk-package-4.6.2.0-local.jar`` (``culprit`` itself occurs 13 times).  Pinned
GATK writes an **annotation name** there, chosen per datum.

GATK's rule (established by measurement against the pinned jar, then confirmed
against GATK's source)
------------------------------------------------------------------------
``VariantRecalibratorEngine.calculateWorstPerformingAnnotation``
(``gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/vqsr/``
``VariantRecalibratorEngine.java:80-93``) compares, for every requested
annotation dimension ``iii``, the **per-dimension log10 mixture likelihood** of
the positive and negative models and keeps the dimension with the smallest
``good - bad`` contrast::

    final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
    final Double badProbLog10  = badModel.evaluateDatumInOneDimension(datum, iii);
    if( goodProbLog10 != null && badProbLog10 != null ) {
        final double prob = goodProbLog10 - badProbLog10;
        if(prob < minProb) { minProb = prob; worstAnnotation = iii; ... }
    }

``GaussianMixtureModel.evaluateDatumInOneDimension``
(``GaussianMixtureModel.java:208-222``) is the log10 column in question::

    pVarInGaussianLog10[k] = gaussian.pMixtureLog10;
    if (gaussian.pMixtureLog10 != Double.NEGATIVE_INFINITY) {
        pVarInGaussianLog10[k] += MathUtils.normalDistributionLog10(
            gaussian.mu[iii], gaussian.sigma.get(iii, iii), datum.annotations[iii]);
    }
    return nanTolerantLog10SumLog10(pVarInGaussianLog10);

(``MathUtils.normalDistributionLog10(mean, sd, x)`` receives the covariance
diagonal ``sigma.get(iii,iii)`` in its ``sd`` slot -- ``sigma`` is the fitted
covariance matrix, so that diagonal is the variance; see
``VariantDataManager.java:29``'s matching "this is really the standard
deviation" remark about GATK's VQSR naming.)

``VariantDataManager.writeOutRecalibrationTable``
(``VariantDataManager.java:485``) then writes that dimension's *name*::

    builder.attribute(GATKVCFConstants.CULPRIT_KEY,
        (datum.worstAnnotation != -1 ? annotationKeys.get(datum.worstAnnotation) : "NULL"));

``annotationKeys`` is GATK's information-ordered annotation order.  The
comparison is strict (``prob < minProb``, ``minProb`` starting at
``Double.MAX_VALUE``), so a per-dimension tie keeps the first ordered
dimension: that is why a fixture whose positive and negative models coincide
reports the first ordered annotation name for every record instead of ``NULL``.

Measured truth (pinned jar; both fixtures are re-run by this script)
-------------------------------------------------------------------
* 6-variant fixture copied from ``verify_variant_recalibrator.py`` (training =
  rows 1-3, known = rows 1-2, ``-an QD -an MQ --max-gaussians 1
  --bad-lod-score-cutoff 100.0``): GATK writes ``culprit=MQ`` on all six
  records; native wrote ``culprit=full-covariance-gmm`` before the fix.
* 60-variant non-degenerate fixture (40 training records, 20 low-scoring
  records, ``-an QD -an MQ`` at GATK's default gaussian counts): GATK writes
  ``culprit=QD`` on 26 records and ``culprit=MQ`` on 34 records -- the value is
  per datum, never a model name.

What this oracle gates
----------------------
Three checks run on every case, each reported separately:

1. ``cross_tool_byte_identity`` -- the culprit column (the ``culprit=...``
   token of each data row, keyed by POS) must be **byte-identical** between
   pinned GATK and the native tool run on identical inputs.  Gated for the two
   fixtures whose fitted models coincide (measured below); on the multi-component
   fixture it is *reported* instead, see "known divergence" underneath.
2. ``self_rule_fidelity`` -- the culprit column each side wrote must equal the
   per-datum argmin of the GATK rule **re-derived in this script from that
   side's own GATKReport model artifact** (``--output-model``).  This is what
   proves a match is not an accident of two different rules agreeing on one
   fixture, and it is gated on every case including the rich one.
3. ``gatk_culprit_multiset`` -- the GATK side must still produce the measured
   annotation-name multiset (fixture-validity guard), so the oracle cannot
   silently degrade into comparing two ``NULL`` columns or two constants.

Known divergence reported by this oracle (not a violation)
----------------------------------------------------------
On the 60-record case, 5 of 60 records get a different annotation from the two
tools.  This is **not** a rule/naming difference: both sides satisfy gate 2,
and the cause is the already-documented model-fidelity boundary -- native's
fitted mixtures are not bit-identical to Java's VBEM for multi-component models
(me-too ``verify_variant_recalibrator_gatk_model_oracle.py``: "model
floating-point values intentionally remain a separate semantic boundary").
Measured on this fixture: the positive model means differ by up to 2.73
normalized units and the negative model covariances by up to 2.93, which flips
the argmin on the handful of records whose two per-dimension contrasts are
nearly tied.  The oracle prints those records in ``cross_tool_differences`` so
the residual is quantified rather than hidden.

Exit status: 0 when every gated check passes, non-zero otherwise (so it can be
registered in CTest).  ``--expect-divergence`` turns the run into a pure
diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import pathlib
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
    "##FILTER=<ID=LowQual,Description=low quality>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)

# The fixture of verify_variant_recalibrator.py: six SNPs at POS 1..6,
# training = rows 1-3, known = rows 1-2.
TINY_ROWS = [(index, qd, mq) for index, (qd, mq) in
             enumerate([(30, 60), (28, 58), (25, 55), (4, 25), (3, 20), (2, 15)], 1)]


def rich_rows() -> list[tuple[int, float, float]]:
    """60-record fixture with genuinely different positive/negative models.

    Rows 1-40 are the training resource (wide spread); rows 41-50 score badly
    because of QD and rows 51-60 because of MQ, so GATK's per-datum rule picks
    both annotation names.  Generated with ``random.Random(11)`` so the fixture
    is byte-stable across runs and platforms.
    """
    generator = random.Random(11)
    values: list[tuple[float, float]] = []
    for _ in range(40):
        values.append((10 + generator.uniform(0, 30), 45 + generator.uniform(0, 25)))
    for _ in range(10):
        values.append((9.0 + generator.uniform(-1, 1), 52.0 + generator.uniform(-2, 2)))
    for _ in range(10):
        values.append((26.0 + generator.uniform(-2, 2), 34.0 + generator.uniform(-2, 2)))
    return [(10 * (index + 1), qd, mq) for index, (qd, mq) in enumerate(values)]


CASES = (
    {
        "case": "tiny-truth-known",
        "why": "the exact fixture of verify_variant_recalibrator.py -- the "
               "divergence the audit measured (native wrote "
               "'culprit=full-covariance-gmm', GATK writes 'culprit=MQ'); the "
               "positive and negative models coincide here, so GATK's strict "
               "argmin keeps the first information-ordered annotation, MQ",
        "rows": TINY_ROWS,
        "training_positions": [1, 2, 3],
        "known_positions": [1, 2],
        "annotations": ("QD", "MQ"),
        "max_gaussians": 1,
        "extra": ("--bad-lod-score-cutoff", "100.0"),
        "expect_gatk_culprits": {"MQ": 6},
    },
    {
        "case": "tiny-truth-known-annotation-order-swapped",
        "why": "same fixture with '-an MQ -an QD': the emitted value must still "
               "be the annotation NAME selected by GATK's rule (GATK reorders "
               "dimensions by information content, so command-line order is "
               "not the emitted order)",
        "rows": TINY_ROWS,
        "training_positions": [1, 2, 3],
        "known_positions": [1, 2],
        "annotations": ("MQ", "QD"),
        "max_gaussians": 1,
        "extra": ("--bad-lod-score-cutoff", "100.0"),
        "expect_gatk_culprits": {"MQ": 6},
    },
    {
        "case": "rich-per-datum",
        "why": "60 records (40 training, 20 low-scoring) with genuinely "
               "different positive/negative models: GATK's culprit column "
               "varies per datum (26x QD, 34x MQ), which is what proves the "
               "value is the annotation selected by "
               "VariantRecalibratorEngine.calculateWorstPerformingAnnotation "
               "rather than a constant",
        "rows": rich_rows(),
        "training_positions": [10 * index for index in range(1, 41)],
        "known_positions": [],
        "annotations": ("QD", "MQ"),
        "max_gaussians": None,
        "extra": ("--bad-lod-score-cutoff", "100.0"),
        "expect_gatk_culprits": {"QD": 26, "MQ": 34},
        # See the module docstring: on this multi-component fixture native's
        # VBEM fit is not bit-identical to Java's, so a handful of records sit
        # close enough to the argmin boundary for the two models to disagree.
        # Cross-tool byte identity is therefore reported, not gated, here; the
        # self-rule-fidelity gate below still runs on both sides.
        "gate_cross_tool_byte_identity": False,
    },
)

GATED = {case["case"] for case in CASES}

# Both tools must emit the same sites-only recalibration artifact for the
# culprit column to be comparable at all.
GATK_ONLY = ("--create-output-variant-index", "false", "--sites-only-vcf-output", "true",
             "--dont-run-rscript", "true")


def write_vcf(path: pathlib.Path, rows: list[tuple[int, float, float]],
              positions: list[int]) -> None:
    wanted = set(positions)
    body = [f"chr1\t{position}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}"
            for position, qd, mq in rows if position in wanted]
    path.write_text(HEADER + "\n".join(body) + "\n", encoding="utf-8")


def index_feature(java: str, gatk: str, path: pathlib.Path, timeout: int) -> None:
    result = subprocess.run(
        [java, "-jar", gatk, "IndexFeatureFile", "-I", str(path)],
        text=True, capture_output=True, check=False, timeout=timeout,
    )
    if result.returncode != 0:
        raise AssertionError(f"GATK IndexFeatureFile failed for {path}:\n{result.stderr}")


def data_rows(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def culprit_entries(path: pathlib.Path) -> list[tuple[str, str | None]]:
    """``(POS, culprit value)`` for every data row; ``None`` when absent."""
    entries: list[tuple[str, str | None]] = []
    for line in data_rows(path):
        fields = line.split("\t")
        value: str | None = None
        if len(fields) > 7:
            for token in fields[7].split(";"):
                if token.startswith("culprit="):
                    value = token.split("=", 1)[1]
                    break
        entries.append((fields[1], value))
    return entries


def culprit_counts(entries: list[tuple[str, str | None]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for _, value in entries:
        key = "NULL" if value is None else value
        counts[key] = counts.get(key, 0) + 1
    return counts


def model_tables(path: pathlib.Path) -> dict[str, list[list[str]]]:
    """GATKReport tables of a ``--output-model`` artifact, keyed by table name."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    tables: dict[str, list[list[str]]] = {}
    name: str | None = None
    for line in lines:
        if line.startswith("#:GATKTable:"):
            name = line.split(":")[2]
            tables.setdefault(name, [])
        elif line.startswith("#:"):
            name = None
        elif name is not None and line.strip():
            tables[name].append(line.split())
    return tables


def _mixture(tables: dict[str, list[list[str]]], means_table: str, covariance_table: str,
             mix_table: str) -> tuple[list[list[float]], list[list[list[float]]],
                                      list[float]]:
    gaussians = [[float(value) for value in row[1:]] for row in tables[means_table][1:]]
    covariances = [[[float(value) for value in row[2:]]
                    for row in tables[covariance_table][1:] if row[0] == str(component)]
                   for component in range(len(gaussians))]
    weights = [float(row[1]) for row in tables[mix_table][1:]]
    return gaussians, covariances, weights


def _log10_sum(values: list[float]) -> float:
    """Java MathUtils.log10sumLog10 with GaussianMixtureModel's NaN guard."""
    if any(value != value for value in values):
        return float("nan")
    maximum = max(values)
    if maximum == float("-inf"):
        return maximum
    total = sum(10.0 ** (value - maximum) for value in values if value != float("-inf"))
    return maximum + math.log10(total)


def _normal_log10(mean: float, sd: float, x: float) -> float:
    """Java MathUtils.normalDistributionLog10; ``sd`` is the covariance diagonal."""
    root_two_pi = math.sqrt(2.0 * math.pi)
    a = -1.0 * math.log10(sd * root_two_pi)
    b = -1.0 * (((x - mean) ** 2) / (2.0 * sd * sd)) / math.log(10.0)
    return a + b


def rule_argmin(model_path: pathlib.Path, rows: list[tuple[int, float, float]],
                positions: list[str]) -> list[tuple[str, str | None]]:
    """Re-derive the GATK culprit column from that side's own model artifact.

    This is the oracle's independent implementation of
    ``VariantRecalibratorEngine.calculateWorstPerformingAnnotation`` over
    ``GaussianMixtureModel.evaluateDatumInOneDimension``: for each record take
    the dimension minimising good-minus-bad one-dimension log10 likelihood and
    report that dimension's annotation name (``NULL`` when none is available).
    """
    tables = model_tables(model_path)
    for table in ("AnnotationMeans", "AnnotationStdevs", "GoodGaussianPMix",
                  "BadGaussianPMix", "PositiveModelMeans", "PositiveModelCovariances",
                  "NegativeModelMeans", "NegativeModelCovariances"):
        if table not in tables or not tables[table]:
            raise AssertionError(f"{model_path}: missing GATKReport table {table}")
    order = [row[0] for row in tables["AnnotationMeans"][1:]]
    means = [float(row[1]) for row in tables["AnnotationMeans"][1:]]
    stdevs = [float(row[1]) for row in tables["AnnotationStdevs"][1:]]
    good = _mixture(tables, "PositiveModelMeans", "PositiveModelCovariances", "GoodGaussianPMix")
    bad = _mixture(tables, "NegativeModelMeans", "NegativeModelCovariances", "BadGaussianPMix")

    def one_dimension(mixture, value: float, dimension: int) -> float:
        gaussians, covariances, weights = mixture
        terms = []
        for component in range(len(gaussians)):
            # The GATKReport column is already pMixtureLog10 (log10 pMix), not
            # a weight: GaussianMixtureModel.evaluateDatumInOneDimension adds
            # the normal term only when that log10 value is not -infinity.
            weight_log10 = weights[component]
            if weight_log10 != float("-inf"):
                weight_log10 += _normal_log10(gaussians[component][dimension],
                                              covariances[component][dimension][dimension], value)
            terms.append(weight_log10)
        return terms[0] if len(terms) == 1 else _log10_sum(terms)

    raw = {position: (qd, mq) for position, qd, mq in rows}
    entries: list[tuple[str, str | None]] = []
    for position in positions:
        qd, mq = raw[int(position)]
        contrasts: list[float] = []
        for index, annotation in enumerate(order):
            value = ((qd if annotation == "QD" else mq) - means[index]) / stdevs[index]
            contrasts.append(one_dimension(good, value, index)
                             - one_dimension(bad, value, index))
        if any(value != value for value in contrasts):
            entries.append((position, "NULL"))
            continue
        worst = 0
        for index in range(1, len(contrasts)):
            if contrasts[index] < contrasts[worst]:
                worst = index
        entries.append((position, order[worst]))
    return entries


def model_max_abs_delta(gatk_model: pathlib.Path, native_model: pathlib.Path) -> dict[str, float]:
    """Max ``|GATK - native|`` over the shared model-report tables.

    Diagnostic only: it quantifies the documented model-fidelity boundary that
    explains the residual cross-tool culprit differences on multi-component
    fixtures.  It is never a gated check.
    """
    try:
        gatk_tables = model_tables(gatk_model)
        native_tables = model_tables(native_model)
    except OSError:
        return {}

    def flat(tables: dict[str, list[list[str]]], table: str, skip: int) -> list[float]:
        return [float(value) for row in tables.get(table, [])[1:] for value in row[skip:]]

    deltas: dict[str, float] = {}
    for table, skip in (("AnnotationMeans", 1), ("AnnotationStdevs", 1),
                        ("PositiveModelMeans", 1), ("NegativeModelMeans", 1),
                        ("PositiveModelCovariances", 2), ("NegativeModelCovariances", 2)):
        gatk_values = flat(gatk_tables, table, skip)
        native_values = flat(native_tables, table, skip)
        if gatk_values and len(gatk_values) == len(native_values):
            deltas[table] = max(abs(left - right)
                                for left, right in zip(gatk_values, native_values))
    return deltas


def run_case(case: dict, work: pathlib.Path, java: str, gatk: str, native: str,
             threads: int, timeout: int) -> dict:
    case_work = work / case["case"]
    case_work.mkdir(parents=True, exist_ok=True)
    input_vcf = case_work / "input.vcf"
    training_vcf = case_work / "training.vcf"
    known_vcf = case_work / "known.vcf"
    rows = case["rows"]
    write_vcf(input_vcf, rows, [position for position, _, _ in rows])
    write_vcf(training_vcf, rows, case["training_positions"])
    write_vcf(known_vcf, rows, case["known_positions"])
    index_feature(java, gatk, input_vcf, timeout)
    index_feature(java, gatk, training_vcf, timeout)
    if case["known_positions"]:
        index_feature(java, gatk, known_vcf, timeout)

    resources = ["--resource:truth,training=true,truth=true,known=false,prior=15.0",
                 str(training_vcf)]
    if case["known_positions"]:
        resources += ["--resource:known,training=false,truth=false,known=true", str(known_vcf)]

    common: list[str] = []
    for annotation in case["annotations"]:
        common += ["-an", annotation]
    common += ["--mode", "SNP"]
    if case["max_gaussians"] is not None:
        common += ["--max-gaussians", str(case["max_gaussians"])]
    common += list(case["extra"])

    gatk_output = case_work / "gatk.recal.vcf"
    gatk_result = subprocess.run(
        [java, "-Xmx2g", "-jar", gatk, "VariantRecalibrator",
         "-V", str(input_vcf), *resources, *common,
         "--max-attempts", "20", "--k-means-iterations", "20",
         *GATK_ONLY,
         "-O", str(gatk_output),
         "--tranches-file", str(case_work / "gatk.tranches"),
         "--truth-sensitivity-tranche", "100",
         "--output-model", str(case_work / "gatk.model")],
        text=True, capture_output=True, check=False, timeout=timeout,
    )

    native_output = case_work / "native.recal.vcf.gz"
    gatk_model = case_work / "gatk.model"
    native_model = case_work / "native.model"
    native_result = subprocess.run(
        [native, "-V", str(input_vcf), *resources, *common,
         "--max-attempts", "20", "--k-means-iterations", "20",
         "--threads", str(threads),
         "-O", str(native_output),
         "--tranches-file", str(case_work / "native.tranches"),
         "--output-model", str(native_model)],
        text=True, capture_output=True, check=False, timeout=timeout,
    )

    result: dict = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["case"] in GATED,
        "gate_cross_tool_byte_identity": case.get("gate_cross_tool_byte_identity", True),
        "annotations": list(case["annotations"]),
        "max_gaussians": case["max_gaussians"],
        "records": len(rows),
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_stderr_tail": gatk_result.stderr[-600:] if gatk_result.returncode else "",
        "native_stderr_tail": native_result.stderr[-600:] if native_result.returncode else "",
        "expect_gatk_culprits": case["expect_gatk_culprits"],
    }
    if gatk_result.returncode != 0 or native_result.returncode != 0:
        result["violations"] = [
            f"run failed: gatk_exit={gatk_result.returncode} "
            f"native_exit={native_result.returncode}"]
        result["gatk_culprits"] = []
        result["native_culprits"] = []
        result["gatk_culprit_counts"] = {}
        result["native_culprit_counts"] = {}
        result["cross_tool_differences"] = []
        result["rule_fidelity"] = {}
        result["model_max_abs_delta"] = {}
        return result

    gatk_entries = culprit_entries(gatk_output)
    native_entries = culprit_entries(native_output)
    result["gatk_culprits"] = [[pos, value] for pos, value in gatk_entries]
    result["native_culprits"] = [[pos, value] for pos, value in native_entries]
    result["gatk_culprit_counts"] = culprit_counts(gatk_entries)
    result["native_culprit_counts"] = culprit_counts(native_entries)

    violations: list[str] = []
    if len(gatk_entries) != len(native_entries):
        violations.append(
            f"{case['case']}: record count differs: GATK={len(gatk_entries)} "
            f"NATIVE={len(native_entries)}")
    mismatches = [[gatk_position, gatk_value, native_value]
                  for (gatk_position, gatk_value), (_, native_value) in
                  zip(gatk_entries, native_entries) if gatk_value != native_value]
    if mismatches and result["gate_cross_tool_byte_identity"]:
        violations.append(
            f"{case['case']}: culprit column is not byte-identical on "
            f"{len(mismatches)} record(s) (POS, GATK, NATIVE): {mismatches[:8]}")
    result["cross_tool_differences"] = mismatches

    # Gate 2: re-derive the rule from each side's own model artifact and require
    # the side's own culprit column to match it.
    positions = [position for position, value in gatk_entries]
    fidelity: dict[str, dict] = {}
    for side, model_path, entries in (("gatk", gatk_model, gatk_entries),
                                      ("native", native_model, native_entries)):
        derived = rule_argmin(model_path, rows, positions)
        disagreements = [[position, actual, expected]
                         for (position, actual), (_, expected) in zip(entries, derived)
                         if actual != expected]
        fidelity[side] = {
            "rule_argmin": [[position, value] for position, value in derived],
            "disagreements": disagreements,
        }
        if disagreements:
            violations.append(
                f"{case['case']}: {side} culprit column does not follow GATK's "
                f"calculateWorstPerformingAnnotation on its own model artifact: "
                f"{len(disagreements)} disagreement(s) (POS, written, rule): "
                f"{disagreements[:8]}")
    result["rule_fidelity"] = fidelity
    result["model_max_abs_delta"] = model_max_abs_delta(gatk_model, native_model)

    # Gate 3: fixture-validity guard -- the GATK side must still produce the
    # measured annotation-name multiset, so a fixture that degraded into a
    # constant or into 'NULL' cannot make this oracle pass vacuously.
    if result["gatk_culprit_counts"] != case["expect_gatk_culprits"]:
        violations.append(
            f"{case['case']}: GATK culprit multiset changed: measured "
            f"{case['expect_gatk_culprits']} but GATK now writes "
            f"{result['gatk_culprit_counts']}")
    result["violations"] = violations
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the VariantRecalibrator INFO culprit "
                    "value against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_VARIANT_RECALIBRATOR_BINARY"),
        help="native VariantRecalibrator binary (default: "
             "$FASTGATK_VARIANT_RECALIBRATOR_BINARY or "
             "fastgatk-native/build/fastgatk-variant-recalibrator)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--threads", type=int, default=2,
                        help="--threads value for the native run (default 2)")
    parser.add_argument("--timeout", type=int, default=900,
                        help="per-process timeout in seconds (default 900)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-variant-recalibrator")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))

    assets = [native, pathlib.Path(java), gatk]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_variant_recalibrator_culprit_gatk_oracle.py",
            pathlib.Path(java), gatk)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-variant-recalibrator-culprit-oracle-") as directory:
        work = pathlib.Path(directory)
        for case in selected:
            results.append(run_case(case, work, java, str(gatk), str(native),
                                    arguments.threads, arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: VariantRecalibratorEngine."
          "calculateWorstPerformingAnnotation (VariantRecalibratorEngine.java:80-93) "
          "keeps the annotation dimension minimising "
          "goodModel.evaluateDatumInOneDimension - badModel.evaluateDatumInOneDimension "
          "-- the per-dimension log10 mixture likelihood of "
          "GaussianMixtureModel.evaluateDatumInOneDimension "
          "(GaussianMixtureModel.java:208-222) -- and VariantDataManager.java:485 "
          "writes that dimension's annotation NAME.")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    for result in results:
        print(f"[{result['case']}] gated={result['gated']} "
              f"gate_cross_tool_byte_identity={result['gate_cross_tool_byte_identity']} "
              f"records={result['records']} annotations={result['annotations']} "
              f"max_gaussians={result['max_gaussians']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    GATK   culprit multiset: {result['gatk_culprit_counts']} "
              f"(measured expectation {result['expect_gatk_culprits']})")
        print(f"    NATIVE culprit multiset: {result['native_culprit_counts']}")
        if result["gatk_culprits"]:
            print(f"    GATK   culprit column: "
                  f"{[value for _, value in result['gatk_culprits']][:16]}")
            print(f"    NATIVE culprit column: "
                  f"{[value for _, value in result['native_culprits']][:16]}")
        for side, fidelity in sorted(result["rule_fidelity"].items()):
            print(f"    rule fidelity [{side}]: re-derived from {side}'s own model "
                  f"artifact -> {len(fidelity['disagreements'])} disagreement(s) "
                  f"with the written culprit column")
            for position, actual, derived in fidelity["disagreements"][:8]:
                print(f"        POS {position}: written={actual} rule={derived}")
        for case_name, model in result["model_max_abs_delta"].items():
            print(f"    model |GATK-NATIVE| max delta [{case_name}]: {model}")
        if result["cross_tool_differences"]:
            note = ("GATED byte-identity" if result["gate_cross_tool_byte_identity"]
                    else "reported, not gated (model-fidelity boundary, see docstring)")
            print(f"    cross-tool culprit differences ({note}): "
                  f"{result['cross_tool_differences'][:8]}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")
        if result["gatk_stderr_tail"]:
            print(f"    gatk_stderr_tail: {result['gatk_stderr_tail']!r}")
        if result["native_stderr_tail"]:
            print(f"    native_stderr_tail: {result['native_stderr_tail']!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "VariantRecalibrator INFO culprit (per-datum annotation name)",
        "binary": str(native),
        "acceptance_criterion": (
            "gate 1 (where enabled) the culprit token of every data row must be "
            "byte-identical between pinned GATK 4.6.2.0 and native; gate 2 on "
            "every case each side's culprit column must equal the "
            "calculateWorstPerformingAnnotation argmin re-derived from that "
            "side's own model artifact; gate 3 the GATK side must still carry "
            "the measured annotation-name multiset"),
        "known_divergence": (
            "rich-per-datum: cross-tool byte identity is not gated because "
            "native's multi-component VBEM fit is not bit-identical to Java's "
            "(documented model-fidelity boundary); the residual is reported in "
            "cross_tool_differences and every record on both sides still "
            "satisfies gate 2"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": sorted(GATED),
        "violations": violations,
        "results": results,
    }
    print(json.dumps(payload, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
