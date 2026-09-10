#!/usr/bin/env python3
"""Compare MarkDuplicates with the pinned GATK/Picard implementation.

``verify_mark_duplicates.py`` is the broad native contract (spill, restart,
index and malformed-option checks).  This verifier is deliberately separate
and only asks the release-pinned Java tool for the semantic slice that can be
compared without depending on BAM compression or dynamic command provenance:
all SAM core fields, duplicate flags, optional tags (including ``DT``/``PG``),
Picard DuplicationMetrics rows, and the duplicate-set histogram.

The fixture intentionally contains stale duplicate state and secondary data.
That catches a particularly dangerous replacement bug: a native pass that
marks new duplicates correctly but leaves an old BAM_FDUP or DT tag on a
representative/read excluded from duplicate-key processing.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
@RG\tID:rg1\tSM:S1\tLB:lib1
@PG\tID:upstream\tPN:upstream\tVN:1
PAIR_A:1:1101:100:100\t99\tchr1\t11\t60\t4M\t=\t31\t24\tACGT\tIIII\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_B:1:1101:101:101\t99\tchr1\t11\t20\t4M\t=\t31\t24\tACGT\t!!!!\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
PAIR_C:1:1101:102:102\t0\tchr1\t21\t30\t4M\t*\t0\t0\tAAAA\tIIII\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_A:1:1101:100:100\t147\tchr1\t31\t60\t4M\t=\t11\t-24\tTGCA\t!!!!\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_B:1:1101:101:101\t1171\tchr1\t31\t20\t4M\t=\t11\t-24\tTGCA\t!!!!\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
SECONDARY:1:1101:103:103\t2304\tchr1\t41\t60\t4M\t*\t0\t0\tCCCC\tIIII\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_X:1:1101:999:999\t99\tchr1\t51\t60\t4M\t=\t71\t24\tACGT\tIIII\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_Y:1:1101:888:888\t99\tchr1\t51\t20\t4M\t=\t71\t24\tACGT\t!!!!\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
PAIR_X:1:1101:999:999\t147\tchr1\t71\t60\t4M\t=\t51\t-24\tTGCA\t!!!!\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
PAIR_Y:1:1101:888:888\t147\tchr1\t71\t20\t4M\t=\t51\t-24\tTGCA\tIIII\tMC:Z:4M\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
single-high\t0\tchr1\t81\t60\t4M\t*\t0\t0\tGGGG\tIIII\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
single-low\t0\tchr1\t81\t20\t4M\t*\t0\t0\tGGGG\t!!!!\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
qcfail-peer\t0\tchr1\t91\t60\t4M\t*\t0\t0\tTTTT\tIIII\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:SQ
qcfail\t512\tchr1\t91\t60\t4M\t*\t0\t0\tTTTT\tIIII\tPG:Z:upstream\tRG:Z:rg1\tDT:Z:LB
"""


def fail(command: list[str], result: subprocess.CompletedProcess[str]) -> None:
    raise AssertionError(
        "command failed ({}): {}\nstdout:\n{}\nstderr:\n{}".format(
            " ".join(command), result.returncode, result.stdout, result.stderr))


def run(command: list[str], *, capture: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=capture, check=False)
    if result.returncode != 0:
        fail(command, result)
    return result


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")] 


def tags(row: list[str]) -> dict[str, str]:
    # Compare tags as a map so a writer-specific optional-field order does not
    # obscure a semantic mismatch.  The duplicate-tag and PG values remain
    # part of the comparison; no provenance field is silently discarded.
    return {field[:2]: field[3:] for field in row[11:] if len(field) >= 4}


def semantic_records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    for row in records(path):
        if len(row) < 11:
            raise AssertionError(f"malformed SAM record in {path}: {row!r}")
        result.append(tuple(row[:11]) + (json.dumps(tags(row), sort_keys=True),))
    return result


def pg_header(path: Path) -> list[dict[str, str]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("@PG\t"):
            continue
        result.append({field[:2]: field[3:] for field in line.split("\t")[1:]
                       if len(field) >= 4})
    return result


def metric_rows(path: Path) -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("##"):
            continue
        fields = line.split("\t")
        if fields[0] == "LIBRARY" or len(fields) < 10:
            continue
        rows[fields[0]] = fields
    return rows


def metric_semantics(path: Path) -> dict[str, tuple[object, ...]]:
    result: dict[str, tuple[object, ...]] = {}
    for library, row in metric_rows(path).items():
        values: list[object] = list(row[:8])
        try:
            values.append(float(row[8]))
        except ValueError as error:
            raise AssertionError(f"invalid PERCENT_DUPLICATION for {library}: {row[8]}") from error
        # Picard writes a blank estimated library size for an unidentifiable
        # population; otherwise it is an integer-valued decimal.
        values.append(None if not row[9] else int(float(row[9])))
        result[library] = tuple(values)
    return result


def histogram_semantics(path: Path) -> list[tuple[float, int, int, int]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        # Picard 4.6 emits BIN/CoverageMult before the three count columns;
        # the native writer emits the compact four-column form.  Both encode
        # the same duplicate-set histogram, so discard the redundant
        # coverage-multiplier column and empty bins.
        start = next(index for index, line in enumerate(lines)
                     if line in ("set_size\tall_sets\toptical_sets\tnon_optical_sets",
                                 "BIN\tCoverageMult\tall_sets\toptical_sets\tnon_optical_sets")) + 1
    except StopIteration:
        return []
    result = []
    for line in lines[start:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) == 5:
            fields = [fields[0], fields[2], fields[3], fields[4]]
        if len(fields) != 4:
            continue
        counts = tuple(int(value) for value in fields[1:])
        if any(counts):
            result.append((float(fields[0]), *counts))
    return result


def compare(java_output: Path, native_output: Path,
            java_metrics: Path, native_metrics: Path) -> None:
    java_records = semantic_records(java_output)
    native_records = semantic_records(native_output)
    if java_records != native_records:
        for index, (expected, actual) in enumerate(zip(java_records, native_records)):
            if expected != actual:
                raise AssertionError({"record_index": index, "picard": expected, "native": actual})
        raise AssertionError({"record_count": (len(java_records), len(native_records))})
    if metric_semantics(java_metrics) != metric_semantics(native_metrics):
        raise AssertionError({"metrics": (metric_semantics(java_metrics), metric_semantics(native_metrics))})
    if histogram_semantics(java_metrics) != histogram_semantics(native_metrics):
        raise AssertionError({"histogram": (histogram_semantics(java_metrics), histogram_semantics(native_metrics))})
    if not any(header.get("ID") == "MarkDuplicates" for header in pg_header(native_output)):
        raise AssertionError("native output is missing the MarkDuplicates @PG header")


def invoke_pair(root: Path, source: Path, work: Path, label: str,
                *, remove_duplicates: bool = False,
                remove_sequencing_duplicates: bool = False,
                tagging_policy: str = "All") -> tuple[int, int]:
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get("FASTGATK_MARK_DUPLICATES_BINARY",
                                os.environ.get("FASTGATK_NATIVE_BUILD",
                                               root / "fastgatk-native/build")))
    if native.is_dir():
        native /= "fastgatk-mark-duplicates"
    java_output = work / f"{label}.picard.sam"
    native_output = work / f"{label}.native.sam"
    java_metrics = work / f"{label}.picard.metrics"
    native_metrics = work / f"{label}.native.metrics"
    java_command = [str(java), "-jar", str(gatk), "MarkDuplicates",
                    "-I", str(source), "-O", str(java_output), "-M", str(java_metrics),
                    "--ASSUME_SORTED", "true", "--CREATE_INDEX", "false",
                    "--TAGGING_POLICY", tagging_policy]
    if remove_duplicates:
        java_command += ["--REMOVE_DUPLICATES", "true"]
    if remove_sequencing_duplicates:
        java_command += ["--REMOVE_SEQUENCING_DUPLICATES", "true"]
    run(java_command)
    native_command = [str(native), "-I", str(source), "-O", str(native_output),
                      "--metrics-file", str(native_metrics),
                      "--create-output-bam-index=false", "--tagging-policy", tagging_policy]
    if remove_duplicates:
        native_command.append("--remove-duplicates")
    if remove_sequencing_duplicates:
        native_command.append("--remove-sequencing-duplicates")
    result = run(native_command)
    summary = json.loads(result.stdout.strip().splitlines()[-1])
    compare(java_output, native_output, java_metrics, native_metrics)
    expected_records = len(records(java_output))
    actual_records = len(records(native_output))
    if expected_records != actual_records:
        raise AssertionError({"output_records": (expected_records, actual_records)})
    return expected_records, summary["duplicate_records"]


def verify_create_index_default(root: Path, source: Path, work: Path) -> None:
    """Pin Picard's CREATE_INDEX=false default and its explicit true path.

    This is intentionally a side-effect assertion: both engines already
    compare record/metrics semantics above, but an unexpected index changes
    the output set observed by Nextflow and downstream scatter tasks.
    """
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get("FASTGATK_MARK_DUPLICATES_BINARY",
                                os.environ.get("FASTGATK_NATIVE_BUILD",
                                               root / "fastgatk-native/build")))
    if native.is_dir():
        native /= "fastgatk-mark-duplicates"

    java_default = work / "index-default.picard.bam"
    native_default = work / "index-default.native.bam"
    java_default_metrics = work / "index-default.picard.metrics"
    native_default_metrics = work / "index-default.native.metrics"
    run([str(java), "-jar", str(gatk), "MarkDuplicates", "-I", str(source),
         "-O", str(java_default), "-M", str(java_default_metrics),
         "--ASSUME_SORTED", "true"])
    run([str(native), "-I", str(source), "-O", str(native_default),
         "--metrics-file", str(native_default_metrics)])
    java_index = java_default.with_suffix(".bai")
    native_index = Path(f"{native_default}.bai")
    if java_index.exists() or native_index.exists():
        raise AssertionError({"create_index_default": False,
                              "picard_index": java_index.exists(),
                              "native_index": native_index.exists()})

    java_explicit = work / "index-explicit.picard.bam"
    native_explicit = work / "index-explicit.native.bam"
    java_explicit_metrics = work / "index-explicit.picard.metrics"
    native_explicit_metrics = work / "index-explicit.native.metrics"
    run([str(java), "-jar", str(gatk), "MarkDuplicates", "-I", str(source),
         "-O", str(java_explicit), "-M", str(java_explicit_metrics),
         "--ASSUME_SORTED", "true", "--CREATE_INDEX", "true"])
    run([str(native), "-I", str(source), "-O", str(native_explicit),
         "--metrics-file", str(native_explicit_metrics),
         "--create-output-bam-index", "true"])
    java_index = java_explicit.with_suffix(".bai")
    native_index = Path(f"{native_explicit}.bai")
    if not java_index.exists() or not native_index.exists():
        raise AssertionError({"create_index_explicit": True,
                              "picard_index": java_index.exists(),
                              "native_index": native_index.exists()})


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get("FASTGATK_MARK_DUPLICATES_BINARY",
                                os.environ.get("FASTGATK_NATIVE_BUILD",
                                               root / "fastgatk-native/build")))
    if native.is_dir():
        native /= "fastgatk-mark-duplicates"
    required = (java, gatk, native)
    if not all(path.is_file() and os.access(path, os.X_OK) for path in (java, native)) or not gatk.is_file():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned GATK/Picard MarkDuplicates oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned GATK/Picard oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-mark-duplicates-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        source.write_text(SAM, encoding="utf-8")
        ordinary_records, duplicate_records = invoke_pair(root, source, work, "ordinary")
        removed_records, removed_duplicates = invoke_pair(
            root, source, work, "removed", remove_duplicates=True, tagging_policy="DontTag")
        sequencing_records, sequencing_duplicates = invoke_pair(
            root, source, work, "sequencing-removed",
            remove_sequencing_duplicates=True, tagging_policy="All")
        verify_create_index_default(root, source, work)
        if ordinary_records != 14 or duplicate_records != 6:
            raise AssertionError({"ordinary_records": ordinary_records,
                                  "duplicate_records": duplicate_records})
        if removed_records != 8 or removed_duplicates != 6:
            raise AssertionError({"removed_records": removed_records,
                                  "removed_duplicates": removed_duplicates})
        if sequencing_records != 12 or sequencing_duplicates != 6:
            raise AssertionError({"sequencing_records": sequencing_records,
                                  "sequencing_duplicates": sequencing_duplicates})
    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "records_compared": 14, "duplicate_records": 6,
                      "remove_duplicates_records_compared": 8,
                      "remove_sequencing_duplicates_records_compared": 12,
                      "create_index_default": False,
                      "create_index_explicit": True,
                      "metrics_and_histogram_semantically_identical": True,
                      "stale_duplicate_state_cleared": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
