#!/usr/bin/env python3
"""Acceptance oracle for the D2 defect: ``-L``-window-dependent HC annotations.

Defect under test (proven, see
``fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md`` and
``fastgatk-native/evidence/2026-09-10-parallel-audit/track-d-streaming-rootcause.md``):
with **no streaming at all**, on

    -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L <WINDOW>
    --emit-ref-confidence GVCF --max-mnp-distance 1 --threads 1
    --add-output-vcf-command-line false

the native HaplotypeCaller annotation evidence at ``20:10020680``
(``ALT=AT,*,<NON_REF>``) depends on the ``-L`` window start, while the pinned
GATK 4.6.2.0 baseline is invariant there (``RAW_MQandDP=28800,8``,
``SB=0,0,3,3``) for every window in which the row set is otherwise identical.

What this oracle does
---------------------
For every window in the pinned set it runs the pinned GATK 4.6.2.0 jar and the
native binary, then

  (a) asserts the native data lines are byte-identical to the GATK data lines
      (``#``-comment/provenance lines and header content are ignored), and
  (b) reports the ``POS 10020680`` ``RAW_MQandDP`` / ``SB`` values for both
      sides, per window.

It exits non-zero on any mismatch so it can be registered in CTest, and it
prints the exact values it compares for each window, so a failure is
self-explanatory.

``--expect-divergence`` turns the same run into a pure diagnostic report and
exits 0: that is how the *current, unfixed* binary can be observed without
turning a test suite red.  The default is the strict assertion.

Window set
----------
``GATING_WINDOWS`` are the windows verified against the pinned GATK (live,
2026-09-10) whose only native divergence is the ``10020680`` annotation row
(``10020381/10020391/10020401/10020411`` currently diverge; the remaining
entries are controls that must stay byte-identical).  ``STRUCTURAL_WINDOWS``
(``10020421``, ``10020431``) were checked because a prior audit reported a third
wrong state there: that report is **refuted** - the pinned GATK itself emits
``RAW_MQandDP=97200,27`` / ``SB=1,2,6,16`` at ``10020680`` for those two windows
and native matches it exactly.  The residual mismatch there is reference-block
granularity (the D3 family: native collapses many ``<NON_REF>`` blocks into one
degenerate ``GQ=0 / PL=0,0,0`` block), which is a different defect, so those
windows are reported and their ``10020680`` row is still gated, but their full
row set is only gated with ``--gate-structural-windows``.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

D2_POSITION = 10020680
D2_ALLELES = "AT,*,<NON_REF>"
WINDOW_END = 10020710

# Window starts whose native divergence against the pinned GATK is (today)
# either nothing or exactly the POS 10020680 annotation row.
GATING_WINDOWS: tuple[int, ...] = (
    10019901,  # window A: 48 rows, byte-identical today (the reference baseline)
    10020201,
    10020301,
    10020351,
    10020381,  # window B: the documented D2 divergence
    10020391,
    10020401,
    10020411,
    10020441,
    10020501,
)

# Windows that additionally expose reference-block granularity divergence (D3).
STRUCTURAL_WINDOWS: tuple[int, ...] = (10020421, 10020431)

# Values observed from the pinned GATK 4.6.2.0 jar on 2026-09-10 for
# POS 10020680 in each window.  Used as a drift guard on the GATK baseline so a
# wrong GATK invocation cannot silently "agree" with a broken native build.
GATK_D2_EXPECTATION: dict[int, dict[str, str]] = {
    10019901: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020201: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020301: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020351: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020381: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020391: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020401: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020411: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020421: {"RAW_MQandDP": "97200,27", "SB": "1,2,6,16"},
    10020431: {"RAW_MQandDP": "97200,27", "SB": "1,2,6,16"},
    10020441: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020451: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020501: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020601: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
    10020661: {"RAW_MQandDP": "28800,8", "SB": "0,0,3,3"},
}


def records(path: Path) -> list[list[str]]:
    """Data rows only: header/provenance lines are ignored by construction."""
    return [line.rstrip("\n").split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def row_at(rows: list[list[str]], position: int) -> list[str] | None:
    for row in rows:
        if len(row) > 1 and row[1].isdigit() and int(row[1]) == position:
            return row
    return None


def info_value(row: list[str] | None, key: str) -> str | None:
    if row is None or len(row) < 8:
        return None
    for field in row[7].split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def format_value(row: list[str] | None, key: str) -> str | None:
    if row is None or len(row) < 10:
        return None
    keys = row[8].split(":")
    if key not in keys:
        return None
    values = row[9].split(":")
    index = keys.index(key)
    return values[index] if index < len(values) else None


def d2_probe(row: list[list[str]]) -> dict[str, object]:
    """The two fields the D2 defect moves, plus the ALT they must belong to."""
    record = row_at(row, D2_POSITION)
    return {
        "present": record is not None,
        "ALT": record[4] if record is not None else None,
        "RAW_MQandDP": info_value(record, "RAW_MQandDP"),
        "SB": format_value(record, "SB"),
    }


def row_field_differences(gatk_row: list[str], native_row: list[str]) -> list[str]:
    """Field-level diff of one VCF row (INFO sub-fields and FORMAT sub-fields)."""
    differences: list[str] = []
    for index in (3, 4, 6):
        if index < len(gatk_row) and index < len(native_row) and gatk_row[index] != native_row[index]:
            label = {3: "REF", 4: "ALT", 6: "QUAL"}[index]
            differences.append(f"{label} {gatk_row[index]!r} -> {native_row[index]!r}")
    gatk_info = {} if len(gatk_row) < 8 else dict(
        field.split("=", 1) if "=" in field else (field, "")
        for field in gatk_row[7].split(";"))
    native_info = {} if len(native_row) < 8 else dict(
        field.split("=", 1) if "=" in field else (field, "")
        for field in native_row[7].split(";"))
    for key in sorted(set(gatk_info) | set(native_info)):
        if gatk_info.get(key) != native_info.get(key):
            differences.append(
                f"INFO {key} {gatk_info.get(key)!r} -> {native_info.get(key)!r}")
    if len(gatk_row) > 9 and len(native_row) > 9:
        if gatk_row[8] != native_row[8]:
            differences.append(f"FORMAT keys {gatk_row[8]!r} -> {native_row[8]!r}")
        keys = gatk_row[8].split(":")
        gatk_values = gatk_row[9].split(":")
        native_values = native_row[9].split(":")
        for index, key in enumerate(keys):
            gatk_value = gatk_values[index] if index < len(gatk_values) else None
            native_value = native_values[index] if index < len(native_values) else None
            if gatk_value != native_value:
                differences.append(f"FORMAT {key} {gatk_value!r} -> {native_value!r}")
    return differences


def divergent_positions(gatk_rows: list[list[str]],
                        native_rows: list[list[str]]) -> list[str]:
    positions: list[str] = []
    for index in range(max(len(gatk_rows), len(native_rows))):
        gatk_row = gatk_rows[index] if index < len(gatk_rows) else None
        native_row = native_rows[index] if index < len(native_rows) else None
        if gatk_row is None:
            positions.append(f"missing-in-gatk@{index}")
        elif native_row is None:
            positions.append(f"missing-in-native@{index}:POS={gatk_row[1]}")
        elif gatk_row != native_row:
            positions.append(native_row[1] if native_row[1] == gatk_row[1]
                             else f"{gatk_row[1]}/{native_row[1]}")
    return positions


def run_gatk(java: str, gatk: Path, reference: Path, bam: Path, window: str,
             output: Path, threads: int) -> None:
    subprocess.run([
        java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
        "-R", str(reference), "-I", str(bam), "-L", window,
        "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1",
        "--native-pair-hmm-threads", str(threads),
        "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
        "-O", str(output),
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_native(native: Path, reference: Path, bam: Path, window: str,
               output: Path, threads: int) -> None:
    subprocess.run([
        str(native),
        "-R", str(reference), "-I", str(bam), "-L", window,
        "--emit-ref-confidence", "GVCF", "--max-mnp-distance", "1",
        "--threads", str(threads),
        "--add-output-vcf-command-line", "false",
        "-O", str(output),
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def compare(label: str, gatk_rows: list[list[str]],
            native_rows: list[list[str]]) -> dict:
    positions = divergent_positions(gatk_rows, native_rows)
    gatk_d2 = d2_probe(gatk_rows)
    native_d2 = d2_probe(native_rows)
    first_row_diff: list[str] = []
    for index in range(min(len(gatk_rows), len(native_rows))):
        if gatk_rows[index] != native_rows[index]:
            first_row_diff = row_field_differences(gatk_rows[index], native_rows[index])
            break
    return {
        "backend": label,
        "gatk_rows": len(gatk_rows),
        "native_rows": len(native_rows),
        "data_rows_byte_identical": gatk_rows == native_rows,
        "divergent_positions": positions,
        "first_divergent_row_field_diff": first_row_diff,
        "gatk_pos_10020680": gatk_d2,
        "native_pos_10020680": native_d2,
        "pos_10020680_identical": gatk_d2 == native_d2,
    }


def print_window_report(comparison: dict, expect_rows: dict | None) -> None:
    window = comparison["window"]
    primary = comparison["backends"][0]
    print(f"[window {window}] gatk_rows={primary['gatk_rows']} "
          f"native_rows={primary['native_rows']} "
          f"data_rows_byte_identical={primary['data_rows_byte_identical']} "
          f"d2_row_identical={primary['pos_10020680_identical']} "
          f"gating={comparison['gating']}")
    for backend in comparison["backends"]:
        positions = backend["divergent_positions"]
        print(f"  [{backend['backend']}] divergent positions: {positions[:10]}"
              + (f" (+{len(positions) - 10} more)" if len(positions) > 10 else ""))
    if expect_rows is not None:
        print(f"  pinned GATK expectation at POS {D2_POSITION}: {expect_rows}")
    for backend in comparison["backends"]:
        gatk_d2 = backend["gatk_pos_10020680"]
        native_d2 = backend["native_pos_10020680"]
        print(f"  [{backend['backend']}] POS {D2_POSITION} ALT={gatk_d2['ALT']}")
        print(f"    GATK   RAW_MQandDP={gatk_d2['RAW_MQandDP']} "
              f"SB={gatk_d2['SB']}  rows={backend['gatk_rows']}")
        print(f"    NATIVE RAW_MQandDP={native_d2['RAW_MQandDP']} "
              f"SB={native_d2['SB']}  rows={backend['native_rows']}")
        if backend["first_divergent_row_field_diff"]:
            print("    first divergent row field diff (GATK -> NATIVE): "
                  + "; ".join(backend["first_divergent_row_field_diff"]))
        if backend["gatk_pos_10020680"]["present"] and gatk_d2["ALT"] != D2_ALLELES:
            print(f"    WARNING: expected ALT {D2_ALLELES!r}, saw {gatk_d2['ALT']!r}")


def parse_windows(specification: str | None) -> tuple[int, ...]:
    if specification is None:
        return GATING_WINDOWS
    starts: list[int] = []
    for token in specification.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            first, last = token.split("-", 1)
            starts.extend(range(int(first), int(last) + 1))
        else:
            starts.append(int(token))
    return tuple(starts)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Acceptance oracle for the -L-window-dependent HC annotation "
                    "(D2) at 20:10020680, against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_HC_BINARY"),
        help="native OpenMP HaplotypeCaller binary "
             "(default: $FASTGATK_HC_BINARY or fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument(
        "--serial-native", default=os.environ.get("FASTGATK_HC_SERIAL_BINARY"),
        help="optional second native backend (build-serial/fastgatk-hc-call); when "
             "given it is compared against GATK with the same rules")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument(
        "--gate-structural-windows", action="store_true",
        help="also assert full-row byte parity for the windows that expose the "
             "separate reference-block-granularity divergence (10020421, 10020431)")
    parser.add_argument(
        "--windows", default=None,
        help="override the window start list, e.g. '10019901,10020381' or "
             "'10020381-10020411'; default is the pinned set below")
    parser.add_argument("--threads", type=int, default=1,
                        help="--threads / --native-pair-hmm-threads value (default 1)")
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(arguments.native) if arguments.native else (
        root / "fastgatk-native/build/fastgatk-hc-call")
    serial = Path(arguments.serial_native) if arguments.serial_native else None
    reference = root / "fixtures/chr20/ref20mnp.fasta"
    bam = root / "fixtures/chr20/mnp.bam"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    gatk = root / ("third_party/gatk-package/gatk-4.6.2.0/"
                   "gatk-package-4.6.2.0-local.jar")

    assets = [native, reference, bam, Path(java), gatk] + ([serial] if serial else [])
    if not all(path.is_file() for path in assets):
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}))
        return 0
    if not Path(java).is_file():
        raise SystemExit(f"missing java: {java}")

    starts = parse_windows(arguments.windows)
    if arguments.windows is None:
        starts = GATING_WINDOWS + STRUCTURAL_WINDOWS
    structural = set(STRUCTURAL_WINDOWS)

    backends: list[tuple[str, Path]] = [("OpenMP", native)]
    if serial is not None:
        backends.append(("Serial", serial))

    results: list[dict] = []
    strict_mode = not arguments.expect_divergence
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-window-oracle-") as directory:
        work = Path(directory)
        for start in starts:
            window = f"20:{start}-{WINDOW_END}"
            gating = (start not in structural) or arguments.gate_structural_windows
            gatk_vcf = work / f"gatk.{start}.vcf"
            run_gatk(java, gatk, reference, bam, window, gatk_vcf, arguments.threads)
            gatk_rows = records(gatk_vcf)
            comparisons = []
            for label, binary in backends:
                native_vcf = work / f"native.{label.lower()}.{start}.vcf"
                run_native(binary, reference, bam, window, native_vcf, arguments.threads)
                comparison = compare(label, gatk_rows, records(native_vcf))
                comparisons.append(comparison)
            expected = GATK_D2_EXPECTATION.get(start)
            gatk_observed = comparisons[0]["gatk_pos_10020680"]
            gatk_drift = expected is not None and (
                gatk_observed["RAW_MQandDP"] != expected["RAW_MQandDP"]
                or gatk_observed["SB"] != expected["SB"])
            result = {
                "window": window,
                "window_start": start,
                "window_end": WINDOW_END,
                "window_class": "structural" if start in structural else "gating",
                "gating": gating,
                "pinned_gatk_expectation_pos_10020680": expected,
                "gatk_baseline_drift": gatk_drift,
                "backends": comparisons,
            }
            results.append(result)

    violations: list[str] = []
    for result in results:
        for backend in result["backends"]:
            label = f"{result['window']} [{backend['backend']}]"
            if not backend["pos_10020680_identical"]:
                violations.append(
                    f"{label}: POS {D2_POSITION} annotation differs "
                    f"(GATK RAW_MQandDP={backend['gatk_pos_10020680']['RAW_MQandDP']} "
                    f"SB={backend['gatk_pos_10020680']['SB']} vs native "
                    f"RAW_MQandDP={backend['native_pos_10020680']['RAW_MQandDP']} "
                    f"SB={backend['native_pos_10020680']['SB']})")
            if result["gating"] and not backend["data_rows_byte_identical"]:
                violations.append(
                    f"{label}: data rows differ at {backend['divergent_positions'][:10]} "
                    f"(gatk_rows={backend['gatk_rows']} "
                    f"native_rows={backend['native_rows']})")
            if result["gatk_baseline_drift"]:
                violations.append(
                    f"{label}: GATK baseline drifted from the pinned expectation "
                    f"{result['pinned_gatk_expectation_pos_10020680']} "
                    f"(observed RAW_MQandDP={backend['gatk_pos_10020680']['RAW_MQandDP']} "
                    f"SB={backend['gatk_pos_10020680']['SB']})")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native, "
          f"fixture fixtures/chr20/mnp.bam, windows '-L 20:<start>-{WINDOW_END}'")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; gate-structural-windows={arguments.gate_structural_windows}; "
          f"threads={arguments.threads}")
    print(f"# backends: {', '.join(label for label, _ in backends)}")
    for result in results:
        print_window_report(result, result["pinned_gatk_expectation_pos_10020680"])
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    divergence_summary = {
        result["window"]: {
            backend["backend"]: {
                "pos_10020680_identical": backend["pos_10020680_identical"],
                "data_rows_byte_identical": backend["data_rows_byte_identical"],
                "gatk_pos_10020680": backend["gatk_pos_10020680"],
                "native_pos_10020680": backend["native_pos_10020680"],
            }
            for backend in result["backends"]
        }
        for result in results
    }
    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "hc -L-window invariance of POS 10020680 annotation evidence (D2)",
        "fixture": "fixtures/chr20/ref20mnp.fasta + fixtures/chr20/mnp.bam",
        "d2_position": D2_POSITION,
        "d2_alleles": D2_ALLELES,
        "acceptance_criterion": (
            f"POS {D2_POSITION} must be byte-identical to pinned GATK in every "
            "window; data rows must be byte-identical for every gating window"),
        "windows": [result["window"] for result in results],
        "gating_windows": [result["window"] for result in results if result["gating"]],
        "strict_mode": strict_mode,
        "violations": violations,
        "per_window": divergence_summary,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
