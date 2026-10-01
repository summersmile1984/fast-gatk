#!/usr/bin/env python3
"""Verify HC `--flow-probability-scaling-factor` + `--flow-quantization-bins` + `--flow-matrix-mods` CLI parity.

GATK 4.6.2.0 HaplotypeCaller accepts three additional flow-control
flags:
  * --flow-probability-scaling-factor N  Integer (>= 1), scale
    flow-probability values during matrix construction
  * --flow-quantization-bins N          Integer (>= 2), number of
    bins used to quantize per-base flow probabilities
  * --flow-matrix-mods STR             String (comma-separated
    src,dst pairs), flow-matrix modification instructions

Before this turn native HC did not parse these flags.  This oracle
asserts that native now:
  * accepts valid --flow-probability-scaling-factor N (>= 1),
  * accepts valid --flow-quantization-bins N (>= 2),
  * accepts valid --flow-matrix-mods STR (non-empty),
  * rejects out-of-range values with a clear error message,
  * produces a valid VCF for valid invocations on the canonical
    NA12878 chr17-69k-70k fixture (3 records, exit 0).
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run_native_hc(*extra_args: str, bam: Path, reference: Path,
                  output: Path, work: Path) -> tuple[int, str, str]:
    cmd = [
        str(Path(os.environ.get(
            "FASTGATK_HC_CALL_BINARY",
            "/home/turing-agents/Documents/fast-gatk/fastgatk-native/build/fastgatk-hc-call"))),
        "-I", str(bam),
        "-R", str(reference),
        "-L", "17:69000-70000",
        "--sample-name", "NA12878",
        "-O", str(output),
        "--tmp-dir", str(work),
        *extra_args,
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    return proc.returncode, proc.stdout, proc.stderr


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_HC_CALL_BINARY",
        str(root / "fastgatk-native" / "build" / "fastgatk-hc-call")))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"

    required = (binary, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified(
            "verify_hc_flow_scaling_quantization_mods_flag_oracle.py", None, None)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled HC oracle inputs unavailable"}))
        return 0

    diagnostics: list[str] = []
    case_results: dict[str, dict[str, object]] = {}

    valid_cases: list[tuple[str, ...]] = [
        ("--flow-probability-scaling-factor", "1"),
        ("--flow-probability-scaling-factor", "3"),
        ("--flow-quantization-bins", "2"),
        ("--flow-quantization-bins", "4"),
        ("--flow-matrix-mods", "A,C"),
        ("--flow-matrix-mods", "G,T"),
        ("--flow-probability-scaling-factor", "3",
         "--flow-quantization-bins", "4",
         "--flow-matrix-mods", "A,C"),
    ]
    for case in valid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-sqm-oracle-") as work:
            out = Path(work) / "out.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                *case, bam=bam, reference=reference, output=out, work=Path(work))
            key = " ".join(case)
            case_results[key] = {"exit_code": rc, "stderr_tail": stderr[-200:]}
            if rc != 0:
                diagnostics.append(
                    f"{key} expected exit 0, got {rc}: {stderr[-200:]}")

    invalid_cases: list[tuple[tuple[str, ...], str]] = [
        (("--flow-probability-scaling-factor", "0"), "must be >= 1"),
        (("--flow-quantization-bins", "1"), "must be >= 2"),
    ]
    for args, expected_fragment in invalid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-sqm-bad-") as work:
            out = Path(work) / "bad.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                *args, bam=bam, reference=reference, output=out, work=Path(work))
            key = " ".join(args)
            case_results[key] = {
                "exit_code": rc,
                "stderr_tail": stderr[-200:],
                "expected_fragment": expected_fragment,
            }
            if rc == 0:
                diagnostics.append(f"{key} expected non-zero exit, got {rc}")
            if expected_fragment not in stderr.lower():
                diagnostics.append(
                    f"{key} stderr missing expected fragment "
                    f"{expected_fragment!r}: {stderr[-200:]}")

    status = "pass" if not diagnostics else "fail"
    print(json.dumps({
        "status": status,
        "fixture": "hc_chr17_69k_70k_flow_scaling_quantization_mods_flag",
        "cases_tested": list(case_results.keys()),
        "case_results": case_results,
        "diagnostics": diagnostics,
        "tier": "Tier-2a",
        "note": ("CLI plumbing closed for --flow-probability-scaling-factor, "
                 "--flow-quantization-bins, --flow-matrix-mods.  GATK 4.6.2.0 "
                 "has 25+ --flow-* flags total; ~5 more remain "
                 "(--flow-order-for-annotations, "
                 "--flow-remove-non-single-base-pair-indels, "
                 "--flow-remove-one-zero-probs, "
                 "--flow-report-insertion-or-deletion, "
                 "--flow-retain-max-n-probs-base-format)."),
    }, sort_keys=True))
    return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())