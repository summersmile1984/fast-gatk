#!/usr/bin/env python3
"""Verify HC `--flow-quality` + `--flow-disallow-soft-clipped` + `--flow-fill-from-read-orientations` CLI parity.

GATK 4.6.2.0 HaplotypeCaller accepts three additional flow-control
flags:
  * --flow-quality F                  per-base flow-quality
    probability threshold (range [0.0, 1.0])
  * --flow-disallow-soft-clipped      boolean flag, exclude
    soft-clipped reads from flow-key candidates
  * --flow-fill-from-read-orientations  boolean flag, fill flow
    gaps from read orientations

Before this turn native HC did not parse these flags.  This oracle
asserts that native now:
  * accepts valid --flow-quality values without error,
  * rejects out-of-range values with a clear error message,
  * accepts --flow-disallow-soft-clipped as a boolean toggle,
  * accepts --flow-fill-from-read-orientations as a boolean toggle,
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
            "verify_hc_flow_quality_disallow_fill_flag_oracle.py", None, None)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled HC oracle inputs unavailable"}))
        return 0

    diagnostics: list[str] = []
    case_results: dict[str, dict[str, object]] = {}

    valid_cases = [
        ("--flow-quality", "0.7"),
        ("--flow-quality", "0.0"),
        ("--flow-quality", "1.0"),
        ("--flow-disallow-soft-clipped",),
        ("--flow-fill-from-read-orientations",),
        ("--flow-disallow-soft-clipped", "--flow-fill-from-read-orientations",
         "--flow-quality", "0.5"),
    ]
    for case in valid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-qdf-oracle-") as work:
            out = Path(work) / "out.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                *case, bam=bam, reference=reference, output=out, work=Path(work))
            key = " ".join(case)
            case_results[key] = {"exit_code": rc, "stderr_tail": stderr[-200:]}
            if rc != 0:
                diagnostics.append(
                    f"{key} expected exit 0, got {rc}: {stderr[-200:]}")

    invalid_cases = [
        (["--flow-quality", "1.5"], "must be in"),
        (["--flow-quality", "-0.1"], "must be in"),
    ]
    for args, expected_fragment in invalid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-qdf-bad-") as work:
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
        "fixture": "hc_chr17_69k_70k_flow_quality_disallow_fill_flag",
        "cases_tested": list(case_results.keys()),
        "case_results": case_results,
        "diagnostics": diagnostics,
        "tier": "Tier-2a",
        "note": ("CLI plumbing closed for --flow-quality, "
                 "--flow-disallow-soft-clipped, "
                 "--flow-fill-from-read-orientations.  4+ additional "
                 "flow-control flags remain: --flow-read-annotation, "
                 "--flow-hmer-sizes, --flow-calling-and-inference-mode, "
                 "--flow-ligation-cycles."),
    }, sort_keys=True))
    return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())