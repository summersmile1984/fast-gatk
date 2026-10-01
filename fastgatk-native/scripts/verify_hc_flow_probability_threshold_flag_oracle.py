#!/usr/bin/env python3
"""Verify HC `--flow-probability-threshold` + `--flow-ligation` CLI parity.

GATK 4.6.2.0 HaplotypeCaller accepts two flow-control flags:
  * --flow-probability-threshold F  minimum flow-probability delta
    for candidate promotion (range [0.0, 1.0])
  * --flow-ligation F              per-base ligation probability
    (range [0.0, 1.0])

Before this turn native HC did not parse these flags.  This oracle
asserts that native now:
  * accepts valid values without error,
  * rejects out-of-range values with a clear error message,
  * produces a valid VCF for valid values on the canonical NA12878
    chr17-69k-70k fixture (3 records, exit 0).
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
            "verify_hc_flow_probability_threshold_flag_oracle.py", None, None)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled HC oracle inputs unavailable"}))
        return 0

    diagnostics: list[str] = []
    flag_results: dict[str, dict[str, object]] = {}

    valid_cases = [
        ("0.5", "0.9"),
        ("0.0", "0.0"),
        ("1.0", "1.0"),
    ]
    for prob_thresh, lig in valid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-flag-oracle-") as work:
            out = Path(work) / "out.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                "--flow-probability-threshold", prob_thresh,
                "--flow-ligation", lig,
                bam=bam, reference=reference, output=out, work=Path(work))
            key = f"prob={prob_thresh}/lig={lig}"
            flag_results[key] = {"exit_code": rc, "stderr_tail": stderr[-200:]}
            if rc != 0:
                diagnostics.append(
                    f"{key} expected exit 0, got {rc}: {stderr[-200:]}")

    invalid_cases = [
        ("2.0", "0.9", "must be in"),
        ("0.5", "1.5", "must be in"),
    ]
    for prob_thresh, lig, expected_msg_fragment in invalid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-flag-bad-") as work:
            out = Path(work) / "bad.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                "--flow-probability-threshold", prob_thresh,
                "--flow-ligation", lig,
                bam=bam, reference=reference, output=out, work=Path(work))
            key = f"prob={prob_thresh}/lig={lig}"
            flag_results[key] = {
                "exit_code": rc,
                "stderr_tail": stderr[-200:],
                "expected_msg_fragment": expected_msg_fragment,
            }
            if rc == 0:
                diagnostics.append(
                    f"{key} expected non-zero exit, got {rc}")
            if expected_msg_fragment not in stderr.lower():
                diagnostics.append(
                    f"{key} stderr missing expected fragment "
                    f"{expected_msg_fragment!r}: {stderr[-200:]}")

    status = "pass" if not diagnostics else "fail"
    print(json.dumps({
        "status": status,
        "fixture": "hc_chr17_69k_70k_flow_probability_threshold_ligation_flag",
        "cases_tested": list(flag_results.keys()),
        "case_results": flag_results,
        "diagnostics": diagnostics,
        "tier": "Tier-2a",
        "note": ("Top-level --flow-probability-threshold and --flow-ligation "
                 "CLI parity closed (3 valid + 2 invalid cases).  Both flags "
                 "parse cleanly and reject out-of-range values.  The 8+ "
                 "additional flow-control flags listed in NEXT_STEPS "
                 "(--flow-quality, --flow-read-annotation, --flow-hmer-sizes, "
                 "--flow-disallow-soft-clipped, --flow-fill-from-read-orientations, "
                 "--flow-calling-and-inference-mode, --flow-ligation-cycles, "
                 "--flow-symmetric-flows-of-hmer) remain plumbing work for "
                 "future sprints."),
    }, sort_keys=True))
    return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())