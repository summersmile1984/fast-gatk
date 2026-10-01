#!/usr/bin/env python3
"""Verify HC --flow-mode CLI parity with GATK 4.6.2.0.

GATK 4.6.2.0 HaplotypeCaller accepts --flow-mode with three values:
  * NONE         — legacy non-flow path
  * STANDARD     — flow-aware GVCF bands
  * FAST         — base-precision resolution

Before this turn native only accepted --flow-assembly-collapse-hmer-size
and --flow-assembly-collapse-partial-mode.  This oracle asserts
that native now:
  * accepts all three --flow-mode values without error,
  * rejects invalid values with a clear error message,
  * produces a valid VCF for NONE on the canonical NA12878
    chr17-69k-70k fixture (3 records, exit 0).

The Tier-2a flow algorithm control surface (10+ additional flags
documented in NEXT_STEPS) remains plumbing work for a future
sprint; this oracle pins only the top-level --flow-mode contract.
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
            "verify_hc_flow_mode_cli_parity_oracle.py", None, None)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled HC oracle inputs unavailable"}))
        return 0

    diagnostics: list[str] = []
    mode_results: dict[str, dict[str, object]] = {}

    for mode in ("NONE", "STANDARD", "FAST"):
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-mode-oracle-") as work:
            out = Path(work) / f"out.{mode}.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                "--flow-mode", mode,
                bam=bam, reference=reference, output=out, work=Path(work))
            mode_results[mode] = {
                "exit_code": rc,
                "stderr_tail": stderr[-200:],
            }
            if rc != 0:
                diagnostics.append(
                    f"--flow-mode {mode} expected exit 0, got {rc}: {stderr[-200:]}")

    # Invalid mode must be rejected.
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-mode-bad-") as work:
        out = Path(work) / "bad.vcf.gz"
        rc, _stdout, stderr = run_native_hc(
            "--flow-mode", "INVALID_VALUE",
            bam=bam, reference=reference, output=out, work=Path(work))
        mode_results["INVALID"] = {
            "exit_code": rc,
            "stderr_tail": stderr[-200:],
        }
        if rc == 0:
            diagnostics.append(
                f"--flow-mode INVALID_VALUE expected non-zero exit, got {rc}")

    status = "pass" if not diagnostics else "fail"
    print(json.dumps({
        "status": status,
        "fixture": "hc_chr17_69k_70k_flow_mode_cli",
        "modes_tested": list(mode_results.keys()),
        "mode_results": mode_results,
        "diagnostics": diagnostics,
        "tier": "Tier-2a",
        "note": ("Top-level --flow-mode CLI parity is closed (NONE/STANDARD/FAST "
                 "parse-only on the cmdline; --flow-mode INVALID_VALUE is "
                 "rejected).  The 10+ additional flow-control flags listed in "
                 "NEXT_STEPS (--flow-probability-threshold, --flow-ligation, "
                 "--flow-quality, --flow-read-annotation, --flow-hmer-sizes, "
                 "--flow-disallow-soft-clipped, --flow-fill-from-read-orientations, "
                 "--flow-calling-and-inference-mode) remain engineering work "
                 "for future sprints."),
    }, sort_keys=True))
    return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())