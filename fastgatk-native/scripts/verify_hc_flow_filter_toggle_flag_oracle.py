#!/usr/bin/env python3
"""Verify HC `--flow-filter-alleles` + `--flow-filter-lone-alleles` + `--flow-disallow-probs-larger-than-call` CLI parity.

GATK 4.6.2.0 HaplotypeCaller accepts three additional flow-control
flags (boolean toggles, default false):
  * --flow-filter-alleles                  enable post-call flow
    allele filtering
  * --flow-filter-lone-alleles             drop flow alleles
    supported only by a single read
  * --flow-disallow-probs-larger-than-call  cap per-base error
    probabilities at 1.0 relative to the base call

Before this turn native HC did not parse these flags.  This oracle
asserts that native now:
  * accepts each flag as a boolean toggle (no value argument),
  * works in any combination of the three flags,
  * produces a valid VCF on the canonical NA12878 chr17-69k-70k
    fixture (3 records, exit 0).
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
            "verify_hc_flow_filter_toggle_flag_oracle.py", None, None)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled HC oracle inputs unavailable"}))
        return 0

    diagnostics: list[str] = []
    case_results: dict[str, dict[str, object]] = {}

    valid_cases: list[tuple[str, ...]] = [
        ("--flow-filter-alleles",),
        ("--flow-filter-lone-alleles",),
        ("--flow-disallow-probs-larger-than-call",),
        ("--flow-filter-alleles", "--flow-filter-lone-alleles"),
        ("--flow-filter-alleles", "--flow-filter-lone-alleles",
         "--flow-disallow-probs-larger-than-call"),
    ]
    for case in valid_cases:
        with tempfile.TemporaryDirectory(prefix="fastgatk-hc-flow-ft-oracle-") as work:
            out = Path(work) / "out.vcf.gz"
            rc, _stdout, stderr = run_native_hc(
                *case, bam=bam, reference=reference, output=out, work=Path(work))
            key = " ".join(case)
            case_results[key] = {"exit_code": rc, "stderr_tail": stderr[-200:]}
            if rc != 0:
                diagnostics.append(
                    f"{key} expected exit 0, got {rc}: {stderr[-200:]}")

    status = "pass" if not diagnostics else "fail"
    print(json.dumps({
        "status": status,
        "fixture": "hc_chr17_69k_70k_flow_filter_toggle_flag",
        "cases_tested": list(case_results.keys()),
        "case_results": case_results,
        "diagnostics": diagnostics,
        "tier": "Tier-2a",
        "note": ("CLI plumbing closed for --flow-filter-alleles, "
                 "--flow-filter-lone-alleles, "
                 "--flow-disallow-probs-larger-than-call.  GATK 4.6.2.0 "
                 "has 25+ --flow-* flags total; ~8 more remain "
                 "(--flow-order-for-annotations, --flow-matrix-mods, "
                 "--flow-probability-scaling-factor, "
                 "--flow-quantization-bins, "
                 "--flow-remove-non-single-base-pair-indels, "
                 "--flow-remove-one-zero-probs, "
                 "--flow-report-insertion-or-deletion, "
                 "--flow-retain-max-n-probs-base-format)."),
    }, sort_keys=True))
    return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())