#!/usr/bin/env python3
"""Pinned GATK oracle for symbolic-ALT allele-usage validation.

GATK's ``VariantContext.validateAlternateAlleles`` deliberately excludes
symbolic alternates from the "observed in a called genotype" check.  This is
important for valid ``<DEL>``/``<CNV>`` records that are 0/0 or no-call.  The
native validator used to treat those strings as ordinary concrete alleles and
rejected them.  Keep this oracle focused on return-code parity; Java exception
wording is release-specific and not part of the replacement contract.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=DEL,Description=deletion>
##ALT=<ID=CNV,Description=copy-number event>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VALIDATE_VARIANTS_BINARY",
        root / "fastgatk-native/build/fastgatk-validate-variants",
    ))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, gatk)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_validate_variants_symbolic_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ValidateVariants symbolic-ALT oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    cases = {
        # A symbolic alternate is not a concrete allele and need not be called.
        "unobserved_symbolic_del": "chr1\t10\t.\tA\t<DEL>\t50\tPASS\t.\tGT\t0/0\n",
        "unobserved_symbolic_cnv_nocall": "chr1\t20\t.\tA\t<CNV>\t50\tPASS\t.\tGT\t./.\n",
        # A concrete alternate still must be present even when a symbolic ALT
        # occurs beside it; this prevents the fix from weakening the check.
        "symbolic_plus_observed_concrete": "chr1\t30\t.\tA\t<DEL>,C\t50\tPASS\t.\tGT\t0/2\n",
        "unobserved_concrete_still_fails": "chr1\t40\t.\tA\t<DEL>,C\t50\tPASS\t.\tGT\t0/0\n",
    }

    with tempfile.TemporaryDirectory(prefix="fastgatk-validate-variants-symbolic-oracle-") as directory:
        work = Path(directory)
        results: dict[str, dict[str, int]] = {}
        for name, row in cases.items():
            path = work / f"{name}.vcf"
            path.write_text(HEADER + row, encoding="utf-8")
            gatk_result = run([
                str(java), "-Xmx1g", "-jar", str(gatk), "ValidateVariants",
                "-V", str(path),
            ])
            native_result = run([str(binary), "-V", str(path)])
            java_ok = int(gatk_result.returncode == 0)
            native_ok = int(native_result.returncode == 0)
            if java_ok != native_ok:
                raise AssertionError(
                    f"symbolic ALT case {name} diverged: "
                    f"GATK={gatk_result.returncode} native={native_result.returncode}; "
                    f"native stderr={native_result.stderr[-1000:]}")
            if name != "unobserved_concrete_still_fails" and not java_ok:
                raise AssertionError(f"GATK unexpectedly rejected valid symbolic case {name}")
            if name == "unobserved_concrete_still_fails" and java_ok:
                raise AssertionError("GATK unexpectedly accepted an unobserved concrete ALT")
            results[name] = {"gatk_returncode": gatk_result.returncode,
                             "native_returncode": native_result.returncode}

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "symbolic_alts_excluded_from_usage_check": True,
        "concrete_alt_check_preserved": True,
        "cases": results,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
