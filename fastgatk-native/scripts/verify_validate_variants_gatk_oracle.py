#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for ValidateVariants argument/validation gates.

The existing contract covers the native validation counters and GVCF coverage
implementation.  This oracle locks the command-line boundary that matters to
the GATK+SLURM+Nextflow replacement path: Barclay Boolean spellings, warning
versus fail-closed behavior, and the two documented incompatible option pairs.
It deliberately compares outcomes, not Java log wording or timestamps.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


VALID_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
chr1\t2\t.\tC\tT\t50\tLowQual\t.\tGT\t0/1
"""

INVALID_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/1
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
        oracle_guard.oracle_not_verified('verify_validate_variants_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ValidateVariants Java oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-validate-variants-gatk-oracle-") as directory:
        work = Path(directory)
        valid = work / "valid.vcf"
        invalid = work / "invalid.vcf"
        dbsnp = work / "dbsnp.vcf"
        valid.write_text(VALID_VCF, encoding="utf-8")
        invalid.write_text(INVALID_VCF, encoding="utf-8")
        dbsnp.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\trs1\tA\tG\t50\tPASS\t.\n",
            encoding="utf-8",
        )

        def gatk_run(path: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
            return run([str(java), "-Xmx1g", "-jar", str(gatk), "ValidateVariants",
                        "-V", str(path), *args])

        def native_run(path: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
            return run([str(binary), "-V", str(path), *args])

        # Explicit Boolean values are what Barclay emits in many generated
        # Nextflow command lines.  Native must accept the same forms even when
        # no reference is supplied and all record validations are excluded.
        forms = {
            "validate_gvcf_false": ["--validate-GVCF", "false",
                                    "--validation-type-to-exclude", "ALL"],
            "skip_filtered_true": ["--do-not-validate-filtered-records", "true",
                                    "--validation-type-to-exclude", "ALL"],
            "fail_overlap_false": ["--fail-gvcf-on-overlap", "false",
                                    "--validation-type-to-exclude", "ALL"],
            "dictionary_validation_true": ["--disable-sequence-dictionary-validation", "true",
                                            "--validation-type-to-exclude", "ALL"],
        }
        form_results: dict[str, tuple[int, int]] = {}
        for label, args in forms.items():
            gatk_result = gatk_run(valid, args)
            native_result = native_run(valid, args)
            if gatk_result.returncode != 0 or native_result.returncode != 0:
                raise AssertionError(
                    f"Boolean form {label} diverged: GATK={gatk_result.returncode} "
                    f"native={native_result.returncode}; native stderr={native_result.stderr[-1000:]}")
            form_results[label] = (gatk_result.returncode, native_result.returncode)

        # ValidateVariants exposes several single-dash Barclay synonyms.  A
        # launcher that forwards GATK command lines must accept these exact
        # spellings as well as the long names above.  Index the dbSNP fixture
        # because GATK requires random access for FeatureInput tracks; native
        # intentionally accepts the same small VCF sequentially.
        index_result = run([
            str(java), "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile",
            "-I", str(dbsnp),
        ])
        if index_result.returncode != 0:
            raise AssertionError(f"failed to index dbSNP fixture: {index_result.stderr[-1000:]}")
        short_forms = {
            "dbsnp_short": ["-D", str(dbsnp), "--validation-type-to-exclude", "ALL"],
            "interval_set_rule_short": ["-isr", "UNION", "--validation-type-to-exclude", "ALL"],
            "skip_filtered_short": ["-do-not-validate-filtered-records", "true",
                                     "--validation-type-to-exclude", "ALL"],
            "warn_short": ["-warn-on-errors", "true", "--validation-type-to-exclude", "ALL"],
            "dictionary_validation_short": ["-disable-sequence-dictionary-validation", "true",
                                             "--validation-type-to-exclude", "ALL"],
        }
        for label, args in short_forms.items():
            gatk_result = gatk_run(valid, args)
            native_result = native_run(valid, args)
            if gatk_result.returncode != 0 or native_result.returncode != 0:
                raise AssertionError(
                    f"short alias {label} diverged: GATK={gatk_result.returncode} "
                    f"native={native_result.returncode}; native stderr={native_result.stderr[-1000:]}")
            form_results[label] = (gatk_result.returncode, native_result.returncode)

        # GATK's warning switch converts a strict validation exception into a
        # successful run.  The same malformed unused ALT must fail without
        # the switch and succeed with the explicit true form.
        strict_args = ["--validation-type-to-exclude", "REF"]
        for args in (strict_args,
                     ["--validation-type-to-exclude", "REF", "--warn-on-errors", "false"],):
            gatk_result = gatk_run(invalid, args)
            native_result = native_run(invalid, args)
            if gatk_result.returncode == 0 or native_result.returncode == 0:
                raise AssertionError("strict invalid ALT did not fail in both implementations")
        warning_results = []
        for warning_flag in (["--warn-on-errors", "true"],):
            args = ["--validation-type-to-exclude", "REF", *warning_flag]
            gatk_result = gatk_run(invalid, args)
            native_result = native_run(invalid, args)
            if gatk_result.returncode != 0 or native_result.returncode != 0:
                raise AssertionError(
                    f"warn-on-errors diverged: GATK={gatk_result.returncode} "
                    f"native={native_result.returncode}; native stderr={native_result.stderr[-1000:]}")
            warning_results.append(warning_flag[0])

        # These are explicit Barclay argument constraints.  Rejecting before
        # opening/decoding input is important for deterministic dispatcher
        # fallback and avoids silently changing GVCF semantics.
        incompatible = [
            ["--validate-GVCF", "true", "--do-not-validate-filtered-records", "true"],
            ["--fail-gvcf-on-overlap", "true", "--do-not-validate-filtered-records", "true"],
        ]
        for args in incompatible:
            gatk_result = gatk_run(valid, args)
            native_result = native_run(valid, args)
            if gatk_result.returncode == 0 or native_result.returncode == 0:
                raise AssertionError(
                    f"incompatible options were accepted: {args}; "
                    f"GATK={gatk_result.returncode} native={native_result.returncode}")

        # Invalid Boolean literals must fail closed at the owning option in
        # both implementations, without being interpreted as a bare switch.
        invalid_boolean = ["--validate-GVCF", "maybe"]
        gatk_result = gatk_run(valid, invalid_boolean)
        native_result = native_run(valid, invalid_boolean)
        if gatk_result.returncode == 0 or native_result.returncode == 0:
            raise AssertionError("invalid Boolean literal was accepted")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "optional_boolean_forms_exact": sorted(form_results),
            "warning_forms_exact": sorted(warning_results),
            "strict_fail_closed": True,
            "incompatible_gvcf_switches_fail_closed": True,
            "invalid_boolean_fail_closed": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
