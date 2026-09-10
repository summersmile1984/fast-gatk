#!/usr/bin/env python3
"""Verify GATK optional-boolean parsing before any input/output is touched."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


# Every entry uses the same three invalid spellings, but the option names cover
# all native tools migrated to the shared parser.  The command line is ordered
# so the invalid value is encountered before the required input/output checks.
TOOLS = {
    "fastgatk-hc-call": ["--native-pair-hmm-use-double-precision"],
    "fastgatk-bqsr": ["--use-original-qualities"],
    "fastgatk-apply-bqsr": ["--use-original-qualities"],
    "fastgatk-mutect2": ["--mitochondria-mode"],
    "fastgatk-filter-mutect-calls": ["--mitochondria-mode"],
    "fastgatk-apply-vqsr": ["--AS"],
    "fastgatk-variant-recalibrator": ["--AS", "--create-output-variant-index"],
    "fastgatk-genotype-gvcf": [
        "--create-output-variant-index", "--use-new-qual-calculator", "--new-qual",
        "--use-posteriors-to-calculate-qual", "--gp-qual"],
    "fastgatk-reblock-gvcf": ["--create-output-variant-index"],
    "fastgatk-left-align-trim": ["--create-output-variant-index"],
    "fastgatk-combine-gvcfs": ["--create-output-variant-index"],
    "fastgatk-variant-filtration": ["--create-output-variant-index"],
    "fastgatk-select-variants": [
        "--create-output-variant-index", "--exclude-non-variants",
        "--exclude-filtered", "--exclude-filtered-variants",
        "--set-filtered-gt-to-nocall",
    ],
    "fastgatk-sort-sam": ["--add-pg-tag", "--create-output-bam-index"],
    "fastgatk-mark-duplicates": ["--add-pg-tag", "--create-output-bam-index"],
    "fastgatk-gather-vcfs": [
        "--reorder-input-by-first-variant", "--REORDER_INPUT_BY_FIRST_VARIANT",
        "--create-output-variant-index", "--CREATE_INDEX"],
    "fastgatk-depth-of-coverage": [
        "--omit-locus-table", "--omit-depth-output-at-each-base",
        "--omit-interval-statistics", "--omit-per-sample-statistics",
        "--omit-intervals", "--omit-sample-summary"],
}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    binaries = {name: build / name for name in TOOLS}
    for name, binary in binaries.items():
        assert binary.is_file() and os.access(binary, os.X_OK), binary

    with tempfile.TemporaryDirectory(prefix="fastgatk-optional-boolean-") as temp:
        work = Path(temp)
        failures: dict[str, dict[str, object]] = {}
        for name, binary in binaries.items():
            for option in TOOLS[name]:
                option_key = option.lstrip("-").replace("-", "_")
                for suffix, option_tokens in (
                    ("inline", [f"{option}=maybe"]),
                    ("empty-inline", [f"{option}="]),
                    ("next-token", [option, "maybe"]),
                ):
                    output = work / f"{name}-{option_key}-{suffix}.vcf.gz"
                # The parser must reject the literal before opening the missing
                # VCF, creating output, or staging any dispatcher sidecar.
                    result = subprocess.run(
                        [str(binary), *option_tokens, "-V", str(work / "missing.vcf.gz"),
                         "-O", str(output)],
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    key = f"{name}:{option}:{suffix}"
                    assert result.returncode != 0, (key, result.stdout, result.stderr)
                    assert "invalid boolean" in result.stderr, (key, result.stderr)
                    assert not output.exists(), (key, output)
                    failures[key] = {
                        "returncode": result.returncode,
                        "stderr": result.stderr.strip(),
                    }

    print(json.dumps({"status": "pass", "tools": failures}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
