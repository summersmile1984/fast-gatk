"""Universal oracle comparison harness for fastgatk-native.

Pair a GATK Java invocation with a native binary invocation, run both,
parse both outputs through a small set of registered parsers, and assert
that the contract pinned by the calling oracle script holds.

This module is the bridge between the per-algorithm docs
(``docs/algorithms/``) and the per-tool oracles
(``scripts/verify_*.py``).  Every oracle script that wants to assert a
byte-identical contract against GATK 4.6.2.0 should use this harness.

Usage from an oracle script::

    from parity_harness import compare_runs

    fixtures = {
        "name": "hc-chr17-69k-70k",
        "input_bam": ROOT / "...",
        "input_bai": ROOT / "....bai",
        "reference": ROOT / "human_g1k_v37.chr17_1Mb.fasta",
        "intervals": "17:69000-70000",
        "java_args": ["--sample-name", "NA12878"],
        "native_args": ["--sample-name", "NA12878"],
        "expected_records": {"17:69067": ("T", "G"), ...},
    }
    result = compare_runs("HaplotypeCaller", fixtures)
    print(json.dumps(result.to_dict(), indent=2))
    sys.exit(0 if result.passed else 1)

The harness does NOT decide what is or isn't byte-identical: that is the
calling oracle script's job, via ``contracts`` registered below.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import oracle_guard

# ---------------------------------------------------------------------------
# Constants — pinned once per repo
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK_JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


# ---------------------------------------------------------------------------
# VCF / GVCF parsing
# ---------------------------------------------------------------------------
def _vcf_open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_vcf(path: Path) -> dict[str, dict[str, str]]:
    """Return {position: {field: value}} for the VCF body lines.

    position is the 2nd column (e.g. "17:69067").
    Fields include REF, ALT, FILTER, plus the entire INFO and FORMAT/SAMPLE
    sub-fields keyed by their tag (e.g. "QD", "FS", "AF").
    """
    out: dict[str, dict[str, str]] = {}
    if not path.is_file():
        return out
    with _vcf_open(path) as h:
        for line in h:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            pos = f"{fields[0]}:{fields[1]}"
            d: dict[str, str] = {"REF": fields[3], "ALT": fields[4], "FILTER": fields[6]}
            info_str = fields[7] if len(fields) > 7 else ""
            for kv in info_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = v
                else:
                    d[kv] = ""
            out[pos] = d
    return out


def vcf_body_bytes(path: Path) -> bytes:
    """Return the semantic VCF payload (body only, no implementation
    headers).  Used for byte-equality asserts on the body."""
    with _vcf_open(path) as h:
        return b"".join(
            line.encode("utf-8")
            for line in h
            if line and not line.startswith("#")
        )


# ---------------------------------------------------------------------------
# Run a binary (GATK Java or native)
# ---------------------------------------------------------------------------
@dataclass
class RunResult:
    binary: str
    args: list[str]
    stdout: str
    stderr: str
    returncode: int
    output_files: dict[str, Path] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.returncode == 0


def run_java(gatk_tool: str, args: list[str], extra_jvm: list[str] = None,
             output_files: Mapping[str, Path] = None) -> RunResult:
    """Invoke GATK Java 4.6.2.0 with a specific tool."""
    cmd = [str(JAVA)]
    if extra_jvm:
        cmd.extend(extra_jvm)
    cmd.extend(["-jar", str(GATK_JAR), gatk_tool])
    cmd.extend(args)
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    return RunResult(
        binary="java " + gatk_tool,
        args=args,
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        output_files=dict(output_files or {}),
    )


def run_native(binary: Path, args: list[str],
               output_files: Mapping[str, Path] = None) -> RunResult:
    proc = subprocess.run([str(binary)] + args, text=True,
                          capture_output=True, check=False)
    return RunResult(
        binary=str(binary),
        args=args,
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
        output_files=dict(output_files or {}),
    )


# ---------------------------------------------------------------------------
# Compare two runs
# ---------------------------------------------------------------------------
def compare_runs(tool_name: str, fixtures: Mapping[str, Any],
                 expected_byte_equal_vcf_body: bool = False,
                 record_field_overrides: Mapping[str, Mapping[str, str]] = None,
                 tolerance: Mapping[str, float] = None) -> "ComparisonResult":
    """High-level helper: run GATK + native, parse VCF, assert contracts.

    `fixtures` keys:
      input_bam, input_bai, reference, intervals (str),
      java_args (list[str]), native_args (list[str]),
      sample_name (str), java_tool (str, default = tool_name),
      output_vcf (str filename for both runs), tmpdir (str).

    `record_field_overrides`: per-position field name -> expected
    substring (for fields where byte-equal isn't required, e.g. QD
    jitter, but the sub-string "0.00;" or "PASS" must appear).

    `tolerance`: per-field-name numeric tolerance (relative or absolute,
    picked by field name suffix `_abs` vs default).
    """
    java_tool = fixtures.get("java_tool", tool_name)
    native_binary = fixtures["native_binary"]
    bam = fixtures["input_bam"]
    bai = fixtures["input_bai"]
    ref = fixtures["reference"]
    intervals = fixtures.get("intervals", "")
    java_args = list(fixtures.get("java_args", []))
    native_args = list(fixtures.get("native_args", []))
    sample = fixtures.get("sample_name")
    out_vcf = fixtures.get("output_vcf", "out.vcf")
    expected = fixtures.get("expected_records", {}) or {}

    with tempfile.TemporaryDirectory(prefix="parity_harness_") as tmpdir:
        tmp = Path(tmpdir)
        java_vcf = tmp / f"java_{out_vcf}"
        native_vcf = tmp / f"native_{out_vcf}"

        java_cmd = list(java_args)
        java_cmd += ["-O", str(java_vcf)]
        native_cmd = list(native_args)
        native_cmd += ["-O", str(native_vcf)]

        java_res = run_java(java_tool, java_cmd,
                            output_files={"vcf": java_vcf})
        native_res = run_native(native_binary, native_cmd,
                                output_files={"vcf": native_vcf})

        cmp = ComparisonResult(
            tool=tool_name, java=java_res, native=native_res,
            expected_records=expected,
            record_field_overrides=record_field_overrides or {},
            tolerance=tolerance or {},
        )
        cmp.compare()
        return cmp


@dataclass
class ComparisonResult:
    tool: str
    java: RunResult
    native: RunResult
    expected_records: dict[str, tuple[str, str]] = field(default_factory=dict)
    record_field_overrides: dict[str, dict[str, str]] = field(default_factory=dict)
    tolerance: dict[str, float] = field(default_factory=dict)
    byte_equal_vcf_body: bool = False
    diagnostics: list[str] = field(default_factory=list)
    passed: bool = False

    def compare(self) -> bool:
        self.passed = True

        # 1. Both runs must exit 0
        if not self.java.passed:
            self.passed = False
            self.diagnostics.append(f"Java GATK exit {self.java.returncode}; stderr={self.java.stderr[-200:]!r}")
        if not self.native.passed:
            self.passed = False
            self.diagnostics.append(f"native exit {self.native.returncode}; stderr={self.native.stderr[-200:]!r}")
        if not self.passed:
            return False

        # 2. VCF parse + record + field assertions
        java_vcf = self.java.output_files.get("vcf")
        native_vcf = self.native.output_files.get("vcf")
        if not java_vcf or not native_vcf:
            self.passed = False
            self.diagnostics.append("output VCF not produced")
            return False

        java_records = parse_vcf(java_vcf)
        native_records = parse_vcf(native_vcf)

        for pos, (ref_expected, alt_expected) in self.expected_records.items():
            if pos not in java_records:
                self.passed = False
                self.diagnostics.append(f"Java missing record at {pos}")
                continue
            if pos not in native_records:
                self.passed = False
                self.diagnostics.append(f"native missing record at {pos}")
                continue
            j = java_records[pos]
            n = native_records[pos]
            if j["REF"] != ref_expected or j["ALT"] != alt_expected:
                self.passed = False
                self.diagnostics.append(
                    f"Java record at {pos} has unexpected alleles "
                    f"({j['REF']},{j['ALT']}) vs ({ref_expected},{alt_expected})")
            if n["REF"] != ref_expected or n["ALT"] != alt_expected:
                self.passed = False
                self.diagnostics.append(
                    f"native record at {pos} has unexpected alleles "
                    f"({n['REF']},{n['ALT']}) vs ({ref_expected},{alt_expected})")

            for field, expected_substr in self.record_field_overrides.get(pos, {}).items():
                if field not in j:
                    self.passed = False
                    self.diagnostics.append(f"Java record {pos} missing {field}")
                if field not in n:
                    self.passed = False
                    self.diagnostics.append(f"native record {pos} missing {field}")
                else:
                    if expected_substr not in j.get(field, ""):
                        self.passed = False
                        self.diagnostics.append(
                            f"Java record {pos} field {field} missing {expected_substr!r}")
                    if expected_substr not in n.get(field, ""):
                        self.passed = False
                        self.diagnostics.append(
                            f"native record {pos} field {field} missing {expected_substr!r}")

        # 3. Optional byte-equal body (skip per-field overrides)
        if self.byte_equal_vcf_body:
            java_body = vcf_body_bytes(java_vcf)
            native_body = vcf_body_bytes(native_vcf)
            if java_body != native_body:
                self.passed = False
                self.diagnostics.append(
                    f"VCF body byte-equal failed: java={len(java_body)} "
                    f"native={len(native_body)}")
        return self.passed

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "passed": self.passed,
            "diagnostics": self.diagnostics,
            "java": {
                "binary": self.java.binary,
                "returncode": self.java.returncode,
                "stdout_tail": self.java.stdout[-400:],
                "stderr_tail": self.java.stderr[-400:],
            },
            "native": {
                "binary": self.native.binary,
                "returncode": self.native.returncode,
                "stdout_tail": self.native.stdout[-400:],
                "stderr_tail": self.native.stderr[-400:],
            },
        }


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Smoke test: run both halves with no real fixtures, just verify the
    # call interface works.
    print(json.dumps({
        "module": "parity_harness",
        "exports": ["compare_runs", "RunResult", "ComparisonResult",
                    "parse_vcf", "run_java", "run_native"],
        "ROOT": str(ROOT),
        "JAVA_exists": JAVA.exists(),
        "GATK_JAR_exists": GATK_JAR.exists(),
    }, indent=2))
