#!/usr/bin/env python3
"""P1 CLI surface: GATK-public options enter native; unknown/cloud fail closed.

Java GATK is an oracle for --help (exit 0) and unknown-option rejection only.
The dispatcher never launches gatk.jar.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "fastgatk-native/scripts"))
from oracle_guard import oracle_ready  # noqa: E402

DISPATCHER = ROOT / "fastgatk-native/dispatcher/fastgatk"
REGISTRY = ROOT / "fastgatk-native/dispatcher/tool_registry.json"
COVERAGE = ROOT / "fastgatk-native/evidence/cli-coverage-48.json"
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
FASTA = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run_dispatcher(*args: str, expect: int = 0,
                   env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    result = subprocess.run([str(DISPATCHER), *args], text=True, capture_output=True, env=merged)
    assert result.returncode == expect, (args, result.returncode, result.stdout, result.stderr)
    return result


def required_args(entry: dict) -> list[str]:
    args: list[str] = []
    for required in entry.get("required_inputs", []):
        alts = [p.strip() for p in str(required).split("|") if p.strip()]
        option = alts[0]
        if option in {"-I", "--input"}:
            args += [option, str(BAM)]
        elif option in {"-V", "--variant"}:
            args += [option, str(BAM)]  # path presence only; dispatcher dry-run
        elif option in {"-R", "--reference"}:
            args += [option, str(FASTA)]
        else:
            args += [option, "dummy"]
    for required in entry.get("required_outputs", []):
        args += [str(required).split("|")[0].strip(), "/tmp/fastgatk-cli-parity.out"]
    return args


def main() -> int:
    assert DISPATCHER.is_file()
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    coverage = json.loads(COVERAGE.read_text(encoding="utf-8"))
    assert registry.get("runtime_java") is False
    assert coverage["tool_count"] == len(registry["tools"])
    assert coverage["summary"].get("must-implement", 0) == 0
    tools = registry["tools"]
    checked = 0
    for name, entry in tools.items():
        args = required_args(entry)
        extra: list[str] = []
        flags = {str(x) for x in entry.get("flag_options", [])}
        values = {str(x) for x in entry.get("value_options", [])}
        if "--QUIET" in flags or "--QUIET" in values:
            extra.append("--QUIET")
        elif "--quiet" in flags or "--quiet" in values:
            extra.append("--quiet")
        if "--verbosity" in values:
            extra += ["--verbosity", "INFO"]
        elif "-verbosity" in values:
            extra += ["-verbosity", "INFO"]
        dry = run_dispatcher("--dry-run", name, *args, *extra)
        plan = json.loads(dry.stdout.splitlines()[-1])
        assert plan["execution_mode"] in {"native", "adapter"}, (name, plan)
        assert "fallback" not in json.dumps(plan).lower()
        if entry.get("status") != "adapter":
            unknown = run_dispatcher(name, *args, "--not-a-gatk-p1-option", expect=2)
            assert unknown.stderr.strip(), (name, unknown.stdout, unknown.stderr)
            error = json.loads(unknown.stderr.splitlines()[-1])
            assert error["category"] == "UNSUPPORTED_PARAMETER", (name, error)
            cloud = run_dispatcher(name, *args, "--cloud-prefetch-buffer", "40", expect=2)
            assert cloud.stderr.strip(), (name, cloud.stdout, cloud.stderr)
            cloud_error = json.loads(cloud.stderr.splitlines()[-1])
            assert cloud_error["category"] == "UNSUPPORTED_PARAMETER", (name, cloud_error)
        checked += 1
    # Oracle: GATK help succeeds; unknown option is rejected. Native --help is 0.
    if oracle_ready("verify_cli_parity_gatk_oracle.py", JAVA, GATK):
        help_java = subprocess.run(
            [str(JAVA), "-jar", str(GATK), "HaplotypeCaller", "--help"],
            text=True, capture_output=True,
        )
        assert help_java.returncode == 0, help_java.stderr[-500:]
        bad_java = subprocess.run(
            [str(JAVA), "-jar", str(GATK), "HaplotypeCaller", "--not-a-gatk-p1-option"],
            text=True, capture_output=True,
        )
        assert bad_java.returncode != 0
        native_help = run_dispatcher("HaplotypeCaller", "--help")
        assert native_help.returncode == 0
    count_reads = ROOT / "fastgatk-native/build/fastgatk-count-reads"
    native_binary = os.environ.get("FASTGATK_COUNT_READS_BINARY", str(count_reads))
    if Path(native_binary).is_file() and BAM.is_file():
        native = subprocess.run(
            [native_binary, "-I", str(BAM), "--QUIET", "--verbosity", "INFO", "--tmp-dir", "/tmp"],
            text=True, capture_output=True,
        )
        assert native.returncode == 0, native.stderr[-800:]
        native_bad = subprocess.run(
            [native_binary, "-I", str(BAM), "--not-a-gatk-p1-option"],
            text=True, capture_output=True,
        )
        assert native_bad.returncode != 0
        assert "unknown option" in native_bad.stderr
    print(json.dumps({
        "status": "pass",
        "tools": checked,
        "must-implement": 0,
        "runtime_java": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
