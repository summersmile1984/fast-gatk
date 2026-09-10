#!/usr/bin/env python3
"""Verify the GATK launcher compatibility boundary.

This gate deliberately uses a capture-only native executable.  It proves that
launcher controls are consumed by the dispatcher, while tool arguments remain
ordered and unknown native parameters cannot accidentally reach a binary.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DISPATCHER = ROOT / "fastgatk-native/dispatcher/fastgatk"
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"


def run(*args: str, env: dict[str, str] | None = None,
        expect: int = 0) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    result = subprocess.run([str(DISPATCHER), *args], text=True,
                            capture_output=True, env=merged, check=False)
    assert result.returncode == expect, (args, result.returncode,
                                         result.stdout, result.stderr)
    return result


def json_stdout(result: subprocess.CompletedProcess[str]) -> dict[str, object]:
    return json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    assert DISPATCHER.is_file() and os.access(DISPATCHER, os.X_OK)
    assert BAM.is_file()
    with tempfile.TemporaryDirectory(prefix="fastgatk-launcher-contract-") as directory:
        work = Path(directory)
        capture = work / "native-argv.json"
        fake_native = work / "capture-native"
        fake_native.write_text(
            "#!" + sys.executable + "\n"
            "import json, os, sys\n"
            "from pathlib import Path\n"
            "Path(os.environ['FASTGATK_CAPTURE']).write_text("
            "json.dumps(sys.argv[1:]), encoding='utf-8')\n",
            encoding="utf-8",
        )
        fake_native.chmod(fake_native.stat().st_mode | stat.S_IXUSR)

        # Both launcher spellings are accepted before and after the tool.  The
        # options are recorded at the dispatcher boundary and never forwarded
        # to a native executable.
        native = run(
            "--java-options", "-Xmx2g", "HaplotypeCaller",
            "-I", str(BAM), "-O", str(work / "native.vcf"),
            "--java-options=-Dfoo=bar", "--dry-run",
            env={"FASTGATK_HC_BINARY": str(fake_native)},
        )
        native_plan = json_stdout(native)
        assert native_plan["execution_mode"] == "native"
        assert native_plan["launcher_args"] == [
            "--java-options", "-Xmx2g", "--java-options=-Dfoo=bar"
        ]
        assert all("java-options" not in token for token in native_plan["argv"])

        # A real native invocation verifies that the stripping is not merely a
        # dry-run presentation detail.
        run(
            "-java-options", "-Xms1g", "HaplotypeCaller",
            "-I", str(BAM), "-O", str(work / "native-real.vcf"),
            "--java-options", "-Dbar=baz",
            env={"FASTGATK_HC_BINARY": str(fake_native),
                 "FASTGATK_CAPTURE": str(capture)},
        )
        captured = json.loads(capture.read_text(encoding="utf-8"))
        assert "--java-options" not in captured
        assert "-java-options" not in captured
        assert str(BAM) in captured

        # Nested @args expansion is bounded, deterministic, and works across
        # the launcher/tool separator.  @@ remains a literal @ escape.
        nested = work / "nested.args"
        nested.write_text(
            f"HaplotypeCaller -I {BAM} -O {work / 'nested.vcf'}\n",
            encoding="utf-8",
        )
        outer = work / "outer.args"
        outer.write_text(f"-- @{nested}\n", encoding="utf-8")
        nested_plan = json_stdout(run(
            f"@{outer}", "--dry-run",
            env={"FASTGATK_HC_BINARY": str(fake_native)},
        ))
        assert nested_plan["execution_mode"] == "native"
        assert nested_plan["launcher_separator"] is True
        assert str(BAM) in nested_plan["argv"]
        sys.path.insert(0, str(DISPATCHER.parent))
        import fastgatk as dispatcher_module  # pylint: disable=import-outside-toplevel
        assert dispatcher_module.expand_argument_files(["@@literal"], cwd=work) == ["@literal"]

        # GATK properties can alter Java/codec/tool defaults outside the
        # native registry.  The dispatcher therefore routes config-bearing
        # calls to Java and preserves both launcher controls in fallback argv.
        config = work / "GATKConfig.properties"
        config.write_text("samjdk.compression_level = 2\n", encoding="utf-8")
        config_plan = json_stdout(run(
            "--dry-run", "--java-options=-Xmx3g",
            "HaplotypeCaller", "-I", str(BAM), "-O", str(work / "config.vcf"),
            "--gatk-config-file", str(config),
            env={"FASTGATK_GATK_BINARY": "/opt/gatk"},
        ))
        assert config_plan["execution_mode"] == "fallback"
        assert config_plan["fallback_reason"].startswith("--gatk-config-file")
        fallback_argv = config_plan["argv"]
        assert fallback_argv[:2] == ["/opt/gatk", "--java-options=-Xmx3g"]
        assert "--gatk-config-file" in fallback_argv
        assert str(config) in fallback_argv
        assert fallback_argv.index("HaplotypeCaller") > fallback_argv.index(str(config))

        # The post-tool -- group is retained for explicit Java fallback.  A
        # native call without fallback fails closed instead of dropping the
        # trailing argument or passing an unknown option to the binary.
        separator_plan = json_stdout(run(
            "--dry-run", "--fallback", "HaplotypeCaller",
            "-I", str(BAM), "-O", str(work / "spark.vcf"), "--",
            "--spark-runner", "LOCAL", "--unknown-spark-option", "kept",
            env={"FASTGATK_GATK_BINARY": "/opt/gatk"},
        ))
        assert separator_plan["execution_mode"] == "fallback"
        separator_argv = separator_plan["argv"]
        assert "--" in separator_argv
        assert "--unknown-spark-option" in separator_argv
        fail_closed = run(
            "HaplotypeCaller", "-I", str(BAM), "-O", str(work / "bad.vcf"),
            "--", "--unknown-spark-option", "kept",
            env={"FASTGATK_HC_BINARY": str(fake_native)}, expect=2,
        )
        assert json.loads(fail_closed.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_PARAMETER"

        # Unknown options remain fail-closed before the capture-only native
        # process is launched.
        unknown = run(
            "HaplotypeCaller", "-I", str(BAM), "-O", str(work / "bad2.vcf"),
            "--definitely-unknown-launcher-or-tool-option",
            env={"FASTGATK_HC_BINARY": str(fake_native)}, expect=2,
        )
        assert json.loads(unknown.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_PARAMETER"
        unknown_launcher = run(
            "--definitely-unknown-launcher-option", "HaplotypeCaller",
            "-I", str(BAM), "-O", str(work / "bad3.vcf"),
            env={"FASTGATK_HC_BINARY": str(fake_native)}, expect=2,
        )
        assert json.loads(unknown_launcher.stderr.splitlines()[-1])["category"] == "UNSUPPORTED_PARAMETER"
        print(json.dumps({
            "status": "pass",
            "checks": 10,
            "launcher": "java-options+gatk-config-file+recursive-args+separator",
            "unknown_option_policy": "fail-closed",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
