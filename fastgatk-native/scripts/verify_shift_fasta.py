#!/usr/bin/env python3
"""Contract and GATK-oracle checks for the native ShiftFasta walker."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/mutect/mito/Homo_sapiens_assembly38.mt_only.fasta"
EXAMPLE_REF = ROOT / "gatk-source/src/test/resources/exampleFASTA.fasta"
EXPECTED_INTERVALS = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/shifted8000.intervals"
EXPECTED_SHIFTED_INTERVALS = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/shifted8000.shifted.intervals"
SHIFT_BIN = Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build")) / "fastgatk-shift-fasta"
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK_JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout}\nstderr={result.stderr}"
        )
    return result


def sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


def compare_bytes(lhs: Path, rhs: Path, label: str) -> None:
    if lhs.read_bytes() != rhs.read_bytes():
        raise AssertionError(f"{label} differs: {lhs} != {rhs}")


def main() -> int:
    assert REF.exists() and REF.with_suffix(REF.suffix + ".fai").exists()
    assert EXAMPLE_REF.exists() and EXAMPLE_REF.with_suffix(EXAMPLE_REF.suffix + ".fai").exists()
    assert EXPECTED_INTERVALS.exists() and EXPECTED_SHIFTED_INTERVALS.exists()
    assert SHIFT_BIN.exists()
    java_oracle = oracle_guard.oracle_ready('verify_shift_fasta.py', JAVA, GATK_JAR)

    with tempfile.TemporaryDirectory(prefix="fastgatk-shift-fasta-") as directory:
        work = Path(directory)
        native = work / "native.fasta"
        native_chain = work / "native.chain"
        native_base = work / "native-mito"
        native_manifest = work / "native.json"
        run([str(SHIFT_BIN), "-R", str(REF), "-O", str(native),
             "--shift-back-output", str(native_chain), "--shift-offset-list", "8000",
             "--interval-file-name", str(native_base), "--line-width", "61",
             "--threads", "2", "--output-manifest", str(native_manifest)])
        assert native.with_suffix(native.suffix + ".fai").exists()
        assert native.with_suffix(".dict").exists()
        payload = json.loads(native_manifest.read_text())
        assert payload["tool"] == "ShiftFasta"
        assert payload["shifted_contigs"] == 1
        assert payload["records"] == 16569
        assert payload["execute_seconds"] >= 0.0
        assert Path(payload["dictionary"]) == native.with_suffix(".dict")
        compare_bytes(native_base.with_suffix(".intervals"), EXPECTED_INTERVALS,
                      "regular interval sidecar")
        compare_bytes(native_base.with_suffix(".shifted.intervals"), EXPECTED_SHIFTED_INTERVALS,
                      "shifted interval sidecar")

        if java_oracle:
            oracle = work / "oracle.fasta"
            oracle_chain = work / "oracle.chain"
            oracle_base = work / "oracle-mito"
            run([str(JAVA), "-jar", str(GATK_JAR), "ShiftFasta", "-R", str(REF),
                 "-O", str(oracle), "--shift-back-output", str(oracle_chain),
                 "--shift-offset-list", "8000", "--interval-file-name", str(oracle_base),
                 "--line-width", "61"])
            compare_bytes(native, oracle, "ShiftFasta FASTA")
            compare_bytes(native.with_suffix(".dict"), oracle.with_suffix(".dict"),
                          "ShiftFasta sequence dictionary")
            compare_bytes(native_chain, oracle_chain, "ShiftFasta chain")
            compare_bytes(native_base.with_suffix(".intervals"), oracle_base.with_suffix(".intervals"),
                          "ShiftFasta regular sidecar")
            compare_bytes(native_base.with_suffix(".shifted.intervals"),
                          oracle_base.with_suffix(".shifted.intervals"),
                          "ShiftFasta shifted sidecar")

        # The default offset is half of each contig.  Applying the native
        # transform twice therefore returns the original sequence (modulo
        # FASTA wrapping and header/index metadata).
        first = work / "example.shifted.fasta"
        second = work / "example.reshifted.fasta"
        run([str(SHIFT_BIN), "-R", str(EXAMPLE_REF), "-O", str(first),
             "--shift-back-output", str(work / "example.first.chain")])
        run([str(SHIFT_BIN), "-R", str(first), "-O", str(second),
             "--shift-back-output", str(work / "example.second.chain")])
        if sequence(second).upper() != sequence(EXAMPLE_REF).upper():
            raise AssertionError("default ShiftFasta round-trip changed the sequence")

    print(json.dumps({"status": "pass", "java_oracle": java_oracle,
                      "tool": "ShiftFasta"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
