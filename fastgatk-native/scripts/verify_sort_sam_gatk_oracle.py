#!/usr/bin/env python3
"""Pinned Picard/GATK 4.6.2.0 oracle for the SortSam replacement slice.

The comparison is deliberately limited to stable SortSam semantics: HTSJDK
coordinate/queryname ordering, the output @HD sort-order declaration, and the
CREATE_INDEX default/explicit contract.  It does not claim Picard byte
identity, compression parity, cloud I/O, or duplicate ordering, which has its
own release-pinned comparator oracle.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


SAM = """@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:100
@SQ\tSN:chr2\tLN:100
@RG\tID:rg1\tSM:S1
zeta\t0\tchr2\t5\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:rg1
tie\t256\tchr1\t10\t20\t4M\t*\t0\t0\tAAAA\t!!!!\tRG:Z:rg1\tHI:i:2
tie\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tHI:i:1
tie\t16\tchr1\t10\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:rg1\tHI:i:1
chr1-late\t0\tchr1\t30\t60\t4M\t*\t0\t0\tCCCC\tIIII\tRG:Z:rg1
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t!!!!\tRG:Z:rg1
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def fail(command: list[str], result: subprocess.CompletedProcess[str]) -> None:
    raise AssertionError(
        f"command failed ({result.returncode}): {' '.join(command)}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def sam_records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise AssertionError(f"malformed SAM record in {path}: {line}")
        # Optional tag order is writer-level.  Keep values and types in the
        # semantic comparison while making their order irrelevant.
        tags = tuple(sorted(fields[11:]))
        result.append(tuple(fields[:11]) + tags)
    return result


def hd_sort_order(path: Path) -> str | None:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("@HD\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:]
                      if ":" in field)
        return fields.get("SO")
    return None


def invoke_pair(root: Path, source: Path, work: Path, label: str,
                sort_order: str, create_index: bool | None) -> None:
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_SORT_SAM_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-sort-sam"
    picard_output = work / f"{label}.picard.sam"
    native_output = work / f"{label}.native.sam"
    picard = [str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
              "-I", str(source), "-O", str(picard_output), "-SO", sort_order,
              "--CREATE_INDEX", "false"]
    native_command = [str(native), "-I", str(source), "-O", str(native_output),
                      "--sort-order", sort_order, "--create-output-bam-index=false"]
    if create_index is not None:
        picard[-1] = "true" if create_index else "false"
        native_command[-1] = "--create-output-bam-index=true" if create_index else "--create-output-bam-index=false"
    picard_result = run(picard)
    if picard_result.returncode != 0:
        fail(picard, picard_result)
    native_result = run(native_command)
    if native_result.returncode != 0:
        fail(native_command, native_result)
    if sam_records(picard_output) != sam_records(native_output):
        raise AssertionError({
            "label": label,
            "picard": sam_records(picard_output),
            "native": sam_records(native_output),
        })
    if hd_sort_order(picard_output) != sort_order or hd_sort_order(native_output) != sort_order:
        raise AssertionError({"label": label, "picard_so": hd_sort_order(picard_output),
                              "native_so": hd_sort_order(native_output), "expected": sort_order})


def invoke_index_pair(root: Path, source: Path, work: Path, label: str,
                      create_index: bool | None) -> None:
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_SORT_SAM_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-sort-sam"
    picard_output = work / f"{label}.picard.bam"
    native_output = work / f"{label}.native.bam"
    value = None if create_index is None else ("true" if create_index else "false")
    picard = [str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
              "-I", str(source), "-O", str(picard_output), "-SO", "coordinate"]
    native_command = [str(native), "-I", str(source), "-O", str(native_output),
                      "--sort-order", "coordinate"]
    if value is not None:
        picard += ["--CREATE_INDEX", value]
        native_command += ["--create-output-bam-index", value]
    picard_result = run(picard)
    if picard_result.returncode != 0:
        fail(picard, picard_result)
    native_result = run(native_command)
    if native_result.returncode != 0:
        fail(native_command, native_result)
    # Picard replaces the .bam suffix (out.bam -> out.bai), while the native
    # HTSlib writer follows its existing sidecar contract (out.bam.bai).
    picard_index = picard_output.with_suffix(".bai")
    native_index = Path(f"{native_output}.bai")
    expected = bool(create_index)
    if picard_index.exists() != expected or native_index.exists() != expected:
        raise AssertionError({"label": label, "expected_index": expected,
                              "picard_index": picard_index.exists(),
                              "native_index": native_index.exists()})


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_SORT_SAM_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-sort-sam"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not (native.is_file() and os.access(native, os.X_OK) and java.is_file() and gatk.is_file()):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned SortSam GATK/Picard oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned GATK/Picard oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-sort-sam-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        source.write_text(SAM, encoding="utf-8")
        invoke_pair(root, source, work, "coordinate", "coordinate", False)
        invoke_pair(root, source, work, "queryname", "queryname", False)
        # Picard's default CREATE_INDEX=false is pinned separately from the
        # explicit true path, which remains required for indexed BAM tasks.
        invoke_index_pair(root, source, work, "index-default", None)
        invoke_index_pair(root, source, work, "index-explicit", True)
    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "sort_orders": ["coordinate", "queryname"],
                      "create_index_default": False,
                      "create_index_explicit": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
