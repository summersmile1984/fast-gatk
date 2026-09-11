#!/usr/bin/env python3
"""Pinned Picard/GATK 4.6.2.0 oracle for SortSam's Picard CLI boundary.

This slice is intentionally semantic: it exercises the long Picard aliases,
optional Boolean switches, compression/validation values, bounded spill
controls, and the coordinate output writer, then compares HTSJDK/native SAM
records and sort-order headers.  It does not claim byte-identical compression
or support Picard-only sort orders and metadata side effects outside this
bounded contract.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


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


def require_success(command: list[str]) -> None:
    result = run(command)
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise AssertionError(f"malformed SAM record in {path}: {line}")
        result.append(tuple(fields[:11]) + tuple(sorted(fields[11:])))
    return result


def header(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("@")] 


def hd_sort_order(path: Path) -> str | None:
    for line in header(path):
        if line.startswith("@HD\t"):
            fields = dict(field.split(":", 1) for field in line.split("\t")[1:]
                          if ":" in field)
            return fields.get("SO")
    return None


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_SORT_SAM_BINARY",
        os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"),
    ))
    if native.is_dir():
        native /= "fastgatk-sort-sam"
    if not (native.is_file() and os.access(native, os.X_OK) and
            java.is_file() and gatk.is_file()):
        oracle_guard.oracle_not_verified('verify_sort_sam_cli_boundary_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("pinned SortSam CLI oracle assets are required")
        print(json.dumps({"status": "skip", "reason": "pinned GATK/Picard oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-sort-sam-cli-oracle-") as directory:
        work = Path(directory)
        source = work / "input.sam"
        reference = work / "reference.fa"
        source.write_text(SAM, encoding="utf-8")
        reference.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 100 + "\n",
                            encoding="utf-8")
        java_tmp = work / "java-tmp"
        native_tmp = work / "native-tmp"
        java_tmp.mkdir()
        native_tmp.mkdir()
        java_output = work / "picard.sam"
        native_output = work / "native.sam"
        native_manifest = work / "native.manifest.json"
        java_command = [str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
                        "--INPUT", str(source), "--OUTPUT", str(java_output),
                        "--REFERENCE_SEQUENCE", str(reference),
                        "--SORT_ORDER", "coordinate", "--CREATE_INDEX", "false",
                        "--COMPRESSION_LEVEL", "0", "--MAX_RECORDS_IN_RAM", "1",
                        "--TMP_DIR", str(java_tmp), "--VALIDATION_STRINGENCY", "SILENT",
                        "--QUIET", "true", "--VERBOSITY", "WARNING",
                        "--USE_JDK_DEFLATER", "false", "--USE_JDK_INFLATER", "false"]
        native_command = [str(native), "--INPUT", str(source), "--OUTPUT", str(native_output),
                          "--REFERENCE_SEQUENCE", str(reference), "--SORT_ORDER", "coordinate",
                          "--CREATE_INDEX", "false", "--COMPRESSION_LEVEL", "0",
                          "--MAX_RECORDS_IN_RAM", "1", "--TMP_DIR", str(native_tmp),
                          "--VALIDATION_STRINGENCY", "SILENT", "--QUIET", "true",
                          "--VERBOSITY", "WARNING", "--USE_JDK_DEFLATER", "false",
                          "--USE_JDK_INFLATER", "false", "--output-manifest", str(native_manifest)]
        require_success(java_command)
        require_success(native_command)
        if records(java_output) != records(native_output):
            raise AssertionError({"picard_records": records(java_output),
                                  "native_records": records(native_output)})
        if hd_sort_order(java_output) != "coordinate" or hd_sort_order(native_output) != "coordinate":
            raise AssertionError({"picard_so": hd_sort_order(java_output),
                                  "native_so": hd_sort_order(native_output)})
        if header(java_output) != header(native_output):
            raise AssertionError({"picard_header": header(java_output),
                                  "native_header": header(native_output)})
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        compatibility = manifest.get("compatibility", {})
        telemetry = manifest.get("telemetry", {})
        if compatibility.get("compression_level") != 0 or telemetry.get("compression_level") != 0:
            raise AssertionError({"compatibility": compatibility, "telemetry": telemetry})
        if compatibility.get("validation_stringency") != "SILENT" or telemetry.get("validation_stringency") != "SILENT":
            raise AssertionError({"compatibility": compatibility, "telemetry": telemetry})
        if compatibility.get("sort_order") != "coordinate" or manifest.get("outputs", [{}])[0].get("complete") is not True:
            raise AssertionError({"compatibility": compatibility, "outputs": manifest.get("outputs")})
        if (work / "native.sam.bai").exists() or (work / "picard.sam.bai").exists():
            raise AssertionError("CREATE_INDEX=false unexpectedly produced an index for SAM output")
    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "picard_aliases": True, "sort_order": "coordinate",
                      "compression_level": 0, "validation_stringency": "SILENT",
                      "create_index": False, "record_order_exact": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
