#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for ApplyVQSR -XL/-ixp traversal semantics."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_APPLY_VQSR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-apply-vqsr"),
))


HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def run(command: list[str], label: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if check and result.returncode != 0:
        raise AssertionError(f"{label} failed ({result.returncode}):\n{result.stderr}")
    return result


def records(path: pathlib.Path) -> list[list[str]]:
    text = path.read_text(encoding="utf-8")
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


def write_fixture(work: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    input_vcf = work / "input.vcf"
    recal_vcf = work / "recal.vcf"
    input_rows = "".join(
        f"chr1\t{position}\t.\tA\tG\t50\tPASS\tQD=10\tGT:DP\t0/1:{20 + position}\n"
        for position in range(1, 8)
    )
    input_vcf.write_text(HEADER + input_rows, encoding="utf-8")
    recal_header = HEADER.replace(
        "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n", "")
    recal_header = recal_header.replace(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>\n",
        "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=VQSLOD>\n",
    ).replace("\tFORMAT\tSAMPLE\n", "\n")
    recal_rows = "".join(
        f"chr1\t{position}\t.\tN\t<VQSR>\t.\tPASS\tEND={position};VQSLOD={'2.0' if position % 2 else '-2.0'};culprit=QD\n"
        for position in range(1, 8)
    )
    recal_vcf.write_text(recal_header + recal_rows, encoding="utf-8")
    return input_vcf, recal_vcf


def index_feature(path: pathlib.Path) -> None:
    run([str(JAVA), "-jar", str(GATK), "IndexFeatureFile", "-I", str(path)], "IndexFeatureFile")


def invoke(java: bool, input_vcf: pathlib.Path, recal_vcf: pathlib.Path,
           output: pathlib.Path, manifest: pathlib.Path, extra: list[str]) -> None:
    if java:
        command = [str(JAVA), "-jar", str(GATK), "ApplyVQSR", "-V", str(input_vcf),
                   "--recal-file", str(recal_vcf), "--lod-score-cutoff", "0", "--mode", "SNP",
                   "--create-output-variant-index", "false", "--QUIET", "true", *extra,
                   "-O", str(output)]
    else:
        command = [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
                   "--lod-score-cutoff", "0", "--mode", "SNP",
                   "--create-output-variant-index", "false", *extra, "-O", str(output),
                   "--output-manifest", str(manifest)]
    run(command, "GATK ApplyVQSR" if java else "native ApplyVQSR")


def comparable(row: list[str]) -> tuple[str, str, str, str, str, str, str, str]:
    # Java and native use different shortest float spellings; compare the
    # score numerically while requiring exact site/filter/INFO semantics.
    info = {item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in row[7].split(";") if "=" in item}
    return (row[0], row[1], row[3], row[4], row[6], info.get("VQSLOD", ""),
            info.get("culprit", ""), row[8] if len(row) > 8 else "")


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not JAVA.exists() or not GATK.exists():
        oracle_guard.oracle_not_verified('verify_apply_vqsr_exclude_intervals_gatk_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ApplyVQSR interval oracle requires the pinned GATK 4.6.2.0 jar")
        print(json.dumps({"status": "skip", "reason": "pinned GATK unavailable"}))
        return
    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-xl-oracle-") as directory:
        work = pathlib.Path(directory)
        input_vcf, recal_vcf = write_fixture(work)
        index_feature(input_vcf)
        index_feature(recal_vcf)

        # -XL is applied after the include set.  Without -L, it still removes
        # loci from the implicit whole-reference traversal.
        cases = {
            "include_exclude": ["-L", "chr1:1-7", "-XL", "chr1:2-3"],
            "exclude_only": ["-XL", "chr1:2-3"],
            # -ixp expands each excluded span in 1-based inclusive space.
            "exclude_padding": ["-XL", "chr1:4-4", "-ixp", "1"],
        }
        expected_positions = {
            "include_exclude": [1, 4, 5, 6, 7],
            "exclude_only": [1, 4, 5, 6, 7],
            "exclude_padding": [1, 2, 6, 7],
        }
        for name, extra in cases.items():
            java_output = work / f"{name}.java.vcf"
            native_output = work / f"{name}.native.vcf"
            native_manifest = work / f"{name}.native.json"
            invoke(True, input_vcf, recal_vcf, java_output, work / "unused.json", extra)
            invoke(False, input_vcf, recal_vcf, native_output, native_manifest, extra)
            java_rows = records(java_output)
            native_rows = records(native_output)
            assert [int(row[1]) for row in java_rows] == expected_positions[name], (name, java_rows)
            assert [int(row[1]) for row in native_rows] == expected_positions[name], (name, native_rows)
            assert len(java_rows) == len(native_rows)
            for java_row, native_row in zip(java_rows, native_rows):
                assert comparable(java_row)[:5] == comparable(native_row)[:5], (name, java_row, native_row)
                assert abs(float(comparable(java_row)[5]) - float(comparable(native_row)[5])) < 1e-6
                assert comparable(java_row)[6:] == comparable(native_row)[6:]
            manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
            assert manifest["interval_excluded_records"] == 2 if name != "exclude_padding" else manifest["interval_excluded_records"] == 3
            assert manifest["interval_exclusion_padding"] == (1 if name == "exclude_padding" else 0)
        print(json.dumps({
            "status": "pass",
            "cases": len(cases),
            "exclude_positions_exact": True,
            "padding_exact": True,
            "java_native_vqslod_filter_exact": True,
            "manifest_telemetry": True,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
