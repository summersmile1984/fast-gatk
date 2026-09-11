#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for SelectVariants filtered-genotype no-call.

The important ordering is the same as GATK: sample subsetting, then
``--set-filtered-gt-to-nocall``, then ``--remove-unused-alternates`` and the
resulting AC/AN/AF refresh.  A filtered genotype that is the sole ALT carrier
therefore disappears from the selected record when non-variants are excluded.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t1\trs1\tA\tG,T\t50\tPASS\tAC=3,2;AN=6;AF=0.5,0.333333;DP=20\tGT:FT:AD:PL:GQ\t0/1:LowDP:5,7,0:50,0,70,80,90,100:50\t0/0:PASS:10,0,0:0,50,80,90,100,110:20\t1|2:LowQual:1,4,5:100,90,80,70,60,0:60
chr1\t2\tno_call_filtered\tC\tT\t40\tPASS\tAC=1;AN=4;AF=0.25;DP=12\tGT:FT:AD:PL:GQ\t./.:LowDP:6,0:0,50,80:.\t0/1:PASS:4,8:50,0,50:50\t0/0:PASS:10,0:0,50,80:20
chr1\t3\tsole_filtered_alt\tG\tA\t60\tPASS\tAC=2;AN=6;AF=0.333333;DP=9\tGT:FT:AD:PL:GQ\t1/1:LowQual:0,9:80,40,0:60\t0/0:PASS:9,0:0,40,80:20\t0/0:PASS:9,0:0,40,80:20
"""


def invoke(java: Path, jar: Path, args: list[str]) -> None:
    result = subprocess.run([str(java), "-jar", str(jar), "SelectVariants", *args],
                            text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"GATK SelectVariants failed: {result.stderr[-3000:]}")


def native(binary: Path, args: list[str]) -> dict:
    result = subprocess.run([str(binary), *args], text=True,
                            capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"native SelectVariants failed: {result.stderr[-3000:]}")
    return json.loads(result.stdout.splitlines()[-1])


def parse_atom(value: str):
    if value in {".", ""}:
        return None
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def records(path: Path) -> list[dict]:
    result: list[dict] = []
    with path.open(encoding="utf-8") as handle:
        samples: list[str] = []
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, object] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, equals, value = item.partition("=")
                    info[key] = (True if not equals else
                                 ([parse_atom(token) for token in value.split(",")]
                                  if "," in value else parse_atom(value)))
            fmt = fields[8].split(":") if len(fields) > 8 else []
            sample_payload = {}
            for name, payload in zip(samples, fields[9:]):
                values = payload.split(":")
                sample_payload[name] = {
                    key: values[index] if index < len(values) else "."
                    for index, key in enumerate(fmt)
                }
            result.append({
                "site": tuple(fields[index] for index in (0, 1, 2, 3, 4, 5, 6)),
                "info": info,
                "samples": sample_payload,
            })
    return result


def projected_rows(rows: list[dict]) -> list[tuple]:
    """Compare the semantic fields affected by filtered-genotype no-call."""
    projected = []
    for row in rows:
        info = row["info"]
        sample_values = tuple(
            (name, payload.get("GT"), payload.get("FT"), payload.get("AD"),
             payload.get("PL"), payload.get("GQ"))
            for name, payload in sorted(row["samples"].items())
        )
        projected.append((row["site"],
                          info.get("AC"), info.get("AN"), info.get("AF"),
                          sample_values))
    return projected


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SELECT_VARIANTS_BINARY",
        root / "fastgatk-native/build/fastgatk-select-variants"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_select_variants_filtered_nocall_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("SelectVariants filtered-genotype oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-select-variants-filtered-nocall-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")

        java_plain = work / "java.vcf"
        invoke(java, jar, ["-V", str(source), "--set-filtered-gt-to-nocall", "true",
                           "-O", str(java_plain)])
        java_rows = records(java_plain)
        if [row["site"][1] for row in java_rows] != ["1", "2", "3"]:
            raise AssertionError(f"unexpected Java no-call rows: {java_rows}")

        native_plain = work / "native.vcf"
        native_summary = native(binary, ["-V", str(source), "--set-filtered-gt-to-nocall",
                                         "true", "--create-output-variant-index=false",
                                         "-O", str(native_plain)])
        native_rows = records(native_plain)
        if projected_rows(native_rows) != projected_rows(java_rows):
            raise AssertionError(("native no-call projection differs from GATK",
                                  projected_rows(native_rows), projected_rows(java_rows)))
        if native_summary.get("filtered_gt_nocall") != 3:
            raise AssertionError(f"unexpected native changed-genotype count: {native_summary}")

        # The same operation must remain sample-index safe after HTSlib
        # subset/translate.  This specifically exercises the phased filtered
        # S3 genotype and its AC/AN/AF refresh in a one-sample output.
        java_subset = work / "java-subset.vcf"
        invoke(java, jar, ["-V", str(source), "--sample-name", "S3",
                           "--set-filtered-gt-to-nocall", "true", "-O", str(java_subset)])
        native_subset = work / "native-subset.vcf"
        subset_summary = native(binary, ["-V", str(source), "--sample-name", "S3",
                                        "--set-filtered-gt-to-nocall", "true",
                                        "--create-output-variant-index=false",
                                        "-O", str(native_subset)])
        if projected_rows(records(native_subset)) != projected_rows(records(java_subset)):
            raise AssertionError(("sample-subset no-call projection differs from GATK",
                                  projected_rows(records(native_subset)),
                                  projected_rows(records(java_subset))))
        if subset_summary.get("filtered_gt_nocall") != 1:
            raise AssertionError(f"unexpected sample-subset changed-genotype count: {subset_summary}")

        # A no-call with a non-PASS FT is not counted a second time.  The
        # false form is also a direct replacement guard: it must preserve the
        # input genotype payload and must not trigger an AC/AN/AF refresh.
        native_false = work / "native-false.vcf"
        false_summary = native(binary, ["-V", str(source),
                                        "--set-filtered-gt-to-nocall=false",
                                        "--create-output-variant-index=false",
                                        "-O", str(native_false)])
        if projected_rows(records(native_false)) != projected_rows(records(source)):
            raise AssertionError("explicit false unexpectedly changed SelectVariants genotypes")
        if false_summary.get("filtered_gt_nocall") != 0:
            raise AssertionError(f"false form changed genotypes: {false_summary}")

        # GATK performs the no-call conversion before removing unused ALTs.
        # Consequently rows 1 and 3 become non-variant and are removed; row 2
        # retains its called T allele and its AC/AN/AF reflect the remaining
        # samples.  This catches an implementation that compacts first.
        java_compacted = work / "java-compacted.vcf"
        invoke(java, jar, ["-V", str(source), "--set-filtered-gt-to-nocall", "true",
                           "--remove-unused-alternates", "true",
                           "--exclude-non-variants", "true", "-O", str(java_compacted)])
        native_compacted = work / "native-compacted.vcf"
        compact_summary = native(binary, ["-V", str(source), "--set-filtered-gt-to-nocall",
                                         "--remove-unused-alternates",
                                         "--exclude-non-variants", "--create-output-variant-index=false",
                                         "-O", str(native_compacted)])
        java_compact_rows = records(java_compacted)
        native_compact_rows = records(native_compacted)
        if projected_rows(native_compact_rows) != projected_rows(java_compact_rows):
            raise AssertionError(("native compacted no-call projection differs from GATK",
                                  projected_rows(native_compact_rows),
                                  projected_rows(java_compact_rows)))
        if [row["site"][1] for row in native_compact_rows] != ["2"]:
            raise AssertionError(f"unexpected compacted positions: {native_compact_rows}")
        if compact_summary.get("filtered_gt_nocall") != 3:
            raise AssertionError(f"unexpected compacted changed-genotype count: {compact_summary}")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(java_rows),
            "filtered_gt_nocall": native_summary.get("filtered_gt_nocall"),
            "sample_subset_filtered_gt_nocall": subset_summary.get("filtered_gt_nocall"),
            "compacted_records": len(native_compact_rows),
            "compacted_positions": [row["site"][1] for row in native_compact_rows],
            "false_form_unchanged": false_summary.get("filtered_gt_nocall") == 0,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
