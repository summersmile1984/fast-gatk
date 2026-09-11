#!/usr/bin/env python3
"""Compare the supported LeftAlignAndTrimVariants surface with GATK 4.6.2.0.

The comparison is deliberately record-semantic rather than header-byte based:
HTSlib and HTSJDK order FORMAT/INFO keys differently and format AF with a
different number of trailing zeroes.  It still checks every emitted record,
site field, INFO annotation, FORMAT value, and genotype value.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[list[str]]:
    opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
    with opener as handle:
        return [line.rstrip("\n").split("\t") for line in handle
                if line and not line.startswith("#")]


def number(value: str) -> tuple[str, object]:
    if value == ".":
        return ("missing", None)
    try:
        return ("number", float(value))
    except ValueError:
        return ("text", value)


def info_map(field: str) -> dict[str, tuple[str, object] | list[tuple[str, object]]]:
    if field == ".":
        return {}
    result: dict[str, tuple[str, object] | list[tuple[str, object]]] = {}
    for item in field.split(";"):
        if "=" not in item:
            result[item] = ("flag", True)
            continue
        key, value = item.split("=", 1)
        values = [number(part) for part in value.split(",")]
        result[key] = values[0] if len(values) == 1 else values
    return result


def format_map(record: list[str]) -> list[dict[str, str]]:
    if len(record) < 10 or record[8] == ".":
        return []
    names = record[8].split(":")
    return [{name: values[index] if index < len(values) else "."
             for index, name in enumerate(names)}
            for values in (sample.split(":") for sample in record[9:])]


def assert_same(native: list[list[str]], gatk: list[list[str]], label: str) -> None:
    assert len(native) == len(gatk), f"{label}: record count {len(native)} != {len(gatk)}"
    for index, (left, right) in enumerate(zip(native, gatk)):
        assert left[:7] == right[:7], f"{label}/{index}: site fields differ: {left[:7]} != {right[:7]}"
        left_info = info_map(left[7])
        right_info = info_map(right[7])
        assert left_info.keys() == right_info.keys(), (
            f"{label}/{index}: INFO keys differ: {left_info} != {right_info}")
        for key in left_info:
            lv, rv = left_info[key], right_info[key]
            lvalues = lv if isinstance(lv, list) else [lv]
            rvalues = rv if isinstance(rv, list) else [rv]
            assert len(lvalues) == len(rvalues), f"{label}/{index}/{key}: INFO width differs"
            for lvalue, rvalue in zip(lvalues, rvalues):
                if lvalue[0] == rvalue[0] == "number":
                    assert abs(float(lvalue[1]) - float(rvalue[1])) < 1e-5, (
                        f"{label}/{index}/{key}: {lvalue} != {rvalue}")
                else:
                    assert lvalue == rvalue, f"{label}/{index}/{key}: {lvalue} != {rvalue}"
        left_formats = format_map(left)
        right_formats = format_map(right)
        assert len(left_formats) == len(right_formats), f"{label}/{index}: sample count differs"
        for sample, (left_fields, right_fields) in enumerate(zip(left_formats, right_formats)):
            assert left_fields.keys() == right_fields.keys(), (
                f"{label}/{index}/sample{sample}: FORMAT keys differ: "
                f"{left_fields} != {right_fields}")
            assert left_fields == right_fields, (
                f"{label}/{index}/sample{sample}: FORMAT differs: "
                f"{left_fields} != {right_fields}")


def run_case(binary: Path, java: Path, jar: Path, source: Path, reference: Path,
             work: Path, label: str, *options: str) -> dict:
    native = work / f"{label}.native.vcf.gz"
    gatk = work / f"{label}.gatk.vcf"
    native_result = subprocess.run(
        [str(binary), "-V", str(source), "-R", str(reference), "-O", str(native), *options],
        text=True, capture_output=True, check=False)
    assert native_result.returncode == 0, f"{label}: native failed: {native_result.stderr}"
    gatk_result = subprocess.run(
        [str(java), "-Xmx1g", "-jar", str(jar), "LeftAlignAndTrimVariants",
         "-V", str(source), "-R", str(reference), "-O", str(gatk), *options],
        text=True, capture_output=True, check=False)
    assert gatk_result.returncode == 0, f"{label}: GATK failed: {gatk_result.stderr}"
    assert_same(records(native), records(gatk), label)
    return json.loads(native_result.stdout.splitlines()[-1])


def index_for_gatk(java: Path, jar: Path, source: Path) -> None:
    result = subprocess.run(
        [str(java), "-Xmx1g", "-jar", str(jar), "IndexFeatureFile", "-I", str(source)],
        text=True, capture_output=True, check=False)
    assert result.returncode == 0, f"index {source}: GATK failed: {result.stderr}"


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-left-align-trim"
    java = root / "third_party/jdk17/bin/java"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not binary.is_file() or not java.is_file() or not jar.is_file():
        oracle_guard.oracle_not_verified('verify_left_align_gatk_oracle.py', java, jar)
        print(json.dumps({"status": "skip", "reason": "native or pinned GATK assets unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-left-align-oracle-") as temp:
        work = Path(temp)
        reference = work / "reference.fa"
        reference.write_text(">chr1\n" + "AC" * 30 + "\n", encoding="utf-8")
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:60\n", encoding="utf-8")
        reference.with_suffix(".fa.fai").write_text(
            "chr1\t60\t6\t60\t61\n", encoding="utf-8")
        header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=60>
##FILTER=<ID=q10,Description=Low quality>
##ALT=<ID=DEL,Description=Deletion>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""
        source = work / "source.vcf"
        source.write_text(
            header
            + "chr1\t4\tmnp\tAC\tCA\t60\tq10\tDP=12\tGT\t0/1\t0/0\n"
            + "chr1\t9\tmulti\tAC\tA,ACAC\t50\tPASS\tDP=20;AC=1,1;AN=4;AF=0.25,0.25\t"
              "GT:AD:PL:GQ:DP\t0/2:8,2,10:80,50,90,60,70,0:1:9\t"
              "0/1:10,4,3:80,50,90,60,70,0:1:10\n"
            + "chr1\t20\thet\tAC\tA,ACAC\t50\tPASS\tDP=20;AC=1,1;AN=4;AF=0.25,0.25\t"
              "GT:AD:PL:GQ:DP\t1/2:8,2,10:80,50,90,60,70,0:1:9\t"
              "0/1:10,4,3:80,50,90,60,70,0:1:10\n"
            + "chr1\t30\tsym\tC\tG,<DEL>\t50\tPASS\tDP=20;AC=1,1;AN=4;AF=0.25,0.25\t"
              "GT:AD:PL:GQ:DP\t0/1:8,2,10:80,50,90,60,70,0:1:9\t"
              "0/0:10,4,3:80,50,90,60,70,0:1:10\n",
            encoding="utf-8")

        split_summary = run_case(binary, java, jar, source, reference, work, "split",
                                 "--split-multi-allelics")
        assert split_summary["interval_skipped"] == 0
        run_case(binary, java, jar, source, reference, work, "split-original",
                 "--split-multi-allelics", "--keep-original-ac")
        run_case(binary, java, jar, source, reference, work, "split-no-trim",
                 "--split-multi-allelics", "--dont-trim-alleles")

        # Interval selection happens before normalization.  The MNP at 4 is
        # excluded, while the two split records at 9 are retained and still
        # use GATK's previous-written-record window.
        index_for_gatk(java, jar, source)
        run_case(binary, java, jar, source, reference, work, "interval",
                 "--split-multi-allelics", "-L", "chr1:9-9")

        # An END span is selected by overlap, not only by POS.  This fixture
        # also exercises the empty-output branch for a non-variant split.
        span_header = header.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate allele>\n", 1)
        span = work / "span.g.vcf"
        span.write_text(span_header + "chr1\t1\tblock\tA\t<NON_REF>\t.\tPASS\tEND=6\tGT\t0/0\t0/0\n",
                        encoding="utf-8")
        index_for_gatk(java, jar, span)
        run_case(binary, java, jar, span, reference, work, "span", "-L", "chr1:4-4")

    print(json.dumps({"status": "pass", "oracle": "GATK-4.6.2.0",
                      "surface": ["split", "orig-counts", "no-trim", "interval", "END-overlap"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
