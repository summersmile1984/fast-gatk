#!/usr/bin/env python3
"""Pinned CLI/interval boundary oracle for LeftAlignAndTrimVariants.

This fixture covers the VariantWalker controls that are easy for a native
implementation to accidentally parse without applying: include/exclude
intervals and their padding, interval set/merging rules, and Barclay Boolean
spellings for the normalization and index switches.
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


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def canonical(record: list[str]) -> list[str]:
    """Ignore legal INFO key ordering differences between HTSlib/htsjdk."""
    if record[7] == ".":
        info = "."
    else:
        info = ";".join(sorted(record[7].split(";")))
    return record[:7] + [info] + record[8:]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-left-align-trim"
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not all(path.is_file() for path in (binary, java, gatk)):
        oracle_guard.oracle_not_verified('verify_left_align_cli_boundary_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("LeftAlignAndTrimVariants CLI oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "native or pinned GATK assets unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-left-align-cli-oracle-") as directory:
        work = Path(directory)
        reference = work / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 40 + "\n", encoding="utf-8")
        reference.with_suffix(".dict").write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n", encoding="utf-8")
        reference.with_suffix(".fa.fai").write_text(
            "chr1\t40\t6\t40\t41\n", encoding="utf-8")
        source = work / "input.vcf"
        source.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=40>\n"
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>\n"
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t4\tindel\tAA\tA\t50\tPASS\t.\tGT\t0/1\n"
            "chr1\t9\tsnp-left\tA\tC\t51\tPASS\t.\tGT\t0/1\n"
            "chr1\t12\tsnp-right\tA\tG\t52\tPASS\t.\tGT\t0/1\n"
            "chr1\t15\tmulti\tA\tC,G\t53\tPASS\tAC=1,1;AN=2;AF=0.5,0.5\tGT:AD:PL:GQ\t0/2:8,2,10:80,50,90,60,70,0:1\n"
            "chr1\t20\tsnp-last\tA\tT\t54\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8")
        indexed = run([str(java), "-Xmx1g", "-jar", str(gatk), "IndexFeatureFile",
                       "-I", str(source)])
        assert indexed.returncode == 0, indexed.stderr

        def compare(label: str, *options: str, check_index: bool = False) -> dict:
            native_output = work / f"{label}.native.vcf.gz"
            gatk_output = work / f"{label}.gatk.vcf"
            manifest = work / f"{label}.manifest.json"
            native = run([
                str(binary), "-V", str(source), "-R", str(reference), "-O", str(native_output),
                "--output-manifest", str(manifest), *options])
            assert native.returncode == 0, f"{label}: native failed: {native.stderr}"
            java_run = run([
                str(java), "-Xmx1g", "-jar", str(gatk), "LeftAlignAndTrimVariants",
                "-V", str(source), "-R", str(reference), "-O", str(gatk_output), *options])
            assert java_run.returncode == 0, f"{label}: GATK failed: {java_run.stderr}"
            native_records = records(native_output)
            gatk_records = records(gatk_output)
            assert [canonical(record) for record in native_records] == [canonical(record) for record in gatk_records], (
                f"{label}: record mismatch\nnative={native_records}\ngatk={gatk_records}")
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if check_index:
                assert not Path(f"{native_output}.tbi").exists(), f"{label}: native index unexpectedly exists"
            return {
                "records": len(native_records),
                "positions": [int(record[1]) for record in native_records],
                "manifest": metadata["compatibility"],
            }

        cases: dict[str, object] = {}
        # -ip expands the include interval before traversal: chr1:10-10 plus
        # two bases selects the records at 9 and 12.
        cases["include-padding"] = compare(
            "include-padding", "-L", "chr1:10-10", "-ip", "2")
        assert cases["include-padding"]["positions"] == [9, 12]  # type: ignore[index]

        # Exclusions are applied after the include set and ixp pads the
        # excluded interval on both sides; the sites at 9/12 disappear.
        cases["exclude-padding"] = compare(
            "exclude-padding", "-L", "chr1:1-30", "-XL", "chr1:10-10", "-ixp", "2")
        assert cases["exclude-padding"]["positions"] == [1, 15, 20]  # type: ignore[index]

        # Repeated include selectors are UNION by default, and -isr is the
        # Barclay short alias.  Abutting selectors with either imr spelling
        # must retain the same selected sites.
        cases["union"] = compare(
            "union", "-L", "chr1:9-9", "-L", "chr1:12-12", "-isr", "UNION",
            "-imr", "OVERLAPPING_ONLY")
        assert cases["union"]["positions"] == [9, 12]  # type: ignore[index]

        # These switches are Boolean-valued Barclay arguments, not bare-only
        # flags.  Exercise separated false forms while comparing output to
        # Java.
        cases["split-false"] = compare("split-false", "--split-multi-allelics", "false")
        cases["dont-trim-false"] = compare("dont-trim-false", "--dont-trim-alleles", "false")
        cases["keep-original-false"] = compare("keep-original-false", "--keep-original-ac", "false")
        cases["index-short-false"] = compare("index-short-false", "-OVI", "false", check_index=True)

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "cases": cases,
            "interval_exclusion_padding_exact": True,
            "optional_boolean_forms_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
