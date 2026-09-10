#!/usr/bin/env python3
"""Release-pinned GATK oracle for VariantsToTable row-shape boundaries.

This deliberately covers behaviors that are easy to get subtly wrong while
still looking plausible in a wide TSV: split multi-allelic List distribution,
the special Number=R -ASF projection, the fact that GATK's molten writer does
not emit -ASF/-ASGF fields, and -EMD's distinction between a missing value and
a literal String value of ``NA``.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get(
    "FASTGATK_VARIANTS_TO_TABLE_BINARY",
    Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
    / "fastgatk-variants-to-table",
))
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>
##INFO=<ID=RINFO,Number=R,Type=Integer,Description=Ref and alt values>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=TXT,Number=1,Type=String,Description=Literal text>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tZETA\tALPHA
"""


def run(command: list[str], *, expected: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if expected and result.returncode != 0:
        raise AssertionError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}")
    return result


def table_bytes(path: Path) -> bytes:
    content = path.read_bytes()
    assert content.startswith(b"RecordID\t") or content.startswith(b"CHROM\t"), content[:120]
    return content


def main() -> int:
    if not BINARY.is_file() or not os.access(BINARY, os.X_OK):
        raise SystemExit(f"missing native VariantsToTable binary: {BINARY}")
    if not JAVA.is_file() or not GATK_JAR.is_file():
        raise SystemExit("pinned GATK 4.6.2.0 runtime is required for this oracle")

    with tempfile.TemporaryDirectory(prefix="fastgatk-variants-to-table-oracle-") as directory:
        work = Path(directory)
        input_vcf = work / "input.vcf"
        input_vcf.write_text(
            HEADER
            + "chr1\t1\trs1\tA\tG,T\t50\tPASS\tAC=2,3;RINFO=10,2,3\tGT:AD:TXT\t0/2:8,1,5:NA\t0/1:8,4,2:present\n"
            + "chr1\t2\trs2\tC\tT\t50\tPASS\tAC=1;RINFO=7,1\tGT:AD:TXT\t0/0:5,0:NA\t0/1:3,2:present\n",
            encoding="utf-8",
        )

        def native_output(name: str, arguments: list[str]) -> Path:
            output = work / f"native-{name}.tsv"
            run([str(BINARY), "-V", str(input_vcf), "-O", str(output), *arguments])
            return output

        def gatk_output(name: str, arguments: list[str]) -> Path:
            output = work / f"gatk-{name}.tsv"
            run([str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "VariantsToTable",
                 "-V", str(input_vcf), "-O", str(output), *arguments])
            return output

        # On normal TSV output, Number=A fields are distributed by ALT and
        # -ASF Number=R first removes REF then distributes the remaining ALT
        # entries.  The exact Java table also fixes sample sort/order.
        split_args = ["-SMA", "-F", "CHROM", "-F", "ALT", "-F", "AC",
                      "-ASF", "RINFO", "-GF", "GT", "-ASGF", "AD"]
        native_split = native_output("split", split_args)
        gatk_split = gatk_output("split", split_args)
        assert table_bytes(native_split) == table_bytes(gatk_split), (
            gatk_split.read_text(encoding="utf-8"), native_split.read_text(encoding="utf-8"))
        split_lines = native_split.read_text(encoding="utf-8").splitlines()
        assert split_lines[0] == "CHROM\tALT\tAC\tRINFO\tALPHA.GT\tALPHA.AD\tZETA.GT\tZETA.AD"
        assert split_lines[1:3] == [
            "chr1\tG\t2\t 2\tA/G\t8,4\tA/T\t8,1",
            "chr1\tT\t3\t 3\tA/G\t8,2\tA/T\t8,5",
        ]

        # In the pinned Java implementation emitMoltenizedOutput iterates
        # fieldsToTake and genotypeFieldsToTake, not asFields/asGenotypeFields.
        # Keep that slightly surprising public behavior byte-identical.
        molten_args = ["-SMA", "--moltenize", "-F", "CHROM", "-F", "ALT",
                       "-ASF", "RINFO", "-GF", "GT", "-ASGF", "AD"]
        native_molten = native_output("molten", molten_args)
        gatk_molten = gatk_output("molten", molten_args)
        assert table_bytes(native_molten) == table_bytes(gatk_molten), (
            gatk_molten.read_text(encoding="utf-8"), native_molten.read_text(encoding="utf-8"))
        molten_text = native_molten.read_text(encoding="utf-8")
        assert "\tRINFO\t" not in molten_text and "\tAD\t" not in molten_text

        # `NA` may be a real String FORMAT value.  -EMD is an absence check,
        # so it must not reject this input; byte equality proves native tracks
        # value presence independently of its display token.
        literal_args = ["-F", "CHROM", "-GF", "TXT", "-EMD"]
        native_literal = native_output("literal-na", literal_args)
        gatk_literal = gatk_output("literal-na", literal_args)
        assert table_bytes(native_literal) == table_bytes(gatk_literal), (
            gatk_literal.read_text(encoding="utf-8"), native_literal.read_text(encoding="utf-8"))
        assert any(line.endswith("\tNA") for line in
                   native_literal.read_text(encoding="utf-8").splitlines()[1:])

        # A genuinely absent annotation remains fail-closed for both tools.
        missing_args = ["-F", "CHROM", "-F", "DOES_NOT_EXIST", "-EMD"]
        native_missing = run([str(BINARY), "-V", str(input_vcf), "-O", str(work / "native-missing.tsv"),
                              *missing_args], expected=False)
        gatk_missing = run([str(JAVA), "-Xmx1g", "-jar", str(GATK_JAR), "VariantsToTable",
                            "-V", str(input_vcf), "-O", str(work / "gatk-missing.tsv"),
                            *missing_args], expected=False)
        assert native_missing.returncode != 0, native_missing.stderr
        assert gatk_missing.returncode != 0, gatk_missing.stderr

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "split_multi_allelic_exact": True,
            "moltenize_exact": True,
            "moltenize_excludes_as_fields": True,
            "error_if_missing_literal_na_exact": True,
            "error_if_missing_absent_field_fail_closed": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
