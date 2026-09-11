#!/usr/bin/env python3
"""Pinned GATK oracle for homozygous ``--use-iupac-sample`` calls.

FastaAlternateReferenceMaker calls GATK's ``getIUPACBase`` for every SNP.
That helper returns the genotype-selected base for homozygous calls, including
hom-ref and hom-ALT.  A native implementation that only emits IUPAC codes for
heterozygotes can therefore silently choose ALT[0] for a hom-ALT genotype on a
multi-ALT record.  This bounded oracle fixes that direct-replacement boundary.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if check and result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout}\nstderr={result.stderr}")
    return result


def fasta_bytes(path: Path) -> bytes:
    return path.read_bytes()


def fasta_sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines()
                   if not line.startswith(">"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    native = Path(os.environ.get(
        "FASTGATK_FASTA_ALTERNATE_REFERENCE_MAKER_BINARY",
        build / "fastgatk-fasta-alternate-reference-maker"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_fasta_alternate_iupac_hom_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"FastaAlternateReferenceMaker IUPAC oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    sequence = "".join(line.strip() for line in reference.read_text().splitlines()
                         if not line.startswith(">"))
    def ref_base(position: int) -> str:
        return sequence[position - 1].upper()

    # Use one multi-ALT SNP record for each homozygous genotype class.  The
    # 2/2 case is the regression: Java must choose ALT[1], not ALT[0].
    positions = (180, 190, 200)
    rows: list[str] = []
    expected_bases: list[str] = []
    for position, genotype in zip(positions, ("2/2", "1/1", "0/0")):
        reference_base = ref_base(position)
        alt1, alt2 = {"A": ("C", "G"), "C": ("G", "T"),
                      "G": ("T", "A"), "T": ("A", "C")}[reference_base]
        rows.append(f"17\t{position}\t.\t{reference_base}\t{alt1},{alt2}\t50\tPASS\t.\tGT\t{genotype}\n")
        expected_bases.append({"2/2": alt2, "1/1": alt1, "0/0": reference_base}[genotype])

    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=17,length=1000000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    with tempfile.TemporaryDirectory(prefix="fastgatk-fasta-iupac-hom-oracle-") as directory:
        work = Path(directory)
        variants = work / "homozygous-multialt.vcf"
        variants.write_text(header + "".join(rows), encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(variants)])
        java_output = work / "java.fasta"
        native_output = work / "native.fasta"
        common = ["-R", str(reference), "-V", str(variants),
                  "--use-iupac-sample", "SAMPLE", "-L", "17:175-205", "-O"]
        run([str(java), "-jar", str(gatk), "FastaAlternateReferenceMaker", *common, str(java_output)])
        run([str(native), *common, str(native_output)])
        if fasta_bytes(java_output) != fasta_bytes(native_output):
            raise AssertionError(
                "native homozygous IUPAC FASTA differs from pinned GATK:\n"
                f"java={java_output.read_text()}native={native_output.read_text()}")
        java_sequence = fasta_sequence(java_output)
        for position, expected in zip(positions, expected_bases):
            observed = java_sequence[position - 175]
            if observed != expected:
                raise AssertionError(
                    f"pinned GATK did not select expected homozygous base at {position}: "
                    f"expected={expected}, observed={observed}")

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "homozygous_iupac_selected_allele_exact": True,
        "multi_alt_hom_alt2_regression": True,
        "genotypes": ["2/2", "1/1", "0/0"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
