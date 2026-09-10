#!/usr/bin/env python3
"""Pin Mutect2 gVCF concrete-allele ownership to GATK's active cores.

``--debug-assembly-variants-out`` shows EventMap events from padded assembly
regions.  Those halo events are not automatically VariantContexts: only the
``AssemblyResultSet.trimTo()`` active-core owner may hand an allele to
SomaticGenotypingEngine.  This real fixture contains two such right-halo
events (69803 and 69807), making it a compact oracle for both the three
low-TLOD concrete records that must be retained and the two that must not be
emitted.
"""

from __future__ import annotations

import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def concrete_alleles(path: Path) -> list[tuple[str, int, str, str]]:
    rows: list[tuple[str, int, str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[4] == "<NON_REF>":
                continue
            alleles = fields[4].split(",")
            assert alleles[-1] == "<NON_REF>", fields
            rows.append((fields[0], int(fields[1]), fields[3], ",".join(alleles[:-1])))
    return rows


def concrete_data_rows(path: Path) -> list[str]:
    """Return full concrete gVCF records, preserving serialized precision."""
    rows: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[4] != "<NON_REF>":
                rows.append(line.rstrip("\n"))
    return rows


def non_reference_tlods(path: Path) -> dict[tuple[str, int, str, str], str]:
    """Return the symbolic <NON_REF> TLOD from each concrete gVCF row.

    This value does *not* come from SomaticReferenceConfidenceModel's pileup
    LOD.  SomaticGenotypingEngine adds a NON_REF likelihood row to the local
    fragment matrix and evaluates it with the same all-alleles-vs-without-alt
    reduction as concrete EventMap alleles.  Keep the serialized result as a
    compact oracle for that otherwise easy-to-conflate boundary.
    """
    rows: dict[tuple[str, int, str, str], str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[4] == "<NON_REF>":
                continue
            alleles = fields[4].split(",")
            assert alleles[-1] == "<NON_REF>", fields
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            tlods = info["TLOD"].split(",")
            assert len(tlods) == len(alleles), fields
            rows[(fields[0], int(fields[1]), fields[3], ",".join(alleles[:-1]))] = tlods[-1]
    return rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 gVCF oracle inputs are required")
        print('{"status":"skip","reason":"bundled GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-gvcf-eventmap-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.g.vcf.gz"
        native_output = work / "native.g.vcf.gz"
        common = ["-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
                  "-L", "17:69000-70000", "-ERC", "GVCF"]
        gatk_run = run([str(java), "-jar", str(gatk), "Mutect2", *common, "-O", str(gatk_output)])
        assert gatk_run.returncode == 0, gatk_run.stderr[-5000:]
        native_run = run([str(binary), *common, "-O", str(native_output)])
        assert native_run.returncode == 0, native_run.stderr[-5000:]
        expected = concrete_alleles(gatk_output)
        observed = concrete_alleles(native_output)
        assert observed == expected, (expected, observed)
        assert concrete_data_rows(native_output) == concrete_data_rows(gatk_output)
        assert len(observed) == 9, observed
        observed_positions = {position for _contig, position, _ref, _alt in observed}
        assert {69124, 69280, 69657}.issubset(observed_positions), observed
        assert not {69803, 69807}.intersection(observed_positions), observed
        expected_non_reference = non_reference_tlods(gatk_output)
        observed_non_reference = non_reference_tlods(native_output)
        assert observed_non_reference == expected_non_reference, (
            expected_non_reference, observed_non_reference)
        print('{"status":"pass","mode":"Mutect2 -ERC GVCF",'
              '"oracle":"GATK-4.6.2.0","concrete_eventmap_alleles":9,'
              '"concrete_data_rows":"exact","non_reference_tlod":"exact"}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
