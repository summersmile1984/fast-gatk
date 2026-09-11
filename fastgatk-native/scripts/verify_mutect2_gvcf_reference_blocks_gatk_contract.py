#!/usr/bin/env python3
"""Pin Mutect2 ``-ERC GVCF`` reference blocks to GATK 4.6.2.0.

This deliberately uses a no-variant part of the bundled chr17 fixture.  It
isolates SomaticReferenceConfidenceModel's per-locus TLOD calculation and
SomaticGVCFBlockCombiner's TLOD-band aggregation from the still-separate
EventMap/variant-row parity work.  The two runs cover GATK's default bands and
an explicit custom band list, comparing compressed data records byte for byte.
"""

from __future__ import annotations

import gzip
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def vcf_data(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line and not line.startswith("#")]


def vcf_header(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line.startswith("#")]


def invoke(java: Path, gatk: Path, binary: Path, reference: Path, bam: Path,
           output: Path, custom_bands: bool) -> tuple[Path, Path]:
    gatk_output = output.with_name(f"gatk-{output.name}")
    native_output = output.with_name(f"native-{output.name}")
    common = ["-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
              "-L", "17:69000-69066", "-ERC", "GVCF"]
    if custom_bands:
        common += ["-LODB", "-1.0", "-LODB", "0.0", "-LODB", "1.0"]
    gatk_run = run([str(java), "-jar", str(gatk), "Mutect2", *common, "-O", str(gatk_output)])
    assert gatk_run.returncode == 0, gatk_run.stderr[-5000:]
    native_run = run([str(binary), *common, "-O", str(native_output)])
    assert native_run.returncode == 0, native_run.stderr[-5000:]
    return gatk_output, native_output


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
        oracle_guard.oracle_not_verified('verify_mutect2_gvcf_reference_blocks_gatk_contract.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 gVCF oracle inputs are required")
        print('{"status":"skip","reason":"bundled GATK oracle unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-gvcf-blocks-") as directory:
        work = Path(directory)
        for name, custom_bands in (("default.g.vcf.gz", False), ("custom.g.vcf.gz", True)):
            gatk_output, native_output = invoke(
                java, gatk, binary, reference, bam, work / name, custom_bands)
            gatk_rows, native_rows = vcf_data(gatk_output), vcf_data(native_output)
            assert native_rows == gatk_rows, (name, gatk_rows, native_rows)
            if custom_bands:
                # All values in this window are in the [-1.0, 0.0) band, so
                # SomaticGVCFBlockCombiner emits one block under the custom
                # list.  The values below are the source-format median/min
                # depth and the minimum raw TLOD, respectively.
                assert native_rows == [
                    "17\t69000\t.\tT\t<NON_REF>\t.\t.\tEND=69066\t"
                    "GT:DP:MIN_DP:TLOD\t0/0:3:2:-7.782e-01"
                ], native_rows
            else:
                # GATK's 0.5-wide default bands split this interval at the
                # same four TLOD boundaries; only the explicit broad custom
                # list above coalesces it into one record.
                assert len(native_rows) == 4, native_rows
                assert native_rows[0].endswith("GT:DP:MIN_DP:TLOD\t0/0:5:3:-7.782e-01"), native_rows
                assert native_rows[-1].startswith("17\t69065\t.\tC\t<NON_REF>\t.\t.\tEND=69066\t"), native_rows
            native_schema = vcf_header(native_output)
            for prefix in (
                "##ALT=<ID=NON_REF,",
                "##FORMAT=<ID=TLOD,",
                "##FORMAT=<ID=MIN_DP,",
                "##INFO=<ID=END,",
                "##GVCFBlock",
            ):
                assert any(line.startswith(prefix) for line in native_schema), (name, prefix)
        print('{"status":"pass","mode":"Mutect2 -ERC GVCF",'
              '"oracle":"GATK-4.6.2.0","windows":2}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
