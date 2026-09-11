#!/usr/bin/env python3
"""Pin a large, real no-call HC traversal to the bundled GATK oracle.

The chr17 and dense chr20 fixtures exercise candidate calling in detail.  This
separate 100 kb CEUTrio window contains 90k real reads but no HC calls under
GATK 4.6.2.0's default filters.  It therefore makes the high-volume Host
boundaries observable without conflating them with candidate annotations:
ReadFilter, activity-region construction, assembly, PairHMM handoff and the
empty VCF writer must all agree with GATK's result.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stderr={result.stderr[-1600:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    reference = root / "testdata/chr20/reference/GRCh37.chr20.fa"
    bam = root / "testdata/real/ceutrio/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"
    bam_index = Path(f"{bam}.bai")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    required = (native, reference, Path(f"{reference}.fai"), bam, bam_index, gatk, java)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_hc_chr20_100k_nocall_gatk_oracle.py', gatk, java)
        raise SystemExit("missing real chr20 no-call oracle inputs")

    interval = "20:10000000-10100000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr20-100k-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        native_manifest = work / "native.json"
        run([
            str(java), "-Xmx2g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", interval,
            "-O", str(gatk_vcf), "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ])
        run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", interval,
            "-O", str(native_vcf), "--output-manifest", str(native_manifest),
            "--threads", "1", "--add-output-vcf-command-line", "false",
        ])

        # Command-line provenance is disabled for both calls, so this checks
        # the complete ordered schema and every emitted record, not a field
        # subset.  A regression that leaks a candidate or changes empty-VCF
        # output cannot be hidden behind set comparison.
        assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
            "100 kb real no-call VCF differs from GATK")
        gatk_rows = [line for line in gatk_vcf.read_text(encoding="utf-8").splitlines()
                     if line and not line.startswith("#")]
        native_rows = [line for line in native_vcf.read_text(encoding="utf-8").splitlines()
                       if line and not line.startswith("#")]
        assert not gatk_rows and not native_rows

        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        telemetry = manifest["telemetry"]
        assert manifest["reads"] >= 90000
        assert manifest["compatibility"]["read_filter_gatk_defaults"] is True
        assert telemetry["assembly_regions"] > 0
        assert telemetry["variant_calls"] == 0
        print(json.dumps({
            "status": "pass",
            "interval": interval,
            "reads": manifest["reads"],
            "assembly_regions": telemetry["assembly_regions"],
            "variant_calls": 0,
            "full_vcf_bit_identical": True,
        }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
