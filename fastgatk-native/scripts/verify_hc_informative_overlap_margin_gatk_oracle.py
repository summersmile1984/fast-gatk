#!/usr/bin/env python3
"""Release-pinned HC oracle for informative-read overlap retention.

GATK's ``--allele-informative-reads-overlap-margin`` is a Host-side
``retainEvidence`` interval, not a graph or PairHMM-window control.  Exercise
the default, an explicit default, and a zero-width request against the pinned
chr17 assembly fixture.  The VCF comparison protects every affected annotation
and FORMAT field, while the manifest makes option ownership observable.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout[-2000:]}\nstderr={result.stderr[-4000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.is_file() for path in (native, java, gatk, bam, reference)):
        oracle_guard.oracle_not_verified('verify_hc_informative_overlap_margin_gatk_oracle.py', java, gatk)
        raise SystemExit("missing HC informative-overlap-margin oracle assets")

    shared_java = [
        str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
        "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
        "--native-pair-hmm-threads", "2",
        "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
        "--seconds-between-progress-updates", "1",
    ]
    shared_native = [
        str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
        "--threads", "2", "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
    ]
    cases: dict[str, tuple[int, list[str]]] = {
        "default": (2, []),
        "explicit-default": (2, ["--allele-informative-reads-overlap-margin", "2"]),
        "zero": (0, ["--allele-informative-reads-overlap-margin", "0"]),
    }
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-overlap-margin-") as directory:
        work = Path(directory)
        gatk_bytes: dict[str, bytes] = {}
        native_bytes: dict[str, bytes] = {}
        for label, (margin, arguments) in cases.items():
            gatk_vcf = work / f"gatk-{label}.vcf"
            native_vcf = work / f"native-{label}.vcf"
            manifest = work / f"native-{label}.json"
            run([*shared_java, *arguments, "-O", str(gatk_vcf)])
            run([*shared_native, *arguments, "-O", str(native_vcf),
                 "--output-manifest", str(manifest)])
            gatk_bytes[label] = gatk_vcf.read_bytes()
            native_bytes[label] = native_vcf.read_bytes()
            if gatk_bytes[label] != native_bytes[label]:
                raise AssertionError(f"{label}: native VCF differs from GATK")
            manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            if manifest_data["compatibility"].get("informative_read_overlap_margin") != margin:
                raise AssertionError(f"{label}: missing compatibility margin")
            if manifest_data["telemetry"].get("informative_read_overlap_margin") != margin:
                raise AssertionError(f"{label}: missing telemetry margin")

        if gatk_bytes["default"] != gatk_bytes["explicit-default"]:
            raise AssertionError("GATK default and explicit 2bp margins differ")
        if native_bytes["default"] != native_bytes["explicit-default"]:
            raise AssertionError("native default and explicit 2bp margins differ")
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": "17:69000-70000",
            "cases": {label: {"margin": margin, "vcf_byte_identical": True}
                      for label, (margin, _) in cases.items()},
            "default_equals_explicit_default": True,
            "scope": "Host retain-evidence overlap ownership",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
