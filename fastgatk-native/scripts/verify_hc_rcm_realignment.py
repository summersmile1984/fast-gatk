#!/usr/bin/env python3
"""Regression for HC's default RCM read-to-haplotype realignment path.

The default must exercise real SW/CIGAR projection and report non-zero
telemetry.  The explicit switch is retained as a backwards-compatible CLI
alias and must remain output-neutral on the pinned fixture.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def decompressed_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        print(json.dumps({"status": "skipped", "reason": "fixture or native binary missing"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-rcm-realign-") as directory:
        work = Path(directory)
        outputs: list[Path] = []
        manifests: list[Path] = []
        for name, extra in (("default", []), ("realigned", [
            "--use-haplotype-realignment-for-rcm"])):
            output = work / f"{name}.g.vcf.gz"
            manifest = work / f"{name}.json"
            subprocess.run([
                str(native), "-I", str(bam), "-R", str(reference), "-L", "17:69000-70000",
                "-O", str(output), "-ERC", "GVCF", "--min-depth", "1",
                "--min-alt-support", "1", "--threads", "2",
                "--create-output-variant-index", "false", "--output-manifest", str(manifest),
                *extra,
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            outputs.append(output)
            manifests.append(manifest)
        baseline = decompressed_lines(outputs[0])
        candidate = decompressed_lines(outputs[1])
        assert baseline == candidate, "RCM realignment changed the pinned VCF/GVCF output"
        default_telemetry = json.loads(manifests[0].read_text(encoding="utf-8"))["telemetry"]
        telemetry = json.loads(manifests[1].read_text(encoding="utf-8"))["telemetry"]
        assert default_telemetry["rcm_haplotype_realignment_used"] is True
        assert default_telemetry["rcm_realigned_observations"] > 0
        assert telemetry["rcm_haplotype_realignment_used"] is True
        assert telemetry["rcm_realigned_observations"] > 0
        assert telemetry["rcm_realignment_fallback_observations"] > 0
        print(json.dumps({
            "status": "pass",
            "records": sum(not line.startswith("#") for line in candidate),
            "rcm_realigned_observations": telemetry["rcm_realigned_observations"],
            "rcm_realignment_fallback_observations": telemetry[
                "rcm_realignment_fallback_observations"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
