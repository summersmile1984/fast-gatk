#!/usr/bin/env python3
"""Verify that ActivityProfile/AssemblyRegion controls cross the tool boundary."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    with tempfile.TemporaryDirectory(prefix="fastgatk-activity-region-") as directory:
        work = Path(directory)
        output = work / "calls.vcf.gz"
        manifest = work / "calls.manifest.json"
        command = [
            str(native), "-I", str(bam), "-R", str(reference), "-L", "17:69000-70000",
            "-O", str(output), "--min-depth", "1", "--min-alt-support", "1",
            "--active-probability-threshold", "0.002",
            "--assembly-region-padding", "8",
            "--max-assembly-region-size", "100",
            "--max-prob-propagation-distance", "50",
            "--output-manifest", str(manifest),
        ]
        result = json.loads(subprocess.check_output(command, text=True).splitlines()[-1])
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        assert result["status"] == "prototype"
        assert telemetry["active_probability_threshold"] == 0.002
        assert telemetry["assembly_region_padding"] == 8
        assert telemetry["max_assembly_region_size"] == 100
        assert telemetry["max_probability_propagation_distance"] == 50
        assert telemetry["activity_filter_size"] == 50
        assert telemetry["activity_effective_max_probability_propagation_distance"] == 100
        assert telemetry["assembly_regions"] > 0
        assert metadata["compatibility"]["activity_profile_controls"] is True
        assert metadata["compatibility"]["activity_bandpass"] is True
        print(json.dumps({
            "status": "pass",
            "assembly_regions": telemetry["assembly_regions"],
            "active_probability_threshold": telemetry["active_probability_threshold"],
            "assembly_region_padding": telemetry["assembly_region_padding"],
            "max_assembly_region_size": telemetry["max_assembly_region_size"],
            "max_probability_propagation_distance": telemetry["max_probability_propagation_distance"],
            "activity_filter_size": telemetry["activity_filter_size"],
            "activity_effective_max_probability_propagation_distance": telemetry[
                "activity_effective_max_probability_propagation_distance"
            ],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
