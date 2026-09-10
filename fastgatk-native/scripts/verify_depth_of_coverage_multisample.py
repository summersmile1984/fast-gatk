#!/usr/bin/env python3
"""Regression for DepthOfCoverage @RG/SM partitioning.

The fixture is deliberately a SAM stream so the contract does not depend on
an external samtools binary or an index.  It exercises both the all-sample
path and an explicit --sample selection.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_DEPTH_OF_COVERAGE_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-depth-of-coverage"),
))
REFERENCE = ROOT / "gatk-source" / "src" / "test" / "resources" / "human_g1k_v37.chr17_1Mb.fasta"


def main() -> None:
    if not BINARY.exists() or not REFERENCE.exists():
        raise SystemExit("missing DepthOfCoverage multisample assets")
    with tempfile.TemporaryDirectory(prefix="fastgatk-depth-of-coverage-multisample-") as temporary:
        work = pathlib.Path(temporary)
        sam = work / "multi.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:17\tLN:1000000\n"
            "@RG\tID:rg1\tSM:S1\n"
            "@RG\tID:rg2\tSM:S2\n"
            "s1a\t0\t17\t69000\t60\t6M\t*\t0\t0\tAAAAAA\tIIIIII\tRG:Z:rg1\n"
            "s1b\t0\t17\t69000\t60\t6M\t*\t0\t0\tCCCCCC\tIIIIII\tRG:Z:rg1\n"
            "s2a\t0\t17\t69000\t60\t6M\t*\t0\t0\tGGGGGG\tIIIIII\tRG:Z:rg2\n"
            "s2b\t0\t17\t69001\t60\t6M\t*\t0\t0\tTTTTTT\tIIIIII\tRG:Z:rg2\n",
            encoding="utf-8",
        )
        all_out = work / "all.out"
        all_manifest = work / "all.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-I", str(sam), "-L", "17:69000-69005",
             "-O", str(all_out), "--output-manifest", str(all_manifest)],
            check=True, text=True, capture_output=True,
        )
        lines = all_out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "Locus,Total_Depth,Average_Depth_S1,Depth_for_S1,Average_Depth_S2,Depth_for_S2"
        assert lines[1:] == [
            "17:69000,3,2.00,2,1.00,1",
            "17:69001,4,2.00,2,2.00,2",
            "17:69002,4,2.00,2,2.00,2",
            "17:69003,4,2.00,2,2.00,2",
            "17:69004,4,2.00,2,2.00,2",
            "17:69005,4,2.00,2,2.00,2",
        ]
        manifest = json.loads(all_manifest.read_text(encoding="utf-8"))
        assert manifest["samples"] == ["S1", "S2"]
        assert manifest["compatibility"]["multi_sample_read_group_partition"] is True
        selected_out = work / "selected.out"
        subprocess.run(
            [str(BINARY), "-R", str(REFERENCE), "-I", str(sam), "-L", "17:69000-69005",
             "--sample", "S2", "-O", str(selected_out)],
            check=True, text=True, capture_output=True,
        )
        selected_lines = selected_out.read_text(encoding="utf-8").splitlines()
        assert selected_lines[0] == "Locus,Total_Depth,Average_Depth_sample,Depth_for_S2"
        assert selected_lines[1] == "17:69000,1,1.00,1"
        assert selected_lines[2] == "17:69001,2,2.00,2"
        print(json.dumps({
            "status": "pass",
            "tool": "DepthOfCoverage",
            "samples": ["S1", "S2"],
            "selected_sample": "S2",
            "partition": "@RG-to-SM",
        }, sort_keys=True))


if __name__ == "__main__":
    main()
