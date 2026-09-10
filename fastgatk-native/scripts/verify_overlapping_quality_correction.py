#!/usr/bin/env python3
"""Check the GATK overlapping-mate base-quality correction contract."""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def run(binary: Path, fixture: Path, *extra: str) -> dict[str, object]:
    completed = subprocess.run(
        [str(binary), "-I", str(fixture), "-O", "-", "--min-depth", "100", *extra],
        check=True, text=True, capture_output=True,
    )
    return json.loads(completed.stdout.strip().splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    fixture = root / "fastgatk-native/tests/fixtures/overlap_quality.sam"
    if not binary.exists() or not fixture.exists():
        raise SystemExit("missing HC binary or overlap fixture")

    corrected = run(binary, fixture)
    assert corrected["overlapping_quality_correction_used"] is True
    assert corrected["overlapping_quality_correction_metadata_available"] is True
    assert corrected["overlapping_pairs"] == 2
    assert corrected["overlapping_bases"] == 12
    assert corrected["overlapping_conflicting_bases"] == 6
    # Six concordant bases per mate are capped to Phred 20 (12 bytes); six
    # conflicting bases per mate are zeroed (12 bytes).
    assert corrected["overlapping_quality_caps"] == 24

    disabled = run(binary, fixture, "--do-not-correct-overlapping-quality")
    assert disabled["overlapping_quality_correction_used"] is False
    assert disabled["overlapping_pairs"] == 0
    assert disabled["do_not_correct_overlapping_base_qualities"] is True
    print(json.dumps({"status": "pass", "corrected_pairs": corrected["overlapping_pairs"],
                      "corrected_bases": corrected["overlapping_bases"],
                      "conflicting_bases": corrected["overlapping_conflicting_bases"],
                      "quality_caps": corrected["overlapping_quality_caps"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
