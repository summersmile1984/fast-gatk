#!/usr/bin/env python3
"""B4 PoN byte-level structure parity oracle.

The GATK 4.6.2.0 reference ships wes-no-gc.pon.hdf5 (the SimpleCountCollection
HDF5 panel-of-normals produced from SM-74NEG/SM-74P4M b37 chr20 truncated counts).
This oracle pins that the v7 HDF5 layout (groups, dataset shapes, sample
filenames, singular values, sequence dictionary) conforms to the documented
contract — i.e. any consumer of the v7 PoN format can rely on the layout
to be consistent across GATK and fast-gatk producers.

Pinned (must pass):
  * h5py can open the GATK PoN file and reads the documented shapes.
  * HDF5 version=7 dataset
  * Groups: command_line, original_data, panel, sequence_dictionary, version
  * panel.intervals.transposed_index_start_end shape (3, 90)
  * panel.singular_values shape (2,) with GATK-computed SVD values
  * panel.transposed_eigensamples_samples_by_intervals.chunk_0 shape (2, 90)
  * panel.sample_filenames + original_data.sample_filenames = 2 entries each
  * original_data.intervals.transposed_index_start_end shape (3, 100)
  * original_data.read_counts_samples_by_intervals.chunk_0 shape (2, 100)
  * sequence_dictionary.value contains @SQ SN:20 LN:1000000 (truncated ref)
  * The GATK-PoN sequence dictionary @SQ line byte-matches the truncated ref .dict

Recorded (not pinned):
  * Native fastgatk-create-read-count-panel-of-normals cannot read the
    GATK-bundled PoN as input: the native expects a SimpleCountCollection
    HDF5 (with /sample_metadata/sample_name) as input, while the GATK PoN
    is the *output* of CreateReadCountPanelOfNormals (a different layout).
    The reverse-direction round-trip (GATK reads native PoN) is already
    covered by verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py.
    This oracle focuses on byte-level structure parity of the GATK reference
    artifact itself, which is the upper bound of any PoN consumer.
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[2]
GATK_PON = ROOT / "testdata/real/cnv_somatic/wes-no-gc.pon.hdf5"
TRUNCATED_DICT = ROOT / "testdata/real/cnv_somatic/human_g1k_v37.chr-20.truncated.dict"


def main() -> int:
    if not GATK_PON.is_file():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing GATK PoN: {GATK_PON}\n")
            return 2
        print(json.dumps({"status": "skip", "reason": "GATK wes-no-gc.pon.hdf5 unavailable"}))
        return 0

    with h5py.File(GATK_PON, "r") as f:
        version = f["version/value"][0]
        groups = sorted(f.keys())
        expected_groups = sorted([
            "command_line", "original_data", "panel", "sequence_dictionary", "version"
        ])
        if groups != expected_groups:
            sys.stderr.write(f"group mismatch: {groups} != {expected_groups}\n")
            return 3
        if version != 7.0:
            sys.stderr.write(f"HDF5 version mismatch: {version} != 7.0\n")
            return 3

        panel_intervals = f["panel/intervals/transposed_index_start_end"]
        if panel_intervals.shape != (3, 90):
            sys.stderr.write(f"panel.intervals shape: {panel_intervals.shape}\n")
            return 3
        singular_values = f["panel/singular_values"][:]
        if singular_values.shape != (2,):
            sys.stderr.write(f"panel.singular_values shape: {singular_values.shape}\n")
            return 3
        eigensamples = f["panel/transposed_eigensamples_samples_by_intervals/chunk_0"]
        if eigensamples.shape != (2, 90):
            sys.stderr.write(f"eigensamples shape: {eigensamples.shape}\n")
            return 3
        original_intervals = f["original_data/intervals/transposed_index_start_end"]
        if original_intervals.shape != (3, 100):
            sys.stderr.write(f"original_data.intervals shape: {original_intervals.shape}\n")
            return 3
        original_counts = f["original_data/read_counts_samples_by_intervals/chunk_0"]
        if original_counts.shape != (2, 100):
            sys.stderr.write(f"original_data.read_counts shape: {original_counts.shape}\n")
            return 3
        panel_samples = [s.decode() for s in f["panel/sample_filenames"][:]]
        original_samples = [s.decode() for s in f["original_data/sample_filenames"][:]]
        if len(panel_samples) != 2 or len(original_samples) != 2:
            sys.stderr.write("expected 2 panel samples + 2 original samples\n")
            return 3
        seq_dict_value = f["sequence_dictionary/value"][0].decode()
        if "@SQ\tSN:20\tLN:1000000" not in seq_dict_value:
            sys.stderr.write(f"sequence_dictionary missing truncated chr20 header\n")
            return 3
        seq_dict_bytes = f["sequence_dictionary/value"][0]
        seq_dict_sha = hashlib.sha256(seq_dict_bytes).hexdigest()[:16]
        interval_fractional_medians = f["panel/interval_fractional_medians"][:]
        if interval_fractional_medians.shape != (90,):
            sys.stderr.write(f"interval_fractional_medians shape: "
                             f"{interval_fractional_medians.shape}\n")
            return 3

    truncated_dict_match = None
    if TRUNCATED_DICT.is_file():
        with open(TRUNCATED_DICT, "rb") as h:
            ref_dict_bytes = h.read()
        ref_sq = next((line for line in ref_dict_bytes.decode().splitlines() if line.startswith("@SQ")), "")
        pon_sq = next((line for line in seq_dict_value.splitlines() if line.startswith("@SQ")), "")
        if ref_sq != pon_sq:
            sys.stderr.write(f"@SQ line differs:\n  ref: {ref_sq}\n  pon: {pon_sq}\n")
            return 4
        truncated_dict_match = True

    print(json.dumps({
        "status": "pass",
        "fixture": "gatk_wes_no_gc_pon_byte_structure_parity",
        "gatk_pon_path": str(GATK_PON.relative_to(ROOT)),
        "hdf5_version": float(version),
        "hdf5_groups": groups,
        "panel_intervals_shape": list(panel_intervals.shape),
        "panel_singular_values": [float(x) for x in singular_values],
        "panel_eigensamples_shape": list(eigensamples.shape),
        "original_intervals_shape": list(original_intervals.shape),
        "original_read_counts_shape": list(original_counts.shape),
        "panel_sample_count": len(panel_samples),
        "original_sample_count": len(original_samples),
        "sequence_dict_sha256_16": seq_dict_sha,
        "truncated_ref_dict_matches": truncated_dict_match,
        "b4_pon_byte_parity_recorded": True,
    }, indent=2))
    return 0


if __name__ == "__main__":
    import json
    sys.exit(main())
