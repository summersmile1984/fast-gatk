#!/usr/bin/env python3
"""Contract checks for the Host VQSR GatherTranches boundary."""

from __future__ import annotations

import csv
import gzip
import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_GATHER_TRANCHES_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-gather-tranches"),
))
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
HEADER = (
    "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,"
    "filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
)


def write_shard(path: pathlib.Path, known_offset: int) -> None:
    path.write_text(
        "# Variant quality score tranches file\n# Version number 6\n"
        + HEADER
        + f"100.0,{10 + known_offset},20,2.0000,1.5000,5.0,VQSRTrancheSNP0.00to100.00,SNP,100,100,1.0000\n"
        + f"99.0,{5 + known_offset},10,1.0000,1.0000,3.0,VQSRTrancheSNP0.00to99.00,SNP,100,90,0.9000\n"
        + f"90.0,{1 + known_offset},5,1.0000,1.0000,1.0,VQSRTrancheSNP0.00to90.00,SNP,100,50,0.5000\n",
        encoding="utf-8",
    )


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-tranches-") as temporary:
        work = pathlib.Path(temporary)
        shard_a = work / "a.tranches"
        shard_b = work / "b.tranches"
        output = work / "gathered.tranches"
        manifest = work / "gathered.json"
        write_shard(shard_a, 0)
        write_shard(shard_b, 2)
        result = subprocess.run(
            [str(BINARY), "-I", str(shard_a), "-I", str(shard_b),
             "--truth-sensitivity-tranche", "100", "-tranche", "90",
             "--mode", "SNP", "-O", str(output), "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=True,
        )
        rows = list(csv.DictReader(
            line for line in output.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ))
        assert [row["targetTruthSensitivity"] for row in rows] == ["100.00", "90.00"], rows
        assert rows[0]["minVQSLod"] == "1.0000", rows
        assert rows[0]["numKnown"] == "4", rows
        assert rows[1]["numKnown"] == "12", rows
        assert rows[1]["callsAtTruthSites"] == "180", rows
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        assert manifest_payload["tool"] == "GatherTranches"
        assert manifest_payload["implementation"] == "fastgatk-gather-tranches"
        assert manifest_payload["status"] == "prototype"
        assert manifest_payload["execution_space"] == "Host"
        assert manifest_payload["determinism"] == "strict"
        assert manifest_payload["primary_output"] == str(output)
        assert manifest_payload["primary_output_kind"] == "tranches"
        assert manifest_payload["compatibility"] == {
            "merge_scattered_tranches": True,
            "stateful_target_selection": True,
            "version6_input": True,
            "gatk_version5_output": True,
            "cloud_uri_staging": False,
            "compressed_output": False,
            "bit_identical_to_gatk": True,
        }
        assert manifest_payload["outputs"] == [{
            "path": str(output), "kind": "tranches", "complete": True,
        }]
        telemetry = manifest_payload["telemetry"]
        assert telemetry["input_files"] == 2
        assert telemetry["merged_tranches"] == 3
        assert telemetry["selected_tranches"] == 2
        assert telemetry["output_bytes"] == output.stat().st_size
        assert telemetry["wall_seconds"] >= 0.0
        assert isinstance(telemetry["resources"], dict)
        assert manifest_payload["merged_tranches"] == 3
        assert json.loads(result.stdout)["selected_tranches"] == 2

        # Real GATK 4.6.2.0 oracle: compare both the version/header and every
        # formatted data row.  This catches the non-obvious stateful target
        # walk and Tranche.tranchesString call-count ordering.
        if not JAVA.exists() or not GATK_JAR.exists():
            raise AssertionError(f"missing GATK oracle runtime: {JAVA} / {GATK_JAR}")
        gatk_output = work / "gatk.tranches"
        subprocess.run(
            [str(JAVA), "-jar", str(GATK_JAR), "GatherTranches",
             "-I", str(shard_a), "-I", str(shard_b),
             "--truth-sensitivity-tranche", "100", "-tranche", "90",
             "--mode", "SNP", "-O", str(gatk_output)],
            text=True, capture_output=True, check=True,
        )
        assert output.read_text(encoding="utf-8") == gatk_output.read_text(encoding="utf-8"), (
            output.read_text(encoding="utf-8"), gatk_output.read_text(encoding="utf-8")
        )

        # VariantRecalibrator's scatter mode can emit an otherwise-valid
        # version-6 tranche header with zero rows for a requested VQSLOD
        # slice.  GATK ignores that slice during gather; accepting it is
        # required for direct scatter/gather replacement rather than merely a
        # hand-crafted non-empty fixture.
        empty_shard = work / "empty.tranches"
        empty_shard.write_text(
            "# Variant quality score tranches file\n# Version number 6\n" + HEADER,
            encoding="utf-8",
        )
        empty_native_output = work / "empty-native.tranches"
        empty_gatk_output = work / "empty-gatk.tranches"
        empty_args = ["-I", str(empty_shard), "-I", str(shard_a),
                      "--mode", "SNP", "-tranche", "100", "-tranche", "90"]
        subprocess.run(
            [str(BINARY), *empty_args, "-O", str(empty_native_output)],
            text=True, capture_output=True, check=True,
        )
        subprocess.run(
            [str(JAVA), "-jar", str(GATK_JAR), "GatherTranches",
             *empty_args, "-O", str(empty_gatk_output)],
            text=True, capture_output=True, check=True,
        )
        assert empty_native_output.read_bytes() == empty_gatk_output.read_bytes(), (
            empty_native_output.read_text(encoding="utf-8"),
            empty_gatk_output.read_text(encoding="utf-8"),
        )

        # --mode is a GATK-required argument; native must fail closed instead
        # of silently selecting SNP.  Check both implementations rather than
        # relying on their intentionally different error text.
        missing_mode_native = subprocess.run(
            [str(BINARY), "-I", str(shard_a), "-O", str(work / "no-mode-native.tranches")],
            text=True, capture_output=True,
        )
        missing_mode_gatk = subprocess.run(
            [str(JAVA), "-jar", str(GATK_JAR), "GatherTranches",
             "-I", str(shard_a), "-O", str(work / "no-mode-gatk.tranches")],
            text=True, capture_output=True,
        )
        assert missing_mode_native.returncode != 0
        assert "--mode is required" in missing_mode_native.stderr
        assert missing_mode_gatk.returncode != 0

        compressed_shard = work / "a.tranches.gz"
        with gzip.open(compressed_shard, "wt", encoding="utf-8") as stream:
            stream.write(shard_a.read_text(encoding="utf-8"))
        compressed_input_output = work / "gathered-from-gz.tranches"
        subprocess.run(
            [str(BINARY), "-I", str(compressed_shard), "-I", str(shard_b),
             "--truth-sensitivity-tranche", "100", "-tranche", "90",
             "--mode", "SNP", "-O", str(compressed_input_output)],
            text=True, capture_output=True, check=True,
        )
        assert compressed_input_output.read_bytes() == output.read_bytes()

        # The output boundary is also stream-safe: a .gz suffix selects zlib
        # compression while preserving the exact GATK tranche text after
        # decompression and recording the choice in the manifest.
        compressed_output = work / "gathered.tranches.gz"
        compressed_manifest = work / "gathered-compressed.json"
        compressed_result = subprocess.run(
            [str(BINARY), "-I", str(shard_a), "-I", str(shard_b),
             "--truth-sensitivity-tranche", "100", "-tranche", "90",
             "--mode", "SNP", "-O", str(compressed_output),
             "--output-manifest", str(compressed_manifest)],
            text=True, capture_output=True, check=True,
        )
        assert gzip.open(compressed_output, "rt", encoding="utf-8").read() == output.read_text(encoding="utf-8")
        compressed_metadata = json.loads(compressed_manifest.read_text(encoding="utf-8"))
        assert compressed_metadata["compatibility"]["compressed_output"] is True
        assert compressed_metadata["telemetry"]["output_compressed"] is True
        assert compressed_metadata["telemetry"]["output_bytes"] == compressed_output.stat().st_size
        assert json.loads(compressed_result.stdout)["selected_tranches"] == 2

        bad_version = work / "bad-version.tranches"
        bad_version.write_text(shard_a.read_text(encoding="utf-8").replace("Version number 6", "Version number 5"), encoding="utf-8")
        invalid_version = subprocess.run(
            [str(BINARY), "-I", str(bad_version), "--mode", "SNP",
             "-O", str(work / "bad-version.out")],
            text=True, capture_output=True,
        )
        assert invalid_version.returncode != 0
        assert "unsupported VQSLOD tranche version" in invalid_version.stderr

        invalid = subprocess.run(
            [str(BINARY), "-I", str(shard_a), "-O", str(work / "invalid.csv"), "--mode", "BAD"],
            text=True, capture_output=True,
        )
        assert invalid.returncode != 0
        assert "--mode must be SNP" in invalid.stderr
        print(json.dumps({
            "status": "pass",
            "input_shards": 2,
            "merged_tranches": 3,
            "selected_tranches": 2,
            "mode_validation": True,
            "gatk_oracle_exact": True,
            "empty_scatter_slice_oracle_exact": True,
            "mode_required": True,
            "short_tranche_alias": True,
            "version_validation": True,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
