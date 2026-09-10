#!/usr/bin/env python3
"""Verify indexed HC region streaming bounds Host staging and preserves loci."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, args: list[str], extra_env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    if extra_env:
        environment.update(extra_env)
    return subprocess.run([str(binary), *args], text=True, capture_output=True,
                          check=False, env=environment)


def records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def keys(rows: list[str]) -> list[tuple[str, str, str, str]]:
    return [(fields[0], fields[1], fields[3], fields[4])
            for fields in (row.split("\t") for row in rows)]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    assert binary.is_file() and os.access(binary, os.X_OK)
    assert bam.is_file() and Path(f"{bam}.bai").is_file()
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    assert reference.is_file() and Path(f"{reference}.fai").is_file()
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-region-streaming-") as directory:
        work = Path(directory)
        common = ["-I", str(bam), "--min-depth", "1", "--min-alt-support", "1"]
        normal = work / "normal.vcf"
        result = run(binary, [*common, "-O", str(normal)])
        assert result.returncode == 0, result.stderr

        one_tile = work / "one-tile.vcf"
        one_manifest = work / "one-tile.manifest.json"
        result = run(binary, [*common, "--stream-by-region", "1000000",
                              "-O", str(one_tile), "--output-manifest", str(one_manifest)])
        assert result.returncode == 0, result.stderr
        assert normal.read_text(encoding="utf-8") == one_tile.read_text(encoding="utf-8")
        one_metadata = json.loads(one_manifest.read_text(encoding="utf-8"))
        one_telemetry = one_metadata["telemetry"]
        assert one_telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert one_telemetry["pipeline_decoded_items"] == one_telemetry["pipeline_computed_items"]
        assert one_telemetry["pipeline_computed_items"] == one_telemetry["pipeline_encoded_items"] == 1
        assert one_telemetry["pipeline_peak_decoded_bytes"] > 0
        assert one_telemetry["pipeline_peak_computed_bytes"] > 0
        assert one_telemetry["pipeline_peak_encoded_bytes"] > 0

        tiled = work / "tiled.vcf"
        tiled_manifest = work / "tiled.manifest.json"
        result = run(binary, [*common, "--stream-by-region", "500",
                              "-O", str(tiled), "--output-manifest", str(tiled_manifest)])
        assert result.returncode == 0, result.stderr
        # SCOPE OF THIS ASSERTION - read before relying on it.
        # This holds for THIS fixture, but it is NOT a general property of
        # --stream-by-region, and the general form is provably unattainable:
        # GATK's own gVCF output depends on where the -L interval starts
        # (reference-block granularity, annotations, and even phasing change
        # with the window), while the streamed path evaluates each tile window
        # separately and stitches the results.  A tile therefore faithfully
        # reproduces a non-streamed GATK run ON ITS OWN WINDOW, which is not the
        # same as a single run on the user's interval.
        # Measured counterexample: fixtures/chr20/mnp.bam with
        # -L 20:10019901-10020710 (48 records non-streamed and in GATK) yields
        # only 20 records at --stream-by-region 500, because 29 <NON_REF> blocks
        # collapse into one; the tile's output is byte-identical to a
        # non-streamed run on that tile's halo window, and pinned GATK produces
        # the same collapsed block for that window.
        # Full evidence:
        # fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md
        # (sections "Wave 0" and "Wave 1").  Treat this as a fixture-scoped
        # regression guard, not as a parity claim for the flag.
        assert keys(records(tiled)) == keys(records(normal))
        metadata = json.loads(tiled_manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        assert telemetry["stream_by_region"] is True
        assert telemetry["stream_indexed"] is True
        assert telemetry["streamed_regions"] > 1
        assert telemetry["streamed_peak_host_bytes"] > 0
        assert metadata["compatibility"]["stream_by_region"] is True
        assert metadata["compatibility"]["stream_indexed"] is True
        assert telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] == telemetry["streamed_regions"]
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0

        # Force a small SLURM allocation so the adaptive queue must split the
        # first indexed core.  The reference-backed interval keeps the
        # assembly domain fixed; output must remain byte-identical to the
        # unsplit normal run while telemetry proves that the budget was
        # actually enforced rather than merely reported.
        bounded = work / "bounded.vcf"
        bounded_manifest = work / "bounded.manifest.json"
        bounded_env = {"SLURM_MEM_PER_NODE": "155K", "SLURM_CPUS_PER_TASK": "1"}
        bounded_common = [*common, "-R", str(reference), "-L", "17:69000-70000"]
        bounded_normal = work / "bounded-normal.vcf"
        result = run(binary, [*bounded_common, "-O", str(bounded_normal)])
        assert result.returncode == 0, result.stderr
        # A fixed stream core is only a serialization boundary.  Its halo
        # must own the full local activity/AssemblyRegion graph before calls
        # are cropped back to the core; otherwise an arbitrary 500bp split
        # can create a native-only indel near 17:69333.
        bounded_tiled = work / "bounded-tiled.vcf"
        result = run(binary, [*bounded_common, "--stream-by-region", "500",
                              "-O", str(bounded_tiled)])
        assert result.returncode == 0, result.stderr
        assert bounded_normal.read_text(encoding="utf-8") == \
            bounded_tiled.read_text(encoding="utf-8")
        result = run(binary, [*bounded_common, "--stream-by-region", "1000000",
                             "-O", str(bounded), "--output-manifest", str(bounded_manifest)],
                     bounded_env)
        assert result.returncode == 0, result.stderr
        bounded_metadata = json.loads(bounded_manifest.read_text(encoding="utf-8"))
        bounded_telemetry = bounded_metadata["telemetry"]
        # Under an intentionally tiny allocation the caller must split the
        # indexed tile.  Each split recomputes activity with its own halo, so
        # exact per-site evidence can differ from the unbounded reference run;
        # enforce the bounded-path invariants instead of claiming byte parity.
        bounded_rows = records(bounded)
        assert bounded_rows
        bounded_keys = keys(bounded_rows)
        assert bounded_keys == sorted(
            bounded_keys, key=lambda key: (key[0], int(key[1]), key[2], key[3]))
        assert len(bounded_keys) == len(set(bounded_keys))
        assert all(69000 <= int(key[1]) <= 70000 for key in bounded_keys)
        assert bounded_telemetry["streamed_region_splits"] > 0
        assert bounded_telemetry["streamed_peak_host_bytes"] <= 155 * 1024
        assert bounded_telemetry["pipeline_decoded_items"] == bounded_telemetry["pipeline_computed_items"]
        assert bounded_telemetry["pipeline_computed_items"] == bounded_telemetry["pipeline_encoded_items"]

        # A tile boundary is an execution detail, not a GVCF block boundary.
        # Compare a bounded normal GVCF to the same interval split into three
        # indexed tiles.  The stitcher must recover the exact block END/DP/
        # MIN_DP/PL representation, including the rounded median depth built
        # from per-site metadata rather than from tile-level medians.
        gvcf_common = [*common, "-R", str(reference), "-L", "17:69000-69100", "-ERC", "GVCF"]
        normal_gvcf = work / "normal.g.vcf"
        result = run(binary, [*gvcf_common, "-O", str(normal_gvcf)])
        assert result.returncode == 0, result.stderr
        tiled_gvcf = work / "tiled.g.vcf"
        result = run(binary, [*gvcf_common, "--stream-by-region", "50",
                              "-O", str(tiled_gvcf)])
        assert result.returncode == 0, result.stderr
        assert normal_gvcf.read_text(encoding="utf-8") == tiled_gvcf.read_text(encoding="utf-8")

        # Exclusions must be applied at projected-locus granularity.  Compare
        # the native key set with pinned GATK after excluding one of the three
        # broad-fixture calls; a read crossing that locus must not cause the
        # neighboring included calls to disappear.
        excluded_native = work / "excluded.vcf"
        excluded_result = run(binary, [*bounded_common, "-XL", "17:69368-69368",
                                       "-O", str(excluded_native)])
        assert excluded_result.returncode == 0, excluded_result.stderr
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if java.is_file() and gatk_jar.is_file():
            excluded_gatk = work / "excluded-gatk.vcf"
            excluded_gatk_result = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
                "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
                "-XL", "17:69368-69368", "-O", str(excluded_gatk),
                "--create-output-variant-index", "false",
                "--seconds-between-progress-updates", "1",
            ], text=True, capture_output=True, check=False)
            assert excluded_gatk_result.returncode == 0, excluded_gatk_result.stderr
            assert keys(records(excluded_native)) == keys(records(excluded_gatk))
            assert all(position != "69368" for _, position, _, _ in keys(records(excluded_native)))

    print(json.dumps({"status": "pass", "one_tile_exact": True,
                      "tiled_loci_exact": True, "adaptive_split_bounded": True,
                      "gvcf_cross_tile_exact": True,
                      "indexed": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
