#!/usr/bin/env python3
"""Assert chr17 1 kb HC walks GATK ActiveRegions independently.

GATK HaplotypeCaller genotypes each ActiveRegion as its own
assemble→likelihood→genotype unit.  The native host must recurse one
assemble→likelihood→genotype pass per padded ActivityRegion on
17:69000-70000.  A region with no retained EventMap event has no PairHMM
allele matrix, so groups can be fewer than scheduled regions; locally
assembled haplotypes must reach PairHMM and VCF alleles still include GATK
4.6.2.0 17:69067 T>G and (with
--dont-use-soft-clipped-bases) 17:69298 A>T.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def records(path: Path) -> list[tuple[str, int, str, str]]:
    return [
        (fields[0], int(fields[1]), fields[3], fields[4])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
        for fields in [line.split("\t")]
    ]


def run(binary: Path, bam: Path, reference: Path, output: Path, manifest: Path,
        extra: list[str]) -> None:
    command = [
        str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-70000",
        "-O", str(output), "--output-manifest", str(manifest),
        "--threads", "2", "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
        *extra,
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise AssertionError(f"native HC failed: {result.stderr[-2000:]}")


def assert_independent(manifest: Path) -> dict:
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    telemetry = payload["telemetry"]
    regions = telemetry["assembly_regions"]
    groups = telemetry["pairhmm_assembly_region_groups"]
    union = telemetry["assembly_region_union_count"]
    assert regions >= 2, {"assembly_regions": regions}
    assert union == regions, {
        "assembly_region_union_count": union,
        "assembly_regions": regions,
    }
    assert 0 < groups <= regions, {
        "pairhmm_assembly_region_groups": groups,
        "assembly_regions": regions,
    }
    assert telemetry["pairhmm_assembly_region_partitioned"] is True
    # `candidate_sites` is the sum of the per-region pre-genotyping sets.
    # Its exact count depends on duplicate EventMap keys which GATK folds into
    # a single site, but every emitted call must have originated in that set.
    assert telemetry["candidate_sites"] >= telemetry["variant_calls"], telemetry
    return telemetry


def assert_haplotypes_reached_pairhmm(telemetry: dict) -> None:
    assert telemetry["graph_haplotype_paths"] > 0, telemetry["graph_haplotype_paths"]
    assert telemetry["pairhmm_graph_haplotypes"] > 0, {
        "pairhmm_graph_haplotypes": telemetry["pairhmm_graph_haplotypes"],
        "graph_haplotype_paths": telemetry["graph_haplotype_paths"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not binary.is_file() or not bam.is_file() or not reference.is_file():
        raise SystemExit("fastgatk-hc-chr17-independent-regions inputs are required")

    with tempfile.TemporaryDirectory(prefix="hc-chr17-regions-") as directory:
        work = Path(directory)
        default_vcf = work / "default.vcf"
        default_json = work / "default.json"
        run(binary, bam, reference, default_vcf, default_json, [])
        default_calls = records(default_vcf)
        assert ("17", 69067, "T", "G") in default_calls, default_calls
        default_tel = assert_independent(default_json)
        assert_haplotypes_reached_pairhmm(default_tel)

        default_b = work / "default_b.vcf"
        run(binary, bam, reference, default_b, work / "default_b.json", [])
        assert default_vcf.read_bytes() == default_b.read_bytes(), (
            "chr17 default VCF is not deterministic")

        soft_vcf = work / "softclip.vcf"
        soft_json = work / "softclip.json"
        run(binary, bam, reference, soft_vcf, soft_json, ["--dont-use-soft-clipped-bases"])
        soft_calls = records(soft_vcf)
        assert ("17", 69067, "T", "G") in soft_calls, soft_calls
        assert ("17", 69298, "A", "T") in soft_calls, soft_calls
        soft_tel = assert_independent(soft_json)

        soft_b = work / "softclip_b.vcf"
        run(binary, bam, reference, soft_b, work / "softclip_b.json",
            ["--dont-use-soft-clipped-bases"])
        assert soft_vcf.read_bytes() == soft_b.read_bytes(), (
            "chr17 soft-clip VCF is not deterministic")

        print(json.dumps({
            "status": "pass",
            "region": "17:69000-70000",
            "assembly_regions": default_tel["assembly_regions"],
            "pairhmm_assembly_region_groups": default_tel["pairhmm_assembly_region_groups"],
            "independent_units": True,
            "allele_69067_TG": True,
            "allele_69298_AT_without_softclips": True,
            "deterministic_bytes": True,
            "graph_haplotype_paths": default_tel["graph_haplotype_paths"],
            "pairhmm_graph_haplotypes": default_tel["pairhmm_graph_haplotypes"],
            "graph_variant_candidates": default_tel["graph_variant_candidates"],
            "softclip_assembly_regions": soft_tel["assembly_regions"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
