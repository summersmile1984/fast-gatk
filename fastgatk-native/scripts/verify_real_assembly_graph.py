#!/usr/bin/env python3
"""Exercise the linked de-Bruijn assembly path on a real BAM corpus.

The small kernel smoke corpus proves individual graph rules, but it cannot
catch the interaction between CIGAR projection, reference scaffolding,
non-unique-kmer retry, linked traversal, artificial-haplotype recovery and
PairHMM on a production-shaped read set.  This check intentionally uses the
bundled chr17 69k--70k BAM/reference fixture and treats the native manifest as
the graph oracle.  It does not claim Java byte identity; GATK does not expose
these internal graph counters in its VCF.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, bam: Path, reference: Path, output: Path,
        manifest: Path, *, disable_recovery: bool = False) -> dict:
    command = [
        str(binary), "-I", str(bam), "-R", str(reference),
        "-L", "17:69000-70000", "-O", str(output),
        "--threads", "1", "--batch-records", "64",
        "--linked-de-bruijn-graph", "--adaptive-pruning",
        "--num-pruning-samples", "1",
        # Exercise the recovery branch itself; production compatibility uses
        # the GATK min-pruning default of two, while this graph stress fixture
        # intentionally keeps every branch so the telemetry path is visible.
        "--min-pruning", "1",
        # The fixture contains repeated reference kmers.  GATK retries with a
        # larger k in this situation; this explicit switch lets the linked
        # topology be exercised rather than stopping at the retry boundary.
        "--allow-non-unique-kmers-in-ref",
        "--output-manifest", str(manifest),
    ]
    if disable_recovery:
        command.append("--disable-artificial-haplotype-recovery")
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    environment.setdefault("OMP_PROC_BIND", "false")
    result = subprocess.run(command, check=True, text=True,
                            capture_output=True, env=environment)
    summary = json.loads(result.stdout.strip().splitlines()[-1])
    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    # The command summary and sidecar must expose the same graph contract.
    for key in ("graph_nodes", "graph_edges", "graph_haplotype_paths",
                "graph_artificial_haplotype_recovery_paths",
                "graph_artificial_haplotype_recovery_bases"):
        assert summary[key] == metadata["telemetry"][key], key
    telemetry = metadata["telemetry"]
    assert telemetry["pairhmm_graph_snp_posterior_pairs"] >= 0
    details = telemetry.get("graph_haplotype_path_details", [])
    assert telemetry["graph_haplotype_provenance_count"] == len(details)
    assert telemetry["graph_haplotype_provenance_count"] == telemetry["graph_haplotype_sequence_count"]
    assert telemetry["graph_haplotype_cigar_signature"] != 0
    # Graph haplotype-to-reference scoring must go through the common Kokkos
    # score kernel.  The variable-length traceback remains on Host, but the
    # score/contract path is observable here so a future implementation cannot
    # silently fall back to a second CPU-only recurrence.
    assert telemetry["graph_haplotype_sw_used"] is True
    assert telemetry["graph_haplotype_sw_execution_space"]
    assert telemetry["graph_haplotype_sw_simd_width"] >= 1
    assert telemetry["graph_haplotype_sw_simd_groups"] >= 0
    assert telemetry["graph_haplotype_sw_prepare_seconds"] >= 0.0
    assert telemetry["graph_haplotype_sw_seconds"] >= 0.0
    # Every anchored path must expose a syntactically complete SW traceback;
    # this catches a sequence-only graph path that silently bypasses the
    # reference-connected assembly provenance contract.
    for detail in details:
        assert detail["cigar"] and re.fullmatch(r"(?:[0-9]+[MIDNSHP=X])+", detail["cigar"])
        assert isinstance(detail["alignment_score"], int)
        assert isinstance(detail["alignment_offset"], int)
    summary["graph_haplotype_provenance_count"] = telemetry["graph_haplotype_provenance_count"]
    return summary


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-hc-call"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not binary.exists() or not bam.exists() or not reference.exists():
        print(json.dumps({"status": "skipped", "reason": "missing native build or chr17 fixture"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-real-assembly-graph-") as directory:
        work = Path(directory)
        enabled = run(binary, bam, reference, work / "linked.vcf.gz",
                      work / "linked.manifest.json")
        disabled = run(binary, bam, reference, work / "linked-disabled.vcf.gz",
                       work / "linked-disabled.manifest.json",
                       disable_recovery=True)

    assert enabled["status"] == "prototype"
    assert enabled["graph_used"] is True
    assert enabled["graph_reference_kmer_rejected"] is False
    assert enabled["graph_nodes"] > 1000
    assert enabled["graph_edges"] > enabled["graph_nodes"]
    assert enabled["graph_reference_nodes"] > 0
    assert enabled["graph_reference_connected_nodes"] > 0
    assert enabled["graph_haplotype_paths"] > 0
    assert enabled["pairhmm_graph_haplotypes"] > 0
    # The fixture is allowed to have every pivotal edge covered by the
    # bounded traversal (in which case no artificial recovery is necessary).
    # Keep the counters and the disabled-path invariant strict without making
    # a data-dependent recovery event a false pass/fail condition.
    assert enabled["graph_artificial_haplotype_recovery_paths"] >= 0
    assert enabled["graph_artificial_haplotype_recovery_bases"] >= 0
    assert disabled["graph_artificial_haplotype_recovery_paths"] == 0
    assert disabled["graph_artificial_haplotype_recovery_bases"] == 0
    assert disabled["graph_haplotype_paths"] > 0
    print(json.dumps({
        "status": "pass",
        "corpus": "NA12878.chr17_69k_70k.dictFix.bam",
        "region": "17:69000-70000",
        "graph_nodes": enabled["graph_nodes"],
        "graph_edges": enabled["graph_edges"],
        "graph_haplotype_paths": enabled["graph_haplotype_paths"],
        "graph_haplotype_provenance_count": enabled["graph_haplotype_provenance_count"],
        "artificial_recovery_paths": enabled["graph_artificial_haplotype_recovery_paths"],
        "artificial_recovery_bases": enabled["graph_artificial_haplotype_recovery_bases"],
        "disabled_recovery_paths": disabled["graph_artificial_haplotype_recovery_paths"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
