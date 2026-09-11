#!/usr/bin/env python3
"""Run the pinned Mutect2 oracle with strict biological-output comparison.

The bundled broad fixture requires exact ordered VCF records, including TLOD,
and exact schema metadata.  ``GATKCommandLine`` is execution provenance: a
native executable must not impersonate a Java invocation, so it is the only
excluded header line in the schema comparison.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import tarfile
import tempfile
from pathlib import Path
import oracle_guard


def run(
    command: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False, cwd=cwd, env=env)


_GATK_DOT_EDGE = re.compile(r"(?P<source>\S+)\s+->\s+(?P<target>\S+)\s+\[label=")
_GATK_DOT_KMER = re.compile(r"-sequenceGraph\.(?P<kmer>\d+)\.0\.2\.cleaned_readthreading_graph\.dot$")
_NATIVE_KMER_ATTEMPT = re.compile(
    r"^\[FASTGATK_KMER_ATTEMPT\] k=(?P<kmer>\d+) "
    r"reference_rejected=(?P<rejected>[01]) .* nodes=(?P<nodes>\d+) paths=(?P<paths>\d+)$",
    re.MULTILINE,
)


def gatk_cleaned_graph_topologies(work: Path) -> list[tuple[int, int, int, int]]:
    """Return (k-mer, vertices, edges, source-to-sink paths) from GATK DOT.

    `--debug-graph-transformations` emits each graph edge twice: once with
    its multiplicity label and once with a styling attribute.  Retaining only
    the labelled edge is therefore the stable graph topology, independent of
    DOT colours or Java identity-hash vertex IDs.
    """
    topologies: list[tuple[int, int, int, int]] = []
    for path in sorted(work.glob("*-sequenceGraph.*.0.2.cleaned_readthreading_graph.dot")):
        match = _GATK_DOT_KMER.search(path.name)
        assert match, path
        edges = [
            (edge.group("source"), edge.group("target"))
            for edge in _GATK_DOT_EDGE.finditer(path.read_text(encoding="utf-8"))
        ]
        assert edges, path
        vertices = {vertex for edge in edges for vertex in edge}
        outgoing: dict[str, list[str]] = {vertex: [] for vertex in vertices}
        incoming: dict[str, int] = {vertex: 0 for vertex in vertices}
        for source, target in edges:
            outgoing[source].append(target)
            incoming[target] += 1
        sources = [vertex for vertex, degree in incoming.items() if degree == 0]
        sinks = {vertex for vertex, targets in outgoing.items() if not targets}
        assert len(sources) == 1 and sinks, (path, sources, sinks)

        visiting: set[str] = set()
        paths_from: dict[str, int] = {}

        def count_paths(vertex: str) -> int:
            if vertex in paths_from:
                return paths_from[vertex]
            assert vertex not in visiting, f"unexpected cycle in cleaned GATK graph: {path}"
            if vertex in sinks:
                paths_from[vertex] = 1
                return 1
            visiting.add(vertex)
            paths_from[vertex] = sum(count_paths(target) for target in outgoing[vertex])
            visiting.remove(vertex)
            return paths_from[vertex]

        topologies.append((int(match.group("kmer")), len(vertices), len(edges), count_paths(sources[0])))
    return topologies


def native_cleaned_graph_topologies(stderr: str) -> list[tuple[int, int, int]]:
    """Extract the chosen, non-rejected native graph attempt per region."""
    return [
        (int(match.group("kmer")), int(match.group("nodes")), int(match.group("paths")))
        for match in _NATIVE_KMER_ATTEMPT.finditer(stderr)
        if match.group("rejected") == "0"
    ]


def gatk_pairhmm_haplotypes(path: Path) -> set[str]:
    return {
        line.split(maxsplit=1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    }


def native_pairhmm_haplotypes(stderr: str) -> set[str]:
    prefix = "[FASTGATK_PAIRHMM_REQUEST] "
    return {
        line[len(prefix):].split("\t", maxsplit=1)[0]
        for line in stderr.splitlines()
        if line.startswith(prefix)
    }


def read_vcf(path: Path) -> tuple[list[str], dict[tuple[str, int, str, str], dict[str, object]]]:
    records: dict[tuple[str, int, str, str], dict[str, object]] = {}
    samples: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            sample_values = fields[9].split(":") if len(fields) > 9 else []
            info: dict[str, str] = {}
            if fields[7] not in {"", "."}:
                for token in fields[7].split(";"):
                    if "=" in token:
                        name, value = token.split("=", 1)
                        info[name] = value
            records[key] = {
                "qual": fields[5],
                "filter": fields[6],
                "info": info,
                "format": dict(zip(format_keys, sample_values)),
            }
    return samples, records


def vcf_data_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def vcf_schema_lines(path: Path) -> list[str]:
    """Return all VCF schema lines except execution-specific provenance."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line.startswith("#") and not line.startswith("##GATKCommandLine=")]


def read_assembly_region_igv(path: Path) -> list[tuple[str, int, int, str, str]]:
    """Read GATK AssemblyRegionWalker's IGV rows without normalizing values."""
    rows: list[tuple[str, int, int, str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("Chromosome\t"):
            continue
        fields = line.split("\t")
        assert len(fields) == 5, fields
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]))
    return rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_mutect2_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-gatk-oracle-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        gatk_f1r2 = work / "gatk.f1r2.tar.gz"
        native_f1r2 = work / "native.f1r2.tar.gz"
        gatk_high_callable = work / "gatk-high-callable.vcf.gz"
        native_high_callable = work / "native-high-callable.vcf.gz"
        native_stats = work / "native.stats.json"
        native_manifest = work / "native.manifest.json"
        gatk_sites_only = work / "gatk-sites-only.vcf.gz"
        native_sites_only = work / "native-sites-only.vcf.gz"
        native_sites_only_stats = work / "native-sites-only.stats.json"
        native_sites_only_manifest = work / "native-sites-only.manifest.json"
        gatk_pairhmm_results = work / "gatk.pairhmm.txt"
        gatk_assembly_regions = work / "gatk.assembly-regions.igv"
        native_assembly_regions = work / "native.assembly-regions.igv"
        native_high_depth_assembly_regions = work / "native.high-depth.assembly-regions.igv"
        gatk_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(gatk_output), "--f1r2-tar-gz", str(gatk_f1r2),
            "--pair-hmm-results-file", str(gatk_pairhmm_results),
            "--assembly-region-out", str(gatk_assembly_regions),
            "--debug-graph-transformations", "true",
        ], cwd=work)
        assert gatk_run.returncode == 0, gatk_run.stderr
        assert gatk_pairhmm_results.stat().st_size > 0
        gatk_pairhmm_rows = sum(
            1 for line in gatk_pairhmm_results.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        )
        gatk_pairhmm_haplotype_set = gatk_pairhmm_haplotypes(gatk_pairhmm_results)
        native_debug_env = os.environ.copy()
        native_debug_env["FASTGATK_DEBUG_KMER_ATTEMPTS"] = "1"
        native_debug_env["FASTGATK_DEBUG_PAIRHMM_REQUESTS"] = "1"
        native_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(native_output), "--f1r2-tar-gz", str(native_f1r2),
            "--stats", str(native_stats), "--output-manifest", str(native_manifest),
            "--assembly-region-out", str(native_assembly_regions),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ], env=native_debug_env)
        assert native_run.returncode == 0, native_run.stderr
        # Mutect2's default companion is a two-column stats table, not a
        # native-only diagnostic document. FilterMutectCalls derives this
        # exact sibling path when --stats is omitted, so retain the byte-level
        # GATK sidecar contract independently of the JSON diagnostic supplied
        # explicitly above.
        gatk_stats_table = Path(f"{gatk_output}.stats")
        native_stats_table = Path(f"{native_output}.stats")
        assert gatk_stats_table.is_file() and native_stats_table.is_file()
        assert native_stats_table.read_text(encoding="utf-8") == \
            gatk_stats_table.read_text(encoding="utf-8")
        # Exercise the source's second callable branch explicitly. At an
        # intentionally unreachable depth it reports only raw isActive()
        # loci, rather than treating the band-pass-expanded AssemblyRegions
        # as callable. This protects both the option and the low-depth rule.
        gatk_high_callable_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "--callable-depth", "1000", "-O", str(gatk_high_callable),
        ])
        assert gatk_high_callable_run.returncode == 0, gatk_high_callable_run.stderr
        native_high_callable_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "--callable-depth", "1000", "-O", str(native_high_callable),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        assert native_high_callable_run.returncode == 0, native_high_callable_run.stderr
        assert Path(f"{native_high_callable}.stats").read_text(encoding="utf-8") == \
            Path(f"{gatk_high_callable}.stats").read_text(encoding="utf-8")
        # Compare the source graph after the exact phase where GATK has run
        # AdaptiveChainPruner and dangling-end recovery.  Final VCF equality
        # alone cannot catch a retained error chain or a missing rejoin when
        # it happens not to alter one of the six emitted sites.
        gatk_graph_topologies = gatk_cleaned_graph_topologies(work)
        native_graph_topologies = native_cleaned_graph_topologies(native_run.stderr)
        assert len(gatk_graph_topologies) == 5, gatk_graph_topologies
        assert native_graph_topologies == [
            (kmer, vertices, paths)
            for kmer, vertices, _edges, paths in gatk_graph_topologies
        ], {
            "gatk_cleaned_graphs": gatk_graph_topologies,
            "native_cleaned_graphs": native_graph_topologies,
            "native_debug": native_run.stderr,
        }
        native_pairhmm_haplotype_set = native_pairhmm_haplotypes(native_run.stderr)
        assert native_pairhmm_haplotype_set == gatk_pairhmm_haplotype_set, {
            "missing": sorted(gatk_pairhmm_haplotype_set - native_pairhmm_haplotype_set),
            "extra": sorted(native_pairhmm_haplotype_set - gatk_pairhmm_haplotype_set),
        }
        assert gatk_assembly_regions.is_file() and native_assembly_regions.is_file()
        gatk_assembly_rows = read_assembly_region_igv(gatk_assembly_regions)
        native_assembly_rows = read_assembly_region_igv(native_assembly_regions)
        # This is a strict AssemblyRegion ownership oracle, independent of the
        # release-specific somatic posterior/TLOD path.  It catches dropped
        # low-depth activity states, incorrect band-pass boundaries, and
        # off-by-one/padding errors before any PairHMM tolerance is considered.
        assert native_assembly_rows == gatk_assembly_rows, {
            "gatk_assembly_rows": gatk_assembly_rows,
            "native_assembly_rows": native_assembly_rows,
        }
        assert native_assembly_rows
        # `--min-depth` is a native calling/HC control, not a Mutect2
        # activity predicate.  Re-run the same fixture with an intentionally
        # impossible calling depth and require the Java AssemblyRegion profile
        # to remain unchanged.  This guards the quality-aware Kokkos branch
        # against reintroducing the old min-depth gate.
        native_high_depth_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(work / "native-high-depth.vcf.gz"),
            "--assembly-region-out", str(native_high_depth_assembly_regions),
            "--min-depth", "100", "--min-alt-support", "1",
        ])
        assert native_high_depth_run.returncode == 0, native_high_depth_run.stderr
        assert read_assembly_region_igv(native_high_depth_assembly_regions) == gatk_assembly_rows
        # The pinned fixture exercises overlapping AssemblyRegion halos.  The
        # native path must report the same stable clipping contract used by
        # GATK's AssemblyBasedCallerUtils: clipped fragments shorter than ten
        # bases are not sent to PairHMM.  Also require active-core ownership
        # to produce one bounded group per activity island; without this,
        # candidate evidence is duplicated across overlapping halos.
        native_stats_values = json.loads(native_stats.read_text(encoding="utf-8"))
        native_manifest_values = json.loads(native_manifest.read_text(encoding="utf-8"))
        native_telemetry = native_manifest_values["telemetry"]
        # The pinned primary oracle is fully reference-backed and its exact
        # GATK AssemblyRegion/PairHMM handoff was verified above.  Do not let
        # a future implementation obtain the same VCF through the native-only
        # pileup recovery path or an unassigned contig-wide candidate group.
        assert native_stats_values["pairhmm_skip_reason"] == "executed"
        assert native_stats_values["somatic_pileup_fallback_candidates"] == 0
        assert native_stats_values["pairhmm_unassigned_candidates"] == 0
        assert native_stats_values["assembly_unassigned_candidates"] == 0
        assert native_stats_values["pairhmm_minimum_read_length_after_trimming"] == 10
        assert native_stats_values["pairhmm_assembly_region_groups"] == 4
        assert native_stats_values["pairhmm_assembly_region_partitioned"] is True
        # This fixture is also the assembly/PairHMM boundary audit.  Keep the
        # native request matrix exactly equal to GATK's debug rows: a change
        # in TLOD tolerances must not conceal a clipped-read or haplotype-phase
        # difference at the AssemblyRegion handoff.
        assert gatk_pairhmm_rows == 706, gatk_pairhmm_rows
        assert native_stats_values["pairhmm_pairs"] == gatk_pairhmm_rows
        assert native_stats_values["pairhmm_haplotypes"] == 18
        # PairHMMLikelihoodCalculationEngine consumes the assembled
        # AssemblyResult haplotypes.  The exact request-row and haplotype-set
        # checks above are the GATK oracle; this native provenance counter
        # must consequently show that at least one graph path reached the
        # Kokkos PairHMM, rather than forcing an equivalent-looking local
        # Cartesian reconstruction to make the counter zero.
        assert native_stats_values["pairhmm_graph_haplotypes"] > 0
        # ReadThreadingAssembler retries the default [10,25] list in +10
        # steps for repeat-rich regions.  This fixture requires a long k-mer;
        # retain that source-level boundary so a packed-64-bit 31-mer ceiling
        # cannot silently return.
        assert native_telemetry["graph_kmer_size_selected"] > 31, native_telemetry
        assert native_telemetry["graph_kmer_iterations"] >= 4, native_telemetry
        gatk_sites_only_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(gatk_sites_only), "--sites-only-vcf-output", "true",
        ])
        assert gatk_sites_only_run.returncode == 0, gatk_sites_only_run.stderr
        native_sites_only_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(native_sites_only), "--sites-only-vcf-output", "true",
            "--stats", str(native_sites_only_stats),
            "--output-manifest", str(native_sites_only_manifest),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        assert native_sites_only_run.returncode == 0, native_sites_only_run.stderr
        assert Path(f"{native_sites_only}.stats").read_text(encoding="utf-8") == \
            Path(f"{gatk_sites_only}.stats").read_text(encoding="utf-8")
        assert native_output.stat().st_size > 0 and Path(f"{native_output}.tbi").is_file()
        assert gatk_output.stat().st_size > 0 and Path(f"{gatk_output}.tbi").is_file()
        assert native_f1r2.stat().st_size > 0 and gatk_f1r2.stat().st_size > 0
        with tarfile.open(native_f1r2, "r:gz") as archive:
            native_members = sorted(member.name for member in archive.getmembers())
        with tarfile.open(gatk_f1r2, "r:gz") as archive:
            gatk_members = sorted(member.name for member in archive.getmembers())
        gatk_samples, gatk_records = read_vcf(gatk_output)
        native_samples, native_records = read_vcf(native_output)
        gatk_sites_only_samples, gatk_sites_only_records = read_vcf(gatk_sites_only)
        native_sites_only_samples, native_sites_only_records = read_vcf(native_sites_only)
        gatk_data = vcf_data_lines(gatk_output)
        native_data = vcf_data_lines(native_output)
        gatk_sites_only_data = vcf_data_lines(gatk_sites_only)
        native_sites_only_data = vcf_data_lines(native_sites_only)
        gatk_schema = vcf_schema_lines(gatk_output)
        native_schema = vcf_schema_lines(native_output)
        assert gatk_sites_only_samples == [] and native_sites_only_samples == []
        assert set(native_sites_only_records) == set(native_records)
        assert all(len(record) == 8 for record in gzip.open(native_sites_only, "rt", encoding="utf-8")
                   if record and not record.startswith("#") for record in [record.rstrip("\n").split("\t")])
        assert all(len(record) == 8 for record in gzip.open(gatk_sites_only, "rt", encoding="utf-8")
                   if record and not record.startswith("#") for record in [record.rstrip("\n").split("\t")])
        shared = set(gatk_records) & set(native_records)
        comparable = ["GT", "DP", "AD"]
        exact = {field: 0 for field in comparable}
        tlod_deltas: list[float] = []
        af_deltas: list[float] = []
        for key in shared:
            gatk_format = gatk_records[key]["format"]
            native_format = native_records[key]["format"]
            for field in comparable:
                if gatk_format.get(field) == native_format.get(field):
                    exact[field] += 1
            gatk_info = gatk_records[key]["info"]
            native_info = native_records[key]["info"]
            if "TLOD" in gatk_info and "TLOD" in native_info:
                tlod_deltas.append(abs(float(gatk_info["TLOD"]) - float(native_info["TLOD"])))
            if "AF" in gatk_info and "AF" in native_info:
                af_deltas.append(abs(float(gatk_info["AF"]) - float(native_info["AF"])))
        assert gatk_samples == ["NA12878"] and native_samples == ["NA12878"]
        # The pinned broad fixture is now a call-set oracle, not merely an
        # overlap smoke test.  Keep the release-specific posterior/GT
        # differences measurable below, but do not allow the low-support
        # candidate boundary to regress to silent site loss.
        assert set(gatk_records) == set(native_records), {
            "gatk_sites": len(gatk_records), "native_sites": len(native_records),
            "missing_from_native": sorted(set(gatk_records) - set(native_records)),
            "extra_in_native": sorted(set(native_records) - set(gatk_records)),
        }
        assert native_data == gatk_data, {
            "native": native_data,
            "gatk": gatk_data,
        }
        assert native_sites_only_data == gatk_sites_only_data, {
            "native_sites_only": native_sites_only_data,
            "gatk_sites_only": gatk_sites_only_data,
        }
        assert native_schema == gatk_schema, {
            "native_schema": native_schema,
            "gatk_schema": gatk_schema,
        }
        # These fields are part of the stable FORMAT contract for the pinned
        # fixture.  A regression here usually means annotation data was
        # updated in one in-memory view but not the view serialized to VCF.
        assert exact["GT"] == len(shared), exact
        assert exact["DP"] == len(shared), exact
        assert exact["AD"] == len(shared), exact
        assert all(69000 <= key[1] <= 70000 for key in native_records), native_records
        assert native_members == ["native.f1r2.ref_histogram", "native.f1r2.alt_histogram", "native.f1r2.alt_table"] or len(native_members) == 3
        assert len(gatk_members) == 3
        # `native_data == gatk_data` above has already made the serialized
        # TLOD contract exact. Keep the numeric reduction explicit as a
        # diagnostic against a formatting-only coincidence.
        if tlod_deltas:
            assert max(tlod_deltas) == 0.0, tlod_deltas
            assert sum(tlod_deltas) == 0.0, tlod_deltas
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "gatk_sites": len(gatk_records),
            "native_sites": len(native_records),
            "shared_sites": len(shared),
            "site_set_exact": set(gatk_records) == set(native_records),
            "site_recall": round(len(shared) / len(gatk_records), 6),
            "site_precision": round(len(shared) / len(native_records), 6),
            "comparable_fields": exact,
            "shared_tlod_records": len(tlod_deltas),
            "shared_tlod_max_abs_delta": round(max(tlod_deltas), 6) if tlod_deltas else None,
            "shared_tlod_mean_abs_delta": round(sum(tlod_deltas) / len(tlod_deltas), 6) if tlod_deltas else None,
            "shared_af_records": len(af_deltas),
            "shared_af_max_abs_delta": round(max(af_deltas), 8) if af_deltas else None,
            "gatk_pairhmm_debug_rows": gatk_pairhmm_rows,
            "pairhmm_haplotype_set_exact": True,
            "pairhmm_haplotype_count": len(gatk_pairhmm_haplotype_set),
            "assembly_region_rows_exact": True,
            "assembly_region_row_count": len(native_assembly_rows),
            "activity_independent_of_min_depth": True,
            "native_pairhmm_request_pairs": native_stats_values["pairhmm_pairs"],
            "native_pairhmm_haplotypes": native_stats_values["pairhmm_haplotypes"],
            "native_graph_kmer_size_selected": native_telemetry["graph_kmer_size_selected"],
            "native_graph_kmer_iterations": native_telemetry["graph_kmer_iterations"],
            "assembly_graph_topology_exact": True,
            "assembly_graph_topologies": [
                {"kmer": kmer, "nodes": nodes, "edges": edges, "paths": paths}
                for kmer, nodes, edges, paths in gatk_graph_topologies
            ],
            "pairhmm_request_set_delta": native_stats_values["pairhmm_pairs"] - gatk_pairhmm_rows,
            "pairhmm_recurrence_oracle": "fastgatk-pairhmm-results-gatk-oracle",
            "sample_contract_exact": True,
            "f1r2_archive_members_exact_count": len(gatk_members) == len(native_members),
            "vcf_data_records_bit_identical": True,
            "sites_only_data_records_bit_identical": True,
            "schema_header_bit_identical": True,
            "mutect2_stats_table_bit_identical": True,
            "mutect2_high_callable_stats_table_bit_identical": True,
            "execution_provenance_header_bit_identical": False,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
