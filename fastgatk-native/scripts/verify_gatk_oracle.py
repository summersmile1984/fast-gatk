#!/usr/bin/env python3
"""Compare the native HaplotypeCaller core with a local GATK oracle.

This integration oracle requires exact, ordered VCF and gVCF data records on
the bundled chr17 fixture.  The ``GATKCommandLine`` header is deliberately
implementation-specific provenance and is checked separately from biological
output records.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import gzip
from pathlib import Path


def records(path: Path) -> set[tuple[str, int, str, str]]:
    result: set[tuple[str, int, str, str]] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        result.add((fields[0], int(fields[1]), fields[3], fields[4]))
    return result


def variant_fields(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def data_lines(path: Path) -> list[str]:
    """Return ordered VCF/gVCF data records exactly as written."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def call_fields(fields: list[list[str]]) -> dict[tuple[str, int, str, str], tuple[str, int]]:
    result: dict[tuple[str, int, str, str], tuple[str, int]] = {}
    for record in fields:
        format_keys = record[8].split(":")
        sample_values = dict(zip(format_keys, record[9].split(":")))
        result[(record[0], int(record[1]), record[3], record[4])] = (
            sample_values["GT"], int(sample_values["GQ"]))
    return result


def sample_values(fields: list[list[str]], key: tuple[str, int, str, str]) -> dict[str, str]:
    for record in fields:
        if (record[0], int(record[1]), record[3], record[4]) != key:
            continue
        return dict(zip(record[8].split(":"), record[9].split(":")))
    raise AssertionError(f"missing variant record: {key}")


def validate_gvcf(fields: list[list[str]], start: int, end: int) -> None:
    assert fields, "gVCF contains no records"
    previous_end = start - 1
    for record in fields:
        position = int(record[1])
        info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in record[7].split(";") if "=" in item}
        record_end = int(info.get("END", position))
        assert position > previous_end, f"overlapping gVCF records at {record[0]}:{position}"
        previous_end = record_end
        sample = dict(zip(record[8].split(":"), record[9].split(":")))
        assert 0 <= int(sample["GQ"]) <= 99
        pl = [int(value) for value in sample["PL"].split(",")]
        if record[4] == "<NON_REF>":
            assert len(pl) == 3 and pl[0] == 0
        elif "<NON_REF>" in record[4].split(","):
            concrete_alt_count = len(record[4].split(",")) - 1
            concrete_alleles = concrete_alt_count + 1
            assert len(pl) == (concrete_alleles + 1) * (concrete_alleles + 2) // 2
            assert all(value < 999 for value in pl)
    assert int(fields[0][1]) <= start
    assert previous_end >= end


def reference_confidence_pls(fields: list[list[str]]) -> dict[int, tuple[int, tuple[int, ...]]]:
    """Return POS -> (GQ, PL) for simple REF/<NON_REF> blocks."""
    result = {}
    for record in fields:
        if record[4] != "<NON_REF>":
            continue
        format_keys = record[8].split(":")
        sample_values = record[9].split(":")
        values = dict(zip(format_keys, sample_values))
        if "PL" not in values or "GQ" not in values:
            continue
        result[int(record[1])] = (int(values["GQ"]),
                                  tuple(int(value) for value in values["PL"].split(",")))
    return result


def reference_confidence_ends(fields: list[list[str]]) -> dict[int, int]:
    """Return one-based POS -> END for simple reference-confidence blocks."""
    result: dict[int, int] = {}
    for record in fields:
        if record[4] != "<NON_REF>":
            continue
        info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in record[7].split(";") if "=" in item}
        result[int(record[1])] = int(info.get("END", record[1]))
    return result


def reference_confidence_samples(
    fields: list[list[str]],
) -> dict[int, tuple[int, str, str, str, str, tuple[int, ...]]]:
    """Return the complete simple REF/<NON_REF> block sample contract."""
    result: dict[int, tuple[int, str, str, str, str, tuple[int, ...]]] = {}
    for record in fields:
        if record[4] != "<NON_REF>":
            continue
        info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in record[7].split(";") if "=" in item}
        values = dict(zip(record[8].split(":"), record[9].split(":")))
        result[int(record[1])] = (
            int(info.get("END", record[1])),
            values.get("GT", ""),
            values.get("DP", ""),
            values.get("GQ", ""),
            values.get("MIN_DP", ""),
            tuple(int(value) for value in values.get("PL", "").split(",")
                  if value != ""),
        )
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    region = "17:69000-69100"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or oracle fixtures; run build_native.sh first")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        gatk_gvcf = work / "gatk.g.vcf.gz"
        native_gvcf = work / "native.g.vcf.gz"
        manifest = work / "native.manifest.json"
        gvcf_manifest = work / "native.g.vcf.gz.manifest.json"
        gatk_command = [
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_vcf),
            "-L", region, "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ]
        subprocess.run(gatk_command, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        gatk_gvcf_command = [
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_gvcf),
            "-L", region, "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ]
        subprocess.run(gatk_gvcf_command, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        native_command = [
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_vcf), "--output-manifest", str(manifest),
            "--threads", "2", "--min-depth", "1", "--min-alt-support", "1",
            # GATK's generic Boolean argument accepts an explicit false;
            # it must leave the standard HC read-filter list intact.
            "--disable-tool-default-read-filters=false",
            # These AssemblyBasedCaller booleans must also retain their
            # default behavior when passed explicitly as false.
            "--dont-use-soft-clipped-bases=false",
            "--do-not-correct-overlapping-quality=false",
            "--dont-increase-kmer-sizes-for-cycles=false",
            "--recover-all-dangling-branches=false",
            "--error-correct-reads=false",
        ]
        subprocess.run(native_command, check=True, stdout=subprocess.DEVNULL)
        native_gvcf_command = [
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_gvcf), "-ERC", "GVCF", "--min-depth", "1",
            "--min-alt-support", "1", "--output-manifest", str(gvcf_manifest),
        ]
        subprocess.run(native_gvcf_command, check=True, stdout=subprocess.DEVNULL)
        gatk_calls = records(gatk_vcf)
        native_calls = records(native_vcf)
        gatk_vcf_fields = variant_fields(gatk_vcf)
        native_vcf_fields = variant_fields(native_vcf)
        gatk_vcf_data = data_lines(gatk_vcf)
        native_vcf_data = data_lines(native_vcf)
        native_call_fields = call_fields(native_vcf_fields)
        gatk_gvcf_fields = variant_fields(gatk_gvcf)
        native_gvcf_fields = variant_fields(native_gvcf)
        gatk_gvcf_data = data_lines(gatk_gvcf)
        native_gvcf_data = data_lines(native_gvcf)
        manifest_data = json.loads(manifest.read_text())
        gvcf_manifest_data = json.loads(gvcf_manifest.read_text())

    sentinel = ("17", 69067, "T", "G")
    assert sentinel in gatk_calls, f"GATK oracle lost expected sentinel call: {sentinel}"
    assert sentinel in native_calls, f"native pipeline lost expected sentinel call: {sentinel}"
    assert native_vcf_data == gatk_vcf_data, (
        f"VCF data-record mismatch: native={native_vcf_data} gatk={gatk_vcf_data}")
    gatk_sentinel_values = sample_values(gatk_vcf_fields, sentinel)
    native_sentinel_values = sample_values(native_vcf_fields, sentinel)
    for field in ("GT", "GQ", "DP", "AD"):
        assert native_sentinel_values.get(field) == gatk_sentinel_values.get(field), (
            f"sentinel {field} mismatch: native={native_sentinel_values.get(field)} "
            f"gatk={gatk_sentinel_values.get(field)}")
    # The call list consumed by telemetry and the VCF writer must agree on
    # every reference-backed PairHMM decision.  The manifest uses 0-based
    # coordinates/tid indices; VCF uses contig names and 1-based POS.
    for call in manifest_data["telemetry"]["calls"]:
        key = ("17", int(call["position"]) + 1, call["ref"], call["alt"])
        assert key in native_call_fields, f"manifest call missing from VCF: {key}"
        assert native_call_fields[key] == (call["genotype"], int(call["gq"])), (
            f"manifest/VCF genotype disagreement at {key}: "
            f"{native_call_fields[key]} != {(call['genotype'], call['gq'])}")
    # The manifest's process-wide bit-identical flag remains conservative:
    # the output GATKCommandLine provenance necessarily identifies the native
    # executable.  The fixture oracle above proves exact biological records.
    assert manifest_data["compatibility"]["bit_identical_to_gatk"] is False
    assert gvcf_manifest_data["compatibility"]["gvcf_standard_fields"] is True
    assert gvcf_manifest_data["compatibility"]["gvcf_candidate_standard_fields"] is True
    assert gvcf_manifest_data["compatibility"]["gvcf_candidate_rebuild_after_genotyping"] is True
    assert manifest_data["compatibility"]["pairhmm_haplotype_likelihoods"] is True
    assert manifest_data["compatibility"]["pairhmm_likelihood_normalization"] is True
    assert manifest_data["telemetry"]["pairhmm_normalization_execution_space"]
    assert manifest_data["compatibility"]["pairhmm_allele_marginalization"] is True
    assert manifest_data["telemetry"]["pairhmm_marginalization_execution_space"]
    assert manifest_data["compatibility"]["smith_waterman_scores"] is True
    assert manifest_data["compatibility"]["smith_waterman_traceback"] is True
    assert manifest_data["compatibility"]["kmer_graph"] is True
    assert manifest_data["telemetry"]["pairhmm_used"] is True
    assert manifest_data["telemetry"]["pairhmm_skip_reason"] == "executed"
    assert manifest_data["compatibility"]["pairhmm_fallback_reason"] == "executed"
    assert manifest_data["telemetry"]["pairhmm_pairs"] > 0
    assert manifest_data["telemetry"]["pairhmm_request_links"] >= \
        manifest_data["telemetry"]["pairhmm_pairs"]
    # HC now consumes the complete K-best EventMap population from the
    # assembled graph, as GATK's PairHMMLikelihoodCalculationEngine does.
    # A nonzero local Cartesian-combination block here would reintroduce
    # unsupported read×haplotype states; the graph-path telemetry below is
    # the positive proof that the assembled population was used.
    assert manifest_data["telemetry"]["pairhmm_haplotype_combination_blocks"] == 0
    assert "haplotype_combination_chunking" in manifest_data["compatibility"]
    assert manifest_data["telemetry"]["pairhmm_haplotypes"] >= 2
    assert manifest_data["telemetry"]["pairhmm_graph_haplotypes"] > 0
    assert manifest_data["compatibility"]["graph_haplotype_likelihoods"] is True
    assert manifest_data["compatibility"]["likelihood_based_haplotype_pruning"] is True
    assert manifest_data["telemetry"]["haplotype_pruning_log10"] == 20
    assert manifest_data["telemetry"]["pairhmm_graph_haplotypes_considered"] >= \
        manifest_data["telemetry"]["pairhmm_graph_haplotypes"] > 0
    # A cyclic graph retry can legitimately reduce the materialized graph set
    # before PairHMM (there is then nothing left for likelihood pruning).  A
    # retry can also be triggered by GATK's reference-k-mer rejection guard;
    # in that case the selected graph is acyclic but the larger selected k is
    # still the observable evidence that the bounded retry path ran.  Keep
    # the strict pruning oracle for the single-k attempt while accepting both
    # GATK retry causes.
    if manifest_data["telemetry"]["graph_kmer_iterations"] == 1:
        assert manifest_data["telemetry"]["pairhmm_graph_haplotypes_pruned"] > 0
    else:
        assert manifest_data["telemetry"]["graph_kmer_iterations"] > 1
        assert (
            manifest_data["telemetry"]["graph_has_non_reference_cycles"] is True
            or manifest_data["telemetry"]["graph_kmer_size_selected"] >
            manifest_data["telemetry"]["graph_kmer_size"]
        ), "graph retry has neither cycle nor selected larger k evidence"
    assert manifest_data["telemetry"]["pairhmm_haplotypes_pruned"] >= \
        manifest_data["telemetry"]["pairhmm_graph_haplotypes_pruned"]
    assert manifest_data["compatibility"]["pairhmm_haplotype_set_marginalization"] is True
    assert manifest_data["compatibility"]["assembly_region_pairhmm_partitioning"] is True
    assert manifest_data["telemetry"]["pairhmm_assembly_region_groups"] >= 1
    assert manifest_data["telemetry"]["pairhmm_unassigned_candidates"] == 0
    assert manifest_data["telemetry"]["pairhmm_assembly_region_partitioned"] is True
    assert manifest_data["telemetry"]["assembly_region_union_count"] >= 1
    assert manifest_data["telemetry"]["assembly_region_union_merges"] >= 0
    assert manifest_data["telemetry"]["assembly_cross_region_rescued"] >= 0
    # AssemblyRegion hard-clipping is a Host semantic boundary before the
    # Kokkos PairHMM call. The chr17 sentinel has one request extending beyond
    # its local window, so the clipping path must be observable in telemetry.
    assert manifest_data["telemetry"]["pairhmm_reads_clipped"] >= 1
    assert manifest_data["telemetry"]["pairhmm_reads_dropped_after_clipping"] == 0
    assert manifest_data["telemetry"]["pairhmm_reads_disqualified"] >= 1
    assert manifest_data["compatibility"]["assembly_region_candidate_partitioning"] is True
    assert manifest_data["telemetry"]["assembly_cross_region_candidates"] == 0
    assert manifest_data["telemetry"]["assembly_unassigned_candidates"] == 0
    assert manifest_data["telemetry"]["sw_used"] is True
    assert manifest_data["telemetry"]["sw_pairs"] == manifest_data["telemetry"]["pairhmm_pairs"]
    assert manifest_data["telemetry"]["sw_traceback_pairs"] == manifest_data["telemetry"]["sw_pairs"]
    # The SW request matrix is exactly the read×haplotype matrix staged for
    # PairHMM.  Keep an explicit provenance/uncertainty contract in the
    # manifest so a caller does not mistake a generic assembly traceback for
    # read-likelihood evidence.
    assert manifest_data["compatibility"]["read_haplotype_cigar_provenance"] is True
    assert manifest_data["telemetry"]["read_haplotype_cigar_pairs"] == \
        manifest_data["telemetry"]["sw_pairs"]
    assert manifest_data["telemetry"]["read_haplotype_cigar_valid_pairs"] == \
        manifest_data["telemetry"]["read_haplotype_cigar_pairs"]
    assert manifest_data["telemetry"]["read_haplotype_cigar_filtered_pairs"] <= \
        manifest_data["telemetry"]["read_haplotype_cigar_pairs"]
    assert manifest_data["telemetry"]["read_haplotype_softclip_filtered_pairs"] <= \
        manifest_data["telemetry"]["read_haplotype_cigar_filtered_pairs"]
    assert manifest_data["telemetry"]["pairhmm_flow_reads_clipped"] <= \
        manifest_data["telemetry"]["pairhmm_flow_reads"]
    assert manifest_data["telemetry"]["read_haplotype_informative_reads"] >= \
        manifest_data["telemetry"]["read_haplotype_uncertain_reads"]
    assert manifest_data["telemetry"]["read_haplotype_informative_reads"] > 0
    assert manifest_data["telemetry"]["read_haplotype_margin_sum"] >= 0.0
    assert manifest_data["telemetry"]["read_haplotype_uncertainty_signature"] != 0
    assert manifest_data["telemetry"]["graph_used"] is True
    assert manifest_data["telemetry"]["graph_nodes"] > 0
    assert manifest_data["compatibility"]["reference_connected_graph"] is True
    assert manifest_data["compatibility"]["activity_profile"] is True
    assert manifest_data["compatibility"]["activity_reference_projection"] is True
    assert manifest_data["telemetry"]["activity_reference_aware"] is True
    assert manifest_data["compatibility"]["reference_confidence_kokkos"] is True
    assert manifest_data["compatibility"]["gvcf_non_ref_likelihood_model"] in {
        "median-qualified-concrete-alleles", "concrete-envelope"
    }
    assert manifest_data["telemetry"]["reference_confidence_execution_space"]
    assert manifest_data["telemetry"]["reference_confidence_seconds"] > 0.0
    assert manifest_data["compatibility"]["read_filter_order"] is True
    assert manifest_data["compatibility"]["read_filter_applied"] is True
    assert manifest_data["compatibility"]["read_filter_gatk_defaults"] is True
    assert manifest_data["telemetry"]["read_filter_min_mapping_quality"] == 20
    assert manifest_data["telemetry"]["read_filter_exclude_duplicates"] is True
    # HaplotypeCallerEngine.makeStandardHCReadFilters() does not filter
    # supplementary alignments by default.
    assert manifest_data["telemetry"]["read_filter_exclude_supplementary"] is False
    assert manifest_data["telemetry"]["read_filter_require_no_n_cigar"] is True
    assert manifest_data["telemetry"]["read_filter_require_read_group"] is True
    assert manifest_data["compatibility"]["genotype_priors"] is True
    assert manifest_data["telemetry"]["graph_reference_connected_nodes"] > 0
    assert manifest_data["telemetry"]["indel_size_to_eliminate_in_ref_model"] == 10
    assert manifest_data["compatibility"]["haplotype_traversal"] is True
    assert manifest_data["telemetry"]["graph_haplotype_paths"] > 0
    assert manifest_data["telemetry"]["graph_haplotype_sequence_count"] > 0
    assert manifest_data["telemetry"]["graph_snp_candidates"] >= 0
    assert "graph_snp_candidate_extraction" in manifest_data["compatibility"]
    assert manifest_data["telemetry"]["assembly_regions"] > 0
    validate_gvcf(gatk_gvcf_fields, 69000, 69100)
    validate_gvcf(native_gvcf_fields, 69000, 69100)
    candidate_gvcf = next(field for field in native_gvcf_fields
                          if field[0] == "17" and int(field[1]) == 69067)
    gatk_candidate_gvcf = next(field for field in gatk_gvcf_fields
                               if field[0] == "17" and int(field[1]) == 69067)
    gatk_candidate_values = dict(zip(gatk_candidate_gvcf[8].split(":"),
                                     gatk_candidate_gvcf[9].split(":")))
    native_candidate_values = dict(zip(candidate_gvcf[8].split(":"),
                                       candidate_gvcf[9].split(":")))
    for field in ("GT", "GQ", "AD"):
        assert native_candidate_values.get(field) == gatk_candidate_values.get(field), (
            f"candidate-site gVCF {field} mismatch: native={native_candidate_values.get(field)} "
            f"gatk={gatk_candidate_values.get(field)}")
    candidate_pl = [int(value) for value in candidate_gvcf[9].split(":")[
        candidate_gvcf[8].split(":").index("PL")].split(",")]
    gatk_candidate_pl = [int(value) for value in gatk_candidate_gvcf[9].split(":")[
        gatk_candidate_gvcf[8].split(":").index("PL")].split(",")]
    assert len(candidate_pl) == 6 and all(value < 999 for value in candidate_pl)
    assert manifest_data["compatibility"]["gvcf_non_ref_likelihood_model"] == \
        "median-qualified-concrete-alleles"
    # Candidate-site rows use GATK's normal VariantContext contract (site
    # QUAL/FILTER, INFO ordering, raw MQ/DP accumulator and FORMAT/SB) rather
    # than the simpler reference-block shape.  Keep a strict text oracle for
    # the pinned sentinel now that the native writer preserves uncapped PL.
    assert candidate_gvcf == gatk_candidate_gvcf, (
        f"candidate-site gVCF row mismatch: native={candidate_gvcf} "
        f"gatk={gatk_candidate_gvcf}")
    gatk_rcm = reference_confidence_pls(gatk_gvcf_fields)
    native_rcm = reference_confidence_pls(native_gvcf_fields)
    gatk_rcm_ends = reference_confidence_ends(gatk_gvcf_fields)
    native_rcm_ends = reference_confidence_ends(native_gvcf_fields)
    # The block writer now follows HomRefBlock exactly: END/GT/DP/GQ/MIN_DP/PL
    # must agree for every simple reference-confidence record, not only the
    # selected sentinel starts below.  Candidate-site gVCF records are checked
    # separately because they carry concrete ALT alleles.
    assert reference_confidence_samples(gatk_gvcf_fields) == \
        reference_confidence_samples(native_gvcf_fields), (
            "simple reference-confidence block mismatch")
    # These loci are outside the candidate span and exercise the RCM SNP plus
    # indel-informative model. The selected block-start triples are required to
    # match GATK 4.6.2.0; the final terminal-soft-clip boundary is tracked as
    # a separate compatibility gap rather than being silently treated as
    # bit-identical.
    rcm_positions = (69000, 69001, 69031, 69032, 69038, 69054, 69068, 69071, 69072, 69080, 69093)
    rcm_exact = 0
    for position in rcm_positions:
        if position in gatk_rcm and position in native_rcm:
            assert native_rcm[position] == gatk_rcm[position], (
                f"RCM mismatch at {position}: {native_rcm[position]} != {gatk_rcm[position]}")
            rcm_exact += 1
    assert rcm_exact >= 10
    # Keep the terminal soft-clip boundary explicit.  The Host CIGAR boundary
    # now matches GATK on this fixture; the boolean is intentionally scoped to
    # the two terminal records and does not claim whole-file bit identity.
    expected_gatk_terminal = {69080: 69092, 69093: 69100}
    assert {position: gatk_rcm_ends[position] for position in expected_gatk_terminal
            if position in gatk_rcm_ends} == expected_gatk_terminal
    terminal_boundary_exact = all(
        native_rcm_ends.get(position) == end
        for position, end in expected_gatk_terminal.items())
    assert terminal_boundary_exact, "terminal reference-confidence boundary mismatch"
    gatk_gvcf_variants = {
        (field[0], int(field[1]), field[3], field[4].split(",")[0]) for field in gatk_gvcf_fields
        if field[4] != "<NON_REF>"
    }
    native_gvcf_variants = {
        (field[0], int(field[1]), field[3], field[4].split(",")[0]) for field in native_gvcf_fields
        if field[4] != "<NON_REF>"
    }
    assert sentinel in gatk_gvcf_variants and sentinel in native_gvcf_variants
    # Candidate-site presence is now exact on the bundled fixture: uncalled
    # SNP/MNP candidates are folded back into RCM blocks, while called sites
    # (and any indel span that must reserve reference bases) remain explicit.
    assert gatk_gvcf_variants == native_gvcf_variants, (
        f"GATK/native gVCF candidate-site mismatch: "
        f"gatk-only={sorted(gatk_gvcf_variants - native_gvcf_variants)}, "
        f"native-only={sorted(native_gvcf_variants - gatk_gvcf_variants)}")
    assert native_gvcf_data == gatk_gvcf_data, (
        f"gVCF data-record mismatch: native={native_gvcf_data} gatk={gatk_gvcf_data}")
    intersection = gatk_calls & native_calls
    union = gatk_calls | native_calls
    print(json.dumps({
        "status": "pass",
        "gatk_calls": len(gatk_calls),
        "native_calls": len(native_calls),
        "intersection": len(intersection),
        "jaccard": len(intersection) / len(union) if union else 1.0,
        "sentinel_call": sentinel,
        "gatk_gvcf_records": len(gatk_gvcf_fields),
        "native_gvcf_records": len(native_gvcf_fields),
        "rcm_block_exact": rcm_exact,
        "rcm_terminal_boundary_exact": terminal_boundary_exact,
        "rcm_terminal_boundary_gatk": {
            str(position): gatk_rcm_ends.get(position) for position in (69080, 69093)},
        "rcm_terminal_boundary_native": {
            str(position): native_rcm_ends.get(position) for position in (69080, 69093)},
        "gvcf_sentinel": sentinel,
        "candidate_gvcf_gt_gq_ad_exact": True,
        "candidate_gvcf_pl_bit_identical": candidate_pl == gatk_candidate_pl,
        "candidate_gvcf_record_text_exact": True,
        "vcf_data_records_bit_identical": True,
        "gvcf_data_records_bit_identical": True,
        "provenance_header_bit_identical": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
