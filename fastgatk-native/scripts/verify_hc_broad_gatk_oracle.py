#!/usr/bin/env python3
"""Broad HaplotypeCaller oracle over the pinned chr17 69k-70k fixture.

The narrow oracle protects the sentinel/gVCF boundary.  This check protects
the assembly evidence boundary: soft-clipped homopolymer noise and a
pair-only artifact must not be emitted as variants, while the three GATK
calls in the 1 kb window must remain present.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def variant_records(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n").split("\t") for line in stream
                if line and not line.startswith("#")]


def vcf_payload_bytes(path: Path) -> bytes:
    """Compare the VCF stream, not BGZF container timestamps/block layout."""
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as stream:
            return stream.read()
    return path.read_bytes()


def header_definitions(path: Path) -> set[str]:
    """Return schema-bearing VCF header lines, excluding provenance/dictionary."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return {
            line.rstrip("\n")
            for line in stream
            if line.startswith(("##FILTER=", "##FORMAT=", "##INFO=", "##ALT="))
        }


def stable_header(path: Path) -> list[str]:
    """Return ordered header lines except GATK's dynamic command provenance."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [
            line.rstrip("\n")
            for line in stream
            if line.startswith("#") and not line.startswith("##GATKCommandLine=")
        ]


def key(record: list[str]) -> tuple[str, int, str, str]:
    return record[0], int(record[1]), record[3], record[4]


def sample(record: list[str]) -> dict[str, str]:
    return dict(zip(record[8].split(":"), record[9].split(":")))


def info(record: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in record[7].split(";"):
        if "=" in item:
            name, value = item.split("=", 1)
            values[name] = value
    return values


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    cram = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.cram"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    region = "17:69000-70000"
    if not native.exists() or not bam.exists() or not cram.exists() or not reference.exists():
        raise SystemExit("missing native build or broad HC fixture")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-broad-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        gatk_float32_vcf = work / "gatk.float32.vcf"
        native_float32_vcf = work / "native.float32.vcf"
        gatk_cram_vcf = work / "gatk.cram.vcf"
        native_cram_vcf = work / "native.cram.vcf"
        gatk_gvcf = work / "gatk.g.vcf.gz"
        native_gvcf = work / "native.g.vcf.gz"
        gatk_float32_gvcf = work / "gatk.float32.g.vcf.gz"
        native_float32_gvcf = work / "native.float32.g.vcf.gz"
        native_streamed_gvcf = work / "native.streamed.g.vcf.gz"
        gatk_default_gvcf = work / "gatk.default.g.vcf"
        native_default_gvcf = work / "native.default.g.vcf"
        manifest = work / "native.json"
        gvcf_manifest = work / "native.g.json"
        streamed_gvcf_manifest = work / "native.streamed.g.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_vcf),
            "-L", region, "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_vcf), "--output-manifest", str(manifest),
            "--threads", "2", "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # The default native compatibility path deliberately uses Float64,
        # but the public PairHMM precision toggle must preserve GATK's
        # Float32/GKL output contract too.  This fixture crosses both simple
        # SNP and assembled graph haplotypes, so compare its whole VCF rather
        # than only accepting that the Float32 kernel executed.
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_float32_vcf),
            "-L", region, "--native-pair-hmm-threads", "2",
            "--native-pair-hmm-use-double-precision", "false",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_float32_vcf), "--threads", "2",
            "--native-pair-hmm-use-double-precision", "false",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # CRAM is a first-class HC input, not a fallback format.  Compare the
        # complete VCF emitted from the same reference-backed decoded reads;
        # this protects the C++ Host CRAM decoder before its batch reaches the
        # Kokkos graph/PairHMM kernels.
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(cram), "-O", str(gatk_cram_vcf),
            "-L", region, "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(cram), "-R", str(reference), "-L", region,
            "-O", str(native_cram_vcf), "--threads", "2",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_gvcf),
            "-L", region, "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_float32_gvcf),
            "-L", region, "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--native-pair-hmm-use-double-precision", "false",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_gvcf), "-ERC", "GVCF", "--output-manifest", str(gvcf_manifest),
            "--threads", "2",
            "--min-depth", "1", "--min-alt-support", "1",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_float32_gvcf), "-ERC", "GVCF", "--threads", "2",
            "--native-pair-hmm-use-double-precision", "false",
            "--min-depth", "1", "--min-alt-support", "1",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Streaming is a Host-side scheduling mode: every tile must publish
        # the same aggregate Kokkos RCM result after its core crop.  Compare
        # it directly to the GATK full interval so an incorrect GQ-band merge
        # or tile-boundary END/DP/GQ/PL reduction cannot hide behind a native
        # aggregate comparison.
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_streamed_gvcf), "-ERC", "GVCF",
            "--stream-by-region", "400", "--output-manifest", str(streamed_gvcf_manifest),
            "--threads", "2",
            "--min-depth", "1", "--min-alt-support", "1",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Default GATK output includes a stable ##source line together with a
        # dynamic command record.  Exercise that real default explicitly: all
        # headers other than the provenance record must be verbatim-identical,
        # and the native writer must not smuggle private status fields into a
        # gVCF consumed by GenotypeGVCFs.
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-O", str(gatk_default_gvcf),
            "-L", region, "-ERC", "GVCF", "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_default_gvcf), "-ERC", "GVCF", "--threads", "2",
            "--min-depth", "1", "--min-alt-support", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gatk_records = variant_records(gatk_vcf)
        native_records = variant_records(native_vcf)
        gatk_gvcf_records = variant_records(gatk_gvcf)
        native_gvcf_records = variant_records(native_gvcf)
        assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
            "normal VCF is not byte-identical when command-line provenance is disabled")
        assert gatk_float32_vcf.read_bytes() == native_float32_vcf.read_bytes(), (
            "Float32 PairHMM normal VCF is not byte-identical when command-line provenance is disabled")
        assert gatk_cram_vcf.read_bytes() == native_cram_vcf.read_bytes(), (
            "CRAM-decoded normal VCF is not byte-identical when command-line provenance is disabled")
        assert vcf_payload_bytes(gatk_gvcf) == vcf_payload_bytes(native_gvcf), (
            "decompressed gVCF is not byte-identical when command-line provenance is disabled")
        assert vcf_payload_bytes(gatk_float32_gvcf) == vcf_payload_bytes(native_float32_gvcf), (
            "Float32 PairHMM decompressed gVCF is not byte-identical when command-line provenance is disabled")
        assert vcf_payload_bytes(gatk_gvcf) == vcf_payload_bytes(native_streamed_gvcf), (
            "streamed decompressed gVCF is not byte-identical when command-line provenance is disabled")
        gatk_gvcf_candidates = [record for record in gatk_gvcf_records
                                if record[4] != "<NON_REF>"]
        native_gvcf_candidates = [record for record in native_gvcf_records
                                  if record[4] != "<NON_REF>"]
        assert gatk_gvcf_candidates == native_gvcf_candidates, (
            "broad HC candidate-site gVCF rows differ: "
            f"gatk={gatk_gvcf_candidates} native={native_gvcf_candidates}")
        # The chr17 fixture contains an isolated singleton mismatch at
        # 17:69124.  GATK's default min-pruning=2 keeps it in the pileup
        # likelihood matrix but does not emit it as an assembled candidate;
        # the native path must fold the same site back into RCM without
        # perturbing neighboring PairHMM calls.
        assert all(record[1] != "69124" for record in native_gvcf_candidates)
        assert header_definitions(gatk_gvcf) == header_definitions(native_gvcf), (
            "candidate-site gVCF schema header mismatch: "
            f"gatk-only={sorted(header_definitions(gatk_gvcf) - header_definitions(native_gvcf))} "
            f"native-only={sorted(header_definitions(native_gvcf) - header_definitions(gatk_gvcf))}")
        assert stable_header(gatk_gvcf) == stable_header(native_gvcf), (
            "candidate-site gVCF stable header mismatch after excluding provenance")
        assert stable_header(gatk_default_gvcf) == stable_header(native_default_gvcf), (
            "default gVCF header mismatch after excluding dynamic GATKCommandLine")
        native_default_headers = stable_header(native_default_gvcf)
        assert "##source=HaplotypeCaller" in native_default_headers
        assert not any(line.startswith("##fastgatk_") for line in native_default_headers)
        gatk_header_definitions = header_definitions(gatk_vcf)
        native_header_definitions = header_definitions(native_vcf)
        assert gatk_header_definitions == native_header_definitions, (
            "normal VCF schema header mismatch: "
            f"gatk-only={sorted(gatk_header_definitions - native_header_definitions)} "
            f"native-only={sorted(native_header_definitions - gatk_header_definitions)}")
        assert stable_header(gatk_vcf) == stable_header(native_vcf), (
            "normal VCF stable header mismatch after excluding dynamic GATKCommandLine")
        gatk_keys = {key(record) for record in gatk_records}
        native_keys = {key(record) for record in native_records}
        assert gatk_keys == native_keys, (
            f"broad HC call-set mismatch: gatk-only={sorted(gatk_keys - native_keys)} "
            f"native-only={sorted(native_keys - gatk_keys)}")
        expected = {
            ("17", 69067, "T", "G"),
            ("17", 69368, "G", "C"),
            ("17", 69631, "C", "T"),
        }
        assert gatk_keys == expected, f"unexpected pinned GATK call set: {sorted(gatk_keys)}"

        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        assert manifest_data["compatibility"]["contig_assembly_metadata"] is True
        # The broad run's GVCF candidate rows use the standard GATK field
        # contract (including SB and uncapped Number=G PL).
        gvcf_manifest_data = json.loads(gvcf_manifest.read_text(encoding="utf-8"))
        assert gvcf_manifest_data["compatibility"]["gvcf_candidate_standard_fields"] is True
        streamed_gvcf_manifest_data = json.loads(streamed_gvcf_manifest.read_text(encoding="utf-8"))
        streamed_gvcf_telemetry = streamed_gvcf_manifest_data["telemetry"]
        assert streamed_gvcf_telemetry["stream_by_region"] is True
        assert streamed_gvcf_telemetry["streamed_regions"] == 3
        # All runs above are reference-backed active regions.  A matching VCF
        # must therefore come from the normal graph/PairHMM ownership path,
        # not from the compatibility path used only by reference-less smoke
        # inputs or candidates outside every AssemblyRegion.
        for output_kind, output_manifest in (
            ("vcf", manifest_data),
            ("gvcf", gvcf_manifest_data),
        ):
            output_telemetry = output_manifest["telemetry"]
            assert output_telemetry["pairhmm_skip_reason"] == "executed", (
                output_kind, output_telemetry["pairhmm_skip_reason"])
            assert output_telemetry["pairhmm_unassigned_candidates"] == 0, (
                output_kind, output_telemetry["pairhmm_unassigned_candidates"])
            assert output_telemetry["assembly_unassigned_candidates"] == 0, (
                output_kind, output_telemetry["assembly_unassigned_candidates"])
        assert telemetry["candidate_softclip_suppressed"] >= 1
        assert telemetry["candidate_fragment_suppressed"] >= 1
        # The broad interval may prune the weak site before final HC
        # genotyping (rather than in the explicit low-support fold-back), so
        # the counter is diagnostic and is allowed to be zero here.  The
        # candidate-set equality above is the normative assertion.
        # The evidence gate is intentionally before PairHMM.  Confirm that
        # the two rejected artifacts are absent from the native output, not
        # merely filtered by the final QUAL threshold.
        assert all(position not in {69298, 69929}
                   for _, position, _, _ in native_keys)

        gatk_by_key = {key(record): record for record in gatk_records}
        native_by_key = {key(record): record for record in native_records}
        qual_delta = {}
        pl_exact = {}
        format_exact = {}
        info_exact = {}
        annotation_exact = {}
        record_text_exact = {}
        for variant in sorted(gatk_keys):
            gatk = gatk_by_key[variant]
            native_record = native_by_key[variant]
            qual_delta[f"{variant[0]}:{variant[1]}"] = round(
                abs(float(gatk[5]) - float(native_record[5])), 6)
            gatk_sample = sample(gatk)
            native_sample = sample(native_record)
            key_name = f"{variant[0]}:{variant[1]}"
            pl_exact[key_name] = gatk_sample.get("PL") == native_sample.get("PL")
            info_exact[key_name] = info(gatk)["DP"] == info(native_record).get("DP")
            gatk_info = info(gatk)
            native_info = info(native_record)
            annotation_fields = (
                "MLEAC", "MLEAF", "MQ", "QD", "FS", "SOR", "MQRankSum",
                "ReadPosRankSum", "BaseQRankSum")
            annotation_exact[key_name] = all(
                gatk_info.get(field) == native_info.get(field)
                for field in annotation_fields
                if field in gatk_info)
            # The normal VCF body is now emitted in GATK's field order and
            # formatting.  Keep this stricter than dictionary-level checks so
            # a future writer change cannot silently reintroduce textual drift.
            record_text_exact[key_name] = "\t".join(gatk) == "\t".join(native_record)
            format_exact[key_name] = all(
                gatk_sample.get(field) == native_sample.get(field)
                for field in ("GT", "AD", "DP", "GQ"))

        # The assembly-pruning and PairHMM input gates are intended to make
        # the release-pinned core call fields exact.  Keep that boundary
        # The normal VCF records are compared as text above.  Whole-file byte
        # identity still includes headers, provenance and command-line/date
        # metadata, which this focused oracle deliberately does not claim.
        assert all(delta == 0 for delta in qual_delta.values()), qual_delta
        assert all(info_exact.values()), info_exact
        assert all(format_exact.values()), format_exact
        assert all(pl_exact.values()), pl_exact
        assert all(annotation_exact.values()), annotation_exact
        assert all(record_text_exact.values()), record_text_exact
        core_fields_bit_identical = (
            all(delta == 0 for delta in qual_delta.values()) and
            all(info_exact.values()) and
            all(format_exact.values()) and
            all(pl_exact.values()))
        core_plus_annotations_bit_identical = core_fields_bit_identical and all(annotation_exact.values())

        print(json.dumps({
            "status": "pass",
            "region": region,
            "gatk_calls": len(gatk_keys),
            "native_calls": len(native_keys),
            "call_set_exact": True,
            "header_definitions_exact": True,
            "stable_header_exact_excluding_command_line": True,
            "gvcf_candidate_records_exact": len(gatk_gvcf_candidates),
            "gvcf_decompressed_full_file_bit_identical_without_provenance": True,
            "float32_gvcf_decompressed_full_file_bit_identical_without_provenance": True,
            "gvcf_streamed_decompressed_full_file_bit_identical_without_provenance": True,
            "gvcf_streamed_regions": streamed_gvcf_telemetry["streamed_regions"],
            "gvcf_header_definitions_exact": True,
            "gvcf_stable_header_exact_excluding_command_line": True,
            "gvcf_default_header_exact_excluding_command_line": True,
            "gvcf_default_private_headers": 0,
            "candidate_softclip_suppressed": telemetry["candidate_softclip_suppressed"],
            "candidate_fragment_suppressed": telemetry["candidate_fragment_suppressed"],
            "candidate_low_support_suppressed": telemetry["candidate_low_support_suppressed"],
            "contig_assembly_metadata": True,
            "qual_abs_delta": qual_delta,
            "format_exact": format_exact,
            "info_dp_exact": info_exact,
            "pl_exact": pl_exact,
            "annotation_exact": annotation_exact,
            "record_text_exact": record_text_exact,
            "normal_vcf_bit_identical_without_provenance": True,
            "float32_normal_vcf_bit_identical_without_provenance": True,
            "cram_normal_vcf_bit_identical_without_provenance": True,
            "core_fields_bit_identical": core_fields_bit_identical,
            "core_plus_annotations_bit_identical": core_plus_annotations_bit_identical,
            "bit_identical": False,
            "bit_identical_scope": "normal VCF and decompressed gVCF are byte-identical with --add-output-vcf-command-line false; BGZF container bytes and dynamic command provenance remain outside this oracle",
        }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
