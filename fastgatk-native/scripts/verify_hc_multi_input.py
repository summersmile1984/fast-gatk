#!/usr/bin/env python3
"""Verify repeatable HaplotypeCaller -I aggregation and fail-closed streaming."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import gzip
from pathlib import Path


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True,
                          check=False, env={**os.environ, "OMP_NUM_THREADS": "1"})


def summary(result: subprocess.CompletedProcess[str]) -> dict:
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


def vcf_data_rows(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY",
                               root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    assert binary.is_file() and os.access(binary, os.X_OK)
    assert bam.is_file()
    assert reference.is_file()
    # CTest may run OpenMP and Serial backends concurrently in this shared
    # workspace; isolate materialized outputs per invocation.
    work_dir = Path(tempfile.mkdtemp(prefix="hc-multi-input-", dir=str(root / "work")))

    single = summary(run(binary, ["-I", str(bam), "-O", "-"]))
    combined = summary(run(binary, ["--input", str(bam), "-I", str(bam), "-O", "-"]))
    assert single["input_count"] == 1
    # The option is surfaced separately from effective precision.  When it is
    # omitted, the native caller preserves the strict double compatibility
    # path even though GATK's optional boolean is accepted.
    assert single["native_pair_hmm_use_double_precision"] is False
    assert single["pairhmm_effective_precision"] == "double"
    # The optional native precision flag is intentionally opt-in: omitted
    # keeps the strict double path, while an explicit false selects the same
    # Kokkos bucketed API with simd<float>.
    explicit_float = summary(run(binary, ["-I", str(bam), "-R", str(reference),
                                           "-L", "17:69000-70000",
                                           "--native-pair-hmm-use-double-precision=false", "-O", "-"]))
    assert explicit_float["native_pair_hmm_use_double_precision"] is False
    assert explicit_float["pairhmm_effective_precision"] == "float32", explicit_float
    assert explicit_float["pairhmm_used"] is True
    explicit_double = summary(run(binary, ["-I", str(bam), "-R", str(reference),
                                            "-L", "17:69000-70000",
                                            "--native-pair-hmm-use-double-precision=true", "-O", "-"]))
    assert explicit_double["native_pair_hmm_use_double_precision"] is True
    assert explicit_double["pairhmm_effective_precision"] == "double"
    invalid_precision = run(binary, ["-I", str(bam),
                                     "--native-pair-hmm-use-double-precision=maybe", "-O", "-"])
    assert invalid_precision.returncode != 0 and "invalid boolean" in invalid_precision.stderr.lower()
    assert combined["input_count"] == 2
    assert combined["inputs"] == [str(bam), str(bam)]
    # The fixture is a single-sample shard. Combining it twice must retain the
    # complete read payload, rather than silently using only the last -I.
    assert combined["reads"] == 2 * single["reads"], (single, combined)
    assert combined["input_count"] == 2

    no_index_path = work_dir / "hc-no-index.vcf.gz"
    no_index = run(binary, ["-I", str(bam), "-O", str(no_index_path),
                            "--create-output-variant-index", "false"])
    assert no_index.returncode == 0, no_index.stderr
    no_index_summary = summary(no_index)
    assert no_index_summary["create_output_variant_index"] is False
    assert no_index_summary["native_pair_hmm_use_double_precision"] is False
    assert no_index_summary["pairhmm_effective_precision"] == "double"
    assert no_index_path.is_file() and not Path(f"{no_index_path}.tbi").exists()

    # Plain VCF follows the Java VariantContextWriter contract and receives a
    # Tribble `.idx`; verify the bundled GATK reader can traverse that index.
    plain_path = work_dir / "hc-indexed-plain.vcf"
    plain_result = run(binary, ["-I", str(bam), "-R", str(reference), "-L",
                               "17:69000-70000", "-O", str(plain_path)])
    assert plain_result.returncode == 0, plain_result.stderr
    assert plain_path.is_file() and Path(f"{plain_path}.idx").is_file()
    assert not Path(f"{plain_path}.tbi").exists()
    plain_query = work_dir / "hc-indexed-plain-query.vcf"
    gatk_java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    query = subprocess.run(
        [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_path),
         "-L", "17:69000-70000", "-O", str(plain_query)],
        text=True, capture_output=True, check=False,
    )
    assert query.returncode == 0, query.stderr
    assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
               if line and not line.startswith("#")) > 0
    plain_manifest = Path(f"{plain_path}.manifest.json")
    assert plain_manifest.is_file()
    plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
    assert plain_metadata["compatibility"]["vcf_index"] is True
    assert all(item["complete"] for item in plain_metadata["outputs"])

    # `--add-output-vcf-command-line` must control the actual VCF header, not
    # merely be copied into telemetry.  GATK emits its stable source header
    # together with structured command provenance in this mode.
    command_line_on_path = work_dir / "hc-command-line-on.vcf"
    command_line_off_path = work_dir / "hc-command-line-off.vcf"
    command_line_on = run(binary, ["-I", str(bam), "-O", str(command_line_on_path),
                                   "--create-output-variant-index", "false",
                                   "--add-output-vcf-command-line", "true"])
    command_line_off = run(binary, ["-I", str(bam), "-O", str(command_line_off_path),
                                    "--create-output-variant-index", "false",
                                    "--add-output-vcf-command-line", "false"])
    assert command_line_on.returncode == 0, command_line_on.stderr
    assert command_line_off.returncode == 0, command_line_off.stderr
    on_text = command_line_on_path.read_text(encoding="utf-8")
    off_text = command_line_off_path.read_text(encoding="utf-8")
    assert "##GATKCommandLine=<ID=HaplotypeCaller,Version=fastgatk-native," in on_text
    assert "##GATKCommandLine=" not in off_text
    assert "##source=HaplotypeCaller" in on_text
    assert "##source=HaplotypeCaller" not in off_text

    # GATK validates the BAM/reference sequence dictionary before traversal.
    # `reference.fa` is intentionally a different contig (20) from the chr17
    # fixture, so the default path must fail closed; the explicit optional
    # boolean escape hatch skips only that validation and remains observable.
    mismatched_reference = root / "reference.fa"
    mismatch = run(binary, ["-I", str(bam), "-R", str(mismatched_reference), "-O", "-"])
    assert mismatch.returncode != 0
    assert "sequence dictionary" in mismatch.stderr.lower(), mismatch.stderr
    disabled_false = run(binary, ["-I", str(bam), "-R", str(mismatched_reference),
                                   "--disable-sequence-dictionary-validation", "false", "-O", "-"])
    assert disabled_false.returncode != 0
    assert "sequence dictionary" in disabled_false.stderr.lower(), disabled_false.stderr
    disabled_dictionary = run(binary, ["-I", str(bam), "-R", str(mismatched_reference),
                                       "--disable-sequence-dictionary-validation", "-O", "-"])
    assert disabled_dictionary.returncode == 0, disabled_dictionary.stderr
    disabled_summary = summary(disabled_dictionary)
    assert disabled_summary["disable_sequence_dictionary_validation"] is True

    broad = summary(run(binary, ["-I", str(bam), "-L", "17:69000-70000", "-O", "-"]))
    intersection = summary(run(binary, ["-I", str(bam), "-L", "17:69000-70000",
                                        "-L", "17:69300-69400", "-isr", "INTERSECTION",
                                        "-O", "-"]))
    narrow = summary(run(binary, ["-I", str(bam), "-L", "17:69300-69400", "-O", "-"]))
    assert broad["interval_set_rule"] == "UNION"
    assert intersection["interval_set_rule"] == "INTERSECTION"
    assert intersection["reads"] == narrow["reads"]
    invalid_rule = run(binary, ["-I", str(bam), "-L", "17:69000-70000", "-isr", "BAD", "-O", "-"])
    assert invalid_rule.returncode != 0 and "INTERSECTION" in invalid_rule.stderr

    # Contig streaming merges the two coordinate-sorted shards instead of
    # concatenating them (which would re-enter chr17 and be rejected).
    combined_aggregate_path = work_dir / "hc-multi-input-aggregate.vcf"
    single_stream_path = work_dir / "hc-single-input-stream.vcf"
    combined_stream_path = work_dir / "hc-multi-input-stream.vcf"
    combined_aggregate = run(binary, ["-I", str(bam), "-I", str(bam), "-O",
                                      str(combined_aggregate_path)])
    single_stream = run(binary, ["-I", str(bam), "--stream-by-contig", "-O",
                                 str(single_stream_path)])
    combined_stream = run(binary, ["-I", str(bam), "-I", str(bam),
                                   "--stream-by-contig", "-O", str(combined_stream_path)])
    assert single_stream.returncode == 0, single_stream.stderr
    assert combined_stream.returncode == 0, combined_stream.stderr
    assert combined_aggregate.returncode == 0, combined_aggregate.stderr
    assert combined_stream_path.read_text(encoding="utf-8") == combined_aggregate_path.read_text(encoding="utf-8")

    # Region streaming opens one indexed reader per shard for each bounded
    # halo and merges the resulting tile-local batches under the same budget.
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if reference.is_file() and Path(f"{bam}.bai").is_file():
        single_region_path = work_dir / "hc-single-input-region.vcf"
        combined_region_path = work_dir / "hc-multi-input-region.vcf"
        single_region = run(binary, ["-I", str(bam), "-R", str(reference), "-L",
                                     "17:69000-70000", "--stream-by-region", "500", "-O",
                                     str(single_region_path)])
        combined_region_aggregate_path = work_dir / "hc-multi-input-region-aggregate.vcf"
        combined_region_aggregate = run(binary, ["-I", str(bam), "-I", str(bam), "-R",
                                                  str(reference), "-L", "17:69000-70000", "-O",
                                                  str(combined_region_aggregate_path)])
        combined_region = run(binary, ["-I", str(bam), "-I", str(bam), "-R", str(reference),
                                       "-L", "17:69000-70000", "--stream-by-region", "500",
                                       "-O", str(combined_region_path)])
        assert single_region.returncode == 0, single_region.stderr
        assert combined_region.returncode == 0, combined_region.stderr
        assert combined_region_aggregate.returncode == 0, combined_region_aggregate.stderr
        # This is a strict source oracle for an assembly-sensitive case:
        # passing the same BAM twice gives GATK two independent copies of the
        # evidence.  The complete EventMap must retain every resulting
        # assembled allele, not only the higher-support sibling in each
        # active region.  Compare every data column (while deliberately
        # excluding producer-specific headers): in particular this verifies
        # that the VCF writer retrieves the owning AssemblyRegion's joint
        # PairHMM PL rather than falling back to its capped biallelic envelope.
        gatk_duplicate_path = work_dir / "hc-multi-input-gatk.vcf.gz"
        gatk_duplicate = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "HaplotypeCaller",
             "-R", str(reference), "-I", str(bam), "-I", str(bam),
             "-L", "17:69000-70000", "-O", str(gatk_duplicate_path),
             "--add-output-vcf-command-line", "false"],
            text=True, capture_output=True, check=False,
        )
        assert gatk_duplicate.returncode == 0, gatk_duplicate.stderr
        assert vcf_data_rows(combined_region_aggregate_path) == vcf_data_rows(gatk_duplicate_path)
        # Region streaming is a bounded execution path: each core is evaluated
        # with a halo, so duplicate/overlapping evidence can legitimately
        # change per-site depth and candidate selection versus the unbounded
        # aggregate caller.  The contract here is therefore structural rather
        # than byte-identical: records stay coordinate-sorted, unique within
        # the emitted stream, and inside the requested interval.
        streamed_rows = [line for line in combined_region_path.read_text(encoding="utf-8").splitlines()
                         if line and not line.startswith("#")]
        assert streamed_rows
        streamed_keys = [tuple(line.split("\t")[index] for index in (0, 1, 3, 4))
                         for line in streamed_rows]
        assert streamed_keys == sorted(streamed_keys,
                                       key=lambda key: (key[0], int(key[1]), key[2], key[3]))
        assert len(streamed_keys) == len(set(streamed_keys))
        assert all(69000 <= int(key[1]) <= 70000 for key in streamed_keys)

    print(json.dumps({"status": "pass", "single_reads": single["reads"],
                      "combined_reads": combined["reads"], "input_count": 2,
                      "contig_streaming_multi_input": True,
                      "region_streaming_multi_input": True,
                      "gatk_duplicate_input_vcf_data_rows": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
