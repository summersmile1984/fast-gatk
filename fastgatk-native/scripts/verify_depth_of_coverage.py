#!/usr/bin/env python3
"""Contract and pinned GATK-oracle checks for native DepthOfCoverage."""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_DEPTH_OF_COVERAGE_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-depth-of-coverage"),
))
BAM = ROOT / "gatk-source" / "src" / "test" / "resources" / "NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source" / "src" / "test" / "resources" / "human_g1k_v37.chr17_1Mb.fasta"
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"


def digest(path: pathlib.Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def run_native(output: pathlib.Path, manifest: pathlib.Path) -> None:
    subprocess.run(
        [str(BINARY), "-R", str(REFERENCE), "-I", str(BAM), "-L", "17:69000-69005",
         "-O", str(output), "--output-manifest", str(manifest), "--batch-records", "2"],
        check=True, text=True, capture_output=True,
    )


def write_reference(path: pathlib.Path, contig: str, sequence: str) -> None:
    path.write_text(f">{contig}\n" + sequence + "\n", encoding="ascii")
    path.with_suffix(path.suffix + ".fai").write_text(
        f"{contig}\t{len(sequence)}\t{len(contig) + 2}\t{len(sequence)}\t{len(sequence) + 1}\n",
        encoding="ascii",
    )
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:{contig}\tLN:{len(sequence)}\n", encoding="ascii",
    )


def write_n_reference(path: pathlib.Path) -> None:
    write_reference(path, "chrN", "ACNNTA")


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not BAM.exists() or not REFERENCE.exists():
        raise SystemExit("missing pinned DepthOfCoverage fixture")
    with tempfile.TemporaryDirectory(prefix="fastgatk-depth-of-coverage-") as temporary:
        work = pathlib.Path(temporary)
        native = work / "native.out"
        native_manifest = work / "native.manifest.json"
        run_native(native, native_manifest)
        assert native.read_text(encoding="utf-8") == (
            "Locus,Total_Depth,Average_Depth_sample,Depth_for_NA12878\n"
            "17:69000,5,5.00,5\n17:69001,5,5.00,5\n17:69002,5,5.00,5\n"
            "17:69003,5,5.00,5\n17:69004,5,5.00,5\n17:69005,5,5.00,5\n"
        )
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert manifest["tool"] == "DepthOfCoverage"
        assert manifest["execution_space"] in {"OpenMP", "Serial"}
        assert manifest["determinism"] == "strict"
        assert manifest["compatibility"]["count_type"] == "COUNT_READS"
        assert manifest["telemetry"]["reads_seen"] == 5
        assert manifest["telemetry"]["observations"] == 30
        telemetry = manifest["telemetry"]
        assert telemetry["count_kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["count_kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["count_kernel_execution_policy"] == "RangePolicy"
        assert telemetry["count_kernel_batches"] >= 1
        assert telemetry["count_kernel_observations"] == telemetry["observations"]
        assert telemetry["count_kernel_prepare_seconds"] >= 0.0
        assert telemetry["count_kernel_execute_seconds"] >= 0.0
        assert telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] >= 1
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0
        assert manifest["telemetry"]["wall_seconds"] >= 0.0
        assert all(entry["complete"] for entry in manifest["outputs"])

        expected_sidecars = [
            native,
            pathlib.Path(str(native) + ".sample_interval_summary"),
            pathlib.Path(str(native) + ".sample_interval_statistics"),
            pathlib.Path(str(native) + ".sample_summary"),
            pathlib.Path(str(native) + ".sample_statistics"),
            pathlib.Path(str(native) + ".sample_cumulative_coverage_counts"),
            pathlib.Path(str(native) + ".sample_cumulative_coverage_proportions"),
        ]
        assert all(path.exists() and path.stat().st_size > 0 for path in expected_sidecars)

        gatk_ready = oracle_guard.oracle_ready('verify_depth_of_coverage.py', JAVA, GATK_JAR)
        oracle = False
        if gatk_ready:
            gatk_prefix = work / "gatk.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(REFERENCE),
                 "-I", str(BAM), "-L", "17:69000-69005", "-O", str(gatk_prefix),
                 "--omit-locus-table", "false", "--omit-interval-statistics", "false",
                 "--omit-per-sample-statistics", "false", "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            for native_path in expected_sidecars:
                gatk_path = pathlib.Path(str(gatk_prefix) + native_path.name[len(native.name):])
                assert gatk_path.exists(), gatk_path
                assert digest(native_path) == digest(gatk_path), (native_path.name, gatk_path.name)
            oracle = True

        # Reference-N handling is an important direct-replacement boundary:
        # GATK excludes N loci unless --include-ref-n-sites is explicit.  Use
        # a tiny SAM/FASTA fixture so this assertion is independent of the
        # chr17 resource's reference composition.
        n_reference = work / "n-reference.fasta"
        write_n_reference(n_reference)
        n_sam = work / "n-reference.sam"
        n_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chrN\tLN:6\n"
            "@RG\tID:rgN\tSM:N_SAMPLE\n"
            "nread\t0\tchrN\t1\t60\t6M\t*\t0\t0\tAAAAAA\tIIIIII\tRG:Z:rgN\n",
            encoding="ascii",
        )
        n_default = work / "n-default.out"
        n_default_manifest = work / "n-default.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam), "-L", "chrN:1-6",
             "-O", str(n_default), "--output-manifest", str(n_default_manifest)],
            check=True, text=True, capture_output=True,
        )
        n_default_lines = n_default.read_text(encoding="utf-8").splitlines()
        assert n_default_lines[1:] == [
            "chrN:1,1,1.00,1", "chrN:2,1,1.00,1",
            "chrN:5,1,1.00,1", "chrN:6,1,1.00,1",
        ]
        n_default_metadata = json.loads(n_default_manifest.read_text(encoding="utf-8"))
        assert n_default_metadata["telemetry"]["loci"] == 6
        assert n_default_metadata["telemetry"]["counted_loci"] == 4
        assert n_default_metadata["telemetry"]["include_ref_n_sites"] is False

        n_include = work / "n-include.out"
        n_include_manifest = work / "n-include.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam), "-L", "chrN:1-6",
             "--include-ref-n-sites", "-O", str(n_include),
             "--output-manifest", str(n_include_manifest)],
            check=True, text=True, capture_output=True,
        )
        n_include_lines = n_include.read_text(encoding="utf-8").splitlines()
        assert len(n_include_lines) == 7 and n_include_lines[1].startswith("chrN:1,")
        n_include_metadata = json.loads(n_include_manifest.read_text(encoding="utf-8"))
        assert n_include_metadata["telemetry"]["counted_loci"] == 6
        assert n_include_metadata["telemetry"]["include_ref_n_sites"] is True

        # GATK applies include intervals, padding, then exclude intervals
        # before constructing any depth summaries.  Keep the same boundary in
        # the native path so excluded coordinates do not remain as zero-depth
        # rows or leak into quantile/cumulative denominators.
        excluded = work / "excluded.out"
        excluded_manifest = work / "excluded.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam),
             "-L", "chrN:1-6", "-XL", "chrN:2-2", "--include-ref-n-sites",
             "-O", str(excluded), "--output-manifest", str(excluded_manifest)],
            check=True, text=True, capture_output=True,
        )
        assert excluded.read_text(encoding="utf-8").splitlines()[1:] == [
            "chrN:1,1,1.00,1", "chrN:3,1,1.00,1", "chrN:4,1,1.00,1",
            "chrN:5,1,1.00,1", "chrN:6,1,1.00,1",
        ]
        excluded_metadata = json.loads(excluded_manifest.read_text(encoding="utf-8"))
        assert excluded_metadata["telemetry"]["excluded_intervals"] == 1
        assert excluded_metadata["telemetry"]["counted_loci"] == 5

        excluded_padded = work / "excluded-padded.out"
        excluded_padded_manifest = work / "excluded-padded.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam),
             "-L", "chrN:1-6", "-XL", "chrN:3-3",
             "--interval-exclusion-padding", "1", "--include-ref-n-sites",
             "-O", str(excluded_padded), "--output-manifest", str(excluded_padded_manifest)],
            check=True, text=True, capture_output=True,
        )
        assert excluded_padded.read_text(encoding="utf-8").splitlines()[1:] == [
            "chrN:1,1,1.00,1", "chrN:5,1,1.00,1", "chrN:6,1,1.00,1",
        ]
        excluded_padded_metadata = json.loads(excluded_padded_manifest.read_text(encoding="utf-8"))
        assert excluded_padded_metadata["compatibility"]["interval_exclusion_padding"] == 1

        padded = work / "padded.out"
        padded_manifest = work / "padded.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam),
             "-L", "chrN:2-2", "--interval-padding", "1",
             "--include-ref-n-sites", "-O", str(padded),
             "--output-manifest", str(padded_manifest)],
            check=True, text=True, capture_output=True,
        )
        assert padded.read_text(encoding="utf-8").splitlines()[1:] == [
            "chrN:1,1,1.00,1", "chrN:2,1,1.00,1", "chrN:3,1,1.00,1",
        ]
        padded_metadata = json.loads(padded_manifest.read_text(encoding="utf-8"))
        assert padded_metadata["compatibility"]["interval_padding"] == 1

        intersected = work / "intersected.out"
        intersected_manifest = work / "intersected.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(n_reference), "-I", str(n_sam),
             "-L", "chrN:1-4", "-L", "chrN:3-6", "--interval-set-rule",
             "INTERSECTION", "--include-ref-n-sites", "-O", str(intersected),
             "--output-manifest", str(intersected_manifest)],
            check=True, text=True, capture_output=True,
        )
        assert intersected.read_text(encoding="utf-8").splitlines()[1:] == [
            "chrN:3,1,1.00,1", "chrN:4,1,1.00,1",
        ]
        intersected_metadata = json.loads(intersected_manifest.read_text(encoding="utf-8"))
        assert intersected_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"

        def run_output_omission(prefix: pathlib.Path, *flags: str) -> set[str]:
            subprocess.run(
                [str(BINARY), "-R", str(n_reference), "-I", str(n_sam), "-L", "chrN:1-6",
                 "-O", str(prefix), *flags],
                check=True, text=True, capture_output=True,
            )
            return {path.name for path in prefix.parent.glob(prefix.name + "*")}

        interval_omitted = run_output_omission(
            work / "omit-intervals.out", "--omit-interval-statistics"
        )
        assert "omit-intervals.out.sample_interval_summary" not in interval_omitted
        assert "omit-intervals.out.sample_interval_statistics" not in interval_omitted
        assert "omit-intervals.out.sample_summary" in interval_omitted
        assert "omit-intervals.out.sample_statistics" in interval_omitted
        sample_omitted = run_output_omission(
            work / "omit-samples.out", "--omit-per-sample-statistics"
        )
        assert "omit-samples.out.sample_interval_summary" in sample_omitted
        assert "omit-samples.out.sample_interval_statistics" in sample_omitted
        assert "omit-samples.out.sample_summary" not in sample_omitted
        assert "omit-samples.out.sample_statistics" not in sample_omitted
        aliases_omitted = run_output_omission(
            work / "omit-aliases.out", "--omit-intervals", "--omit-sample-summary"
        )
        assert not any(name.endswith((".sample_interval_summary", ".sample_interval_statistics",
                                      ".sample_summary", ".sample_statistics"))
                       for name in aliases_omitted)

        n_java_oracle = False
        if gatk_ready:
            n_bam = work / "n-reference.bam"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "SortSam", "-I", str(n_sam),
                 "-O", str(n_bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                check=True, text=True, capture_output=True,
            )
            n_java_prefix = work / "n-java.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(n_reference),
                 "-I", str(n_bam), "-L", "chrN:1-6", "-O", str(n_java_prefix),
                 "--omit-locus-table", "false", "--omit-interval-statistics", "false",
                 "--omit-per-sample-statistics", "false", "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            assert n_default.read_bytes() == n_java_prefix.read_bytes()
            n_java_include_prefix = work / "n-java-include.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(n_reference),
                 "-I", str(n_bam), "-L", "chrN:1-6", "--include-ref-n-sites",
                 "-O", str(n_java_include_prefix), "--omit-locus-table", "false",
                 "--omit-interval-statistics", "false", "--omit-per-sample-statistics", "false",
                 "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            assert n_include.read_bytes() == n_java_include_prefix.read_bytes()
            n_java_oracle = True

            n_java_excluded_prefix = work / "n-java-excluded.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(n_reference),
                 "-I", str(n_bam), "-L", "chrN:1-6", "-XL", "chrN:2-2",
                 "--include-ref-n-sites", "-O", str(n_java_excluded_prefix),
                 "--omit-locus-table", "false", "--omit-interval-statistics", "false",
                 "--omit-per-sample-statistics", "false", "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            assert excluded.read_bytes() == n_java_excluded_prefix.read_bytes()

            n_java_excluded_padded_prefix = work / "n-java-excluded-padded.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(n_reference),
                 "-I", str(n_bam), "-L", "chrN:1-6", "-XL", "chrN:3-3",
                 "--interval-exclusion-padding", "1", "--include-ref-n-sites",
                 "-O", str(n_java_excluded_padded_prefix), "--omit-locus-table", "false",
                 "--omit-interval-statistics", "false", "--omit-per-sample-statistics", "false",
                 "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            assert excluded_padded.read_bytes() == n_java_excluded_padded_prefix.read_bytes()

            n_java_padded_prefix = work / "n-java-padded.out"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(n_reference),
                 "-I", str(n_bam), "-L", "chrN:2-2", "--interval-padding", "1",
                 "--include-ref-n-sites", "-O", str(n_java_padded_prefix),
                 "--omit-locus-table", "false", "--omit-interval-statistics", "false",
                 "--omit-per-sample-statistics", "false", "--omit-depth-output-at-each-base", "false"],
                check=True, text=True, capture_output=True,
            )
            assert padded.read_bytes() == n_java_padded_prefix.read_bytes()

        # GATK 4.6.2.0 exposes COUNT_FRAGMENTS in the command-line enum but
        # deliberately rejects it in CoverageUtils.  Keep the native mode as
        # an explicit extension (stable read-name+locus de-duplication), and
        # pin that boundary so a future GATK change cannot silently turn this
        # into an unverified compatibility claim.
        fragment_reference = work / "fragment-reference.fasta"
        write_reference(fragment_reference, "chrF", "ACGTAC")
        fragment_sam = work / "fragments.sam"
        fragment_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chrF\tLN:6\n"
            "@RG\tID:rgF\tSM:F_SAMPLE\n"
            "pair\t99\tchrF\t1\t60\t6M\t=\t1\t6\tAAAAAA\tIIIIII\tRG:Z:rgF\n"
            "pair\t147\tchrF\t1\t60\t6M\t=\t1\t-6\tCCCCCC\tIIIIII\tRG:Z:rgF\n",
            encoding="ascii",
        )
        fragment_reads = work / "fragment-reads.out"
        fragment_reads_manifest = work / "fragment-reads.manifest.json"
        subprocess.run(
            [str(BINARY), "-R", str(fragment_reference), "-I", str(fragment_sam),
             "-L", "chrF:1-6", "--count-type", "COUNT_READS", "-O", str(fragment_reads),
             "--output-manifest", str(fragment_reads_manifest)],
            check=True, text=True, capture_output=True,
        )
        assert all(line.endswith(",2,2.00,2") for line in fragment_reads.read_text(encoding="ascii").splitlines()[1:])
        fragment_output = work / "fragment-counts.out"
        fragment_manifest = work / "fragment-counts.manifest.json"
        fragment_run = subprocess.run(
            [str(BINARY), "-R", str(fragment_reference), "-I", str(fragment_sam),
             "-L", "chrF:1-6", "--count-type", "COUNT_FRAGMENTS", "-O", str(fragment_output),
             "--output-manifest", str(fragment_manifest), "--batch-records", "1"],
            check=True, text=True, capture_output=True,
        )
        assert all(line.endswith(",1,1.00,1") for line in fragment_output.read_text(encoding="ascii").splitlines()[1:])
        fragment_metadata = json.loads(fragment_manifest.read_text(encoding="utf-8"))
        assert fragment_metadata["compatibility"]["count_type"] == "COUNT_FRAGMENTS"
        assert fragment_metadata["telemetry"]["count_type"] == "COUNT_FRAGMENTS"
        assert fragment_metadata["telemetry"]["fragment_keys"] == 6
        quality_fragment_sam = work / "fragments-quality.sam"
        quality_fragment_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chrF\tLN:6\n"
            "@RG\tID:rgF\tSM:F_SAMPLE\n"
            "pair\t99\tchrF\t1\t60\t6M\t=\t1\t6\tAAAAAA\t!!!!!!\tRG:Z:rgF\n"
            "pair\t147\tchrF\t1\t60\t6M\t=\t1\t-6\tCCCCCC\tIIIIII\tRG:Z:rgF\n",
            encoding="ascii",
        )
        quality_fragment_output = work / "fragments-quality.out"
        subprocess.run(
            [str(BINARY), "-R", str(fragment_reference), "-I", str(quality_fragment_sam),
             "-L", "chrF:1-6", "--count-type", "COUNT_FRAGMENTS", "--min-base-quality", "20",
             "--batch-records", "1", "-O", str(quality_fragment_output)],
            check=True, text=True, capture_output=True,
        )
        assert all(line.endswith(",1,1.00,1") for line in quality_fragment_output.read_text(encoding="ascii").splitlines()[1:])
        fragment_gatk_supported = False
        if gatk_ready:
            fragment_bam = work / "fragments.bam"
            subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "SortSam", "-I", str(fragment_sam),
                 "-O", str(fragment_bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                check=True, text=True, capture_output=True,
            )
            java_fragment_prefix = work / "gatk-fragments"
            java_fragment = subprocess.run(
                [str(JAVA), "-jar", str(GATK_JAR), "DepthOfCoverage", "-R", str(fragment_reference),
                 "-I", str(fragment_bam), "-L", "chrF:1-6", "--count-type", "COUNT_FRAGMENTS",
                 "-O", str(java_fragment_prefix), "--omit-locus-table", "false",
                 "--omit-interval-statistics", "false", "--omit-per-sample-statistics", "false",
                 "--omit-depth-output-at-each-base", "false"],
                check=False, text=True, capture_output=True,
            )
            if java_fragment.returncode == 0:
                fragment_gatk_supported = True
                assert fragment_output.read_bytes() == java_fragment_prefix.read_bytes()
            else:
                assert "Fragment based counting is currently unsupported" in java_fragment.stderr
        print(json.dumps({
            "status": "pass",
            "tool": "DepthOfCoverage",
            "records": 6,
            "sidecars": len(expected_sidecars) - 1,
            "gatk_oracle": oracle,
            "reference_n_oracle": n_java_oracle,
            "count_fragments_extension": True,
            "gatk_count_fragments_supported": fragment_gatk_supported,
            "bit_identical_to_gatk": oracle,
            "manifest": manifest,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
