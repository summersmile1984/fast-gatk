#!/usr/bin/env python3
"""Verify deterministic duplicate marking/removal, metrics and index contracts."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
INST:1:1101:100:100\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII
INST:1:1101:101:101\t0\tchr1\t11\t20\t4M\t*\t0\t0\tACGT\t!!!!
unique\t0\tchr1\t21\t60\t4M\t*\t0\t0\tTGCA\tIIII
"""

PAIRED_SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
@RG\tID:rg1\tSM:S1\tLB:lib1
PAIR_A:1:1101:100:100\t99\tchr1\t11\t60\t4M\t=\t31\t24\tACGT\tIIII\tMC:Z:4M\tRG:Z:rg1
PAIR_B:1:1101:101:101\t99\tchr1\t11\t20\t4M\t=\t31\t24\tACGT\t!!!!\tMC:Z:4M\tRG:Z:rg1
PAIR_C:1:1101:102:102\t0\tchr1\t21\t30\t4M\t*\t0\t0\tAAAA\tIIII\tRG:Z:rg1
PAIR_A:1:1101:100:100\t147\tchr1\t31\t60\t4M\t=\t11\t-24\tTGCA\t!!!!\tMC:Z:4M\tRG:Z:rg1
PAIR_B:1:1101:101:101\t147\tchr1\t31\t20\t4M\t=\t11\t-24\tTGCA\tIIII\tMC:Z:4M\tRG:Z:rg1
"""

MULTI_LIBRARY_SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
@RG\tID:rgA\tSM:S1\tLB:libA
@RG\tID:rgB\tSM:S1\tLB:libB
a-high\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rgA
a-low\t0\tchr1\t11\t20\t4M\t*\t0\t0\tACGT\t!!!!\tRG:Z:rgA
b-high\t0\tchr1\t21\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:rgB
b-low\t0\tchr1\t21\t20\t4M\t*\t0\t0\tTGCA\t!!!!\tRG:Z:rgB
"""

GIANT_GROUP_SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
g0\t0\tchr1\t31\t60\t4M\t*\t0\t0\tACGT\tIIII
g1\t0\tchr1\t31\t50\t4M\t*\t0\t0\tACGT\tIIII
g2\t0\tchr1\t31\t40\t4M\t*\t0\t0\tACGT\tIIII
g3\t0\tchr1\t31\t30\t4M\t*\t0\t0\tACGT\tIIII
g4\t0\tchr1\t31\t20\t4M\t*\t0\t0\tACGT\tIIII
g5\t0\tchr1\t31\t10\t4M\t*\t0\t0\tACGT\tIIII
g6\t0\tchr1\t31\t5\t4M\t*\t0\t0\tACGT\tIIII
"""


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def semantic_records(path: Path) -> list[tuple[str, int, int, str]]:
    """Return the Picard-visible duplicate semantics, ignoring PG provenance."""
    result = []
    for row in records(path):
        dt = next((field for field in row[11:] if field.startswith("DT:Z:")), "")
        result.append((row[0], int(row[1]), int(row[3]), dt))
    return result


def metric_table(path: Path) -> dict[str, list[str]]:
    """Parse the standard Picard DuplicationMetrics table."""
    rows: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("##"):
            continue
        fields = line.split("\t")
        if fields[0] == "LIBRARY" or len(fields) < 10:
            continue
        rows[fields[0]] = fields
    return rows


def run(binary: Path, source: Path, output: Path, manifest: Path, metrics: Path, *extra: str) -> dict:
    result = subprocess.run(
        [str(binary), "-I", str(source), "-O", str(output),
         "--metrics-file", str(metrics), "--output-manifest", str(manifest), *extra],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-mark-duplicates"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-mark-duplicates-") as temp:
        work = Path(temp)
        source = work / "input.sam"
        source.write_text(SAM, encoding="utf-8")

        marked = work / "marked.sam"
        marked_manifest = work / "marked.manifest.json"
        marked_metrics = work / "marked.metrics.txt"
        scratch = work / "scratch"
        summary = run(binary, source, marked, marked_manifest, marked_metrics,
                      "--tmp-dir", str(scratch), "--max-records-in-memory", "2",
                      "--tagging-policy", "All")
        assert summary["status"] == "prototype"
        assert summary["input_records"] == 3 and summary["output_records"] == 3
        assert summary["duplicate_records"] == 1
        # Picard only classifies paired fragments as optical duplicates.  The
        # two unpaired reads use sequencer-like names but remain a normal
        # library duplicate and receive DT:Z:LB when tagging is enabled.
        assert summary["optical_duplicate_records"] == 0
        assert marked.exists()
        marked_lines = marked.read_text(encoding="utf-8").splitlines()
        assert any(line.startswith("@PG\t") and "ID:MarkDuplicates" in line
                   for line in marked_lines)
        rows = records(marked)
        assert len(rows) == 3
        assert int(rows[1][1]) & 0x400
        assert "DT:Z:LB" in rows[1][11:]
        metrics_text = marked_metrics.read_text(encoding="utf-8")
        unpaired_metrics = metric_table(marked_metrics)["Unknown Library"]
        assert unpaired_metrics[5] == "1"  # UNPAIRED_READ_DUPLICATES
        assert unpaired_metrics[6] == "0"  # READ_PAIR_DUPLICATES
        assert unpaired_metrics[7] == "0"  # READ_PAIR_OPTICAL_DUPLICATES
        metadata = json.loads(marked_manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["mark_or_remove"] is True
        assert metadata["compatibility"]["bounded_group"] is False
        assert metadata["compatibility"]["spill_runs"] is True
        assert metadata["telemetry"]["spill_run_count"] >= 1
        assert metadata["compatibility"]["optical_duplicate_detection"] is True
        assert metadata["compatibility"]["max_records_in_memory"] == 2
        assert scratch.is_dir()
        assert not list(scratch.glob("fastgatk-markduplicates-run-*.bin"))
        assert all(item["complete"] for item in metadata["outputs"])

        pg_marked = work / "marked-pg.sam"
        pg_marked_manifest = work / "marked-pg.manifest.json"
        pg_marked_metrics = work / "marked-pg.metrics.txt"
        pg_result = run(binary, source, pg_marked, pg_marked_manifest, pg_marked_metrics,
                         "--add-pg-tag=true", "--program-record-id", "md-test",
                         "--program-group-name", "MarkDuplicates", "--program-group-version", "2.1",
                         "--program-group-command-line", "fastgatk MarkDuplicates test")
        assert pg_result["status"] == "prototype"
        pg_lines = pg_marked.read_text(encoding="utf-8").splitlines()
        assert any(line.startswith("@PG\t") and "ID:md-test" in line and
                   "PN:MarkDuplicates" in line and "VN:2.1" in line and
                   "CL:fastgatk MarkDuplicates test" in line for line in pg_lines)
        pg_metadata = json.loads(pg_marked_manifest.read_text(encoding="utf-8"))
        assert pg_metadata["compatibility"]["pg_line"] is True

        invalid_pg = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(work / "invalid-pg.sam"),
             "--add-pg-tag=maybe"], text=True, capture_output=True, check=False)
        assert invalid_pg.returncode != 0 and "invalid boolean" in invalid_pg.stderr
        invalid_index = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(work / "invalid-index.bam"),
             "--create-output-bam-index=maybe"], text=True, capture_output=True, check=False)
        assert invalid_index.returncode != 0 and "invalid boolean" in invalid_index.stderr

        # Without scratch, a strict bound remains fail-closed rather than
        # silently growing an in-memory duplicate group.
        bounded = subprocess.run([
            str(binary), "-I", str(source), "-O", str(work / "bounded.sam"),
            "--metrics-file", str(work / "bounded.metrics.txt"),
            "--max-records-in-memory", "1",
        ], text=True, capture_output=True, check=False)
        assert bounded.returncode != 0 and "RESOURCE_EXHAUSTED" in bounded.stderr

        paired_source = work / "paired.sam"
        paired_source.write_text(PAIRED_SAM, encoding="utf-8")
        paired_output = work / "paired-marked.sam"
        paired_manifest = work / "paired.manifest.json"
        paired_metrics = work / "paired.metrics.txt"
        paired_summary = run(binary, paired_source, paired_output, paired_manifest, paired_metrics,
                             "--tagging-policy", "All")
        assert paired_summary["duplicate_records"] == 2
        paired_rows = records(paired_output)
        assert sum(bool(int(row[1]) & 0x400) for row in paired_rows) == 2
        paired_flags = {(row[0], int(bool(int(row[1]) & 0x40))): bool(int(row[1]) & 0x400)
                        for row in paired_rows if row[0].startswith("PAIR_") and row[0] != "PAIR_C:1:1101:102:102"}
        assert paired_flags["PAIR_A:1:1101:100:100", 1] is False
        assert paired_flags["PAIR_A:1:1101:100:100", 0] is False
        assert paired_flags["PAIR_B:1:1101:101:101", 1] is True
        assert paired_flags["PAIR_B:1:1101:101:101", 0] is True
        paired_metadata = json.loads(paired_manifest.read_text(encoding="utf-8"))
        assert paired_metadata["compatibility"]["paired_fragment_consensus"] is True
        paired_metrics_text = paired_metrics.read_text(encoding="utf-8")
        paired_metric_row = metric_table(paired_metrics)["lib1"]
        assert paired_metric_row[6] == "1"  # READ_PAIR_DUPLICATES
        assert paired_metric_row[7] == "1"  # READ_PAIR_OPTICAL_DUPLICATES
        assert "## HISTOGRAM\tjava.lang.Double" in paired_metrics_text
        assert any(line.split("\t") == ["1", "0", "0", "1"]
                   for line in paired_metrics_text.splitlines())
        assert any(line.split("\t") == ["2", "1", "1", "0"]
                   for line in paired_metrics_text.splitlines())
        # The only duplicate pair in this fixture is optical, so Picard's
        # estimator receives an empty non-optical population and emits a
        # blank ESTIMATED_LIBRARY_SIZE.
        assert paired_summary["estimated_library_size"] == 0

        multi_source = work / "multi-library.sam"
        multi_source.write_text(MULTI_LIBRARY_SAM, encoding="utf-8")
        multi_output = work / "multi-library-marked.sam"
        multi_manifest = work / "multi-library.manifest.json"
        multi_metrics = work / "multi-library.metrics.txt"
        multi_summary = run(binary, multi_source, multi_output, multi_manifest, multi_metrics,
                            "--tmp-dir", str(scratch), "--max-records-in-memory", "2")
        assert multi_summary["duplicate_records"] == 2
        metric_rows = [line.split("\t") for line in multi_metrics.read_text(encoding="utf-8").splitlines()
                       if line.startswith("lib")]
        assert {row[0] for row in metric_rows} == {"libA", "libB"}
        assert all(row[5] == "1" for row in metric_rows)
        multi_metadata = json.loads(multi_manifest.read_text(encoding="utf-8"))
        assert multi_metadata["compatibility"]["per_library_metrics"] is True
        assert multi_metadata["telemetry"]["library_count"] == 2

        giant_source = work / "giant-group.sam"
        giant_source.write_text(GIANT_GROUP_SAM, encoding="utf-8")
        giant_output = work / "giant-group-marked.sam"
        giant_manifest = work / "giant-group.manifest.json"
        giant_metrics = work / "giant-group.metrics.txt"
        giant_summary = run(binary, giant_source, giant_output, giant_manifest, giant_metrics,
                            "--tmp-dir", str(scratch), "--max-records-in-memory", "2")
        assert giant_summary["duplicate_records"] == 6
        assert giant_summary["spill_run_count"] >= 4
        giant_rows = records(giant_output)
        assert sum(bool(int(row[1]) & 0x400) for row in giant_rows) == 6
        giant_metadata = json.loads(giant_manifest.read_text(encoding="utf-8"))
        assert giant_metadata["compatibility"]["chunked_group_spill"] is True
        assert giant_metadata["telemetry"]["external_group_record_limit"] == 2

        # A preempted second pass must leave a durable first-pass checkpoint
        # and spill runs that a retry can consume without rereading the input.
        resume_output = work / "resumed.sam"
        resume_manifest = work / "resumed.manifest.json"
        resume_metrics = work / "resumed.metrics.txt"
        resume_scratch = work / "resume-scratch"
        interrupted = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(resume_output),
             "--metrics-file", str(resume_metrics),
             "--output-manifest", str(resume_manifest),
             "--tmp-dir", str(resume_scratch), "--max-records-in-memory", "2",
             "--tagging-policy", "All"],
            text=True, capture_output=True, check=False,
            env={**os.environ, "FASTGATK_MARKDUPLICATES_FAIL_AFTER_CHECKPOINT": "1"},
        )
        assert interrupted.returncode != 0
        assert "checkpoint interruption" in interrupted.stderr
        checkpoints = list(resume_scratch.glob("fastgatk-markduplicates-checkpoint-*.txt"))
        assert len(checkpoints) == 1
        assert list(resume_scratch.glob("fastgatk-markduplicates-run-*.bin"))
        resumed_summary = run(binary, source, resume_output, resume_manifest, resume_metrics,
                              "--tmp-dir", str(resume_scratch), "--max-records-in-memory", "2",
                              "--tagging-policy", "All", "--resume-spill")
        assert resumed_summary["resumed"] is True
        assert semantic_records(marked) == semantic_records(resume_output)
        assert metric_table(marked_metrics) == metric_table(resume_metrics)
        resumed_metadata = json.loads(resume_manifest.read_text(encoding="utf-8"))
        assert resumed_metadata["compatibility"]["spill_checkpoint"] is True
        assert resumed_metadata["compatibility"]["resumed"] is True
        assert resumed_metadata["telemetry"]["resumed"] is True
        assert not list(resume_scratch.glob("fastgatk-markduplicates-checkpoint-*.txt"))
        assert not list(resume_scratch.glob("fastgatk-markduplicates-run-*.bin"))

        indexed = work / "indexed.bam"
        indexed_manifest = work / "indexed.manifest.json"
        indexed_metrics = work / "indexed.metrics.txt"
        indexed_summary = run(binary, source, indexed, indexed_manifest, indexed_metrics,
                              "--create-output-bam-index=true")
        assert indexed_summary["output_records"] == 3
        assert indexed.exists() and Path(f"{indexed}.bai").exists()

        removed = work / "removed.sam"
        removed_manifest = work / "removed.manifest.json"
        removed_metrics = work / "removed.metrics.txt"
        removed_summary = run(binary, source, removed, removed_manifest, removed_metrics,
                              "--remove-duplicates")
        assert removed_summary["output_records"] == 2
        assert removed.exists()

        # When the pinned GATK/Picard oracle is available, compare the actual
        # duplicate flag and DT tag semantics on both an unpaired and a paired
        # fixture.  Header command-line provenance and compression are not
        # part of this semantic check.
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if java.is_file() and gatk_jar.is_file():
            oracle_output = work / "picard-oracle.sam"
            oracle_metrics = work / "picard-oracle.metrics"
            oracle = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "MarkDuplicates",
                "-I", str(source), "-O", str(oracle_output), "-M", str(oracle_metrics),
                "--ASSUME_SORTED", "true", "--CREATE_INDEX", "false",
                "--TAGGING_POLICY", "All",
            ], text=True, capture_output=True, check=False)
            assert oracle.returncode == 0, oracle.stderr
            oracle_native = work / "native-oracle.sam"
            oracle_native_metrics = work / "native-oracle.metrics"
            native_oracle = run(binary, source, oracle_native, work / "native-oracle.manifest.json",
                                 oracle_native_metrics, "--tagging-policy", "All",
                                 "--add-pg-tag=true", "--create-output-bam-index=false")
            assert native_oracle["duplicate_records"] == 1
            assert semantic_records(oracle_output) == semantic_records(oracle_native)
            oracle_metric = metric_table(oracle_metrics)["Unknown Library"]
            native_metric = metric_table(oracle_native_metrics)["Unknown Library"]
            assert oracle_metric[:9] == native_metric[:9]
            assert bool(oracle_metric[9]) == bool(native_metric[9])

            paired_oracle_output = work / "picard-paired-oracle.sam"
            paired_oracle_metrics = work / "picard-paired-oracle.metrics"
            paired_oracle = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "MarkDuplicates",
                "-I", str(paired_source), "-O", str(paired_oracle_output), "-M", str(paired_oracle_metrics),
                "--ASSUME_SORTED", "true", "--CREATE_INDEX", "false",
                "--TAGGING_POLICY", "All",
            ], text=True, capture_output=True, check=False)
            assert paired_oracle.returncode == 0, paired_oracle.stderr
            paired_native_oracle = work / "native-paired-oracle.sam"
            paired_native_oracle_metrics = work / "native-paired-oracle.metrics"
            run(binary, paired_source, paired_native_oracle,
                work / "native-paired-oracle.manifest.json", paired_native_oracle_metrics,
                "--tagging-policy", "All", "--add-pg-tag=true",
                "--create-output-bam-index=false")
            assert semantic_records(paired_oracle_output) == semantic_records(paired_native_oracle)
            oracle_paired_metric = metric_table(paired_oracle_metrics)["lib1"]
            native_paired_metric = metric_table(paired_native_oracle_metrics)["lib1"]
            assert oracle_paired_metric[:9] == native_paired_metric[:9]
            assert bool(oracle_paired_metric[9]) == bool(native_paired_metric[9])

    print(json.dumps({"status": "pass", "duplicate_records": summary["duplicate_records"],
                      "optical_duplicate_records": summary["optical_duplicate_records"],
                      "removed_output_records": removed_summary["output_records"], "indexed": True,
                      "picard_oracle": (java.is_file() and gatk_jar.is_file())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
