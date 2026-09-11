#!/usr/bin/env python3
"""Verify constrained-memory SortSam spill/merge, order and BAM index contracts."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


SAM = """@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:100
@RG\tID:rg1\tSM:S1
r3\t0\tchr1\t30\t60\t4M\t*\t0\t0\tTGCA\tIIII
r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII
r4\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t!!!!
r2\t0\tchr1\t20\t60\t4M\t*\t0\t0\tGATC\tIIII
"""


def records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def semantic_records(path: Path) -> list[tuple[str, int, str, int]]:
    return [(row[0], int(row[1]), row[2], int(row[3]))
            for row in records(path)]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-sort-sam"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-sort-sam-") as temp:
        work = Path(temp)
        source = work / "input.sam"
        output = work / "sorted.sam"
        manifest = work / "sorted.manifest.json"
        spill = work / "spill"
        source.write_text(SAM, encoding="utf-8")
        result = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(output),
             "--sort-order", "coordinate", "--max-records-in-memory", "2",
             "--tmp-dir", str(spill), "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "prototype"
        assert summary["input_records"] == 4 and summary["output_records"] == 4
        assert summary["spill_runs"] == 2
        lines = output.read_text(encoding="utf-8").splitlines()
        assert any(line == "@HD\tVN:1.6\tSO:coordinate" for line in lines)
        assert [int(record[3]) for record in records(output)[:3]] == [10, 20, 30]
        assert records(output)[3][0] == "r4"
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["external_memory_sort"] is True
        assert metadata["compatibility"]["spill_runs"] is True
        assert metadata["telemetry"]["spill_runs"] == 2
        assert all(item["complete"] for item in metadata["outputs"])

        bam_output = work / "sorted.bam"
        bam_manifest = work / "sorted.bam.manifest.json"
        bam_result = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(bam_output),
             "--max-records-in-memory", "1", "--tmp-dir", str(spill),
             "--create-output-bam-index=true",
             "--output-manifest", str(bam_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert bam_result.returncode == 0, bam_result.stderr
        assert bam_output.exists() and Path(f"{bam_output}.bai").exists()
        bam_metadata = json.loads(bam_manifest.read_text(encoding="utf-8"))
        assert any(item["kind"] == "alignment-index" and item["complete"] for item in bam_metadata["outputs"])

        queryname = work / "queryname.sam"
        query_result = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(queryname), "--sort-order", "queryname",
             "--max-records-in-memory", "2"], text=True, capture_output=True, check=False,
        )
        assert query_result.returncode == 0, query_result.stderr
        assert [record[0] for record in records(queryname)] == ["r1", "r2", "r3", "r4"]

        pg_output = work / "with-pg.sam"
        pg_manifest = work / "with-pg.manifest.json"
        pg_result = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(pg_output),
             "--max-records-in-memory", "2", "--add-pg-tag=true",
             "--program-record-id", "sort-test", "--program-group-name", "SortSam",
             "--program-group-version", "1.2", "--program-group-command-line", "fastgatk SortSam test",
             "--output-manifest", str(pg_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert pg_result.returncode == 0, pg_result.stderr
        pg_lines = pg_output.read_text(encoding="utf-8").splitlines()
        assert any(line.startswith("@PG\t") and "ID:sort-test" in line and
                   "PN:SortSam" in line and "VN:1.2" in line and
                   "CL:fastgatk SortSam test" in line for line in pg_lines)
        pg_metadata = json.loads(pg_manifest.read_text(encoding="utf-8"))
        assert pg_metadata["compatibility"]["pg_line"] is True

        invalid_pg = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(work / "invalid-pg.sam"),
             "--add-pg-tag=maybe"], text=True, capture_output=True, check=False)
        assert invalid_pg.returncode != 0 and "invalid boolean" in invalid_pg.stderr
        invalid_index = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(work / "invalid-index.bam"),
             "--create-output-bam-index=maybe"], text=True, capture_output=True, check=False)
        assert invalid_index.returncode != 0 and "invalid boolean" in invalid_index.stderr

        # An output-stage failure leaves a complete input-pass checkpoint and
        # spill runs when --resume-spill is enabled.  A retry with the same
        # contract reuses those runs instead of decoding/staging the input a
        # second time.
        resume_spill = work / "resume-spill"
        blocked_output = work / "blocked-output"
        blocked_output.mkdir()
        resume_manifest = work / "resume.manifest.json"
        failed = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(blocked_output),
             "--max-records-in-memory", "2", "--tmp-dir", str(resume_spill),
             "--resume-spill", "--output-manifest", str(resume_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert failed.returncode != 0
        checkpoint_files = list(resume_spill.glob("fastgatk-sort-checkpoint-*.txt"))
        assert len(checkpoint_files) == 1
        assert len(list(resume_spill.glob("fastgatk-sort-run-*.bam"))) == 2
        blocked_output.rmdir()
        resumed_result = subprocess.run(
            [str(binary), "-I", str(source), "-O", str(blocked_output),
             "--max-records-in-memory", "2", "--tmp-dir", str(resume_spill),
             "--resume-spill", "--output-manifest", str(resume_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert resumed_result.returncode == 0, resumed_result.stderr
        resumed_summary = json.loads(resumed_result.stdout.splitlines()[-1])
        assert resumed_summary["resumed"] is True
        resumed_metadata = json.loads(resume_manifest.read_text(encoding="utf-8"))
        assert resumed_metadata["compatibility"]["spill_checkpoint"] is True
        assert resumed_metadata["compatibility"]["resumed"] is True
        assert not list(resume_spill.glob("fastgatk-sort-checkpoint-*.txt"))
        assert not list(resume_spill.glob("fastgatk-sort-run-*.bam"))

        # Compare HTSJDK's coordinate/queryname tie-break semantics whenever
        # the pinned GATK jar is available.  This catches subtle same-locus
        # ordering differences (strand, flags, MAPQ and mate fields) that a
        # simple position-only assertion cannot detect.
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if oracle_guard.oracle_ready('verify_sort_sam.py', java, gatk_jar):
            oracle_coordinate = work / "picard-coordinate.sam"
            oracle = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "SortSam",
                "-I", str(source), "-O", str(oracle_coordinate), "-SO", "coordinate",
                "--CREATE_INDEX", "false",
            ], text=True, capture_output=True, check=False)
            assert oracle.returncode == 0, oracle.stderr
            assert semantic_records(oracle_coordinate) == semantic_records(output)
            oracle_queryname = work / "picard-queryname.sam"
            oracle = subprocess.run([
                str(java), "-Xmx1g", "-jar", str(gatk_jar), "SortSam",
                "-I", str(source), "-O", str(oracle_queryname), "-SO", "queryname",
                "--CREATE_INDEX", "false",
            ], text=True, capture_output=True, check=False)
            assert oracle.returncode == 0, oracle.stderr
            assert semantic_records(oracle_queryname) == semantic_records(queryname)

    print(json.dumps({"status": "pass", "records": summary["output_records"],
                      "spill_runs": summary["spill_runs"], "indexed": True,
                      "picard_oracle": (java.is_file() and gatk_jar.is_file())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
