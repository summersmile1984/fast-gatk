#!/usr/bin/env python3
"""Verify the native LearnReadOrientationModel strand-mixture contract."""

from __future__ import annotations

import json
import gzip
import io
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
import oracle_guard


def prior_rows(path: Path) -> dict[str, list[str]]:
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == 1
        text = archive.extractfile(members[0]).read().decode("utf-8")
    lines = [line for line in text.splitlines() if line]
    assert lines[0].startswith("#<METADATA>SAMPLE=")
    assert lines[1].split("\t")[:4] == ["context", "rev_comp", "f1r2_a", "f1r2_c"]
    rows = [line.split("\t") for line in lines[2:]]
    assert len(rows) == 64
    assert all(len(row) == 16 for row in rows)
    return {row[0]: row for row in rows}


def prior_context_order(path: Path) -> list[str]:
    """Return the serialized row order (GATK writes a HashMap value view)."""
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == 1
        text = archive.extractfile(members[0]).read().decode("utf-8")
    return [line.split("\t", 1)[0] for line in text.splitlines()[2:] if line]


def prior_lines(path: Path) -> list[str]:
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == 1
        text = archive.extractfile(members[0]).read().decode("utf-8")
    return [line for line in text.splitlines() if line][2:]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_LEARN_ORIENTATION_BINARY",
        str(Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-learn-read-orientation-model")))
    with tempfile.TemporaryDirectory(prefix="fastgatk-orientation-") as directory:
        work = Path(directory)
        shard1 = work / "shard-1.f1r2.tsv"
        shard2 = work / "shard-2.f1r2.tsv"
        header = "# FASTGATK-MUTECT2-F1R2 v2\n# contig\tposition\talt\tf1r2\tr1f2\ttotal\n"
        shard1.write_text(header + "17\t100\tA\t9\t1\t10\n17\t101\tC\t0\t4\t4\n", encoding="utf-8")
        shard2.write_text(header + "17\t102\tA\t3\t0\t3\n17\t103\tT\t1\t7\t8\n", encoding="utf-8")
        output = work / "artifact-prior.tar.gz"
        manifest = work / "artifact-prior.manifest.json"
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        result = json.loads(subprocess.check_output([
            str(native), "-I", str(shard1), "-I", str(shard2), "-O", str(output),
            "--sample", "TUMOR", "--threads", "2", "--output-manifest", str(manifest)
        ], text=True, env=env).splitlines()[-1])
        assert result == {
            "tool": "LearnReadOrientationModel", "status": "prototype", "sample": "TUMOR",
            "records": 4, "contexts": 64, "f1r2_alt_support": 13, "r1f2_alt_support": 12,
        }
        assert output.exists() and output.stat().st_size > 0
        with tarfile.open(output, "r:gz") as archive:
            members = archive.getmembers()
            assert len(members) == 1 and members[0].name == "TUMOR.orientation_priors"
            text = archive.extractfile(members[0]).read().decode("utf-8")
        lines = [line for line in text.splitlines() if line]
        assert lines[0] == "#<METADATA>SAMPLE=TUMOR"
        assert lines[1].split("\t")[:4] == ["context", "rev_comp", "f1r2_a", "f1r2_c"]
        rows = [line.split("\t") for line in lines[2:]]
        assert len(rows) == 64
        for row in rows:
            assert len(row) == 16
            probabilities = [float(value) for value in row[2:14]]
            assert abs(sum(probabilities) - 1.0) < 1e-9
            assert row[14:] == ["4", "4"]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["tar_gz"] is True
        assert metadata["compatibility"]["artifact_prior_table"] is True
        assert metadata["compatibility"]["full_gatk_em"] is False
        assert metadata["telemetry"]["f1r2_alt_support"] == 13
        telemetry = metadata["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_execution_space"]
        assert telemetry["kernel_batches"] > 0
        assert telemetry["kernel_observations"] == 4
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0

        # Legacy Mutect2 sidecars may be compressed while they are staged
        # between workflow tasks.  The compressed reader must feed the same
        # fixed-order aggregation as the plain TSV reader.
        compressed_shard = work / "shard-1.f1r2.tsv.gz"
        with gzip.open(compressed_shard, "wt", encoding="utf-8") as stream:
            stream.write(shard1.read_text(encoding="utf-8"))
        compressed_output = work / "artifact-prior-compressed.tsv.tar.gz"
        compressed_result = json.loads(subprocess.check_output([
            str(native), "-I", str(compressed_shard), "-I", str(shard2),
            "-O", str(compressed_output), "--sample", "TUMOR", "--threads", "2",
        ], text=True, env=env).splitlines()[-1])
        assert compressed_result == result
        assert prior_lines(compressed_output) == prior_lines(output)

        input_archive = work / "collect-f1r2.tar.gz"
        with tarfile.open(input_archive, "w:gz") as archive:
            payload = (header + "17\t104\tG\t5\t2\t7\n").encode("utf-8")
            info = tarfile.TarInfo("TUMOR.f1r2")
            info.size = len(payload)
            archive.addfile(info, fileobj=io.BytesIO(payload))
        tar_output = work / "artifact-prior-from-tar.tar.gz"
        tar_result = json.loads(subprocess.check_output([
            str(native), "-I", str(input_archive), "-O", str(tar_output),
            "--sample", "TUMOR", "--output-manifest", str(work / "tar.manifest.json")
        ], text=True, env=env).splitlines()[-1])
        assert tar_result["records"] == 1
        assert tar_output.exists() and tar_output.stat().st_size > 0

        # Exercise the real GATK CollectF1R2Counts -> LearnReadOrientationModel
        # contract when the pinned local oracle is available.  This is the
        # standard three-member tar.gz path, not the legacy FASTGATK TSV shim.
        gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk"
        java = root / "third_party/jdk17/bin/java"
        standard_oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
        if oracle_guard.oracle_ready('verify_learn_read_orientation_model.py', gatk, java):
            oracle_env = env.copy()
            oracle_env["JAVA_HOME"] = str(java.parent.parent)
            oracle_env["PATH"] = str(java.parent) + os.pathsep + oracle_env.get("PATH", "")
            collect_archive = work / "gatk-collect-f1r2.tar.gz"
            collect = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "CollectF1R2Counts",
                "-R", str(root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
                "-I", str(root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"),
                "-L", "17:69000-70000", "-O", str(collect_archive),
            ], env=oracle_env, text=True, capture_output=True)
            assert collect.returncode == 0, collect.stderr[-4000:]
            gatk_prior = work / "gatk-standard-prior.tar.gz"
            native_prior = work / "native-standard-prior.tar.gz"
            gatk_run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "LearnReadOrientationModel",
                "-I", str(collect_archive), "-O", str(gatk_prior), "--max-depth", "200",
            ], env=oracle_env, text=True, capture_output=True)
            assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
            native_run = subprocess.run([
                str(native), "-I", str(collect_archive), "-O", str(native_prior),
                "--sample", "NA12878", "--threads", "2",
                "--output-manifest", str(work / "native-standard.manifest.json"),
            ], env=env, text=True, capture_output=True)
            assert native_run.returncode == 0, native_run.stderr[-4000:]
            native_summary = json.loads(native_run.stdout.splitlines()[-1])
            assert native_summary["contexts"] == 64 and native_summary["contexts_with_data"] > 0
            native_manifest = json.loads((work / "native-standard.manifest.json").read_text(encoding="utf-8"))
            assert native_manifest["compatibility"]["full_gatk_em"] is True
            assert native_manifest["compatibility"]["em_kokkos"] is True
            assert native_manifest["telemetry"]["em_kokkos_observations"] > 0
            assert native_manifest["telemetry"]["em_kokkos_batches"] > 0
            assert native_manifest["telemetry"]["em_execution_space"]
            assert native_manifest["telemetry"]["kernel_lifecycle"] == (
                "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
            )
            assert native_manifest["telemetry"]["kernel_execution_policy"] == "RangePolicy"
            assert native_manifest["telemetry"]["kernel_batches"] > 0
            assert native_manifest["telemetry"]["kernel_observations"] > 0
            gatk_rows = prior_rows(gatk_prior)
            native_rows = prior_rows(native_prior)
            # ArtifactPriorCollection serializes HashMap.values(), whose
            # bucket order is observable in the GATK artifact.  Keep this
            # ordering contract explicit even though downstream lookup is by
            # context.  Rows with no observations should also be byte-stable;
            # any remaining differences are then confined to the EM math.
            assert prior_context_order(gatk_prior) == prior_context_order(native_prior)
            for expected_line, actual_line in zip(prior_lines(gatk_prior), prior_lines(native_prior)):
                if expected_line.split("\t")[-2:] == ["0", "0"]:
                    assert expected_line == actual_line
            max_delta = 0.0
            count_mismatches = 0
            for context, expected in gatk_rows.items():
                actual = native_rows[context]
                for index in range(2, 14):
                    max_delta = max(max_delta, abs(float(expected[index]) - float(actual[index])))
                if expected[14:] != actual[14:]:
                    count_mismatches += 1
            assert max_delta <= 1e-10 and count_mismatches == 0
            standard_oracle = {
                "status": "pass", "contexts": 64, "contexts_with_data": native_summary["contexts_with_data"],
                "max_probability_delta": max_delta, "count_mismatches": count_mismatches,
            }
        print(json.dumps({"status": "pass", "records": 4, "contexts": 64,
                          "archive_members": 1, "f1r2_alt_support": 13,
                          "r1f2_alt_support": 12, "tar_input": True,
                          "standard_gatk_oracle": standard_oracle}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
