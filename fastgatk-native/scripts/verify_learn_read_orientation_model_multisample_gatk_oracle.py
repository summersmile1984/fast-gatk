#!/usr/bin/env python3
"""Pin LearnReadOrientationModel's multi-sample tar writer to GATK 4.6.2.0."""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
import oracle_guard


def clone_samples(source: Path, destination: Path, samples: tuple[str, ...]) -> None:
    """Clone a valid CollectF1R2Counts archive with distinct sample headers."""
    with tarfile.open(source, "r:gz") as archive, tarfile.open(destination, "w:gz") as output:
        members = archive.getmembers()
        for sample in samples:
            for member in members:
                payload = archive.extractfile(member).read()
                payload = payload.replace(b"NA12878", sample.encode("ascii"))
                copied = tarfile.TarInfo(member.name.replace("NA12878", sample))
                copied.mode = member.mode
                copied.uid = member.uid
                copied.gid = member.gid
                copied.mtime = member.mtime
                copied.size = len(payload)
                output.addfile(copied, io.BytesIO(payload))


def members_by_sample(path: Path) -> dict[str, bytes]:
    with tarfile.open(path, "r:gz") as archive:
        result: dict[str, bytes] = {}
        for member in archive.getmembers():
            name = member.name.removeprefix("./")
            assert name.endswith(".orientation_priors"), name
            result[name] = archive.extractfile(member).read()
        return result


def assert_prior_equivalent(expected: bytes, actual: bytes) -> float:
    expected_lines = expected.decode("utf-8").splitlines()
    actual_lines = actual.decode("utf-8").splitlines()
    assert expected_lines[:2] == actual_lines[:2]
    assert len(expected_lines) == len(actual_lines) == 66
    max_delta = 0.0
    for expected_line, actual_line in zip(expected_lines[2:], actual_lines[2:]):
        left, right = expected_line.split("\t"), actual_line.split("\t")
        assert left[:2] == right[:2] and left[14:] == right[14:]
        for index in range(2, 14):
            max_delta = max(max_delta, abs(float(left[index]) - float(right[index])))
    # GATK and native use binary64 special functions; the pinned contract is
    # the same tolerance as the existing single-sample oracle, not byte text.
    assert max_delta <= 1e-10, max_delta
    return max_delta


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_LEARN_ORIENTATION_BINARY",
        str(root / "fastgatk-native/build/fastgatk-learn-read-orientation-model")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk"
    java = root / "third_party/jdk17/bin/java"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    assert native.exists() and reference.exists() and bam.exists()
    oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
    max_depth_delta = None
    with tempfile.TemporaryDirectory(prefix="fastgatk-orientation-multisample-") as directory:
        work = Path(directory)
        collect = work / "collect.tar.gz"
        java_tmp = work / "java-tmp"
        native_tmp = work / "native-tmp"
        java_tmp.mkdir()
        native_tmp.mkdir()
        samples = ("SAMPLE1", "SAMPLE2")
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        if oracle_guard.oracle_ready('verify_learn_read_orientation_model_multisample_gatk_oracle.py', gatk, java):
            oracle_env = env.copy()
            oracle_env["JAVA_HOME"] = str(java.parent.parent)
            oracle_env["PATH"] = str(java.parent) + os.pathsep + oracle_env.get("PATH", "")
            collect_one = work / "collect-one.tar.gz"
            run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "CollectF1R2Counts",
                "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000", "-O", str(collect_one),
            ], env=oracle_env, text=True, capture_output=True)
            assert run.returncode == 0, run.stderr[-4000:]
            clone_samples(collect_one, collect, samples)
            gatk_output = work / "gatk-priors.tar.gz"
            run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "LearnReadOrientationModel",
                "-I", str(collect), "-O", str(gatk_output), "--max-depth", "200",
                "--QUIET", "true", "--tmp-dir", str(java_tmp),
                "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
                "--verbosity", "WARNING",
            ], env=oracle_env, text=True, capture_output=True)
            assert run.returncode == 0, run.stderr[-4000:]
            oracle = {"status": "pass", "gatk_members": sorted(members_by_sample(gatk_output))}
        native_output = work / "native-priors.tar.gz"
        manifest = work / "native.manifest.json"
        run = subprocess.run([
            str(native), "-I", str(collect), "-O", str(native_output),
            "--output-manifest", str(manifest), "--threads", "2",
            "--QUIET", "true", "--tmp-dir", str(native_tmp),
            "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
            "--verbosity", "WARNING",
        ], env=env, text=True, capture_output=True)
        assert run.returncode == 0, run.stderr[-4000:]
        summary = json.loads(run.stdout.splitlines()[-1])
        native_members = members_by_sample(native_output)
        assert sorted(native_members) == [f"{sample}.orientation_priors" for sample in samples]
        assert oracle["status"] == "skip" or native_members.keys() == set(oracle["gatk_members"])
        if oracle["status"] == "pass":
            gatk_members = members_by_sample(gatk_output)
            for sample in samples:
                assert_prior_equivalent(gatk_members[f"{sample}.orientation_priors"],
                                        native_members[f"{sample}.orientation_priors"])

            # GATK's max-depth is an EM observation bound.  It does not
            # re-bin higher-depth rows already present in a CollectF1R2Counts
            # histogram: those rows remain part of numExamples but are not
            # visited by the per-iteration E/M update.  This fixture has
            # source bins through depth 200, so max-depth=2 protects this
            # otherwise easy-to-miss data/algorithm boundary.
            gatk_bounded_output = work / "gatk-priors-max-depth-2.tar.gz"
            run = subprocess.run([
                sys.executable, str(gatk), "--java-options", "-Xmx2g", "LearnReadOrientationModel",
                "-I", str(collect), "-O", str(gatk_bounded_output), "--max-depth", "2",
                "--QUIET", "true", "--tmp-dir", str(java_tmp),
                "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
                "--verbosity", "WARNING",
            ], env=oracle_env, text=True, capture_output=True)
            assert run.returncode == 0, run.stderr[-4000:]
            native_bounded_output = work / "native-priors-max-depth-2.tar.gz"
            run = subprocess.run([
                str(native), "-I", str(collect), "-O", str(native_bounded_output),
                "--max-depth", "2", "--threads", "2", "--QUIET", "true",
                "--tmp-dir", str(native_tmp), "--use-jdk-deflater", "false",
                "--use-jdk-inflater", "false", "--verbosity", "WARNING",
            ], env=env, text=True, capture_output=True)
            assert run.returncode == 0, run.stderr[-4000:]
            bounded_gatk_members = members_by_sample(gatk_bounded_output)
            bounded_native_members = members_by_sample(native_bounded_output)
            max_depth_delta = 0.0
            for sample in samples:
                max_depth_delta = max(
                    max_depth_delta,
                    assert_prior_equivalent(
                        bounded_gatk_members[f"{sample}.orientation_priors"],
                        bounded_native_members[f"{sample}.orientation_priors"],
                    ),
                )
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["samples"] == 2
        assert metadata["compatibility"]["collect_f1r2_tar_input"] is True
        assert metadata["compatibility"]["quiet"] is True
        assert metadata["telemetry"]["quiet"] is True
        assert metadata["compatibility"]["verbosity"] == "WARNING"
        assert metadata["telemetry"]["verbosity"] == "WARNING"
        assert len(native_members) == 2
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            assert oracle["status"] == "pass"
            assert max_depth_delta is not None
        print(json.dumps({"status": "pass", "samples": 2, "members": sorted(native_members),
                          "member_content_exact": oracle["status"] == "pass",
                          "max_depth_em_boundary": max_depth_delta is not None,
                          "max_depth": 2 if max_depth_delta is not None else None,
                          "max_depth_max_probability_delta": max_depth_delta,
                          "gatk_oracle": oracle}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
