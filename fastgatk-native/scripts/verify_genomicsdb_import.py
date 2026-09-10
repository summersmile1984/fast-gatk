#!/usr/bin/env python3
"""Verify the resource-aware GenomicsDBImport external-backend adapter."""
from __future__ import annotations

import json
import gzip
import os
import stat
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    binary = build / "fastgatk-genomicsdb-import"
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-import-") as temp:
        work = Path(temp)
        fake_backend = work / "fake-gatk"
        fake_backend.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, sys\n"
            "args = sys.argv[1:]\n"
            "assert args and args[0] == 'GenomicsDBImport'\n"
            "workspace = pathlib.Path(args[args.index('--genomicsdb-workspace-path') + 1])\n"
            "workspace.mkdir(parents=True, exist_ok=True)\n"
            "(workspace / 'callset.json').write_text('{}\\n')\n"
            "(workspace / 'argv.txt').write_text('\\n'.join(args))\n"
            "sys.exit(0)\n",
            encoding="utf-8",
        )
        fake_backend.chmod(fake_backend.stat().st_mode | stat.S_IXUSR)
        scratch = work / "scratch"
        scratch.mkdir()
        input_vcf = work / "sample-1.g.vcf.gz"
        with gzip.open(input_vcf, "wt", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##contig=<ID=chr1,length=100>\n")
            handle.write("##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n")
            handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
            handle.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n")
            handle.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n")
            handle.write("##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum depth>\n")
            handle.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=1\tGT:DP:GQ:MIN_DP:PL\t0/0:10:99:10:0,99,990\n")
        input_vcf_chr2 = work / "sample-2.g.vcf.gz"
        with gzip.open(input_vcf_chr2, "wt", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##contig=<ID=chr2,length=100>\n")
            handle.write("##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n")
            handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
            handle.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n")
            handle.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n")
            handle.write("##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum depth>\n")
            handle.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n")
            handle.write("chr2\t1\t.\tC\t<NON_REF>\t.\t.\tEND=1\tGT:DP:GQ:MIN_DP:PL\t0/0:7:99:7:0,99,990\n")
        workspace = work / "workspace"
        manifest = work / "workspace.manifest.json"
        result = subprocess.run(
            [
                str(binary),
                "-V", str(input_vcf),
                "--genomicsdb-workspace-path", str(workspace),
                "--batch-size", "999",
                "--reader-threads", "999",
                "--output-manifest", str(manifest),
                "--release-specific-option", "preserved",
            ],
            text=True,
            capture_output=True,
            env={
                **os.environ,
                "FASTGATK_GATK_BINARY": str(fake_backend),
                "SLURM_JOB_ID": "adapter-test",
                "SLURM_CPUS_PER_TASK": "4",
                "SLURM_MEM_PER_NODE": "64",
                "SLURM_TMPDIR": str(scratch),
            },
            check=False,
        )
        assert result.returncode == 0, (result.stdout, result.stderr)
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["status"] == "adapter"
        assert 0 < summary["effective_batch_size"] < 999
        assert 0 < summary["effective_reader_threads"] <= 4
        assert workspace.is_dir() and (workspace / "callset.json").stat().st_size > 0
        assert (workspace / "fastgatk-inputs.tsv").is_file()
        assert str(input_vcf.resolve()) in (workspace / "fastgatk-inputs.tsv").read_text(encoding="utf-8")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["status"] == "adapter"
        assert metadata["execution_mode"] == "external-genomicsdb"
        assert metadata["outputs"][0]["complete"] is True
        assert metadata["compatibility"]["external_backend"] is True
        assert metadata["compatibility"]["workspace_input_index"] is True
        assert metadata["telemetry"]["workspace_input_count"] == 1
        assert metadata["telemetry"]["tmp_dir"] == str(scratch)
        backend_argv = (workspace / "argv.txt").read_text(encoding="utf-8").splitlines()
        assert "--release-specific-option" in backend_argv
        assert "preserved" in backend_argv
        assert str(scratch) in backend_argv
        assert "--output-manifest" not in backend_argv
        assert int(backend_argv[backend_argv.index("--batch-size") + 1]) == summary["effective_batch_size"]
        assert int(backend_argv[backend_argv.index("--reader-threads") + 1]) == summary["effective_reader_threads"]

        # Reusing a populated workspace must be explicit.  This prevents a
        # failed/retried Nextflow task from silently appending into a stale
        # TileDB workspace; the GATK option remains available when overwrite
        # is intentional.
        refused = subprocess.run(
            [str(binary), "-V", str(input_vcf),
             "--genomicsdb-workspace-path", str(workspace)],
            text=True, capture_output=True,
            env={**os.environ, "FASTGATK_GATK_BINARY": str(fake_backend)},
            check=False,
        )
        assert refused.returncode != 0 and "non-empty" in refused.stderr
        overwritten = subprocess.run(
            [str(binary), "-V", str(input_vcf),
             "--genomicsdb-workspace-path", str(workspace),
             "--overwrite-existing-genomicsdb-workspace"],
            text=True, capture_output=True,
            env={**os.environ, "FASTGATK_GATK_BINARY": str(fake_backend)},
            check=False,
        )
        assert overwritten.returncode == 0, overwritten.stderr

        # Two retried tasks targeting the same workspace must not write at
        # once.  The first process holds the atomic sibling-directory lock
        # while its backend is running; the second fails closed with a
        # resource error, and the lock is removed after the first completes.
        slow_backend = work / "slow-gatk"
        slow_backend.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, sys, time\n"
            "args = sys.argv[1:]\n"
            "workspace = pathlib.Path(args[args.index('--genomicsdb-workspace-path') + 1])\n"
            "time.sleep(0.8)\n"
            "workspace.mkdir(parents=True, exist_ok=True)\n"
            "(workspace / 'callset.json').write_text('{}\\n')\n"
            "sys.exit(0)\n",
            encoding="utf-8",
        )
        slow_backend.chmod(slow_backend.stat().st_mode | stat.S_IXUSR)
        concurrent_workspace = work / "concurrent-workspace"
        concurrent_manifest = work / "concurrent.manifest.json"
        first = subprocess.Popen(
            [str(binary), "-V", str(input_vcf),
             "--genomicsdb-workspace-path", str(concurrent_workspace),
             "--output-manifest", str(concurrent_manifest)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "FASTGATK_GATK_BINARY": str(slow_backend)},
        )
        lock_path = Path(f"{concurrent_workspace}.fastgatk.lock")
        for _ in range(80):
            if lock_path.is_dir():
                break
            time.sleep(0.02)
        assert lock_path.is_dir(), "first GenomicsDB writer did not acquire its lock"
        second = subprocess.run(
            [str(binary), "-V", str(input_vcf),
             "--genomicsdb-workspace-path", str(concurrent_workspace)],
            text=True, capture_output=True,
            env={**os.environ, "FASTGATK_GATK_BINARY": str(slow_backend)},
            check=False,
        )
        assert second.returncode != 0 and "workspace is locked" in second.stderr
        first_stdout, first_stderr = first.communicate(timeout=10)
        assert first.returncode == 0, (first_stdout, first_stderr)
        assert not lock_path.exists(), "workspace lock was not released after publish"
        concurrent_metadata = json.loads(concurrent_manifest.read_text(encoding="utf-8"))
        assert concurrent_metadata["compatibility"]["workspace_lock"] is True
        assert concurrent_metadata["telemetry"]["workspace_lock_path"] == str(lock_path)

        # The native sparse-index backend is an explicit no-Java mode for
        # native GenotypeGVCFs.  It publishes the same portable input index
        # consumed by gendb:// without pretending to be a TileDB workspace.
        native_workspace = work / "native-workspace"
        native_manifest = work / "native-workspace.manifest.json"
        native_result = subprocess.run([
            str(binary), "-V", str(input_vcf),
            "--genomicsdb-workspace-path", str(native_workspace),
            "--output-manifest", str(native_manifest),
            "--fastgatk-native-workspace",
        ], text=True, capture_output=True, env=os.environ.copy(), check=False)
        assert native_result.returncode == 0, native_result.stderr
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        assert native_summary["execution_mode"] == "native-sparse-index"
        native_metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert native_metadata["execution_mode"] == "native-sparse-index"
        assert native_metadata["compatibility"]["native_workspace_index"] is True
        assert native_metadata["compatibility"]["native_workspace_materialized"] is True
        assert native_metadata["compatibility"]["external_backend"] is False
        assert native_metadata["compatibility"]["bit_identical_to_gatk"] is False
        assert native_metadata["compatibility"]["tiledb_storage"] is False
        assert native_metadata["compatibility"]["genomicsdb_query_compatible"] is False
        assert native_metadata["compatibility"]["storage_format"] == "fastgatk-portable-sparse-v1"
        assert (native_workspace / "fastgatk-workspace.json").is_file()
        assert (native_workspace / "fastgatk-native-inputs.tsv").is_file()
        record_index = native_workspace / "fastgatk-record-index.tsv"
        assert record_index.is_file()
        record_rows = [line for line in record_index.read_text(encoding="utf-8").splitlines()
                       if line and not line.startswith("#")]
        assert len(record_rows) == 1 and "\tchr1\t0\t1\t1" in record_rows[0]
        native_workspace_metadata = json.loads(
            (native_workspace / "fastgatk-workspace.json").read_text(encoding="utf-8"))
        assert native_workspace_metadata["record_index"] == "fastgatk-record-index.tsv"
        assert native_workspace_metadata["record_index_schema"] == 1
        assert native_workspace_metadata["storage_format"] == "fastgatk-portable-sparse-v1"
        assert native_workspace_metadata["tiledb_schema"] is False
        assert native_workspace_metadata["genomicsdb_query_compatible"] is False
        assert native_workspace_metadata["query_engine"] == "fastgatk-native-genotype-gvcf"
        assert list((native_workspace / "native-inputs").glob("input-*"))

        # GATK preserves -V/sample-map order in the GenomicsDB callset map.
        # A filename sort would be invisible for the usual sample-1/sample-2
        # fixture but would reorder downstream genotype columns for a reverse
        # invocation.  The native workspace must retain the explicit order
        # while still deduplicating repeated paths.
        native_workspace_reversed = work / "native-workspace-reversed"
        reversed_result = subprocess.run([
            str(binary), "-V", str(input_vcf_chr2), "-V", str(input_vcf),
            "-V", str(input_vcf_chr2),
            "--genomicsdb-workspace-path", str(native_workspace_reversed),
            "--output-manifest", str(work / "native-workspace-reversed.manifest.json"),
            "--fastgatk-native-workspace",
        ], text=True, capture_output=True, env=os.environ.copy(), check=False)
        assert reversed_result.returncode == 0, reversed_result.stderr
        reversed_index = [line for line in (
            native_workspace_reversed / "fastgatk-inputs.tsv").read_text(
                encoding="utf-8").splitlines() if line and not line.startswith("#")]
        assert reversed_index == [str(input_vcf_chr2.resolve()), str(input_vcf.resolve())]
        reversed_metadata = json.loads((native_workspace_reversed / "fastgatk-workspace.json").read_text(
            encoding="utf-8"))
        assert reversed_metadata["inputs"] == reversed_index

        # A preempted native import must leave an atomic per-input checkpoint
        # that a retried SLURM/Nextflow task can validate and reuse.  Inject a
        # failure after the first shard, then resume with the same input set;
        # the completed copy is not rewritten and the checkpoint is removed
        # only after the record index and metadata publish successfully.
        native_resume_workspace = work / "native-resume-workspace"
        native_resume_manifest = work / "native-resume.manifest.json"
        interrupted = subprocess.run([
            str(binary), "-V", str(input_vcf), "-V", str(input_vcf_chr2),
            "--genomicsdb-workspace-path", str(native_resume_workspace),
            "--fastgatk-native-workspace",
        ], text=True, capture_output=True, env={
            **os.environ, "FASTGATK_GENOMICSDB_FAIL_AFTER_INPUTS": "1"
        }, check=False)
        assert interrupted.returncode != 0
        assert "checkpoint interruption" in interrupted.stderr
        checkpoint = native_resume_workspace / "fastgatk-native-workspace.checkpoint"
        assert checkpoint.is_file() and checkpoint.stat().st_size > 0
        assert (native_resume_workspace / "native-inputs" / "input-00000000.gz").is_file()
        assert not (native_resume_workspace / "native-inputs" / "input-00000001.gz").exists()
        resumed = subprocess.run([
            str(binary), "-V", str(input_vcf), "-V", str(input_vcf_chr2),
            "--genomicsdb-workspace-path", str(native_resume_workspace),
            "--output-manifest", str(native_resume_manifest),
            "--resume-native-workspace",
        ], text=True, capture_output=True, env={
            **os.environ, "FASTGATK_GENOMICSDB_FAIL_AFTER_INPUTS": ""
        }, check=False)
        assert resumed.returncode == 0, (resumed.stdout, resumed.stderr)
        resumed_summary = json.loads(resumed.stdout.splitlines()[-1])
        assert resumed_summary["checkpoint_resumed"] is True
        resumed_metadata = json.loads(native_resume_manifest.read_text(encoding="utf-8"))
        assert resumed_metadata["compatibility"]["native_workspace_checkpoint"] is True
        assert resumed_metadata["compatibility"]["checkpoint_resumed"] is True
        assert resumed_metadata["telemetry"]["checkpoint_inputs_reused"] == 1
        assert resumed_metadata["telemetry"]["checkpoint_inputs_completed"] == 2
        assert not checkpoint.exists()
        assert len(list((native_resume_workspace / "native-inputs").glob("input-*.gz"))) == 2

        # Consume the workspace through the native GenotypeGVCFs reader.  A
        # generated sidecar is therefore tested end-to-end, not only as a
        # metadata-writing operation.  Remove the original source first to
        # prove the materialized workspace is self-contained.
        genotype = Path(os.environ.get(
            "FASTGATK_GENOTYPE_BINARY", str(build / "fastgatk-genotype-gvcf")))
        # A multi-shard native workspace must retain every sample header while
        # using the sparse record-span index to skip a shard that cannot overlap
        # the requested interval.  This is the key distinction from filtering
        # after decoding: chr2's body is never parsed for a chr1 query.
        native_workspace_multi = work / "native-workspace-multi"
        native_multi_manifest = work / "native-workspace-multi.manifest.json"
        native_multi = subprocess.run([
            str(binary), "-V", str(input_vcf), "-V", str(input_vcf_chr2),
            "--genomicsdb-workspace-path", str(native_workspace_multi),
            "--output-manifest", str(native_multi_manifest),
            "--fastgatk-native-workspace",
        ], text=True, capture_output=True, env=os.environ.copy(), check=False)
        assert native_multi.returncode == 0, native_multi.stderr
        input_vcf.unlink()
        input_vcf_chr2.unlink()

        reversed_output = work / "native-reversed-genotyped.vcf.gz"
        reversed_genotyped = subprocess.run([
            str(genotype), "-V", f"gendb://{native_workspace_reversed}", "-L", "chr1:1-1",
            "-O", str(reversed_output), "--output-manifest",
            str(work / "native-reversed-genotyped.manifest.json"),
            "--include-non-variant-sites",
        ], text=True, capture_output=True, check=False)
        assert reversed_genotyped.returncode == 0, reversed_genotyped.stderr
        with gzip.open(reversed_output, "rt", encoding="utf-8") as stream:
            reversed_lines = stream.read().splitlines()
        reversed_header = next(line for line in reversed_lines if line.startswith("#CHROM"))
        assert reversed_header.endswith("\tS2\tS1"), reversed_header
        reversed_rows = [line for line in reversed_lines if line and not line.startswith("#")]
        assert len(reversed_rows) == 1
        reversed_fields = reversed_rows[0].split("\t")
        # S2 is the chr2-only shard and must remain a missing genotype for a
        # chr1 query; S1 is the first (chr1) row in the original workspace.
        assert reversed_fields[9].startswith("./.") and reversed_fields[10].startswith("0/0"), reversed_rows[0]

        native_output = work / "native-genotyped.vcf.gz"
        native_genotype_manifest = work / "native-genotyped.manifest.json"
        native_genotyped = subprocess.run([
            str(genotype), "-V", f"gendb://{native_workspace}", "-L", "chr1:1-1",
            "-O", str(native_output), "--output-manifest", str(native_genotype_manifest),
        ], text=True, capture_output=True, check=False)
        assert native_genotyped.returncode == 0, native_genotyped.stderr
        assert native_output.is_file() and Path(f"{native_output}.tbi").is_file()
        native_genotype_metadata = json.loads(native_genotype_manifest.read_text(encoding="utf-8"))
        assert native_genotype_metadata["compatibility"]["native_record_index"] is True
        assert native_genotype_metadata["telemetry"]["native_record_indexed_inputs"] == 1
        assert native_genotype_metadata["telemetry"]["native_record_index_skipped_inputs"] == 0

        multi_output = work / "native-multi-genotyped.vcf.gz"
        multi_manifest = work / "native-multi-genotyped.manifest.json"
        multi_genotyped = subprocess.run([
            str(genotype), "-V", f"gendb://{native_workspace_multi}", "-L", "chr1:1-1",
            "-O", str(multi_output), "--output-manifest", str(multi_manifest),
            "--include-non-variant-sites",
        ], text=True, capture_output=True, check=False)
        assert multi_genotyped.returncode == 0, multi_genotyped.stderr
        multi_metadata = json.loads(multi_manifest.read_text(encoding="utf-8"))
        assert multi_metadata["telemetry"]["native_record_indexed_inputs"] == 2
        assert multi_metadata["telemetry"]["native_record_index_skipped_inputs"] == 1
        with gzip.open(multi_output, "rt", encoding="utf-8") as stream:
            multi_lines = stream.read().splitlines()
        chrom_header = next(line for line in multi_lines if line.startswith("#CHROM"))
        assert chrom_header.endswith("\tS1\tS2"), chrom_header
        assert [line for line in multi_lines if line and not line.startswith("#")]
        assert all(not line.startswith("chr2\t") for line in multi_lines)

        # The inverse query exercises global RID remapping: the first shard
        # is skipped, while the chr2-only shard is emitted under the merged
        # output dictionary with its S2 sample column intact.
        multi_chr2_output = work / "native-multi-chr2.vcf.gz"
        multi_chr2_manifest = work / "native-multi-chr2.manifest.json"
        multi_chr2 = subprocess.run([
            str(genotype), "-V", f"gendb://{native_workspace_multi}", "-L", "chr2:1-1",
            "-O", str(multi_chr2_output), "--output-manifest", str(multi_chr2_manifest),
            "--include-non-variant-sites",
        ], text=True, capture_output=True, check=False)
        assert multi_chr2.returncode == 0, multi_chr2.stderr
        multi_chr2_metadata = json.loads(multi_chr2_manifest.read_text(encoding="utf-8"))
        assert multi_chr2_metadata["telemetry"]["native_record_indexed_inputs"] == 2
        assert multi_chr2_metadata["telemetry"]["native_record_index_skipped_inputs"] == 1
        with gzip.open(multi_chr2_output, "rt", encoding="utf-8") as stream:
            multi_chr2_lines = stream.read().splitlines()
        assert any(line.startswith("##contig=<ID=chr2") for line in multi_chr2_lines)
        chr2_rows = [line for line in multi_chr2_lines if line.startswith("chr2\t")]
        assert len(chr2_rows) == 1
        row_fields = chr2_rows[0].split("\t")
        assert len(row_fields) == 11 and row_fields[0] == "chr2"
        assert row_fields[9].startswith("./.") and row_fields[9].split(":")[3] == "."
        assert row_fields[10].startswith("0/0") and row_fields[10].split(":")[3] == "7"

        remote_native = subprocess.run([
            str(binary), "-V", "https://example.invalid/sample.g.vcf.gz",
            "--genomicsdb-workspace-path", str(work / "remote-native-workspace"),
            "--fastgatk-native-workspace",
        ], text=True, capture_output=True, env=os.environ.copy(), check=False)
        assert remote_native.returncode != 0
        assert "BACKEND_UNAVAILABLE" in remote_native.stderr

        dry = subprocess.run(
            [
                str(root / "fastgatk-native/dispatcher/fastgatk"),
                "--dry-run", "GenomicsDBImport", "-V", "sample.g.vcf.gz",
                "--genomicsdb-workspace-path", str(workspace),
                "--release-specific-option", "preserved",
            ],
            text=True,
            capture_output=True,
            check=False,
            env={**os.environ, "FASTGATK_GENOMICSDBIMPORT_BINARY": str(binary)},
        )
        assert dry.returncode == 0, dry.stderr
        dry_plan = json.loads(dry.stdout.splitlines()[-1])
        assert dry_plan["execution_mode"] == "adapter"
        assert "--release-specific-option" in dry_plan["argv"]

    print(json.dumps({"status": "pass", "effective_batch": summary["effective_batch_size"],
                      "effective_reader_threads": summary["effective_reader_threads"],
                      "native_input_order_preserved": True,
                      "native_reversed_query_sample_order": ["S2", "S1"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
