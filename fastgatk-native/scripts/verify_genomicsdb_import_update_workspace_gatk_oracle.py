#!/usr/bin/env python3
"""Pin GenomicsDBImport incremental-workspace semantics to GATK 4.6.2.0.

The adapter is allowed to delegate the real TileDB operation to GATK, but its
workspace path/preflight and input-index metadata must preserve the same
incremental boundary.  GATK ignores intervals on an update and reuses the
interval set from the initial import. The native sparse workspace implements
the same safe input-index boundary by atomically rebuilding its materialized
VCF/index set; it remains explicitly non-TileDB.
"""

from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run_checked(command: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False, **kwargs)
    assert result.returncode == 0, ("\n".join(command), result.stdout, result.stderr[-4000:])
    return result


def samples(path: Path) -> list[str]:
    header = next(line for line in path.read_text(encoding="utf-8").splitlines()
                  if line.startswith("#CHROM"))
    return header.split("\t")[9:]


def data_rows(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def callset_samples(workspace: Path) -> list[str]:
    callset = json.loads((workspace / "callset.json").read_text(encoding="utf-8"))
    return [entry["sample_name"] for entry in callset["callsets"]]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    importer = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_IMPORT_BINARY", build / "fastgatk-genomicsdb-import"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    gvcf = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/tiny.g.vcf"
    gvcf_index = Path(f"{gvcf}.idx")
    reference = root / "reference.fa"
    required = (importer, java, jar, gvcf, gvcf_index, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genomicsdb_import_update_workspace_gatk_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled GenomicsDB update oracle assets unavailable: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    initial_interval = "20:69491-69512"
    update_interval = "20:69762-69772"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-update-oracle-") as directory:
        work = Path(directory)
        gvcf2 = work / "sample-2.g.vcf"
        gvcf2.write_text(gvcf.read_text(encoding="utf-8").replace("SAMPLE_1", "SAMPLE_2"),
                         encoding="utf-8")
        shutil.copyfile(gvcf_index, Path(f"{gvcf2}.idx"))

        java_workspace = work / "java-workspace"
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenomicsDBImport",
            "-V", str(gvcf), "--genomicsdb-workspace-path", str(java_workspace),
            "-L", initial_interval, "--create-output-variant-index", "false",
        ])
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenomicsDBImport",
            "-V", str(gvcf2), "--genomicsdb-update-workspace-path", str(java_workspace),
            # GATK warns that update uses the initial interval set.  Keep a
            # different interval here so an adapter that accidentally strips
            # or applies it cannot pass this oracle.
            "-L", update_interval, "--create-output-variant-index", "false",
        ])
        assert callset_samples(java_workspace) == ["SAMPLE_1", "SAMPLE_2"]
        java_output = work / "java-joint.vcf"
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{java_workspace}",
            "-L", initial_interval, "-O", str(java_output),
            "--create-output-variant-index", "false",
        ])

        wrapper = work / "gatk-wrapper"
        wrapper.write_text(
            "#!/usr/bin/env python3\n"
            "import subprocess, sys\n"
            f"sys.exit(subprocess.run([{str(java)!r}, '-Xmx1g', '-jar', {str(jar)!r}, *sys.argv[1:]]).returncode)\n",
            encoding="utf-8",
        )
        wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR)

        adapter_workspace = work / "adapter-workspace"
        first_manifest = work / "adapter-first.manifest.json"
        run_checked([
            str(importer), "-V", str(gvcf),
            "--genomicsdb-workspace-path", str(adapter_workspace),
            "-L", initial_interval, "--batch-size", "1", "--reader-threads", "1",
            "--output-manifest", str(first_manifest),
        ], env={**os.environ, "FASTGATK_GATK_BINARY": str(wrapper)})

        update_manifest = work / "adapter-update.manifest.json"
        update_result = run_checked([
            str(importer), "-V", str(gvcf2),
            "--genomicsdb-update-workspace-path", str(adapter_workspace),
            "-L", update_interval, "--batch-size", "1", "--reader-threads", "1",
            "--output-manifest", str(update_manifest),
        ], env={**os.environ, "FASTGATK_GATK_BINARY": str(wrapper)})
        assert json.loads(update_result.stdout.splitlines()[-1])["external_update_workspace"] is True
        assert callset_samples(adapter_workspace) == callset_samples(java_workspace)
        indexed = [line for line in (adapter_workspace / "fastgatk-inputs.tsv").read_text(
            encoding="utf-8").splitlines() if line and not line.startswith("#")]
        assert indexed == [str(gvcf.resolve()), str(gvcf2.resolve())]
        metadata = json.loads(update_manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["external_update_workspace"] is True
        assert metadata["telemetry"]["workspace_input_count"] == 2
        assert metadata["outputs"][0]["complete"] is True

        adapter_output = work / "adapter-joint.vcf"
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{adapter_workspace}",
            "-L", initial_interval, "-O", str(adapter_output),
            "--create-output-variant-index", "false",
        ])
        assert samples(adapter_output) == samples(java_output) == ["SAMPLE_1", "SAMPLE_2"]
        # The command-line header embeds the temporary workspace path, so
        # compare the data rows byte-for-byte and the sample header separately.
        assert data_rows(adapter_output) == data_rows(java_output)

        # Native sparse update is a deterministic append/rebuild of the
        # materialized input index. It must preserve input order and remain
        # explicitly non-TileDB.
        native_workspace = work / "native-workspace"
        native_map = work / "native-sample-map.tsv"
        native_map.write_text(f"NATIVE_SAMPLE_1 {gvcf.resolve()}\n", encoding="utf-8")
        run_checked([
            str(importer), "-V", str(gvcf),
            "--sample-name-map", str(native_map),
            "--genomicsdb-workspace-path", str(native_workspace),
            "--fastgatk-native-workspace",
        ])
        native_update_manifest = work / "native-update.manifest.json"
        native_update = run_checked([
            str(importer), "-V", str(gvcf2),
            "--genomicsdb-update-workspace-path", str(native_workspace),
            "--fastgatk-native-workspace", "-L", update_interval,
            "--output-manifest", str(native_update_manifest),
        ])
        native_summary = json.loads(native_update.stdout.splitlines()[-1])
        assert native_summary["native_incremental_update"] is True
        native_index = [line for line in (native_workspace / "fastgatk-inputs.tsv").read_text(
            encoding="utf-8").splitlines() if line and not line.startswith("#")]
        assert native_index == [str(gvcf.resolve()), str(gvcf2.resolve())]
        native_metadata = json.loads((native_workspace / "fastgatk-workspace.json").read_text(
            encoding="utf-8"))
        assert native_metadata["input_count"] == 2
        assert native_metadata["incremental_update"] is True
        assert native_metadata["sample_name_map"] == [{
            "input": str(gvcf.resolve()), "sample_name": "NATIVE_SAMPLE_1"
        }]
        native_manifest = json.loads(native_update_manifest.read_text(encoding="utf-8"))
        assert native_manifest["compatibility"]["native_incremental_update"] is True
        assert native_manifest["compatibility"]["bit_identical_to_gatk"] is False

        # The no-Java GenotypeGVCFs path consumes the rebuilt index directly.
        native_genotype = Path(os.environ.get(
            "FASTGATK_GENOTYPE_BINARY", build / "fastgatk-genotype-gvcf"))
        if native_genotype.is_file():
            native_output = work / "native-joint.vcf.gz"
            native_genotype_manifest = work / "native-genotype.manifest.json"
            native_result = run_checked([
                str(native_genotype), "-R", str(reference),
                "-V", f"gendb://{native_workspace}",
                "-O", str(native_output), "--create-output-variant-index", "false",
                "--output-manifest", str(native_genotype_manifest),
            ])
            native_metadata = json.loads(native_genotype_manifest.read_text(encoding="utf-8"))
            assert native_metadata["telemetry"]["genomicsdb_workspace_input_index"] is True
            assert native_output.is_file()

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "incremental_update": True,
        "initial_interval_reused": True,
        "adapter_callset_samples_exact": True,
        "input_index_preserved": True,
        "native_sparse_incremental_rebuild": True,
        "native_workspace_remains_non_tiledb": True,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
