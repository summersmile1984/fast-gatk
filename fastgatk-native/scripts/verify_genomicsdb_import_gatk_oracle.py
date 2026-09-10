#!/usr/bin/env python3
"""Verify the GenomicsDBImport adapter against a pinned Java GATK backend.

This is intentionally a narrow boundary oracle: the native executable must
preserve the complete GenomicsDBImport argument surface while delegating the
actual TileDB write to GATK, and the resulting workspace must remain readable
by GATK GenotypeGVCFs.  It does not claim that the native sparse workspace is
TileDB-compatible (that mode is covered by verify_genomicsdb_import.py).
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    binary = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_IMPORT_BINARY", build / "fastgatk-genomicsdb-import"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    gvcf = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/tiny.g.vcf"
    reference = root / "reference.fa"
    required = (binary, java, jar, gvcf, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled GATK GenomicsDBImport oracle assets unavailable: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    region = "20:69491-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-import-gatk-oracle-") as directory:
        work = Path(directory)
        wrapper = work / "gatk-wrapper"
        # The adapter intentionally invokes an executable beginning with the
        # GATK tool name.  Keep this wrapper in the test so the production
        # binary remains a normal direct-replacement executable.
        wrapper.write_text(
            "#!/usr/bin/env python3\n"
            "import subprocess, sys\n"
            f"sys.exit(subprocess.run([{str(java)!r}, '-Xmx1g', '-jar', {str(jar)!r}, *sys.argv[1:]]).returncode)\n",
            encoding="utf-8",
        )
        wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR)
        workspace = work / "workspace"
        manifest = work / "adapter.manifest.json"
        result = subprocess.run([
            str(binary), "-V", str(gvcf),
            "--genomicsdb-workspace-path", str(workspace),
            "--intervals", region,
            "--batch-size", "1", "--reader-threads", "1",
            "--output-manifest", str(manifest),
        ], text=True, capture_output=True, check=False,
            env={**os.environ, "FASTGATK_GATK_BINARY": str(wrapper),
                 "FASTGATK_NATIVE_BUILD": str(build)})
        assert result.returncode == 0, result.stderr[-4000:]
        summary = json.loads(result.stdout.splitlines()[-1])
        assert summary["tool"] == "GenomicsDBImport"
        assert summary["execution_mode"] == "external-genomicsdb"
        assert summary["effective_batch_size"] == 1
        assert summary["effective_reader_threads"] == 1
        assert workspace.is_dir()
        assert (workspace / "callset.json").is_file()
        assert (workspace / "vidmap.json").is_file()
        assert (workspace / "vcfheader.vcf").is_file()
        assert any(path.is_dir() and (path / "__array_schema.tdb").is_file()
                   for path in workspace.iterdir())
        index = workspace / "fastgatk-inputs.tsv"
        assert index.is_file() and str(gvcf.resolve()) in index.read_text(encoding="utf-8")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["status"] == "adapter"
        assert metadata["compatibility"]["external_backend"] is True
        assert metadata["compatibility"]["workspace_input_index"] is True
        assert metadata["compatibility"]["bit_identical_to_gatk"] is True
        assert metadata["compatibility"]["tiledb_storage"] is True
        assert metadata["compatibility"]["genomicsdb_query_compatible"] is True
        assert metadata["compatibility"]["storage_format"] == "external-genomicsdb"
        assert metadata["telemetry"]["workspace_input_count"] == 1

        # Re-open the workspace with the same pinned GATK release.  This
        # catches adapter-side workspace publication/race regressions that a
        # fake backend contract cannot observe.
        output = work / "gatk-joint.vcf"
        reopened = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{workspace}",
            "-L", region, "-O", str(output),
            "--create-output-variant-index", "false",
        ], text=True, capture_output=True, check=False)
        assert reopened.returncode == 0, reopened.stderr[-4000:]
        rows = [line for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")]
        assert rows, "GATK could not read the adapter-published workspace"

    print(json.dumps({"status": "pass", "gatk_version": "4.6.2.0",
                      "execution_mode": "external-genomicsdb",
                      "workspace_reopen": True, "genotyped_records": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
