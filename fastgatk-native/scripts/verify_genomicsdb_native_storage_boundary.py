#!/usr/bin/env python3
"""Pin the native sparse workspace's non-TileDB storage boundary.

The bundled GATK library provides a release-coupled query API and JNI importer,
but no matching public importer/schema headers.  The no-Java workspace must
therefore remain a portable VCF materialization plus sparse span index.  This
oracle proves that its metadata says so and that neither the isolated native
GenomicsDB bridge nor GATK 4.6.2.0 mistake it for a real TileDB workspace.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    importer = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_IMPORT_BINARY", build / "fastgatk-genomicsdb-import"))
    bridge = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_BRIDGE_BINARY", build / "fastgatk-genomicsdb-export"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    gvcf = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/tiny.g.vcf"
    reference = root / "reference.fa"
    required = (importer, bridge, java, jar, gvcf, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genomicsdb_native_storage_boundary.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GENOMICSDB_BRIDGE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"missing native storage-boundary assets: {missing}")
        print(json.dumps({"status": "skipped", "reason": "storage-boundary assets unavailable"}))
        return 0

    region = "20:69491-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-native-boundary-") as directory:
        work = Path(directory)
        workspace = work / "native-workspace"
        manifest = work / "native.manifest.json"
        imported = subprocess.run([
            str(importer), "-V", str(gvcf),
            "--genomicsdb-workspace-path", str(workspace),
            "--intervals", region, "--fastgatk-native-workspace",
            "--output-manifest", str(manifest),
        ], text=True, capture_output=True, check=False)
        assert imported.returncode == 0, imported.stderr[-4000:]

        adapter_metadata = json.loads(manifest.read_text(encoding="utf-8"))
        workspace_metadata = json.loads(
            (workspace / "fastgatk-workspace.json").read_text(encoding="utf-8"))
        compatibility = adapter_metadata["compatibility"]
        assert compatibility["storage_format"] == "fastgatk-portable-sparse-v1"
        assert compatibility["tiledb_storage"] is False
        assert compatibility["genomicsdb_query_compatible"] is False
        assert workspace_metadata["tiledb_schema"] is False
        assert workspace_metadata["genomicsdb_query_compatible"] is False
        assert not (workspace / "callset.json").exists()
        assert not (workspace / "vidmap.json").exists()
        assert not any(workspace.glob("*/__array_schema.tdb"))

        bridge_result = subprocess.run(
            [str(bridge), str(workspace), str(work / "bridge.tsv")],
            text=True, capture_output=True, check=False)
        assert bridge_result.returncode != 0
        assert "no queryable arrays" in bridge_result.stderr

        gatk_result = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{workspace}", "-L", region,
            "-O", str(work / "gatk.vcf"), "--create-output-variant-index", "false",
        ], text=True, capture_output=True, check=False)
        assert gatk_result.returncode != 0, "GATK unexpectedly accepted a non-TileDB native workspace"

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "storage_format": "fastgatk-portable-sparse-v1",
        "bridge_rejected": True,
        "gatk_rejected": True,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
