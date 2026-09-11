#!/usr/bin/env python3
"""Pin GenomicsDBImport sample-name-map semantics to GATK 4.6.2.0.

The native sparse workspace is deliberately not claimed to be TileDB.  This
oracle checks the compatibility boundary that matters to a direct replacement
in a SLURM/Nextflow graph: a two-column sample map must rename the sample in
the published workspace and that name must survive a downstream
GenotypeGVCFs invocation.  The Java side supplies the release-pinned expected
callset/header shape; native output rows remain an explicitly bounded model
boundary.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run_checked(command: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False, **kwargs)
    assert result.returncode == 0, ("\n".join(command), result.stdout, result.stderr[-4000:])
    return result


def sample_from_header(path: Path) -> list[str]:
    header = next(line for line in path.read_text(encoding="utf-8").splitlines()
                  if line.startswith("#CHROM"))
    fields = header.split("\t")
    return fields[9:]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    importer = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_IMPORT_BINARY", build / "fastgatk-genomicsdb-import"))
    genotype = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY", build / "fastgatk-genotype-gvcf"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    gvcf = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/tiny.g.vcf"
    reference = root / "reference.fa"
    required = (importer, genotype, java, jar, gvcf, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genomicsdb_import_sample_map_gatk_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled sample-map oracle assets unavailable: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    interval = "20:69763-69773"
    requested_sample = "FASTGATK_SAMPLE_MAP_RENAMED"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-sample-map-oracle-") as directory:
        work = Path(directory)
        sample_map = work / "sample-map.tsv"
        sample_map.write_text(f"{requested_sample}\t{gvcf.resolve()}\n", encoding="utf-8")

        java_workspace = work / "java-workspace"
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenomicsDBImport",
            "--sample-name-map", str(sample_map),
            "--genomicsdb-workspace-path", str(java_workspace),
            "-L", interval, "--create-output-variant-index", "false",
        ])
        java_callset = json.loads((java_workspace / "callset.json").read_text(encoding="utf-8"))
        assert java_callset["callsets"][0]["sample_name"] == requested_sample
        java_output = work / "java-genotyped.vcf"
        run_checked([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{java_workspace}",
            "-L", interval, "-O", str(java_output),
            "--create-output-variant-index", "false",
        ])

        native_workspace = work / "native-workspace"
        native_manifest = work / "native.manifest.json"
        native_summary = json.loads(run_checked([
            str(importer), "--sample-name-map", str(sample_map),
            "--genomicsdb-workspace-path", str(native_workspace),
            "--fastgatk-native-workspace", "--output-manifest", str(native_manifest),
        ], env={**os.environ, "FASTGATK_NATIVE_BUILD": str(build)}).stdout.splitlines()[-1])
        assert native_summary["sample_name_map_entries"] == 1
        metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["sample_name_map"] is True
        assert metadata["telemetry"]["sample_name_map_entries"] == 1
        workspace_metadata = json.loads(
            (native_workspace / "fastgatk-workspace.json").read_text(encoding="utf-8"))
        assert workspace_metadata["sample_name_map_entries"] == 1
        assert workspace_metadata["sample_name_map"][0]["sample_name"] == requested_sample
        mapped_input = next((native_workspace / "native-inputs").glob("input-*"))
        assert sample_from_header(mapped_input) == [requested_sample]

        native_output = work / "native-genotyped.vcf"
        run_checked([
            str(genotype), "-R", str(reference), "-V", f"gendb://{native_workspace}",
            "-L", interval, "-O", str(native_output),
            "--create-output-variant-index", "false",
        ])
        assert sample_from_header(java_output) == [requested_sample]
        assert sample_from_header(native_output) == [requested_sample]

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "sample_name_map": True,
        "requested_sample": requested_sample,
        "java_and_native_headers_match": True,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
