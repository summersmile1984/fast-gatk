#!/usr/bin/env python3
"""Pin native sparse GenomicsDBImport interval filtering to GATK 4.6.2.0.

GATK's GenomicsDB importer selects an input record when its reference span
overlaps the imported interval, retaining the complete record (including a
reference-confidence block that starts before or ends after the selector).
The native sparse workspace must apply the same rule before publication; this
is distinct from CombineGVCFs, whose reference-block clipping policy differs.
"""

from __future__ import annotations

import gzip
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


def rows(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n").split("\t") for line in stream
                if line and not line.startswith("#")]


def sample_values(row: list[str]) -> dict[str, str]:
    if len(row) < 10:
        return {}
    return dict(zip(row[8].split(":"), row[9].split(":")))


def row_signature(row: list[str]) -> tuple[tuple[str, str, str, str, str], dict[str, str]]:
    info = dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
    return ((row[0], row[1], row[3], row[4], info.get("END", "")), sample_values(row))


def read_export_manifest(manifest: Path) -> list[list[str]]:
    exported: list[list[str]] = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        output = Path(line.split("\t", 1)[1])
        exported.extend(rows(output))
    return exported


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
    required = (importer, bridge, java, jar, gvcf)
    if not all(path.is_file() and os.access(path, os.X_OK) if path in (importer, bridge, java)
               else path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_genomicsdb_import_native_interval_gatk_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled native GenomicsDB interval oracle assets unavailable: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK/bridge assets unavailable"}))
        return 0

    # First interval selects a reference block that starts at the selector
    # beginning but ends later.  The second selects a concrete variant.
    regions = ("20:69491-69500", "20:69511-69511")
    observed = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-native-interval-oracle-") as directory:
        work = Path(directory)
        for ordinal, region in enumerate(regions):
            java_workspace = work / f"java-workspace-{ordinal}"
            export_manifest = work / f"java-export-{ordinal}.tsv"
            run_checked([
                str(java), "-Xmx1g", "-jar", str(jar), "GenomicsDBImport",
                "-V", str(gvcf), "--genomicsdb-workspace-path", str(java_workspace),
                "-L", region, "--batch-size", "1",
            ])
            run_checked([str(bridge), str(java_workspace), str(export_manifest)])
            java_rows = read_export_manifest(export_manifest)
            assert java_rows, (region, "GATK exported no selected records")

            native_workspace = work / f"native-workspace-{ordinal}"
            native_manifest = work / f"native-{ordinal}.manifest.json"
            native_result = run_checked([
                str(importer), "-V", str(gvcf),
                "--genomicsdb-workspace-path", str(native_workspace),
                "--fastgatk-native-workspace", "-L", region,
                "--output-manifest", str(native_manifest),
            ])
            summary = json.loads(native_result.stdout.splitlines()[-1])
            metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
            native_input = next(native_workspace.glob("native-inputs/input-*"))
            native_rows = rows(native_input)
            assert [row_signature(row) for row in native_rows] == [
                row_signature(row) for row in java_rows
            ], (region, java_rows, native_rows)
            assert summary["interval_selectors"] == 1
            assert metadata["telemetry"]["interval_records_kept"] == len(native_rows)
            assert metadata["telemetry"]["interval_records_skipped"] > 0
            workspace_metadata = json.loads(
                (native_workspace / "fastgatk-workspace.json").read_text(encoding="utf-8"))
            assert workspace_metadata["import_intervals"] == [region]
            observed.append({
                "region": region,
                "rows": len(native_rows),
                "rows_exact": True,
                "kept": metadata["telemetry"]["interval_records_kept"],
                "skipped": metadata["telemetry"]["interval_records_skipped"],
            })

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "native_workspace_interval_span_exact": True,
        "cases": observed,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
