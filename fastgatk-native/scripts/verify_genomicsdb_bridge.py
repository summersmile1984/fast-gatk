#!/usr/bin/env python3
"""Exercise direct replacement of a real GATK TileDB/GenomicsDB workspace.

The native sparse workspace tests intentionally do not claim TileDB
compatibility. This verifier creates a genuine workspace with the pinned GATK
jar, exports it through the isolated legacy-ABI bridge, and then consumes it
through native GenotypeGVCFs. GT/AD/PL are compared with GATK's direct gendb
query; annotation formatting remains covered by the broader GenotypeGVCFs
oracle.
"""

from __future__ import annotations

import json
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def records(path: Path) -> dict[tuple[str, int, str, str], list[str]]:
    result: dict[tuple[str, int, str, str], list[str]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            result[(fields[0], int(fields[1]), fields[3], fields[4])] = fields
    return result


def sample_fields(row: list[str]) -> dict[str, str]:
    return dict(zip(row[8].split(":"), row[9].split(":")))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    helper = Path(os.environ.get(
        "FASTGATK_GENOMICSDB_BRIDGE_BINARY", build / "fastgatk-genomicsdb-export"))
    genotype = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY", build / "fastgatk-genotype-gvcf"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gvcf = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/engine/GenomicsDBIntegration/tiny.g.vcf"
    reference = root / "reference.fa"
    required = (helper, genotype, jar, java, gvcf, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GENOMICSDB_BRIDGE"):
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"missing GenomicsDB bridge assets: {missing}")
        print(json.dumps({"status": "skipped", "reason": "GenomicsDB bridge assets unavailable"}))
        return 0

    region = "20:69491-70000"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-bridge-oracle-") as directory:
        work = Path(directory)
        workspace = work / "gatk-workspace"
        exported_manifest = work / "export.tsv"
        gatk_import = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "GenomicsDBImport",
            "-V", str(gvcf), "--genomicsdb-workspace-path", str(workspace),
            "--intervals", region, "--batch-size", "1",
        ], text=True, capture_output=True, check=False,
            env={**os.environ, "PATH": f"{java.parent}:{os.environ.get('PATH', '')}"})
        assert gatk_import.returncode == 0, gatk_import.stderr[-4000:]
        assert (workspace / "callset.json").is_file()
        assert (workspace / "vidmap.json").is_file()
        assert any(path.is_dir() and (path / "__array_schema.tdb").is_file()
                   for path in workspace.iterdir())

        exported = subprocess.run(
            [str(helper), str(workspace), str(exported_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert exported.returncode == 0, exported.stderr[-4000:]
        output_paths = []
        for line in exported_manifest.read_text(encoding="utf-8").splitlines():
            if line and not line.startswith("#"):
                output_paths.append(Path(line.split("\t", 1)[1]))
        assert output_paths and all(path.is_file() and path.stat().st_size > 0
                                    for path in output_paths)
        exported_records = {}
        for path in output_paths:
            exported_records.update(records(path))
        assert exported_records, "bridge exported no records"
        exported_sites = {
            (contig, pos, ref, alt)
            for (contig, pos, ref, alt) in exported_records
            for alt in alt.split(",")
            if alt != "<NON_REF>"
        }

        # Exercise TileDB column pushdown at a single-base boundary.  The
        # pinned tiny fixture has adjacent records at 20:69511 and 20:69512;
        # an inclusive-end conversion bug would therefore leak the latter.
        narrow_manifest = work / "export-narrow.tsv"
        narrow = subprocess.run(
            [str(helper), str(workspace), str(narrow_manifest),
             "--interval", "20:69511-69511"],
            text=True, capture_output=True, check=False,
        )
        assert narrow.returncode == 0, narrow.stderr[-4000:]
        narrow_paths = [
            Path(line.split("\t", 1)[1])
            for line in narrow_manifest.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        assert narrow_paths and all(path.is_file() and path.stat().st_size > 0
                                    for path in narrow_paths)
        narrow_records = {}
        for path in narrow_paths:
            narrow_records.update(records(path))
        assert narrow_records, "narrow bridge export no records"
        assert all(contig == "20" and pos == 69511
                   for contig, pos, _ref, _alt in narrow_records), narrow_records

        java_output = work / "gatk-joint.vcf"
        java_joint = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", f"gendb://{workspace}",
            "-L", region, "-O", str(java_output),
            "--create-output-variant-index", "false",
        ], text=True, capture_output=True, check=False)
        assert java_joint.returncode == 0, java_joint.stderr[-4000:]

        native_output = work / "native-joint.vcf"
        native_manifest = work / "native-joint.manifest.json"
        native_joint = subprocess.run([
            str(genotype), "-R", str(reference), "-V", f"gendb://{workspace}",
            "-L", region, "-O", str(native_output),
            "--create-output-variant-index", "false",
            "--output-manifest", str(native_manifest),
        ], text=True, capture_output=True, check=False,
            env={**os.environ, "FASTGATK_GENOMICSDB_BRIDGE": str(helper)})
        assert native_joint.returncode == 0, native_joint.stderr[-4000:]
        native_summary = json.loads(native_joint.stdout.splitlines()[-1])
        metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert native_summary["genomicsdb_bridge"] is True
        assert metadata["compatibility"]["genomicsdb_bridge"] is True
        assert metadata["telemetry"]["genomicsdb_bridge"] is True
        assert metadata["telemetry"]["indexed_inputs"] >= 1
        assert metadata["telemetry"]["indexed_interval_queries"] >= 1
        gatk_records = records(java_output)
        native_records = records(native_output)
        assert set(gatk_records) == set(native_records), (set(gatk_records), set(native_records))
        assert set(native_records).issubset(exported_sites)
        for key in gatk_records:
            gatk_sample = sample_fields(gatk_records[key])
            native_sample = sample_fields(native_records[key])
            for field in ("GT", "AD", "PL"):
                assert native_sample.get(field) == gatk_sample.get(field), (key, field)

    print(json.dumps({"status": "pass", "arrays": len(output_paths),
                      "exported_records": len(exported_records),
                      "genotyped_records": len(native_records)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
