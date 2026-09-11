#!/usr/bin/env python3
"""Pinned GATK oracle for the Picard-style GatherVcfs CLI boundary."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def data_records(path: Path) -> list[str]:
    opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
    with opener as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-gather-vcfs"
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not all(path.is_file() for path in (binary, java, gatk)):
        oracle_guard.oracle_not_verified('verify_gather_vcfs_cli_boundary_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("GatherVcfs CLI oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "native or pinned GATK assets unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-cli-oracle-") as directory:
        work = Path(directory)
        shard1 = work / "shard-1.vcf"
        shard2 = work / "shard-2.vcf"
        shard1.write_text(HEADER + "chr1\t1\tfirst\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n", encoding="utf-8")
        shard2.write_text(HEADER + "chr1\t2\tsecond\tA\tT\t60\tPASS\tDP=12\tGT\t0/0\n", encoding="utf-8")

        native_output = work / "native.vcf.gz"
        native_manifest = work / "native.manifest.json"
        native = run([
            str(binary), "--INPUT", str(shard1), "--INPUT", str(shard2),
            "--OUTPUT", str(native_output), "--COMMENT", "cli-boundary",
            "--CREATE_INDEX", "false", "--COMPRESSION_LEVEL", "0",
            "--QUIET", "true", "--VERBOSITY", "WARNING",
            "--output-manifest", str(native_manifest)])
        assert native.returncode == 0, native.stderr
        assert native_output.is_file() and not Path(f"{native_output}.tbi").exists()
        native_text = gzip.open(native_output, "rt", encoding="utf-8").read()
        assert "##GatherVcfs.comment=cli-boundary" in native_text

        gatk_output = work / "gatk.vcf.gz"
        gatk_result = run([
            str(java), "-Xmx1g", "-jar", str(gatk), "GatherVcfs",
            "--INPUT", str(shard1), "--INPUT", str(shard2),
            "--OUTPUT", str(gatk_output), "--COMMENT", "cli-boundary",
            "--CREATE_INDEX", "false", "--COMPRESSION_LEVEL", "0",
            "--QUIET", "true", "--VERBOSITY", "WARNING",
        ])
        assert gatk_result.returncode == 0, gatk_result.stderr
        assert not Path(f"{gatk_output}.tbi").exists()
        assert data_records(native_output) == data_records(gatk_output)
        metadata = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert metadata["compatibility"]["comments"] == 1
        assert metadata["compatibility"]["compression_level"] == 0
        assert metadata["compatibility"]["vcf_index"] is False
        assert metadata["compatibility"]["quiet"] is True

        # The short Picard aliases are independently covered by the existing
        # gather contract; this oracle also checks that the mixed uppercase
        # spelling can be used in a direct GATK replacement command.
        summary = json.loads(native.stdout.splitlines()[-1])
        assert summary["input_files"] == 2 and summary["input_records"] == 2
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(data_records(native_output)),
            "uppercase_input_output": True,
            "comment_exact": True,
            "compression_level": 0,
            "create_index_false": True,
            "quiet_verbosity_accepted": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
