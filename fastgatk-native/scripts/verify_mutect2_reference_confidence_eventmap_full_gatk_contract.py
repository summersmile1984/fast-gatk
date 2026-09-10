#!/usr/bin/env python3
"""Pin full EventMap reference-confidence output to GATK 4.6.2.0.

The short BP-resolution contract exercises the symbolic writer shape.  This
fixture covers the full chr17 EventMap window, where the source alternates
between no-variation pileups and PairHMM-qualified, realigned evidence inside
five active transactions.  Both reference-confidence modes must be byte exact
after excluding only execution provenance.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def without_provenance(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith("##GATKCommandLine=")]


def data_rows(lines: list[str]) -> list[str]:
    return [line for line in lines if not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party" / "jdk17" / "bin" / "java"
    gatk = root / "third_party" / "gatk-package" / "gatk-4.6.2.0" / \
        "gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled full EventMap reference-confidence oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "oracle inputs unavailable"}))
        return 0

    rows_by_fixture_mode: dict[str, int] = {}
    # The full window exercises the partition-and-merge path.  This narrow
    # interval deliberately contains one EventMap transaction, so it follows
    # the direct (one active AssemblyRegion) return path instead.  Keep both
    # under the same strict VCF oracle: reference-confidence output cannot
    # silently fall back to the pre-finalize activity pileup when there is
    # only one active region.
    fixtures = (
        ("full", "17:69000-70000", 9, ()),
        # Each output core reuses the source-sized -L assembly input but
        # releases its Host read batch before the next core.  This is the
        # accelerator's bounded-memory execution path; its emitted records
        # must remain identical to GATK's one-pass traversal.
        ("streamed-full", "17:69000-70000", 9, ("--stream-by-region", "400")),
        ("single-eventmap", "17:69460-69490", None, ()),
    )
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-reference-confidence-eventmap-") as directory:
        work = Path(directory)
        for fixture_name, interval, expected_concrete_records, native_extra in fixtures:
            for mode in ("GVCF", "BP_RESOLUTION"):
                gatk_output = work / f"gatk.{fixture_name}.{mode}.vcf.gz"
                native_output = work / f"native.{fixture_name}.{mode}.vcf.gz"
                native_stats = work / f"native.{fixture_name}.{mode}.stats.json"
                common = ["-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
                          "-L", interval, "-ERC", mode,
                          "--create-output-variant-index", "false",
                          "--add-output-vcf-command-line", "false"]
                gatk_run = run([str(java), "-Xmx2g", "-jar", str(gatk), "Mutect2",
                                *common, "-O", str(gatk_output)])
                assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
                native_run = run([str(binary), *common, *native_extra, "-O", str(native_output),
                                  "--stats", str(native_stats)])
                assert native_run.returncode == 0, native_run.stderr[-4000:]
                expected = without_provenance(gatk_output)
                observed = without_provenance(native_output)
                assert observed == expected, (
                    f"{fixture_name} EventMap {mode} VCF differs from GATK")
                rows = data_rows(observed)
                rows_by_fixture_mode[f"{fixture_name}:{mode}"] = len(rows)
                if expected_concrete_records is not None:
                    concrete = [row for row in rows if row.split("\t")[4] != "<NON_REF>"]
                    assert len(concrete) == expected_concrete_records, (mode, concrete)
                if interval == "17:69000-70000" and mode == "BP_RESOLUTION":
                    assert len(rows) == 1001, len(rows)
                if fixture_name == "single-eventmap":
                    stats = json.loads(native_stats.read_text(encoding="utf-8"))
                    assert stats["pairhmm_assembly_region_groups"] == 1, stats
                if fixture_name == "streamed-full":
                    stats = json.loads(native_stats.read_text(encoding="utf-8"))
                    assert stats["stream_by_region"] is True, stats
                    assert stats["streamed_regions"] == 3, stats

    print(json.dumps({
        "status": "pass",
        "oracle": "GATK-4.6.2.0",
        "regions": ["17:69000-70000", "17:69460-69490"],
        "streamed_core_size": 400,
        "modes": ["GVCF", "BP_RESOLUTION"],
        "data_rows": rows_by_fixture_mode,
        "vcf_exact_without_execution_provenance": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
