#!/usr/bin/env python3
"""Pin Mutect2's repeatable ReadThreadingAssembler k-mer contract to GATK.

GATK 4.6.2.0 defaults to ``--kmer-size 10 --kmer-size 25`` and normalizes
explicit repeated values to the same sorted retry list.  A single native k=10
attempt can silently diverge on cyclic active regions, so this checks both the
CLI/list handoff and a real tumor/normal call window.  It also pins a k=35
request: GATK's retry series extends beyond a 64-bit two-bit k-mer, while the
native Host front end must retain the full string key before Kokkos receives
its compact node IDs.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", ROOT / "fastgatk-native/build/fastgatk-mutect2"))
REFERENCE = ROOT / "testdata/downloads/reference/hs37d5.fa.gz"
TUMOR = ROOT / "testdata/real/dream_synthetic/chr20/tumor.bam"
NORMAL = ROOT / "testdata/real/dream_synthetic/chr20/normal.bam"
TUMOR_SM = "synthetic.challenge.set1.tumor"
NORMAL_SM = "synthetic.challenge.set1.normal"
WINDOW = "20:10022000-10023500"


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout[-2000:]}\nstderr={result.stderr[-4000:]}")


def records(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream if not line.startswith("#")]


def allele_identities(rows: list[str]) -> list[tuple[str, str, str, str]]:
    return [tuple(row.split("\t")[field] for field in (0, 1, 3, 4)) for row in rows]


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL)
    if not all(path.is_file() for path in required):
        raise RuntimeError(f"missing oracle inputs: {[str(path) for path in required if not path.is_file()]}")

    common = [
        "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
        "--tumor-sample", TUMOR_SM, "--normal-sample", NORMAL_SM,
        "-L", WINDOW,
    ]
    repeated = ["--kmer-size", "25", "--kmer-size", "10"]
    long_kmer = ["--kmer-size", "35"]
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-kmer-list-") as directory:
        work = Path(directory)
        gatk_default = work / "gatk-default.vcf.gz"
        gatk_repeated = work / "gatk-repeated.vcf.gz"
        native_default = work / "native-default.vcf.gz"
        native_repeated = work / "native-repeated.vcf.gz"
        gatk_long_kmer = work / "gatk-long-kmer.vcf.gz"
        native_long_kmer = work / "native-long-kmer.vcf.gz"
        native_default_manifest = work / "native-default.json"
        native_repeated_manifest = work / "native-repeated.json"
        native_long_kmer_manifest = work / "native-long-kmer.json"
        run([
            str(JAVA), "-jar", str(GATK), "Mutect2", *common,
            "-O", str(gatk_default),
        ])
        run([
            str(JAVA), "-jar", str(GATK), "Mutect2", *common, *repeated,
            "-O", str(gatk_repeated),
        ])
        run([
            str(JAVA), "-jar", str(GATK), "Mutect2", *common, *long_kmer,
            "-O", str(gatk_long_kmer),
        ])
        run([
            str(NATIVE), *common, "-O", str(native_default),
            "--output-manifest", str(native_default_manifest),
        ])
        run([
            str(NATIVE), *common, *repeated, "-O", str(native_repeated),
            "--output-manifest", str(native_repeated_manifest),
        ])
        run([
            str(NATIVE), *common, *long_kmer, "-O", str(native_long_kmer),
            "--output-manifest", str(native_long_kmer_manifest),
        ])

        gatk_default_records = records(gatk_default)
        gatk_repeated_records = records(gatk_repeated)
        native_default_records = records(native_default)
        native_repeated_records = records(native_repeated)
        gatk_long_kmer_records = records(gatk_long_kmer)
        native_long_kmer_records = records(native_long_kmer)
        assert gatk_default_records == gatk_repeated_records, (
            "GATK default/repeated k-mer records differ",
            gatk_default_records, gatk_repeated_records)
        assert native_default_records == gatk_default_records, (
            "native default k-mer records differ from GATK",
            native_default_records, gatk_default_records)
        assert native_repeated_records == gatk_repeated_records, (
            "native repeated k-mer records differ from GATK",
            native_repeated_records, gatk_repeated_records)
        assert native_default_records == native_repeated_records, (
            "native default/repeated k-mer records differ",
            native_default_records, native_repeated_records)
        assert allele_identities(native_long_kmer_records) == allele_identities(gatk_long_kmer_records), (
            "native k=35 allele identities differ from GATK",
            native_long_kmer_records, gatk_long_kmer_records)

        default_manifest = json.loads(native_default_manifest.read_text(encoding="utf-8"))
        repeated_manifest = json.loads(native_repeated_manifest.read_text(encoding="utf-8"))
        long_kmer_manifest = json.loads(native_long_kmer_manifest.read_text(encoding="utf-8"))
        default_telemetry = default_manifest["telemetry"]
        repeated_telemetry = repeated_manifest["telemetry"]
        assert default_telemetry["graph_kmer_sizes"] == [10, 25], default_telemetry
        assert repeated_telemetry["graph_kmer_sizes"] == [10, 25], repeated_telemetry
        assert default_telemetry["graph_kmer_size"] == 10, default_telemetry
        assert repeated_telemetry["graph_kmer_size"] == 10, repeated_telemetry
        long_kmer_telemetry = long_kmer_manifest["telemetry"]
        assert long_kmer_telemetry["graph_kmer_sizes"] == [35], long_kmer_telemetry
        assert long_kmer_telemetry["graph_kmer_size"] == 35, long_kmer_telemetry

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": WINDOW,
            "default_kmer_sizes": [10, 25],
            "reordered_explicit_kmer_sizes": [25, 10],
            "normalized_kmer_sizes": [10, 25],
            "gatk_default_repeated_records_exact": True,
            "native_default_gatk_records_exact": True,
            "native_repeated_gatk_records_exact": True,
            "long_kmer": 35,
            "native_long_kmer_gatk_alleles_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
