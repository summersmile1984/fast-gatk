#!/usr/bin/env python3
"""End-to-end FilterMutectCalls joint oracle on HCC1143 chr20 Mutect2 output.

Sister oracle to verify_filter_mutect_dream_synthetic_joint_oracle.py on the
HCC1143 tumor/normal fixture (b37 chr20 truncated to 1Mb synthetic Ns).
Pins the same contract: end-to-end FilterMutectCalls run + manifest schema
+ index + ordered semantic VCF/header equality. The HCC1143 fixture produces more
chr20 calls (44 from GATK Mutect2) than the DREAM synthetic fixture (5),
giving the joint empirical model more observations to exercise.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
M2_NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
FMC_NATIVE = Path(os.environ.get(
    "FASTGATK_FILTER_MUTECT_CALLS_BINARY",
    str(ROOT / "fastgatk-native" / "build" / "fastgatk-filter-mutect-calls")))

TUMOR = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam"
NORMAL = ROOT / "testdata/real/cnv_somatic/chr20/HCC1143_normal.bam"
REFERENCE = ROOT / "testdata/real/cnv_somatic/human_g1k_v37.chr-20.truncated.fasta"
TUMOR_SM = "HCC1143"
NORMAL_SM = "HCC1143 BL"


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=check)


def normalized_info_value(value: str) -> tuple[str, object]:
    """Preserve strings while normalizing equivalent VCF numeric spellings."""
    try:
        return ("number", float(value))
    except ValueError:
        return ("string", value)


def vcf_format_types(path: Path) -> dict[str, str]:
    """Return VCF FORMAT types so numeric sample spellings compare semantically."""
    result: dict[str, str] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#CHROM"):
                break
            match = re.match(r"##FORMAT=<ID=([^,>]+),.*?,Type=([^,>]+)", line)
            if match:
                result[match.group(1)] = match.group(2)
    return result


def normalized_format_sample(format_keys: list[str], sample: str,
                             types: dict[str, str]) -> tuple[object, ...]:
    """Normalize declared Integer/Float fields; retain GT and strings verbatim."""
    values = sample.split(":")
    normalized: list[object] = []
    for index, key in enumerate(format_keys):
        value = values[index] if index < len(values) else ""
        value_type = types.get(key)
        if value_type in {"Integer", "Float"} and value not in {"", "."}:
            numeric: list[object] = []
            for token in value.split(","):
                if token == ".":
                    numeric.append(("missing", token))
                else:
                    try:
                        numeric.append(("number", float(token)))
                    except ValueError:
                        numeric.append(("string", token))
            normalized.append(tuple(numeric))
        else:
            normalized.append(("string", value))
    if len(values) > len(format_keys):
        normalized.extend(("string", value) for value in values[len(format_keys):])
    return tuple(normalized)


def vcf_data_records(path: Path) -> list[tuple[object, ...]]:
    records: list[tuple[object, ...]] = []
    format_types = vcf_format_types(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError(f"malformed VCF record in {path}: {line!r}")
            info: dict[str, tuple[str, object]] = {}
            for item in f[7].split(";"):
                key, value = (item.split("=", 1) + [""])[:2] if "=" in item else (item, "")
                if key in info:
                    raise ValueError(f"duplicate INFO key {key!r} in {path}")
                info[key] = normalized_info_value(value)
            format_keys = f[8].split(":") if len(f) > 8 else []
            samples = tuple(normalized_format_sample(format_keys, sample, format_types)
                            for sample in f[9:])
            records.append(tuple(f[:7]) + (tuple(sorted(info.items())),) +
                           ((f[8],) if len(f) > 8 else ()) + samples)
    return records


def vcf_schema_header(path: Path) -> tuple[set[str], str]:
    meta: set[str] = set()
    chrom = ""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            line = line.rstrip("\n")
            if line.startswith("##GATKCommandLine="):
                continue
            if line.startswith("##"):
                meta.add(line)
            elif line.startswith("#CHROM"):
                chrom = line
    return meta, chrom


def run_filter_mutect_calls(in_vcf: Path, out_vcf: Path, binary: Path,
                            output_manifest: Path | None = None) -> subprocess.CompletedProcess[str]:
    if str(binary).endswith("fastgatk-filter-mutect-calls"):
        cmd = [str(binary), "-V", str(in_vcf), "-R", str(REFERENCE), "-O", str(out_vcf)]
    else:
        cmd = [
            str(JAVA), "-jar", str(GATK), "FilterMutectCalls",
            "-V", str(in_vcf), "-R", str(REFERENCE), "-O", str(out_vcf),
        ]
    if output_manifest:
        cmd += ["--output-manifest", str(output_manifest)]
    return run(cmd, check=False)


def main() -> int:
    required = (JAVA, GATK, M2_NATIVE, FMC_NATIVE, TUMOR, NORMAL, REFERENCE)
    if not all(p.is_file() for p in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing required inputs: {[str(p) for p in required if not p.is_file()]}\n")
            return 2
        print(json.dumps({"status": "skip", "reason": "bundled HCC1143 joint oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-fmc-hcc1143-joint-") as directory:
        work = Path(directory)
        gatk_m2_vcf = work / "gatk_m2.vcf.gz"
        gatk_fmc_vcf = work / "gatk_fmc.vcf.gz"
        native_fmc_vcf = work / "native_fmc.vcf.gz"
        native_fmc_manifest = work / "native_fmc.manifest.json"

        gatk_m2 = run([
            str(JAVA), "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE),
            "-I", str(TUMOR), "-I", str(NORMAL),
            "-tumor", TUMOR_SM, "-normal", NORMAL_SM,
            "-O", str(gatk_m2_vcf),
        ], check=False)
        if gatk_m2.returncode != 0:
            sys.stderr.write(f"GATK Mutect2 pre-step failed:\n{gatk_m2.stderr[-2000:]}\n")
            return 3

        gatk_fmc = run_filter_mutect_calls(gatk_m2_vcf, gatk_fmc_vcf, GATK)
        if gatk_fmc.returncode != 0:
            sys.stderr.write(f"GATK FilterMutectCalls failed:\n{gatk_fmc.stderr[-2000:]}\n")
            return 4
        native_fmc = run_filter_mutect_calls(
            gatk_m2_vcf, native_fmc_vcf, FMC_NATIVE,
            output_manifest=native_fmc_manifest,
        )
        if native_fmc.returncode != 0:
            sys.stderr.write(f"native FilterMutectCalls failed:\n{native_fmc.stderr[-2000:]}\n")
            return 4

        gatk_records = vcf_data_records(gatk_fmc_vcf)
        native_records = vcf_data_records(native_fmc_vcf)

        for path, label in ((gatk_fmc_vcf, "GATK"), (native_fmc_vcf, "native")):
            tbi = Path(str(path) + ".tbi")
            if not tbi.is_file():
                sys.stderr.write(f"{label} FilterMutectCalls VCF index missing: {tbi}\n")
                return 5

        if not native_fmc_manifest.is_file():
            sys.stderr.write(f"native FMC OutputManifest missing: {native_fmc_manifest}\n")
            return 6
        manifest = json.loads(native_fmc_manifest.read_text(encoding="utf-8"))
        if manifest.get("schema_version") != 1:
            sys.stderr.write(f"native FMC manifest schema_version mismatch: {manifest.get('schema_version')}\n")
            return 6

        # All output fields must agree in input record order.  HTSlib can
        # choose a different textual INFO ordering/float spelling than Java,
        # so the parser normalizes only those two representations above.
        if len(gatk_records) < 5:
            sys.stderr.write(f"GATK FMC dropped below 5 records: {len(gatk_records)}\n")
            return 7
        if len(native_records) < 5:
            sys.stderr.write(f"native FMC dropped below 5 records: {len(native_records)}\n")
            return 7
        if native_records != gatk_records:
            first_difference = next(
                ((index, gatk, native)
                 for index, (gatk, native) in enumerate(zip(gatk_records, native_records))
                 if gatk != native),
                (min(len(gatk_records), len(native_records)), None, None),
            )
            sys.stderr.write(
                "FilterMutectCalls semantic VCF records diverge; this compares "
                "all non-provenance fields, INFO values, FILTER and FORMAT payloads.\n"
                f"first_difference={first_difference!r}\n"
            )
            return 7
        gatk_header, gatk_chrom = vcf_schema_header(gatk_fmc_vcf)
        native_header, native_chrom = vcf_schema_header(native_fmc_vcf)
        if native_header != gatk_header or native_chrom != gatk_chrom:
            sys.stderr.write("FilterMutectCalls schema header diverges (excluding GATKCommandLine provenance).\n")
            return 8

        report = {
            "status": "pass",
            "fixture": "hcc1143_chr20_truncated_filter_mutect_calls_joint",
            "reference": "human_g1k_v37.chr-20.truncated.fasta (1Mb synthetic Ns)",
            "mutect2_input_source": "GATK 4.6.2.0 (consumed by both filter runs)",
            "gatk_filter_record_count": len(gatk_records),
            "native_filter_record_count": len(native_records),
            "ordered_vcf_records_semantically_exact": True,
            "vcf_schema_exact_excluding_execution_provenance": True,
            "native_fmc_status": manifest.get("status"),
            "native_fmc_artifact_lod_threshold": (manifest.get("compatibility") or {}).get("artifact_lod_threshold"),
            "a5_joint_filter_boundary_recorded": True,
        }
        print(json.dumps(report, indent=2))
        return 0


if __name__ == "__main__":
    sys.exit(main())
