#!/usr/bin/env python3
"""Run GenotypeGVCFs on the same GATK gVCF and compare native joint output.

This is deliberately a diagnostic oracle: it gates the shared record/FORMAT
contract and reports exact GT/GQ/PL/QUAL/AF agreement without pretending that
the native implementation is already bit-identical for every cohort model.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def fields(path: Path) -> list[list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def semantic_header(path: Path) -> list[str]:
    """Return stable header definitions, excluding run/provenance metadata."""
    opener = gzip.open if path.suffix == ".gz" else open
    ignored = ("##GATKCommandLine=", "##source=", "##fileDate=",
               "##contig=", "##GVCFBlock", "##fastgatk_genotype_gvcfs_status=")
    return sorted({line.rstrip() for line in opener(path, "rt", encoding="utf-8")
                   if line.startswith("##") and not line.startswith(ignored)})


def record_map(path: Path) -> dict[tuple[str, int, str, str], list[str]]:
    result: dict[tuple[str, int, str, str], list[str]] = {}
    for value in fields(path):
        result[(value[0], int(value[1]), value[3], value[4])] = value
    return result


def sample_format(record: list[str]) -> dict[str, str]:
    names = record[8].split(":")
    values = record[9].split(":") if len(record) > 9 else []
    return dict(zip(names, values))


def info_map(record: list[str]) -> dict[str, str]:
    return {item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in record[7].split(";") if "=" in item}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or GATK fixtures; run build_native.sh first")
    if not gatk_jar.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-oracle-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf.gz"
        gatk_joint = work / "gatk.joint.vcf.gz"
        native_joint = work / "native.joint.vcf.gz"
        native_manifest = work / "native.joint.manifest.json"
        gatk_nonvariant = work / "gatk.nonvariant.vcf"
        native_nonvariant = work / "native.nonvariant.vcf.gz"
        native_nonvariant_manifest = work / "native.nonvariant.manifest.json"
        gatk_gp_qual = work / "gatk.gp-qual.vcf.gz"
        native_gp_qual = work / "native.gp-qual.vcf.gz"
        native_gp_qual_manifest = work / "native.gp-qual.manifest.json"
        # Use the complete bundled chr17 window.  It contains both a simple
        # SNP and records with rank-sum annotations, so strict text parity is
        # exercised across multiple records rather than a single coincidence.
        region = "17:69000-70000"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-O", str(gatk_gvcf), "-ERC", "GVCF",
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "true",
            "--seconds-between-progress-updates", "1",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        gatk_genotype = subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-L", region,
            "-O", str(gatk_joint), "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True)
        if gatk_genotype.returncode != 0:
            raise RuntimeError(f"GATK GenotypeGVCFs failed: {gatk_genotype.stderr[-4000:]}")
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(gatk_gvcf),
            "-L", region, "-O", str(native_joint),
            "--gatk-compatible-annotations",
            "--output-manifest", str(native_manifest),
        ], check=True, text=True, capture_output=True)
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        gatk_header = semantic_header(gatk_joint)
        native_header = semantic_header(native_joint)
        gatk_records = record_map(gatk_joint)
        native_records = record_map(native_joint)
        gatk_nonvariant_result = subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-L", region,
            "-O", str(gatk_nonvariant), "--include-non-variant-sites",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True)
        if gatk_nonvariant_result.returncode != 0:
            raise RuntimeError(f"GATK include-non-variant-sites failed: {gatk_nonvariant_result.stderr[-4000:]}")
        native_nonvariant_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(gatk_gvcf),
            "-L", region, "-O", str(native_nonvariant),
            "--include-non-variant-sites", "--output-manifest", str(native_nonvariant_manifest),
        ], check=True, text=True, capture_output=True)
        native_nonvariant_manifest_data = json.loads(
            native_nonvariant_manifest.read_text(encoding="utf-8"))
        gatk_nonvariant_records = record_map(gatk_nonvariant)
        native_nonvariant_records = record_map(native_nonvariant)
        # GATK's --use-posteriors-to-calculate-qual is conditional: it only
        # rewrites QUAL when an input genotype carries a GP posterior field.
        # The HC gVCF generated above has PL/GQ but no GP, so this invocation
        # is a useful regression guard against accidentally treating PL as a
        # posterior and changing QUAL.
        gatk_gp_result = subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk_jar), "GenotypeGVCFs",
            "-R", str(reference), "-V", str(gatk_gvcf), "-L", region,
            "-O", str(gatk_gp_qual), "--use-posteriors-to-calculate-qual",
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True)
        if gatk_gp_result.returncode != 0:
            raise RuntimeError(f"GATK GenotypeGVCFs --gp-qual failed: {gatk_gp_result.stderr[-4000:]}")
        native_gp_result = subprocess.run([
            str(native), "-R", str(reference), "-V", str(gatk_gvcf),
            "-L", region, "-O", str(native_gp_qual), "--gp-qual",
            "--gatk-compatible-annotations", "--output-manifest",
            str(native_gp_qual_manifest),
        ], check=True, text=True, capture_output=True)
        native_gp_manifest = json.loads(
            native_gp_qual_manifest.read_text(encoding="utf-8"))
        gatk_gp_records = record_map(gatk_gp_qual)
        native_gp_records = record_map(native_gp_qual)

    intersection = sorted(set(gatk_records) & set(native_records))
    assert intersection, "GATK and native GenotypeGVCFs produced no common records"
    common_gt = common_gq = common_pl = common_qual = common_af = 0
    common_record_text = common_info_fields = common_format_fields = 0
    qual_deltas: list[float] = []
    af_pairs: list[tuple[str, str]] = []
    af_numeric_deltas: list[float] = []
    mle_ac_exact = mle_af_exact = 0
    info_field_diffs: list[dict[str, object]] = []
    format_field_diffs: list[dict[str, object]] = []
    for key in intersection:
        gatk = gatk_records[key]
        native_record = native_records[key]
        gatk_sample = sample_format(gatk)
        native_sample = sample_format(native_record)
        if gatk_sample.get("GT") == native_sample.get("GT"):
            common_gt += 1
        if gatk_sample.get("GQ") == native_sample.get("GQ"):
            common_gq += 1
        if gatk_sample.get("PL") == native_sample.get("PL"):
            common_pl += 1
        if gatk[5] == native_record[5]:
            common_qual += 1
        else:
            try:
                qual_deltas.append(abs(float(gatk[5]) - float(native_record[5])))
            except ValueError:
                pass
        gatk_af = info_map(gatk).get("AF", "")
        native_af = info_map(native_record).get("AF", "")
        if gatk_af == native_af:
            common_af += 1
        else:
            af_pairs.append((gatk_af, native_af))
        try:
            af_numeric_deltas.append(abs(float(gatk_af) - float(native_af)))
        except ValueError:
            pass
        gatk_info = info_map(gatk)
        native_info = info_map(native_record)
        if gatk == native_record:
            common_record_text += 1
        gatk_info_keys = sorted(gatk_info)
        native_info_keys = sorted(native_info)
        if gatk_info_keys == native_info_keys:
            common_info_fields += 1
        elif len(info_field_diffs) < 8:
            info_field_diffs.append({
                "record": f"{key[0]}:{key[1]}:{key[2]}>{key[3]}",
                "gatk_only": sorted(set(gatk_info_keys) - set(native_info_keys)),
                "native_only": sorted(set(native_info_keys) - set(gatk_info_keys)),
            })
        gatk_format_keys = sorted(gatk_sample)
        native_format_keys = sorted(native_sample)
        if gatk_format_keys == native_format_keys:
            common_format_fields += 1
        elif len(format_field_diffs) < 8:
            format_field_diffs.append({
                "record": f"{key[0]}:{key[1]}:{key[2]}>{key[3]}",
                "gatk_only": sorted(set(gatk_format_keys) - set(native_format_keys)),
                "native_only": sorted(set(native_format_keys) - set(gatk_format_keys)),
            })
        if gatk_info.get("MLEAC") == native_info.get("MLEAC"):
            mle_ac_exact += 1
        try:
            if abs(float(gatk_info.get("MLEAF", "nan")) -
                   float(native_info.get("MLEAF", "nan"))) <= 1.0e-6:
                mle_af_exact += 1
        except ValueError:
            pass

    assert native_summary["status"] == "contract-compatible"
    assert manifest["compatibility"]["allele_frequency_calculator_em"] is True
    assert manifest["telemetry"]["cohort_af_kernel_calls"] > 0
    assert manifest["telemetry"]["cohort_af_converged"] > 0
    assert manifest["compatibility"]["gatk_annotation_compatibility"] is True
    header_gatk_only = sorted(set(gatk_header) - set(native_header))
    header_native_only = sorted(set(native_header) - set(gatk_header))
    assert not header_gatk_only and not header_native_only
    assert max(qual_deltas, default=0.0) <= 0.01
    assert max(af_numeric_deltas, default=0.0) <= 1.0e-6
    assert mle_ac_exact == len(intersection)
    assert mle_af_exact == len(intersection)
    assert common_info_fields == len(intersection)
    assert common_format_fields == len(intersection)
    # The strict annotation profile now canonicalizes INFO order and the
    # GATK standard float precisions, making the complete data row raw-text
    # identical for the bundled chr17 multi-record oracle window.
    assert common_record_text == len(intersection)
    nonvariant_keys = {
        key for key in gatk_nonvariant_records
        if key[3] == "."
    }
    native_nonvariant_keys = {
        key for key in native_nonvariant_records
        if key[3] == "."
    }
    assert nonvariant_keys == native_nonvariant_keys
    nonvariant_gt_dp_gq_exact = 0
    nonvariant_mismatches: list[dict[str, object]] = []
    for key in sorted(nonvariant_keys):
        gatk = sample_format(gatk_nonvariant_records[key])
        native_record = sample_format(native_nonvariant_records[key])
        gatk_values = {
            "GT": gatk.get("GT"),
            "DP": gatk.get("DP"),
            "RGQ": gatk.get("RGQ"),
            "INFO_DP": info_map(gatk_nonvariant_records[key]).get("DP"),
        }
        native_values = {
            "GT": native_record.get("GT"),
            "DP": native_record.get("DP"),
            "RGQ": native_record.get("GQ"),
            "INFO_DP": info_map(native_nonvariant_records[key]).get("DP"),
        }
        if gatk_values == native_values:
            nonvariant_gt_dp_gq_exact += 1
        elif len(nonvariant_mismatches) < 16:
            nonvariant_mismatches.append({
                "record": f"{key[0]}:{key[1]}:{key[2]}>{key[3]}",
                "gatk": gatk_values,
                "native": native_values,
            })
    assert native_nonvariant_manifest_data["compatibility"]["include_non_variant_sites"] is True
    gp_intersection = sorted(set(gatk_gp_records) & set(native_gp_records))
    assert gp_intersection == sorted(set(gatk_records) & set(native_records))
    gp_qual_exact = sum(
        gatk_gp_records[key][5] == native_gp_records[key][5]
        for key in gp_intersection
    )
    gp_record_text_exact = sum(
        gatk_gp_records[key] == native_gp_records[key]
        for key in gp_intersection
    )
    assert gp_qual_exact == len(gp_intersection)
    assert gp_record_text_exact == len(gp_intersection)
    assert native_gp_manifest["telemetry"]["use_posteriors_to_calculate_qual"] is True
    assert native_gp_manifest["telemetry"]["posterior_kernel_calls"] == 0
    assert native_gp_manifest["telemetry"]["posterior_samples"] == 0
    print(json.dumps({
        "status": "pass",
        "gatk_records": len(gatk_records),
        "native_records": len(native_records),
        "intersection": len(intersection),
        "gt_exact": common_gt,
        "gq_exact": common_gq,
        "pl_exact": common_pl,
        "qual_exact": common_qual,
        "af_exact": common_af,
        "record_text_exact": common_record_text,
        "info_field_set_exact": common_info_fields,
        "format_field_set_exact": common_format_fields,
        "header_semantic_exact": not header_gatk_only and not header_native_only,
        "header_gatk_only": header_gatk_only,
        "header_native_only": header_native_only,
        "max_qual_delta": max(qual_deltas) if qual_deltas else 0.0,
        "af_pairs": af_pairs,
        "info_field_diffs": info_field_diffs,
        "format_field_diffs": format_field_diffs,
        "max_af_numeric_delta": max(af_numeric_deltas) if af_numeric_deltas else 0.0,
        "mle_ac_exact": mle_ac_exact,
        "mle_af_exact": mle_af_exact,
        "cohort_iterations": manifest["telemetry"]["cohort_af_iterations"],
        "include_nonvariant_total_gatk_records": len(gatk_nonvariant_records),
        "include_nonvariant_total_native_records": len(native_nonvariant_records),
        "nonvariant_gatk_records": len(nonvariant_keys),
        "nonvariant_native_records": len(native_nonvariant_keys),
        "nonvariant_gt_dp_gq_exact": nonvariant_gt_dp_gq_exact,
        "nonvariant_mismatches": nonvariant_mismatches,
        "gp_qual_without_posteriors_exact": gp_record_text_exact == len(gp_intersection),
        "gp_qual_records": len(gp_intersection),
        "gp_qual_exact": gp_qual_exact,
        "gp_qual_posterior_kernel_calls": native_gp_manifest["telemetry"]["posterior_kernel_calls"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
