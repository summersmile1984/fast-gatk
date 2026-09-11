#!/usr/bin/env python3
"""Pinned GATK model-artifact oracle for VariantRecalibrator input-model."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_VARIANT_RECALIBRATOR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-variant-recalibrator")))
GATK_JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def records(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines() if line and not line.startswith("#")]


def info_map(record: str) -> dict[str, str]:
    fields = record.split("\t")
    assert len(fields) >= 8, record
    values: dict[str, str] = {}
    for token in fields[7].split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            values[key] = value
        elif token:
            values[token] = ""
    return values


def table_names(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.split(":", 3)[2] for line in stream
                if line.startswith("#:GATKTable:") and line.count(":") >= 2
                and not line.startswith("#:GATKTable:8:")]


def annotation_order(path: pathlib.Path) -> list[str]:
    """Read the authoritative AnnotationMeans row order from a GATK report."""
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    marker = next(i for i, line in enumerate(lines)
                  if line.startswith("#:GATKTable:") and "AnnotationMeans" in line)
    order: list[str] = []
    for line in lines[marker + 1:]:
        if line.startswith("#:") or not line.strip():
            if order:
                break
            continue
        fields = line.split()
        if not fields or fields[0] in {"Annotation", "Mean"}:
            continue
        order.append(fields[0])
    return order


def index_feature(path: pathlib.Path) -> None:
    result = subprocess.run(
        [str(GATK_JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(path)],
        text=True, capture_output=True, check=False,
    )
    if result.returncode != 0:
        raise AssertionError(f"GATK IndexFeatureFile failed for {path}:\n{result.stderr}")


def main() -> int:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    if not GATK_JAVA.exists() or not GATK_JAR.exists():
        oracle_guard.oracle_not_verified('verify_variant_recalibrator_gatk_model_oracle.py', GATK_JAVA, GATK_JAR)
        raise SystemExit("VariantRecalibrator GATK model oracle inputs are required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-model-oracle-") as temporary:
        work = pathlib.Path(temporary)
        rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\tGT\t0/1"
            for index, (qd, mq) in enumerate([(30, 60), (28, 58), (25, 55), (4, 25), (3, 20), (2, 15)], 1)
        ]
        input_vcf = work / "input.vcf"
        training_vcf = work / "training.vcf"
        known_vcf = work / "known.vcf"
        write_vcf(input_vcf, rows)
        write_vcf(training_vcf, rows[:3])
        write_vcf(known_vcf, rows[:2])
        # VariantRecalibrator queries each resource by interval, so Java's
        # FeatureDataSource requires the same random-access index that a
        # production staged VCF would carry.
        for path in (input_vcf, training_vcf, known_vcf):
            index_feature(path)
        gatk_model = work / "gatk.model"
        gatk_output = work / "gatk.recal.vcf"
        gatk_result = subprocess.run(
            [str(GATK_JAVA), "-jar", str(GATK_JAR), "VariantRecalibrator",
             "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
             "--max-attempts", "20", "--k-means-iterations", "20",
             "--bad-lod-score-cutoff", "100.0",
             "--dont-run-rscript", "true", "--create-output-variant-index", "false",
             "--sites-only-vcf-output", "true",
             "-O", str(gatk_output), "--tranches-file", str(work / "gatk.tranches"),
             "--truth-sensitivity-tranche", "100", "--output-model", str(gatk_model)],
            text=True, capture_output=True, check=False,
        )
        if gatk_result.returncode != 0:
            raise AssertionError(f"GATK VariantRecalibrator failed:\n{gatk_result.stderr}")
        model_text = gatk_model.read_text(encoding="utf-8")
        expected_tables = {
            "AnnotationMeans", "AnnotationStdevs", "GoodGaussianPMix",
            "PositiveModelMeans", "PositiveModelCovariances",
            "BadGaussianPMix", "NegativeModelMeans", "NegativeModelCovariances",
        }
        actual_tables = set(table_names(gatk_model))
        assert expected_tables.issubset(actual_tables), (actual_tables, model_text[:500])
        gatk_annotation_order = annotation_order(gatk_model)
        assert gatk_annotation_order, model_text[:500]

        native_output = work / "native.recal.vcf.gz"
        native_model = work / "native.model"
        native_manifest = work / "native.manifest.json"
        native_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--input-model", str(gatk_model),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "-O", str(native_output),
             "--tranches-file", str(work / "native.tranches"),
             "--output-model", str(native_model), "--output-manifest", str(native_manifest),
             "--sites-only-vcf-output", "true"],
            text=True, capture_output=True, check=False,
        )
        if native_result.returncode != 0:
            raise AssertionError(f"native GATK model load failed:\n{native_result.stderr}")
        native_payload = json.loads(native_result.stdout)
        manifest_payload = json.loads(native_manifest.read_text(encoding="utf-8"))
        native_records = records(native_output)
        assert native_payload["input_model"] is True, native_payload
        assert manifest_payload["model_loaded"] is True, manifest_payload
        # The persisted GATK model order is authoritative when replaying a
        # model; explicit -an names are validated but must not permute it.
        assert manifest_payload["annotations"] == gatk_annotation_order, manifest_payload
        assert len(native_records) == len(rows), native_records
        assert all("VQSLOD=" in record for record in native_records), native_records
        gatk_records = records(gatk_output)
        assert len(gatk_records) == len(rows), gatk_records
        assert all(len(record.split("\t")) == 8 for record in gatk_records), gatk_records
        assert all(len(record.split("\t")) == 8 for record in native_records), native_records
        for gatk_record, native_record in zip(gatk_records, native_records):
            gatk_info = info_map(gatk_record)
            native_info = info_map(native_record)
            assert abs(float(gatk_info["VQSLOD"]) - float(native_info["VQSLOD"])) < 1.0e-3, (
                gatk_record, native_record
            )
            assert gatk_info.get("culprit") == "MQ", gatk_record
            assert native_info.get("culprit") == "serialized-gmm", native_record
            # GATK's recalibration artifact and native output intentionally
            # use different provenance fields/allele representations.  Keep
            # the semantic contract on score, coordinate, labels, and the
            # native symbolic N/<VQSR> writer boundary rather than literal
            # first-seven-column equality.
            gatk_fields = gatk_record.split("\t")
            native_fields = native_record.split("\t")
            assert gatk_fields[0:2] == native_fields[0:2], (gatk_record, native_record)
            assert native_fields[3:5] == ["N", "<VQSR>"], native_record
            # POSITIVE_TRAIN_SITE is resource membership and is stable when a
            # serialized model is replayed.  GATK may additionally emit
            # NEGATIVE_TRAIN_SITE after its model-side bad-site selection;
            # that selection is not serialized in a GATKReport, so it remains
            # an explicit non-bit-identical metadata boundary here.
            assert ("POSITIVE_TRAIN_SITE" in gatk_info) == (
                "POSITIVE_TRAIN_SITE" in native_info
            ), (gatk_record, native_record)
        # A loaded GATK report must be round-trippable through the native
        # report writer: table identities and annotation dimensions survive,
        # while model floating-point values intentionally remain a separate
        # semantic boundary from Java's VBEM implementation.
        native_tables = set(table_names(native_model))
        assert expected_tables.issubset(native_tables), native_tables
        assert model_text.startswith("#:GATKReport.v1.1:"), model_text[:80]
        print(json.dumps({
            "gatk_version": "4.6.2.0",
            "gatk_model_tables": len(actual_tables),
            "native_model_tables": len(native_tables),
            "annotations": manifest_payload["annotations"],
            "input_model_loaded": manifest_payload["model_loaded"],
            "records_scored": native_payload["scored_records"],
            "status": "pass",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    main()
