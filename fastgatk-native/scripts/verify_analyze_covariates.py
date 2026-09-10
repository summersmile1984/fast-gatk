#!/usr/bin/env python3
"""Verify AnalyzeCovariates CSV/PDF compatibility against GATK 4.6.2.0."""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, **kwargs)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_ANALYZE_COVARIATES_BINARY",
        root / "fastgatk-native/build/fastgatk-analyze-covariates"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    known_sites = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/BQSR/"
        "dbsnp_132.b37.excluding_sites_after_129.chr17_69k_70k.vcf")
    if not native.exists() or not bam.exists() or not reference.exists():
        print(json.dumps({"status": "skipped", "reason": "missing native fixture"}))
        return 0
    if not jar.exists() or not java.exists() or not known_sites.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit("missing GATK/JDK/known-sites oracle fixture")
        print(json.dumps({"status": "skipped", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-analyze-covariates-") as directory:
        work = Path(directory)
        report = work / "before.table"
        gatk_csv = work / "gatk.csv"
        native_csv = work / "native.csv"
        native_csv_pair = work / "native-pair.csv"
        native_pdf = work / "native.pdf"
        manifest = work / "native.manifest.json"
        run([str(java), "-Xmx1g", "-jar", str(jar), "BaseRecalibrator",
             "-R", str(reference), "-I", str(bam), "-O", str(report),
             "-L", "17:69000-69100", "--known-sites", str(known_sites)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        run([str(java), "-Xmx1g", "-jar", str(jar), "AnalyzeCovariates",
             "-before", str(report), "-csv", str(gatk_csv)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        native_result = run([str(native), "-before", str(report), "-csv", str(native_csv),
                             "-plots", str(native_pdf), "--output-manifest", str(manifest)],
                            capture_output=True)
        assert native_result.stdout.strip()
        assert native_csv.read_bytes() == gatk_csv.read_bytes(), "CSV is not bit-identical"
        assert native_pdf.read_bytes().startswith(b"%PDF-1.4")
        # A repeated before/after report exercises role ordering and the same
        # delta aggregation that GATK uses when two reports are supplied.
        run([str(native), "-before", str(report), "-after", str(report),
             "-csv", str(native_csv_pair)], stdout=subprocess.DEVNULL)
        pair_lines = native_csv_pair.read_text(encoding="utf-8").splitlines()
        assert len(pair_lines) == 2 * len(gatk_csv.read_text(encoding="utf-8").splitlines()) - 1
        mismatched_report = work / "mismatched.table"
        mismatched_report.write_text(
            re.sub(r"(quantizing_levels\s+)16", r"\g<1>17",
                   report.read_text(encoding="utf-8"), count=1),
            encoding="utf-8")
        mismatch_result = subprocess.run(
            [str(native), "-before", str(report), "-after", str(mismatched_report),
             "-csv", str(work / "mismatch.csv")],
            text=True, capture_output=True)
        assert mismatch_result.returncode != 0
        assert "incompatible AnalyzeCovariates report arguments" in mismatch_result.stderr

        # GATK's RecalibrationArgumentCollection.compareReportArguments()
        # intentionally ignores indels_context_size (AnalyzeCovariates only
        # combines already materialized report rows).  Keep native behavior
        # aligned with that bounded Java compatibility rule: a pair that Java
        # accepts must not be rejected before CSV generation.
        ignored_argument_report = work / "ignored-indels-context.table"
        report_lines = report.read_text(encoding="utf-8").splitlines()
        in_arguments = False
        changed = False
        for line_index, line in enumerate(report_lines):
            if line.startswith("#:GATKTable:"):
                in_arguments = line.startswith("#:GATKTable:Arguments:")
            if in_arguments and line.startswith("indels_context_size"):
                match = re.search(r"(indels_context_size\s+)(3)(\s*)$", line)
                assert match, "unexpected fixed-width indels_context_size row"
                report_lines[line_index] = (
                    line[:match.start(2)] + "4" + line[match.end(2):])
                changed = True
                break
        assert changed
        ignored_argument_report.write_text(
            "\n".join(report_lines) + "\n", encoding="utf-8")
        gatk_pair_csv = work / "gatk-ignored-indels.csv"
        native_pair_csv = work / "native-ignored-indels.csv"
        run([str(java), "-Xmx1g", "-jar", str(jar), "AnalyzeCovariates",
             "-before", str(report), "-after", str(ignored_argument_report),
             "-csv", str(gatk_pair_csv)], stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
        run([str(native), "-before", str(report), "-after", str(ignored_argument_report),
             "-csv", str(native_pair_csv)], stdout=subprocess.DEVNULL)
        assert native_pair_csv.read_bytes() == gatk_pair_csv.read_bytes(), (
            "native rejected or changed a Java-compatible ignored argument pair")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["tool"] == "AnalyzeCovariates"
        assert metadata["execution_space"] in {"OpenMP", "Serial"}
        assert metadata["determinism"] == "strict"
        assert metadata["compatibility"]["gatk_report_v1_1"] is True
        assert metadata["compatibility"]["csv"] is True
        assert metadata["telemetry"]["execution_space"]
        assert metadata["telemetry"]["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert metadata["telemetry"]["kernel_execution_policy"] == "RangePolicy"
        assert metadata["telemetry"]["kernel_batches"] == 1
        assert metadata["telemetry"]["kernel_observations"] == metadata["telemetry"]["rows"]
        assert metadata["telemetry"]["csv_bytes"] == native_csv.stat().st_size
        assert metadata["telemetry"]["pdf_bytes"] == native_pdf.stat().st_size
        assert metadata["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in metadata["outputs"])
        print(json.dumps({
            "status": "pass",
            "csv_bit_identical": True,
            "csv_rows": len(gatk_csv.read_text(encoding="utf-8").splitlines()) - 1,
            "before_after_rows": len(pair_lines) - 1,
            "pdf_header": "%PDF-1.4",
            "execution_space": metadata["telemetry"]["execution_space"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
