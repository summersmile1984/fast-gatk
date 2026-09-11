#!/usr/bin/env python3
"""GATK oracle for HC post-PairHMM read disqualification and PLs."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def record(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "17" and fields[1] == "69067":
            values = dict(zip(fields[8].split(":"), fields[9].split(":")))
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            return values, info
    raise AssertionError("missing chr17 sentinel call")


def record_line(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "17" and fields[1] == "69067":
            return line
    raise AssertionError("missing chr17 sentinel call")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, java, gatk, bam, reference)
    if not all(path.exists() for path in required):
        oracle_guard.oracle_not_verified('verify_hc_likelihood_filter.py', java, gatk)
        missing = [str(path) for path in required if not path.exists()]
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit(f"missing required HC GATK oracle inputs: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-likelihood-filter-") as directory:
        work = Path(directory)
        # The production HC default is the GATK 30-phred call threshold.  The
        # reference-backed PairHMM marginalization now places this bounded
        # singleton above that gate, so the default path must retain the call
        # rather than silently applying the old approximate-PL expectation.
        default_output = work / "native-default.vcf"
        default_manifest = work / "native-default.manifest.json"
        default_command = [
            str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(default_output), "--threads", "2", "--min-depth", "1",
            "--min-alt-support", "1", "--output-manifest", str(default_manifest),
        ]
        subprocess.run(default_command, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        default_values, default_info = record(default_output)
        assert default_values["GT"] == "1/1"
        assert default_values["GQ"] == "3"
        assert default_values["DP"] == "1"
        assert default_values["AD"] == "0,1"
        assert default_info["DP"] == "1"
        default_telemetry = json.loads(default_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert default_telemetry["standard_confidence_for_calling"] == 30
        output = work / "native.vcf"
        manifest = work / "native.manifest.json"
        gatk_output = work / "gatk.vcf"
        gatk_command = [
            str(java), "-jar", str(gatk), "HaplotypeCaller",
            "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(gatk_output), "--min-pruning", "1",
            "--standard-min-confidence-threshold-for-calling", "0",
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ]
        gatk_run = subprocess.run(gatk_command, text=True, capture_output=True)
        assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
        command = [
            str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(output), "--threads", "2", "--min-depth", "1",
            "--min-alt-support", "1", "--standard-min-confidence-threshold-for-calling", "0",
            "--output-manifest", str(manifest),
        ]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        values, info = record(output)
        gatk_values, gatk_info = record(gatk_output)
        native_line = record_line(output)
        gatk_line = record_line(gatk_output)
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        # The sentinel has PL 43,3,0.  GATK HaplotypeCaller defaults to
        # USE_PLS_TO_ASSIGN, so the native default must remain 1/1 even though
        # the assuming-HW prior is available to the AF/model stage.
        flat_output = work / "native-flat.vcf"
        flat_manifest = work / "native-flat.manifest.json"
        flat_command = [
            str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(flat_output), "--threads", "2", "--min-depth", "1",
            "--min-alt-support", "1", "--standard-min-confidence-threshold-for-calling", "0",
            "--no-genotype-priors",
            "--output-manifest", str(flat_manifest),
        ]
        subprocess.run(flat_command,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        flat_values, _ = record(flat_output)
    assert values["GT"] == "1/1"
    assert values["GQ"] == "3"
    assert values["DP"] == "1"
    assert values["AD"] == "0,1"
    assert info["DP"] == "1"
    assert native_line == gatk_line
    assert values == gatk_values
    assert info == gatk_info
    assert flat_values["GT"] == values["GT"] == "1/1"
    assert flat_values["GQ"] == values["GQ"]
    assert telemetry["genotype_priors"] is True
    assert telemetry["genotype_assignment_method"] == "USE_PLS_TO_ASSIGN"
    assert telemetry["joint_genotype_priors_used"] is False
    assert telemetry["pairhmm_reads_disqualified"] >= 1
    assert telemetry["pairhmm_reads_clipped"] >= 1
    assert telemetry["pairhmm_read_disqualification_threshold"] > 0
    assert telemetry["requested_batch_records"] >= telemetry["effective_batch_records"] >= 1
    assert telemetry["adaptive_batch_reductions"] >= 0
    print(json.dumps({"status": "pass", "tool": "HaplotypeCaller",
                      "pairhmm_reads_disqualified": telemetry["pairhmm_reads_disqualified"],
                      "pairhmm_reads_clipped": telemetry["pairhmm_reads_clipped"],
                      "sentinel": {key: values[key] for key in ("GT", "GQ", "DP", "AD")},
                      "pl": values["PL"], "source_record_bit_identical": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
