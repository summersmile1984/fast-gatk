#!/usr/bin/env python3
"""P2 recover-all: multi-fork dangling tails vs GATK 4.6.2.0.

Two incomplete alternate tails share a prefix, then fork.  GATK
``--recover-all-dangling-branches`` recovers each dangling terminal (the
walk still uses the heaviest incoming edge *within* a terminal's chain).
This oracle drives the shipped HaplotypeCaller, not a reimplemented walker.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


# Same topology as the kernel port of ReadThreadingGraphUnitTest.testForkedDanglingEnds.
PREFIX = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
REF_SUFFIX = "GCTAGCTAATCG"
ALT1_SUFFIX = "ACTAGCTAATCG"
ALT2_SUFFIX = "ACTAGATAATCG"
LEFT = "ATGCTAGCTAGATCGTAC" * 5
RIGHT = "CGATACGTAGCTATGCTA" * 5


def records(path: Path) -> list[tuple[str, int, str, str]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append((fields[0], int(fields[1]), fields[3], fields[4]))
    return rows


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout[-1500:]}\nstderr={result.stderr[-3000:]}"
        )
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    if not all(path.exists() for path in (native, java, gatk)):
        oracle_guard.oracle_not_verified("verify_hc_recover_all_gatk_oracle.py", java, gatk)
        raise SystemExit("missing recover-all oracle assets")

    reference_seq = LEFT + PREFIX + REF_SUFFIX + RIGHT
    alt1 = LEFT + PREFIX + ALT1_SUFFIX + RIGHT
    alt2 = LEFT + PREFIX + ALT2_SUFFIX + RIGHT
    contig_len = len(reference_seq)

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-recover-all-") as directory:
        work = Path(directory)
        reference = work / "ref.fa"
        reference.write_text(">chr1\n" + reference_seq + "\n", encoding="utf-8")
        (work / "ref.fa.fai").write_text(
            f"chr1\t{contig_len}\t6\t{contig_len}\t{contig_len + 1}\n", encoding="utf-8")
        header = ["@HD\tVN:1.6\tSO:coordinate",
                  f"@SQ\tSN:chr1\tLN:{contig_len}",
                  "@RG\tID:rg\tSM:S1"]
        lines = list(header)
        read_id = 0
        prefix_start = len(LEFT) + 1
        for sequence in (reference_seq,):
            for copy in range(8):
                lines.append(
                    f"r{read_id}\t0\tchr1\t1\t60\t{contig_len}M\t*\t0\t0\t"
                    f"{sequence}\t{'I' * contig_len}\tRG:Z:rg")
                read_id += 1
        # Incomplete alternate tails: they share PREFIX then fork, and never
        # reach RIGHT, so they dangle instead of rejoining the reference sink.
        for suffix in (ALT1_SUFFIX, ALT2_SUFFIX):
            sequence = PREFIX + suffix
            cigar = f"{len(sequence)}M"
            for copy in range(8):
                lines.append(
                    f"r{read_id}\t0\tchr1\t{prefix_start}\t60\t{cigar}\t*\t0\t0\t"
                    f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg")
                read_id += 1
        sam = work / "input.sam"
        sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
        run([str(java), "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
             "-R", str(reference), "-O", str(work / "ref.dict"),
             "--TRUNCATE_NAMES_AT_WHITESPACE", "true"])
        bam = work / "input.bam"
        run([str(java), "-Xmx1g", "-jar", str(gatk), "SortSam",
             "-I", str(sam), "-O", str(bam), "-SO", "coordinate", "--CREATE_INDEX", "true"])

        common = [
            "-R", str(reference), "-I", str(bam), "-L", f"chr1:1-{contig_len}",
            "--kmer-size", "15", "--dont-increase-kmer-sizes-for-cycles",
            "--min-pruning", "1", "--min-dangling-branch-length", "4",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--allow-non-unique-kmers-in-ref",
        ]
        recover = ["--recover-all-dangling-branches"]

        gatk_default = work / "gatk.default.vcf"
        gatk_all = work / "gatk.recover.vcf"
        native_default = work / "native.default.vcf"
        native_all = work / "native.recover.vcf"
        native_all_manifest = work / "native.recover.json"
        native_default_manifest = work / "native.default.json"

        run([str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
             "-O", str(gatk_default)])
        run([str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common, *recover,
             "-O", str(gatk_all)])
        run([str(native), *common, "--threads", "1", "-O", str(native_default),
             "--output-manifest", str(native_default_manifest)])
        run([str(native), *common, *recover, "--threads", "1", "-O", str(native_all),
             "--output-manifest", str(native_all_manifest)])

        gatk_default_keys = records(gatk_default)
        gatk_all_keys = records(gatk_all)
        native_default_keys = records(native_default)
        native_all_keys = records(native_all)
        if native_all_keys != gatk_all_keys:
            raise AssertionError(
                f"recover-all VCF keys {native_all_keys} != GATK {gatk_all_keys}")
        if native_default_keys != gatk_default_keys:
            raise AssertionError(
                f"default dangling VCF keys {native_default_keys} != GATK {gatk_default_keys}")

        recover_tel = json.loads(native_all_manifest.read_text(encoding="utf-8"))["telemetry"]
        default_tel = json.loads(native_default_manifest.read_text(encoding="utf-8"))["telemetry"]
        if not recover_tel.get("graph_recover_all_dangling_branches"):
            raise AssertionError("recover-all flag not recorded in native manifest")
        if default_tel.get("graph_recover_all_dangling_branches"):
            raise AssertionError("default run unexpectedly recorded recover-all")

        recover_haps = {item["sequence"] for item in recover_tel.get("graph_haplotype_path_details") or []
                        if item.get("sequence")}
        default_haps = {item["sequence"] for item in default_tel.get("graph_haplotype_path_details") or []
                        if item.get("sequence")}
        extra_haps = recover_haps - default_haps
        extra_vcf = set(gatk_all_keys) - set(gatk_default_keys)
        recovered = int(recover_tel.get("graph_dangling_recovered_paths") or 0)
        # The fixture must exercise GATK recover-all on more than one dangling
        # fork.  Matching an empty/default VCF, or recovering a single
        # heaviest-edge walk, is not this contract.
        if not extra_vcf:
            raise AssertionError(
                "GATK recover-all VCF equals default; fixture did not exercise dangling forks")
        if recovered < 2:
            raise AssertionError(
                "recover-all recovered fewer than two dangling forks: "
                f"recovered_paths={recovered} recover_haps={len(recover_haps)} "
                f"default_haps={len(default_haps)} extra_vcf={sorted(extra_vcf)}")
        if len(extra_haps) < 2:
            raise AssertionError(
                "native recover-all did not add two extra haplotypes vs default: "
                f"recover={len(recover_haps)} default={len(default_haps)} "
                f"extra={len(extra_haps)}")

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "gatk_default_vcf": [f"{c}:{p}:{r}>{a}" for c, p, r, a in gatk_default_keys],
            "gatk_recover_all_vcf": [f"{c}:{p}:{r}>{a}" for c, p, r, a in gatk_all_keys],
            "native_recover_all_vcf": [f"{c}:{p}:{r}>{a}" for c, p, r, a in native_all_keys],
            "native_dangling_recovered_paths": recover_tel.get("graph_dangling_recovered_paths"),
            "native_recover_haplotypes": len(recover_haps),
            "native_default_haplotypes": len(default_haps),
            "extra_recover_haplotypes": len(extra_haps),
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
