#!/usr/bin/env python3
"""Local regression for the two-pass VQSR tranche joint (scatter -> gather).

Scenario (all-native, no Java oracle):
  1. Train a small VariantRecalibrator model on a synthetic 10-variant input.
  2. Split the input into two shards; score each shard against the SAME model
     with --output-tranches-for-scatter and identical requested VQSLOD slices
     -> two v6 scatter tranche CSVs.
  3. GatherTranches merges the two CSVs into one v5 joint tranche file;
     the merge must be deterministic (byte-identical on rerun).
  4. ApplyVQSR assigns FILTERs using the gathered joint tranches; the
     application must be deterministic.

This pins the native scatter/gather/apply wiring (Track B2 local scope) and
the tranche-walk determinism without asserting GATK byte parity.
"""
from __future__ import annotations

import gzip
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
VR = ROOT / "fastgatk-native" / "build" / "fastgatk-variant-recalibrator"
GATHER = ROOT / "fastgatk-native" / "build" / "fastgatk-gather-tranches"
APPLY = ROOT / "fastgatk-native" / "build" / "fastgatk-apply-vqsr"

HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)

APPLY_HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=x>\n"
    "##INFO=<ID=culprit,Number=1,Type=String,Description=x>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"command failed ({result.returncode}):\n"
                             f"{' '.join(command)}\n{result.stderr[-1200:]}")
    return result


def main() -> int:
    if not all(path.is_file() for path in (VR, GATHER, APPLY)):
        raise SystemExit("native VQSR binaries are required")
    import tempfile
    with tempfile.TemporaryDirectory(prefix="vqsr-scatter-joint-") as temporary:
        work = pathlib.Path(temporary)
        base = [
            f"chr1\t{i}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\tGT\t0/1"
            for i, (qd, mq) in enumerate(
                [(30, 60), (28, 58), (25, 55), (24, 54), (23, 53), (22, 52),
                 (20, 50), (18, 48), (15, 45), (12, 42)], 1)
        ]
        write_vcf(work / "input.vcf", base)
        write_vcf(work / "shard1.vcf", base[:5])
        write_vcf(work / "shard2.vcf", base[5:])
        write_vcf(work / "training.vcf", base[:5])
        write_vcf(work / "known.vcf", base[:3])
        resources = [
            "--resource:truth,training=true,truth=true,known=false,prior=15.0",
            str(work / "training.vcf"),
            "--resource:known,training=false,truth=false,known=true",
            str(work / "known.vcf"),
        ]
        common = ["-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "2",
                  "--k-means-iterations", "20", "--bad-lod-score-cutoff", "100.0",
                  "--create-output-variant-index", "false"]
        slices = ["--vqslod-tranche", "10.0", "--vqslod-tranche", "5.0",
                  "--vqslod-tranche", "0.0", "--vqslod-tranche", "-5.0",
                  "--vqslod-tranche", "-10.0"]
        model = work / "model.report"
        run([str(VR), "-V", str(work / "input.vcf"), *resources, *common,
             "-O", str(work / "recal.vcf.gz"), "--tranches-file",
             str(work / "tranches.txt"), "--output-model", str(model)])

        scatter_paths = []
        for shard in ("shard1", "shard2"):
            scatter = work / f"{shard}_scatter.csv"
            run([str(VR), "-V", str(work / f"{shard}.vcf"), "--input-model",
                 str(model), *resources, *common, "-O",
                 str(work / f"{shard}_recal.vcf.gz"), "--tranches-file",
                 str(scatter), "--output-tranches-for-scatter", *slices])
            text = scatter.read_text(encoding="utf-8")
            assert "Version number 6" in text, scatter
            assert "requestedVQSLOD" in text, scatter
            scatter_paths.append(str(scatter))

        gathered = work / "gathered.csv"
        run([str(GATHER), "-I", scatter_paths[0], "-I", scatter_paths[1],
             "-O", str(gathered), "--mode", "SNP"])
        rerun = work / "gathered2.csv"
        run([str(GATHER), "-I", scatter_paths[0], "-I", scatter_paths[1],
             "-O", str(rerun), "--mode", "SNP"])
        assert gathered.read_bytes() == rerun.read_bytes(), "gather not deterministic"
        joint = gathered.read_text(encoding="utf-8")
        assert "Version number 5" in joint or "requestedVQSLOD" in joint

        # ApplyVQSR with the gathered joint tranches on a small VQSLOD-tagged
        # input; application must be deterministic and every row filtered.
        apply_input = work / "apply_in.vcf"
        apply_input.write_text(
            APPLY_HEADER + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.5;culprit=QD\n"
            "chr1\t2\t.\tA\tG\t.\tPASS\tVQSLOD=0.8;culprit=QD\n"
            "chr1\t3\t.\tA\tG\t.\tPASS\tVQSLOD=-2.0;culprit=MQ\n"
            "chr1\t4\t.\tA\tG\t.\tPASS\tVQSLOD=-30.0;culprit=MQ\n",
            encoding="utf-8")
        run([str(APPLY), "-V", str(apply_input), "-O", str(work / "a1.vcf.gz"),
             "--recal-file", str(work / "recal.vcf.gz"), "--tranches-file",
             str(gathered), "--mode", "SNP", "--lod-score-cutoff", "0",
             "--create-output-variant-index", "false"])
        run([str(APPLY), "-V", str(apply_input), "-O", str(work / "a2.vcf.gz"),
             "--recal-file", str(work / "recal.vcf.gz"), "--tranches-file",
             str(gathered), "--mode", "SNP", "--lod-score-cutoff", "0",
             "--create-output-variant-index", "false"])
        with gzip.open(work / "a1.vcf.gz", "rt") as handle:
            a1 = handle.read()
        with gzip.open(work / "a2.vcf.gz", "rt") as handle:
            a2 = handle.read()
        assert a1 == a2, "apply-vqsr not deterministic"
        assert a1.count("\tPASS\t") + a1.count("\tLOW_VQSLOD\t") >= 4, a1
        print('{"suite":"vqsr-scatter-joint","scatter_shards":2,'
              '"gather_deterministic":true,"apply_deterministic":true,'
              '"status":"pass"}')
    return 0


if __name__ == "__main__":
    sys.exit(main())
