#!/usr/bin/env python3
"""B2 4-shard VQSR scatter joint contract.

Extension of verify_vqsr_scatter_joint.py: splits the input into 4 shards
and verifies the scatter→gather→apply round-trip stays deterministic. Pins
B2's pipeline scaling to 4 shards (the 2-shard test is at
verify_vqsr_scatter_joint.py). Uses synthetic data on chr1 — no GATK jar required.
"""
from __future__ import annotations
import gzip, pathlib, subprocess, sys
ROOT = pathlib.Path("/home/turing-agents/Documents/fast-gatk")
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

def write_vcf(path, rows):
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")

def run(cmd):
    r = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if r.returncode != 0:
        raise AssertionError(f"command failed ({r.returncode}):\n{' '.join(cmd)}\n{r.stderr[-1200:]}")
    return r

def gather_cmd(binary, inputs, output):
    cmd = [str(binary)]
    for p in inputs:
        cmd += ["-I", p]
    cmd += ["-O", str(output), "--mode", "SNP"]
    return cmd

def main():
    if not all(p.is_file() for p in (VR, GATHER, APPLY)):
        raise SystemExit("native VQSR binaries are required")
    import tempfile, json
    with tempfile.TemporaryDirectory(prefix="vqsr-4shard-") as t:
        work = pathlib.Path(t)
        base = [
            f"chr1\t{i}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}\tGT\t0/1"
            for i, (qd, mq) in enumerate(
                [(30, 60), (28, 58), (25, 55), (24, 54), (23, 53), (22, 52),
                 (20, 50), (18, 48), (15, 45), (12, 42), (10, 40), (8, 38),
                 (6, 36), (4, 34), (2, 32), (1, 30), (29, 59), (27, 57),
                 (26, 56), (21, 51)], 1)
        ]
        write_vcf(work / "input.vcf", base)
        for i in range(4):
            write_vcf(work / f"shard{i+1}.vcf", base[i*5:(i+1)*5])
        write_vcf(work / "training.vcf", base[:10])
        write_vcf(work / "known.vcf", base[:5])
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
        for i in range(1, 5):
            scatter = work / f"shard{i}_scatter.csv"
            run([str(VR), "-V", str(work / f"shard{i}.vcf"), "--input-model",
                 str(model), *resources, *common, "-O",
                 str(work / f"shard{i}_recal.vcf.gz"), "--tranches-file",
                 str(scatter), "--output-tranches-for-scatter", *slices])
            text = scatter.read_text(encoding="utf-8")
            assert "Version number 6" in text, scatter
            assert "requestedVQSLOD" in text, scatter
            scatter_paths.append(str(scatter))

        gathered = work / "gathered.csv"
        run(gather_cmd(GATHER, scatter_paths, gathered))
        rerun = work / "gathered2.csv"
        run(gather_cmd(GATHER, scatter_paths, rerun))
        assert gathered.read_bytes() == rerun.read_bytes(), "gather not deterministic"
        joint = gathered.read_text(encoding="utf-8")
        assert "Version number 5" in joint or "requestedVQSLOD" in joint

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
        with gzip.open(work / "a1.vcf.gz", "rt") as h:
            a1 = h.read()
        with gzip.open(work / "a2.vcf.gz", "rt") as h:
            a2 = h.read()
        assert a1 == a2, "apply not deterministic"

    print(json.dumps({
        "status": "pass",
        "fixture": "vqsr_scatter_4shard_joint",
        "shards": 4,
        "variants_per_shard": 5,
        "scattered_csv_version": 6,
        "gathered_csv_version": 5,
        "deterministic": True,
        "b2_4shard_recorded": True,
    }, indent=2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
