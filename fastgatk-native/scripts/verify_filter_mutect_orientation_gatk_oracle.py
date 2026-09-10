#!/usr/bin/env python3
"""Pin the FilterMutectCalls orientation-quality writer to GATK 4.6.2.0.

The orientation posterior itself is a release-specific model, but its VCF
quality encoding is a stable GATK boundary: QualityUtils.errorProbToQual()
bounds the phred quality to [1, 93].  In particular, a zero artifact
posterior is written as ROQ=93, not an unbounded sentinel such as 1000.
This oracle also exercises a reverse-complement prior and an MNV so the
orientation filter is proven to consume the same prior-table wire format.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def read_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, str] = {}
            if fields[7] not in {"", "."}:
                for token in fields[7].split(";"):
                    if "=" in token:
                        key, value = token.split("=", 1)
                        info[key] = value
            records.append({"pos": fields[1], "filter": fields[6], "roq": info.get("ROQ", "")})
    return records


def reverse_complement(context: str) -> str:
    return context.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def prior_table() -> str:
    columns = [
        "context", "rev_comp", "f1r2_a", "f1r2_c", "f1r2_g", "f1r2_t",
        "f2r1_a", "f2r1_c", "f2r1_g", "f2r1_t", "hom_ref", "germline_het",
        "somatic_het", "hom_var", "num_examples", "num_alt_examples",
    ]
    mapping = [7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10, 11]
    rows: list[str] = []
    seen: set[str] = set()
    for left in "ACGT":
        for middle in "ACGT":
            for right in "ACGT":
                context = left + middle + right
                reverse = reverse_complement(context)
                canonical = min(context, reverse)
                if canonical in seen:
                    continue
                seen.add(canonical)

                # Strong F1R2_C prior in canonical AAA.  GATK's
                # ArtifactPrior.getReverseComplement() maps this to F2R1_G
                # in TTT, the context used by the pinned reference fixture.
                canonical_prior = [0.0] * 12
                canonical_prior[1] = 0.90
                canonical_prior[5] = 0.01
                canonical_prior[8:] = [0.01, 0.01, 0.03, 0.01]
                total = sum(canonical_prior)
                canonical_prior = [value / total for value in canonical_prior]
                reverse_prior = [0.0] * 12
                for source, target in enumerate(mapping):
                    reverse_prior[target] = canonical_prior[source]

                rows.append("\t".join([
                    canonical, reverse,
                    *[f"{value:.12g}" for value in canonical_prior], "10", "10",
                ]))
                rows.append("\t".join([
                    reverse, canonical,
                    *[f"{value:.12g}" for value in reverse_prior], "10", "10",
                ]))
    assert len(rows) == 64
    return "#<METADATA>SAMPLE=TUMOR\n" + "\t".join(columns) + "\n" + "\n".join(rows) + "\n"


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    required = (java, gatk, reference, binary)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK orientation oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK orientation oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-orientation-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=Count of reads in F1R2 pair orientation supporting each allele>
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=Count of reads in F2R1 pair orientation supporting each allele>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
        body = (
            # Reference is T at this pinned all-T fixture.  The reverse
            # complement prior has artifact mass for G, not C, so these two
            # records are the zero-posterior/ROQ=93 guard.
            "17\t69000\t.\tT\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,10:20,0\n"
            "17\t69000\t.\tT\tC\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,5:20,5\n"
            # The second base is T>G and therefore exercises the reverse-
            # complement F2R1_G artifact state in the same prior table.
            "17\t69000\t.\tTT\tCG\t.\tPASS\tTLOD=20\tGT:AD:F1R2:F2R1\t0/1:20,10:20,10:20,0\n"
        )
        input_path.write_text(header + body, encoding="utf-8")
        index_run = run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)])
        assert index_run.returncode == 0, index_run.stderr
        stats = work / "stats.table"
        stats.write_text(
            "#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n",
            encoding="utf-8",
        )
        prior = work / "TUMOR.orientation_priors.tar.gz"
        with tarfile.open(prior, "w:gz") as archive:
            payload = prior_table().encode("utf-8")
            member = tarfile.TarInfo("TUMOR.orientation_priors")
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))

        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        common = [
            "--orientation-bias-artifact-priors", str(prior),
            "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.5",
        ]
        gatk_run = run([
            str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
            "-V", str(input_path), "-O", str(gatk_output), "--stats", str(stats),
            "--lenient", "true", *common,
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = run([
            str(binary), "-R", str(reference), "-V", str(input_path),
            "-O", str(native_output), "--stats", str(stats), *common,
        ])
        assert native_run.returncode == 0, native_run.stderr

        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)
        assert len(gatk_records) == len(native_records) == 3
        assert [record["pos"] for record in gatk_records] == [record["pos"] for record in native_records]
        assert [record["roq"] for record in gatk_records[:2]] == ["93", "93"], gatk_records
        assert [record["roq"] for record in native_records[:2]] == ["93", "93"], native_records
        assert "orientation" not in gatk_records[0]["filter"]
        assert "orientation" not in gatk_records[1]["filter"]
        assert "orientation" in gatk_records[2]["filter"]
        assert "orientation" not in native_records[0]["filter"]
        assert "orientation" not in native_records[1]["filter"]
        assert "orientation" in native_records[2]["filter"]
        # No native orientation quality may exceed GATK's maximum phred
        # quality, including the nonzero MNV posterior.
        native_quality = [float(record["roq"]) for record in native_records]
        assert all(1.0 <= quality <= 93.0 for quality in native_quality), native_records
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": len(gatk_records),
            "gatk_roq": [record["roq"] for record in gatk_records],
            "native_roq": [record["roq"] for record in native_records],
            "zero_posterior_roq_93_exact": [record["roq"] for record in native_records[:2]] == ["93", "93"],
            "orientation_filter_pattern_exact": True,
            "orientation_quality_bounded_1_93": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
