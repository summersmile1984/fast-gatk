#!/usr/bin/env python3
"""Pinned GATK/native oracle for joint GenotypeGVCFs spanning-deletion cleanup.

GATK treats ``*`` as a structural non-variant only when a concrete deletion
record owns the locus.  An orphan ``*`` is removed at the final allele-subset
boundary and any genotype containing it becomes a no-call.  This fixture keeps
both the orphan and supported cases in scope, while comparing the aggregate
and stream-by-locus native merge paths against GATK 4.6.2.0.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def rows(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def run(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")


HEADER = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled genotype likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""

ORPHAN_ROW = (
    "17\t69067\t.\tT\tC,*,<NON_REF>\t.\tPASS\tDP=30\t"
    "GT:DP:AD:PL\t0/1:15:8,7,0,0:100,0,100,120,120,120,150,150,150,150\t"
    "0/2:15:8,0,7,0:100,120,100,120,0,120,100,120,100,120\n"
)

SUPPORTED_ROWS = (
    "17\t69066\t.\tAT\tA,<NON_REF>\t.\tPASS\tDP=30\t"
    "GT:DP:AD:PL\t0/1:15:8,7,0:100,0,100,120,120,120\t"
    "0/0:15:15,0,0:0,100,100,120,120,120\n"
    "17\t69067\t.\tT\tC,*,<NON_REF>\t.\tPASS\tDP=30\t"
    "GT:DP:AD:PL\t0/0:15:15,0,0,0:0,100,100,120,120,120,150,150,150,150\t"
    "0/2:15:8,0,7,0:100,120,100,120,0,120,100,120,100,120\n"
)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY",
        str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, native, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK spanning-deletion oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-spanning-deletion-oracle-") as directory:
        work = Path(directory)
        case_rows = {"orphan": ORPHAN_ROW, "supported": SUPPORTED_ROWS}
        case_summary: dict[str, object] = {}
        for case, body in case_rows.items():
            source = work / f"{case}.g.vcf"
            source.write_text(HEADER + body, encoding="utf-8")
            run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(source)],
                f"GATK IndexFeatureFile ({case})")
            gatk_output = work / f"gatk-{case}.vcf.gz"
            run([str(java), "-Xmx1g", "-jar", str(gatk), "GenotypeGVCFs",
                 "-R", str(reference), "-V", str(source),
                 "--create-output-variant-index", "false", "-O", str(gatk_output)],
                f"GATK GenotypeGVCFs ({case})")
            expected = rows(gatk_output)
            assert expected, {"case": case, "reason": "GATK emitted no rows"}
            if case == "orphan":
                fields = expected[0].split("\t")
                assert fields[4] == "C", expected[0]
                assert fields[9].startswith("0/1:"), expected[0]
                assert fields[10].startswith("./.:"), expected[0]
                assert "AC=1;AF=0.500;AN=2" in fields[7], expected[0]
            else:
                assert len(expected) == 2, expected
                assert "C,*" in expected[1].split("\t")[4], expected[1]

            traversals: dict[str, bool] = {}
            telemetry: dict[str, int] = {}
            for traversal in ("aggregate", "stream-by-locus"):
                output = work / f"native-{case}-{traversal}.vcf.gz"
                manifest = work / f"native-{case}-{traversal}.manifest.json"
                command = [str(native), "-R", str(reference), "-V", str(source),
                           "--gatk-compatible-annotations",
                           "--create-output-variant-index=false", "-O", str(output),
                           "--output-manifest", str(manifest)]
                if traversal == "stream-by-locus":
                    command.append("--stream-by-locus")
                run(command, f"native GenotypeGVCFs ({case}/{traversal})")
                actual = rows(output)
                assert actual == expected, {
                    "case": case, "traversal": traversal,
                    "gatk": expected, "native": actual}
                metadata = json.loads(manifest.read_text(encoding="utf-8"))
                assert metadata["compatibility"]["spanning_deletion_orphan_cleanup"] is True
                telemetry[traversal] = metadata["telemetry"]["orphan_spanning_deletion_loci"]
                traversals[traversal] = True
            assert telemetry["aggregate"] == telemetry["stream-by-locus"]
            assert telemetry["aggregate"] == (1 if case == "orphan" else 0)
            case_summary[case] = {
                "rows": len(expected),
                "aggregate_exact": traversals["aggregate"],
                "stream_by_locus_exact": traversals["stream-by-locus"],
                "orphan_spanning_deletion_loci": telemetry["aggregate"],
            }

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "orphan_star_removed_and_genotype_nocall": True,
        "supported_star_preserved": True,
        "aggregate_stream_exact": True,
        "cases": case_summary,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
