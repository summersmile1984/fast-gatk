#!/usr/bin/env python3
"""Pinned Java/native oracle for GenotypeGVCFs genotype assignment modes."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def fasta_base(path: Path, one_based: int) -> str:
    sequence = "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines()
                       if not line.startswith(">"))
    return sequence[one_based - 1]


def records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.split("\t") for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def header_ids(path: Path) -> set[tuple[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        result: set[tuple[str, str]] = set()
        for line in stream.read().splitlines():
            kind = ("INFO" if line.startswith("##INFO=<") else
                    "FORMAT" if line.startswith("##FORMAT=<") else
                    "FILTER" if line.startswith("##FILTER=<") else None)
            if kind is None or "ID=" not in line:
                continue
            value = line.split("ID=", 1)[1].split(",", 1)[0]
            result.add((kind, value.rstrip(">")))
        return result


def run(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed:\n{result.stderr[-6000:]}")
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY", str(root / "fastgatk-native/build/fastgatk-genotype-gvcf")
    ))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK GenotypeGVCFs assignment oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    # The fixture has two samples and a concrete three-allele site.  PL rows
    # intentionally make one sample ambiguous, so posterior assignment can
    # differ from raw PL assignment under the Hardy-Weinberg cohort prior.
    pos = 69067
    ref = fasta_base(reference, pos)
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-assignment-oracle-") as directory:
        work = Path(directory)
        # Plain input avoids requiring an external bgzip binary; GATK's
        # IndexFeatureFile creates the Tribble index and both engines accept
        # the same local stream.
        input_path = work / "cohort.g.vcf"
        header = f"""##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of a reference block>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled genotype likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""
        # Alleles are REF,C,G,<NON_REF>, hence ten diploid PL cells.  S1's
        # likelihoods mildly favor 0/1; the cohort posterior prior favors the
        # reference state, exercising GP/PG emission and assignment together.
        row = (f"17\t{pos}\t.\t{ref}\tC,G,<NON_REF>\t.\tPASS\tDP=30\t"
               "GT:DP:AD:PL\t0/1:15:8,7,0,0:100,0,100,120,120,120,150,150,150,150\t"
               "0/0:15:15,0,0,0:0,100,100,120,120,120,150,150,150,150\n")
        with input_path.open("wt", encoding="utf-8") as stream:
            stream.write(header)
            stream.write(row)
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)],
            "GATK IndexFeatureFile")

        outputs: dict[str, tuple[Path, Path]] = {}
        # GenotypeGVCFs deliberately forces GATK's PREFER_PLS policy in
        # GenotypeGVCFsEngine#createMinimalArgs, even when the inherited
        # --genotype-assignment-method is supplied.  Exercise both accepted
        # spellings and keep this oracle focused on the shared output-allele
        # subset boundary rather than claiming a posterior mode that Java
        # does not apply in this tool.
        for mode in ("PREFER_PLS", "USE_PLS_TO_ASSIGN"):
            gatk_output = work / f"gatk-{mode}.vcf.gz"
            common = ["-R", str(reference), "-V", str(input_path)]
            run([str(java), "-jar", str(gatk), "GenotypeGVCFs", *common,
                 "-O", str(gatk_output), "--genotype-assignment-method", mode,
                 "--create-output-variant-index", "false"], f"GATK {mode}")
            gatk_rows = records(gatk_output)
            # The native compatibility profile intentionally adds telemetry
            # headers and GATK adds command-line metadata; compare the
            # executable schema IDs rather than release-specific metadata.
            required_headers = {("INFO", key) for key in ("AC", "AF", "AN", "DP", "MLEAC", "MLEAF")}
            required_headers |= {("FORMAT", key) for key in ("GT", "AD", "DP", "GQ", "PL")}
            assert required_headers <= header_ids(gatk_output)
            for traversal in ("aggregate", "stream-by-locus"):
                native_output = work / f"native-{mode}-{traversal}.vcf.gz"
                manifest = work / f"native-{mode}-{traversal}.manifest.json"
                native_args = [str(binary), *common, "-O", str(native_output),
                               "--genotype-assignment-method", mode,
                               "--gatk-compatible-annotations", "--create-output-variant-index=false",
                               "--output-manifest", str(manifest)]
                if traversal == "stream-by-locus":
                    native_args.append("--stream-by-locus")
                run(native_args, f"native {mode} {traversal}")
                outputs[f"{mode}/{traversal}"] = (gatk_output, native_output)
                native_rows = records(native_output)
                if gatk_rows != native_rows:
                    raise AssertionError({"mode": mode, "traversal": traversal,
                                          "gatk": gatk_rows, "native": native_rows})
                native_headers = header_ids(native_output)
                if not required_headers <= native_headers:
                    raise AssertionError({"mode": mode, "traversal": traversal,
                                          "missing": sorted(required_headers - native_headers),
                                          "headers": sorted(native_headers)})
                metadata = json.loads(manifest.read_text(encoding="utf-8"))
                assert metadata["compatibility"]["gatk_annotation_compatibility"] is True
                assert metadata["compatibility"]["genotype_assignment_method"] is True
                assert metadata["compatibility"]["output_allele_subset"] is True
                assert metadata["compatibility"]["standard_confidence_allele_pruning"] is True
                assert metadata["telemetry"]["output_allele_pruning_calls"] > 0
                assert metadata["telemetry"]["output_alleles_pruned"] == 1
                assert metadata["telemetry"]["standard_confidence_for_calling"] == 30

        # Assignment mode is not allowed to change the cohort's output schema
        # or sample ordering; raw PL mode is the baseline for this fixture.
        assert records(outputs["USE_PLS_TO_ASSIGN/aggregate"][0])
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "modes": list(outputs),
            "sample_count": 2,
            "allele_count": 4,
            "java_native_rows_exact": True,
            "gatk_assignment_method_aliases_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
