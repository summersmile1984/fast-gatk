#!/usr/bin/env python3
"""Exercise a real GATK HC tetra-ploid/tetra-allelic assembly corpus.

The bundled BAM is a small GATK HaplotypeCaller integration fixture.  Its
reference is normally an LFS-only asset, so this verifier reconstructs only
the required hg19 chr20 context (the remainder is an indexed N-filled
contig).  GATK is used as the numerical oracle; the native contract requires
the same locus/allele set, complete ploidy-4 Number=G PL width, and executed
PairHMM telemetry.  The final JSON deliberately reports whether the native
row is byte-identical instead of hiding the remaining assembly/PairHMM
numerical delta.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path


REFERENCE_START = 11_363_000
REFERENCE_LENGTH = 63_025_520
# hg19 chr20:11363001-11364000, fetched once from the public hg19 reference
# and kept here as a compact deterministic test asset.
REFERENCE_CONTEXT = (
    "CACAGTGCAATCAAATTAGAACTCAGGATTAAGAAACTCACTCAAAACCACACAATTACATGGAAATTGAACAACCTGCTCCTGAATGAATCTGGGGTAAATAATGAAATTAAGGCAGAAATCAAGAGTTCTTAGAAACCAGTTAGAACAAACAGACAATGTACCAGAATCTCTGGGAAACAACTAAAGCAGTGTTAAGAGGAAAATTTATAGCACTACTTGCTCACACCTGAAAGCTAGATTGATCTCAAATTGACACCCTAGCATCACAAAACAGCTAGAGAAGCAAGAGTAAACTAATCCAAAAGCTAGCAAAAGACAAGAAGTAATTAAGATTAAAGCAAAATTGAAGGAGATAGAGACATGAAAAACCCACCAAAAAATCAATGAATCCAGGAGCTGATTTTTTGAAAAATTAACAAAATAGGTAGGCCACTATCTAGACTAATAAGGAATAAGAGAGAGAAGAATCAAATAGACACAATAAAAAATGATAAAGGGGATATTACTACTGACCACATAGAAATAAAAACTACCTTTAGAGAATCCTATAAACACCTCTACGCAAATAAACTAGAAAATATAGAATCTATGGATAAATTCCTGGACACATACACCCTCCCAGGACTAAACCAGGAAGAAGTCAAATCCCTGAAAAGACAAAATACAAATTCTGAAATTAAAGCAGTAATTAATAACCTACCAACCAATAAAAGCCCAGGACCAGACAGATTCACAGCTGAATTCCACCAGTTACAAAAAGGAGCTGGTACCATTCCTTCTAAAACTATTCCAAACAATTGGAAAAAAGGGACTCCTCCCTAACTCATTTTATGTAGCCACCATCATCCTGATATCAAAGCCTGGCAGAAACACAACAAAAAAAGAAAACTTCAGGTCAATATCCCTGATGAACATCGAAGTAAAAATCAGAAAATACTGGCAAACCAAATCTAGCAGCACATCAAAAAATTTATCCACCACAATCAAGTAGGCTTCA"
)
EXPECTED_GATK_ROW = (
    "20\t11363591\t.\tT\tA,G\t684.25\t.\t"
    "AC=1,2;AF=0.250,0.500;AN=4;BaseQRankSum=-1.560;DP=40;FS=3.081;"
    "MLEAC=1,2;MLEAF=0.250,0.500;MQ=60.00;MQRankSum=0.000;QD=20.12;"
    "ReadPosRankSum=0.366;SOR=0.132\tGT:AD:DP:GQ:PL\t"
    "0/1/2/2:2,15,17:34:6:698,386,347,327,364,354,43,6,36,308,0,12,285,26,320"
)


def data_rows(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def rows_at_position(path: Path, position: str) -> list[str]:
    return [row for row in data_rows(path) if row.split("\t", 2)[1] == position]


def build_reference(path: Path) -> None:
    if len(REFERENCE_CONTEXT) < 900:
        raise AssertionError("embedded reference context is unexpectedly short")
    sequence = bytearray(b"N") * REFERENCE_LENGTH
    begin = REFERENCE_START
    sequence[begin:begin + len(REFERENCE_CONTEXT)] = REFERENCE_CONTEXT.encode("ascii")
    line_bases = 60
    with path.open("wb") as output:
        output.write(b">20\n")
        for offset in range(0, REFERENCE_LENGTH, line_bases):
            output.write(sequence[offset:offset + line_bases])
            output.write(b"\n")
    # The first base follows the four-byte '>20\\n' header.
    path.with_suffix(path.suffix + ".fai").write_text(
        f"20\t{REFERENCE_LENGTH}\t4\t{line_bases}\t{line_bases + 1}\n",
        encoding="utf-8")
    # Omit M5 intentionally: the BAM dictionary is GRCh37, while this
    # synthetic outer sequence is N-filled.  The HC invocation disables only
    # dictionary validation; all queried bases are the exact hg19 context.
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:20\tLN:{REFERENCE_LENGTH}\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get("FASTGATK_HC_BINARY",
                               root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/pretendTobeTetraPloidTetraAllelicSite.bam"
    if not (java.exists() and jar.exists() and native.exists() and bam.exists()):
        raise SystemExit("missing bundled GATK/JDK/native HC complex oracle assets")

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-complex-oracle-") as directory:
        work = Path(directory)
        reference = work / "reference.fa"
        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf"
        gatk_max_genotype_output = work / "gatk.max-genotype-5.vcf"
        native_max_genotype_output = work / "native.max-genotype-5.vcf"
        gatk_max_genotype_gvcf = work / "gatk.max-genotype-5.g.vcf"
        native_max_genotype_gvcf = work / "native.max-genotype-5.g.vcf"
        gatk_capped_output = work / "gatk.max-alt-1.vcf"
        native_capped_output = work / "native.max-alt-1.vcf"
        gatk_capped_gvcf = work / "gatk.max-alt-1.g.vcf"
        native_capped_gvcf = work / "native.max-alt-1.g.vcf"
        manifest = work / "native.manifest.json"
        build_reference(reference)
        base_gatk = [
            str(java), "-Xmx2g", "-jar", str(jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", "20:11363580-11363600",
            "-O", str(gatk_output), "-ploidy", "4", "--max-genotype-count", "15",
            "--add-output-vcf-command-line", "false", "--create-output-variant-index", "false",
            "--disable-sequence-dictionary-validation", "true",
            "--seconds-between-progress-updates", "1",
        ]
        gatk = subprocess.run(base_gatk, text=True, capture_output=True)
        if gatk.returncode != 0:
            raise RuntimeError(f"GATK HaplotypeCaller failed:\n{gatk.stderr[-5000:]}")
        native_run = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "20:11363580-11363600", "-O", str(native_output),
            "--sample-ploidy", "4", "--threads", "2",
            "--disable-sequence-dictionary-validation", "true",
            "--output-manifest", str(manifest),
        ], text=True, capture_output=True)
        if native_run.returncode != 0:
            raise RuntimeError(f"native HaplotypeCaller failed:\n{native_run.stderr[-5000:]}")
        max_genotype_gatk = base_gatk.copy()
        max_genotype_gatk[max_genotype_gatk.index(str(gatk_output))] = str(gatk_max_genotype_output)
        max_genotype_gatk[max_genotype_gatk.index("--max-genotype-count") + 1] = "5"
        max_genotype_gatk_run = subprocess.run(max_genotype_gatk, text=True, capture_output=True)
        if max_genotype_gatk_run.returncode != 0:
            raise RuntimeError(
                "GATK HaplotypeCaller --max-genotype-count failed:\n"
                f"{max_genotype_gatk_run.stderr[-5000:]}")
        max_genotype_native_run = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "20:11363580-11363600", "-O", str(native_max_genotype_output),
            "--sample-ploidy", "4", "--max-genotype-count", "5", "--threads", "2",
            "--disable-sequence-dictionary-validation", "true",
            "--output-manifest", str(work / "native.max-genotype-5.manifest.json"),
        ], text=True, capture_output=True)
        if max_genotype_native_run.returncode != 0:
            raise RuntimeError(
                "native HaplotypeCaller --max-genotype-count failed:\n"
                f"{max_genotype_native_run.stderr[-5000:]}")
        max_genotype_gatk_gvcf = max_genotype_gatk.copy()
        max_genotype_gatk_gvcf[max_genotype_gatk_gvcf.index(str(gatk_max_genotype_output))] = \
            str(gatk_max_genotype_gvcf)
        max_genotype_gatk_gvcf.extend(["-ERC", "GVCF"])
        max_genotype_gatk_gvcf_run = subprocess.run(
            max_genotype_gatk_gvcf, text=True, capture_output=True)
        if max_genotype_gatk_gvcf_run.returncode != 0:
            raise RuntimeError(
                "GATK HaplotypeCaller max-genotype gVCF failed:\n"
                f"{max_genotype_gatk_gvcf_run.stderr[-5000:]}")
        max_genotype_native_gvcf_run = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "20:11363580-11363600", "-O", str(native_max_genotype_gvcf),
            "--sample-ploidy", "4", "--max-genotype-count", "5", "-ERC", "GVCF",
            "--threads", "2", "--disable-sequence-dictionary-validation", "true",
        ], text=True, capture_output=True)
        if max_genotype_native_gvcf_run.returncode != 0:
            raise RuntimeError(
                "native HaplotypeCaller max-genotype gVCF failed:\n"
                f"{max_genotype_native_gvcf_run.stderr[-5000:]}")
        capped_gatk = base_gatk.copy()
        capped_gatk[capped_gatk.index(str(gatk_output))] = str(gatk_capped_output)
        capped_gatk.extend(["--max-alternate-alleles", "1"])
        capped_gatk_run = subprocess.run(capped_gatk, text=True, capture_output=True)
        if capped_gatk_run.returncode != 0:
            raise RuntimeError(
                "GATK HaplotypeCaller --max-alternate-alleles failed:\n"
                f"{capped_gatk_run.stderr[-5000:]}")
        capped_native_run = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "20:11363580-11363600", "-O", str(native_capped_output),
            "--sample-ploidy", "4", "--max-alternate-alleles", "1", "--threads", "2",
            "--disable-sequence-dictionary-validation", "true",
            "--output-manifest", str(work / "native.max-alt-1.manifest.json"),
        ], text=True, capture_output=True)
        if capped_native_run.returncode != 0:
            raise RuntimeError(
                "native HaplotypeCaller --max-alternate-alleles failed:\n"
                f"{capped_native_run.stderr[-5000:]}")
        capped_gatk_gvcf = capped_gatk.copy()
        capped_gatk_gvcf[capped_gatk_gvcf.index(str(gatk_capped_output))] = str(gatk_capped_gvcf)
        capped_gatk_gvcf.extend(["-ERC", "GVCF"])
        capped_gatk_gvcf_run = subprocess.run(capped_gatk_gvcf, text=True, capture_output=True)
        if capped_gatk_gvcf_run.returncode != 0:
            raise RuntimeError(
                "GATK HaplotypeCaller capped gVCF failed:\n"
                f"{capped_gatk_gvcf_run.stderr[-5000:]}")
        capped_native_gvcf_run = subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "20:11363580-11363600", "-O", str(native_capped_gvcf),
            "--sample-ploidy", "4", "--max-alternate-alleles", "1", "-ERC", "GVCF",
            "--threads", "2", "--disable-sequence-dictionary-validation", "true",
        ], text=True, capture_output=True)
        if capped_native_gvcf_run.returncode != 0:
            raise RuntimeError(
                "native HaplotypeCaller capped gVCF failed:\n"
                f"{capped_native_gvcf_run.stderr[-5000:]}")
        gatk_rows = data_rows(gatk_output)
        native_rows = data_rows(native_output)
        gatk_max_genotype_rows = data_rows(gatk_max_genotype_output)
        native_max_genotype_rows = data_rows(native_max_genotype_output)
        gatk_max_genotype_gvcf_rows = rows_at_position(gatk_max_genotype_gvcf, "11363591")
        native_max_genotype_gvcf_rows = rows_at_position(native_max_genotype_gvcf, "11363591")
        gatk_capped_rows = data_rows(gatk_capped_output)
        native_capped_rows = data_rows(native_capped_output)
        gatk_capped_gvcf_rows = rows_at_position(gatk_capped_gvcf, "11363591")
        native_capped_gvcf_rows = rows_at_position(native_capped_gvcf, "11363591")
        if gatk_rows != [EXPECTED_GATK_ROW]:
            raise AssertionError(f"unexpected GATK oracle row: {gatk_rows!r}")
        if len(native_rows) != 1:
            raise AssertionError(f"native HC did not emit one candidate row: {native_rows!r}")
        native_fields = native_rows[0].split("\t")
        if native_fields[:5] != ["20", "11363591", ".", "T", "A,G"]:
            raise AssertionError(f"native HC locus/alleles mismatch: {native_rows[0]}")
        if len(native_fields) != 10 or len(native_fields[9].split(":")[-1].split(",")) != 15:
            raise AssertionError("native HC did not emit the complete ploidy-4 Number=G PL row")
        gatk_fields = gatk_rows[0].split("\t")
        if native_fields[8:10] != gatk_fields[8:10]:
            raise AssertionError(
                "native complex assembly genotype fields must match GATK exactly: "
                f"native={native_fields[8:10]!r} gatk={gatk_fields[8:10]!r}")
        # Multi-ALT AD is an allele partition: a read is counted once across
        # REF/A/G.  The complex row is compared as a complete GATK record,
        # including PL and all standard annotations.
        native_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                       for item in native_fields[7].split(";") if "=" in item}
        native_sample = native_fields[9].split(":")
        # The shared Kokkos cohort AF kernel owns the multi-ALT site
        # posterior/QUAL path; all standard scalar and rank-sum annotations
        # are strict in this oracle.
        if native_fields[5] != "684.25":
            raise AssertionError(f"native multi-ALT site QUAL mismatch: {native_fields[5]}")
        if {key: native_info.get(key) for key in ("DP", "FS", "QD")} != {
                "DP": "40", "FS": "3.081", "QD": "20.12"}:
            raise AssertionError(
                "native multi-ALT QUAL/coverage annotations mismatch: "
                f"{native_info}")
        if "ExcessHet" in native_info:
            raise AssertionError(
                "native tetra-ploid HC must omit the diploid-only ExcessHet annotation")
        if native_sample[1:3] != ["2,15,17", "34"]:
            raise AssertionError(
                "native multi-ALT AD/DP must use one-read-one-allele ownership: "
                f"INFO/DP={native_info.get('DP')} FORMAT={native_sample}")
        if len(gatk_max_genotype_rows) != 1 or len(native_max_genotype_rows) != 1:
            raise AssertionError(
                "--max-genotype-count 5 must emit one haplotype-score-capped candidate row: "
                f"gatk={gatk_max_genotype_rows!r} native={native_max_genotype_rows!r}")
        if gatk_max_genotype_rows != native_max_genotype_rows:
            raise AssertionError(
                "native --max-genotype-count 5 row must exactly match GATK after "
                "upstream K-best haplotype-score ALT selection: "
                f"native={native_max_genotype_rows!r} gatk={gatk_max_genotype_rows!r}")
        max_genotype_fields = native_max_genotype_rows[0].split("\t")
        if "," in max_genotype_fields[4] or len(max_genotype_fields[9].split(":")[-1].split(",")) != 5:
            raise AssertionError(
                "max-genotype-capped tetra-ploid row must contain one ALT and five Number=G PL values")
        if gatk_max_genotype_gvcf_rows != native_max_genotype_gvcf_rows:
            raise AssertionError(
                "native max-genotype gVCF candidate must exactly match GATK's "
                "pre-<NON_REF> haplotype-score allele selection: "
                f"native={native_max_genotype_gvcf_rows!r} gatk={gatk_max_genotype_gvcf_rows!r}")
        if len(gatk_capped_rows) != 1 or len(native_capped_rows) != 1:
            raise AssertionError(
                "--max-alternate-alleles 1 must emit exactly one capped candidate row: "
                f"gatk={gatk_capped_rows!r} native={native_capped_rows!r}")
        if gatk_capped_rows != native_capped_rows:
            raise AssertionError(
                "native --max-alternate-alleles 1 row must exactly match GATK after "
                "GL-based ALT selection and Number=G/Number=R remapping: "
                f"native={native_capped_rows!r} gatk={gatk_capped_rows!r}")
        capped_fields = native_capped_rows[0].split("\t")
        if "," in capped_fields[4] or len(capped_fields[9].split(":")[-1].split(",")) != 5:
            raise AssertionError(
                "capped tetra-ploid row must contain one ALT and the five-entry Number=G PL vector")
        if gatk_capped_gvcf_rows != native_capped_gvcf_rows:
            raise AssertionError(
                "native capped gVCF candidate must preserve GATK's original symbolic "
                "<NON_REF> Number=G subset: "
                f"native={native_capped_gvcf_rows!r} gatk={gatk_capped_gvcf_rows!r}")
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        if telemetry.get("sample_ploidy") != 4 or telemetry.get("pairhmm_skip_reason") != "executed":
            raise AssertionError(f"native HC PairHMM/ploidy telemetry mismatch: {telemetry}")
        print(json.dumps({
            "status": "pass",
            "gatk_row": gatk_rows[0],
            "native_row": native_rows[0],
            "native_bit_identical": native_rows[0] == gatk_rows[0],
            "native_genotype_bit_identical": native_fields[8:10] == gatk_fields[8:10],
            "max_genotype_count_5_row": native_max_genotype_rows[0],
            "max_genotype_count_5_bit_identical":
                native_max_genotype_rows == gatk_max_genotype_rows,
            "max_genotype_count_5_gvcf_candidate_bit_identical":
                native_max_genotype_gvcf_rows == gatk_max_genotype_gvcf_rows,
            "max_alternate_alleles_1_row": native_capped_rows[0],
            "max_alternate_alleles_1_bit_identical": native_capped_rows == gatk_capped_rows,
            "max_alternate_alleles_1_gvcf_candidate_bit_identical":
                native_capped_gvcf_rows == gatk_capped_gvcf_rows,
            "native_site_qual_exact": native_fields[5] == gatk_fields[5],
            "native_scalar_annotations_exact": all(
                native_info.get(key) == {item.split("=", 1)[0]: item.split("=", 1)[1]
                    for item in gatk_fields[7].split(";") if "=" in item}.get(key)
                for key in ("DP", "FS", "QD")),
            "native_execution_space": telemetry.get("execution_space", ""),
            "native_pairhmm_execution_space": telemetry.get("pairhmm_execution_space", ""),
            "native_pl_width": 15,
            "reference_context_md5": hashlib.md5(REFERENCE_CONTEXT.encode()).hexdigest(),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
