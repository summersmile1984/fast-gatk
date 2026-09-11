#!/usr/bin/env python3
"""Contract and GATK-oracle checks for the native FASTA walkers."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
BUILD = Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
FASTA_BIN = BUILD / "fastgatk-fasta-reference-maker"
ALTERNATE_BIN = BUILD / "fastgatk-fasta-alternate-reference-maker"
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK_JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if check and result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout}\nstderr={result.stderr}"
        )
    return result


def reference_sequence() -> str:
    return "".join(line.strip() for line in REF.read_text().splitlines() if not line.startswith(">"))


def write_vcf(path: Path, sequence: str, *, offset: int = 0) -> None:
    def base(position: int) -> str:
        return sequence[position - 1]

    # One SNP, one anchored insertion and one anchored deletion.  The REF
    # alleles are taken from the actual reference so both Java and native paths
    # exercise the same simple-variant branch without reference mismatch noise.
    snp_pos = 100 + offset
    ins_pos = 200 + offset
    del_pos = 300 + offset
    snp_ref = base(snp_pos)
    snp_alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[snp_ref.upper()]
    ins_ref = base(ins_pos)
    del_ref = base(del_pos) + base(del_pos + 1)
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=17,length=1000000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"17\t{snp_pos}\t.\t{snp_ref}\t{snp_alt}\t50\tPASS\t.\tGT\t0/1\n"
        f"17\t{ins_pos}\t.\t{ins_ref}\t{ins_ref}T\t50\tPASS\t.\tGT\t1/1\n"
        f"17\t{del_pos}\t.\t{del_ref}\t{del_ref[0]}\t50\tPASS\t.\tGT\t1/1\n"
    )


def read_manifest(path: Path) -> dict:
    return json.loads(path.read_text())


def compare_outputs(native: Path, oracle: Path) -> None:
    native_bytes = native.read_bytes()
    oracle_bytes = oracle.read_bytes()
    if native_bytes != oracle_bytes:
        limit = min(len(native_bytes), len(oracle_bytes))
        first = next((index for index in range(limit)
                      if native_bytes[index] != oracle_bytes[index]), limit)
        native_seq = "".join(line.strip() for line in native.read_text().splitlines() if not line.startswith(">"))
        oracle_seq = "".join(line.strip() for line in oracle.read_text().splitlines() if not line.startswith(">"))
        sequence_mismatches = [(index + 1, lhs, rhs)
                               for index, (lhs, rhs) in enumerate(zip(native_seq, oracle_seq))
                               if lhs != rhs][:12]
        raise AssertionError(
            f"FASTA output differs at byte {first} (native_len={len(native_bytes)}, "
            f"oracle_len={len(oracle_bytes)}):\n native={native_bytes[max(0, first-32):first+64]!r}\n"
            f" oracle={oracle_bytes[max(0, first-32):first+64]!r}\n"
            f" sequence_mismatches={sequence_mismatches}"
        )


def main() -> int:
    assert REF.exists() and (REF.with_suffix(REF.suffix + ".fai")).exists()
    assert FASTA_BIN.exists() and ALTERNATE_BIN.exists()
    java_oracle = oracle_guard.oracle_ready('verify_fasta_reference_tools.py', JAVA, GATK_JAR)
    with tempfile.TemporaryDirectory(prefix="fastgatk-fasta-") as directory:
        work = Path(directory)
        native = work / "reference.native.fasta"
        manifest = work / "reference.native.json"
        run([str(FASTA_BIN), "-R", str(REF), "-L", "17:100-260", "-L", "17:261-360",
             "-O", str(native), "--line-width", "37", "--threads", "2",
             "--output-manifest", str(manifest)])
        assert native.with_suffix(native.suffix + ".fai").exists()
        assert native.with_suffix(".dict").exists()
        payload = read_manifest(manifest)
        assert payload["tool"] == "FastaReferenceMaker"
        assert payload["records"] == 261
        assert payload["execution_space"]
        assert payload["execute_seconds"] >= 0.0

        if java_oracle:
            oracle = work / "reference.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaReferenceMaker", "-R", str(REF),
                 "-L", "17:100-260", "-L", "17:261-360", "-O", str(oracle),
                 "--line-width", "37"])
            compare_outputs(native, oracle)
            compare_outputs(native.with_suffix(".dict"), oracle.with_suffix(".dict"))

        sequence = reference_sequence()
        variants = work / "variants.vcf"
        write_vcf(variants, sequence)
        if java_oracle:
            run([str(JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(variants)])
        alternate = work / "alternate.native.fasta"
        alternate_manifest = work / "alternate.native.json"
        run([str(ALTERNATE_BIN), "-R", str(REF), "-V", str(variants), "-L", "17:1-1000",
             "-O", str(alternate), "--line-width", "41", "--output-manifest",
             str(alternate_manifest)])
        assert alternate.with_suffix(alternate.suffix + ".fai").exists()
        assert alternate.with_suffix(".dict").exists()
        alternate_payload = read_manifest(alternate_manifest)
        assert alternate_payload["tool"] == "FastaAlternateReferenceMaker"
        assert alternate_payload["records"] == 1000
        assert alternate_payload["execute_seconds"] >= 0.0

        if java_oracle:
            alternate_oracle = work / "alternate.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaAlternateReferenceMaker", "-R", str(REF),
                 "-V", str(variants), "-L", "17:1-1000", "-O", str(alternate_oracle),
                 "--line-width", "41"])
            compare_outputs(alternate, alternate_oracle)
            compare_outputs(alternate.with_suffix(".dict"), alternate_oracle.with_suffix(".dict"))

        # Verify the explicit mask-priority branch and IUPAC sample branch on
        # the same small interval.  The mask VCF has a SNP at a different
        # locus, while the input VCF has a heterozygous SNP at snp_pos.
        mask = work / "mask.vcf"
        mask_pos = 150
        mask_ref = sequence[mask_pos - 1]
        mask_alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[mask_ref.upper()]
        mask.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"17\t{mask_pos}\t.\t{mask_ref}\t{mask_alt}\t50\tPASS\t.\n"
        )
        if java_oracle:
            run([str(JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(mask)])
        masked = work / "alternate.masked.native.fasta"
        run([str(ALTERNATE_BIN), "-R", str(REF), "-V", str(variants), "--snp-mask", str(mask),
             "--snp-mask-priority", "-L", "17:1-1000", "-O", str(masked)])
        assert masked.exists()
        if java_oracle:
            masked_oracle = work / "alternate.masked.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaAlternateReferenceMaker", "-R", str(REF),
                 "-V", str(variants), "--snp-mask", str(mask), "--snp-mask-priority",
                 "-L", "17:1-1000", "-O", str(masked_oracle)])
            compare_outputs(masked, masked_oracle)

        # GATK's SNP mask is a feature mask, not a callset filter: its
        # isMasked() helper checks VariantContext::isSNP() and intentionally
        # does not drop FILTERed records.  Keep this regression separate from
        # the PASS mask above so the native reader cannot accidentally reuse
        # the source-call filtering policy.
        filtered_mask = work / "mask.filtered.vcf"
        filtered_pos = 160
        filtered_ref = sequence[filtered_pos - 1]
        filtered_alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[filtered_ref.upper()]
        filtered_mask.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "##FILTER=<ID=LowQual,Description=synthetic filtered SNP>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"17\t{filtered_pos}\t.\t{filtered_ref}\t{filtered_alt}\t10\tLowQual\t.\n"
        )
        if java_oracle:
            run([str(JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(filtered_mask)])
        filtered_native = work / "alternate.filtered-mask.native.fasta"
        run([str(ALTERNATE_BIN), "-R", str(REF), "-V", str(variants),
             "--snp-mask", str(filtered_mask), "--snp-mask-priority",
             "-L", "17:1-1000", "-O", str(filtered_native)])
        if java_oracle:
            filtered_oracle = work / "alternate.filtered-mask.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaAlternateReferenceMaker", "-R", str(REF),
                 "-V", str(variants), "--snp-mask", str(filtered_mask), "--snp-mask-priority",
                 "-L", "17:1-1000", "-O", str(filtered_oracle)])
            compare_outputs(filtered_native, filtered_oracle)

        # A mixed SNP+indel record is VariantContext.Type.MIXED, so neither
        # the simple SNP nor simple-indel branch in GATK applies.  This guards
        # the record-level classification retained by the native ALT flattener.
        mixed_pos = 170
        mixed_ref = sequence[mixed_pos - 1]
        mixed_snp = {"A": "C", "C": "G", "G": "T", "T": "A"}[mixed_ref.upper()]
        mixed = work / "mixed.vcf"
        mixed.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            f"17\t{mixed_pos}\t.\t{mixed_ref}\t{mixed_snp},{mixed_ref}T\t50\tPASS\t.\tGT\t1/2\n"
        )
        if java_oracle:
            run([str(JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(mixed)])
        mixed_native = work / "alternate.mixed.native.fasta"
        run([str(ALTERNATE_BIN), "-R", str(REF), "-V", str(mixed),
             "-L", "17:1-1000", "-O", str(mixed_native)])
        if java_oracle:
            mixed_oracle = work / "alternate.mixed.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaAlternateReferenceMaker", "-R", str(REF),
                 "-V", str(mixed), "-L", "17:1-1000", "-O", str(mixed_oracle)])
            compare_outputs(mixed_native, mixed_oracle)

        # IUPAC selection is diploid, but it is not limited to biallelic VCF
        # records.  A 1/2 genotype in a pure multi-ALT SNP record must encode
        # the two selected concrete bases (rather than silently falling back
        # to the first ALT).
        multi_pos = 180
        multi_ref = sequence[multi_pos - 1].upper()
        multi_alts = {
            "A": ("C", "G"), "C": ("G", "T"),
            "G": ("T", "A"), "T": ("A", "C"),
        }[multi_ref]
        multi = work / "multi-alt-iupac.vcf"
        multi.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            f"17\t{multi_pos}\t.\t{multi_ref}\t{multi_alts[0]},{multi_alts[1]}\t50\tPASS\t.\tGT\t1/2\n"
        )
        if java_oracle:
            run([str(JAVA), "-jar", str(GATK_JAR), "IndexFeatureFile", "-I", str(multi)])
        multi_native = work / "alternate.multi-alt-iupac.native.fasta"
        run([str(ALTERNATE_BIN), "-R", str(REF), "-V", str(multi),
             "--use-iupac-sample", "SAMPLE", "-L", "17:1-1000", "-O", str(multi_native)])
        if java_oracle:
            multi_oracle = work / "alternate.multi-alt-iupac.oracle.fasta"
            run([str(JAVA), "-jar", str(GATK_JAR), "FastaAlternateReferenceMaker", "-R", str(REF),
                 "-V", str(multi), "--use-iupac-sample", "SAMPLE",
                 "-L", "17:1-1000", "-O", str(multi_oracle)])
            compare_outputs(multi_native, multi_oracle)

    print(json.dumps({"status": "pass", "java_oracle": java_oracle,
                      "tools": ["FastaReferenceMaker", "FastaAlternateReferenceMaker"]},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
