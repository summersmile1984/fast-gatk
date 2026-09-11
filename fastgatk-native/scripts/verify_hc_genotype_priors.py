#!/usr/bin/env python3
"""Verify GATK assumingHW priors are consumed for a multi-ALT HC record."""
from __future__ import annotations

import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def record(path: Path) -> tuple[list[str], dict[str, str]]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "chr1" and fields[1] == "11":
            values = dict(zip(fields[8].split(":"), fields[9].split(":")))
            return fields, values
    raise AssertionError("missing synthetic multi-ALT sentinel")


def java_assuming_hw_oracle(root: Path, work: Path) -> dict[str, list[float]]:
    """Read GenotypePriorCalculator.assumingHW directly from GATK 4.6.2.0."""
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    if not jar.exists() or not java.exists():
        oracle_guard.oracle_not_verified('verify_hc_genotype_priors.py', jar, java)
        raise AssertionError("GATK/JDK is required for the genotype-prior oracle")
    source = work / "HCPriorOracle.java"
    source.write_text(
        """import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.GenotypePriorCalculator;
import java.util.Arrays;
public class HCPriorOracle {
  public static void main(String[] args) {
    final var calc = GenotypePriorCalculator.assumingHW(Math.log10(1.0e-3), Math.log10(1.0e-4));
    final var snps = Arrays.asList(Allele.create(\"A\", true), Allele.create(\"C\"), Allele.create(\"G\"));
    final var indel = Arrays.asList(Allele.create(\"A\", true), Allele.create(\"AT\"), Allele.create(\"A\"));
    for (int ploidy : new int[]{2,3}) {
      System.out.println(\"SNP\" + ploidy + \"=\" + Arrays.toString(calc.getLog10Priors(ploidy, snps)));
      System.out.println(\"INDEL\" + ploidy + \"=\" + Arrays.toString(calc.getLog10Priors(ploidy, indel)));
    }
  }
}
""",
        encoding="utf-8")
    classes = work / "java-classes"
    classes.mkdir()
    subprocess.run([str(java.with_name("javac")), "-cp", str(jar), "-d", str(classes), str(source)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(java), "-cp", f"{classes}:{jar}", "HCPriorOracle"],
                            check=True, text=True, capture_output=True)
    vectors: dict[str, list[float]] = {}
    for line in result.stdout.splitlines():
        name, values = line.split("=", 1)
        vectors[name] = [float(value.strip()) for value in values.strip("[]").split(",")]
    return vectors


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    if not binary.exists():
        raise SystemExit("missing HC binary")

    sam_header = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:40\n"
        "@RG\tID:rg1\tSM:genotype-prior\n"
    )
    # Use a non-repetitive sequence so ordinary default HC produces the
    # candidate via its retained assembly EventMap.
    reference_sequence = "CGTATCGGCTAGCTTACGATCGGTACCTGATCAGTCCGTA"
    reference_text = ">chr1\n" + reference_sequence + "\n"
    reads = []
    for index, alternate in enumerate(("C", "G")):
        for copy in range(4):
            sequence = reference_sequence[:10] + alternate + reference_sequence[11:21]
            reads.append(
                f"{alternate.lower()}{index}{copy}\t0\tchr1\t1\t60\t21M\t*\t0\t0\t"
                f"{sequence}\t{'5' * 21}\tRG:Z:rg1\n"
            )

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-genotype-priors-") as directory:
        work = Path(directory)
        java_vectors = java_assuming_hw_oracle(root, work)
        sam = work / "multi-alt.sam"
        reference = work / "reference.fa"
        sam.write_text(sam_header + "".join(reads), encoding="utf-8")
        reference.write_text(reference_text, encoding="utf-8")

        def run(label: str, extra: list[str], gvcf: bool = False, ploidy: int = 2) -> tuple[list[str], dict[str, str], dict[str, object]]:
            output = work / f"{label}.vcf"
            manifest = work / f"{label}.manifest.json"
            subprocess.run([
                str(binary), "-I", str(sam), "-R", str(reference), "-L", "chr1:1-40",
                "-O", str(output), "--threads", "2", "--min-depth", "1",
                "--min-alt-support", "1", "--sample-ploidy", str(ploidy),
                *( ["--emit-ref-confidence", "GVCF"] if gvcf else []),
                *extra, "--output-manifest", str(manifest),
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            fields, values = record(output)
            return fields, values, json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]

        prior_fields, prior, prior_telemetry = run(
            "prior", ["--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES"])
        default_fields, default_values, default_telemetry = run("default", [])
        flat_fields, flat, flat_telemetry = run("flat", ["--no-genotype-priors"])
        legacy_fields, legacy, legacy_telemetry = run(
            "legacy-q45", ["--no-genotype-priors", "--pcr-indel-model", "NONE"])
        gvcf_prior_fields, gvcf_prior, gvcf_prior_telemetry = run(
            "prior-gvcf", ["--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES"], True)
        gvcf_default_fields, gvcf_default, gvcf_default_telemetry = run(
            "default-gvcf", [], True)
        gvcf_flat_fields, gvcf_flat, gvcf_flat_telemetry = run(
            "flat-gvcf", ["--no-genotype-priors"], True)
        gvcf_poly_prior_fields, gvcf_poly_prior, gvcf_poly_prior_telemetry = run(
            "prior-gvcf-poly", ["--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES"], True, 3)
        gvcf_poly_flat_fields, gvcf_poly_flat, gvcf_poly_flat_telemetry = run(
            "flat-gvcf-poly", ["--no-genotype-priors"], True, 3)

    assert prior_fields[4] == "C,G" and flat_fields[4] == "C,G"
    # PairHMM concrete-ALT marginalization uses one (locus,read) row.  The
    # former biallelic row contract could treat the sibling ALT as REF and
    # produced a spurious 0/1 call; the concrete matrix correctly preserves
    # the 1/2 likelihood and the cross-ALT PL entries.
    assert prior["PL"] == flat["PL"] == "174,87,75,87,0,75"
    assert legacy["PL"] == "174,87,75,87,0,75"
    assert prior["GT"] == "1/2" and prior["GQ"] == "49"
    prior_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                  for item in prior_fields[7].split(";") if "=" in item}
    default_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                    for item in default_fields[7].split(";") if "=" in item}
    flat_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                 for item in flat_fields[7].split(";") if "=" in item}
    assert prior_info["AC"] == "1,1"
    assert default_values["GT"] == flat["GT"] == "1/2"
    assert default_values["GQ"] == flat["GQ"] == "75"
    assert default_info["AC"] == flat_info["AC"] == "1,1"
    assert default_telemetry["genotype_priors"] is True
    assert default_telemetry["genotype_assignment_method"] == "USE_PLS_TO_ASSIGN"
    assert default_telemetry["joint_genotype_priors_used"] is False
    assert flat["GT"] == "1/2" and flat["GQ"] == "75"
    assert flat_info["AC"] == "1,1"
    assert prior_telemetry["genotype_priors"] is True
    assert prior_telemetry["genotype_assignment_method"] == "USE_POSTERIOR_PROBABILITIES"
    assert prior_telemetry["joint_genotype_priors_used"] is True
    assert prior_telemetry["genotype_prior_kernel_calls"] > 0
    assert prior_telemetry["genotype_prior_execution_space"]
    assert flat_telemetry["genotype_priors"] is False
    assert flat_telemetry["genotype_assignment_method"] == "USE_PLS_TO_ASSIGN"
    assert flat_telemetry["joint_genotype_priors_used"] is False
    assert flat_telemetry["genotype_prior_kernel_calls"] == 0
    assert gvcf_prior_fields[4] == "C,G,<NON_REF>"
    assert gvcf_flat_fields[4] == "C,G,<NON_REF>"
    # The symbolic row is retained in the emitted matrix but is not a
    # candidate-site GT.  Prior-aware assignment and flat PLS assignment use
    # the same concrete PL prefix but differ in GQ as in GATK's assignment
    # boundary.
    assert gvcf_prior["GT"] == "1/2" and gvcf_flat["GT"] == "1/2"
    assert gvcf_prior["PL"] == gvcf_flat["PL"]
    assert gvcf_prior["PL"] == "174,87,75,87,0,75,174,87,87,174"
    assert gvcf_prior["GQ"] == "49" and gvcf_flat["GQ"] == "75"
    assert gvcf_default["GT"] == gvcf_flat["GT"] == "1/2"
    assert gvcf_default["GQ"] == gvcf_flat["GQ"] == "75"
    assert gvcf_default_telemetry["genotype_priors"] is True
    assert gvcf_default_telemetry["genotype_assignment_method"] == "USE_PLS_TO_ASSIGN"
    assert gvcf_default_telemetry["joint_genotype_priors_used"] is False
    assert gvcf_prior_telemetry["genotype_priors"] is True
    assert gvcf_prior_telemetry["genotype_assignment_method"] == "USE_POSTERIOR_PROBABILITIES"
    assert gvcf_prior_telemetry["joint_genotype_priors_used"] is True
    assert gvcf_prior_telemetry["genotype_prior_kernel_calls"] > 0
    assert gvcf_prior_telemetry["genotype_prior_execution_space"]
    assert gvcf_flat_telemetry["genotype_priors"] is False
    assert gvcf_flat_telemetry["genotype_assignment_method"] == "USE_PLS_TO_ASSIGN"
    assert gvcf_flat_telemetry["joint_genotype_priors_used"] is False
    assert gvcf_flat_telemetry["genotype_prior_kernel_calls"] == 0
    assert gvcf_poly_prior_fields[4] == "C,G,<NON_REF>"
    assert gvcf_poly_flat_fields[4] == "C,G,<NON_REF>"
    assert gvcf_poly_prior["GT"] == "0/1/2" and gvcf_poly_flat["GT"] == "1/1/2"
    assert gvcf_poly_prior["PL"] == gvcf_poly_flat["PL"]
    assert gvcf_poly_prior["PL"] == "172,92,80,73,92,12,0,80,0,73,172,92,80,92,12,80,172,92,92,172"
    assert gvcf_poly_prior["GQ"] == "15" and gvcf_poly_flat["GQ"] == "0"
    assert gvcf_poly_prior_telemetry["genotype_priors"] is True
    assert gvcf_poly_prior_telemetry["genotype_assignment_method"] == "USE_POSTERIOR_PROBABILITIES"
    assert gvcf_poly_prior_telemetry["joint_genotype_priors_used"] is True
    assert gvcf_poly_prior_telemetry["genotype_prior_kernel_calls"] > 0
    assert gvcf_poly_prior_telemetry["genotype_prior_execution_space"]
    assert gvcf_poly_flat_telemetry["genotype_priors"] is False
    assert gvcf_poly_flat_telemetry["joint_genotype_priors_used"] is False
    assert gvcf_poly_flat_telemetry["genotype_prior_kernel_calls"] == 0
    assert java_vectors["SNP2"][1] == -3.4771212547196626
    assert java_vectors["SNP2"][2] == -6.477121254719663
    assert prior_telemetry["genotype_prior_snp_normalization"] == "log10(3)"
    native_prior_candidates = {
        item["alt"]: item for item in prior_telemetry["genotype_prior_candidates"]
    }
    assert set(native_prior_candidates) >= {"C", "G"}
    for alt in ("C", "G"):
        native = native_prior_candidates[alt]
        assert math.isclose(native["het"], java_vectors["SNP2"][1], rel_tol=0.0, abs_tol=1e-12)
        expected_hom_alt = java_vectors["SNP2"][2 if alt == "C" else 5]
        assert math.isclose(native["hom_alt"], expected_hom_alt, rel_tol=0.0, abs_tol=1e-12)
    print(json.dumps({
        "status": "pass",
        "tool": "HaplotypeCaller",
        "multi_alt": "C,G",
        "prior_gt": prior["GT"],
        "flat_gt": flat["GT"],
        "gvcf_prior_gq": gvcf_prior["GQ"],
        "gvcf_flat_gq": gvcf_flat["GQ"],
        "gvcf_poly_prior_gq": gvcf_poly_prior["GQ"],
        "gvcf_poly_flat_gq": gvcf_poly_flat["GQ"],
        "pl": prior["PL"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
