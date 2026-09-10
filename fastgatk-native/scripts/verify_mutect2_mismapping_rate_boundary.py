#!/usr/bin/env python3
"""Pin Mutect2's global-mismapping-rate boundary to GATK 4.6.2.0.

GATK accepts positive Q values (the per-read cap) and negative values (the
explicitly disabled cap).  A zero value is not a third "uncapped" mode: it
becomes a zero log-likelihood floor and GATK rejects it in
AlleleLikelihoods.normalizeLikelihoods.  Native must fail closed as well;
silently flattening every read would change the downstream TLOD model.  The
small direct Java check also guards the less visible disabled-cap behavior:
GATK leaves a valid ``-Infinity`` likelihood cell unchanged when Q is -1.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def verify_gatk_disabled_cap(java: Path, gatk: Path) -> None:
    """Run GATK's AlleleLikelihoods API on a zero-probability cell.

    Mutect2 passes Double.NEGATIVE_INFINITY to normalizeLikelihoods for a
    negative CLI Q.  Calling the actual 4.6.2.0 class prevents this boundary
    from being reduced to a Python reimplementation of the Java contract.
    """
    source = r'''
import java.util.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import htsjdk.variant.variantcontext.*;
var alleles = new IndexedAlleleList<Allele>(List.of(Allele.REF_A, Allele.ALT_C));
var samples = new IndexedSampleList("S");
var evidence = List.of(List.of(new SimpleInterval("1",1,1), new SimpleInterval("1",2,2), new SimpleInterval("1",3,3)));
var filtered = List.of(List.<SimpleInterval>of());
var values = new double[][][] {{{0.0, Double.NEGATIVE_INFINITY, -3.0}, {-1.0,-1.0,-1.0}}};
var likelihoods = AlleleLikelihoods.createAlleleLikelihoods(alleles, samples, evidence, filtered, values);
likelihoods.normalizeLikelihoods(Double.NEGATIVE_INFINITY, true);
System.out.println(Arrays.toString(new double[]{likelihoods.sampleMatrix(0).get(0,0), likelihoods.sampleMatrix(0).get(0,1), likelihoods.sampleMatrix(0).get(0,2)}));
'''
    prefs_dir = tempfile.mkdtemp(prefix="fastgatk-jshell-prefs-")
    check = subprocess.run(
        [str(java).replace("/java", "/jshell"),
         f"-J-Djava.util.prefs.userRoot={prefs_dir}",
         "--class-path", str(gatk)],
        input=source, text=True, capture_output=True, check=False)
    if check.returncode != 0:
        raise AssertionError(f"GATK disabled-cap probe failed: {check.stderr[-2000:]}")
    if "[0.0, -Infinity, -3.0]" not in check.stdout:
        raise AssertionError(
            "GATK disabled-cap probe changed a -Infinity likelihood: "
            f"{check.stdout[-2000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 boundary inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-mismapping-rate-") as directory:
        work = Path(directory)
        gatk_run = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(work / "gatk.vcf.gz"),
            "--phred-scaled-global-read-mismapping-rate", "0",
        ])
        native_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(work / "native.vcf.gz"),
            "--phred-scaled-global-read-mismapping-rate", "0",
            "--min-depth", "1", "--min-alt-support", "1",
        ])
        if gatk_run.returncode == 0:
            raise AssertionError("GATK unexpectedly accepted Q=0")
        if native_run.returncode == 0:
            raise AssertionError("native unexpectedly accepted Q=0")
        gatk_text = f"{gatk_run.stdout}\n{gatk_run.stderr}"
        native_text = f"{native_run.stdout}\n{native_run.stderr}"
        if "minimum reference likelihood fall must be negative" not in gatk_text:
            raise AssertionError(f"unexpected GATK Q=0 failure: {gatk_text[-2000:]}")
        if "calling thresholds must be positive" not in native_text:
            raise AssertionError(f"unexpected native Q=0 failure: {native_text[-2000:]}")
        verify_gatk_disabled_cap(java, gatk)

        # Negative Q is the real GATK disabled-cap path.  Both implementations
        # must accept it on the pinned fixture; exact TLOD equality is covered
        # by the broad oracle and is intentionally not claimed here because
        # assembly/PairHMM inputs remain release-specific.
        gatk_uncapped = run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference),
            "-I", str(bam), "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(work / "gatk-uncapped.vcf.gz"),
            "--phred-scaled-global-read-mismapping-rate", "-1",
        ])
        native_uncapped = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(work / "native-uncapped.vcf.gz"),
            "--phred-scaled-global-read-mismapping-rate", "-1",
            "--min-depth", "1", "--min-alt-support", "1",
        ])
        if gatk_uncapped.returncode != 0:
            raise AssertionError(f"GATK unexpectedly rejected Q=-1: {gatk_uncapped.stderr[-2000:]}")
        if native_uncapped.returncode != 0:
            raise AssertionError(f"native unexpectedly rejected Q=-1: {native_uncapped.stderr[-2000:]}")
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "phred_scaled_global_read_mismapping_rate": 0,
            "gatk_rejected": True,
            "native_rejected": True,
            "negative_cap": "gatk_java_and_native_accepted",
            "gatk_disabled_inf_probe": "[0.0, -Infinity, -3.0]",
            "semantic": "zero-invalid-negative-disables-positive-caps-and-preserves-inf",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
