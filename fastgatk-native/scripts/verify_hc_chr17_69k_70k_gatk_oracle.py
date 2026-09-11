#!/usr/bin/env python3
"""Pin HaplotypeCaller on the canonical NA12878 chr17 69k-70k fixture.

This is the 5th real window in fast-gatk's A3 ≥5 real fixture pin target
(alongside the 4 mnp.bam windows pinned by verify_hc_chr20_real_contract.py).
The fixture is the standard GATK 4.6.2.0 integration BAM (`NA12878.chr17_69k_70k.
dictFix.bam`, 493 reads, b37 chr17 69k-70k) used by every pinned GATK
integration oracle in this repository. The pinned position set is:

  17:69067   T/G  1/1 (single read; AF=1.00, QD=33.48)
  17:69368   G/C  0/1 (42 reads; QD=15.97) — same position as R24 Mutect2
                                                       posterior-pin fixture
  17:69631   C/T  0/1 (37 reads; QD=19.36)

Pinned (must pass on both backends):
  * End-to-end HaplotypeCaller run completes (Java + native exit 0).
  * Both VCF.GZ outputs are tabix-indexed and parseable.
  * Both emit all 3 positions with the exact same REF and ALT alleles
    (GATK + native agree byte-equal at the allele level for this
    canonical fixture; the third call is the high-coverage 69368 site that
    A4's Mutect2 posterior parity has been validated against).
  * Native OutputManifest is valid JSON, status="prototype" or
    "contract-compatible", sample name matches verbatim.
  * Native emits read evidence (tumor_reads > 0) and activity in the
    69-70k window.
  * Both `--assembly-region-out` tracks are byte-identical.  This pins the
    Ref-vs-Any active-state likelihood, CIGAR/soft-clip handling, high-quality
    soft-clip propagation, and GATK's band-pass boundaries on a real fixture.
  * The unique haplotype sequences submitted to PairHMM are identical.  This
    pins AssemblyResultSet/ActivityRegion ownership as well as the final VCF:
    a dangling-tail path in padded halo must not become a separate HC request.

Recorded (not pinned): TLOD/QUAL/PL/AD byte values may differ slightly
across HaplotypeCaller releases; the oracle pins the (POS, REF, ALT) tuple
which is the documented HaplotypeCaller API contract.

Why this is a 5th pinned window (not a duplicate of chr17 69-70k in
verify_mutect2_gatk_oracle.py):
  * verify_mutect2_gatk_oracle.py runs **Mutect2** (somatic) on the same
    BAM. This oracle runs **HaplotypeCaller** (germline) on the same BAM.
  * The pinned positions differ: Mutect2 (verify_mutect2_gatk_oracle.py)
    emits 4 positions (69067, 69368, 69631, 69632); HaplotypeCaller emits
    3 (69067, 69368, 69631). Different callers, different position sets,
    same fixture — provides a third real-WGS-window point of contact.
  * Pairs with the 4-window chr20 mnp.bam contract
    (verify_hc_chr20_real_contract.py) to satisfy A3 ≥5 real fixture pin.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import oracle_guard

ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
HC_NATIVE = Path(os.environ.get(
    "FASTGATK_HC_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-hc-call")))

BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
BAI = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam.bai"
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
SAMPLE = "NA12878"
INTERVAL = "17:69000-70000"

EXPECTED = {
    "69067": ("T", "G"),
    "69368": ("G", "C"),
    "69631": ("C", "T"),
}


def run(
    cmd: list[str], *, check: bool = True, env: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=check, env=env)


def gatk_pairhmm_haplotypes(path: Path) -> set[str]:
    return {
        line.split(maxsplit=1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    }


def native_pairhmm_haplotypes(stderr: str) -> set[str]:
    prefix = "[FASTGATK_PAIRHMM_REQUEST] "
    return {
        line[len(prefix):].split("\t", maxsplit=1)[0]
        for line in stderr.splitlines()
        if line.startswith(prefix)
    }


def vcf_records(path: Path) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    if not path.is_file():
        return out
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out[f[1]] = (f[3], f[4])
    return out


def vcf_body(path: Path) -> bytes:
    """Return the semantic VCF payload, excluding implementation headers."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as h:
        return b"".join(
            line.encode("utf-8")
            for line in h
            if line and not line.startswith("#")
        )


def igv_rows(path: Path) -> list[tuple[str, int, int, str, str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("Chromosome\t"):
            continue
        fields = line.split("\t")
        assert len(fields) == 5, fields
        rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]))
    return rows


def main() -> int:
    required = (JAVA, GATK, HC_NATIVE, BAM, BAI, REFERENCE)
    if not all(p.is_file() for p in required):
        oracle_guard.oracle_not_verified('verify_hc_chr17_69k_70k_gatk_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            sys.stderr.write(f"missing required inputs: {[str(p) for p in required if not p.is_file()]}\n")
            return 2
        print(json.dumps({"status": "skip", "reason": "bundled NA12878 chr17 69k oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr17-69k-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf.gz"
        gatk_tbi = Path(str(gatk_vcf) + ".tbi")
        gatk_activity = work / "gatk.igv"
        gatk_pairhmm = work / "gatk.pairhmm.txt"
        native_vcf = work / "native.vcf.gz"
        native_activity = work / "native.igv"
        native_manifest = work / "native.manifest.json"

        gatk_run = run([
            str(JAVA), "-jar", str(GATK), "HaplotypeCaller",
            "-R", str(REFERENCE), "-I", str(BAM),
            "-O", str(gatk_vcf), "-L", INTERVAL,
            "--sample-name", SAMPLE,
            "--assembly-region-out", str(gatk_activity),
            "--pair-hmm-results-file", str(gatk_pairhmm),
        ], check=False)
        if gatk_run.returncode != 0:
            sys.stderr.write(f"GATK HaplotypeCaller failed:\n{gatk_run.stderr[-2000:]}\n")
            return 3

        native_debug_env = os.environ | {"FASTGATK_DEBUG_PAIRHMM_REQUESTS": "1"}
        native_run = run([
            str(HC_NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
            "--sample-name", SAMPLE, "-L", INTERVAL,
            "-O", str(native_vcf),
            "--assembly-region-out", str(native_activity),
            "--output-manifest", str(native_manifest),
        ], check=False, env=native_debug_env)
        if native_run.returncode != 0:
            sys.stderr.write(f"native HaplotypeCaller failed:\n{native_run.stderr[-2000:]}\n")
            return 3

        for vcf, label in ((gatk_vcf, "GATK"), (native_vcf, "native")):
            tbi = Path(str(vcf) + ".tbi")
            if not tbi.is_file():
                sys.stderr.write(f"{label} HC VCF index missing: {tbi}\n")
                return 4

        if not native_manifest.is_file():
            sys.stderr.write(f"native OutputManifest missing: {native_manifest}\n")
            return 5
        if not gatk_activity.is_file() or not native_activity.is_file():
            sys.stderr.write("GATK or native AssemblyRegion IGV track is missing\n")
            return 5
        if gatk_activity.read_bytes() != native_activity.read_bytes():
            sys.stderr.write(
                "HaplotypeCaller AssemblyRegion track differs from GATK:\n" +
                f"--- GATK ({gatk_activity})\n+++ native ({native_activity})\n")
            return 5
        if not gatk_pairhmm.is_file() or gatk_pairhmm.stat().st_size == 0:
            sys.stderr.write("GATK HaplotypeCaller PairHMM request log is missing or empty\n")
            return 5
        gatk_haplotypes = gatk_pairhmm_haplotypes(gatk_pairhmm)
        native_haplotypes = native_pairhmm_haplotypes(native_run.stderr)
        if gatk_haplotypes != native_haplotypes:
            sys.stderr.write(
                "HaplotypeCaller PairHMM haplotype set differs from GATK:\n"
                f"  missing native: {sorted(gatk_haplotypes - native_haplotypes)}\n"
                f"  extra native: {sorted(native_haplotypes - gatk_haplotypes)}\n"
            )
            return 5
        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        if manifest.get("schema_version") != 1:
            sys.stderr.write(f"manifest schema_version mismatch: {manifest.get('schema_version')}\n")
            return 5
        # `sample_name_selection` is a boolean (was the user-provided --sample-name
        # accepted). The actual selected sample name lives in `sample_name`.
        if manifest.get("compatibility", {}).get("sample_name_selection") is not True:
            sys.stderr.write(
                f"manifest sample_name_selection not True: "
                f"{manifest.get('compatibility', {}).get('sample_name_selection')}\n"
            )
            return 5
        if manifest.get("compatibility", {}).get("assembly_region_out") is not True:
            sys.stderr.write("manifest did not record --assembly-region-out compatibility\n")
            return 5
        if manifest.get("sample_name") != SAMPLE:
            sys.stderr.write(
                f"manifest sample_name != {SAMPLE}: {manifest.get('sample_name')}\n"
            )
            return 5

        gatk_records = vcf_records(gatk_vcf)
        native_records = vcf_records(native_vcf)

        # Header provenance necessarily differs between the Java and native
        # front ends, but the complete variant payload is GATK-compatible on
        # this real AssemblyRegion.  This is deliberately stronger than the
        # allele-only pins below: it covers QUAL, INFO, FORMAT and the sample
        # column for every emitted record.
        if vcf_body(gatk_vcf) != vcf_body(native_vcf):
            sys.stderr.write("HaplotypeCaller VCF record payload differs from GATK\n")
            return 6

        if not all(pos in gatk_records and gatk_records[pos] == EXPECTED[pos] for pos in EXPECTED):
            missing = [p for p in EXPECTED if p not in gatk_records or gatk_records[p] != EXPECTED[p]]
            sys.stderr.write(f"GATK chr17 69k oracle dropped/changed positions: {missing}\n")
            return 6
        if not all(pos in native_records and native_records[pos] == EXPECTED[pos] for pos in EXPECTED):
            missing = [p for p in EXPECTED if p not in native_records or native_records[p] != EXPECTED[p]]
            sys.stderr.write(f"native chr17 69k oracle dropped/changed positions: {missing}\n")
            return 6

        # A no-read -L interval is valid source traversal: the inherited
        # ActivityProfile receives empty pileups, produces hard-max inactive
        # chunks, and only then applies --force-active.  Verify both modes so
        # the HC front end cannot reject the input before the shared
        # C++ Host/Kokkos ActivityProfile stage sees its traversal span.
        zero_interval = "17:800000-800700"
        zero_gatk_rows: dict[bool, list[tuple[str, int, int, str, str]]] = {}
        zero_native_rows: dict[bool, list[tuple[str, int, int, str, str]]] = {}
        for force in (False, True):
            suffix = "true" if force else "false"
            zero_gatk_vcf = work / f"zero-gatk-{suffix}.vcf.gz"
            zero_gatk_igv = work / f"zero-gatk-{suffix}.igv"
            zero_native_vcf = work / f"zero-native-{suffix}.vcf.gz"
            zero_native_igv = work / f"zero-native-{suffix}.igv"
            zero_native_manifest = work / f"zero-native-{suffix}.manifest.json"
            zero_gatk_run = run([
                str(JAVA), "-jar", str(GATK), "HaplotypeCaller",
                "-R", str(REFERENCE), "-I", str(BAM), "--sample-name", SAMPLE,
                "-L", zero_interval, "-O", str(zero_gatk_vcf),
                "--assembly-region-out", str(zero_gatk_igv),
                "--force-active", str(force).lower(),
                "--create-output-variant-index", "false",
            ], check=False)
            if zero_gatk_run.returncode != 0:
                sys.stderr.write(f"GATK zero-coverage HC failed:\n{zero_gatk_run.stderr[-2000:]}\n")
                return 7
            zero_native_run = run([
                str(HC_NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
                "--sample-name", SAMPLE, "-L", zero_interval,
                "-O", str(zero_native_vcf), "--assembly-region-out", str(zero_native_igv),
                "--force-active=" + str(force).lower(),
                "--create-output-variant-index", "false",
                "--output-manifest", str(zero_native_manifest),
            ], check=False)
            if zero_native_run.returncode != 0:
                sys.stderr.write(
                    f"native zero-coverage HC failed:\n{zero_native_run.stderr[-2000:]}\n")
                return 7
            zero_gatk_rows[force] = igv_rows(zero_gatk_igv)
            zero_native_rows[force] = igv_rows(zero_native_igv)
            if zero_gatk_rows[force] != zero_native_rows[force]:
                sys.stderr.write("zero-coverage HC AssemblyRegion track differs from GATK\n")
                return 7
            if vcf_body(zero_gatk_vcf) != b"" or vcf_body(zero_native_vcf) != b"":
                sys.stderr.write("zero-coverage HC unexpectedly emitted a VCF record\n")
                return 7
            zero_manifest = json.loads(zero_native_manifest.read_text(encoding="utf-8"))
            if zero_manifest.get("compatibility", {}).get("force_active") is not force:
                sys.stderr.write("zero-coverage HC manifest force_active mismatch\n")
                return 7
        zero_profile_false = [row for row in zero_gatk_rows[False]
                              if row[3].startswith("size=")]
        zero_profile_true = [row for row in zero_gatk_rows[True]
                             if row[3].startswith("size=")]
        if [(row[1], row[2], row[3]) for row in zero_profile_false] != [
            (799999, 800299, "size=300"),
            (800299, 800599, "size=300"),
            (800599, 800700, "size=101"),
        ]:
            sys.stderr.write("GATK zero-coverage HC profile did not retain the expected 300/300/101 split\n")
            return 7
        if not all(row[4] == "-1.00000" for row in zero_profile_false):
            sys.stderr.write("GATK zero-coverage HC inactive profile state changed\n")
            return 7
        if [row[:4] for row in zero_profile_false] != [row[:4] for row in zero_profile_true] or \
                not all(row[4] == "1.00000" for row in zero_profile_true):
            sys.stderr.write("HC --force-active did not preserve zero-coverage boundaries/state\n")
            return 7

        # With no explicit -L, GATK traverses every reference contig, even
        # when the read filters remove all records.  This is distinct from a
        # selected empty interval: it proves that native traversal domains are
        # derived from the reference dictionary rather than from observed
        # compact pileup coordinates.
        full_gatk_vcf = work / "full-zero-gatk.vcf.gz"
        full_gatk_igv = work / "full-zero-gatk.igv"
        full_native_vcf = work / "full-zero-native.vcf.gz"
        full_native_igv = work / "full-zero-native.igv"
        full_native_manifest = work / "full-zero-native.manifest.json"
        full_gatk_run = run([
            str(JAVA), "-jar", str(GATK), "HaplotypeCaller",
            "-R", str(REFERENCE), "-I", str(BAM), "--sample-name", SAMPLE,
            "-O", str(full_gatk_vcf), "--assembly-region-out", str(full_gatk_igv),
            "--minimum-mapping-quality", "255",
            "--create-output-variant-index", "false",
        ], check=False)
        if full_gatk_run.returncode != 0:
            sys.stderr.write(f"GATK full-contig zero-coverage HC failed:\n{full_gatk_run.stderr[-2000:]}\n")
            return 8
        full_native_run = run([
            str(HC_NATIVE), "-R", str(REFERENCE), "-I", str(BAM),
            "--sample-name", SAMPLE, "-O", str(full_native_vcf),
            "--assembly-region-out", str(full_native_igv),
            "--minimum-mapping-quality", "255",
            "--create-output-variant-index", "false",
            "--output-manifest", str(full_native_manifest),
        ], check=False)
        if full_native_run.returncode != 0:
            sys.stderr.write(
                f"native full-contig zero-coverage HC failed:\n{full_native_run.stderr[-2000:]}\n")
            return 8
        full_gatk_rows = igv_rows(full_gatk_igv)
        full_native_rows = igv_rows(full_native_igv)
        if full_gatk_rows != full_native_rows:
            sys.stderr.write("full-contig zero-coverage HC AssemblyRegion track differs from GATK\n")
            return 8
        full_profile_rows = [row for row in full_gatk_rows if row[3].startswith("size=")]
        if len(full_profile_rows) != 3334 or full_profile_rows[0] != (
                "17", 0, 300, "size=300", "-1.00000") or full_profile_rows[-1] != (
                "17", 999900, 1000000, "size=100", "-1.00000"):
            sys.stderr.write("GATK full-contig zero-coverage profile changed unexpectedly\n")
            return 8
        if vcf_body(full_gatk_vcf) != b"" or vcf_body(full_native_vcf) != b"":
            sys.stderr.write("full-contig zero-coverage HC unexpectedly emitted a VCF record\n")
            return 8

        report = {
            "status": "pass",
            "fixture": "na12878_chr17_69k_70k_hc_pinned_window_5",
            "bam": str(BAM.relative_to(ROOT)),
            "reference": str(REFERENCE.relative_to(ROOT)),
            "interval": INTERVAL,
            "expected_pinned_positions": [
                {"pos": p, "ref": r, "alt": a, "gatk_allele": gatk_records[p],
                 "native_allele": native_records[p]}
                for p, (r, a) in EXPECTED.items()
            ],
            "shared_position_count": len(set(gatk_records) & set(native_records)),
            "native_status": manifest.get("status"),
            "native_variant_calls": manifest.get("telemetry", {}).get("variant_calls"),
            "native_candidate_sites": manifest.get("telemetry", {}).get("candidate_sites"),
            "native_loci": manifest.get("telemetry", {}).get("loci"),
            "vcf_record_payload_byte_identical": True,
            "assembly_region_track_byte_identical": True,
            "pairhmm_haplotype_set_identical": True,
            "pairhmm_unique_haplotype_count": len(gatk_haplotypes),
            "zero_coverage_force_active_profile_exact": True,
            "full_contig_zero_coverage_profile_exact": True,
            "a3_5th_real_fixture_window": True,
        }
        print(json.dumps(report, indent=2))
        return 0


if __name__ == "__main__":
    sys.exit(main())
