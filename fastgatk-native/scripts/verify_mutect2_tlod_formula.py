#!/usr/bin/env python3
"""Independent TLOD engine-parity verifier (A4, no GATK needed).

Runs native Mutect2 on the pinned chr17:69000-70000 fixture with the
FASTGATK_DEBUG_TLOD fragment-matrix dump enabled, then recomputes each
single-ALT candidate's TLOD from the dumped candidate x fragment ref/alt
log-likelihood matrix using a faithful Python transcription of GATK 4.6.2.0
SomaticLikelihoodsEngine (Dirichlet posterior fixed point over fragments,
logEvidence with natural-log domain, TLOD = log10(e) * (logEvidence(all) -
logEvidence(no-ALT))).  The native kernel must reproduce the transcribed Java
formula on its own inputs to ~1e-4 log10 at every emitted single-ALT locus.

This isolates the posterior/evidence engine from the upstream assembly and
read/haplotype preparation parity that is tracked separately: any engine
change that shifts the Dirichlet fixed point, xLogx/negligible thresholds,
natural-log scaling or digamma use would fail here regardless of upstream
input differences.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

LN10 = math.log(10.0)
K_MISSING = -1.0e299


def digamma(x: float) -> float:
    """Commons-math3 digamma (recursive below 8, asymptotic above)."""
    value = 0.0
    while x < 8.0:
        value -= 1.0 / x
        x += 1.0
    inverse_square = 1.0 / (x * x)
    return (value + math.log(x) - 0.5 / x
            - inverse_square * (1.0 / 12.0
                                - inverse_square * (1.0 / 120.0
                                                    - inverse_square / 252.0)))


def softmax(logits):
    high = max(logits)
    exp_values = [math.exp(value - high) for value in logits]
    total = sum(exp_values)
    return [value / total for value in exp_values]


def log_dirichlet_normalization(params):
    return (math.lgamma(sum(params))
            - sum(math.lgamma(param) for param in params))


def log_evidence_biallelic(ref_log10, alt_log10):
    """Java SomaticLikelihoodsEngine.logEvidence for a 2-allele matrix.

    ref_log10/alt_log10 are per-fragment log10 likelihoods (None when the
    engine has no likelihood for that allele on the fragment).  Returns the
    natural-log evidence.
    """
    fragments = len(ref_log10)
    prior = (1.0, 1.0)

    # alleleFractionsPosterior fixed point: alpha <- counts + prior.
    alpha = [1.0, 1.0]
    for _ in range(1000):
        alpha_sum = alpha[0] + alpha[1]
        log_weight_ref = digamma(alpha[0]) - digamma(alpha_sum)
        log_weight_alt = digamma(alpha[1]) - digamma(alpha_sum)
        counts = [0.0, 0.0]
        for index in range(fragments):
            if alt_log10[index] is None:
                counts[0] += 1.0
                continue
            logits = [log_weight_ref + ref_log10[index] * LN10,
                      log_weight_alt + alt_log10[index] * LN10]
            responsibilities = softmax(logits)
            counts[0] += responsibilities[0]
            counts[1] += responsibilities[1]
        next_alpha = [prior[0] + counts[0], prior[1] + counts[1]]
        relative_l1 = (abs(next_alpha[0] - alpha[0])
                       + abs(next_alpha[1] - alpha[1])) / sum(next_alpha)
        alpha = next_alpha
        if relative_l1 < 0.001:
            break

    alpha_sum = alpha[0] + alpha[1]
    posterior_normalization = log_dirichlet_normalization(alpha)
    log_weight_ref = digamma(alpha[0]) - digamma(alpha_sum)
    log_weight_alt = digamma(alpha[1]) - digamma(alpha_sum)

    likelihood_entropy = 0.0
    for index in range(fragments):
        if alt_log10[index] is None:
            likelihood_entropy += ref_log10[index] * LN10
            continue
        logits = [log_weight_ref + ref_log10[index] * LN10,
                  log_weight_alt + alt_log10[index] * LN10]
        responsibilities = softmax(logits)
        entropy = sum(x * math.log(x) for x in responsibilities
                      if x >= 1e-8)
        ref_contribution = 0.0 if responsibilities[0] < 1e-10 \
            else ref_log10[index] * LN10 * responsibilities[0]
        alt_contribution = 0.0 if responsibilities[1] < 1e-10 \
            else alt_log10[index] * LN10 * responsibilities[1]
        likelihood_entropy += ref_contribution + alt_contribution - entropy
    return 0.0 - posterior_normalization + likelihood_entropy


def tlod_biallelic(ref_log10, alt_log10):
    """log10(e) * (logEvidence(all) - logEvidence(no-ALT))."""
    fragments = len(ref_log10)
    no_alt_evidence = sum(ref_log10[index] * LN10
                          for index in range(fragments))
    all_evidence = log_evidence_biallelic(ref_log10, alt_log10)
    return (all_evidence - no_alt_evidence) / LN10


def run(command, env=None):
    return subprocess.run(command, text=True, capture_output=True,
                          check=False, env=env)


def parse_dump(path):
    """Candidate position/ref/alt -> [ref_log10 list, alt_log10 list|None]."""
    candidates = {}
    header = None
    refs = []
    alts = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("[FASTGATK_DEBUG_TLOD:") and "candidate" in line:
            if header is not None:
                key = (header[0], header[1], header[2])
                if key not in candidates:
                    candidates[key] = (refs, alts)
            fields = {}
            for token in line.split():
                if "=" in token:
                    name, _, value = token.partition("=")
                    fields[name] = value
            header = (int(fields["pos"]), fields["ref"], fields["alt"])
            refs = []
            alts = []
            continue
        match = re.match(r"\s*g=(\d+) ref=(-?[\d.eE+]+) alt=(-?[\d.eE+]+)",
                         line)
        if match is not None and header is not None:
            ref = float(match.group(2))
            alt = float(match.group(3))
            refs.append(ref if ref > K_MISSING / 2 else None)
            alts.append(alt if alt > K_MISSING / 2 else None)
    if header is not None:
        key = (header[0], header[1], header[2])
        if key not in candidates:
            candidates[key] = (refs, alts)
    return candidates


def read_vcf(path):
    """VCF position (1-based)/ref/alt -> TLOD float."""
    records = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            position = int(fields[1])
            ref, alt = fields[3], fields[4]
            tlod = None
            for token in fields[7].split(";"):
                if token.startswith("TLOD="):
                    tlod = float(token.split("=", 1)[1])
            if tlod is not None:
                records[(position, ref, alt)] = tlod
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")
    ))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, bam, reference, Path(f"{bam}.bai"), Path(f"{reference}.fai"))
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_MUTECT2_FIXTURE") == "1":
            raise SystemExit("bundled Mutect2 fixture inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled Mutect2 fixture unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-tlod-formula-") as directory:
        work = Path(directory)
        native_output = work / "native.vcf.gz"
        dump = work / "native.tlod.txt"
        native_env = dict(os.environ)
        native_env["FASTGATK_DEBUG_TLOD"] = "1"
        native_run = run([
            str(binary), "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000",
            "-O", str(native_output), "--min-depth", "1",
            "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ], native_env)
        if native_run.returncode != 0:
            raise AssertionError(
                f"native Mutect2 run failed: {native_run.stderr[-2000:]}")
        dump.write_text(native_run.stderr, encoding="utf-8")
        full_tlod = {}
        for line in native_run.stderr.splitlines():
            if not line.startswith("[FASTGATK_TLOD_FULL] "):
                continue
            fields = {}
            for token in line.split():
                if "=" in token:
                    name, _, value = token.partition("=")
                    fields[name] = value
            try:
                full_tlod[(int(fields["pos"]), fields["ref"], fields["alt"])] = float(fields["tlod"])
            except (KeyError, ValueError):
                pass
        candidates = parse_dump(dump)
        records = read_vcf(native_output)
        if not candidates or not records:
            raise AssertionError("no candidates or VCF records parsed")

        # Locus with more than one concrete ALT is scored by the merged
        # multiallelic model in calculate_somatic; the biallelic formula does
        # not apply there, so exclude those (keyed by internal pos + ref).
        alt_count = {}
        for (position, ref, _alt) in candidates:
            alt_count[(position, ref)] = alt_count.get((position, ref), 0) + 1

        deltas = []
        skipped = 0
        for (position, ref, alt), (refs, alts) in candidates.items():
            if alt_count.get((position, ref), 0) > 1:
                skipped += 1
                continue
            # Dump coordinates are one less than the serialized VCF position.
            key = (position + 1, ref, alt)
            if key not in records:
                continue
            if len(refs) != len(alts):
                raise AssertionError("matrix width mismatch")
            expected = tlod_biallelic(refs, alts)
            observed = records[key]
            serialized = observed
            full_observed = full_tlod.get((position, ref, alt))
            if full_observed is not None:
                observed = full_observed
            deltas.append(abs(expected - observed))
            print("SITE pos=%d ref=%s alt=%s expected=%.17g vcf=%.17g full=%s diff=%.3g" %
                  (position + 1, ref, alt, expected, serialized,
                   ("%.17g" % full_observed) if full_observed is not None else "-",
                   expected - observed))
        if not deltas:
            raise AssertionError("no comparable single-ALT TLOD loci")
        maximum = max(deltas)
        mean = sum(deltas) / len(deltas)
        # Full-precision engine TLOD (FASTGATK_TLOD_FULL echo) proves the
        # Python/Java SomaticLikelihoodsEngine transcription matches the
        # native engine at machine precision (observed max ~1e-14).  The VCF
        # TLOD serialization is 8 significant digits, so a VCF-only
        # measurement caps at ~1e-7 relative.
        assert maximum < (1e-8 if full_tlod else 1e-4), {
            "maximum_tlod_formula_delta": maximum,
            "full_precision_loci": len(full_tlod),
            "deltas": deltas,
        }
        print(json.dumps({
            "status": "pass",
            "loci_checked": len(deltas),
            "per_locus_deltas": deltas,
            "multiallelic_loci_skipped": skipped,
            "maximum_tlod_formula_delta": maximum,
            "mean_tlod_formula_delta": mean,
            "engine": "python-transcribed SomaticLikelihoodsEngine on "
                      "native fragment matrix",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
