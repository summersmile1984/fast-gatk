#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the sample column of a materialized spanning
locus: GATK inherits the hom-ref call's **source** GQ and DP.

The dense traversal publishes a row at every coordinate an input record's span
covers (``verify_genotype_gvcf_dense_materialize_gatk_oracle.py`` pins the row
synthesis itself).  That row is a projection of its single source record onto
``[REF, *]``, and GATK's ``cleanupGenotypeAnnotations()`` then rewrites the
call as hom-ref **without recomputing GQ**: the RGQ it emits is the source
record's own FORMAT GQ, and DP is the source record's DP.

Measured on GATK's own chr20 gVCF (pinned 4.6.2.0), which is why these exact
values are pinned below:

* coordinate 10008964/10008965 -- covered by the record starting at 10008952
  (``REF=CACACACACACACA``, sample ``GT=2/4, DP=34, GQ=99``); GATK emits
  ``0/0:34:99`` with ``FORMAT=GT:DP:RGQ``.
* coordinate 10076991 -- covered by the record starting at 10076989
  (``REF=CAA``, sample ``GT=2/3, DP=31, GQ=41``); GATK emits ``0/0:31:41``.

A covered coordinate whose spanning record has **no called deletion** is still
visited and published, as a *monomorphic* reference row (``ALT=.``,
``QUAL=Infinity``, ``FILTER`` unset, ``INFO=DP;MLEAC=.;MLEAF=.``) -- measured at
10062936-10062938 (source 10062935, ``GT=0/4`` with ``DP=23, GQ=16``, its called
ALT being longer than its REF), 10087821-10087822 (source 10087820, ``GT=0/2``,
``DP=53, GQ=99``) and 10098309-10098310 (source 10098308, hom-ref ``GT=0/0``,
``DP=2, GQ=6``).  Native used to drop those seven rows entirely because its
materialization required a *called* deletion; the sample column again inherits
the source DP/GQ (23:16, 53:99, 2:6) and the QUAL is the ``+Infinity`` that
``GenotypingEngine`` assigns when the whole posterior mass is on "no variant".

Native used to re-read GQ from its own output record, where the compute stage
had already replaced it with the PL-derived value ``0``: because the hom-ref
conversion was additionally guarded on ``gq > 0``, the sample column came out
``./.:34:0`` / ``./.:31:41``-shaped no-call -- the GT stayed ``./.`` and RGQ was
zeroed, while DP (inherited separately) stayed correct.  That is a whole
(coordinate, sample) disagreement with GATK on every such row.

Gated here (the two source GQ values differ, 99 != 41, so no constant
satisfies both):

1. native publishes a row at each pinned coordinate;
2. each native row is byte-identical to GATK's row at the same coordinate
   (CHROM..sample, full text);
3. GATK's own rows at those coordinates are hom-ref ``0/0`` with the source GQ
   (``0/0:34:99`` / ``0/0:31:41``) and ``FORMAT=GT:DP:RGQ`` -- the measured
   truth, asserted independently of native so the gate cannot be satisfied by
   native and GATK being wrong together;
4. neither pinned row is a no-call and neither RGQ is zero.

Exit status: 0 when every gated case passes, non-zero otherwise.
``--expect-divergence`` turns the run into a diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

# The dense window that materializes the pinned coordinates, and the coordinates
# themselves with the sample column GATK publishes there.  The source GQ values
# (99 and 41) differ on purpose: a fix that hard-codes one of them fails here.
WINDOW = "20:10000000-10099999"
PINNED = {
    10008964: "0/0:34:99",   # source is the record at 10008952 (GQ=99, DP=34)
    10008965: "0/0:34:99",   # same source, same span
    10076991: "0/0:31:41",   # source is the record at 10076989 (GQ=41, DP=31)
}
# Coordinates whose spanning source record has NO called deletion: GATK still
# visits them and publishes a monomorphic reference row (ALT='.', QUAL=Infinity).
PINNED_MONOMORPHIC = {
    10062936: "0/0:23:16",   # source 10062935 (GT 0/4, DP=23, GQ=16)
    10062937: "0/0:23:16",
    10062938: "0/0:23:16",
    10087821: "0/0:53:99",   # source 10087820 (GT 0/2, DP=53, GQ=99)
    10087822: "0/0:53:99",
    10098309: "0/0:2:6",     # source 10098308 (GT 0/0, DP=2, GQ=6)
    10098310: "0/0:2:6",
}
# REF-only rows whose ExcessHet GATK renders through htsjdk's raw-double path
# ("0.00", probed: formatVCFDouble(0.0) == "0.00") instead of the "%.4f" String
# that ExcessHet.java builds for every other row ("0.0000").  Native used to
# print four decimals everywhere, which cost these five rows.
PINNED_REF_ONLY_EXCESS_HET = {
    10041698: "0.00",
    10077008: "0.00",
    10077010: "0.00",
    10098308: "0.00",
    10099270: "0.00",
}
DENSE = "--include-non-variant-sites"


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stderr={result.stderr[-1600:]}")


def read_rows(path: Path) -> dict[int, str]:
    rows: dict[int, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows[int(fields[1])] = line
    return rows


def columns(row: str) -> tuple[str, str]:
    fields = row.split("\t")
    return fields[8], fields[9]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expect-divergence", action="store_true",
                        help="report the measured state and always exit 0")
    arguments = parser.parse_args()
    strict = not arguments.expect_divergence

    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_GENOTYPE_GVCF_BINARY",
        root / "fastgatk-native/build/fastgatk-genotype-gvcf"))
    reference = root / "testdata/chr20/reference/GRCh37.chr20.fa"
    source = root / ("gatk-source/src/test/resources/org/broadinstitute/"
                     "hellbender/tools/haplotypecaller/"
                     "expected.testGVCFMode.gatk4.g.vcf")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    required = (native, reference, Path(f"{reference}.fai"), source, gatk, java)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified(
            "verify_genotype_gvcf_spanning_source_gq_gatk_oracle.py", gatk, java)
        raise SystemExit("missing real chr20 spanning-locus oracle inputs")

    violations: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-span-gq-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        run([str(java), "-Xmx4g", "-jar", str(gatk), "GenotypeGVCFs",
             "-R", str(reference), "-V", str(source), "-L", WINDOW, DENSE,
             "-O", str(gatk_vcf), "--create-output-variant-index", "false"])
        run([str(native), "-R", str(reference), "-V", str(source), "-L", WINDOW,
             DENSE, "--gatk-compatible-annotations", "-O", str(native_vcf)])
        gatk_rows = read_rows(gatk_vcf)
        native_rows = read_rows(native_vcf)

        pinned = {**PINNED, **PINNED_MONOMORPHIC}
        # Control: the pinned rows must come from sources with different GQ, so
        # that a constant cannot satisfy every case below.
        distinct_rgq = {sample.split(":")[2] for sample in pinned.values()}
        if len(distinct_rgq) < 2:
            violations.append("control: the pinned RGQ values are not distinct")

        for position in sorted(pinned):
            expected = pinned[position]
            truth = gatk_rows.get(position)
            row = native_rows.get(position)
            if truth is None:
                violations.append(f"POS {position}: GATK published no row")
                continue
            # 3. the measured truth, asserted against GATK itself.
            gatk_format, gatk_sample = columns(truth)
            if gatk_format != "GT:DP:RGQ" or gatk_sample != expected:
                violations.append(
                    f"POS {position}: GATK sample is {gatk_format} {gatk_sample}, "
                    f"pinned truth is GT:DP:RGQ {expected}")
            # 4. the shape rules, stated independently of the exact numbers.
            if gatk_sample.split(":")[0] != "0/0":
                violations.append(f"POS {position}: GATK call is not hom-ref")
            if gatk_sample.split(":")[2] == "0":
                violations.append(f"POS {position}: GATK RGQ is zero")
            if position in PINNED_MONOMORPHIC:
                # The shape of a covered coordinate with no called deletion: a
                # monomorphic row with an infinite QUAL, asserted against GATK
                # itself rather than inferred from native.
                fields = truth.split("\t")
                if fields[4] != "." or fields[5] != "Infinity":
                    violations.append(
                        f"POS {position}: GATK row is not the measured monomorphic "
                        f"shape (ALT={fields[4]!r}, QUAL={fields[5]!r})")
            # 1. native must publish the row.
            if row is None:
                violations.append(f"POS {position}: native published no row")
                continue
            # 2. byte-identical rows.
            if row != truth:
                violations.append(
                    f"POS {position}: native {columns(row)} != GATK "
                    f"{columns(truth)}")

        # REF-only rows that must carry GATK's raw-double ExcessHet rendering.
        for position, expected in sorted(PINNED_REF_ONLY_EXCESS_HET.items()):
            truth = gatk_rows.get(position)
            row = native_rows.get(position)
            if truth is None or row is None:
                violations.append(
                    f"POS {position}: missing row (GATK={truth is not None}, "
                    f"native={row is not None})")
                continue
            if f"ExcessHet={expected}" not in truth.split("\t")[7]:
                violations.append(
                    f"POS {position}: GATK ExcessHet is not the measured {expected}")
            if row != truth:
                violations.append(
                    f"POS {position}: REF-only ExcessHet row is not byte-identical")

    report = {
        "gate": "genotype-gvcf-spanning-source-gq",
        "gatk_version": "4.6.2.0",
        "intervals": WINDOW,
        "ref_only_excess_het": {str(position): token
                                for position, token
                                in sorted(PINNED_REF_ONLY_EXCESS_HET.items())},
        "pinned": {str(position): sample
                   for position, sample in sorted({**PINNED, **PINNED_MONOMORPHIC}.items())},
        "violations": violations,
        "status": "pass" if not violations else "divergence",
        "mode": "strict" if strict else "expect-divergence",
    }
    for violation in violations:
        print(f"VIOLATION: {violation}")
    print(json.dumps(report, sort_keys=True))
    if not strict:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
