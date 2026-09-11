#!/usr/bin/env python3
"""Strict oracle for the SelectVariants reference-only record (pinned GATK 4.6.2.0).

Observation under test
----------------------
``fastgatk-native/scripts/verify_select_variants.py`` used to assert
``ref_only_records == []`` (line 820): native *dropped* an all-hom-ref record
under ``--remove-unused-alternates``.  An audit measured that pinned GATK
**keeps** that record and emits it with ``ALT='.'``.  This script re-measures
both tools on identical inputs and gates the result.

GATK's rule (measured against the pinned jar, then confirmed in GATK source)
--------------------------------------------------------------------------
1. ``SelectVariants.subsetGenotypesBySampleNames``
   (``SelectVariants.java:1195-1233``) trims the ALT alleles no genotype uses
   and hands the ``VariantContext`` to ``GATKVariantContextUtils.trimAlleles``
   (``SelectVariants.java:1233`` -> ``GATKVariantContextUtils.java:1455``).
   That returns a **single-allele** context (reference only) rather than a
   dropped record, and the writer renders the missing ALT as ``.``.
2. The record is therefore *not* removed by ``--remove-unused-alternates``.  It
   is removed only when ``--exclude-non-variants`` is also set:
   ``SelectVariants.java:708-715`` runs *after* the allele trimming::

       if (excludeNonVariants) {
           final boolean nonVariant = ! result.isPolymorphicInSamples() || GATKVariantContextUtils.isSpanningDeletionOnly(result);
           if (nonVariant) { return; }
       }

   and the flag defaults to ``false`` (``SelectVariants.java:265-266``).
3. For that reference-only subset the genotype payload is rebuilt by
   ``AlleleSubsettingUtils.subsetAlleles``, whose ``:94-100`` says::

       if (newLikelihoods.length > 1) { ...recompute GQ... }
       else {  //if we subset to just ref allele, keep the GQ
           newLog10GQ = g.getGQ()/-10.0;
       }

   i.e. the **original GQ is preserved** (not re-derived from the collapsed PL
   vector) and GQ is emitted only when the input genotype had one
   (``AlleleSubsettingUtils.java:128``).  PL collapses to the single all-ref
   cell, which is old PL index 0 (``subsettedPLIndices``, :426-441), and AD
   keeps the reference count (``getNewAlleleBasedReadCountAnnotation``).
4. Site annotations are then refreshed by htsjdk's chromosome-count
   calculation (``GATKVariantContextUtils.java:1243``), so a reference-only
   record carries ``AN`` (including ``AN=0`` for an all-no-call genotype set)
   and no ``AC``/``AF``.

Measured truth (pinned jar, this script re-runs every case)
-----------------------------------------------------------
The fixture is ``chr1 20 . A C,G 50 PASS <INFO> GT:AD:GQ:PL 0/0:30,0,0:30:0,30,60,40,70,80``
unless a case says otherwise::

    default (no flags)                  -> kept untrimmed: ALT C,G, INFO unchanged
    --remove-unused-alternates          -> kept, ALT '.', INFO AN=2, 0/0:30:30:0
    --exclude-non-variants              -> dropped (0 rows)
    --remove-unused-alternates + the above -> dropped (0 rows)

so the retained/dropped behaviour is option-dependent, and the option that
decides it is ``--exclude-non-variants``, not ``--remove-unused-alternates``.

Ordering matters, and GATK's order is: trim (:696) -> filtered genotypes to
no-call (:698) -> non-variant test (:708-715).  A record whose only ALT carrier
is a filtered genotype therefore survives trimming, becomes non-variant when
that genotype is converted to no-call, and is dropped under
``--exclude-non-variants`` (case ``rua-env-drops-after-nocall-conversion``).

Scope and known, separately-scoped divergences
----------------------------------------------
This gate byte-compares data rows.  Three *pre-existing* divergences are outside
its subject (the reference-only record) and are therefore reported rather than
gated; all three involve key ordering, AF formatting or an ordering difference
that predates this gate, not which reference-only records are emitted:

* FORMAT key order: after allele subsetting GATK re-emits FORMAT keys in its own
  rebuild order (``GT:AD:GQ:PL``), while native preserves the input header
  order.  Case ``diagnostic-input-format-order`` (the literal fixture of
  ``verify_select_variants.py``, which declares ``GT:AD:PL:GQ``) is gated on
  site/INFO byte identity plus per-key FORMAT value equality and prints the key
  order of both sides; the difference is not allowed to hide a value change.
* INFO key order / AF precision on records that still carry ALTs (GATK writes
  ``AC=...;AF=...;AN=...`` with ``%.3g``-style AF; native keeps the input order
  and its own float formatting).  Case ``rua-mixed-file`` gates the
  reference-only row byte-for-byte and reports the sibling row's divergence.
* ``diagnostic-nocall-before-trim`` (reported, not gated): without
  ``--exclude-non-variants`` the same filtered-ALT fixture shows that native
  converts filtered genotypes to no-call *before* trimming, so it trims the
  allele away while GATK keeps ``ALT=C``; both now emit one record.  Fixing that
  ordering is a separate change to the ``--set-filtered-gt-to-nocall`` contract
  (already covered by ``verify_select_variants_filtered_nocall_oracle.py``), so
  this oracle only reports it.

Exit status: 0 when every gated check passes, non-zero otherwise (so it can be
registered in CTest).  ``--expect-divergence`` turns the run into a pure
diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
    "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)

# FORMAT declared in GATK's own post-subsetting key order, so that FORMAT key
# ordering (a separate, pre-existing divergence: see the module docstring) cannot
# mask the record-retention behaviour this oracle is about.
CANONICAL = HEADER
REF_ONLY_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:GQ:PL\t0/0:30,0,0:30:0,30,60,40,70,80\n"
GQ99_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:GQ:PL\t0/0:30,0,0:99:0,30,60,40,70,80\n"
NOCALL_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:GQ:PL\t./.:30,0,0:30:0,30,60,40,70,80\n"
SITE_AC_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\tAC=0,0;AN=2;AF=0,0\tGT:AD:GQ:PL\t0/0:30,0,0:30:0,30,60,40,70,80\n"
ALREADY_REF_ROW = "chr1\t20\t.\tA\t.\t50\tPASS\t.\tGT:AD:GQ:PL\t0/0:30:30:0\n"
HAPLOID_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:GQ:PL\t0:30,0,0:30:0,30,60\n"
# The record carries no GQ value even though the header declares it.
NO_GQ_ROW = "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:PL\t0/0:30,0,0:0,30,60,40,70,80"

# The literal fixture and arguments of verify_select_variants.py:801-815: FORMAT
# declared in input order (GT:AD:PL:GQ), which is also the sequence that file
# asserts.  GATK re-emits the keys in its own order; see the module docstring.
LEGACY_ORDER_FIXTURE = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:PL:GQ\t0/0:30,0,0:0,30,60,40,70,80:30\n"
)

RUA = ["--remove-unused-alternates"]
ENV = ["--exclude-non-variants"]

# A filtered ALT carrier: GATK trims the unused ALTs first
# (SelectVariants.java:696) and only then applies --set-filtered-gt-to-nocall
# (:698), so the record still carries ALT C at trim time; the later no-call
# conversion is what makes the record non-variant, and only then does
# --exclude-non-variants (:708-715) drop it.  A pre-trim non-variant test cannot
# see that, which is why the drop has to be re-evaluated after trimming.
FT_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100>\n"
    "##FILTER=<ID=LowDP,Description=Low depth>\n"
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n"
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n"
    "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)
FT_ROW = ("chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:FT:AD:GQ:PL\t"
          "0/1:LowDP:15,15,0:50:50,0,80,99,99,99\t0/0:PASS:30,0,0:30:0,30,60,40,70,80\n")
NOCALL = ["--set-filtered-gt-to-nocall"]

CASES = (
    {
        "case": "rua-keeps-ref-only",
        "why": "the divergence under audit: native dropped the all-hom-ref record "
               "under --remove-unused-alternates, pinned GATK keeps it with ALT='.'",
        "header": CANONICAL,
        "records": [REF_ONLY_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0"],
    },
    {
        "case": "default-keeps-untrimmed",
        "why": "without --remove-unused-alternates GATK emits the record untouched "
               "(ALT C,G, INFO unchanged); guards against over-trimming",
        "header": CANONICAL,
        "records": [REF_ONLY_ROW],
        "args": [],
        "expect_gatk_rows": ["chr1\t20\t.\tA\tC,G\t50\tPASS\t.\tGT:AD:GQ:PL\t0/0:30,0,0:30:0,30,60,40,70,80"],
    },
    {
        "case": "exclude-non-variants-drops",
        "why": "the option that actually decides retention: --exclude-non-variants "
               "drops the same record (SelectVariants.java:708-715, default false at :265-266)",
        "header": CANONICAL,
        "records": [REF_ONLY_ROW],
        "args": ENV,
        "expect_gatk_rows": [],
    },
    {
        "case": "rua-plus-exclude-non-variants-drops",
        "why": "both flags: still dropped, so the fix must not resurrect records "
               "that GATK removes",
        "header": CANONICAL,
        "records": [REF_ONLY_ROW],
        "args": RUA + ENV,
        "expect_gatk_rows": [],
    },
    {
        "case": "rua-preserves-original-gq",
        "why": "AlleleSubsettingUtils.java:94-100 keeps the *original* GQ when the "
               "subset is the reference allele alone (input GQ=99 here)",
        "header": CANONICAL,
        "records": [GQ99_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:99:0"],
    },
    {
        "case": "rua-all-nocall-an-zero",
        "why": "a reference-only subset of an all-no-call genotype set still writes "
               "AN, at zero (htsjdk chromosome-count refresh)",
        "header": CANONICAL,
        "records": [NOCALL_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=0\tGT:AD:GQ:PL\t./.:30:30:0"],
    },
    {
        "case": "rua-clears-ac-and-af",
        "why": "a reference-only record has no ALT, so AN survives the refresh and "
               "AC/AF do not",
        "header": CANONICAL,
        "records": [SITE_AC_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0"],
    },
    {
        "case": "rua-existing-ref-only-untouched",
        "why": "an input record that is already reference-only takes no trimming "
               "path in GATK and keeps INFO '.' -- no AN is invented",
        "header": CANONICAL,
        "records": [ALREADY_REF_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\t.\tGT:AD:GQ:PL\t0/0:30:30:0"],
    },
    {
        "case": "rua-omits-absent-gq",
        "why": "GQ is emitted only when the input genotype had one "
               "(AlleleSubsettingUtils.java:128)",
        "header": CANONICAL,
        "records": [NO_GQ_ROW + "\n"],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:PL\t0/0:30:0"],
    },
    {
        "case": "rua-haploid",
        "why": "the collapse to the all-ref PL cell is ploidy-aware (haploid GT 0)",
        "header": CANONICAL,
        "records": [HAPLOID_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=1\tGT:AD:GQ:PL\t0:30:30:0"],
    },
    {
        "case": "rua-mixed-file",
        "why": "a reference-only record next to a surviving variant record: the "
               "reference-only row must stay byte-identical",
        "header": CANONICAL,
        "records": ["chr1\t10\t.\tA\tC\t50\tPASS\t.\tGT:AD:GQ:PL\t0/1:15,15:50:50,0,80\n",
                    REF_ONLY_ROW],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t10\t.\tA\tC\t50\tPASS\t.\tGT:AD:GQ:PL\t0/1:15,15:50:50,0,80",
                             "chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0"],
        # The sibling row diverges in INFO only (native refreshes AC/AN/AF of a
        # record that needs no trimming; GATK leaves INFO '.').  That is the
        # pre-existing annotation-refresh divergence described in the docstring,
        # so byte identity is gated for the reference-only POS only.
        "gate_positions": ["20"],
    },
    {
        "case": "rua-env-drops-after-nocall-conversion",
        "why": "the non-variant test must be re-evaluated after trimming: GATK "
               "trims first (SelectVariants.java:696), converts filtered "
               "genotypes to no-call second (:698) and only then applies "
               "--exclude-non-variants (:708-715), so this record is dropped",
        "header": FT_HEADER,
        "records": [FT_ROW],
        "args": NOCALL + RUA + ENV,
        "expect_gatk_rows": [],
    },
    {
        "case": "diagnostic-nocall-before-trim",
        "why": "reported, not gated: the same fixture without "
               "--exclude-non-variants shows a *separate* pre-existing ordering "
               "difference.  GATK keeps ALT C because it trims while the "
               "filtered genotype still calls it (SelectVariants.java:696 runs "
               "before :698); native converts to no-call first and therefore "
               "trims the allele away.  Both tools now emit one record, with "
               "different ALTs and FORMAT/INFO key orders",
        "header": FT_HEADER,
        "records": [FT_ROW],
        "args": NOCALL + RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\tC\t50\tPASS\tAC=0;AF=0.00;AN=2\tGT:AD:FT:GQ:PL\t"
                             "./.:15,15:LowDP:50:50,0,80\t0/0:30,0:PASS:30:0,30,60"],
        "gate": False,
    },
    {
        "case": "diagnostic-input-format-order",
        "why": "the literal verify_select_variants.py fixture (FORMAT declared "
               "GT:AD:PL:GQ): the reference-only record and its values must match "
               "GATK; only the pre-existing FORMAT key-order divergence is excused",
        "header": LEGACY_ORDER_FIXTURE,
        "records": [],
        "args": RUA,
        "expect_gatk_rows": ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0"],
        "gate_format_keys": True,
    },
)

# Cases with ``"gate": False`` are measured and printed but do not fail the run:
# they document divergences that this oracle is not the subject of.
GATED = {case["case"] for case in CASES if case.get("gate", True)}
REPORTED_ONLY = {case["case"] for case in CASES if not case.get("gate", True)}


def data_rows(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream.read().splitlines()
                if line and not line.startswith("#")]


def row_position(row: str) -> str:
    return row.split("\t")[1]


def format_map(row: str) -> dict[str, str]:
    """``FORMAT key -> value`` for a single-sample data row."""
    fields = row.split("\t")
    if len(fields) < 10:
        return {}
    keys = fields[8].split(":")
    values = fields[9].split(":")
    return dict(zip(keys, values))


def index_feature(java: str, gatk: str, path: pathlib.Path, timeout: int) -> None:
    result = subprocess.run(
        [java, "-jar", gatk, "IndexFeatureFile", "-I", str(path)],
        text=True, capture_output=True, check=False, timeout=timeout,
    )
    if result.returncode != 0:
        raise AssertionError(f"GATK IndexFeatureFile failed for {path}:\n{result.stderr}")


def run_case(case: dict, work: pathlib.Path, java: str, gatk: str, native: str,
             timeout: int) -> dict:
    case_work = work / case["case"]
    case_work.mkdir(parents=True, exist_ok=True)
    # Plain (uncompressed) input: GATK IndexFeatureFile can index it with its
    # Tribble index without needing an external bgzip binary, and both engines
    # read the same bytes.
    input_vcf = case_work / "input.vcf"
    input_vcf.write_text(case["header"] + "".join(case["records"]), encoding="utf-8")
    index_feature(java, gatk, input_vcf, timeout)

    gatk_output = case_work / "gatk.vcf"
    gatk_result = subprocess.run(
        [java, "-Xmx1g", "-jar", gatk, "SelectVariants",
         "-V", str(input_vcf), "-O", str(gatk_output), *case["args"]],
        text=True, capture_output=True, check=False, timeout=timeout,
    )
    native_output = case_work / "native.vcf.gz"
    native_result = subprocess.run(
        [native, "-V", str(input_vcf), "-O", str(native_output), *case["args"]],
        text=True, capture_output=True, check=False, timeout=timeout,
    )

    result: dict = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["case"] in GATED,
        "reported_only": case["case"] in REPORTED_ONLY,
        "args": list(case["args"]),
        "gate_positions": case.get("gate_positions"),
        "gate_format_keys": bool(case.get("gate_format_keys")),
        "expect_gatk_rows": case["expect_gatk_rows"],
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_stderr_tail": gatk_result.stderr[-600:] if gatk_result.returncode else "",
        "native_stderr_tail": native_result.stderr[-600:] if native_result.returncode else "",
        "violations": [],
    }
    if gatk_result.returncode != 0 or native_result.returncode != 0:
        result["violations"].append(
            f"run failed: gatk_exit={gatk_result.returncode} "
            f"native_exit={native_result.returncode}")
        result["gatk_rows"] = []
        result["native_rows"] = []
        result["format_key_order"] = {}
        return result

    gatk_rows = data_rows(gatk_output)
    native_rows = data_rows(native_output)
    result["gatk_rows"] = gatk_rows
    result["native_rows"] = native_rows

    violations: list[str] = []
    # Fixture-validity gate: the GATK side must still reproduce the measured
    # contract, so the oracle cannot pass vacuously if GATK's behaviour moves.
    if gatk_rows != case["expect_gatk_rows"]:
        violations.append(
            f"{case['case']}: pinned GATK no longer reproduces the measured rows: "
            f"expected {case['expect_gatk_rows']} but measured {gatk_rows}")

    if case.get("gate_format_keys"):
        result["format_key_order"] = {
            "gatk": [row.split("\t")[8] for row in gatk_rows],
            "native": [row.split("\t")[8] for row in native_rows],
        }
        if result["gated"] and len(gatk_rows) != len(native_rows):
            violations.append(
                f"{case['case']}: record count differs: GATK={len(gatk_rows)} "
                f"NATIVE={len(native_rows)}")
        for gatk_row, native_row in zip(gatk_rows, native_rows):
            gatk_site = gatk_row.split("\t")[:8]
            native_site = native_row.split("\t")[:8]
            if gatk_site != native_site and result["gated"]:
                violations.append(
                    f"{case['case']}: site/INFO columns differ at POS "
                    f"{row_position(gatk_row)}: GATK={gatk_site} NATIVE={native_site}")
            elif format_map(gatk_row) != format_map(native_row) and result["gated"]:
                violations.append(
                    f"{case['case']}: FORMAT values differ at POS {row_position(gatk_row)}: "
                    f"GATK={format_map(gatk_row)} NATIVE={format_map(native_row)}")
    elif case.get("gate_positions"):
        wanted = set(case["gate_positions"])
        gatk_selected = {row_position(row): row for row in gatk_rows if row_position(row) in wanted}
        native_selected = {row_position(row): row for row in native_rows if row_position(row) in wanted}
        result["gated_rows"] = {
            "positions": sorted(wanted),
            "gatk": [gatk_selected[key] for key in sorted(gatk_selected)],
            "native": [native_selected[key] for key in sorted(native_selected)],
        }
        if gatk_selected != native_selected and result["gated"]:
            violations.append(
                f"{case['case']}: gated rows at POS {sorted(wanted)} are not "
                f"byte-identical: GATK={[gatk_selected[key] for key in sorted(gatk_selected)]} "
                f"NATIVE={[native_selected[key] for key in sorted(native_selected)]}")
        other_positions = [row_position(row) for row in gatk_rows if row_position(row) not in wanted]
        result["reported_positions"] = other_positions
    else:
        if gatk_rows != native_rows:
            differences = [[index, gatk_row, native_row]
                           for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows))
                           if gatk_row != native_row]
            result["cross_tool_differences"] = differences[:4]
            if result["gated"] and len(gatk_rows) != len(native_rows):
                violations.append(
                    f"{case['case']}: record count differs: GATK={len(gatk_rows)} "
                    f"NATIVE={len(native_rows)}")
            if result["gated"] and differences:
                violations.append(
                    f"{case['case']}: data rows are not byte-identical "
                    f"(index, GATK, NATIVE): {differences[:4]}")

    result["violations"] = violations
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the SelectVariants reference-only record "
                    "against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_SELECT_VARIANTS_BINARY"),
        help="native SelectVariants binary (default: "
             "$FASTGATK_SELECT_VARIANTS_BINARY or "
             "fastgatk-native/build/fastgatk-select-variants)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the measured GATK/native rows and exit 0 "
             "instead of gating byte parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-select-variants")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))

    assets = [native, pathlib.Path(java), gatk]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_select_variants_refonly_gatk_oracle.py",
            pathlib.Path(java), gatk)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-select-variants-refonly-oracle-") as directory:
        work = pathlib.Path(directory)
        for case in selected:
            results.append(run_case(case, work, java, str(gatk), str(native),
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: SelectVariants.subsetGenotypesBySampleNames "
          "(SelectVariants.java:1195-1233) trims unused ALTs to a reference-only "
          "VariantContext (GATKVariantContextUtils.java:1455) instead of dropping it; "
          "the record is removed only by --exclude-non-variants "
          "(SelectVariants.java:708-715, default false at :265-266).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker} args={result['args']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    GATK   rows: {result['gatk_rows']}")
        print(f"    NATIVE rows: {result['native_rows']}")
        if result.get("cross_tool_differences"):
            print(f"    cross-tool differences (index, GATK, NATIVE): "
                  f"{result['cross_tool_differences']}")
        if result["gatk_rows"] != result["expect_gatk_rows"]:
            print(f"    expected GATK rows (measured): {result['expect_gatk_rows']}")
        for side, keys in sorted(result.get("format_key_order", {}).items()):
            print(f"    FORMAT key order [{side}]: {keys}")
        if result.get("gated_rows"):
            print(f"    gated rows (POS {result['gated_rows']['positions']}): "
                  f"GATK={result['gated_rows']['gatk']} NATIVE={result['gated_rows']['native']}")
        if result.get("reported_positions"):
            print(f"    reported, not gated, POS: {result['reported_positions']}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")
        if result["gatk_stderr_tail"]:
            print(f"    gatk_stderr_tail: {result['gatk_stderr_tail']!r}")
        if result["native_stderr_tail"]:
            print(f"    native_stderr_tail: {result['native_stderr_tail']!r}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "SelectVariants reference-only record retention under --remove-unused-alternates",
        "binary": str(native),
        "acceptance_criterion": (
            "on every gated case the data rows written by pinned GATK 4.6.2.0 and "
            "native must be byte-identical for the same input and arguments "
            "(record retention, ALT='.', GT/AD/PL/GQ values and AN/AC/AF "
            "annotations included), and the GATK side must still reproduce the "
            "rows measured when this oracle was written"),
        "known_divergence": (
            "reported, not gated, and outside this oracle's subject: (1) GATK "
            "re-emits FORMAT keys in its own rebuild order after allele "
            "subsetting while native preserves input header order "
            "(case diagnostic-input-format-order gates values per key, not key "
            "order); (2) on records that need no trimming native refreshes "
            "AC/AN/AF where GATK leaves INFO untouched (case rua-mixed-file "
            "gates the reference-only row only)"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": sorted(GATED),
        "reported_only_cases": sorted(REPORTED_ONLY),
        "violations": violations,
        "results": results,
    }
    print(json.dumps(payload, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
