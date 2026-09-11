#!/usr/bin/env python3
"""Strict oracle for the gVCF `/ <NON_REF>` symbolic-allele AF prior divergence.

Root class under test
---------------------
`AlleleFrequencyCalculator` derives every allele's Dirichlet pseudocount from
the allele *length*, not from a fixed SNP/indel flag::

    // gatk-source/.../afcalc/AlleleFrequencyCalculator.java:175-176
    final double[] priorPseudocounts = alleles.stream()
        .mapToDouble(a -> a.isReference() ? refPseudocount
                        : (a.length() == refLength ? snpPseudocount : indelPseudocount))
        .toArray();

with ``refLength = vc.getReference().length()``.  htsjdk's symbolic alleles have
length 0 while ``Allele.SPAN_DEL`` has length 1 (measured against the pinned
htsjdk 4.2.0: ``SPAN_DEL.length()==1``, ``NON_REF.length()==0``,
``<DEL>.length()==0``).  Consequences for a 1 bp REF record:

* concrete SNV/1-for-1 ALT -> SNP pseudocount;
* concrete indel ALT      -> indel pseudocount;
* symbolic ``*``           -> **SNP** pseudocount (length 1 == refLength);
* symbolic ``<NON_REF>``   -> indel pseudocount (length 0 != refLength).

A previous round fixed that rule for the ordinary-VCF candidate QUAL
(``calling_pipeline.cpp`` ``spanning_pseudocount``).  The gVCF writers construct
their own pseudocount vector and still used ``indel_heterozygosity``
unconditionally for ``*`` (``hc_call.cpp``, the gVCF candidate-site AF block),
so a reference-confidence record that carries ``*`` gets a different
``P(no variant)`` and therefore a different QUAL.

Measured divergence (pre-fix, this repo, ``-ERC BP_RESOLUTION``)
----------------------------------------------------------------
``chr1:698`` of the ``bp-span-del-frameshift`` fixture::

    GATK   chr1 698 . C *,A,<NON_REF> 282.04 ... MLEAC=0,1,0 ...
    NATIVE chr1 698 . C A,*,<NON_REF> 291.07 ... MLEAC=1,0,0 ...

Everything else (including the ``chr1:697 TC>T`` concrete deletion record and
every reference block) is byte-identical.  The same 282.04/291.07 pair appears
in plain ``-ERC GVCF`` mode.

Comparison rules (and what is deliberately *not* gated here)
-----------------------------------------------------------
* Every data row is compared field by field; ``#`` provenance lines and the
  header are ignored by construction.
* ``QUAL`` is additionally re-checked on the **raw, unnormalized** rows for
  every record of every case: that is the quantity this gate exists for.
* Two divergences in this area are *separately tracked* refactors and are
  **not** gated (neither fixed nor hidden): native always renders the symbolic
  spanning deletion as the last concrete ALT (GATK places ``*`` by haplotype
  order) and native writes a redundant ``END=`` INFO on concrete-indel records
  in plain GVCF mode.  Gated cases therefore use ``-ERC BP_RESOLUTION`` (no
  ``END=``) and a canonical allele-order normalization is applied to records
  whose ALT *set* is identical but whose order differs; records whose ALT set
  differs are violations, never normalized.  The normalizations applied are
  printed per case (``allele_order_normalizations``) and are limited to a
  documented whitelist (ALT, ``AC``/``AF``/``MLEAC``/``MLEAF``, ``AD``, ``PL``,
  ``GT``/``PGT``); any other comma-valued key makes the record
  non-normalizable, i.e. it must be byte-identical.
* The plain-GVCF ``END=`` probes are run and reported as
  ``out_of_scope_notes`` (never counted as violations, never fixed here).

Exit status: 0 when every gated case's data rows are byte-identical (modulo the
documented canonicalization), non-zero otherwise, so it can be registered as a
strict test.  ``--expect-divergence`` turns the run into a diagnostic that
always exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path
import oracle_guard

REFERENCE_LENGTH = 1500
HOMOPOLYMER_START = 600          # 1-based, 10 x 'A' -> 600..609
TANDEM_START = 900               # 1-based, 8 x 'CAG' -> 900..923
SNP_POSITION = 698               # 1-based; reference base is 'C' (seed 11)
DEFAULT_INTERVAL = "chr1:500-780"

BP_RESOLUTION = ("-ERC", "BP_RESOLUTION")
GVCF = ("-ERC", "GVCF")

# `-ERC BP_RESOLUTION` keeps every reference site as its own record, which is
# what makes a byte-identical comparison possible here: in plain GVCF mode
# native additionally writes a redundant `END=` on concrete-indel records (a
# separately tracked annotation-emission divergence, reported below).
CASES: tuple[dict, ...] = (
    {
        "case": "bp-span-del-frameshift",
        "why": "reported pileup (30 ref + 12 chr1:698 C>A + 6 one-bp-deletion "
               "reads inside the chr1:698-700 CCC run, deletion encoded as "
               "sequence removal) at -ERC BP_RESOLUTION; this is the case whose "
               "chr1:698 QUAL differs (GATK 282.04 / native 291.07 pre-fix)",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
    },
    {
        "case": "bp-span-del-cigar",
        "why": "same pileup, deletion expressed as a real 138M1D161M CIGAR: "
               "the divergence must not depend on the pathological encoding",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "cigar", "drop_base": 700}},
    },
    {
        "case": "bp-span-del-deep",
        "why": "deeper deletion support (10 reads) at the same locus: the "
               "QUAL divergence must survive a different allele balance",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 10, "start": 560, "step": 4, "length": 300,
                           "mode": "cigar", "drop_base": 700}},
    },
    {
        "case": "bp-snp-only-control",
        "why": "control: the same 30 ref + 12 C>A reads without any deletion, "
               "so no `*` allele exists; rows must be byte-identical with no "
               "normalization at all",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": None},
    },
    {
        "case": "bp-insertion-only-control",
        "why": "control: a concrete 1 bp insertion, i.e. a real INDEL allele "
               "that must keep the indel pseudocount",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": None,
                   "ins": {"count": 10, "start": 440, "step": 3, "length": 300,
                           "insert_at": 698, "insert_base": "A"},
                   "del": None},
    },
    {
        "case": "bp-ref-only-control",
        "why": "control: reference-only pileup, no candidate allele at all",
        "interval": DEFAULT_INTERVAL,
        "extra": BP_RESOLUTION,
        "layout": {"ref": {"count": 20, "start": 440, "step": 3, "length": 300},
                   "snp": None, "del": None},
    },
)

# Same fixtures in plain `-ERC GVCF` mode.  They must show the *same* QUAL pair
# once the rule is fixed, but they also carry the separately tracked native-only
# `END=` INFO on concrete-indel records, so they are probes, not gates.
OUT_OF_SCOPE_CASES: tuple[dict, ...] = (
    {
        "case": "gvcf-span-del-frameshift",
        "why": "reported pileup in plain -ERC GVCF mode; QUAL reports the same "
               "282.04/291.07 pair, but native additionally writes END= on the "
               "chr1:697 concrete deletion record (separately tracked)",
        "interval": DEFAULT_INTERVAL,
        "extra": GVCF,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "frameshift", "drop_base": 700}},
    },
    {
        "case": "gvcf-cigar-del",
        "why": "concrete CIGAR deletion in plain -ERC GVCF mode; probes the "
               "same rule and shows the native-only END= annotation",
        "interval": DEFAULT_INTERVAL,
        "extra": GVCF,
        "layout": {"ref": {"count": 30, "start": 440, "step": 3, "length": 300},
                   "snp": {"count": 12, "start": 440, "step": 3, "length": 300},
                   "del": {"count": 6, "start": 560, "step": 4, "length": 300,
                           "mode": "cigar", "drop_base": 700}},
    },
)

ALL_CASES: tuple[dict, ...] = CASES + OUT_OF_SCOPE_CASES

# INFO/FORMAT keys this gate is willing to permute together with the ALT list.
# Anything else that carries a comma makes a record non-normalizable.
INFO_NUMBER_A = ("AC", "AF", "MLEAC", "MLEAF")
FORMAT_NUMBER_A: tuple[str, ...] = ()
FORMAT_NUMBER_R = ("AD",)
FORMAT_NUMBER_G = ("PL",)
FORMAT_GENOTYPE = ("GT", "PGT")
# Keys that carry commas but are not allele-indexed, so permutation leaves them
# alone: GENOTYPE/RGQ-style scalars plus GATK's Number=2 RAW_MQandDP.
SAFE_KEYS = set(INFO_NUMBER_A) | set(FORMAT_NUMBER_A) | set(FORMAT_NUMBER_R) \
    | set(FORMAT_NUMBER_G) | set(FORMAT_GENOTYPE) | {
        "DP", "GQ", "MIN_DP", "RGQ", "SB", "PS", "PID", "FT", "FILTER", ".",
        "RAW_MQandDP",
    }


def build_reference() -> str:
    random.seed(11)
    bases = list("".join(random.choice("ACGT") for _ in range(REFERENCE_LENGTH)))
    bases[HOMOPOLYMER_START - 1:HOMOPOLYMER_START - 1 + 10] = list("A" * 10)
    bases[TANDEM_START - 1:TANDEM_START - 1 + 24] = list("CAG" * 8)
    return "".join(bases)


def build_reads(layout: dict, reference_text: str) -> list[tuple[int, str, str, str]]:
    """(start, cigar, sequence, tag) for every read of a case, sorted by start."""
    placed: list[tuple[int, str, str, str]] = []
    group = layout["ref"]
    for offset in range(group["count"]):
        start = group["start"] + offset * group["step"]
        placed.append((start, f"{group['length']}M",
                       reference_text[start - 1:start - 1 + group["length"]], "r"))

    group = layout.get("snp")
    if group:
        for offset in range(group["count"]):
            start = group["start"] + offset * group["step"]
            sequence = list(reference_text[start - 1:start - 1 + group["length"]])
            sequence[SNP_POSITION - start] = group.get("alt", "A")
            placed.append((start, f"{group['length']}M", "".join(sequence), "s"))

    group = layout.get("ins")
    if group:
        for offset in range(group["count"]):
            start = group["start"] + offset * group["step"]
            offset_in_read = group["insert_at"] - start
            sequence = (reference_text[start - 1:start - 1 + offset_in_read]
                        + group["insert_base"]
                        + reference_text[start - 1 + offset_in_read:
                                         start - 1 + group["length"]])
            cigar = f"{offset_in_read}M1I{group['length'] - offset_in_read}M"
            placed.append((start, cigar, sequence, "i"))

    group = layout.get("del")
    if group:
        length, drop = group["length"], group["drop_base"]
        for offset in range(group["count"]):
            start = group["start"] + offset * group["step"]
            if group["mode"] == "frameshift":
                sequence = (reference_text[start - 1:drop - 1]
                            + reference_text[drop:start - 1 + length])
                cigar = f"{length - 1}M"
            else:
                lead = drop - start
                sequence = (reference_text[start - 1:drop - 1]
                            + reference_text[drop:start - 1 + length])
                cigar = f"{lead}M1D{length - lead - 1}M"
            placed.append((start, cigar, sequence, "d"))
    return placed


def write_sam(work: Path, name: str, layout: dict, reference_text: str) -> Path:
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             f"@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}", "@RG\tID:rg\tSM:S1"]
    reads = build_reads(layout, reference_text)
    for index, (start, cigar, sequence, tag) in enumerate(
            sorted(reads, key=lambda item: item[0])):
        lines.append(f"{tag}{index}\t0\tchr1\t{start}\t60\t{cigar}\t*\t0\t0\t"
                     f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg")
    path = work / f"{name}.sam"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def records(path: Path) -> list[list[str]]:
    """Data rows only: header/provenance lines are ignored by construction."""
    if not path.is_file():
        return []
    return [line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def split_info(text: str) -> dict:
    if not text or text == ".":
        return {}
    return dict(field.split("=", 1) if "=" in field else (field, "")
                for field in text.split(";"))


def join_info(values: dict) -> str:
    if not values:
        return "."
    return ";".join(f"{key}={value}" if value != "" else key
                    for key, value in values.items())


def genotype_index(left: int, right: int) -> int:
    """Diploid Number=G index for alleles left <= right."""
    if left > right:
        left, right = right, left
    return right * (right + 1) // 2 + left


def canonical_allele_order(alts: list[str]) -> list[str] | None:
    """Sorted concrete ALTs with a trailing <NON_REF> pinned last."""
    concrete = [alt for alt in alts if alt != "<NON_REF>"]
    order = sorted(concrete)
    if len(order) != len(set(order)):
        return None
    if "<NON_REF>" in alts:
        order.append("<NON_REF>")
    return order


def normalize_row(row: list[str]) -> tuple[list[str], bool]:
    """Canonicalize a record: sorted concrete ALTs, ``<NON_REF>`` pinned last,
    every allele-indexed field permuted with it, and phased genotypes written
    in ascending allele order.

    The canonicalization is applied *unconditionally and identically* to both
    files, so a record whose ALT order already matches is still permuted
    (identity permutation) and its phased genotype strings are still ordered
    the same way.  Returns ``(row, alt_order_changed)``.  Raises
    ``ValueError`` when the record carries a comma-valued field outside the
    documented whitelist, i.e. when permutation would be unsafe.
    """
    if len(row) < 8:
        return row, False
    alts = row[4].split(",")
    if alts == ["."] or not alts:
        return row, False
    order = canonical_allele_order(alts)
    if order is None:
        raise ValueError(f"repeated ALT allele in {row[4]!r}")
    changed = order != alts
    index = {alt: position for position, alt in enumerate(order)}
    perm = [index[alt] for alt in alts]
    alleles = len(alts) + 1

    info = split_info(row[7])
    for key, value in info.items():
        if key not in SAFE_KEYS and "," in value:
            raise ValueError(f"non-normalizable INFO field {key}={value}")

    def permute_alt_valued(text: str) -> str:
        values = text.split(",")
        if len(values) != len(alts):
            raise ValueError(f"Number=A arity mismatch: {text}")
        out = [""] * len(alts)
        for old, value in enumerate(values):
            out[perm[old]] = value
        return ",".join(out)

    def permute_ref_valued(text: str) -> str:
        values = text.split(",")
        if len(values) != alleles:
            raise ValueError(f"Number=R arity mismatch: {text}")
        out = [values[0]] * alleles
        for old in range(len(alts)):
            out[1 + perm[old]] = values[1 + old]
        return ",".join(out)

    def remap_index(value: int) -> int:
        return 0 if value == 0 else 1 + perm[value - 1]

    def permute_genotype(text: str) -> str:
        separator = "|" if "|" in text else "/" if "/" in text else None
        if separator is None:
            return str(remap_index(int(text)))
        mapped = sorted(remap_index(int(part))
                        for part in text.split(separator) if part != ".")
        return separator.join(str(value) for value in mapped)

    def permute_number_g(text: str) -> str:
        values = [int(part) for part in text.split(",")]
        if len(values) != alleles * (alleles + 1) // 2:
            raise ValueError(f"Number=G arity mismatch for {alleles} alleles")
        out = [0] * len(values)
        for left in range(alleles):
            for right in range(left, alleles):
                new_left, new_right = sorted((remap_index(left), remap_index(right)))
                out[genotype_index(new_left, new_right)] = \
                    values[genotype_index(left, right)]
        return ",".join(str(value) for value in out)

    for key in INFO_NUMBER_A:
        if key in info and info[key] not in ("", "."):
            info[key] = permute_alt_valued(info[key])

    normalized = list(row)
    normalized[4] = ",".join(order)
    normalized[7] = join_info(info)
    if len(row) > 9:
        keys = row[8].split(":")
        values = row[9].split(":")
        if len(keys) != len(values):
            raise ValueError("FORMAT/sample arity mismatch")
        out_values: list[str] = []
        for key, value in zip(keys, values):
            if value in (".", ""):
                out_values.append(value)
            elif key in FORMAT_NUMBER_A:
                out_values.append(permute_alt_valued(value))
            elif key in FORMAT_NUMBER_R:
                out_values.append(permute_ref_valued(value))
            elif key in FORMAT_NUMBER_G:
                out_values.append(permute_number_g(value))
            elif key in FORMAT_GENOTYPE:
                out_values.append(permute_genotype(value))
            elif key not in SAFE_KEYS and "," in value:
                raise ValueError(f"non-normalizable FORMAT field {key}={value}")
            else:
                out_values.append(value)
        normalized[9] = ":".join(out_values)
    return normalized, changed


def compare_rows(gatk_rows: list[list[str]],
                 native_rows: list[list[str]]) -> dict:
    """Field-level comparison with the documented canonicalization applied."""
    violations: list[str] = []
    normalizations: list[str] = []
    phasing_notes: list[str] = []
    qual_mismatches: list[str] = []
    first_diff: str = ""
    if len(gatk_rows) != len(native_rows):
        violations.append(f"row count differs: gatk={len(gatk_rows)} "
                          f"native={len(native_rows)}")
    for position in range(max(len(gatk_rows), len(native_rows))):
        gatk_raw = gatk_rows[position] if position < len(gatk_rows) else None
        native_raw = native_rows[position] if position < len(native_rows) else None
        if gatk_raw is None or native_raw is None:
            violations.append(
                f"record {position}: missing in "
                f"{'native' if native_raw is None else 'gatk'}")
            continue
        # Raw QUAL identity: the quantity this gate exists for.
        if gatk_raw[5] != native_raw[5]:
            qual_mismatches.append(
                f"{gatk_raw[0]}:{gatk_raw[1]} {gatk_raw[3]}>{gatk_raw[4]} "
                f"{gatk_raw[5]} -> {native_raw[5]}")
        if gatk_raw == native_raw:
            continue
        if set(gatk_raw[4].split(",")) != set(native_raw[4].split(",")):
            violations.append(
                f"record {position} ({gatk_raw[0]}:{gatk_raw[1]}): ALT set "
                f"differs {gatk_raw[4]!r} vs {native_raw[4]!r}")
            continue
        try:
            gatk_norm, gatk_applied = normalize_row(gatk_raw)
            native_norm, native_applied = normalize_row(native_raw)
        except ValueError as error:
            violations.append(
                f"record {position} ({gatk_raw[0]}:{gatk_raw[1]}): "
                f"{error}")
            continue
        if gatk_applied or native_applied:
            normalizations.append(
                f"{gatk_raw[0]}:{gatk_raw[1]} ALT {gatk_raw[4]!r} / "
                f"{native_raw[4]!r} -> {gatk_norm[4]!r}")
        # The phasing triple is only exempted on a record whose symbolic ALT
        # order actually differs; elsewhere it is gated like every other field.
        exclude = PHASING_KEYS if gatk_raw[4] != native_raw[4] else frozenset()
        kept, dropped = field_differences(gatk_norm, native_norm, exclude)
        if dropped:
            phasing_notes.append(
                f"{gatk_raw[0]}:{gatk_raw[1]} (ALT order differs, phasing keys "
                f"reported not gated): {dropped}")
        if kept:
            violations.append(
                f"record {position} ({gatk_raw[0]}:{gatk_raw[1]}): {kept}")
            if not first_diff:
                first_diff = kept
    return {"violations": violations,
            "allele_order_normalizations": normalizations,
            "phasing_notes": phasing_notes,
            "raw_qual_mismatches": qual_mismatches,
            "data_rows_byte_identical": (not violations and not qual_mismatches),
            "first_divergent_row_field_diff": first_diff}


def field_differences(gatk_row: list[str], native_row: list[str],
                      exclude_format: frozenset = frozenset()) -> tuple[str, str]:
    """Return ``(kept_differences, excluded_differences)`` as two strings.

    ``exclude_format`` names FORMAT keys whose differences are reported
    separately instead of being counted: the gate uses it only for the phasing
    triple (``PGT``/``PID``/``PS``) of records whose symbolic ALT order differs
    between the two files, where native's ``PGT`` string encodes the
    pre-reorder allele index of ``*`` (the separately tracked symbolic-ALT
    ordering gap).  Every other field is always gated.
    """
    differences: list[str] = []
    excluded: list[str] = []
    for index, label in ((3, "REF"), (4, "ALT"), (5, "QUAL"), (6, "FILTER")):
        if index < len(gatk_row) and index < len(native_row) \
                and gatk_row[index] != native_row[index]:
            differences.append(f"{label} {gatk_row[index]!r} -> {native_row[index]!r}")
    if len(gatk_row) > 7 and len(native_row) > 7:
        gatk_info = split_info(gatk_row[7])
        native_info = split_info(native_row[7])
        for key in sorted(set(gatk_info) | set(native_info)):
            if gatk_info.get(key) != native_info.get(key):
                differences.append(
                    f"INFO {key} {gatk_info.get(key)!r} -> {native_info.get(key)!r}")
    if len(gatk_row) > 9 and len(native_row) > 9:
        if gatk_row[8] != native_row[8]:
            differences.append(f"FORMAT keys {gatk_row[8]!r} -> {native_row[8]!r}")
        keys = gatk_row[8].split(":")
        gatk_values = gatk_row[9].split(":")
        native_values = native_row[9].split(":")
        for index, key in enumerate(keys):
            gatk_value = gatk_values[index] if index < len(gatk_values) else None
            native_value = native_values[index] if index < len(native_values) else None
            if gatk_value != native_value:
                text = f"FORMAT {key} {gatk_value!r} -> {native_value!r}"
                if key in exclude_format:
                    excluded.append(text)
                else:
                    differences.append(text)
    return ("; ".join(differences), "; ".join(excluded))


PHASING_KEYS = frozenset({"PGT", "PID", "PS"})


def qual_summary(rows: list[list[str]]) -> list[dict]:
    return [{"POS": row[1], "REF": row[3], "ALT": row[4], "QUAL": row[5]}
            for row in rows if len(row) > 5]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the gVCF/<NON_REF> symbolic-allele "
                    "allele-frequency prior divergence against pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_HC_BINARY"),
        help="native HaplotypeCaller binary (default: $FASTGATK_HC_BINARY or "
             "fastgatk-native/build/fastgatk-hc-call)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the divergence and exit 0 instead of "
             "asserting byte parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help=f"run only the named case(s); default is all of "
                             f"{[case['case'] for case in ALL_CASES]}")
    parser.add_argument("--threads", type=int, default=2,
                        help="--threads value for the native run (default 2)")
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(arguments.native) if arguments.native else (
        root / "fastgatk-native/build/fastgatk-hc-call")
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))

    assets = [native, Path(java), gatk]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified('verify_hc_gvcf_symbolic_prior_gatk_oracle.py', Path(java), gatk)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}))
        return 0

    selected = [case for case in ALL_CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    gated = {case["case"] for case in CASES}

    strict_mode = not arguments.expect_divergence
    results: list[dict] = []
    run_errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-gvcf-prior-oracle-") as directory:
        work = Path(directory)
        reference_text = build_reference()
        reference = work / "ref.fa"
        reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
        (work / "ref.fa.fai").write_text(
            f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
            encoding="utf-8")
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=300)

        for case in selected:
            sam = write_sam(work, case["case"], case["layout"], reference_text)
            bam = work / f"{case['case']}.bam"
            subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam",
                            "-I", str(sam), "-O", str(bam), "-SO", "coordinate",
                            "--CREATE_INDEX", "true"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=300)
            common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
                      "--min-pruning", "1",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false"]
            common += list(case["extra"])
            gatk_vcf = work / f"gatk.{case['case']}.vcf"
            native_vcf = work / f"native.{case['case']}.vcf"
            gatk_run = subprocess.run(
                [java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                 "-O", str(gatk_vcf)],
                capture_output=True, text=True, timeout=900)
            native_run = subprocess.run(
                [str(native), *common, "--threads", str(arguments.threads),
                 "-O", str(native_vcf)],
                capture_output=True, text=True, timeout=900)
            if gatk_run.returncode != 0 or native_run.returncode != 0:
                run_errors.append(
                    f"{case['case']}: java_exit={gatk_run.returncode} "
                    f"native_exit={native_run.returncode} "
                    f"java_stderr={gatk_run.stderr[-300:]!r} "
                    f"native_stderr={native_run.stderr[-300:]!r}")
                continue
            gatk_rows = records(gatk_vcf)
            native_rows = records(native_vcf)
            comparison = compare_rows(gatk_rows, native_rows)
            comparison["case"] = case["case"]
            comparison["why"] = case["why"]
            comparison["interval"] = case["interval"]
            comparison["extra_args"] = list(case["extra"])
            comparison["gated"] = case["case"] in gated
            comparison.setdefault("phasing_notes", [])
            comparison["gatk_rows"] = len(gatk_rows)
            comparison["native_rows"] = len(native_rows)
            comparison["gatk_qual"] = qual_summary(gatk_rows)
            comparison["native_qual"] = qual_summary(native_rows)
            comparison["gatk_rows_all"] = gatk_rows
            comparison["native_rows_all"] = native_rows
            results.append(comparison)

    violations: list[str] = list(run_errors)
    out_of_scope_notes: list[str] = []
    for result in results:
        if result["gated"]:
            violations += [f"{result['case']}: {item}"
                           for item in result["violations"]]
            violations += [f"{result['case']}: raw QUAL {item}"
                           for item in result["raw_qual_mismatches"]]
        elif not result["data_rows_byte_identical"]:
            out_of_scope_notes.append(
                f"{result['case']}: raw QUAL mismatches "
                f"{result['raw_qual_mismatches'][:5]}; field diffs "
                f"{result['violations'][:5]}")

    print(f"# {Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# reference-confidence (-ERC) pileups with and without a symbolic "
          "spanning deletion, so `*` and <NON_REF> enter the AF pseudocounts")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}; threads={arguments.threads}")
    print(f"# gated cases (mismatch => exit 1): {[case['case'] for case in CASES]}")
    print(f"# out-of-scope probes (reported, not gated): "
          f"{[case['case'] for case in OUT_OF_SCOPE_CASES]}")
    for result in results:
        print(f"[{result['case']}] interval={result['interval']} "
              f"extra={result['extra_args']} gated={result['gated']}")
        print(f"    why: {result['why']}")
        print(f"    gatk_rows={result['gatk_rows']} native_rows={result['native_rows']} "
              f"data_rows_byte_identical={result['data_rows_byte_identical']} "
              f"allele_order_normalizations={len(result['allele_order_normalizations'])} "
              f"raw_qual_mismatches={len(result['raw_qual_mismatches'])}")
        for item in result["allele_order_normalizations"]:
            print(f"    normalized: {item}")
        for item in result["raw_qual_mismatches"]:
            print(f"    RAW QUAL MISMATCH: {item}")
        for item in result["phasing_notes"]:
            print(f"    PHASING NOTE (reported, not gated): {item}")
        for row in result["gatk_rows_all"]:
            print(f"    GATK   {chr(9).join(row)}")
        for row in result["native_rows_all"]:
            print(f"    NATIVE {chr(9).join(row)}")
        if result["first_divergent_row_field_diff"]:
            print("    first divergent row field diff (GATK -> NATIVE): "
                  + result["first_divergent_row_field_diff"])
    if run_errors:
        print(f"# {len(run_errors)} run error(s):")
        for error in run_errors:
            print(f"#   - {error}")
    if violations:
        print(f"# {len(violations)} violation(s):")
        for violation in violations:
            print(f"#   - {violation}")
    if out_of_scope_notes:
        print(f"# {len(out_of_scope_notes)} out-of-scope note(s) (not counted as "
              f"violations):")
        for note in out_of_scope_notes:
            print(f"#   - {note}")

    payload = {
        "status": ("diagnostic" if not strict_mode else
                   ("fail" if violations else "pass")),
        "release": "GATK 4.6.2.0",
        "oracle": "HaplotypeCaller gVCF/<NON_REF> symbolic-allele AF prior parity",
        "binary": str(native),
        "acceptance_criterion": (
            "for every gated case the data rows (header/provenance ignored) must "
            "be byte-identical to pinned GATK 4.6.2.0 after the documented "
            "allele-order canonicalization, and QUAL must be identical on the "
            "raw rows"),
        "strict_mode": strict_mode,
        "cases": [case["case"] for case in selected],
        "gated_cases": [case["case"] for case in CASES],
        "out_of_scope_cases": [case["case"] for case in OUT_OF_SCOPE_CASES],
        "violations": violations,
        "out_of_scope_notes": out_of_scope_notes,
        "results": results,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
