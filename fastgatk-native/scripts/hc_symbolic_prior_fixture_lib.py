#!/usr/bin/env python3
"""Shared fixture + comparison machinery for the symbolic-allele (`*`) prior fixtures.

This module exists because the seven sites of `.diag/round-length-prior-sweep.md`
section C.2 are all of one class -- GATK selects an allele's prior/pseudocount
from the allele's *length* and *symbolic-ness*, native used a fixed value -- and
they all need the same two things to become observable:

1. a synthetic pileup in which the symbolic spanning deletion ``*`` enters an
   allele-frequency matrix **without** being the called allele (this is the only
   regime in which the ``*`` pseudocount is visible: once ``*`` is called its
   effective allele count is data-dominated and the prior cancels), and
2. a comparison that can canonicalize native's symbolic-ALT ordering
   (``hc_call.cpp`` renders ``*`` as the last concrete ALT, GATK places it by
   haplotype order) identically on both files for arbitrary ploidy.

Fixture recipe (identical reference to the two previous rounds' oracles so the
results are comparable): a 1500 bp ``chr1`` built with ``random.seed(11)``, a
10 x ``A`` homopolymer at 600, 8 x ``CAG`` at 900, and base ``chr1:698`` = ``C``
inside the ``698-700 CCC`` run.  Reads are 300M placed on a stride; the optional
second group carries ``chr1:698 C>A`` and the optional deletion group carries a
one-base deletion of base 700, encoded either as a real ``138M1D161M`` CIGAR or
as sequence removal behind a ``299M`` CIGAR.

Genotype enumeration follows htsjdk's ``GenotypeIndexCalculator`` (the same
closed form as ``unrank_genotype`` in ``fastgatk-kernels/src/genotype.cpp``), so
Number=G fields can be permuted for any ploidy, not just 2.
"""
from __future__ import annotations

import math
import random
import subprocess
from pathlib import Path

REFERENCE_LENGTH = 1500
HOMOPOLYMER_START = 600          # 1-based, 10 x 'A' -> 600..609
TANDEM_START = 900               # 1-based, 8 x 'CAG' -> 900..923
SNP_POSITION = 698               # 1-based; reference base is 'C' (seed 11)
DEFAULT_INTERVAL = "chr1:500-780"

BP_RESOLUTION = ("-ERC", "BP_RESOLUTION")
ORDINARY_VCF: tuple[str, ...] = ()


# --------------------------------------------------------------------------- #
# fixture construction
# --------------------------------------------------------------------------- #
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


def layout(reference_reads: int = 30, snp_reads: int = 12, deletion_reads: int = 0,
           mode: str = "cigar", snp_alt: str = "A") -> dict:
    """The standard span-del fixture layout (see the module docstring)."""
    result = {"ref": {"count": reference_reads, "start": 440, "step": 3, "length": 300}}
    result["snp"] = ({"count": snp_reads, "start": 440, "step": 3, "length": 300,
                      "alt": snp_alt} if snp_reads else None)
    result["del"] = ({"count": deletion_reads, "start": 560, "step": 4, "length": 300,
                      "mode": mode, "drop_base": 700} if deletion_reads else None)
    return result


def write_sam(work: Path, name: str, case_layout: dict, reference_text: str) -> Path:
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             f"@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}", "@RG\tID:rg\tSM:S1"]
    reads = build_reads(case_layout, reference_text)
    for index, (start, cigar, sequence, tag) in enumerate(
            sorted(reads, key=lambda item: item[0])):
        lines.append(f"{tag}{index}\t0\tchr1\t{start}\t60\t{cigar}\t*\t0\t0\t"
                     f"{sequence}\t{'I' * len(sequence)}\tRG:Z:rg")
    path = work / f"{name}.sam"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def prepare_reference(work: Path, java: str, gatk: str) -> Path:
    """Write ``ref.fa`` (+ ``.fai`` + ``.dict``) and return the fasta path."""
    reference_text = build_reference()
    reference = work / "ref.fa"
    reference.write_text(">chr1\n" + reference_text + "\n", encoding="utf-8")
    (work / "ref.fa.fai").write_text(
        f"chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{REFERENCE_LENGTH + 1}\n",
        encoding="utf-8")
    subprocess.run([java, "-Xmx1g", "-jar", gatk, "CreateSequenceDictionary",
                    "-R", str(reference), "-O", str(work / "ref.dict"),
                    "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, timeout=300)
    return reference


def make_bam(work: Path, name: str, case_layout: dict, reference_text: str,
             java: str, gatk: str) -> Path:
    sam = write_sam(work, name, case_layout, reference_text)
    bam = work / f"{name}.bam"
    subprocess.run([java, "-Xmx1g", "-jar", gatk, "SortSam",
                    "-I", str(sam), "-O", str(bam), "-SO", "coordinate",
                    "--CREATE_INDEX", "true"],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, timeout=300)
    return bam


# --------------------------------------------------------------------------- #
# VCF reading / canonicalization
# --------------------------------------------------------------------------- #
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


def genotype_count(allele_count: int, ploidy: int) -> int:
    """C(allele_count + ploidy - 1, ploidy) -- htsjdk's genotype cardinality."""
    return math.comb(allele_count + ploidy - 1, ploidy)


def unrank_genotype(index: int, allele_count: int, ploidy: int) -> list[int]:
    """htsjdk/GenotypeIndexCalculator enumeration order.

    Mirrors ``unrank_genotype`` in ``fastgatk-kernels/src/genotype.cpp`` verbatim
    (itself the closed form of htsjdk's genotype index): the list is ordered by
    the largest allele first, then recursively.
    """
    alleles = [0] * ploidy
    available = allele_count
    for position in range(ploidy - 1, -1, -1):
        selected = available - 1
        for candidate in range(available):
            group = genotype_count(candidate + 1, position)
            if index < group:
                selected = candidate
                break
            index -= group
        alleles[position] = selected
        available = selected + 1
    return alleles


def _number_g_permutation(allele_count: int, ploidy: int,
                          remap) -> dict[int, int]:
    order = [tuple(sorted(unrank_genotype(index, allele_count, ploidy)))
             for index in range(genotype_count(allele_count, ploidy))]
    lookup = {key: index for index, key in enumerate(order)}
    return {index: lookup[tuple(sorted(remap(allele) for allele in key))]
            for index, key in enumerate(order)}


def canonical_allele_order(alts: list[str]) -> list[str] | None:
    """Sorted concrete ALTs with a trailing ``<NON_REF>`` pinned last."""
    concrete = [alt for alt in alts if alt != "<NON_REF>"]
    order = sorted(concrete)
    if len(order) != len(set(order)):
        return None
    if "<NON_REF>" in alts:
        order.append("<NON_REF>")
    return order


INFO_NUMBER_A = ("AC", "AF", "MLEAC", "MLEAF")
FORMAT_NUMBER_A: tuple[str, ...] = ()
FORMAT_NUMBER_R = ("AD",)
# `GP` (Number=G genotype posteriors) is emitted by GATK/native when
# --genotype-assignment-method USE_POSTERIOR_PROBABILITIES is used, so it must be
# permuted with the ALT list exactly like `PL`.
FORMAT_NUMBER_G = ("PL", "GP", "PG")
FORMAT_GENOTYPE = ("GT", "PGT")
SAFE_KEYS = set(INFO_NUMBER_A) | set(FORMAT_NUMBER_A) | set(FORMAT_NUMBER_R) \
    | set(FORMAT_NUMBER_G) | set(FORMAT_GENOTYPE) | {
        "DP", "GQ", "MIN_DP", "RGQ", "SB", "PS", "PID", "FT", "FILTER", ".",
        "RAW_MQandDP",
    }
PHASING_KEYS = frozenset({"PGT", "PID", "PS"})


def normalize_row(row: list[str]) -> tuple[list[str], bool]:
    """Canonicalize one record so both files can be compared field by field.

    Sorted concrete ALTs with ``<NON_REF>`` last, every allele-indexed field
    permuted with them (``AC/AF/MLEAC/MLEAF`` Number=A, ``AD`` Number=R, ``PL``
    Number=G for **any** ploidy, ``GT``/``PGT`` indices remapped), and phased
    genotype strings written in ascending allele order.  The canonicalization is
    applied unconditionally and identically to both files.  Raises ``ValueError``
    when the record carries a comma-valued field outside the documented
    whitelist, i.e. when permutation would be unsafe.
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
    position_of = {alt: index for index, alt in enumerate(order)}
    perm = [position_of[alt] for alt in alts]
    alleles = len(alts) + 1

    def remap_index(value: int) -> int:
        return 0 if value == 0 else 1 + perm[value - 1]

    info = split_info(row[7])
    for key, value in info.items():
        if key not in SAFE_KEYS and "," in value:
            raise ValueError(f"non-normalizable INFO field {key}={value}")

    def permute_number_a(text: str) -> str:
        values = text.split(",")
        if len(values) != len(alts):
            raise ValueError(f"Number=A arity mismatch: {text}")
        out = [""] * len(alts)
        for old, value in enumerate(values):
            out[perm[old]] = value
        return ",".join(out)

    def permute_number_r(text: str) -> str:
        values = text.split(",")
        if len(values) != alleles:
            raise ValueError(f"Number=R arity mismatch: {text}")
        out = [values[0]] * alleles
        for old in range(len(alts)):
            out[1 + perm[old]] = values[1 + old]
        return ",".join(out)

    def permute_number_g(text: str) -> str:
        # Values are position-independent (PL integers, GP phred-scaled
        # posteriors as decimals), so permute the entries themselves.
        values = text.split(",")
        ploidy = 1
        while (ploidy <= 8
               and genotype_count(alleles, ploidy) != len(values)):
            ploidy += 1
        if genotype_count(alleles, ploidy) != len(values):
            raise ValueError(f"Number=G arity mismatch for {alleles} alleles: "
                             f"{len(values)} entries")
        mapping = _number_g_permutation(alleles, ploidy, remap_index)
        out = [""] * len(values)
        for old, value in enumerate(values):
            out[mapping[old]] = value
        return ",".join(out)

    def permute_genotype(text: str) -> str:
        separator = "|" if "|" in text else "/" if "/" in text else None
        if separator is None:
            return str(remap_index(int(text)))
        mapped = sorted(remap_index(int(part))
                        for part in text.split(separator) if part != ".")
        return separator.join(str(value) for value in mapped)

    for key in INFO_NUMBER_A:
        if key in info and info[key] not in ("", "."):
            info[key] = permute_number_a(info[key])

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
                out_values.append(permute_number_a(value))
            elif key in FORMAT_NUMBER_R:
                out_values.append(permute_number_r(value))
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


def field_differences(gatk_row: list[str], native_row: list[str],
                      exclude_format: frozenset = frozenset()) -> tuple[str, str]:
    """Return ``(kept_differences, excluded_differences)`` as two strings."""
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
        # Align FORMAT fields by *name*, not by position: one side may carry
        # extra keys (GATK emits GP/PG under USE_POSTERIOR_PROBABILITIES, native
        # does not), and positional alignment would then compare unrelated
        # values and hide the real per-field differences.
        gatk_keys = gatk_row[8].split(":")
        native_keys = native_row[8].split(":")
        if gatk_keys != native_keys:
            differences.append(f"FORMAT keys {gatk_row[8]!r} -> {native_row[8]!r}")
        gatk_values = gatk_row[9].split(":")
        native_values = native_row[9].split(":")
        gatk_by_key = dict(zip(gatk_keys, gatk_values))
        native_by_key = dict(zip(native_keys, native_values))
        for key in sorted(set(gatk_keys) | set(native_keys)):
            gatk_value = gatk_by_key.get(key)
            native_value = native_by_key.get(key)
            if gatk_value != native_value:
                text = f"FORMAT {key} {gatk_value!r} -> {native_value!r}"
                if key in exclude_format:
                    excluded.append(text)
                else:
                    differences.append(text)
    return ("; ".join(differences), "; ".join(excluded))


def difference_field(item: str) -> str:
    """Field name of one ``field_differences`` item, e.g. ``INFO QD '8.20' ...``
    -> ``QD``, ``FORMAT PL ...`` -> ``PL``, ``QUAL ...`` -> ``QUAL``."""
    parts = item.split(" ")
    if parts[0] in ("INFO", "FORMAT") and len(parts) > 1:
        return parts[1]
    return parts[0]


def compare_rows(gatk_rows: list[list[str]], native_rows: list[list[str]],
                 divergence_positions: frozenset = frozenset(),
                 tolerated: dict[str, frozenset] | None = None,
                 allow_divergence: bool = False,
                 divergence_fields: frozenset = frozenset({"QUAL", "QD"})) -> dict:
    """Field-level comparison of two VCF data-row lists.

    ``divergence_positions`` names the POS values of the records the fixture is
    built to expose.  Those positions are allowed to differ **only** in
    ``divergence_fields`` (``QUAL`` plus the QUAL-derived ``QD = QUAL/DP``
    annotation, which moves with it) and only when ``require_divergence`` is set
    (``allow_divergence``): QUAL is the quantity the defect moves, every other
    field must still agree exactly.  In strict mode a QUAL difference at those
    positions is a violation like any other.  ``tolerated``
    names further (POS -> fields) pairs that are separately tracked divergences:
    they are reported as notes, never counted as violations and never hidden; the
    key ``"*"`` applies to every record of the case.
    """
    tolerated = tolerated or {}
    violations: list[str] = []
    notes: list[str] = []
    normalizations: list[str] = []
    qual_divergences: list[str] = []
    unexpected_qual: list[str] = []
    first_diff = ""
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
        site = gatk_raw[1]
        if gatk_raw == native_raw:
            continue
        if set(gatk_raw[4].split(",")) != set(native_raw[4].split(",")):
            violations.append(
                f"record {position} ({gatk_raw[0]}:{site}): ALT set differs "
                f"{gatk_raw[4]!r} vs {native_raw[4]!r}")
            continue
        try:
            gatk_norm, gatk_applied = normalize_row(gatk_raw)
            native_norm, native_applied = normalize_row(native_raw)
        except ValueError as error:
            violations.append(f"record {position} ({gatk_raw[0]}:{site}): {error}")
            continue
        if gatk_applied or native_applied:
            normalizations.append(
                f"{gatk_raw[0]}:{site} ALT {gatk_raw[4]!r} / {native_raw[4]!r} -> "
                f"{gatk_norm[4]!r}")
        exclude = PHASING_KEYS if gatk_raw[4] != native_raw[4] else frozenset()
        kept, dropped = field_differences(gatk_norm, native_norm, exclude)
        if dropped:
            notes.append(f"{gatk_raw[0]}:{site} (ALT order differs, phasing keys "
                         f"reported not gated): {dropped}")
        if not kept:
            continue
        allowed_here = site in divergence_positions and allow_divergence
        # A "*" entry in `tolerated` applies to every record of the case.
        fields = (tolerated.get(site, frozenset())
                  | tolerated.get("*", frozenset()))
        for item in kept.split("; "):
            field = difference_field(item)
            if field in fields:
                notes.append(f"{gatk_raw[0]}:{site} separately tracked "
                             f"(reported, not gated): {item}")
                continue
            if field in divergence_fields and allowed_here:
                qual_divergences.append(
                    f"{gatk_raw[0]}:{site} {gatk_raw[3]}>{gatk_raw[4]} {item}")
                continue
            if field in divergence_fields:
                unexpected_qual.append(
                    f"{gatk_raw[0]}:{site} {gatk_raw[3]}>{gatk_raw[4]} {item}")
            violations.append(f"record {position} ({gatk_raw[0]}:{site}): {item}")
            if not first_diff:
                first_diff = item
    return {"violations": violations,
            "separately_tracked_notes": notes,
            "allele_order_normalizations": normalizations,
            "qual_divergences": qual_divergences,
            "unexpected_qual_mismatches": unexpected_qual,
            "data_rows_identical": (not violations and not unexpected_qual
                                    and not qual_divergences),
            "first_divergent_row_field_diff": first_diff}


def qual_summary(rows: list[list[str]]) -> list[dict]:
    return [{"POS": row[1], "REF": row[3], "ALT": row[4], "QUAL": row[5]}
            for row in rows if len(row) > 5]
