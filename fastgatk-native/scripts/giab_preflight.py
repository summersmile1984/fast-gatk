#!/usr/bin/env python3
"""Disk-backed GIAB biological preflight; execution is supplied by the CLI.

No caller or evaluator runs here. A successful tool exit alone is never an
eligibility decision. Large transient SAM/depth files live under data-root and
are removed only after their complete audit succeeds; failed evidence remains.
"""
from __future__ import annotations

import collections
import fcntl
import gzip
import hashlib
import json
import re
import shutil
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path

from giab_assets import file_digest, load_lock

UINT32_MAX = (1 << 32) - 1
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
SCOPES = {"chr20": ("chr20",), "autosomes": AUTOSOMES}
CIGAR = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")
QNAME = re.compile(r"[!-?A-~]{1,254}\Z")
SEQUENCE = re.compile(r"[=ACMGRSVTWYHKDBNacmgrsvtwyhkdbn]+\Z")
QUALITY = re.compile(r"[!-~]+\Z")
CONTIG_ATTRIBUTE = re.compile(r'(?:^|,)([^=,]+)=("(?:[^"\\]|\\.)*"|[^,]*)')
REMOVED_TAGS = {"DT", "DI", "DS"}
FLAG_NAMES = {0x4: "unmapped", 0x100: "secondary", 0x200: "qc_fail",
              0x400: "duplicate", 0x800: "supplementary"}


class GateFailure(RuntimeError):
    def __init__(self, message: str, **evidence):
        super().__init__(message)
        self.evidence = evidence


def _require(condition, message: str, **evidence) -> None:
    if not condition:
        raise GateFailure(message, **evidence)


def _json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _failure(exc: Exception) -> dict:
    return {"status": "blocked", "reason": str(exc),
            "error_type": type(exc).__name__,
            "evidence": getattr(exc, "evidence", {})}


def _sam_header(lines: list[str]) -> dict:
    result = {kind: [] for kind in ("HD", "SQ", "RG", "PG", "CO")}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        kind = parts[0][1:]
        if kind == "CO":
            result["CO"].append("\t".join(parts[1:]))
            continue
        fields = {}
        for field in parts[1:]:
            key, sep, value = field.partition(":")
            _require(sep and key not in fields, "Malformed SAM header field", line=line)
            fields[key] = value
        result.setdefault(kind, []).append(fields)
    _require(len(result["HD"]) == 1 and result["HD"][0].get("SO") == "coordinate",
             "BAM must declare coordinate sort order")
    sq = result["SQ"]
    _require(sq and len({s.get("SN") for s in sq}) == len(sq),
             "Missing or duplicate BAM SQ records")
    for entry in sq:
        _require(entry.get("SN") and int(entry.get("LN", 0)) > 0,
                 "Invalid BAM sequence dictionary", entry=entry)
    rgs = result["RG"]
    _require(rgs and len({r.get("ID") for r in rgs}) == len(rgs),
             "Missing or duplicate read groups")
    for rg in rgs:
        _require(all(rg.get(key) for key in ("ID", "SM", "LB", "PL")),
                 "RG requires ID/SM/LB/PL", read_group=rg)
        _require(rg["SM"] == "HG002" and rg["PL"].upper() == "ILLUMINA",
                 "Input must be single-sample HG002 Illumina", read_group=rg)
    pg = result["PG"]
    _require(len({p.get("ID") for p in pg}) == len(pg) and all(p.get("ID") for p in pg),
             "Missing or duplicate PG identifiers")
    programs = {p["ID"]: p for p in pg}
    for program in pg:
        visited = set()
        while program.get("PP"):
            previous = program["PP"]
            _require(previous in programs and previous not in visited,
                     "Broken or cyclic PG ancestry", previous=previous)
            visited.add(previous)
            program = programs[previous]
    return result


def _payload(fields: list[str], tags: dict, rg: dict, cigar: list) -> dict:
    """HtsReader byte spans have no NUL; CIGAR offsets count uint32 ops."""
    length = 0 if fields[9] == "*" else len(fields[9])
    values = {"records": 1, "bases": length, "qualities": length,
              "cigar_ops": len(cigar), "names": len(fields[0]),
              "read_groups": len(tags.get("RG", ("", ""))[1]),
              "samples": len(rg.get("SM", "")), "libraries": len(rg.get("LB", "")),
              "platforms": len(rg.get("PL", "")), "platform_units": len(rg.get("PU", "")),
              "flow_orders": len(rg.get("FO", ""))}
    for tag, name in (("OA", "original_alignments"), ("XM", "mate_contigs"),
                      ("BI", "insertion_qualities"), ("BD", "deletion_qualities"),
                      ("t0", "flow_t0_phred")):
        kind, value = tags.get(tag, ("", ""))
        values[name] = len(value) if kind == "Z" else 0
    kind, value = tags.get("tp", ("", ""))
    values["flow_tp"] = value.count(",") if kind == "B" else 0
    return values


def _audit_sam(path: Path, audit: Path, name: str, *, cleared: bool = False) -> dict:
    header_lines = []
    counts = collections.Counter()
    per_contig = {}
    errors = collections.Counter()
    examples = []
    ordered = hashlib.sha256()
    contig_digests = {}
    last_key = (-1, -1)
    header = None
    sq_order = {}
    lengths = {}
    read_groups = {}

    def error(kind, line_number, **details):
        errors[kind] += 1
        if len(examples) < 100:
            examples.append({"kind": kind, "sam_line": line_number, **details})

    with path.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            if header is None and line.startswith("@"):
                header_lines.append(line)
                continue
            if header is None:
                header = _sam_header(header_lines)
                sq_order = {s["SN"]: i for i, s in enumerate(header["SQ"])}
                lengths = {s["SN"]: int(s["LN"]) for s in header["SQ"]}
                read_groups = {r["ID"]: r for r in header["RG"]}
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                error("short_record", line_number)
                continue
            try:
                flag, pos, mapq, mate_pos, tlen = map(int, (fields[1], fields[3], fields[4], fields[7], fields[8]))
            except ValueError:
                error("invalid_integer", line_number, qname=fields[0])
                continue
            counts["records"] += 1
            chrom = fields[2]
            metrics = per_contig.get(chrom)
            if metrics is None:
                metrics = {"counts": collections.Counter(), "aggregate_payload": collections.Counter()}
                per_contig[chrom] = metrics
            stats = metrics["counts"]
            stats["records"] += 1
            for bit, label in FLAG_NAMES.items():
                if flag & bit:
                    counts[label] += 1
                    stats[label] += 1
            base_count = 0 if fields[9] == "*" else len(fields[9])
            counts["bases"] += base_count
            stats["bases"] += base_count
            if not QNAME.fullmatch(fields[0]) or fields[0] == "*":
                error("invalid_qname", line_number, qname=fields[0])
            if not (0 <= flag <= 65535 and 0 <= mapq <= 255 and pos >= 0 and mate_pos >= 0
                    and -(1 << 31) <= tlen < (1 << 31)):
                error("invalid_core_fields", line_number, qname=fields[0])
            if chrom == "*":
                key = (len(sq_order), 0)
                if pos != 0 or not flag & 0x4:
                    error("invalid_unplaced_record", line_number, qname=fields[0])
            elif chrom in sq_order:
                key = (sq_order[chrom], pos)
                if not 1 <= pos <= lengths[chrom]:
                    error("position_out_of_bounds", line_number, chrom=chrom, pos=pos)
            else:
                error("unknown_contig", line_number, chrom=chrom)
                key = last_key
            if key < last_key:
                error("not_coordinate_sorted", line_number, chrom=chrom, pos=pos)
            last_key = key
            mate_chrom = chrom if fields[6] == "=" else fields[6]
            if mate_chrom != "*" and (mate_chrom not in lengths or
                    not 0 <= mate_pos <= lengths.get(mate_chrom, -1)):
                error("invalid_mate_coordinate", line_number, qname=fields[0])
            if fields[9] == "*":
                counts["missing_sequence"] += 1
                if not flag & (0x100 | 0x800):
                    error("missing_primary_sequence", line_number, qname=fields[0])
            elif not SEQUENCE.fullmatch(fields[9]):
                error("invalid_sequence", line_number, qname=fields[0])
            if fields[10] == "*":
                counts["missing_quality"] += 1
                if fields[9] != "*":
                    error("missing_quality_with_sequence", line_number, qname=fields[0])
            elif (fields[9] == "*" or len(fields[9]) != len(fields[10]) or
                  not QUALITY.fullmatch(fields[10])):
                error("invalid_quality", line_number, qname=fields[0])
            cigar = CIGAR.findall(fields[5]) if fields[5] != "*" else []
            if fields[5] != "*" and "".join(n + op for n, op in cigar) != fields[5]:
                error("invalid_cigar", line_number, qname=fields[0])
            if not flag & 0x4 and not cigar:
                error("mapped_record_without_cigar", line_number, qname=fields[0])
            query_span = sum(int(n) for n, op in cigar if op in "MIS=X")
            reference_span = sum(int(n) for n, op in cigar if op in "MDN=X")
            if cigar and fields[9] != "*" and query_span != len(fields[9]):
                error("cigar_sequence_length", line_number, qname=fields[0])
            if chrom in lengths and cigar and pos + max(1, reference_span) - 1 > lengths[chrom]:
                error("alignment_out_of_bounds", line_number, qname=fields[0])
            tags = {}
            for field in fields[11:]:
                components = field.split(":", 2)
                if len(components) != 3 or len(components[0]) != 2 or components[0] in tags:
                    error("invalid_or_duplicate_tag", line_number, qname=fields[0], tag=field)
                    continue
                tags[components[0]] = (components[1], components[2])
            rg_tag = tags.get("RG", ("", ""))
            if rg_tag[0] != "Z" or rg_tag[1] not in read_groups:
                error("missing_or_unknown_rg", line_number, qname=fields[0], rg=rg_tag[1])
            if "OQ" in tags:
                counts["oq_present"] += 1
                kind, value = tags["OQ"]
                if kind != "Z" or fields[9] == "*" or len(value) != len(fields[9]) or not QUALITY.fullmatch(value):
                    error("invalid_oq", line_number, qname=fields[0])
                else:
                    counts["oq_valid"] += 1
                    if value != fields[10]:
                        counts["oq_differs_from_qual"] += 1
            for tag in ("BI", "BD"):
                if tag in tags:
                    kind, value = tags[tag]
                    if kind != "Z" or len(value) != len(fields[9]) or not QUALITY.fullmatch(value):
                        error("invalid_indel_quality", line_number, qname=fields[0], tag=tag)
            if cleared and (flag & 0x400 or REMOVED_TAGS.intersection(tags)):
                error("duplicate_state_not_cleared", line_number, qname=fields[0])
            if chrom in lengths and pos > 0:
                metrics["aggregate_payload"].update(_payload(fields, tags, read_groups.get(rg_tag[1], {}), cigar))
            canonical = fields[:11]
            canonical[1] = str(flag & ~0x400)
            canonical.extend(sorted(field for field in fields[11:] if field[:2] not in REMOVED_TAGS))
            encoded = ("\t".join(canonical) + "\n").encode("ascii")
            ordered.update(encoded)
            if chrom not in contig_digests:
                contig_digests[chrom] = hashlib.sha256()
            contig_digests[chrom].update(encoded)
    if header is None:
        header = _sam_header(header_lines)
    header_path = audit / f"{name}.header.sam"
    header_path.write_text("".join(header_lines), encoding="ascii")
    result = {"status": "eligible" if not errors and counts["records"] else "blocked",
              "header": header, "header_path": str(header_path), "counts": dict(counts),
              "per_contig": per_contig, "errors": dict(errors), "error_examples": examples,
              "canonical_ordered_sha256": ordered.hexdigest(),
              "contig_canonical_sha256": {k: v.hexdigest() for k, v in contig_digests.items()},
              "canonicalization": "SAM fields, flag 0x400 cleared, DT/DI/DS omitted, tags sorted; all other fields retained. Ordered equality proves multiset conservation and preserves multiplicity.",
              "full_decode_sam": str(path)}
    _json(audit / f"{name}.audit.json", result)
    return result


def _reference(raw: Path, prepared: Path, run) -> tuple[dict, dict]:
    reference = prepared / "reference.fa"
    partial = prepared / "reference.fa.part"
    with gzip.open(raw, "rb") as source, partial.open("wb") as target:
        shutil.copyfileobj(source, target, length=8 << 20)
    partial.replace(reference)
    run("samtools", ["faidx", str(reference)], "reference-faidx")
    dictionary = prepared / "reference.dict"
    run("samtools", ["dict", "-o", str(dictionary), str(reference)], "reference-dict")
    sequences = {}
    with dictionary.open() as stream:
        for line in stream:
            if not line.startswith("@SQ\t"):
                continue
            fields = dict(item.split(":", 1) for item in line.rstrip().split("\t")[1:])
            name = fields["SN"]
            _require(name not in sequences and re.fullmatch(r"[a-fA-F0-9]{32}", fields.get("M5", "")),
                     "Reference dictionary lacks unique SN/M5", entry=fields)
            sequences[name] = {"length": int(fields["LN"]), "md5": fields["M5"].lower()}
    fai = Path(str(reference) + ".fai")
    indexed = {}
    with fai.open() as stream:
        for line in stream:
            fields = line.rstrip().split("\t")
            _require(len(fields) >= 5 and fields[0] not in indexed, "Malformed reference FAI", line=line)
            indexed[fields[0]] = int(fields[1])
    _require(list(indexed.items()) == [(n, s["length"]) for n, s in sequences.items()] and sequences,
             "Reference FAI and M5 dictionary disagree")
    return sequences, {"status": "eligible", "path": str(reference), **file_digest(reference),
                       "fai": {"path": str(fai), **file_digest(fai)},
                       "dictionary": {"path": str(dictionary), **file_digest(dictionary)},
                       "sequences": sequences}


def _bam_reference(audit: dict, reference: dict) -> dict:
    actual = [(s["SN"], int(s["LN"])) for s in audit["header"]["SQ"]]
    expected = [(n, s["length"]) for n, s in reference.items()]
    _require(actual == expected, "Full BAM/reference SN/LN/order mismatch", bam=actual, reference=expected)
    checked, missing = [], []
    for seq in audit["header"]["SQ"]:
        if "M5" in seq:
            _require(seq["M5"].lower() == reference[seq["SN"]]["md5"],
                     "BAM/reference M5 mismatch", bam_sequence=seq,
                     reference_sequence=reference[seq["SN"]])
            checked.append(seq["SN"])
        else:
            missing.append(seq["SN"])
    return {"status": "eligible", "sn_ln_order": "checked_match", "m5_checked": checked,
            "m5_absent": missing,
            "evidence_boundary": "Absent BAM M5 cannot establish every historical alignment reference base; pinned BAM/reference source provenance is required in addition to SN/LN."}


def _provenance(audit: dict, lock: dict) -> dict:
    declared = lock.get("provenance", {}).get("bam", {})
    problems = []
    if not declared.get("evidence"):
        problems.append("No public input-transformation provenance evidence is pinned")
    if declared.get("duplicates_removed") is not False:
        problems.append("Public provenance does not establish duplicates were retained")
    if declared.get("bqsr_applied") is not False or declared.get("quality_processing") != "unrecalibrated":
        problems.append("QUAL processing eligibility is unknown: header absence of ApplyBQSR is not proof of original qualities; review and pin public preparation provenance")
    programs = audit["header"]["PG"]
    saw_bwa = saw_markdup = False
    for program in programs:
        text = " ".join(program.values())
        lower = text.lower()
        if any(word in lower for word in ("applybqsr", "printreads", "bqsr", "recalibrat", "abra", "calmd -e")):
            problems.append("Potential quality/alignment transformation in PG: " + program["ID"])
        if re.search(r"(?:remove_duplicates|remove_sequencing_duplicates)(?:=|\s+)true", lower):
            problems.append("PG declares duplicate removal: " + program["ID"])
        identity = (program.get("PN") or program["ID"]).lower()
        if identity.startswith("bwa"):
            saw_bwa = True
        elif identity.startswith("markduplicates"):
            saw_markdup = True
            for option in ("remove_duplicates", "remove_sequencing_duplicates"):
                if not re.search(option + r"(?:=|\s+)false", lower):
                    problems.append("MarkDuplicates PG does not explicitly retain reads: " + option)
        elif identity.startswith("samtools"):
            command = program.get("CL", "")
            if not re.search(r"\bsamtools\s+(?:sort|view|merge|index)\b", command):
                problems.append("Unreviewed samtools transformation: " + program["ID"])
            elif re.search(r"\bsamtools\s+view\b", command) and re.search(r"(?:^|\s)-(?:f|F|q|s|e|L|r|R|d|D|N|x|B)(?:\s|[0-9])", command):
                problems.append("samtools view may have filtered or modified original reads: " + program["ID"])
        else:
            problems.append("Unreviewed input program: " + program["ID"])
    if not saw_bwa or not saw_markdup:
        problems.append("Expected BWA and MarkDuplicates provenance is incomplete")
    if audit["counts"].get("oq_differs_from_qual", 0):
        problems.append("OQ differs from QUAL; unexplained quality transformation requires review, not OQ restoration")
    return {"status": "blocked" if problems else "eligible", "quality_eligibility": "unknown" if problems else "eligible",
            "declared_provenance": declared, "reasons": problems,
            "policy": "Do not infer original QUAL from a filename or missing ApplyBQSR PG; do not restore OQ."}


def _check_bam_index(bam: Path, index: Path, audit: dict, reference: dict, work: Path, run, name: str) -> dict:
    bed = work / "reference.bed"
    if not bed.exists():
        bed.write_text("".join(f"{n}\t0\t{s['length']}\n" for n, s in reference.items()))
    output = work / f"{name}.index-query.count"
    run("samtools", ["view", "-c", "-M", "-L", str(bed), "-X", str(bam), str(index)],
        f"{name}-index-query", stdout=output)
    actual = int(output.read_text().strip())
    expected = sum(v["counts"]["records"] for n, v in audit["per_contig"].items() if n != "*")
    _require(actual == expected, "Indexed BAM traversal count differs from full sequential decode",
             expected=expected, actual=actual, count_path=str(output))
    stats_file = work / f"{name}.idxstats.tsv"
    run("samtools", ["idxstats", str(bam)], f"{name}-idxstats", stdout=stats_file)
    seen = set()
    with stats_file.open() as stream:
        for line in stream:
            chrom, length, mapped, unmapped = line.rstrip().split("\t")
            _require(chrom not in seen, "Duplicate idxstats contig", chrom=chrom)
            seen.add(chrom)
            stats = audit["per_contig"].get(chrom, {}).get("counts", {})
            _require(int(mapped) + int(unmapped) == stats.get("records", 0) and
                     int(unmapped) == stats.get("unmapped", 0),
                     "BAM index per-contig counts disagree", chrom=chrom, index_line=line, decoded=stats)
            _require(int(length) == (0 if chrom == "*" else reference[chrom]["length"]),
                     "BAM index dictionary differs", chrom=chrom)
    _require(seen == set(reference) | {"*"}, "Incomplete BAM idxstats dictionary")
    return {"status": "eligible", "indexed_placed_records": actual,
            "query_count": str(output), "idxstats": str(stats_file), "index": str(index),
            "index_digest": file_digest(index)}


def _vcf_header(stream) -> tuple[list[str], dict, list[str]]:
    metadata, contigs = [], {}
    for line in stream:
        if line.startswith("##"):
            metadata.append(line)
            if line.startswith("##contig=<"):
                attributes = {m.group(1): m.group(2).strip('"')
                              for m in CONTIG_ATTRIBUTE.finditer(line.strip()[10:-1])}
                name = attributes.get("ID")
                _require(name and name not in contigs, "Invalid or duplicate VCF contig", line=line)
                length = attributes.get("length", attributes.get("Length"))
                contigs[name] = int(length) if length is not None else None
            continue
        _require(line.startswith("#CHROM\t"), "VCF is missing column header", line=line)
        columns = line.rstrip("\n").split("\t")
        _require(columns[:8] == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
                 "Invalid VCF columns", columns=columns)
        return metadata, contigs, columns
    raise GateFailure("Empty VCF or missing #CHROM")


def _vcf_rows(stream, contigs: dict):
    order = {name: i for i, name in enumerate(contigs)}
    previous = (-1, -1)
    for number, line in enumerate(stream, 1):
        fields = line.rstrip("\n").split("\t")
        _require(len(fields) >= 8 and fields[0] in order, "Malformed VCF record or absent dictionary entry",
                 record=number, line=line[:2000])
        position = int(fields[1])
        key = (order[fields[0]], position)
        _require(position > 0 and key >= previous, "VCF is not coordinate sorted", record=number, line=line[:2000])
        previous = key
        _require(re.fullmatch(r"[ACGTNacgtn]+", fields[3]), "Invalid VCF REF", record=number, line=line[:2000])
        length = contigs[fields[0]]
        end = position + len(fields[3]) - 1
        for item in fields[7].split(";"):
            if item.startswith("END="):
                end = max(end, int(item[4:]))
        _require(length is None or end <= length, "VCF record exceeds declared contig", record=number, line=line[:2000])
        yield fields, line


def _known_sites(raw: Path, output: Path, reference: dict, fasta: Path, work: Path, run, name: str) -> dict:
    plain = work / f"{name}.compatible.vcf"
    removed, retained = collections.Counter(), collections.Counter()
    with gzip.open(raw, "rt", encoding="ascii") as source, plain.open("w", encoding="ascii") as target:
        metadata, contigs, columns = _vcf_header(source)
        for chrom, length in contigs.items():
            if chrom in reference:
                _require(length == reference[chrom]["length"], "Known-sites contig LN mismatch or missing LN",
                         resource=name, chrom=chrom, observed=length, expected=reference[chrom]["length"])
        for line in metadata:
            if not line.startswith("##contig=<"):
                target.write(line)
        for chrom, seq in reference.items():
            target.write(f"##contig=<ID={chrom},length={seq['length']}>\n")
        target.write("\t".join(columns) + "\n")
        for fields, line in _vcf_rows(source, contigs):
            _require(len(fields) == len(columns), "Known-sites VCF field count mismatch", resource=name, line=line[:2000])
            if fields[0] in reference:
                retained[fields[0]] += 1
                target.write(line)
            else:
                removed[fields[0]] += 1
    summary = {"status": "pending", "raw": str(raw), "path": str(output),
               "retained_records_by_contig": dict(retained), "removed_records_by_contig": dict(removed),
               "removed_dictionary_contigs": [c for c in contigs if c not in reference],
               "selection_rule": "Contig membership of caller reference only, independent of truth and calls"}
    _json(work / f"{name}.selection.json", summary)
    _require(sum(retained.values()) > 0, "Known-sites contains no compatible records", resource=name)
    sort_tmp = work / f"{name}.sort"
    run("bcftools", ["sort", "-T", str(sort_tmp), "-Oz", "-o", str(output), str(plain)], f"{name}-sort")
    run("bcftools", ["index", "-f", "-t", str(output)], f"{name}-index")
    direct = _validate_vcf(output, reference, sample=None, exact_dictionary=True)
    # In bcftools 1.21, -N bypasses the normalization path and its REF check.
    # Discard normalized output; the locked/derived input is never rewritten.
    run("bcftools", ["norm", "-f", str(fasta), "-c", "e", "-Ou", "-o", "/dev/null", str(output)], f"{name}-ref-check")
    index_count = work / f"{name}.index-count.txt"
    run("bcftools", ["index", "-n", str(output)], f"{name}-index-count", stdout=index_count)
    _require(int(index_count.read_text().strip()) == sum(retained.values()),
             "Known-sites index count mismatch", resource=name)
    # This indexed query exercises TBI offsets as well as index metadata.
    queried = work / f"{name}.indexed.vcf.gz"
    run("bcftools", ["view", "--no-version", "-r", ",".join(reference), "-Oz", "-o", str(queried), str(output)], f"{name}-indexed-query")
    checked = _validate_vcf(queried, reference, sample=None, exact_dictionary=True)
    _require(checked["records_sha256"] == direct["records_sha256"] and checked["records"] == direct["records"],
             "Known-sites indexed traversal differs from sequential decode", resource=name)
    summary.update(status="eligible", digest=file_digest(output), index_digest=file_digest(Path(str(output) + ".tbi")),
                   records=direct["records"], records_sha256=direct["records_sha256"],
                   selection_evidence=str(work / f"{name}.selection.json"))
    _json(work / f"{name}.selection.json", summary)
    plain.unlink()
    queried.unlink()
    return summary


def _validate_vcf(path: Path, reference: dict, *, sample: str | None, exact_dictionary: bool = False,
                  selected: tuple[str, ...] | None = None) -> dict:
    counts = collections.Counter()
    hashes = {}
    digest = hashlib.sha256()
    with gzip.open(path, "rt", encoding="ascii") as stream:
        _, contigs, columns = _vcf_header(stream)
        if sample:
            _require(columns[8:] == ["FORMAT", sample], "Truth must have exactly sample HG002", columns=columns)
        if exact_dictionary:
            _require(list(contigs.items()) == [(n, s["length"]) for n, s in reference.items()],
                     "Derived known-sites dictionary is not exactly the caller reference")
        for fields, line in _vcf_rows(stream, contigs):
            _require(len(fields) == len(columns), "VCF record/header column count mismatch", line=line[:2000])
            if selected is not None:
                _require(fields[0] in selected, "Scope VCF leaked an out-of-scope record", line=line[:2000])
            if sample:
                _require("GT" in fields[8].split(":"), "Truth record lacks GT", line=line[:2000])
                gt_index = fields[8].split(":").index("GT")
                values = fields[9].split(":")
                _require(gt_index < len(values), "Truth record lacks genotype value", line=line[:2000])
                alleles = re.split(r"[/|]", values[gt_index])
                _require(alleles and all(a == "." or (a.isdigit() and int(a) <= len(fields[4].split(","))) for a in alleles),
                         "Truth genotype allele is invalid", line=line[:2000])
            counts[fields[0]] += 1
            # bcftools may render harmless float/INFO values differently. Alleles,
            # genotype and all other VCF fields are nevertheless required intact.
            encoded = line.encode("ascii")
            digest.update(encoded)
            if fields[0] not in hashes:
                hashes[fields[0]] = hashlib.sha256()
            hashes[fields[0]].update(encoded)
    return {"records": sum(counts.values()), "counts": dict(counts), "contigs": contigs,
            "records_sha256": digest.hexdigest(),
            "contig_sha256": {c: h.hexdigest() for c, h in hashes.items()}, "columns": columns}


def _benchmark_bed(raw: Path, reference: dict, scopes: tuple[str, ...], prepared: Path, work: Path) -> dict:
    database = work / "benchmark.sqlite"
    order = {chrom: i for i, chrom in enumerate(reference)}
    result = {}
    with sqlite3.connect(database) as db:
        db.execute("PRAGMA temp_store=FILE")
        db.execute("PRAGMA cache_size=-8192")
        # Maintain the sort B-tree during insertion: no external sort can spill
        # SQLite's temporary files onto an implicit /tmp filesystem.
        db.execute("CREATE TABLE intervals (rank INTEGER, chrom TEXT, start INTEGER, end INTEGER, ordinal INTEGER, PRIMARY KEY(rank,start,end,ordinal)) WITHOUT ROWID")
        count = 0
        with raw.open() as stream:
            for number, line in enumerate(stream, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.split()
                _require(len(fields) >= 3, "Malformed benchmark BED", line_number=number, line=line)
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                _require(chrom in reference and 0 <= start < end <= reference.get(chrom, {}).get("length", 0),
                         "Benchmark BED outside caller reference", line_number=number, line=line)
                db.execute("INSERT INTO intervals VALUES (?,?,?,?,?)", (order[chrom], chrom, start, end, count))
                count += 1
        _require(count > 0, "Benchmark BED has no intervals")
        db.commit()
        for scope in scopes:
            directory = prepared / scope
            directory.mkdir(exist_ok=True)
            bed = directory / "scope.bed"
            evaluation = directory / "evaluation.bed"
            chromosomes = SCOPES[scope]
            for chrom in chromosomes:
                _require(chrom in reference, "Reference lacks requested primary contig", scope=scope, chrom=chrom)
            bed.write_text("".join(f"{chrom}\t0\t{reference[chrom]['length']}\n" for chrom in chromosomes))
            bases = intervals = 0
            with evaluation.open("w") as target:
                for chrom in chromosomes:
                    begin = end = None
                    for start, stop in db.execute("SELECT start,end FROM intervals WHERE rank=? ORDER BY start,end", (order[chrom],)):
                        if end is not None and start > end:
                            target.write(f"{chrom}\t{begin}\t{end}\n")
                            bases += end - begin
                            intervals += 1
                            begin = end = None
                        begin = start if begin is None else begin
                        end = stop if end is None else max(end, stop)
                    if end is not None:
                        target.write(f"{chrom}\t{begin}\t{end}\n")
                        bases += end - begin
                        intervals += 1
            _require(bases > 0, "Scope has no benchmark bases", scope=scope)
            result[scope] = {"scope_bed": str(bed), "evaluation_bed": str(evaluation),
                             "scope_bases": sum(reference[c]["length"] for c in chromosomes),
                             "confident_bases": bases, "merged_intervals": intervals,
                             "scope_digest": file_digest(bed), "evaluation_digest": file_digest(evaluation),
                             "domain_rule": "Whole declared primary chromosomes intersected with original BED; union counts, no coverage-dependent clipping"}
    return result


def _truth_scope(raw: Path, full: dict, scope: str, reference: dict, prepared: Path, work: Path, run) -> dict:
    selected = SCOPES[scope]
    for chrom in selected:
        _require(full["contigs"].get(chrom) == reference[chrom]["length"],
                 "Truth/caller-reference primary contig LN mismatch or missing LN", scope=scope, chrom=chrom,
                 truth_length=full["contigs"].get(chrom), reference_length=reference[chrom]["length"])
    directory = prepared / scope
    truth = directory / "truth.vcf.gz"
    run("bcftools", ["view", "--no-version", "-R", str(directory / "scope.bed"), "--regions-overlap", "1",
                     "-Oz", "-o", str(truth), str(raw)], f"{scope}-truth-select")
    checked = _validate_vcf(truth, reference, sample="HG002", selected=selected)
    _require(checked["counts"] == {c: full["counts"][c] for c in selected if c in full["counts"]},
             "Truth indexed scope selection lost or duplicated records", scope=scope,
             expected={c: full["counts"].get(c, 0) for c in selected}, actual=checked["counts"])
    # Compare against another sequential bcftools rendering, not the raw lexical
    # VCF: BCF encoding may canonicalize numerics without changing a record.
    sequential = work / f"{scope}.truth-sequential.vcf.gz"
    run("bcftools", ["view", "--no-version", "-T", str(directory / "scope.bed"), "--targets-overlap", "1",
                     "-Oz", "-o", str(sequential), str(raw)], f"{scope}-truth-sequential")
    sequential_audit = _validate_vcf(sequential, reference, sample="HG002", selected=selected)
    _require(checked["contig_sha256"] == sequential_audit["contig_sha256"],
             "Truth indexed selection differs from sequential full-scope alleles/genotypes", scope=scope)
    run("bcftools", ["index", "-f", "-t", str(truth)], f"{scope}-truth-index")
    indexed = work / f"{scope}.truth-indexed.vcf.gz"
    run("bcftools", ["view", "--no-version", "-R", str(directory / "scope.bed"), "--regions-overlap", "1",
                     "-Oz", "-o", str(indexed), str(truth)], f"{scope}-truth-derived-index-query")
    index_audit = _validate_vcf(indexed, reference, sample="HG002", selected=selected)
    _require(index_audit["records_sha256"] == checked["records_sha256"] and
             index_audit["records"] == checked["records"], "Derived truth index loses or changes records", scope=scope)
    run("bcftools", ["norm", "-f", str(prepared / "reference.fa"), "-c", "e",
                     "-Ou", "-o", "/dev/null", str(truth)], f"{scope}-truth-ref-check")
    sequential.unlink()
    indexed.unlink()
    return {"status": "eligible", "reference_eligibility": "eligible", "path": str(truth),
            "digest": file_digest(truth), "index_digest": file_digest(Path(str(truth) + ".tbi")),
            "records": checked["records"], "records_by_contig": checked["counts"],
            "selection": "Whole scope, indexed record overlap=1, all alleles retained; not prefiltered by confident BED",
            "ref_check": "bcftools norm -f reference.fa -c e -Ou -o /dev/null; input unchanged, no REF repair",
            "scope_limit": "A failure conservatively blocks this scope, including mismatches outside confident BED; it is not a claim that all no-alt callsets are ineligible."}


def _coverage(bam: Path, reference: dict, scopes: tuple[str, ...], work: Path, run) -> dict:
    chromosomes = tuple(c for c in reference if any(c in SCOPES[s] for s in scopes))
    bed = work / "coverage-domain.bed"
    bed.write_text("".join(f"{c}\t0\t{reference[c]['length']}\n" for c in chromosomes))
    depth = work / "coverage.depth.tsv"
    run("samtools", ["depth", "-aa", "-b", str(bed), "-G", "0xF04", "-q", "0", "-Q", "0", str(bam)],
        "coverage-depth", stdout=depth)
    stats = {c: {"positions": 0, "depth_sum": 0, "covered_bases": 0, "below_10x_bases": 0} for c in chromosomes}
    digest = hashlib.sha256()
    with depth.open() as stream:
        for number, line in enumerate(stream, 1):
            digest.update(line.encode("ascii"))
            chrom, position, coverage = line.rstrip().split("\t")
            position, coverage = int(position), int(coverage)
            _require(chrom in stats and coverage >= 0, "Invalid depth row", line_number=number, line=line)
            current = stats[chrom]
            _require(position == current["positions"] + 1 and position <= reference[chrom]["length"],
                     "Depth output has duplicated, missing or out-of-range positions", line_number=number, line=line)
            current["positions"] += 1
            current["depth_sum"] += coverage
            current["covered_bases"] += coverage > 0
            current["below_10x_bases"] += coverage < 10
    for chrom, current in stats.items():
        _require(current["positions"] == reference[chrom]["length"], "Incomplete full-domain depth output", chrom=chrom, observed=current)
    result = {"status": "eligible", "per_contig": stats, "scopes": {}, "depth_sha256": digest.hexdigest(),
              "definition": "samtools depth -aa, base/map Q>=0, exclude flags 0xF04; deletions and reference skips not counted; overlapping mates not collapsed; common BAM has old duplicate flags cleared",
              "low_coverage_threshold": 10, "bed_unchanged": True}
    for scope in scopes:
        total = {key: sum(stats[c][key] for c in SCOPES[scope]) for key in next(iter(stats.values()))}
        bases = total["positions"]
        result["scopes"][scope] = {**total, "mean_depth": total["depth_sum"] / bases,
                                     "coverage_fraction": total["covered_bases"] / bases,
                                     "low_coverage_fraction": total["below_10x_bases"] / bases}
    _json(work / "coverage.json", result)
    depth.unlink()
    return result


def _capacity(audit: dict, scopes: tuple[str, ...]) -> dict:
    source_root = Path(__file__).resolve().parents[2]
    evidence = {}
    for relative in ("fastgatk-core/src/hts_reader.cpp", "fastgatk-core/include/fastgatk/io/hts_reader.hpp",
                     "fastgatk-native/src/hc_call.cpp"):
        evidence[relative] = file_digest(source_root / relative)
    result = {"status": "eligible", "uint32_max": UINT32_MAX, "source": evidence, "scopes": {},
              "representation": "HtsReader::next emits one byte/base/quality and non-NUL metadata strings; cigar_ops offsets count uint32 elements, not bytes; HC append_batch checks each accumulated offset and record ordinal.",
              "selection": "All coordinate-placed records in HC whole-chromosome -L scopes, before caller read filters; preprocessing preserves sequence/name/RG/mapping and this payload footprint.",
              "not_a_memory_guarantee": True}
    for scope in scopes:
        payload = collections.Counter()
        for chrom in SCOPES[scope]:
            payload.update(audit["per_contig"].get(chrom, {}).get("aggregate_payload", {}))
        exceeded = {key: value for key, value in payload.items()
                    if value > UINT32_MAX + (1 if key == "records" else 0)}
        result["scopes"][scope] = {"status": "not_run_capacity_limit" if exceeded else "within_uint32_capacity",
                                     "payload_elements": dict(payload), "exceeded": exceeded,
                                     "memory_status": "not_measured"}
    return result


def _prepare_inputs(data_root: Path, run_tool, *, scopes, reference_only: bool) -> dict:
    """Prepare shared inputs and persist fail-closed, scope-specific evidence.

    ``run_tool(tool, args, name, stdout=Path|None)`` must run the pinned tool,
    preserve stderr/argv evidence and raise on any failure. No retries or cache
    reuse occur here. An explicit new invocation makes a new audit directory.
    """
    data_root = Path(data_root).resolve()
    data_root.mkdir(parents=True, exist_ok=True)
    scopes = tuple(scopes)
    _require(scopes and len(set(scopes)) == len(scopes) and all(s in SCOPES for s in scopes),
             "Only nonempty, unique chr20/autosomes scopes are supported", scopes=scopes)
    prepared = data_root / "prepared"
    prepared.mkdir(exist_ok=True)
    # flock is released even on exceptions; a leftover pathname is not a stale lock.
    with (prepared / ".preflight.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        attempt = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
        work = prepared / "audit" / attempt
        work.mkdir(parents=True)
        manifest_path = data_root / ("reference-input-manifest.json" if reference_only else "input-manifest.json")
        manifest = {"schema_version": 1, "sample": "HG002", "reference": "GRCh38_no_alt",
                    "attempt": attempt, "input_status": "pending", "reference_eligibility": "pending",
                    "preflight_kind": "reference_inputs" if reference_only else "complete_inputs",
                    "audit_directory": str(work), "data_root": str(data_root), "gates": {}, "tools": [],
                    "scopes": {s: {"input_status": "pending", "reference_eligibility": "pending"} for s in scopes},
                    "accuracy_status": "not_evaluated", "query_ref_check": "required_after_calling_not_yet_evaluated"}
        if manifest_path.exists():
            shutil.copyfile(manifest_path, work / "previous-input-manifest.json")

        def save():
            _json(manifest_path, manifest)
            _json(work / "input-manifest.json", manifest)

        def run(tool, args, name, *, stdout=None):
            entry = {"tool": tool, "argv": args, "name": f"biological-{attempt}-{name}",
                     "stdout": str(stdout) if stdout else None, "status": "pending"}
            manifest["tools"].append(entry)
            save()
            try:
                result = run_tool(tool, args, entry["name"], stdout=stdout)
                entry.update(status="completed", execution=result)
            except Exception as exc:
                entry.update(_failure(exc))
                save()
                raise
            save()
            return result

        def gate(name, operation):
            manifest["gates"][name] = {"status": "pending"}
            save()
            try:
                result = operation()
                manifest["gates"][name] = result
            except Exception as exc:
                result = _failure(exc)
                manifest["gates"][name] = result
                _json(work / f"{name}.failure.json", result)
            save()
            return result

        save()
        try:
            lock = load_lock()
            _require(lock.get("sample") == "HG002" and lock.get("reference") == "GRCh38_no_alt", "Unexpected asset lock cohort/reference")
            _json(work / "assets-lock.json", lock)
            raw = {}
            verified = {}
            required = ("reference", "truth", "truth_index", "benchmark_bed",
                        "dbsnp", "dbsnp_index", "mills", "mills_index")

            def verify_raw(name):
                spec = lock["assets"][name]
                path = data_root / "raw" / spec["filename"]
                _require(path.is_file(), "Required raw asset is missing; fetch before biological preflight", asset=name, path=str(path))
                digest = file_digest(path)
                _require(digest["md5"] == spec["md5"] and
                         ("size_bytes" not in spec or digest["bytes"] == spec["size_bytes"]),
                         "Raw asset fails pinned checksum/length", asset=name, path=str(path), observed=digest)
                raw[name] = path
                verified[name] = {"path": str(path), **digest}
            for name in required:
                verify_raw(name)
            manifest["raw_assets"] = verified
            sequences, ref = _reference(raw["reference"], prepared, run)
            manifest["gates"]["reference"] = ref
            manifest["reference_path"] = str(prepared / "reference.fa")
            save()
            bed_result = gate("benchmark_bed", lambda: {"status": "eligible", "scopes": _benchmark_bed(raw["benchmark_bed"], sequences, scopes, prepared, work)})
            truth_result = gate("truth_full_decode", lambda: {"status": "eligible", **_validate_vcf(raw["truth"], sequences, sample="HG002")})
            for scope in scopes:
                scope_info = manifest["scopes"][scope]
                if bed_result["status"] == "eligible" and truth_result["status"] == "eligible":
                    scope_info.update(bed_result["scopes"][scope])
                    truth = gate(f"truth_{scope}", lambda scope=scope: _truth_scope(raw["truth"], truth_result, scope, sequences, prepared, work, run))
                    scope_info["truth"] = truth
                    scope_info["reference_eligibility"] = truth.get("reference_eligibility", "blocked")
                else:
                    scope_info["reference_eligibility"] = "blocked"
                    scope_info["truth"] = {"status": "blocked", "reason": "Full truth/BED prerequisite failed"}
            for resource in ("dbsnp", "mills"):
                result = gate(resource, lambda resource=resource: _known_sites(raw[resource], prepared / f"{resource}.vcf.gz",
                              sequences, prepared / "reference.fa", work, run, resource))
                manifest.setdefault("known_sites", {})[resource] = result.get("path")
            prerequisites = ("reference", "dbsnp", "mills", "benchmark_bed", "truth_full_decode") + tuple(f"truth_{s}" for s in scopes)
            if reference_only or any(manifest["gates"].get(g, {}).get("status") != "eligible" for g in prerequisites):
                # The finally block publishes incomplete/blocked input status.
                # Do not inspect, hash or decode the 46 GB BAM until REF gates pass.
                return manifest
            for name in ("bam", "bam_index"):
                verify_raw(name)

            def audit_source():
                run("samtools", ["quickcheck", "-v", str(raw["bam"])], "source-quickcheck")
                decoded = work / "source.sam"
                run("samtools", ["view", "--no-PG", "-h", str(raw["bam"])], "source-decode", stdout=decoded)
                result = _audit_sam(decoded, work, "source")
                _require(result["status"] == "eligible", "Full source SAM audit failed", audit_path=str(work / "source.audit.json"), errors=result["errors"])
                manifest["gates"]["bam_reference"] = _bam_reference(result, sequences)
                manifest["gates"]["provenance"] = _provenance(result, lock)
                result["index_check"] = _check_bam_index(raw["bam"], raw["bam_index"], result, sequences, work, run, "source")
                decoded.unlink()
                result["full_decode_sam_removed_after_success"] = True
                return result

            source = gate("source_bam", audit_source)
            common = None
            if source["status"] == "eligible":
                manifest["gates"]["aggregate_capacity"] = _capacity(source, scopes)
            if source["status"] == "eligible" and manifest["gates"].get("provenance", {}).get("status") == "eligible":
                def prepare_common():
                    target = prepared / "common.bam"
                    run("samtools", ["view", "--no-PG", "-b", "--remove-flags", "0x400", "-x", "DT", "-x", "DI", "-x", "DS",
                                     "-o", str(target), str(raw["bam"])], "common-create")
                    run("samtools", ["quickcheck", "-v", str(target)], "common-quickcheck")
                    run("samtools", ["index", str(target)], "common-index")
                    decoded = work / "common.sam"
                    run("samtools", ["view", "--no-PG", "-h", str(target)], "common-decode", stdout=decoded)
                    result = _audit_sam(decoded, work, "common", cleared=True)
                    _require(result["status"] == "eligible" and result["counts"]["records"] == source["counts"]["records"] and
                             result["canonical_ordered_sha256"] == source["canonical_ordered_sha256"],
                             "Common BAM violates declared record conservation", source_audit=str(work / "source.audit.json"), common_audit=str(work / "common.audit.json"))
                    _require(result["header"] == source["header"], "Common BAM changed SQ/RG/historical PG/header")
                    result["index_check"] = _check_bam_index(target, Path(str(target) + ".bai"), result, sequences, work, run, "common")
                    result.update(path=str(target), digest=file_digest(target))
                    decoded.unlink()
                    result["full_decode_sam_removed_after_success"] = True
                    return result

                common = gate("common_bam", prepare_common)
                if common["status"] == "eligible":
                    manifest["common_bam"] = common["path"]
                    if "chr20" in scopes:
                        def prepare_pilot():
                            pilot = prepared / "chr20.bam"
                            run("samtools", ["view", "--no-PG", "-b", "-o", str(pilot), common["path"], "chr20"], "chr20-create")
                            run("samtools", ["quickcheck", "-v", str(pilot)], "chr20-quickcheck")
                            run("samtools", ["index", str(pilot)], "chr20-index")
                            decoded = work / "chr20.sam"
                            run("samtools", ["view", "--no-PG", "-h", str(pilot)], "chr20-decode", stdout=decoded)
                            result = _audit_sam(decoded, work, "chr20", cleared=True)
                            _require(result["status"] == "eligible" and result["header"] == common["header"] and
                                     result["canonical_ordered_sha256"] == common["contig_canonical_sha256"].get("chr20") and
                                     set(result["per_contig"]) == {"chr20"},
                                     "Pilot did not preserve exactly chr20 records and full SQ/RG/history")
                            result["index_check"] = _check_bam_index(pilot, Path(str(pilot) + ".bai"), result, sequences, work, run, "chr20")
                            result.update(path=str(pilot), digest=file_digest(pilot),
                                          interpretation="Local chr20 preprocessing experiment; interchromosomal mates are absent, not WGS MarkDuplicates/BQSR state")
                            decoded.unlink()
                            result["full_decode_sam_removed_after_success"] = True
                            return result
                        pilot = gate("chr20_bam", prepare_pilot)
                        manifest["scopes"]["chr20"]["bam"] = pilot.get("path")
                    if "autosomes" in scopes:
                        manifest["scopes"]["autosomes"]["bam"] = common["path"]
                    gate("coverage", lambda: _coverage(Path(common["path"]), sequences, scopes, work, run))
        except Exception as exc:
            manifest["gates"]["preparation"] = _failure(exc)
            _json(work / "preparation.failure.json", manifest["gates"]["preparation"])
        finally:
            required_common = ("reference", "dbsnp", "mills", "benchmark_bed", "truth_full_decode")
            if not reference_only:
                required_common += ("source_bam", "bam_reference", "provenance", "common_bam", "coverage", "aggregate_capacity")
            for scope in scopes:
                info = manifest["scopes"][scope]
                required_gates = required_common + (f"truth_{scope}",)
                if not reference_only and scope == "chr20":
                    required_gates += ("chr20_bam",)
                blocked = [name for name in required_gates if manifest["gates"].get(name, {}).get("status") != "eligible"]
                if "preparation" in manifest["gates"]:
                    blocked.append("preparation")
                info["blocking_gates"] = blocked
                info["input_status"] = "blocked" if blocked else ("pending" if reference_only else "eligible")
                if info["reference_eligibility"] == "pending" or (not reference_only and manifest["gates"].get("bam_reference", {}).get("status") != "eligible"):
                    info["reference_eligibility"] = "blocked"
            manifest["input_status"] = "eligible" if all(s["input_status"] == "eligible" for s in manifest["scopes"].values()) else "blocked"
            manifest["reference_eligibility"] = "eligible" if all(s["reference_eligibility"] == "eligible" for s in manifest["scopes"].values()) else "blocked"
            if reference_only:
                manifest["reference_inputs_status"] = "blocked" if any(s["blocking_gates"] for s in manifest["scopes"].values()) else "eligible"
                manifest["input_status"] = "blocked" if manifest["reference_inputs_status"] == "blocked" else "pending"
                manifest["remaining_input_proofs"] = ["Full BAM decode/index and reference dictionary", "QUAL/duplicate-retention provenance", "Common BAM conservation", "chr20 materialization", "Coverage and aggregate capacity"]
            manifest["completed_at"] = datetime.now(timezone.utc).isoformat()
            save()
        return manifest


def prepare_inputs(data_root: Path, run_tool, *, scopes=("chr20", "autosomes")) -> dict:
    """Prepare all shared biological inputs, returning/persisting input-manifest.json."""
    return _prepare_inputs(data_root, run_tool, scopes=scopes, reference_only=False)


def prepare_reference_inputs(data_root: Path, run_tool, *, scopes=("chr20", "autosomes")) -> dict:
    """Run identical reference/truth/known-sites gates before BAM acquisition.

    Publishes reference-input-manifest.json. ``reference_inputs_status`` is the
    early gate; ``input_status`` remains pending until full BAM proof is done.
    No proof from this narrower operation is reused without revalidation.
    """
    return _prepare_inputs(data_root, run_tool, scopes=scopes, reference_only=True)
