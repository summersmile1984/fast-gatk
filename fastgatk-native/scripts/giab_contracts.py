#!/usr/bin/env python3
"""Forensic, disk-backed output contracts, independent of biological scoring.

The runner supplies pinned samtools/bcftools. Inputs are never repaired, sorted,
normalized in place, or substituted with another branch's outputs. SAM/VCF
serializations, SQLite canonical multisets, and every difference remain on disk.
"""
from __future__ import annotations

import collections
import hashlib
import json
import re
import sqlite3
import traceback
from contextlib import closing
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path

from giab_assets import file_digest
from giab_preflight import (CIGAR, CONTIG_ATTRIBUTE, QNAME, QUALITY, SEQUENCE, GateFailure,
                            _json, _require, _sam_header, _vcf_header)

STAGES = {"markdup": "MarkDuplicates", "bqsr": "BaseRecalibrator",
          "apply": "ApplyBQSR", "hc": "HaplotypeCaller", "genotype": "GenotypeGVCFs"}
OUTPUTS = {"markdup": "markdup_bam", "bqsr": "recal_table", "apply": "recal_bam",
           "hc": "gvcf", "genotype": "vcf"}
REPORT_KEYS = {
    "Arguments": ("Argument",), "Quantized": ("QualityScore",),
    "RecalTable0": ("ReadGroup", "EventType"),
    "RecalTable1": ("ReadGroup", "QualityScore", "EventType"),
    "RecalTable2": ("ReadGroup", "QualityScore", "CovariateValue", "CovariateName", "EventType"),
}
REPORT_COLUMNS = {
    "Arguments": ("Argument", "Value"),
    "Quantized": ("QualityScore", "Count", "QuantizedScore"),
    "RecalTable0": ("ReadGroup", "EventType", "EmpiricalQuality", "EstimatedQReported", "Observations", "Errors"),
    "RecalTable1": ("ReadGroup", "QualityScore", "EventType", "EmpiricalQuality", "Observations", "Errors"),
    "RecalTable2": ("ReadGroup", "QualityScore", "CovariateValue", "CovariateName", "EventType", "EmpiricalQuality", "Observations", "Errors"),
}
NUMERIC_COLUMNS = {"QualityScore", "Count", "QuantizedScore", "EmpiricalQuality",
                   "EstimatedQReported", "Observations", "Errors"}


def _compact(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _line(stream, value):
    stream.write(_compact(value) + "\n")


def _metadata(path: Path, *, digest=False) -> dict:
    _require(path.is_file(), "Required artifact is absent", path=str(path))
    stat = path.stat()
    result = {"path": str(path.resolve()), "bytes": stat.st_size,
              "mtime_ns": stat.st_mtime_ns, "device": stat.st_dev, "inode": stat.st_ino}
    if digest:
        result.update(file_digest(path))
    return result


def _json_metadata(path: Path):
    """Parse metadata incrementally; HC's potentially huge calls array is hashed.

    The unmodified sidecar remains the authoritative record-level evidence.
    Only one call (or scalar token) is materialized at a time, including in the
    otherwise unbounded aggregate-HC manifest and telemetry JSON.
    Repeated object keys are accepted only when their values agree; conflicting
    duplicates are ambiguous. The source bytes and checksum preserve both.
    """
    decoder = json.JSONDecoder()
    whitespace = re.compile(r"\s*")
    buffer, offset, eof = "", 0, False
    with path.open() as stream:
        def extend():
            nonlocal buffer, offset, eof
            chunk = stream.read(65536)
            buffer, offset = buffer[offset:] + chunk, 0
            eof = not chunk

        def peek():
            nonlocal offset
            while True:
                offset = whitespace.match(buffer, offset).end()
                if offset < len(buffer):
                    return buffer[offset]
                if eof:
                    return ""
                extend()

        def take(expected):
            nonlocal offset
            _require(peek() == expected, "Malformed JSON metadata punctuation", path=str(path), expected=expected)
            offset += 1

        def scalar():
            nonlocal offset
            peek()
            while True:
                try:
                    value, end = decoder.raw_decode(buffer, offset)
                    if not eof and (end == len(buffer) or buffer[end] not in " \t\r\n,:]}"):
                        extend()
                        continue
                    _require(end == len(buffer) or buffer[end] in " \t\r\n,:]}",
                             "Invalid JSON scalar delimiter", path=str(path))
                    _require(not isinstance(value, float) or Decimal(str(value)).is_finite(),
                             "Non-finite JSON metadata number", path=str(path))
                    offset = end
                    return value
                except json.JSONDecodeError:
                    _require(not eof, "Truncated or malformed JSON metadata", path=str(path))
                    extend()

        def parse(key=None):
            char = peek()
            if char == "{":
                take("{")
                result = {}
                if peek() == "}":
                    take("}")
                    return result
                while True:
                    name = scalar()
                    _require(isinstance(name, str), "Invalid JSON metadata key", path=str(path))
                    take(":")
                    value = parse(name)
                    _require(name not in result or result[name] == value,
                             "Conflicting duplicate JSON metadata key", path=str(path), key=name)
                    result[name] = value
                    if peek() == "}":
                        take("}")
                        return result
                    take(",")
            if char == "[":
                take("[")
                projected = key == "calls"
                result, count, digest = [], 0, hashlib.sha256()
                if peek() != "]":
                    while True:
                        value = parse()
                        if projected:
                            digest.update((_compact(value) + "\n").encode())
                        else:
                            result.append(value)
                        count += 1
                        if peek() == "]":
                            break
                        take(",")
                take("]")
                return {"retained_in_source_array": "calls", "records": count, "sha256": digest.hexdigest()} if projected else result
            _require(char != "", "Unexpected end of JSON metadata", path=str(path))
            return scalar()

        result = parse()
        _require(peek() == "", "Trailing JSON metadata content", path=str(path))
        return result


def _database(path: Path):
    _require(not path.exists(), "Refusing to overwrite previous forensic database", path=str(path))
    db = sqlite3.connect(path)
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA cache_size=-8192")
    return db


def _number(text: str) -> str:
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise GateFailure("Invalid numeric report cell", value=text) from exc
    _require(value.is_finite(), "Non-finite numeric report cell", value=text)
    rendered = format(value, "f")
    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")
    return "0" if value == 0 else rendered


def _difference(left, right):
    with localcontext() as context:
        context.prec = max(len(left), len(right), 32) * 2
        return _number(str(Decimal(left) - Decimal(right)))


def _run_callback(tools, result):
    prefix = "contract-" + hashlib.sha256(result["metadata_path"].encode()).hexdigest()[:20]
    def run(tool, argv, name, *, stdout=None):
        entry = {"tool": tool, "argv": list(map(str, argv)), "name": prefix + "-" + name,
                 "stdout_path": str(stdout) if stdout is not None else None, "status": "pending"}
        result.setdefault("commands", []).append(entry)
        try:
            entry["result"] = tools(tool, entry["argv"], entry["name"], stdout=stdout)
            entry["status"] = "completed"
            return entry["result"]
        except Exception as exc:
            entry.update(status="failed", error_type=type(exc).__name__, reason=str(exc))
            raise
    return run


def _reference(reference: Path):
    fai = Path(str(reference) + ".fai")
    entries = {}
    with fai.open() as stream:
        for number, line in enumerate(stream, 1):
            cells = line.rstrip("\n").split("\t")
            _require(len(cells) >= 5 and cells[0] not in entries, "Malformed reference FAI", line=number)
            length, offset, bases, width = map(int, cells[1:5])
            _require(length > 0 and offset >= 0 and 0 < bases <= width,
                     "Invalid reference FAI geometry", line=number)
            entries[cells[0]] = (length, offset, bases, width)
    _require(entries, "Empty reference FAI")
    dictionary = reference.with_suffix(".dict")
    sequences = {}
    with dictionary.open() as stream:
        for line in stream:
            if line.startswith("@SQ\t"):
                fields = dict(field.split(":", 1) for field in line.rstrip("\n").split("\t")[1:])
                _require(fields.get("SN") not in sequences and re.fullmatch(r"[0-9a-fA-F]{32}", fields.get("M5", "")),
                         "Invalid reference dictionary SN/M5", dictionary=str(dictionary))
                sequences[fields["SN"]] = (int(fields["LN"]), fields["M5"].lower())
    _require([(name, value[0]) for name, value in sequences.items()] == [(name, value[0]) for name, value in entries.items()],
             "Reference dictionary and FAI disagree")
    entries = {name: (*value, sequences[name][1]) for name, value in entries.items()}
    return entries


def _index(path: Path, *, bam=False):
    candidates = [Path(str(path) + suffix) for suffix in ((".bai", ".csi") if bam else (".tbi", ".csi"))]
    if bam:
        candidates.append(path.with_suffix(".bai"))
    present = [p for p in candidates if p.is_file() and p.stat().st_size]
    _require(present, "Missing or empty output index", path=str(path), candidates=list(map(str, candidates)))
    return [_metadata(p, digest=True) for p in present]


def _native_bundle(stage: str, output: Path):
    manifest_path = Path(str(output) + ".manifest.json")
    manifest = _json_metadata(manifest_path)
    _require(isinstance(manifest, dict) and manifest.get("schema_version") == 1,
             "Invalid native manifest schema", path=str(manifest_path))
    _require(manifest.get("tool") == STAGES[stage], "Wrong native manifest tool", manifest=manifest)
    implementation = {"markdup": "fastgatk-mark-duplicates", "bqsr": "fastgatk-bqsr",
                      "apply": "fastgatk-bqsr", "hc": "fastgatk-hc-call", "genotype": "fastgatk-genotype-gvcf"}[stage]
    _require(manifest.get("implementation") == implementation, "Unexpected native implementation", path=str(manifest_path))
    _require(Path(manifest.get("primary_output", "")).resolve() == output.resolve(),
             "Manifest refers to a different primary output", path=str(manifest_path))
    declared = manifest.get("outputs")
    _require(isinstance(declared, list) and declared, "Native manifest lacks output bundle")
    by_kind, seen = {}, set()
    for item in declared:
        _require(isinstance(item, dict) and item.get("path") and item.get("kind") and item.get("complete") is True,
                 "Native manifest declares an incomplete artifact", artifact=item)
        path = Path(item["path"]).resolve()
        _require(path.parent == output.resolve().parent and path not in seen,
                 "Manifest artifact is borrowed, duplicated, or outside its branch directory", artifact=item)
        seen.add(path)
        metadata = _metadata(path)
        _require(metadata["bytes"] > 0, "Empty native bundle artifact", artifact=item)
        by_kind.setdefault(item["kind"], []).append(metadata)
    _require(output.resolve() in seen, "Primary output absent from native bundle")
    required = {"markdup": ("duplicate-metrics", "alignment-index"),
                "bqsr": ("bqsr-covariate-table",), "apply": ("bam-or-cram-index",),
                "hc": ("vcf-index",), "genotype": ("vcf-index",)}[stage]
    for kind in required:
        _require(len(by_kind.get(kind, [])) == 1, "Missing or ambiguous mandatory native sidecar", kind=kind)
    if stage == "bqsr":
        _require(Path(by_kind["bqsr-covariate-table"][0]["path"]) == Path(str(output.resolve()) + ".covariates.tsv"),
                 "Native BQSR sidecar path differs from ApplyBQSR's consumed path")
    if stage == "hc":
        _require(manifest.get("sample_name") == "HG002", "Native HC manifest sample differs from HG002")
    return {"metadata": _metadata(manifest_path, digest=True), "content": manifest, "by_kind": by_kind}


def _bam_store(path: Path, directory: Path, name: str, reference: dict, run):
    sam = directory / (name + ".sam")
    database = directory / (name + ".sqlite")
    _require(not database.exists() and not sam.exists(), "Audit evidence directory was already used", directory=str(directory))
    run("samtools", ["quickcheck", "-v", str(path)], name + "-quickcheck")
    run("samtools", ["view", "--no-PG", "-h", str(path)], name + "-decode", stdout=sam)
    counts = collections.Counter()
    per_contig = collections.Counter()
    header_lines, header = [], None
    digest = hashlib.sha256()
    last = (-1, -1)
    with closing(_database(database)) as db:
        db.execute("CREATE TABLE records (identity TEXT, qual TEXT, duplicate INTEGER, n INTEGER, first_line INTEGER, PRIMARY KEY(identity,qual,duplicate)) WITHOUT ROWID")
        db.execute("CREATE INDEX records_apply ON records(identity,duplicate,qual)")
        with sam.open(encoding="ascii") as stream:
            for number, line in enumerate(stream, 1):
                if header is None and line.startswith("@"):
                    header_lines.append(line)
                    continue
                if header is None:
                    header = _sam_header(header_lines)
                    ranks = {sq["SN"]: i for i, sq in enumerate(header["SQ"])}
                    read_groups = {rg["ID"] for rg in header["RG"]}
                fields = line.rstrip("\n").split("\t")
                _require(len(fields) >= 11, "Short SAM record", path=str(sam), line=number)
                flag, position, mapq, mate_position, tlen = map(int, (fields[1], fields[3], fields[4], fields[7], fields[8]))
                _require(0 <= flag <= 65535 and 0 <= mapq <= 255 and position >= 0 and mate_position >= 0,
                         "Invalid SAM integer", path=str(sam), line=number)
                _require(QNAME.fullmatch(fields[0]), "Invalid SAM QNAME", line=number)
                chrom = fields[2]
                _require(chrom == "*" or chrom in reference, "Unknown SAM contig", line=number, chrom=chrom)
                _require(fields[6] in ("*", "=") or fields[6] in reference,
                         "Unknown mate contig", line=number)
                if fields[6] not in ("*", "="):
                    _require(mate_position <= reference[fields[6]][0], "Mate position outside reference", line=number)
                if fields[6] == "=":
                    _require(chrom in reference and mate_position <= reference[chrom][0], "Mate position outside reference", line=number)
                key = (len(ranks), 0) if chrom == "*" else (ranks[chrom], position)
                _require(key >= last, "SAM coordinate order violation", line=number, previous=last, current=key)
                last = key
                sequence, quality = fields[9:11]
                _require(sequence == "*" or SEQUENCE.fullmatch(sequence), "Invalid SAM SEQ", line=number)
                _require(quality == "*" or (QUALITY.fullmatch(quality) and sequence != "*" and len(quality) == len(sequence)),
                         "Illegal or length-mismatched SAM QUAL", line=number)
                cigar = list(CIGAR.finditer(fields[5])) if fields[5] != "*" else []
                _require(fields[5] == "*" or "".join(m.group(0) for m in cigar) == fields[5],
                         "Malformed CIGAR", line=number)
                reference_span = sum(int(m[1]) for m in cigar if m[2] in "MDN=X")
                query_span = sum(int(m[1]) for m in cigar if m[2] in "MIS=X")
                if cigar and sequence != "*":
                    _require(query_span == len(sequence), "CIGAR/SEQ length mismatch", line=number)
                if chrom != "*":
                    _require(position <= reference[chrom][0] and position - 1 + reference_span <= reference[chrom][0],
                             "Alignment exceeds reference", line=number)
                if not flag & 4:
                    _require(chrom != "*" and position > 0 and cigar, "Mapped read lacks mapping", line=number)
                tags = {}
                for raw in fields[11:]:
                    parts = raw.split(":", 2)
                    _require(len(parts) == 3 and parts[0] not in tags, "Malformed/duplicate SAM tag", line=number, tag=raw)
                    tags[parts[0]] = parts[1:]
                _require(tags.get("RG", [None, None])[0] == "Z" and tags["RG"][1] in read_groups,
                         "Record lacks valid read group", line=number)
                for tag in ("OQ", "BI", "BD"):
                    if tag in tags:
                        _require(tags[tag][0] == "Z" and QUALITY.fullmatch(tags[tag][1]) and sequence != "*" and len(tags[tag][1]) == len(sequence),
                                 "Invalid quality tag", line=number, tag=tag)
                counts["records"] += 1
                counts["duplicate_records"] += bool(flag & 0x400)
                counts["missing_qualities"] += quality == "*"
                counts["duplicate_tags"] += any(t in tags for t in ("DT", "DI", "DS"))
                counts["bases"] += 0 if sequence == "*" else len(sequence)
                per_contig[chrom] += 1
                # All mandatory fields and all auxiliary tags except PG survive.
                # Duplication and QUAL remain explicit columns, never discarded.
                core = fields[:10]
                core[1] = str(flag & ~0x400)
                core[3], core[4], core[7], core[8] = map(str, (position, mapq, mate_position, tlen))
                identity = _compact([core, sorted((tag, *value) for tag, value in tags.items() if tag != "PG")])
                db.execute("INSERT INTO records VALUES (?,?,?,?,?) ON CONFLICT(identity,qual,duplicate) DO UPDATE SET n=n+1",
                           (identity, quality, int(bool(flag & 0x400)), 1, number))
                digest.update(line.encode("ascii"))
                if counts["records"] % 10000 == 0:
                    db.commit()
        if header is None:
            header = _sam_header(header_lines)
        _require([(s["SN"], int(s["LN"])) for s in header["SQ"]] == [(n, value[0]) for n, value in reference.items()],
                 "Full BAM/reference dictionary mismatch", path=str(path))
        for sequence in header["SQ"]:
            if "M5" in sequence:
                _require(sequence["M5"].lower() == reference[sequence["SN"]][4],
                         "BAM/reference M5 mismatch", bam_sequence=sequence)
        db.commit()
    result = {"metadata": _metadata(path), "sam": str(sam), "sqlite": str(database), "header": header,
              "counts": dict(counts), "per_contig": dict(per_contig), "records_sha256": digest.hexdigest()}
    _json(directory / (name + ".json"), result)
    return result


def _bam_index(path, summary, directory, name, reference, run):
    indices = _index(path, bam=True)
    bed = directory / (name + ".index-domain.bed")
    with bed.open("w") as stream:
        for chrom, entry in reference.items():
            stream.write(f"{chrom}\t0\t{entry[0]}\n")
    count = directory / (name + ".indexed-count.txt")
    run("samtools", ["view", "--no-PG", "-c", "-M", "-L", str(bed), str(path)], name + "-indexed-query", stdout=count)
    indexed = int(count.read_text().strip())
    expected = sum(value for chrom, value in summary["per_contig"].items() if chrom != "*")
    _require(indexed == expected, "BAM indexed query disagrees with full decode", indexed=indexed, sequential=expected)
    return {"indices": indices, "indexed_placed_records": indexed, "query_evidence": str(count)}


def _merge(left, right):
    """Merge unique ordered (key, payload) iterators, without collecting ties."""
    left, right = iter(left), iter(right)
    a, b = next(left, None), next(right, None)
    while a is not None or b is not None:
        if b is None or (a is not None and a[0] < b[0]):
            yield a[0], a[1], None
            a = next(left, None)
        elif a is None or b[0] < a[0]:
            yield b[0], None, b[1]
            b = next(right, None)
        else:
            yield a[0], a[1], b[1]
            a, b = next(left, None), next(right, None)


def _bam_rows(db, mode):
    columns = {"markdup": "identity,qual", "apply": "identity,duplicate",
               "exact": "identity,qual,duplicate", "duplicates": "identity,qual"}[mode]
    where = " WHERE duplicate=1" if mode == "duplicates" else ""
    for row in db.execute(f"SELECT {columns},SUM(n),MIN(first_line) FROM records{where} GROUP BY {columns} ORDER BY {columns}"):
        yield tuple(row[:-2]), (row[-2], row[-1])


def _bam_difference(left, right, path, mode):
    changed = missing_left = missing_right = 0
    with closing(sqlite3.connect(left["sqlite"])) as a, closing(sqlite3.connect(right["sqlite"])) as b, path.open("w") as stream:
        for key, x, y in _merge(_bam_rows(a, mode), _bam_rows(b, mode)):
            xn, yn = (x[0] if x else 0), (y[0] if y else 0)
            if xn != yn:
                changed += 1
                missing_left += max(0, yn - xn)
                missing_right += max(0, xn - yn)
                _line(stream, {"record_identity": json.loads(key[0]), "projected_fields": key[1:],
                               "left_count": xn, "right_count": yn,
                               "left_sam_line": x[1] if x else None, "right_sam_line": y[1] if y else None})
    return {"different_keys": changed, "left_only_records": missing_right, "right_only_records": missing_left,
            "differences": str(path), "projection": mode,
            "left_sam": left["sam"], "right_sam": right["sam"]}


def _quality_rows(db, keep_duplicate):
    cols = "identity,duplicate,qual" if keep_duplicate else "identity,qual"
    for row in db.execute(f"SELECT {cols},SUM(n),MIN(first_line) FROM records GROUP BY {cols} ORDER BY {cols}"):
        identity = tuple(row[:-3])
        yield (identity, row[-3]), (row[-2], row[-1])


def _quality_difference(left, right, directory, *, keep_duplicate):
    database = directory / "quality-residuals.sqlite"
    evidence = directory / "quality-differences.jsonl"
    counts = collections.Counter()
    histogram = collections.Counter()
    changed_histogram = collections.Counter()
    with closing(_database(database)) as residual, closing(sqlite3.connect(left["sqlite"])) as a, closing(sqlite3.connect(right["sqlite"])) as b:
        residual.execute("CREATE TABLE residual (side INTEGER, identity TEXT, qual TEXT, n INTEGER, line INTEGER, PRIMARY KEY(side,identity,qual)) WITHOUT ROWID")
        for (identity, qual), x, y in _merge(_quality_rows(a, keep_duplicate), _quality_rows(b, keep_duplicate)):
            xn, yn = (x[0] if x else 0), (y[0] if y else 0)
            same = min(xn, yn)
            counts["identical_quality_records"] += same
            if qual != "*":
                histogram[0] += same * len(qual)
                changed_histogram[0] += same
            for side, value, common in ((0, x, yn), (1, y, xn)):
                if value and value[0] > common:
                    residual.execute("INSERT INTO residual VALUES (?,?,?,?,?)", (side, _compact(identity), qual, value[0] - common, value[1]))
        residual.commit()
        # Exact quality strings cancel first. Residual indistinguishable records
        # pair lexicographically, not by arbitrary coordinate-tie iteration order.
        aa = iter(residual.execute("SELECT identity,qual,n,line FROM residual WHERE side=0 ORDER BY identity,qual"))
        bb = iter(residual.execute("SELECT identity,qual,n,line FROM residual WHERE side=1 ORDER BY identity,qual"))
        x, y = next(aa, None), next(bb, None)
        with evidence.open("w") as stream:
            while x is not None or y is not None:
                if y is None or (x is not None and x[0] < y[0]):
                    counts["left_unmatched_records"] += x[2]
                    _line(stream, {"identity": json.loads(x[0]), "left": x[1:], "right": None})
                    x = next(aa, None)
                    continue
                if x is None or y[0] < x[0]:
                    counts["right_unmatched_records"] += y[2]
                    _line(stream, {"identity": json.loads(y[0]), "left": None, "right": y[1:]})
                    y = next(bb, None)
                    continue
                n = min(x[2], y[2])
                counts["different_quality_records"] += n
                delta = collections.Counter()
                changed = None
                if x[1] == "*" or y[1] == "*":
                    counts["missing_quality_pairs"] += n
                else:
                    _require(len(x[1]) == len(y[1]), "Matched reads have different quality lengths", identity=x[0])
                    delta.update(ord(r) - ord(l) for l, r in zip(x[1], y[1]))
                    changed = len(x[1]) - delta[0]
                    changed_histogram[changed] += n
                    for value, amount in delta.items():
                        histogram[value] += amount * n
                _line(stream, {"identity": json.loads(x[0]), "multiplicity": n, "left_qual": x[1], "right_qual": y[1],
                               "left_sam_line": x[3], "right_sam_line": y[3], "changed_bases": changed,
                               "delta_histogram_right_minus_left": dict(delta)})
                x = next(aa, None) if x[2] == n else (x[0], x[1], x[2] - n, x[3])
                y = next(bb, None) if y[2] == n else (y[0], y[1], y[2] - n, y[3])
    return {"counts": dict(counts), "base_delta_histogram_right_minus_left": dict(sorted(histogram.items())),
            "changed_bases_per_read": dict(sorted(changed_histogram.items())), "differences": str(evidence),
            "residuals_sqlite": str(database), "left_sam": left["sam"], "right_sam": right["sam"],
            "pairing": "Exact quality multisets cancel; remaining identical-record copies pair lexicographically. Physical identity of otherwise identical copies is unknowable.",
            "duplicate_flag_in_pairing_key": keep_duplicate}


def _report_db(path):
    db = _database(path)
    db.execute("CREATE TABLE tables (name TEXT PRIMARY KEY, columns_json TEXT, numeric_json TEXT, declared_rows INTEGER) WITHOUT ROWID")
    db.execute("CREATE TABLE cells (name TEXT, key TEXT, values_json TEXT, line INTEGER, PRIMARY KEY(name,key)) WITHOUT ROWID")
    return db


def _insert_report(db, table, keys, values, number):
    _require(all(key in values and values[key] != "" for key in keys), "Empty semantic report key", table=table, line=number)
    key = [_number(values[k]) if k == "QualityScore" else values[k] for k in keys]
    if table == "RecalTable2" and values["CovariateName"] == "Cycle":
        key[keys.index("CovariateValue")] = _number(values["CovariateValue"])
    try:
        db.execute("INSERT INTO cells VALUES (?,?,?,?)", (table, _compact(key), _compact(values), number))
    except sqlite3.IntegrityError as exc:
        raise GateFailure("Duplicate semantic report key", table=table, key=key, line=number) from exc


def _bqsr_report(path, directory, name):
    database = directory / (name + ".sqlite")
    tables = {}
    with closing(_report_db(database)) as db, path.open() as stream:
        lines = enumerate(stream, 1)
        _, first = next(lines)
        match = re.fullmatch(r"#:GATKReport\.v1\.1:([0-9]+)\s*", first)
        _require(match is not None, "Invalid GATK report header", path=str(path))
        declared_tables = int(match[1])
        for number, line in lines:
            if not line.strip():
                continue
            definition = line.rstrip("\n").split(":")
            _require(len(definition) >= 6 and definition[:2] == ["#", "GATKTable"], "Invalid GATK table definition", line=number)
            columns_count, rows_count = int(definition[2]), int(definition[3])
            _require(columns_count > 0 and rows_count >= 0 and len(definition[4:-1]) == columns_count and definition[-1] == ";",
                     "Invalid GATK table dimensions/formats", line=number)
            _, title = next(lines)
            name_parts = title.rstrip("\n").split(":", 3)
            _require(len(name_parts) == 4 and name_parts[:2] == ["#", "GATKTable"], "Malformed GATK table title", line=number + 1)
            table = name_parts[2]
            _require(table in REPORT_KEYS and table not in tables, "Unexpected or duplicate BQSR report table", table=table)
            _, header = next(lines)
            tokens = list(re.finditer(r"\S+", header))
            columns = [token[0] for token in tokens]
            _require(tuple(columns) == REPORT_COLUMNS[table] and len(columns) == columns_count,
                     "BQSR report column schema mismatch", table=table, columns=columns)
            starts = [token.start() for token in tokens]
            numeric = sorted(set(columns) & NUMERIC_COLUMNS)
            db.execute("INSERT INTO tables VALUES (?,?,?,?)", (table, _compact(columns), _compact(numeric), rows_count))
            for _ in range(rows_count):
                row_number, row = next(lines)
                _require(row.strip() and not row.startswith("#"), "Truncated GATK table", table=table, line=row_number)
                if "\t" in header:
                    fields = row.rstrip("\n").split("\t")
                else:
                    text = row.rstrip("\n")
                    fields = [text[start:starts[i + 1] if i + 1 < len(starts) else len(text)].strip()
                              for i, start in enumerate(starts)]
                _require(len(fields) == columns_count, "GATK table width mismatch", table=table, line=row_number)
                values = dict(zip(columns, fields))
                for column in numeric:
                    value = Decimal(_number(values[column]))
                    _require(value >= 0, "Negative BQSR numeric cell", table=table, column=column, line=row_number)
                    if column in ("QualityScore", "QuantizedScore", "Count", "Observations"):
                        _require(value == int(value), "Nonintegral BQSR integer cell", table=table, column=column, line=row_number)
                    if column in ("QualityScore", "QuantizedScore"):
                        _require(value <= 93, "Quality score out of SAM range", table=table, line=row_number)
                if "EventType" in values:
                    _require(values["EventType"] in ("M", "I", "D"), "Unknown BQSR event", line=row_number)
                    _require(Decimal(values["Errors"]) <= Decimal(values["Observations"]), "Errors exceed observations", line=row_number)
                _insert_report(db, table, REPORT_KEYS[table], values, row_number)
            tables[table] = {"rows": rows_count, "columns": columns, "formats": definition[4:-1], "description": name_parts[3]}
            db.commit()
        _require(set(tables) == set(REPORT_KEYS) and len(tables) == declared_tables, "Incomplete BQSR report", tables=tables)
        _require(tables["Quantized"]["rows"] == 94, "Quantization map must cover Q0 through Q93")
    return {"metadata": _metadata(path, digest=True), "sqlite": str(database), "tables": tables}


def _covariates(path, directory):
    database = directory / "covariates.sqlite"
    rows = 0
    with closing(_database(database)) as db, path.open() as stream:
        _require(stream.readline().rstrip("\n") == "# FASTGATK-BQSR-COVARIATES v2", "Unsupported native covariate version", path=str(path))
        expected = "# read_group\tcycle\tcontext\traw_quality\tcount\tmismatches\tempirical_quality\tdelta\tevent_type"
        _require(stream.readline().rstrip("\n") == expected, "Malformed native covariate header")
        db.execute("CREATE TABLE covariates (key TEXT PRIMARY KEY, values_json TEXT, line INTEGER) WITHOUT ROWID")
        for number, line in enumerate(stream, 3):
            fields = line.rstrip("\n").split("\t")
            _require(len(fields) == 9, "Malformed native covariate row", line=number)
            group, cycle, context, quality, count, errors, empirical, delta, event = fields
            cycle, quality, count, errors, delta = map(int, (cycle, quality, count, errors, delta))
            _require(group and context and event in ("M", "I", "D") and 0 <= quality <= 93 and 0 <= errors <= count
                     and Decimal(_number(empirical)) >= 0 and -20 <= delta <= 20,
                     "Invalid native covariate values", line=number)
            key = _compact([group, cycle, context, quality, event])
            try:
                db.execute("INSERT INTO covariates VALUES (?,?,?)", (key, _compact(fields), number))
            except sqlite3.IntegrityError as exc:
                raise GateFailure("Duplicate native covariate key", line=number, key=key) from exc
            rows += 1
            if rows % 10000 == 0:
                db.commit()
        db.commit()
    return {"metadata": _metadata(path, digest=True), "sqlite": str(database), "rows": rows,
            "retained_for_apply": True, "schema_source": "bqsr_tool.cpp:run_base_recalibrator/read_covariates v2"}


def _metrics(path, directory):
    database = directory / "metrics.sqlite"
    sections = {}
    table = None
    columns = None
    pending = None
    histogram_number = 0
    with closing(_report_db(database)) as db, path.open() as stream:
        for number, line in enumerate(stream, 1):
            text = line.rstrip("\n")
            if text.startswith("## METRICS CLASS"):
                _require("DuplicationMetrics" in text, "Unexpected MarkDuplicates metrics class", line=number)
                pending = "metrics"
                columns = None
            elif text.startswith("## HISTOGRAM"):
                histogram_number += 1
                pending = "histogram-" + str(histogram_number)
                columns = None
            elif not text.strip():
                table = None
                columns = None
            elif text.startswith("#"):
                continue
            elif columns is None:
                _require(pending is not None, "Data outside a metrics section", line=number)
                table, pending = pending, None
                columns = text.split("\t")
                _require(len(columns) > 1 and len(set(columns)) == len(columns) and table not in sections,
                         "Invalid metrics header", line=number)
                if table == "metrics":
                    required = {"LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS",
                                "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION"}
                    _require(required <= set(columns) and columns[0] == "LIBRARY", "Incomplete duplicate metrics header")
                numeric = columns[1:] if table == "metrics" else columns
                sections[table] = {"columns": columns, "rows": 0}
                db.execute("INSERT INTO tables VALUES (?,?,?,NULL)", (table, _compact(columns), _compact(numeric)))
            else:
                fields = text.split("\t")
                _require(len(fields) == len(columns), "Metrics row width mismatch", line=number)
                values = dict(zip(columns, fields))
                for column in columns if table != "metrics" else columns[1:]:
                    if values[column] != "":
                        _require(Decimal(_number(values[column])) >= 0, "Negative duplicate metric", column=column, line=number)
                if table == "metrics":
                    _require(Decimal(values["PERCENT_DUPLICATION"]) <= 1, "Invalid duplication fraction", line=number)
                    _require(Decimal(values["READ_PAIR_OPTICAL_DUPLICATES"]) <= Decimal(values["READ_PAIR_DUPLICATES"]),
                             "Optical duplicate count exceeds pair duplicates", line=number)
                else:
                    values[columns[0]] = _number(values[columns[0]])
                _insert_report(db, table, (columns[0],), values, number)
                sections[table]["rows"] += 1
        _require("metrics" in sections and sections["metrics"]["rows"] > 0 and pending is None,
                 "Missing/truncated per-library duplicate metrics")
        db.commit()
    return {"metadata": _metadata(path, digest=True), "sqlite": str(database), "sections": sections}


def _scope(db, bed, reference):
    ranks = {name: rank for rank, name in enumerate(reference)}
    db.execute("CREATE TABLE scope (rank INTEGER, start INTEGER, end INTEGER, PRIMARY KEY(rank,start,end)) WITHOUT ROWID")
    count = 0
    with bed.open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split()
            _require(len(fields) >= 3 and fields[0] in reference, "Invalid scope BED", line=number)
            start, end = map(int, fields[1:3])
            _require(0 <= start < end <= reference[fields[0]][0], "Scope BED outside reference", line=number)
            db.execute("INSERT OR IGNORE INTO scope VALUES (?,?,?)", (ranks[fields[0]], start, end))
            count += 1
    _require(count, "Empty scope BED")


def _interval_union(rows):
    current = None
    for rank, start, end in rows:
        if current is None:
            current = (rank, start, end)
        elif rank == current[0] and start <= current[2]:
            current = (rank, current[1], max(end, current[2]))
        else:
            yield current
            current = (rank, start, end)
    if current is not None:
        yield current


def _non_n(reference_file, entry, start, end):
    """FAI byte arithmetic; no contig-sized FASTA allocation or faidx argument list."""
    _, offset, bases, width = entry[:4]
    non_n = 0
    while start < end:
        stop = min(end, start + (1 << 20))
        begin_byte = offset + start // bases * width + start % bases
        last = stop - 1
        end_byte = offset + last // bases * width + last % bases + 1
        reference_file.seek(begin_byte)
        data = reference_file.read(end_byte - begin_byte).replace(b"\n", b"").replace(b"\r", b"")
        _require(len(data) == stop - start, "FAI disagrees with plain reference FASTA", start=start, end=stop)
        non_n += len(data) - data.count(b"N") - data.count(b"n")
        start = stop
    return non_n


def _coverage_gaps(db, reference, fasta, directory, name):
    evidence = directory / (name + ".coverage-gaps.jsonl")
    contigs = tuple(reference)
    counts = collections.Counter()
    coverage = iter(_interval_union(db.execute("SELECT rank,start,end FROM spans ORDER BY rank,start,end")))
    current = next(coverage, None)
    with fasta.open("rb") as source, evidence.open("w") as target:
        def gap(rank, start, end):
            if start >= end:
                return
            chrom = contigs[rank]
            non_n = _non_n(source, reference[chrom], start, end)
            counts["gap_intervals"] += 1
            counts["uncovered_bases"] += end - start
            counts["uncovered_non_n_bases"] += non_n
            _line(target, {"chrom": chrom, "start": start, "end": end, "non_n_bases": non_n})
        for rank, start, end in _interval_union(db.execute("SELECT rank,start,end FROM scope ORDER BY rank,start,end")):
            counts["scope_bases"] += end - start
            cursor = start
            while current is not None and (current[0] < rank or (current[0] == rank and current[2] <= start)):
                current = next(coverage, None)
            while current is not None and current[0] == rank and current[1] < end:
                gap(rank, cursor, min(end, current[1]))
                cursor = max(cursor, min(end, current[2]))
                if current[2] >= end:
                    break
                current = next(coverage, None)
            gap(rank, cursor, end)
    return {**dict(counts), "uncovered_non_n_bases": counts["uncovered_non_n_bases"], "evidence": str(evidence),
            "definition": "Union of reference blocks and full variant reference spans; only reference N/n may explain uncovered scope bases."}


def _vcf_definitions(metadata):
    definitions = {kind: {} for kind in ("INFO", "FORMAT", "FILTER", "ALT")}
    for line in metadata:
        for kind, entries in definitions.items():
            prefix = "##" + kind + "=<"
            if line.startswith(prefix):
                fields = {match[1]: match[2].strip('"')
                          for match in CONTIG_ATTRIBUTE.finditer(line.rstrip("\n")[len(prefix):-1])}
                identifier = fields.get("ID")
                _require(identifier and identifier not in entries, "Invalid/duplicate VCF header declaration", declaration=line)
                entries[identifier] = {key: value for key, value in fields.items() if key != "Description"}
    return definitions


def _vcf_semantic(fields, definitions):
    # Do not decompose or left-align alleles: retain representation differences.
    # Only order-insensitive serialization and numeric spelling are canonicalized.
    def value(kind, key, raw):
        if raw is not None and definitions[kind].get(key, {}).get("Type") in ("Integer", "Float"):
            return ",".join(_number(part) if part != "." else part for part in raw.split(","))
        return raw
    info = []
    for item in fields[7].split(";") if fields[7] != "." else []:
        key, sep, value_text = item.partition("=")
        info.append((key, value("INFO", key, value_text) if sep else None))
    fmt = fields[8].split(":")
    values = fields[9].split(":")
    sample = dict(zip(fmt, values + ["."] * (len(fmt) - len(values))))
    sample = {key: value("FORMAT", key, raw) for key, raw in sample.items()}
    if "/" in sample["GT"]:
        sample["GT"] = "/".join(sorted(sample["GT"].split("/"), key=lambda a: (-1 if a == "." else int(a))))
    qual = fields[5] if fields[5] == "." else _number(fields[5])
    value = [fields[0], int(fields[1]), fields[2], fields[3], fields[4], qual,
             sorted(fields[6].split(";")), sorted(info), sorted(sample.items())]
    return _compact(value)


def _vcf_audit(path, directory, name, reference_path, reference, scope_bed, run, *, gvcf):
    sequential = directory / (name + ".vcf")
    indexed = directory / (name + ".indexed.vcf")
    database = directory / (name + ".sqlite")
    indices = _index(path)
    # -N bypasses REF validation in bcftools 1.21. Output is deliberately discarded.
    run("bcftools", ["norm", "-f", str(reference_path), "-c", "e", "-Ou", "-o", "/dev/null", str(path)], name + "-ref-check")
    run("bcftools", ["view", "--no-version", "-Ov", "-o", str(sequential), str(path)], name + "-decode")
    counts = collections.Counter()
    digest = hashlib.sha256()
    ranks = {name: rank for rank, name in enumerate(reference)}
    last = (-1, -1)
    last_block = (-1, -1)
    last_variant_anchor = (-1, -1)
    with closing(_database(database)) as db, sequential.open() as stream:
        db.execute("CREATE TABLE records (key TEXT PRIMARY KEY, n INTEGER, first_record INTEGER) WITHOUT ROWID")
        db.execute("CREATE TABLE spans (rank INTEGER, start INTEGER, end INTEGER, PRIMARY KEY(rank,start,end)) WITHOUT ROWID")
        _scope(db, scope_bed, reference)
        metadata, contigs, columns = _vcf_header(stream)
        _require(columns[8:] == ["FORMAT", "HG002"], "VCF must contain exactly HG002", columns=columns)
        _require(list(contigs.items()) == [(n, value[0]) for n, value in reference.items()],
                 "VCF dictionary differs from full reference", observed=contigs)
        _json(directory / (name + ".header.json"), {"metadata": metadata, "columns": columns, "contigs": contigs})
        definitions = _vcf_definitions(metadata)
        for number, line in enumerate(stream, 1):
            fields = line.rstrip("\n").split("\t")
            _require(len(fields) == 10 and fields[0] in ranks, "Invalid VCF row", record=number, path=str(sequential))
            chrom, position, ref, alt = fields[0], int(fields[1]), fields[3], fields[4].split(",")
            rank = ranks[chrom]
            key = (rank, position)
            _require(key >= last and position > 0, "Unsorted VCF", record=number)
            last = key
            _require(re.fullmatch(r"[ACGTNacgtn]+", ref) and len(set(alt)) == len(alt), "Invalid REF or duplicate ALT", record=number)
            _require(all(a == "*" or re.fullmatch(r"[ACGTNacgtn]+|<[A-Za-z0-9_:]+>", a) for a in alt),
                     "Invalid small-variant ALT", record=number, alt=alt)
            _require(all(a.upper() != ref.upper() for a in alt), "ALT equals REF", record=number)
            info = {}
            for entry in fields[7].split(";") if fields[7] != "." else []:
                tag, sep, value = entry.partition("=")
                _require(tag and tag not in info, "Duplicate/invalid INFO tag", record=number)
                info[tag] = value if sep else None
            end = position + len(ref) - 1
            if "END" in info:
                declared_end = int(info["END"])
                _require(declared_end >= end, "END precedes variant reference span", record=number)
                end = declared_end
            _require(end <= reference[chrom][0], "VCF span exceeds reference", record=number)
            _require(db.execute("SELECT 1 FROM scope WHERE rank=? AND start<? AND end>? LIMIT 1",
                                (rank, end, position - 1)).fetchone() is not None,
                     "VCF record lies entirely outside calling scope", record=number, chrom=chrom, position=position)
            if fields[5] != ".":
                _require(Decimal(_number(fields[5])) >= 0, "Negative VCF QUAL", record=number)
            fmt = fields[8].split(":")
            values = fields[9].split(":")
            _require(len(set(fmt)) == len(fmt) and fmt and fmt[0] == "GT" and len(values) <= len(fmt),
                     "Invalid FORMAT or absent first GT", record=number)
            gt = values[0]
            _require(re.fullmatch(r"(?:\.|[0-9]+)(?:[/|](?:\.|[0-9]+))?", gt), "Malformed GT", record=number, gt=gt)
            alleles = re.split(r"[/|]", gt)
            _require((gt == "." or len(alleles) == 2) and all(a == "." or int(a) <= len(alt) for a in alleles),
                     "Non-diploid or out-of-range GT", record=number, gt=gt)
            block = alt == ["<NON_REF>"]
            if gvcf:
                _require("<NON_REF>" in alt, "gVCF record lacks NON_REF likelihood allele", record=number)
                if block:
                    _require(len(ref) == 1 and all(a in ("0", ".") for a in alleles), "Invalid reference-block genotype/REF", record=number)
                    _require(last_block[0] != rank or position > last_block[1], "Overlapping reference blocks", record=number)
                    _require(last_variant_anchor != key, "Reference block overlaps variant anchor", record=number)
                    last_block = (rank, end)
                    counts["reference_blocks"] += 1
                else:
                    _require(last_block[0] != rank or position > last_block[1], "Variant anchor lies inside a reference block", record=number)
                    last_variant_anchor = key
                    counts["variant_records"] += 1
                # A deletion span may overlap later variants or reference blocks
                # on the other haplotype. It is not a row-overlap violation.
                db.execute("INSERT OR IGNORE INTO spans VALUES (?,?,?)", (rank, position - 1, end))
            else:
                _require("<NON_REF>" not in alt and "<*>" not in alt,
                         "Final VCF contains reference-confidence residue", record=number)
                _require(not all(a == "0" for a in alleles), "Final VCF contains a hom-reference record", record=number)
            counts["records"] += 1
            counts["no_call_records"] += all(a == "." for a in alleles)
            counts["filtered_records"] += fields[6] not in ("PASS", ".")
            db.execute("INSERT INTO records VALUES (?,?,?) ON CONFLICT(key) DO UPDATE SET n=n+1", (_vcf_semantic(fields, definitions), 1, number))
            digest.update(line.encode("ascii"))
            if number % 10000 == 0:
                db.commit()
        db.commit()
        coverage = _coverage_gaps(db, reference, reference_path, directory, name) if gvcf else None
    run("bcftools", ["view", "--no-version", "-r", ",".join(reference), "-Ov", "-o", str(indexed), str(path)], name + "-indexed-query")
    indexed_digest = hashlib.sha256()
    indexed_count = 0
    with indexed.open() as stream:
        _vcf_header(stream)
        for line in stream:
            indexed_digest.update(line.encode("ascii"))
            indexed_count += 1
    _require(indexed_count == counts["records"] and indexed_digest.hexdigest() == digest.hexdigest(),
             "Indexed VCF query differs from full decode", sequential=counts["records"], indexed=indexed_count)
    result = {"metadata": _metadata(path), "indices": indices, "sqlite": str(database), "vcf": str(sequential),
              "indexed_vcf": str(indexed), "header": str(directory / (name + ".header.json")),
              "counts": dict(counts), "records_sha256": digest.hexdigest(), "coverage": coverage, "definitions": definitions,
              "ref_check": "bcftools norm -f REF -c e -Ou -o /dev/null; original input unchanged",
              "empty_callset": counts["records"] == 0}
    _json(directory / (name + ".json"), result)
    if coverage:
        _require(coverage["uncovered_non_n_bases"] == 0, "gVCF leaves unexplained non-N scope gaps", **coverage)
    return result


def _header_conservation(left, right):
    # New PG nodes and record PG provenance are allowed; RG/SQ/HD/CO are not
    # silently stripped. Comments and non-SO HD fields are retained as evidence.
    _require(left["header"]["SQ"] == right["header"]["SQ"]
             and {rg["ID"]: rg for rg in left["header"]["RG"]} == {rg["ID"]: rg for rg in right["header"]["RG"]},
             "BAM SQ/RG changed during preprocessing", left=left["header"], right=right["header"])
    old_pg = {p["ID"]: p for p in left["header"]["PG"]}
    new_pg = {p["ID"]: p for p in right["header"]["PG"]}
    _require(all(new_pg.get(key) == value for key, value in old_pg.items()), "Historical BAM PG was lost or modified")


def audit_stage(stage: str, branch: str, *, input_path: Path | None, output_path: Path,
                reference: Path, scope_bed: Path, work: Path, tools) -> dict:
    """Audit one stage; all contract exceptions become persisted failed evidence."""
    directory = Path(work) / stage
    directory.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 1, "stage": stage, "branch": branch, "contract_status": "failed",
              "accuracy_status": "not_evaluated", "metadata_path": str(directory / "audit.json"), "artifacts": {}}
    try:
        _require(stage in STAGES and branch in ("native", "gatk"), "Unsupported contract stage/branch", stage=stage, branch=branch)
        output_path, reference, scope_bed = map(lambda p: Path(p).resolve(), (output_path, reference, scope_bed))
        input_path = Path(input_path).resolve() if input_path is not None else None
        result["consumed"] = {"input": _metadata(input_path) if input_path else None,
                              "output": _metadata(output_path), "reference": _metadata(reference),
                              "reference_fai": _metadata(Path(str(reference) + ".fai"), digest=True),
                              "reference_dictionary": _metadata(reference.with_suffix(".dict"), digest=True),
                              "scope_bed": _metadata(scope_bed, digest=True)}
        _require(input_path != output_path, "Input and output alias the same path")
        if input_path:
            _require(not input_path.samefile(output_path), "Input and output alias the same inode")
        run = _run_callback(tools, result)
        reference_entries = _reference(reference)
        bundle = _native_bundle(stage, output_path) if branch == "native" else None
        result["native_bundle"] = bundle
        if branch == "native" and stage == "hc":
            telemetry_path = output_path.parent / "hc.telemetry.json"
            telemetry = _json_metadata(telemetry_path)
            embedded = bundle["content"].get("telemetry")
            _require(isinstance(telemetry, dict) and isinstance(embedded, dict),
                     "HC telemetry metadata is not an object", path=str(telemetry_path))
            # ResourceSnapshot::budget() resamples filesystem space on each
            # serialization (resource.cpp:402-442), including after file writes.
            sampled = {"scratch_free_bytes", "scratch_hard_bytes"}
            observations = {}
            comparable = []
            for label, document in (("sidecar", telemetry), ("manifest", embedded)):
                resources = document.get("resources")
                _require(isinstance(resources, dict) and
                         all(type(resources.get(key)) is int and resources[key] >= 0 for key in sampled),
                         "Invalid HC scratch-space observations", source=label)
                observations[label] = {key: resources[key] for key in sorted(sampled)}
                comparable.append({**document, "resources": {key: value for key, value in resources.items() if key not in sampled}})
            _require(comparable[0] == comparable[1],
                     "HC telemetry differs beyond independently sampled scratch space", path=str(telemetry_path))
            result["hc_telemetry"] = {"metadata": _metadata(telemetry_path, digest=True), "content": telemetry,
                                      "sampled_resource_values": observations}
        artifacts = result["artifacts"]
        if stage in ("markdup", "apply"):
            _require(input_path is not None, "BAM conservation audit requires its own input")
            before = artifacts["input"] = _bam_store(input_path, directory, "input", reference_entries, run)
            after = artifacts["output"] = _bam_store(output_path, directory, "output", reference_entries, run)
            _header_conservation(before, after)
            after["index"] = _bam_index(output_path, after, directory, "output", reference_entries, run)
            conservation = artifacts["conservation"] = _bam_difference(before, after, directory / "conservation.jsonl", stage)
            if stage == "markdup":
                _require(not after["counts"].get("duplicate_tags"), "DontTag output contains DT/DI/DS tags")
                metrics = Path(bundle["by_kind"]["duplicate-metrics"][0]["path"]) if bundle else output_path.with_suffix(".metrics.txt")
                artifacts["metrics"] = _metrics(metrics, directory)
                artifacts["duplicate_changes"] = _bam_difference(before, after, directory / "duplicate-changes.jsonl", "duplicates")
            else:
                artifacts["quality_changes"] = _quality_difference(before, after, directory, keep_duplicate=True)
                _require(not artifacts["quality_changes"]["counts"].get("missing_quality_pairs"),
                         "ApplyBQSR introduced or removed missing QUAL values")
            _require(conservation["different_keys"] == 0, "Stage changed forbidden SAM fields or record multiplicities", **conservation)
        elif stage == "bqsr":
            artifacts["output"] = _bqsr_report(output_path, directory, "output")
            if bundle:
                artifacts["covariates"] = _covariates(Path(bundle["by_kind"]["bqsr-covariate-table"][0]["path"]), directory)
        else:
            if stage == "genotype":
                _require(input_path is not None, "Final VCF audit requires its branch's gVCF")
                artifacts["input"] = _vcf_audit(input_path, directory, "input", reference, reference_entries, scope_bed, run, gvcf=True)
            artifacts["output"] = _vcf_audit(output_path, directory, "output", reference, reference_entries, scope_bed, run, gvcf=stage == "hc")
            if stage == "genotype":
                result["empty_callset_interpretation"] = "A structurally valid empty final VCF is allowed only with complete input gVCF coverage. Biological sensitivity remains unmeasured."
        result["contract_status"] = "checked_match"
        result["meaning"] = "Stage-specific structural/conservation contract checked; not equality to GATK and not biological accuracy."
    except Exception as exc:
        result.update(contract_status="failed", error_type=type(exc).__name__, reason=str(exc),
                      evidence=getattr(exc, "evidence", {}), traceback=traceback.format_exc())
    _json(directory / "audit.json", result)
    return result


def _table_rows(db):
    for table, key, values, line in db.execute("SELECT name,key,values_json,line FROM cells ORDER BY name,key"):
        yield (table, key), (json.loads(values), line)


def _report_difference(left, right, path):
    changed_rows = changed_cells = 0
    numeric_summary = {}
    with closing(sqlite3.connect(left["sqlite"])) as a, closing(sqlite3.connect(right["sqlite"])) as b, path.open("w") as stream:
        def schema(db):
            return {name: {"columns": json.loads(columns), "numeric": json.loads(numeric)}
                    for name, columns, numeric in db.execute("SELECT name,columns_json,numeric_json FROM tables ORDER BY name")}
        sa, sb = schema(a), schema(b)
        schema_differences = {key: {"left": sa.get(key), "right": sb.get(key)}
                              for key in sorted(sa.keys() | sb.keys()) if sa.get(key) != sb.get(key)}
        for (table, key), x, y in _merge(_table_rows(a), _table_rows(b)):
            differences = {}
            if x is not None and y is not None:
                for column in sorted(x[0].keys() | y[0].keys()):
                    xv, yv = x[0].get(column), y[0].get(column)
                    numeric = column in sa[table]["numeric"] and column in sb[table]["numeric"]
                    if table == "Arguments" and column == "Value" and xv is not None and yv is not None:
                        numeric = all(re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][+-]?[0-9]+)?", v) for v in (xv, yv))
                    if numeric and xv not in (None, "") and yv not in (None, ""):
                        xn, yn = _number(xv), _number(yv)
                        if Decimal(xn) == Decimal(yn):
                            continue
                        delta = _difference(xn, yn)
                        differences[column] = {"left": xv, "right": yv, "delta_left_minus_right": delta}
                        name = table + "." + column
                        summary = numeric_summary.setdefault(name, {"different_cells": 0, "max_absolute_delta": "0"})
                        summary["different_cells"] += 1
                        magnitude = Decimal(delta).copy_abs()
                        if magnitude > Decimal(summary["max_absolute_delta"]):
                            summary["max_absolute_delta"] = str(magnitude)
                    elif xv != yv:
                        differences[column] = {"left": xv, "right": yv}
                if not differences:
                    continue
            changed_rows += 1
            changed_cells += len(differences)
            _line(stream, {"table": table, "semantic_key": json.loads(key), "left": x, "right": y, "cell_differences": differences})
    return {"different_rows": changed_rows, "different_cells": changed_cells, "schema_differences": schema_differences,
            "numeric_summary": numeric_summary, "differences": str(path),
            "numeric_rule": "Exact decimal comparison at published precision; no numerical tolerance and no rounded-away differences."}


def _vcf_difference(left, right, path):
    changed = left_only = right_only = 0
    def rows(db):
        for key, n, first in db.execute("SELECT key,n,first_record FROM records ORDER BY key"):
            yield key, (n, first)
    with closing(sqlite3.connect(left["sqlite"])) as a, closing(sqlite3.connect(right["sqlite"])) as b, path.open("w") as stream:
        for key, x, y in _merge(rows(a), rows(b)):
            xn, yn = x[0] if x else 0, y[0] if y else 0
            if xn != yn:
                changed += 1
                left_only += max(0, xn - yn)
                right_only += max(0, yn - xn)
                _line(stream, {"record": json.loads(key), "native_count": xn, "gatk_count": yn,
                               "native_vcf_record": x[1] if x else None, "gatk_vcf_record": y[1] if y else None})
    return {"different_keys": changed, "native_only_records": left_only, "gatk_only_records": right_only,
            "differences": str(path), "native_vcf": left["vcf"], "gatk_vcf": right["vcf"],
            "meaning": "Representation-preserving record multiset comparison, not haplotype equivalence or GIAB scoring. INFO/FORMAT order and unphased GT allele order are immaterial; allele representation and annotations are retained."}


def _stage_artifacts(branch, stage, branch_name, reference, scope_bed):
    entry = branch.get("stages", {}).get(stage, {})
    contract = entry.get("contract", {})
    _require(contract.get("contract_status") in ("checked_match", "checked_differences"),
             "Comparison requires that branch's successful stage contract", stage=stage,
             execution_status=entry.get("execution_status"), contract_status=contract.get("contract_status"))
    _require(contract.get("branch") == branch_name and contract.get("stage") == stage,
             "Borrowed stage contract", stage=stage, expected_branch=branch_name)
    for key, path in (("reference", reference), ("scope_bed", scope_bed)):
        previous = contract["consumed"][key]
        current = _metadata(Path(path).resolve())
        _require(all(previous.get(field) == value for field, value in current.items()),
                 "Comparison reference/scope differs from audited input", stage=stage, input_kind=key)
    input_key = {"bqsr": "markdup_bam", "apply": "markdup_bam", "hc": "recal_bam", "genotype": "gvcf"}.get(stage)
    if input_key is not None:
        _require(Path(contract["consumed"]["input"]["path"]).resolve() == Path(branch["outputs"][input_key]).resolve(),
                 "Stage consumed another branch's intermediate", stage=stage)
    artifacts = contract["artifacts"]
    output = Path(branch["outputs"][OUTPUTS[stage]]).resolve()
    _require(Path(artifacts["output"]["metadata"]["path"]).resolve() == output,
             "Stage audit describes a different branch output", stage=stage)
    _require(_metadata(output) == {k: v for k, v in artifacts["output"]["metadata"].items()
                                  if k in ("path", "bytes", "mtime_ns", "device", "inode")},
             "Audited stage output changed before comparison", stage=stage)
    return artifacts


def compare_branches(native: dict, gatk: dict, *, reference: Path, scope_bed: Path, work: Path, tools) -> dict:
    """Compare independent, already-audited branches without substituting inputs."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 1, "contract_status": "checked_match", "accuracy_status": "not_evaluated", "stages": {},
              "meaning": "Native/GATK contracts are independent of GIAB biological truth; observed differences are not accuracy failures."}
    for stage in STAGES:
        directory = work / stage
        directory.mkdir(parents=True, exist_ok=True)
        item = {"contract_status": "failed", "metadata_path": str(directory / "comparison.json")}
        result["stages"][stage] = item
        try:
            a = _stage_artifacts(native, stage, "native", reference, scope_bed)
            b = _stage_artifacts(gatk, stage, "gatk", reference, scope_bed)
            ap, bp = Path(a["output"]["metadata"]["path"]), Path(b["output"]["metadata"]["path"])
            _require(not ap.samefile(bp), "Branches borrowed the same output inode", native=str(ap), gatk=str(bp))
            item["native_artifacts"], item["gatk_artifacts"] = a, b
            changed = False
            if stage in ("markdup", "apply"):
                item["records"] = _bam_difference(a["output"], b["output"], directory / "record-differences.jsonl", "exact")
                changed = bool(item["records"]["different_keys"])
                item["duplicate_record_sets"] = _bam_difference(a["output"], b["output"], directory / "duplicate-set-differences.jsonl", "duplicates")
                item["headers"] = {key: {"native": a["output"]["header"][key], "gatk": b["output"]["header"][key]}
                                   for key in ("SQ", "RG") if a["output"]["header"][key] != b["output"]["header"][key]}
                changed |= bool(item["headers"])
                if stage == "markdup":
                    item["metrics"] = _report_difference(a["metrics"], b["metrics"], directory / "metric-differences.jsonl")
                    changed |= bool(item["metrics"]["different_rows"] or item["metrics"]["schema_differences"])
                else:
                    item["qualities"] = _quality_difference(a["output"], b["output"], directory, keep_duplicate=False)
                    item["qualities"]["boundary"] = "Cross-branch quality pairing excludes upstream duplicate bit only; exact record/duplicate-set comparison above still reports every flag difference."
            elif stage == "bqsr":
                item["reports"] = _report_difference(a["output"], b["output"], directory / "report-differences.jsonl")
                changed = bool(item["reports"]["different_rows"] or item["reports"]["schema_differences"])
            else:
                item["records"] = _vcf_difference(a["output"], b["output"], directory / "variant-differences.jsonl")
                changed = bool(item["records"]["different_keys"])
                item["header_evidence"] = {"native": a["output"]["header"], "gatk": b["output"]["header"]}
                item["header_definitions"] = {key: {"native": a["output"]["definitions"].get(key), "gatk": b["output"]["definitions"].get(key)}
                                              for key in ("INFO", "FORMAT", "FILTER", "ALT")
                                              if a["output"]["definitions"].get(key) != b["output"]["definitions"].get(key)}
                changed |= bool(item["header_definitions"])
            item["contract_status"] = "checked_differences" if changed else "checked_match"
        except Exception as exc:
            item.update(contract_status="failed", reason=str(exc), error_type=type(exc).__name__,
                        evidence=getattr(exc, "evidence", {}), traceback=traceback.format_exc())
        _json(directory / "comparison.json", item)
    statuses = {item["contract_status"] for item in result["stages"].values()}
    result["contract_status"] = "failed" if "failed" in statuses else ("checked_differences" if "checked_differences" in statuses else "checked_match")
    return result
