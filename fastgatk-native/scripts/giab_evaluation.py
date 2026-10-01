#!/usr/bin/env python3
"""Independent, fail-closed hap.py/vcfeval scoring; never a biological pass flag."""
from __future__ import annotations

import csv
from contextlib import closing
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import re
import shlex
import sqlite3

from giab_assets import _parse_md5, _stratification_checksums, file_digest, load_lock
from giab_contracts import _reference
from giab_environment import PinnedTools, save_json
from giab_preflight import SCOPES, _require, _vcf_header, _vcf_rows

CHECKED = {"checked_match", "checked_differences"}
COUNTS = {"truth_total": "TRUTH.TOTAL", "truth_tp": "TRUTH.TP", "truth_fn": "TRUTH.FN",
          "query_total": "QUERY.TOTAL", "query_tp": "QUERY.TP", "query_fp": "QUERY.FP",
          "query_unk": "QUERY.UNK", "fp_gt": "FP.gt", "fp_al": "FP.al"}
ERROR_COLUMNS = ("side", "error_class", "chrom", "pos", "ref", "alt", "gt", "decision",
                 "decision_subtype", "variant_type", "filter", "regions", "superlocus", "record_identity")
POLICY = {"engine": "vcfeval", "gender": "male", "threads": 16,
          "primary_filter": "PASS", "retained_filters": ["ALL", "PASS"],
          "genotype_aware": True, "truth_prefiltered_by_confident_bed": False,
          "precision": "QUERY.TP / (QUERY.TP + QUERY.FP)",
          "recall": "TRUTH.TP / (TRUTH.TP + TRUTH.FN)",
          "f1": "2*precision*recall/(precision+recall); null for undefined or zero denominator",
          "unknown_query": "Outside confident BED is UNK, not FP",
          "missing_class": "not_evaluated, never fabricated zero counts",
          "strata": "Overlapping layers cannot be summed",
          "error_identity": "Exact annotated-record [chrom,pos,ref,alt,gt]; FP overlap is diagnostic, not haplotype-aware shared biological errors",
          "accuracy_threshold": None, "timeout_seconds": 12 * 3600}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load(path: Path) -> dict:
    with path.open() as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), "Expected JSON object", path=str(path))
    return value


def _path(value, label: str) -> Path:
    _require(isinstance(value, (str, Path)) and bool(str(value)), "Missing " + label)
    path = Path(value).resolve()
    _require(path.is_file(), "Missing " + label, path=str(path))
    return path


def _error(error: Exception) -> dict:
    return {"error": str(error), "error_type": type(error).__name__,
            "evidence": getattr(error, "evidence", {})}


def _count(value) -> int | None:
    if value is None or str(value).strip().lower() in ("", ".", "na", "nan", "null"):
        return None
    try:
        number = Decimal(str(value))
    except InvalidOperation as error:
        raise ValueError(f"Invalid comparator count: {value!r}") from error
    _require(number.is_finite() and number >= 0 and number == number.to_integral_value(),
             "Invalid nonnegative integer comparator count", value=value)
    return int(number)


def _metric(row: dict, confident_bases: int, *, missing=False) -> dict:
    result = {"variant_type": row["Type"], "filter": row["Filter"],
              "subset": row.get("Subset", "*"), "subtype": row.get("Subtype", "*"),
              "genotype": row.get("Genotype", "*"),
              **{key: _count(row.get(column)) for key, column in COUNTS.items()}}
    result["confident_bases"] = (confident_bases if result["subset"] == "*" else
                                  _count(row.get("Subset.IS_CONF.Size")))
    p_den = None if result["query_tp"] is None or result["query_fp"] is None else result["query_tp"] + result["query_fp"]
    r_den = None if result["truth_tp"] is None or result["truth_fn"] is None else result["truth_tp"] + result["truth_fn"]
    precision = result["query_tp"] / p_den if p_den else None
    recall = result["truth_tp"] / r_den if r_den else None
    result.update(precision=precision, recall=recall,
                  f1=2 * precision * recall / (precision + recall)
                  if precision is not None and recall is not None and precision + recall else None)
    result["undefined_metrics"] = [key for key in ("precision", "recall", "f1") if result[key] is None]
    result["missing_counts"] = [key for key in COUNTS if result[key] is None]
    result["metric_status"] = "not_evaluated" if missing or result["missing_counts"] else "measured"
    result["raw_comparator"] = row
    if missing:
        result["reason"] = "Comparator omitted this class/filter; absence is not evidence of zero counts"
    return result


def parse_happy(prefix: Path, confident_bases: int) -> list[dict]:
    """Parse unthresholded extended rows, preserving unequal truth/query TPs.

    Every raw extended row remains in the original CSV. Structured rows retain
    all SNP/INDEL subsets/subtypes/genotypes at QQ='*'. Missing global ALL/PASS
    rows are explicit null/not_evaluated entries, not inferred empty classes.
    """
    _require(type(confident_bases) is int and confident_bases >= 0, "Invalid confident base count")
    prefix = Path(prefix)
    summary = {}
    with Path(str(prefix) + ".summary.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        _require(reader.fieldnames and {"Type", "Filter", "TRUTH.TP", "QUERY.FP"} <= set(reader.fieldnames),
                 "Malformed hap.py summary CSV")
        for row in reader:
            _require(None not in row, "Malformed summary CSV row")
            key = (row["Type"], row["Filter"])
            _require(key not in summary, "Duplicate summary row", key=key)
            summary[key] = row
    metrics, seen = [], set()
    with Path(str(prefix) + ".extended.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"Type", "Filter", "Subset", "Subtype", "Genotype", "QQ", "QUERY.TP", *COUNTS.values()}
        _require(reader.fieldnames and required <= set(reader.fieldnames), "Malformed hap.py extended CSV")
        for row in reader:
            _require(None not in row, "Malformed extended CSV row")
            if row["QQ"] != "*" or row["Type"] not in ("SNP", "INDEL"):
                continue
            key = tuple(row[k] for k in ("Type", "Filter", "Subset", "Subtype", "Genotype"))
            _require(key not in seen, "Duplicate unthresholded comparator metric", key=key)
            seen.add(key)
            metric = _metric(row, confident_bases)
            if key[2:] == ("*", "*", "*"):
                raw_summary = summary.get(key[:2])
                if raw_summary is None:
                    metric.update(metric_status="not_evaluated", reason="Global extended row lacks corresponding summary row")
                else:
                    metric["raw_summary"] = raw_summary
                    for column in COUNTS.values():
                        if column in raw_summary:
                            _require(_count(raw_summary[column]) == _count(row.get(column)),
                                     "Summary/extended count disagreement", key=key, column=column)
            metrics.append(metric)
    for variant_type in ("SNP", "INDEL"):
        for filtering in ("ALL", "PASS"):
            if (variant_type, filtering, "*", "*", "*") not in seen:
                metric = _metric({"Type": variant_type, "Filter": filtering}, confident_bases, missing=True)
                if (variant_type, filtering) in summary:
                    metric["raw_summary"] = summary[variant_type, filtering]
                    metric["reason"] = "Summary row lacks extended QUERY.TP evidence; no truth-side TP substitution"
                metrics.append(metric)
    return metrics


def _error_records(vcf: Path, output: Path) -> dict:
    counts = {"FN": 0, "FP": 0, "UNK": 0}
    with gzip.open(vcf, "rt") as stream, output.open("x", newline="") as target:
        _, contigs, columns = _vcf_header(stream)
        _require(columns[8:] == ["FORMAT", "TRUTH", "QUERY"], "Unexpected annotated comparator samples", columns=columns)
        writer = csv.DictWriter(target, fieldnames=ERROR_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for fields, _ in _vcf_rows(stream, contigs):
            _require(len(fields) == 11, "Malformed annotated comparator record")
            keys = fields[8].split(":")
            _require({"GT", "BD", "BK", "BVT"} <= set(keys), "Comparator lacks decision/type FORMAT fields")
            info = dict(item.split("=", 1) if "=" in item else (item, "") for item in fields[7].split(";"))
            for side, sample in (("truth", fields[9]), ("query", fields[10])):
                values = dict(zip(keys, sample.split(":")))
                decision = values.get("BD", ".")
                error_class = ("FN" if side == "truth" and decision == "FN" else
                               "FP" if side == "query" and decision == "FP" else
                               "UNK" if side == "query" and decision in ("N", "UNK") else None)
                if error_class is None:
                    continue
                gt = values.get("GT", ".")
                identity = [fields[0], int(fields[1]), fields[3], fields[4], gt]
                writer.writerow({"side": side, "error_class": error_class, "chrom": fields[0], "pos": fields[1],
                                 "ref": fields[3], "alt": fields[4], "gt": gt, "decision": decision,
                                 "decision_subtype": values.get("BK", "."), "variant_type": values.get("BVT", "."),
                                 "filter": values.get("FT", fields[6]), "regions": info.get("Regions", "."),
                                 "superlocus": info.get("BS", "."),
                                 "record_identity": json.dumps(identity, separators=(",", ":"))})
                counts[error_class] += 1
    return {"counts": counts, "identity_semantics": POLICY["error_identity"],
            "count_semantics": "Annotated records, not hap.py allele counts; UNK retained separately"}


def _digest(path: Path, expected: dict, label: str) -> dict:
    actual = file_digest(path)
    _require(expected.get("sha256") == actual["sha256"], "Changed or unfrozen " + label,
             path=str(path), expected=expected.get("sha256"), observed=actual["sha256"])
    return actual


def _tool_eligibility(environment: dict) -> dict:
    lock = load_lock()
    cpus = environment.get("cpus", [])
    _require(environment.get("status") == "eligible" and environment.get("resource_status") == "enforced_within_limit",
             "Evaluator environment has not demonstrated resource enforcement")
    _require(len(cpus) == 16 and len(set(cpus)) == 16 and all(type(c) is int and c >= 0 for c in cpus)
             and environment.get("memory_bytes") == 64 * 1024 ** 3, "Evaluator requires 16 CPUs and 64 GiB")
    for tool in ("samtools", "bcftools", "happy"):
        observed = environment.get("images", {}).get(tool, {})
        requested = lock["images"][tool]
        _require(observed.get("requested") == requested and observed.get("id") and
                 any(d.split("@")[-1] == requested.split("@")[-1] for d in observed.get("repo_digests", [])),
                 "Unattested pinned evaluator image", tool=tool)
    versions = environment.get("versions", {})
    for tool in ("samtools", "bcftools"):
        _require(versions.get(tool, {}).get("output", "").startswith(tool + " 1.21\n"), "Missing pinned tool version", tool=tool)
    proof = versions.get("happy", {}).get("source_attestation", {})
    _require(proof.get("version") == "v0.3.12" and proof.get("image") == lock["images"]["happy"]
             and proof.get("source_sha256") == "1ba8d0245b8c7d6a68ca5b78ab8d0d87104e9a3b36c86169319c5add2395b6a5"
             and proof.get("files") and all(proof.get(k) == [] for k in ("missing", "different", "additional")),
             "Missing exact official hap.py v0.3.12 source attestation")
    _require(versions.get("rtg", {}).get("output", "").startswith("Product: RTG Tools 3.10.1\n"),
             "Missing pinned actual RTG/vcfeval version")
    return {"images": environment["images"], "versions": versions}


def _bed_rows(path: Path, reference: dict):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split()
            _require(len(fields) >= 3, "Malformed BED", path=str(path), line=number)
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            _require(chrom in reference and 0 <= start < end <= reference[chrom][0],
                     "BED interval outside reference", path=str(path), line=number, interval=fields[:3])
            yield chrom, start, end


def _beds(scope_bed: Path, evaluation_bed: Path, reference: dict, work: Path) -> dict:
    scope = {}
    for chrom, start, end in _bed_rows(scope_bed, reference):
        _require(chrom not in scope, "Scope must declare each whole chromosome once")
        _require(start == 0 and end == reference[chrom][0], "Scope cannot shrink the declared chromosomes")
        scope[chrom] = end
    _require(scope, "Empty scope BED")
    database = work / "evaluation-bed.sqlite"
    _require(not database.exists(), "Refusing reused BED audit")
    bases = records = 0
    with closing(sqlite3.connect(database)) as db:
        db.execute("PRAGMA temp_store=FILE")
        db.execute("CREATE TABLE intervals (chrom TEXT, start INTEGER, end INTEGER, PRIMARY KEY(chrom,start,end)) WITHOUT ROWID")
        for chrom, start, end in _bed_rows(evaluation_bed, reference):
            _require(chrom in scope, "Evaluation BED leaks outside calling scope", chrom=chrom)
            db.execute("INSERT OR IGNORE INTO intervals VALUES (?,?,?)", (chrom, start, end))
            records += 1
        db.commit()
        current, begin, end = None, 0, 0
        for chrom, start, stop in db.execute("SELECT chrom,start,end FROM intervals ORDER BY chrom,start,end"):
            if chrom != current or start > end:
                bases += end - begin
                current, begin, end = chrom, start, stop
            else:
                end = max(end, stop)
        bases += end - begin
    _require(bases > 0, "Empty confident evaluation domain")
    return {"confident_bases": bases, "bed_records": records, "chromosomes": list(scope),
            "scope_digest": file_digest(scope_bed), "evaluation_digest": file_digest(evaluation_bed),
            "union_evidence": str(database)}


def _stratifications(tsv: Path | None, reference: dict, *, data_root: Path, synthetic: bool) -> dict:
    if tsv is None:
        _require(synthetic, "Production requires all official stratifications")
        return {"status": "not_requested", "kind": "synthetic", "layers": [], "tsv": None}
    tsv = _path(tsv, "stratification TSV")
    receipt = None
    if not synthetic:
        lock = load_lock()
        directory = data_root / "prepared" / "stratifications"
        _require(tsv == (directory / lock["stratification"]["entrypoint"]).resolve(), "Unpinned official stratification entrypoint")
        receipt_path = data_root / "prepared" / "stratifications-manifest.json"
        receipt = _load(receipt_path)
        _require(receipt.get("status") == "verified", "Stratification acquisition is not verified")
        manifests = {}
        for name in ("stratification_manifest", "stratification_member_manifest"):
            source = data_root / "raw" / lock["assets"][name]["filename"]
            _digest(source, {"sha256": receipt.get("source_sha256", {}).get(name)}, name)
            manifests[name] = _parse_md5(source)
        checksums = _stratification_checksums(lock, manifests)
        _digest(tsv, receipt.get("members", {}).get(lock["stratification"]["entrypoint"], {}), "official TSV")
        _require(file_digest(tsv)["md5"] == lock["stratification"]["entrypoint_md5"], "Official TSV MD5 mismatch")
    layers, seen = [], set()
    with tsv.open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            _require(len(fields) == 2 and fields[0] and fields[0] not in seen, "Invalid stratification TSV", line=number)
            label, relative = fields
            seen.add(label)
            bed = (tsv.parent / relative).resolve()
            _require(bed.is_relative_to(tsv.parent) and bed.is_file(), "Stratification path escapes tree or is missing", path=str(bed))
            digest = file_digest(bed)
            if receipt is not None:
                member = bed.relative_to(directory).as_posix()
                expected = receipt.get("members", {}).get(member, {})
                _require(digest["sha256"] == expected.get("sha256") and digest["md5"] == expected.get("md5")
                         and digest["md5"] == checksums.get(member),
                         "Official BED checksum mismatch", member=member)
            records = sum(1 for _ in _bed_rows(bed, reference))
            other_sample = any(sample != "HG002" for sample in re.findall(r"HG00[1-7]", label + " " + relative))
            xy = bool(re.search(r"(?:^|[/_])(?:XY|chrX|chrY)(?:[/_.]|$)", relative))
            layers.append({"label": label, "path": str(bed), "digest": digest, "records": records,
                           "reference_bounds": "checked", "other_sample_specific": other_sample,
                           "outside_autosomal_biology": xy,
                           "interpretation": "synthetic fixture only" if synthetic else
                           "not HG002/autosomal biology" if other_sample or xy else "overlapping official HG002-applicable layer"})
    _require(layers, "Empty stratification TSV")
    if receipt is not None:
        _require([(x["label"], str(Path(x["path"]).relative_to(directory))) for x in layers] ==
                 [(x["label"], x["member"]) for x in receipt.get("layers", [])], "Official layer set changed")
    return {"status": "eligible", "kind": "synthetic" if synthetic else "GIAB_v3.6",
            "tsv": str(tsv), "layers": layers, "overlap_policy": POLICY["strata"],
            "manifest": str(receipt_path) if receipt is not None else None,
            "source_sha256": receipt.get("source_sha256") if receipt is not None else None}


def _tools(data_root: Path, output: Path, environment: dict, mounts: list[Path]):
    runner = PinnedTools(data_root, load_lock()["images"], environment["cpus"], environment["memory_bytes"], mounts=mounts)
    runner.logs = output / "tool-logs"
    runner.scratch = output / "scratch" / "validation"
    runner.scratch.mkdir(parents=True)
    counter = 0
    def run(tool, args, name, *, stdout=None):
        nonlocal counter
        counter += 1
        return runner(tool, args, f"{counter:04d}-{name}", stdout=stdout)
    return run


def _check_vcf(path: Path, reference: dict, fasta: Path, scope_bed: Path, chromosomes: list[str], work: Path, run, name: str) -> dict:
    _require(path.name.endswith(".vcf.gz"), "Comparator requires bgzip VCF", path=str(path))
    indices = [p for p in (Path(str(path) + ".tbi"), Path(str(path) + ".csi")) if p.is_file() and p.stat().st_size]
    _require(indices, "Missing VCF index", path=str(path))
    run("bcftools", ["norm", "-f", str(fasta), "-c", "e", "-Ou", "-o", "/dev/null", str(path)], name + "-ref-check")
    rendered = []
    for indexed in (False, True):
        target = work / (name + (".indexed" if indexed else ".sequential") + ".vcf.gz")
        args = ["view", "--no-version"]
        if indexed:
            args += ["-R", str(scope_bed), "--regions-overlap", "1"]
        run("bcftools", [*args, "-Oz", "-o", str(target), str(path)], name + ("-index-query" if indexed else "-decode"))
        digest = hashlib.sha256()
        count = filtered = 0
        with gzip.open(target, "rt", encoding="ascii") as stream:
            _, contigs, columns = _vcf_header(stream)
            _require(columns[8:] == ["FORMAT", "HG002"], "VCF must have exactly sample HG002", path=str(path), columns=columns)
            _require(all(contigs.get(c) == reference[c][0] for c in chromosomes), "VCF scope dictionary differs from reference", path=str(path))
            for fields, line in _vcf_rows(stream, contigs):
                _require(len(fields) == 10 and fields[0] in chromosomes, "VCF records escape scope or sample schema", line=line[:1000])
                alt = fields[4].split(",")
                _require(alt and "." not in alt and "<NON_REF>" not in alt and "<*>" not in alt,
                         "Final comparator input contains reference-confidence blocks", line=line[:1000])
                fmt, values = fields[8].split(":"), fields[9].split(":")
                _require(len(set(fmt)) == len(fmt) and "GT" in fmt and fmt.index("GT") < len(values), "Missing/duplicate GT FORMAT")
                gt = values[fmt.index("GT")]
                _require(re.fullmatch(r"(?:\.|[0-9]+)(?:[/|](?:\.|[0-9]+))?", gt), "Malformed GT", gt=gt)
                alleles = re.split(r"[/|]", gt)
                _require((gt == "." or len(alleles) == 2) and all(a == "." or int(a) <= len(alt) for a in alleles), "Invalid diploid GT", gt=gt)
                _require(not all(a == "0" for a in alleles), "Final input contains hom-reference record")
                digest.update(line.encode("ascii"))
                count += 1
                filtered += fields[6] not in ("PASS", ".")
        rendered.append({"path": str(target), "records": count, "filtered_records": filtered, "records_sha256": digest.hexdigest()})
    _require(all(rendered[0][k] == rendered[1][k] for k in ("records", "records_sha256")), "VCF index loses or changes full-scope records", path=str(path))
    return {"status": "eligible", "path": str(path), "digest": file_digest(path), "indices": {str(p): file_digest(p) for p in indices},
            "full_scope_ref_checked": True, "confident_bed_prefilter": False, "decodes": rendered,
            "records": rendered[0]["records"], "filtered_records": rendered[0]["filtered_records"]}


def _score(*, truth: Path, query: Path, reference: Path, scope_bed: Path, evaluation_bed: Path,
           data_root: Path, output: Path, environment: dict, experiment_kind: str,
           common: dict, mounts: list[Path]) -> dict:
    from giab_execution import execute_stage

    result = {"accuracy_status": "not_evaluated", "execution": {"execution_status": "not_run", "resource_status": "not_measured"},
              "metrics": [], "artifacts": {}, "errors": [], "experiment_kind": experiment_kind}
    prefix = output / "happy"
    result["artifacts"] = {key: str(Path(str(prefix) + suffix)) for key, suffix in
                           (("summary_csv", ".summary.csv"), ("extended_csv", ".extended.csv"),
                            ("metrics_json_gz", ".metrics.json.gz"), ("evaluation_vcf", ".vcf.gz"),
                            ("evaluation_vcf_index", ".vcf.gz.tbi"), ("runinfo_json", ".runinfo.json"))}
    result["artifacts"].update(error_records_tsv=str(output / "error-records.tsv"),
                               score_manifest=str(output / "score.json"), tool_logs=str(output / "tool-logs"),
                               execution_manifest=str(output / "execution" / "stage.json"))
    save_json(output / "score.json", result)
    try:
        tools = _tools(data_root, output, environment, mounts)
        work = output / "validation"
        work.mkdir()
        result["query_validation"] = _check_vcf(query, common["reference"], reference, scope_bed,
                                                 common["domain"]["chromosomes"], work, tools, "query")
        save_json(output / "score.json", result)
        scratch = output / "scratch" / "happy"
        scratch.mkdir()
        command = ["/opt/hap.py/bin/hap.py", str(truth), str(query), "-r", str(reference),
                   "-f", str(evaluation_bed), "-o", str(prefix), "--engine", "vcfeval", "--gender", "male",
                   "--threads", "16", "--scratch-prefix", str(scratch), "--no-roc", "--no-fixchr",
                   "--no-adjust-conf-regions", "--write-vcf", "--write-counts"]
        if common["stratification"].get("tsv"):
            command += ["--stratification", common["stratification"]["tsv"]]
        result["reproduction"] = {"argv": command, "command": shlex.join(command), "image": load_lock()["images"]["happy"],
                                  "note": "Use the saved execution Docker envelope to preserve CPU/memory/swap limits"}
        result["execution"] = execute_stage(command, data_root=data_root, stage_dir=output / "execution",
                                             scratch=scratch, environment=environment, timeout_seconds=12 * 3600,
                                             name="happy-" + output.name, mounts=mounts, image=load_lock()["images"]["happy"],
                                             extra_env={"RTG_MEM": "40g", "RTG_JAVA_OPTS": "-XX:ActiveProcessorCount=16"})
        _require(result["execution"].get("execution_status") == "completed", "Comparator workload did not complete")
        _require(result["execution"].get("resource_status") in ("enforced_within_limit", "limit_hit"),
                 "Comparator resource enforcement was not demonstrated")
        for key in ("summary_csv", "extended_csv", "metrics_json_gz", "evaluation_vcf", "evaluation_vcf_index", "runinfo_json"):
            _path(result["artifacts"][key], key)
        with gzip.open(result["artifacts"]["metrics_json_gz"], "rt") as stream:
            raw_metrics = json.load(stream)
        _require(isinstance(raw_metrics, dict) and "metrics" in raw_metrics, "Malformed comparator metrics JSON")
        _load(Path(result["artifacts"]["runinfo_json"]))
        result["metrics"] = parse_happy(prefix, common["domain"]["confident_bases"])
        result["error_records"] = _error_records(Path(result["artifacts"]["evaluation_vcf"]), Path(result["artifacts"]["error_records_tsv"]))
        measured = [row for row in result["metrics"] if row["metric_status"] == "measured" and
                    row["subset"] == row["subtype"] == row["genotype"] == "*" and row["filter"] in ("PASS", "ALL")]
        result["accuracy_status"] = "measured" if measured else "not_evaluated"
        result["metric_completeness"] = {"missing_global_rows": [dict(variant_type=r["variant_type"], filter=r["filter"])
                                                                 for r in result["metrics"] if r["subset"] == r["subtype"] == r["genotype"] == "*"
                                                                 and r["metric_status"] != "measured"]}
        result["artifact_digests"] = {key: file_digest(Path(result["artifacts"][key])) for key in
                                     ("summary_csv", "extended_csv", "metrics_json_gz", "evaluation_vcf", "evaluation_vcf_index", "error_records_tsv", "runinfo_json")}
    except Exception as error:
        result["accuracy_status"] = "blocked"
        result["errors"].append(_error(error))
        stage_manifest = output / "execution" / "stage.json"
        if result["execution"].get("execution_status") == "not_run" and stage_manifest.is_file():
            try:
                result["execution"] = _load(stage_manifest)
            except Exception as evidence_error:
                result["errors"].append(_error(evidence_error))
    result["ended_at"] = _now()
    save_json(output / "score.json", result)
    return result


def _common(*, truth: Path, reference: Path, scope_bed: Path, evaluation_bed: Path, data_root: Path,
            output: Path, environment: dict, experiment_kind: str, stratification_tsv: Path | None,
            mounts: list[Path]) -> dict:
    result = {"tool_eligibility": _tool_eligibility(environment), "reference": _reference(reference)}
    result["domain"] = _beds(scope_bed, evaluation_bed, result["reference"], output)
    result["stratification"] = _stratifications(stratification_tsv, result["reference"], data_root=data_root,
                                               synthetic=experiment_kind == "synthetic_smoke")
    tools = _tools(data_root, output, environment, mounts)
    result["truth_validation"] = _check_vcf(truth, result["reference"], reference, scope_bed,
                                            result["domain"]["chromosomes"], output, tools, "truth")
    save_json(output / "eligibility.json", result)
    return result


def score_callset(*, truth: Path, query: Path, reference: Path, scope_bed: Path, evaluation_bed: Path,
                  data_root: Path, output: Path, environment: dict, experiment_kind: str,
                  stratification_tsv: Path | None = None, mounts: list[Path] | None = None) -> dict:
    """Execute one explicitly synthetic comparator case in a fresh directory.

    This helper cannot confer GIAB eligibility. Production uses evaluate_run,
    which requires a completed independent pipeline and input audit snapshot.
    Truth/query are indexed bgzip VCFs with exactly sample HG002. Whole tiny
    reference chromosomes belong in scope_bed; confident BED may be narrower.
    """
    output, data_root = Path(output).resolve(), Path(data_root).resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {"accuracy_status": "blocked", "execution": {"execution_status": "not_run", "resource_status": "not_measured"},
              "metrics": [], "artifacts": {}, "errors": [], "experiment_kind": experiment_kind}
    try:
        _require(experiment_kind == "synthetic_smoke", "Direct score_callset requires explicit synthetic_smoke; production must use evaluate_run")
        truth, query, reference, scope_bed, evaluation_bed = (_path(p, label) for p, label in
                ((truth, "truth"), (query, "query"), (reference, "reference"), (scope_bed, "scope BED"), (evaluation_bed, "evaluation BED")))
        mounts = list(dict.fromkeys([output, *(Path(p).resolve() for p in mounts or []),
                                     *(p.parent for p in (truth, query, reference, scope_bed, evaluation_bed)),
                                     *([Path(stratification_tsv).resolve().parent] if stratification_tsv else [])]))
        work = output / "eligibility"
        work.mkdir()
        common = _common(truth=truth, reference=reference, scope_bed=scope_bed, evaluation_bed=evaluation_bed,
                         data_root=data_root, output=work, environment=environment, experiment_kind=experiment_kind,
                         stratification_tsv=stratification_tsv, mounts=mounts)
        result = _score(truth=truth, query=query, reference=reference, scope_bed=scope_bed, evaluation_bed=evaluation_bed,
                        data_root=data_root, output=output, environment=environment, experiment_kind=experiment_kind,
                        common=common, mounts=mounts)
        result.update(eligibility=str(work / "eligibility.json"), stratification=common["stratification"],
                      interpretation="Synthetic comparator semantics only; not HG002/GIAB accuracy evidence")
    except Exception as error:
        result["errors"].append(_error(error))
    save_json(output / "score.json", result)
    return result


def _production_gate(manifest: dict, scope: str, paths: dict) -> None:
    info = manifest.get("scopes", {}).get(scope, {})
    _require(manifest.get("sample") == "HG002" and manifest.get("preflight_kind") == "complete_inputs"
             and manifest.get("experiment_kind", "giab") == "giab", "Production requires complete real HG002 input preflight")
    _require(info.get("input_status") == "eligible" and info.get("reference_eligibility") == "eligible"
             and info.get("truth", {}).get("status") == "eligible", "Requested scope is not truth/reference eligible", scope=scope)
    gates = manifest.get("gates", {})
    required = ["reference", "dbsnp", "mills", "benchmark_bed", "truth_full_decode", "truth_" + scope,
                "source_bam", "bam_reference", "provenance", "common_bam", "coverage", "aggregate_capacity"]
    if scope == "chr20":
        required.append("chr20_bam")
    _require(all(gates.get(name, {}).get("status") == "eligible" for name in required),
             "Incomplete requested-scope biological/provenance audit", failed=[n for n in required if gates.get(n, {}).get("status") != "eligible"])
    _require(gates["provenance"].get("quality_eligibility") == "eligible", "Original QUAL provenance unknown or blocked")
    _digest(paths["reference"], gates["reference"], "reference")
    for key in ("fai", "dictionary"):
        expected = gates["reference"].get(key, {})
        actual_path = Path(str(paths["reference"]) + ".fai") if key == "fai" else paths["reference"].with_suffix(".dict")
        _digest(actual_path, expected, key)
    _digest(paths["truth"], info["truth"].get("digest", {}), "scope truth")
    _digest(Path(str(paths["truth"]) + ".tbi"), info["truth"].get("index_digest", {}), "truth index")
    _digest(paths["scope_bed"], info.get("scope_digest", {}), "scope BED")
    _digest(paths["evaluation_bed"], info.get("evaluation_digest", {}), "evaluation BED")
    reference = _reference(paths["reference"])
    _require(list(_bed_rows(paths["scope_bed"], reference)) == [(c, 0, reference[c][0]) for c in SCOPES[scope]],
             "Production scope is not the entire declared chromosomes")


def _candidates(report: dict):
    for repetition in report.get("repetitions", []):
        number = repetition["number"]
        for branch in ("native", "gatk"):
            candidate = repetition.get("branches", {}).get(branch, {})
            yield f"rep-{number:02d}/{branch}", number, branch, report.get("profile", {}), candidate, False
    for auxiliary in report.get("auxiliary", []):
        name = auxiliary.get("name", "")
        _require(re.fullmatch(r"[A-Za-z0-9_-]+", name), "Unsafe or missing auxiliary name")
        yield "auxiliary/" + name, None, auxiliary.get("branch", "native"), auxiliary.get("profile"), auxiliary, True


def _candidate_error(candidate: dict, auxiliary: bool, report: dict) -> str | None:
    if candidate.get("execution_status") != "completed" or not candidate.get("scoring_candidate"):
        return "Independent candidate did not complete; partial or borrowed output cannot be scored"
    if candidate.get("contract_status") not in CHECKED:
        return "Candidate output contracts did not pass"
    required = ("hc", "genotype") if auxiliary else ("markdup", "bqsr", "apply", "hc", "genotype")
    checks = [(name, candidate.get("stages", {}).get(name, {})) for name in required]
    if auxiliary:
        source = next((r for r in report.get("repetitions", []) if r.get("number") == candidate.get("source_repetition")), {})
        upstream = source.get("branches", {}).get("native", {})
        if candidate.get("shared_native_recal_bam") != upstream.get("outputs", {}).get("recal_bam"):
            return "Aggregate does not name the independently audited native recal BAM"
        checks += [("source/" + name, upstream.get("stages", {}).get(name, {})) for name in ("markdup", "bqsr", "apply")]
    failed = [name for name, stage in checks if stage.get("execution_status") != "completed" or stage.get("contract_status") not in CHECKED]
    return "Failed upstream execution/contracts: " + ", ".join(failed) if failed else None


def evaluate_run(run: Path, output: Path | None = None) -> dict:
    """Score every independently successful candidate; preserve prior attempts.

    A terminated run may contain failed siblings. Its ended_at, not overall
    execution success, permits evaluation. No unrelated scope is a veto.
    """
    run = Path(run).resolve()
    run_path = run / "run.json" if run.is_dir() else run
    run_dir = run_path.parent
    output = Path(output).resolve() if output is not None else run_dir / "evaluation"
    with (run_dir / ".run.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        report = _load(run_path)
        _require(report.get("ended_at"), "Cannot evaluate a live or unfinished run")
        output.mkdir(parents=True, exist_ok=False)
        # The snapshot preserves every prior evaluation pointer/status before publication.
        save_json(output / "run-before-evaluation.json", report)
        result = {"schema_version": 1, "run_path": str(run_path), "output": str(output),
                  "experiment_kind": report.get("experiment_kind"), "scope": report.get("scope"),
                  "input_status": report.get("input_status", "pending"), "accuracy_status": "not_evaluated",
                  "errors": [], "scoring_policy": POLICY, "stratification": {"status": "not_evaluated", "layers": []},
                  "evaluations": [], "started_at": _now()}
        candidates = []
        common = None
        try:
            candidates = list(_candidates(report))
            for ident, number, branch, configuration, candidate, auxiliary in candidates:
                reason = _candidate_error(candidate, auxiliary, report)
                result["evaluations"].append({"id": ident, "repetition": number, "branch": branch,
                                              "configuration": configuration, "accuracy_status": "not_evaluated",
                                              "execution": {"execution_status": "not_run", "resource_status": "not_measured"},
                                              "metrics": [], "artifacts": {}, "errors": [],
                                              "candidate_status": {key: candidate.get(key) for key in ("execution_status", "contract_status", "resource_status")},
                                              "not_evaluated_reason": reason})
            manifest = _load(_path(report.get("input_manifest"), "run input manifest snapshot"))
            scope = report.get("scope")
            kind = report.get("experiment_kind")
            _require(scope in SCOPES and kind in ("giab", "synthetic_smoke"), "Unknown scope or experiment kind")
            synthetic = kind == "synthetic_smoke"
            _require(not synthetic or manifest.get("experiment_kind") == "synthetic_smoke", "Synthetic run lacks explicit synthetic input declaration")
            info = manifest.get("scopes", {}).get(scope, {})
            _require(info.get("input_status") == "eligible", "Requested scope input gate is not eligible", scope=scope)
            paths = {"truth": _path(info.get("truth", {}).get("path"), "scope truth"),
                     "reference": _path(manifest.get("reference_path"), "reference"),
                     "scope_bed": _path(info.get("scope_bed"), "scope BED"),
                     "evaluation_bed": _path(info.get("evaluation_bed"), "evaluation BED")}
            data_root = Path(report["data_root"]).resolve()
            if not synthetic:
                _production_gate(manifest, scope, paths)
            environment = _load(_path(report.get("environment_manifest_snapshot", report.get("environment_manifest")), "environment snapshot"))
            stratification = manifest.get("stratification_tsv") if synthetic else data_root / "prepared" / "stratifications" / load_lock()["stratification"]["entrypoint"]
            work = output / "eligibility"
            work.mkdir()
            mounts = list(dict.fromkeys([run_dir, output, *(p.parent for p in paths.values()),
                                         *([Path(stratification).resolve().parent] if stratification else [])]))
            common = _common(**paths, data_root=data_root, output=work, environment=environment,
                             experiment_kind=kind, stratification_tsv=stratification, mounts=mounts)
            result.update(stratification=common["stratification"], eligibility=str(work / "eligibility.json"),
                          confident_bases=common["domain"]["confident_bases"], tool_eligibility=common["tool_eligibility"])
            for entry, (_, _, _, _, candidate, _) in zip(result["evaluations"], candidates):
                if entry["not_evaluated_reason"]:
                    continue
                directory = output / entry["id"]
                directory.mkdir(parents=True)
                try:
                    query = _path(candidate.get("outputs", {}).get("vcf"), "independent query VCF")
                    scored = _score(**paths, query=query, data_root=data_root, output=directory, environment=environment,
                                    experiment_kind=kind, common=common, mounts=list(dict.fromkeys([*mounts, query.parent])))
                    entry.update(scored)
                except Exception as error:
                    entry.update(accuracy_status="blocked", errors=[_error(error)])
                save_json(output / "evaluation.json", result)
        except Exception as error:
            result["errors"].append(_error(error))
            for entry in result["evaluations"]:
                if not entry.get("not_evaluated_reason"):
                    entry.update(accuracy_status="blocked", errors=[_error(error)])
        statuses = [entry["accuracy_status"] for entry in result["evaluations"]]
        result["accuracy_status"] = "measured" if "measured" in statuses else "blocked" if result["errors"] or "blocked" in statuses else "not_evaluated"
        result["measurement_completeness"] = "all_candidates_measured" if statuses and all(s == "measured" for s in statuses) else "partial_or_unmeasured"
        result["interpretation"] = ("Synthetic comparator/pipeline evidence only, not GIAB accuracy" if result["experiment_kind"] == "synthetic_smoke" else
                                    "HG002 v5.0q draft/work-in-progress benchmark, recommended successor to v4.2.1; chr20 is not WGS evidence" if result["scope"] == "chr20" else
                                    "HG002 v5.0q draft/work-in-progress benchmark, recommended successor to v4.2.1; measured is not accuracy_pass")
        result["ended_at"] = _now()
        save_json(output / "evaluation.json", result)
        report["evaluation_manifest"] = str(output / "evaluation.json")
        report["accuracy_status"] = result["accuracy_status"]
        for entry, (_, _, _, _, candidate, _) in zip(result["evaluations"], candidates):
            candidate["accuracy_status"] = entry["accuracy_status"]
        save_json(run_path, report)
        return result
