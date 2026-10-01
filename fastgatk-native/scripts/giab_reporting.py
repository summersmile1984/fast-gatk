#!/usr/bin/env python3
"""Evidence-only GIAB reports: execution, contracts, resources and scores stay independent."""
from __future__ import annotations

import copy
import csv
import itertools
import json
import math
import re
import sqlite3
import statistics
import tempfile
from pathlib import Path

from giab_environment import save_json

STATUSES = ("input_status", "execution_status", "contract_status", "resource_status", "accuracy_status")
METRIC_KEYS = ("variant_type", "filter", "subset", "subtype", "genotype")
SCORES = ("truth_total", "truth_tp", "truth_fn", "query_total", "query_tp", "query_fp", "query_unk",
          "fp_gt", "fp_al", "precision", "recall", "f1", "confident_bases")
TIMINGS = ("pipeline_seconds", "cpu_seconds", "wall_seconds_including_contracts")
RESOURCES = ("seconds", "user_seconds", "system_seconds", "cpu_seconds", "max_rss_kb",
             "memory_peak_bytes", "cgroup_cpu_core_seconds", "cgroup_measurement_window_seconds",
             "allocated_cpu_utilization", "output_bytes", "scratch_peak_sampled_bytes", "io_stat_delta",
             "cpu_stat_delta", "memory_events_delta", "filesystem_inputs", "filesystem_outputs",
             "supervisor_wall_seconds")
ERROR_LIMITATION = (
    "Common truth FN means the same exact truth record identity. Query FP intersections use exact "
    "coordinates, alleles and genotype only; equivalent haplotypes with different representations "
    "may not intersect. These are diagnostic record overlaps, not haplotype-aware shared biological "
    "error counts. UNK remains separate from FP. Counts are unique record identities, not hap.py counts."
)


def _read_json(path: Path) -> dict:
    with path.open() as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _resolve(path: str, base: Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else base / value


def _number(value) -> bool:
    return type(value) in (int, float) and math.isfinite(value)


def _difference(native, gatk):
    return native - gatk if _number(native) and _number(gatk) else None


def _cell(value):
    # Explicit null, not an empty cell or zero, preserves missing-value semantics.
    if value is None or isinstance(value, (dict, list, bool)):
        return json.dumps(value, sort_keys=True, allow_nan=False)
    return value


def _write_tsv(path: Path, rows, fields: tuple | list) -> None:
    temporary = path.with_name(path.name + ".part")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _cell(row.get(key)) for key in fields})
    temporary.replace(path)


def _stratum(subset) -> dict:
    name = str(subset or "").lower()
    categories = []
    for category, tokens in (
        ("low_complexity", ("lowcomplex", "low_complex", "low-complex", "lcr")),
        ("homopolymer_tandem_repeat", ("homopoly", "tandem", "simplerepeat", "simple_repeat", "str_", "_str", "tr_")),
        ("low_mappability", ("mappab", "nonunique", "non_unique")),
        ("segmental_duplication", ("segdup", "segmental")),
        ("gc", ("gccontent", "gc_content", "gc_", "_gc")),
        ("coding_functional", ("coding", "cds", "exon", "functional", "difficultgene")),
        ("difficult_easy_union", ("difficult", "easy", "alldifficult", "all_difficult")),
    ):
        if any(token in name for token in tokens):
            categories.append(category)
    other_samples = sorted(set(re.findall(r"hg\d{3}", name)) - {"hg002"})
    sex_layer = bool(re.search(r"(?:^|[^a-z0-9])(?:chr)?[xy](?:[^a-z0-9]|$)|autosomesandsex|sexchrom", name))
    return {"highlight_categories": categories, "other_sample_layers": other_samples,
            "outside_autosome_biology": sex_layer,
            "interpretation": "not_HG002_autosome_biology" if other_samples or sex_layer else "overlapping_stratum_do_not_sum"}


def _primary(metrics: list[dict]) -> dict:
    result = {filter_name: {kind: None for kind in ("SNP", "INDEL")} for filter_name in ("PASS", "ALL")}
    for row in metrics:
        if all(row.get(key) == "*" for key in ("subset", "subtype", "genotype")):
            if row.get("filter") in result and row.get("variant_type") in result[row["filter"]]:
                result[row["filter"]][row["variant_type"]] = row
    return result


def _score_pairs(repetitions: list[dict], evaluations: dict) -> list[dict]:
    paired = []
    for repetition in repetitions:
        number = repetition["number"]
        native_id, gatk_id = f"rep-{number:02d}/native", f"rep-{number:02d}/gatk"
        native, gatk = evaluations.get(native_id, {}), evaluations.get(gatk_id, {})
        indexes = [{tuple(row.get(key) for key in METRIC_KEYS): row for row in item.get("metrics", [])}
                   for item in (native, gatk)]
        keys = set(indexes[0]) | set(indexes[1])
        keys.update((kind, filter_name, "*", "*", "*") for kind in ("SNP", "INDEL") for filter_name in ("PASS", "ALL"))
        for key in sorted(keys, key=lambda key: tuple(str(value) for value in key)):
            n, g = indexes[0].get(key), indexes[1].get(key)
            measured = native.get("accuracy_status") == gatk.get("accuracy_status") == "measured"
            paired.append({"repetition": number, "native_id": native_id, "gatk_id": gatk_id,
                           **dict(zip(METRIC_KEYS, key)), "native": n, "gatk": g,
                           "pair_status": "measured" if measured and n is not None and g is not None else "not_evaluated",
                           "native_minus_gatk": {field: _difference(n.get(field), g.get(field))
                                                 if measured and n is not None and g is not None else None for field in SCORES},
                           "stratum": _stratum(key[2])})
    return paired


def _paired_performance(repetitions: list[dict]) -> dict:
    samples = []
    for rep in repetitions:
        native, gatk = (rep.get("branches", {}).get(branch, {}) for branch in ("native", "gatk"))
        complete = all(item.get("scoring_candidate") and item.get("metrics_cover_complete_chain") for item in (native, gatk))
        samples.append({"repetition": rep.get("number"), "sample_class": rep.get("sample_class"),
                        "branch_order": rep.get("branch_order"), "complete_pair": complete,
                        "native": {key: native.get(key) for key in TIMINGS},
                        "gatk": {key: gatk.get(key) for key in TIMINGS},
                        "native_minus_gatk": {key: _difference(native.get(key), gatk.get(key)) if complete else None for key in TIMINGS}})
    complete = [sample for sample in samples if sample["complete_pair"]]
    medians = {}
    for key in TIMINGS:
        values = [sample["native_minus_gatk"][key] for sample in complete]
        medians[key] = statistics.median(values) if len(values) > 1 and all(_number(value) for value in values) else None
    return {"raw_samples": samples, "first_sample": samples[0] if samples else None,
            "complete_pair_count": len(complete), "median_native_minus_gatk": medians,
            "median_population": "All complete paired observations, including the separately labelled first round; never an inferred warmup",
            "interpretation": "single_observation_no_stability_claim" if len(complete) == 1 else
                              "no_complete_pairs" if not complete else "paired_observations_no_percentile_claim"}


def _error_identity(row: dict) -> tuple[str, tuple] | None:
    # Column spelling is the evaluator's public error TSV contract.
    side, decision = row.get("side"), row.get("error_class")
    category = "truth_FN" if side == "truth" and decision == "FN" else (
        "query_" + decision if side == "query" and decision in ("FP", "UNK") else None)
    if category is None:
        return None
    fields = ("chrom", "pos", "ref", "alt", "gt")
    if any(row.get(key) in (None, "") for key in fields):
        raise ValueError(f"Error record lacks exact identity fields {fields}")
    return category, tuple(row[key] for key in fields)


def _error_overlap(native: dict, gatk: dict, output: Path, base: Path, repetition: int) -> dict:
    result = {"repetition": repetition, "status": "not_evaluated", "limitation": ERROR_LIMITATION,
              "path": None, "counts": None, "sources": {branch: item.get("artifacts", {}).get("error_records_tsv")
                                                         for branch, item in (("native", native), ("gatk", gatk))}}
    if any(item.get("accuracy_status") != "measured" for item in (native, gatk)) or not all(result["sources"].values()):
        result["reason"] = "Both independently measured error-record artifacts are required; missing evidence is not an empty error set"
        return result
    path = output / f"error-overlap-rep-{repetition:02d}.tsv"
    try:
        # Disk-backed ordered merge bounds RAM even when a poor callset has millions of errors.
        with tempfile.TemporaryDirectory(prefix=".giab-error-overlap-", dir=output) as scratch:
            with sqlite3.connect(str(Path(scratch) / "records.sqlite")) as db:
                db.execute("CREATE TABLE records (category TEXT, identity TEXT, branch TEXT, payload TEXT, PRIMARY KEY(category, identity, branch, payload)) WITHOUT ROWID")
                for branch, source in result["sources"].items():
                    with _resolve(source, base).open(newline="") as stream:
                        reader = csv.DictReader(stream, delimiter="\t")
                        required = {"side", "error_class", "decision", "chrom", "pos", "ref", "alt", "gt"}
                        if not required.issubset(reader.fieldnames or []):
                            raise ValueError(f"Error TSV missing columns {sorted(required - set(reader.fieldnames or []))}: {source}")
                        for row in reader:
                            identity = _error_identity(row)
                            if identity is not None:
                                category, key = identity
                                db.execute("INSERT OR IGNORE INTO records VALUES (?, ?, ?, ?)",
                                           (category, json.dumps(key), branch, json.dumps(row, sort_keys=True)))
                    db.commit()
                counts = {category: {membership: 0 for membership in ("common", "native_only", "gatk_only")}
                          for category in ("truth_FN", "query_FP", "query_UNK")}

                def rows():
                    query = db.execute("SELECT category, identity, branch, payload FROM records ORDER BY category, identity, branch, payload")
                    for (category, identity), group in itertools.groupby(query, key=lambda row: row[:2]):
                        records = {"native": [], "gatk": []}
                        for _, _, branch, payload in group:
                            records[branch].append(json.loads(payload))
                        membership = "common" if all(records.values()) else "native_only" if records["native"] else "gatk_only"
                        counts[category][membership] += 1
                        chrom, pos, ref, alt, gt = json.loads(identity)
                        yield {"category": category, "membership": membership, "chrom": chrom, "pos": pos,
                               "ref": ref, "alt": alt, "gt": gt, "native_records": records["native"], "gatk_records": records["gatk"]}

                _write_tsv(path, rows(), ("category", "membership", "chrom", "pos", "ref", "alt", "gt", "native_records", "gatk_records"))
                result.update(status="compared_exact_records", path=str(path), counts=counts)
    except (OSError, ValueError, csv.Error, sqlite3.Error) as error:
        result.update(status="not_evaluated", error=str(error), error_type=type(error).__name__)
    return result


def _failure_records(run: dict, observations: list[dict], evaluation: dict | None, report_errors: list[dict]) -> list[dict]:
    failures = [{"location": "report", **error} for error in report_errors]
    if run.get("errors") or run.get("blocked") or run.get("execution_status") != "completed" or run.get("contract_status") == "failed":
        failures.append({"location": "run", **{key: run.get(key) for key in STATUSES},
                         "errors": run.get("errors", []), "blocked": run.get("blocked"), "reproduction": run.get("reproduction")})
    for rep in run.get("repetitions", []):
        comparison = rep.get("comparisons", {})
        if comparison.get("contract_status") == "failed":
            failures.append({"location": f"rep-{rep['number']:02d}/comparison", "details": comparison})
    for item in observations:
        if item.get("execution_status") != "completed" or item.get("contract_status") not in ("checked_match", "checked_differences") or item.get("resource_status") != "enforced_within_limit":
            failures.append({"location": item["id"], "details": {key: value for key, value in item.items() if key not in ("stages", "accuracy")}})
        for stage, entry in item.get("stages", {}).items():
            if entry.get("execution_status") != "completed" or entry.get("contract_status") not in ("checked_match", "checked_differences") or entry.get("resource_status") != "enforced_within_limit":
                failures.append({"location": f"{item['id']}/{stage}", "details": entry})
    if evaluation is not None:
        if evaluation.get("errors") or evaluation.get("accuracy_status") != "measured":
            failures.append({"location": "evaluation", "accuracy_status": evaluation.get("accuracy_status"), "errors": evaluation.get("errors", [])})
        for entry in evaluation.get("evaluations", []):
            if entry.get("accuracy_status") != "measured" or entry.get("errors"):
                failures.append({"location": "evaluation/" + entry["id"], "details": entry})
    return failures


def _table_rows(summary: dict):
    yield {"row_kind": "run_status", **{key: summary.get(key) for key in STATUSES}, "details": summary["caveats"]}
    for observation in summary["observations"]:
        common = {"id": observation["id"], "repetition": observation.get("repetition"), "branch": observation.get("branch"),
                  "configuration": observation.get("configuration"), "sample_class": observation.get("sample_class"),
                  **{key: observation.get(key) for key in STATUSES}}
        yield {**common, "row_kind": "auxiliary_resources" if observation["is_auxiliary"] else "branch_resources",
               **{key: observation.get(key) for key in TIMINGS}, "details": {key: observation.get(key) for key in
               ("metrics_cover_complete_chain", "performance_domain", "partial", "failed_stage", "reproduction")}}
        for stage, entry in observation.get("stages", {}).items():
            yield {**common, "row_kind": "stage_resources", "stage": stage,
                   **{key: entry.get(key) for key in STATUSES if key in entry},
                   **{key: entry.get(key) for key in RESOURCES}, "details": entry}
        evaluation = observation.get("accuracy") or {}
        metrics = evaluation.get("metrics", [])
        for metric in metrics:
            primary = all(metric.get(key) == "*" for key in ("subset", "subtype", "genotype"))
            yield {**common, "row_kind": "accuracy_main" if primary else "accuracy_stratum", **metric,
                   "details": {"stratum": _stratum(metric.get("subset")), "artifacts": evaluation.get("artifacts"),
                               "undefined_metrics": metric.get("undefined_metrics")}}
        for filter_name, variants in observation["main_accuracy"].items():
            for kind, metric in variants.items():
                if metric is None:
                    yield {**common, "row_kind": "accuracy_main", "filter": filter_name, "variant_type": kind,
                           "subset": "*", "subtype": "*", "genotype": "*", "metric_status": "not_evaluated",
                           "details": "No measured row available; counts and scores remain null"}
        if evaluation.get("execution"):
            execution = evaluation["execution"]
            yield {**common, "row_kind": "evaluator_resources", **{key: execution.get(key) for key in RESOURCES},
                   "execution_status": execution.get("execution_status"), "resource_status": execution.get("resource_status"),
                   "details": execution}
    for pair in summary["accuracy_differences"]:
        yield {"row_kind": "accuracy_native_minus_gatk", "repetition": pair["repetition"],
               **{key: pair.get(key) for key in METRIC_KEYS}, **pair["native_minus_gatk"],
               "metric_status": pair["pair_status"], "details": {"native": pair["native"], "gatk": pair["gatk"], "stratum": pair["stratum"]}}
    for sample in summary["paired_performance"]["raw_samples"]:
        yield {"row_kind": "performance_native_minus_gatk", "repetition": sample["repetition"],
               "sample_class": sample["sample_class"], **sample["native_minus_gatk"], "details": sample}
    yield {"row_kind": "performance_paired_median", **summary["paired_performance"]["median_native_minus_gatk"],
           "details": {key: value for key, value in summary["paired_performance"].items() if key not in ("raw_samples", "first_sample")}}
    measurements = summary.get("measurement_summary") or {}
    for branch, data in measurements.get("branches", {}).items():
        yield {"row_kind": "performance_branch_median", "branch": branch, **data.get("medians", {}),
               "details": {"complete_sample_count": data.get("complete_sample_count"), "median_population": data.get("median_population")}}
    for overlap in summary["error_overlaps"]:
        yield {"row_kind": "exact_error_overlap", "repetition": overlap["repetition"], "details": overlap}


def report_run(run: Path, output: Path | None = None) -> dict:
    """Write reports for complete, partial or blocked runs, without rerunning any workload."""
    run_path = Path(run).resolve()
    if run_path.is_dir():
        run_path /= "run.json"
    source = _read_json(run_path)
    run_dir = run_path.parent
    output = Path(output).resolve() if output is not None else run_dir
    output.mkdir(parents=True, exist_ok=True)
    report_errors = []
    evaluation = None
    evaluation_path = None
    if source.get("evaluation_manifest"):
        evaluation_path = _resolve(source["evaluation_manifest"], run_dir)
        try:
            evaluation = _read_json(evaluation_path)
            if evaluation.get("schema_version") != 1 or not isinstance(evaluation.get("evaluations"), list):
                raise ValueError("Unsupported or incomplete evaluation manifest schema")
            if evaluation.get("run_path") and _resolve(evaluation["run_path"], evaluation_path.parent).resolve() not in (run_path, run_dir):
                raise ValueError("Evaluation manifest belongs to a different run")
        except (OSError, ValueError) as error:
            evaluation = None
            report_errors.append({"artifact": str(evaluation_path), "error": str(error), "error_type": type(error).__name__})
    evaluations = {}
    for entry in (evaluation or {}).get("evaluations", []):
        if entry["id"] in evaluations:
            raise ValueError(f"Duplicate evaluation identity: {entry['id']}")
        evaluations[entry["id"]] = entry
    observations = []
    for repetition in source.get("repetitions", []):
        for branch, item in repetition.get("branches", {}).items():
            identity = f"rep-{repetition['number']:02d}/{branch}"
            observations.append({**copy.deepcopy(item), "id": identity, "branch": branch,
                                 "repetition": repetition["number"], "sample_class": repetition.get("sample_class"),
                                 "configuration": source.get("profile"), "is_auxiliary": False,
                                 "performance_domain": "five_stage_pipeline"})
    for item in source.get("auxiliary", []):
        name = item.get("name", "native-aggregate")
        identity = "auxiliary/" + name
        observations.append({**copy.deepcopy(item), "id": identity, "repetition": None,
                             "configuration": item.get("profile"), "is_auxiliary": True,
                             "performance_domain": "HC_and_Genotype_only_using_first_native_recal_BAM_not_a_five_stage_pipeline"})
    for observation in observations:
        score = evaluations.get(observation["id"])
        observation.setdefault("input_status", source.get("input_status", "pending"))
        observation.setdefault("resource_status", "not_measured")
        observation["accuracy_status"] = score.get("accuracy_status", "not_evaluated") if score else "not_evaluated"
        observation["accuracy"] = score
        observation["main_accuracy"] = _primary(score.get("metrics", []) if score else [])
    overlaps = []
    for rep in source.get("repetitions", []):
        number = rep["number"]
        overlaps.append(_error_overlap(evaluations.get(f"rep-{number:02d}/native", {}),
                                       evaluations.get(f"rep-{number:02d}/gatk", {}), output,
                                       evaluation_path.parent if evaluation_path else run_dir, number))
    for item in overlaps:
        if item.get("error"):
            report_errors.append({"location": "error_overlap", **item})
    caveats = [
        "GIAB HG002 v5.0q replaces v4.2.1 as the recommended benchmark but remains draft/work-in-progress; the truth release must be stated with every formal score.",
        "GATK is an independent implementation comparator, not biological truth; contract differences, accuracy and resource enforcement are separate conclusions.",
        "Measured accuracy is not an accuracy pass. No biological acceptance or noninferiority thresholds were declared.",
        "Main results are PASS; ALL and every available stratum are retained. Strata overlap and must not be summed. Other-sample and XY strata are not HG002 autosomal biological categories.",
        "TRUTH.TP and QUERY.TP have distinct evaluator semantics. BED-exterior query UNK is not FP. Missing metrics and zero-denominator ratios stay null.",
        "pipeline_seconds is the runner's sum of attempted stage workload seconds; partial chains are labelled. wall_seconds_including_contracts and whole-run wall are separate; evaluator resources are never added to pipeline sums.",
        "Process maximum RSS and cgroup peak have different meanings; cgroup peak includes descendants and charged page cache. I/O block counters are not bytes; scratch peaks are sampled, not exact maxima.",
        "First round is reported separately but not discarded as warmup. Medians use actual completed observations, no p95 or stability claim from a single pair.",
        "Aggregate auxiliary results reuse first-round native recal BAM and measure only HC/Genotype; they are independent of full-chain and streamed configuration success.",
        ERROR_LIMITATION,
    ]
    if source.get("experiment_kind") == "synthetic_smoke":
        caveats.insert(0, "SYNTHETIC SMOKE ONLY: all scores and resource observations are fixture evidence, not HG002/GIAB biological accuracy or production resource evidence.")
    if source.get("scope") == "chr20":
        caveats.append("CHR20 PILOT ONLY: local preprocessing excludes interchromosomal mates; it does not establish whole-genome MarkDuplicates/BQSR state or chr1–chr22 accuracy/performance.")
    elif source.get("scope") == "autosomes":
        caveats.append("Autosome baseline is a single observed pair, not a stable speedup or a percentile estimate. chrX/chrY are out of scope.")
    if evaluation is None:
        caveats.append("No readable attached evaluation exists. Accuracy is not_evaluated, not zero and not a pass; runner execution status remains unchanged.")
    failures = _failure_records(source, observations, evaluation, report_errors)
    summary = {**copy.deepcopy(source), "schema_version": 1, "report_schema_version": 1,
               "run_path": str(run_path), "report_output": str(output),
               "accuracy_status": evaluation.get("accuracy_status", "not_evaluated") if evaluation else "not_evaluated",
               "runner_accuracy_status": source.get("accuracy_status"), "observations": observations,
               "evaluation": evaluation, "evaluation_manifest": str(evaluation_path) if evaluation_path else None,
               "measurement_summary": source.get("measurement_summary"),
               "paired_performance": _paired_performance(source.get("repetitions", [])),
               "accuracy_differences": _score_pairs(source.get("repetitions", []), evaluations),
               "error_overlaps": overlaps, "caveats": caveats, "report_errors": report_errors,
               "failures": failures, "status_semantics": "Five independent dimensions; no accuracy_pass and no successful sibling masking a failed branch",
               "reproduction_policy": "Use the saved runner reproduction and per-stage Docker argv/request/environment to retain the pinned resource envelope. Choose a fresh run output; never reuse partial products.",
               "artifacts": {"summary_json": str(output / "giab_summary.json"), "summary_tsv": str(output / "giab_summary.tsv"),
                             "failure_details_json": str(output / "failure-details.json")}}
    fields = ("row_kind", "id", "repetition", "sample_class", "branch", "configuration", "stage", *STATUSES,
              *TIMINGS, *(key for key in RESOURCES if key not in TIMINGS), *METRIC_KEYS, *SCORES, "metric_status", "details")
    _write_tsv(output / "giab_summary.tsv", _table_rows(summary), fields)
    save_json(output / "failure-details.json", {"schema_version": 1, "run_path": str(run_path), "failures": failures,
                                               "reproduction": source.get("reproduction"), "policy": summary["reproduction_policy"]})
    save_json(output / "giab_summary.json", summary)
    return summary
