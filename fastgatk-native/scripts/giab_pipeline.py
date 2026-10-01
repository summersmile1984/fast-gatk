#!/usr/bin/env python3
"""Serial, provenance-gated HG002 pipelines; calling is not an accuracy verdict.

Only ``run_pipeline`` is public. The production path reads the complete input
manifest; an explicitly supplied synthetic_smoke manifest is an in-process
verification interface, never a production CLI override. Execution, biological
contracts and later truth scoring remain separate evidence channels.
"""
from __future__ import annotations

import fcntl
import json
import math
import os
import shlex
import statistics
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from giab_assets import load_lock
from giab_environment import PinnedTools, ROOT, save_json, sha256, source_fingerprint
from giab_preflight import SCOPES, _audit_sam, _capacity

STAGES = ("markdup", "bqsr", "apply", "hc", "genotype")
BINARIES = dict(zip(STAGES, ("mark-duplicates", "bqsr", "apply-bqsr", "hc-call", "genotype-gvcf")))
OUTPUT_KEYS = dict(zip(STAGES, ("markdup_bam", "recal_table", "recal_bam", "gvcf", "vcf")))
CHECKED = {"checked_match", "checked_differences"}
PROFILE = "region1m-locus-float32"
GIB = 1024 ** 3


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require(condition, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _load(path: Path) -> dict:
    with path.open() as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"Expected JSON object: {path}")
    return value


def _path(value, label: str) -> Path:
    _require(isinstance(value, str) and bool(value), f"Missing {label} path")
    path = Path(value)
    _require(path.is_absolute() and path.is_file(), f"Missing absolute {label} file: {path}")
    return path.resolve()


def _verify_digest(path: Path, expected: str, label: str) -> None:
    _require(isinstance(expected, str) and len(expected) == 64, f"Missing SHA-256 for {label}")
    _require(sha256(path) == expected, f"Changed {label}: {path}; repeat preflight explicitly")


def _verify_environment(environment: dict, data_root: Path) -> dict:
    _require(environment.get("status") == "eligible", "Execution environment preflight is not eligible")
    cpus = environment.get("cpus", [])
    _require(len(cpus) == 16 and all(type(cpu) is int and cpu >= 0 for cpu in cpus)
             and len(set(cpus)) == 16, "Environment must pin 16 distinct actual CPU IDs")
    _require(set(cpus) <= os.sched_getaffinity(0), "Pinned CPUs are no longer permitted to this process")
    _require(environment.get("memory_bytes") == 64 * GIB, "Environment must pin exactly 64 GiB")
    _require(environment.get("resource_status") == "enforced_within_limit",
             "Environment lacks observed CPU/memory/swap enforcement")
    image = environment.get("envelope", {}).get("image", {}).get("id", "")
    _require(image.startswith("sha256:"), "Execution envelope image is not content pinned")
    build = environment["build"]
    _require(_load(data_root / "build-manifest.json") == build,
             "Environment and build-manifest.json do not describe the same build")
    current = source_fingerprint()
    _require(current == build.get("source"), "Native source differs from the pinned build source manifest")
    build_dir = Path(build["build_directory"]).resolve()
    cache = build_dir / "CMakeCache.txt"
    _verify_digest(cache, build.get("cache_sha256"), "native CMake cache")
    cache_text = cache.read_text()
    for required in ("CMAKE_BUILD_TYPE:STRING=Release", "FASTGATK_KOKKOS_BACKEND:STRING=OPENMP",
                     "Kokkos_ENABLE_OPENMP:BOOL=ON"):
        _require(required in cache_text, f"Pinned build is not approved OpenMP Release: {required}")
    for tool in BINARIES.values():
        entry = build["binaries"][tool]
        binary = _path(entry.get("path"), tool)
        _require(binary.parent == build_dir and binary.name == "fastgatk-" + tool,
                 f"Binary is outside the single pinned build: {binary}")
        _verify_digest(binary, entry.get("sha256"), tool)
    oracle = environment["oracle"]
    for key, pinned in (("java", ROOT / "third_party/jdk17/bin/java"),
                        ("gatk", ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar")):
        path = _path(oracle.get(key), "pinned " + key)
        _require(path == pinned.resolve(), f"Unpinned oracle {key}: {path}")
        _verify_digest(path, oracle.get(key + "_sha256"), key)
    runtime = environment["runtime"]
    _require(runtime.get("artifacts"), "Missing pinned execution runtime artifacts")
    for name, entry in runtime["artifacts"].items():
        _verify_digest(_path(entry.get("path"), name), entry.get("sha256"), name)
    images = load_lock()["images"]
    for tool in ("samtools", "bcftools"):
        observed = environment.get("images", {}).get(tool, {})
        _require(observed.get("id") and observed.get("requested") == images[tool]
                 and any(digest.split("@")[-1] == images[tool].split("@")[-1]
                         for digest in observed.get("repo_digests", [])),
                 f"Missing observed pinned {tool} image digest")
    return {"status": "eligible", "checked_at": _now(), "source_sha256": current["sha256"],
            "build_manifest": str(data_root / "build-manifest.json"),
            "binary_sha256": {name: entry["sha256"] for name, entry in build["binaries"].items()},
            "policy": "Rehash existing pinned artifacts; never rebuild or choose a different backend"}


def _input_paths(manifest: dict, scope: str, synthetic: bool) -> dict:
    _require(manifest.get("sample") == "HG002", "Only single-sample HG002 is supported")
    info = manifest.get("scopes", {}).get(scope, {})
    _require(info.get("input_status") == "eligible", f"Requested scope {scope} input is not eligible: {info.get('blocking_gates', [])}")
    if not synthetic:
        _require(manifest.get("preflight_kind") == "complete_inputs", "Production requires complete-input preflight, not reference-only preparation")
        _require(manifest.get("experiment_kind", "giab") == "giab", "Synthetic manifests cannot be loaded by the production path")
        _require(info.get("reference_eligibility") == "eligible"
                 and info.get("truth", {}).get("status") == "eligible", "Requested scope truth/reference gate failed")
        gates = manifest.get("gates", {})
        required = ["reference", "dbsnp", "mills", "benchmark_bed", "truth_full_decode", "truth_" + scope,
                    "source_bam", "bam_reference", "provenance", "common_bam", "coverage", "aggregate_capacity"]
        if scope == "chr20":
            required.append("chr20_bam")
        failed = [name for name in required if gates.get(name, {}).get("status") != "eligible"]
        _require(not failed, f"Requested scope lacks complete input audit: {failed}")
        _require(gates["provenance"].get("quality_eligibility") == "eligible", "Production QUAL provenance is unknown or blocked")
    reference = _path(manifest.get("reference_path"), "reference")
    paths = {"reference": reference, "bam": _path(info.get("bam"), "scope BAM"),
             "scope_bed": _path(info.get("scope_bed"), "scope BED"),
             "evaluation_bed": _path(info.get("evaluation_bed"), "evaluation BED"),
             **{name: _path(manifest.get("known_sites", {}).get(name), name) for name in ("dbsnp", "mills")}}
    if not synthetic:
        expected_bam = manifest.get("common_bam") if scope == "autosomes" else manifest["gates"]["chr20_bam"].get("path")
        _require(paths["bam"] == _path(expected_bam, "audited common/pilot BAM"), "Scope BAM is not the audited shared input")
        fai = {}
        with Path(str(reference) + ".fai").open() as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                fai[fields[0]] = int(fields[1])
        with paths["scope_bed"].open() as stream:
            rows = [line.rstrip("\n").split("\t")[:3] for line in stream if line.strip() and not line.startswith("#")]
        _require(rows == [[chrom, "0", str(fai[chrom])] for chrom in SCOPES[scope]],
                 "Production calling scope must contain the entire declared chromosomes, without shrinking")
        _path(info["truth"].get("path"), "truth")
    return paths


def _outputs(directory: Path) -> dict:
    return {"markdup_bam": str(directory / "markdup.bam"),
            "markdup_metrics": str(directory / "markdup.metrics.txt"),
            "recal_table": str(directory / "recal.table"), "recal_bam": str(directory / "recal.bam"),
            "gvcf": str(directory / "calls.g.vcf.gz"), "vcf": str(directory / "calls.vcf.gz"),
            "hc_telemetry": str(directory / "hc.telemetry.json")}


def _commands(branch: str, outputs: dict, paths: dict, scratch: Path, scope: str,
              environment: dict, *, aggregate: bool = False) -> dict[str, list[str]]:
    p = {key: str(value) for key, value in paths.items()}
    intervals = [value for chrom in SCOPES[scope] for value in ("-L", chrom)]
    shared_bqsr = ["--known-sites", p["dbsnp"], "--known-sites", p["mills"],
                   "--mismatches-context-size", "2", "--indels-context-size", "3", "--maximum-cycle-value", "500"]
    apply_options = ["--bqsr-recal-file", outputs["recal_table"], "--preserve-qscores-less-than", "6",
                     "--quantize-quals", "0", "--use-original-qualities", "false", "--emit-original-quals", "false",
                     "--create-output-bam-index", "true"]
    hc_options = ["-I", outputs["recal_bam"], "-R", p["reference"], *intervals,
                  "-O", outputs["gvcf"], "-ERC", "GVCF", "--sample-ploidy", "2",
                  "--native-pair-hmm-use-double-precision", "false", "--create-output-variant-index", "true"]
    gt_options = ["-R", p["reference"], "-V", outputs["gvcf"], "-O", outputs["vcf"],
                  "--sample-ploidy", "2", "--standard-min-confidence-threshold-for-calling", "30",
                  "--max-alternate-alleles", "6", "--genotype-assignment-method", "PREFER_PLS",
                  "--use-new-qual-calculator", "true", "--use-posteriors-to-calculate-qual", "false",
                  "--create-output-variant-index", "true"]
    if branch == "native":
        binary = {stage: environment["build"]["binaries"][tool]["path"] for stage, tool in BINARIES.items()}
        return {
            "markdup": [binary["markdup"], "-I", p["bam"], "-O", outputs["markdup_bam"],
                        "--metrics-file", outputs["markdup_metrics"], "--tmp-dir", str(scratch / "markdup"),
                        "--max-records-in-memory", "100000", "--optical-duplicate-pixel-distance", "2500",
                        "--tagging-policy", "DontTag", "--create-output-bam-index", "true"],
            "bqsr": [binary["bqsr"], "-I", outputs["markdup_bam"], "-R", p["reference"],
                     "-O", outputs["recal_table"], "--batch-records", "4096", *shared_bqsr],
            "apply": [binary["apply"], "-I", outputs["markdup_bam"], "-R", p["reference"],
                      "-O", outputs["recal_bam"], "--batch-records", "4096", *apply_options],
            "hc": [binary["hc"], *hc_options, "--threads", "16", "--batch-records", "4096",
                   *([] if aggregate else ["--stream-by-region", "1000000"]),
                   "--telemetry", outputs["hc_telemetry"]],
            "genotype": [binary["genotype"], *gt_options, "--gatk-compatible-annotations",
                         *([] if aggregate else ["--stream-by-locus"])]}
    oracle = environment["oracle"]
    def java(stage: str, tool: str) -> list[str]:
        return [oracle["java"], "-XX:ActiveProcessorCount=16", "-Xmx40g",
                "-Djava.io.tmpdir=" + str(scratch / stage), "-jar", oracle["gatk"], tool]
    return {
        "markdup": [*java("markdup", "MarkDuplicates"), "-I", p["bam"], "-O", outputs["markdup_bam"],
                    "-M", outputs["markdup_metrics"], "--TMP_DIR", str(scratch / "markdup"),
                    "--MAX_RECORDS_IN_RAM", "100000", "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", "2500",
                    "--TAGGING_POLICY", "DontTag", "--REMOVE_DUPLICATES", "false",
                    "--REMOVE_SEQUENCING_DUPLICATES", "false", "--CLEAR_DT", "true", "--CREATE_INDEX", "true"],
        "bqsr": [*java("bqsr", "BaseRecalibrator"), "-I", outputs["markdup_bam"], "-R", p["reference"],
                 "-O", outputs["recal_table"], *shared_bqsr, "--use-original-qualities", "false"],
        "apply": [*java("apply", "ApplyBQSR"), "-I", outputs["markdup_bam"], "-R", p["reference"],
                  "-O", outputs["recal_bam"], *apply_options],
        "hc": [*java("hc", "HaplotypeCaller"), *hc_options, "--native-pair-hmm-threads", "16"],
        "genotype": [*java("genotype", "GenotypeGVCFs"), *gt_options]}


def _artifacts(directory: Path) -> list[dict]:
    # Include partial files/sidecars, not scratch, decoded audit SAMs or logs.
    return [{"path": str(path), "bytes": path.stat().st_size} for path in sorted(directory.iterdir())
            if path.is_file() and not path.is_symlink()]


def _resource_status(stages: list[dict]) -> str:
    measured = [s.get("resource_status", "not_measured") for s in stages if s.get("attempted")]
    for value in ("limit_hit", "unenforced", "not_measured"):
        if value in measured:
            return value
    return "enforced_within_limit" if measured else "not_measured"


def _sum_metric(stages: list[dict], name: str) -> float | None:
    values = [entry.get(name) for entry in stages]
    if not values or any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in values):
        return None
    return sum(values)


def _run_branch(branch: str, directory: Path, paths: dict, scope: str, environment: dict,
                data_root: Path, run_output: Path, tools, save, result: dict, *, stages=STAGES, aggregate=False) -> None:
    from giab_contracts import audit_stage
    from giab_execution import execute_stage

    scratch = directory / "scratch"
    for stage in stages:
        (scratch / stage).mkdir(parents=True, exist_ok=True)
        (directory / "stages" / stage).mkdir(parents=True, exist_ok=True)
        (directory / "contracts" / stage).mkdir(parents=True, exist_ok=True)
    outputs = result["outputs"]
    commands = _commands(branch, outputs, paths, scratch, scope, environment, aggregate=aggregate)
    stage_inputs = {"markdup": paths["bam"], "bqsr": Path(outputs["markdup_bam"]),
                    "apply": Path(outputs["markdup_bam"]), "hc": Path(outputs["recal_bam"]),
                    "genotype": Path(outputs["gvcf"])}
    result.update(started_at=_now(), execution_status="failed", contract_status="not_evaluated",
                  resource_status="not_measured", accuracy_status="not_evaluated", stages={
                      stage: {"execution_status": "not_run", "contract_status": "not_evaluated",
                              "resource_status": "not_measured", "attempted": False,
                              "argv": commands[stage], "reproduction_command": shlex.join(commands[stage]),
                              "output_path": outputs[OUTPUT_KEYS[stage]]} for stage in stages})
    began = time.monotonic()
    save()
    blocked_by = None
    for stage in stages:
        entry = result["stages"][stage]
        stage_dir = directory / "stages" / stage
        if blocked_by:
            entry.update(blocked_by=blocked_by, reason="Earlier stage failed; no downstream output is reusable")
            save_json(stage_dir / "stage.json", entry)
            save()
            continue
        entry.update(attempted=True, started_at=_now())
        save()
        try:
            execution = execute_stage(commands[stage], data_root=data_root, stage_dir=stage_dir,
                                      scratch=scratch / stage, environment=environment,
                                      timeout_seconds=(6 if scope == "chr20" else 24) * 3600,
                                      name=directory.name + "-" + stage, mounts=[run_output])
            entry.update(execution)
            entry["attempted"] = True
            entry["reproduction"] = {"workload_argv": commands[stage], "environment": entry.get("env", {}),
                                     "docker_argv": entry.get("docker_argv"),
                                     "request_path": str(stage_dir / "request.json"),
                                     "note": "Use the saved Docker envelope, not the workload command alone, to retain limits and pinned runtime"}
            if entry.get("execution_status") != "completed":
                blocked_by = stage
                result["execution_status"] = entry.get("execution_status", "failed")
            else:
                audit = audit_stage(stage, branch, input_path=stage_inputs[stage],
                                    output_path=Path(outputs[OUTPUT_KEYS[stage]]), reference=paths["reference"],
                                    scope_bed=paths["scope_bed"], work=directory / "contracts" / stage, tools=tools)
                entry["contract"] = audit
                entry["contract_status"] = audit.get("contract_status", "failed")
                if entry["contract_status"] not in CHECKED:
                    entry["contract_status"] = "failed"
                    result["contract_status"] = "failed"
                    blocked_by = stage
                elif entry.get("resource_status") != "enforced_within_limit":
                    entry["reason"] = "Workload completed without demonstrated resource enforcement"
                    result["failure_kind"] = "resource_enforcement"
                    blocked_by = stage
        except Exception as error:
            entry.update(error=str(error), error_type=type(error).__name__)
            if entry.get("execution_status") == "completed":
                entry["contract_status"] = "failed"
                result["contract_status"] = "failed"
            else:
                entry["execution_status"] = "failed"
            blocked_by = stage
        entry.setdefault("ended_at", _now())
        entry["artifacts"] = _artifacts(directory)
        save_json(stage_dir / "stage.json", entry)
        result["artifacts"] = entry["artifacts"]
        save()
    measured = [entry for entry in result["stages"].values() if entry.get("attempted")]
    if not blocked_by:
        result["execution_status"] = "completed"
        result["contract_status"] = ("checked_differences" if any(entry["contract_status"] == "checked_differences"
                                                                 for entry in measured) else "checked_match")
    result.update(ended_at=_now(), wall_seconds_including_contracts=time.monotonic() - began,
                  failed_stage=blocked_by, partial=bool(blocked_by),
                  resource_status=_resource_status(measured),
                  pipeline_seconds=_sum_metric(measured, "seconds"),
                  cpu_seconds=_sum_metric(measured, "cpu_seconds"),
                  metrics_cover_complete_chain=not bool(blocked_by), artifacts=_artifacts(directory))
    result["scoring_candidate"] = result["execution_status"] == "completed" and result["contract_status"] in CHECKED
    save()


def _capacity_for_run(manifest: dict, scope: str, synthetic: bool, paths: dict, work: Path, tools) -> dict:
    if not synthetic:
        capacity = manifest["gates"]["aggregate_capacity"]
        for relative, digest in capacity.get("source", {}).items():
            _verify_digest(ROOT / relative, digest["sha256"], "aggregate capacity source " + relative)
        entry = capacity.get("scopes", {}).get(scope, {})
        _require(entry.get("status") in ("within_uint32_capacity", "not_run_capacity_limit"),
                 "Aggregate payload capacity was not audited for requested scope")
        return entry
    work.mkdir(parents=True, exist_ok=True)
    decoded = work / "input.sam"
    tools("samtools", ["view", "--no-PG", "-h", str(paths["bam"])], "synthetic-capacity", stdout=decoded)
    audit = _audit_sam(decoded, work, "capacity")
    _require(audit["status"] == "eligible", "Synthetic aggregate input SAM audit failed")
    capacity = _capacity(audit, (scope,))
    save_json(work / "capacity.json", capacity)
    decoded.unlink()
    return {**capacity["scopes"][scope], "evidence": str(work / "capacity.json"),
            "full_decode_sam_removed_after_success": True}


def _summaries(repetitions: list[dict], synthetic: bool, scope: str) -> dict:
    result = {"interpretation": "synthetic_smoke_not_GIAB_evidence" if synthetic else
              ("paired_chr20_local_preprocessing" if scope == "chr20" else "single_observation_no_stability_claim"),
              "first_sample_reported_separately": True, "branches": {}, "paired_samples": []}
    for branch in ("native", "gatk"):
        samples = []
        for rep in repetitions:
            observation = rep["branches"][branch]
            samples.append({"repetition": rep["number"], "sample_class": rep["sample_class"],
                            "execution_status": observation.get("execution_status"),
                            "contract_status": observation.get("contract_status"),
                            "complete": observation.get("scoring_candidate", False),
                            **{key: observation.get(key) for key in
                               ("pipeline_seconds", "cpu_seconds", "wall_seconds_including_contracts")}})
        complete = [sample for sample in samples if sample["complete"]]
        medians = {}
        if len(complete) > 1:
            for key in ("pipeline_seconds", "cpu_seconds", "wall_seconds_including_contracts"):
                values = [sample[key] for sample in complete]
                medians[key] = statistics.median(values) if all(isinstance(v, (float, int)) and math.isfinite(v) for v in values) else None
        result["branches"][branch] = {"raw_samples": samples, "first_sample": samples[0] if samples else None,
                                      "complete_sample_count": len(complete), "medians": medians,
                                      "median_population": "all completed samples, including the separately labelled first sample",
                                      "subsequent_samples": samples[1:]}
    for rep in repetitions:
        native, gatk = (rep["branches"][name] for name in ("native", "gatk"))
        if native.get("scoring_candidate") and gatk.get("scoring_candidate"):
            ntime, gtime = native.get("pipeline_seconds"), gatk.get("pipeline_seconds")
            result["paired_samples"].append({"repetition": rep["number"], "branch_order": rep["branch_order"],
                                             "native_seconds": ntime, "gatk_seconds": gtime,
                                             "observed_gatk_over_native": gtime / ntime if ntime and gtime is not None else None})
    return result


def run_pipeline(data_root: Path, output: Path, *, scope="chr20", profile=PROFILE,
                 cpus=16, memory_gib=64, prepared_manifest: dict | None = None,
                 repetitions: int | None = None) -> dict:
    """Run independent serial native/GATK chains and persist every failure.

    Output must be fresh. Production repetitions are fixed to three chr20 pairs
    or one autosome pair. A supplied synthetic_smoke manifest defaults to one
    pair and may request additional pairs; it never inherits GIAB eligibility.
    Missing prerequisites return a blocked run.json rather than starting work.
    A concurrent/reused output is rejected without overwriting its evidence.
    """
    data_root, output = Path(data_root).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".run.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        _require(not any(path.name != ".run.lock" for path in output.iterdir()),
                 f"Output is not fresh; choose a new directory: {output}")
        synthetic = prepared_manifest is not None
        report = {"schema_version": 1, "data_root": str(data_root), "output": str(output),
                  "sample": "HG002", "scope": scope, "experiment_kind": "synthetic_smoke" if synthetic else "giab",
                  "profile": {"name": profile, "configuration_name": profile + "-t16-m64g", "cpus": cpus,
                              "memory_gib": memory_gib, "hc_stream_region_bases": 1000000,
                              "genotype_stream_by_locus": True, "native_pair_hmm_use_double_precision": False,
                              "batch_records": 4096, "optical_duplicate_pixel_distance": 2500,
                              "markdup_max_records_in_memory": 100000, "java_heap_gib": 40},
                  "started_at": _now(), "input_status": "pending", "execution_status": "failed",
                  "contract_status": "not_evaluated", "resource_status": "not_measured",
                  "accuracy_status": "not_evaluated", "environment_manifest": str(data_root / "environment-manifest.json"),
                  "input_manifest": str(output / "input-manifest.json"), "repetitions": [], "auxiliary": [], "errors": []}
        report["reproduction"] = {
            "api": "giab_pipeline.run_pipeline",
            "kwargs": {"data_root": str(data_root), "output": str(output), "scope": scope, "profile": profile,
                       "cpus": cpus, "memory_gib": memory_gib, "repetitions": repetitions},
            "prepared_manifest_snapshot": str(output / "input-manifest.json") if synthetic else None,
            "policy": "Explicit reruns must choose a fresh output directory; no automatic retry or partial reuse"}
        if not synthetic:
            command = [sys.executable, str(Path(__file__).with_name("benchmark_giab.py")), "run",
                       "--data-root", str(data_root), "--scope", scope, "--profile", profile,
                       "--cpus", str(cpus), "--memory-gib", str(memory_gib), "--output", str(output)]
            report["reproduction"]["argv"] = command
            report["reproduction"]["command"] = shlex.join(command)
        def save():
            save_json(output / "run.json", report)
        save()
        began = time.monotonic()
        try:
            _require(scope in SCOPES, "Only chr20 and autosomes scopes are supported")
            _require(profile == PROFILE and type(cpus) is int and cpus == 16
                     and type(memory_gib) is int and memory_gib == 64,
                     "Approved profile is region1m-locus-float32 with exactly 16 CPUs and 64 GiB")
            count = (1 if synthetic or scope == "autosomes" else 3) if repetitions is None else repetitions
            _require(type(count) is int and count > 0, "Repetitions must be a positive integer")
            _require(synthetic or count == (3 if scope == "chr20" else 1), "Production repetition count cannot be changed")
            report["requested_repetitions"] = count
            manifest = json.loads(json.dumps(prepared_manifest)) if synthetic else _load(data_root / "input-manifest.json")
            save_json(output / "input-manifest.json", manifest)
            _require(not synthetic or manifest.get("experiment_kind") == "synthetic_smoke",
                     "prepared_manifest is reserved exclusively for explicitly labelled synthetic_smoke")
            paths = _input_paths(manifest, scope, synthetic)
            report["input_status"] = "eligible"
            report["input_eligibility_boundary"] = ("Synthetic biological eligibility is supplied by the caller; no formal truth REF gate or GIAB scoring claim" if synthetic else
                                                     "Only requested scope gates apply; global status of unrelated scopes is not a veto")
            report["preprocessing_domain"] = ("full shared WGS BAM; chromosome restriction only at HC" if scope == "autosomes" else
                                                 "chr20 local preprocessing, without interchromosomal mates; not WGS MarkDuplicates/BQSR")
            report["paths"] = {key: str(path) for key, path in paths.items()}
            environment = _load(data_root / "environment-manifest.json")
            save_json(output / "environment-manifest.json", environment)
            report["environment_manifest_snapshot"] = str(output / "environment-manifest.json")
            if "build" in environment:
                save_json(output / "build-manifest.json", environment["build"])
            report["provenance_check"] = _verify_environment(environment, data_root)
            report["profile"]["cpu_ids"] = environment["cpus"]
            audit_tools = PinnedTools(data_root, load_lock()["images"], environment["cpus"], environment["memory_bytes"], mounts=[output])
            audit_tools.logs = output / "tool-logs"
            audit_tools.scratch = output / "scratch" / "contracts"
            audit_tools.scratch.mkdir(parents=True, exist_ok=True)
            tool_counter = 0
            run_id = uuid.uuid4().hex[:12]
            def tools(tool, args, name, *, stdout=None):
                nonlocal tool_counter
                tool_counter += 1
                return audit_tools(tool, args, f"{run_id}-{tool_counter:06d}-{name}", stdout=stdout)
            from giab_contracts import compare_branches
            for number in range(1, count + 1):
                order = ["native", "gatk"] if number % 2 else ["gatk", "native"]
                rep_dir = output / f"rep-{number:02d}"
                repetition = {"number": number, "sample_class": "first_sample" if number == 1 else "subsequent_sample",
                              "branch_order": order, "branches": {}, "comparisons": {"contract_status": "not_evaluated"}}
                for branch in ("native", "gatk"):
                    directory = rep_dir / branch
                    directory.mkdir(parents=True)
                    repetition["branches"][branch] = {"branch": branch, "directory": str(directory), "outputs": _outputs(directory),
                                                       "execution_status": "not_run", "contract_status": "not_evaluated"}
                report["repetitions"].append(repetition)
                save()
                for branch in order:
                    result = repetition["branches"][branch]
                    try:
                        _run_branch(branch, Path(result["directory"]), paths, scope, environment, data_root, output, tools, save, result)
                    except Exception as error:
                        result.update(execution_status="failed", contract_status="failed", error=str(error),
                                      error_type=type(error).__name__, partial=True, scoring_candidate=False)
                        save()
                native, gatk = (repetition["branches"][branch] for branch in ("native", "gatk"))
                if all(branch.get("scoring_candidate") for branch in (native, gatk)):
                    try:
                        work = rep_dir / "comparison"
                        work.mkdir()
                        repetition["comparisons"] = compare_branches(native, gatk, reference=paths["reference"],
                                                                      scope_bed=paths["scope_bed"], work=work, tools=tools)
                    except Exception as error:
                        repetition["comparisons"] = {"contract_status": "failed", "error": str(error), "error_type": type(error).__name__}
                else:
                    repetition["comparisons"]["reason"] = "At least one complete branch is unavailable; independent successful output remains a scoring candidate"
                save()
            auxiliary = {"name": "native-aggregate", "branch": "native", "source_repetition": 1,
                         "profile": "aggregate-float32-t16-m64g", "execution_status": "not_run",
                         "contract_status": "not_evaluated", "resource_status": "not_measured", "accuracy_status": "not_evaluated"}
            report["auxiliary"].append(auxiliary)
            save()
            try:
                capacity = _capacity_for_run(manifest, scope, synthetic, paths, output / "aggregate-capacity", tools)
                auxiliary["capacity"] = capacity
                if capacity["status"] == "not_run_capacity_limit":
                    auxiliary.update(execution_status="not_run_capacity_limit", reason="Audited uint32 aggregate payload capacity exceeded; more RAM cannot fix offsets")
                elif scope != "chr20":
                    auxiliary["reason"] = "Aggregate auxiliary calling is predeclared for chr20 only, not the default autosome path"
                else:
                    source = report["repetitions"][0]["branches"]["native"]
                    applied = source.get("stages", {}).get("apply", {})
                    if applied.get("execution_status") != "completed" or applied.get("contract_status") not in CHECKED or applied.get("resource_status") != "enforced_within_limit":
                        auxiliary["reason"] = "First native recal BAM is unavailable or failed its contract/resource gate"
                    else:
                        directory = output / "auxiliary" / "native-aggregate"
                        directory.mkdir(parents=True)
                        auxiliary.update(directory=str(directory), outputs=_outputs(directory),
                                         shared_native_recal_bam=source["outputs"]["recal_bam"])
                        auxiliary["outputs"]["recal_bam"] = source["outputs"]["recal_bam"]
                        _run_branch("native", directory, paths, scope, environment, data_root, output, tools, save,
                                    auxiliary, stages=("hc", "genotype"), aggregate=True)
                        for unused in ("markdup_bam", "markdup_metrics", "recal_table"):
                            auxiliary["outputs"].pop(unused)
            except Exception as error:
                auxiliary.update(execution_status="failed", reason=str(error), error_type=type(error).__name__)
            branches = [branch for rep in report["repetitions"] for branch in rep["branches"].values()]
            report["execution_status"] = "completed" if all(b.get("scoring_candidate") for b in branches) else "failed"
            contracts = [b.get("contract_status") for b in branches] + [r["comparisons"].get("contract_status") for r in report["repetitions"]]
            report["contract_status"] = ("failed" if "failed" in contracts else "not_evaluated" if "not_evaluated" in contracts else
                                          "checked_differences" if "checked_differences" in contracts else "checked_match")
            report["resource_status"] = _resource_status([stage for branch in branches for stage in branch.get("stages", {}).values()])
            report["measurement_summary"] = _summaries(report["repetitions"], synthetic, scope)
            report["auxiliary_status_is_independent"] = True
        except Exception as error:
            if report["input_status"] == "pending":
                report["input_status"] = "blocked"
            report["errors"].append({"error": str(error), "error_type": type(error).__name__})
            report["blocked"] = not bool(report["repetitions"])
            report["execution_status"] = "failed"
        finally:
            report.update(ended_at=_now(), wall_seconds_including_gates_contracts_auxiliary=time.monotonic() - began)
            report["accuracy_status"] = "blocked" if report["input_status"] != "eligible" else "not_evaluated"
            report["accuracy_interpretation"] = "No biological pass/fail threshold; this runner never measures accuracy"
            save()
        return report
