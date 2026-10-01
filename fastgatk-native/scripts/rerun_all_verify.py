#!/usr/bin/env python3
"""Re-run the per-tool verify_*.py scripts and persist fixed-format reports.

For each of the 48 production tools in `fastgatk-native/CMakeLists.txt` (the
`fastgatk-genomicsdb-export` bridge binary has no verify script and is
explicitly skipped), this driver:

  1. enumerates the `verify_*.py` scripts whose filename matches the tool's
     native binary stem (longest-composite-first, simple prefixes last);
  2. runs every script as a subprocess (timeout 600s, parallel up to
     ``min(8, ncpu)`` workers), capturing stdout/stderr/exit/elapsed;
  3. writes ``<tool>-rerun-<YYYYMMDD>.md`` + ``<tool>-rerun-<YYYYMMDD>.json``
     under ``<repo>/fastgatk-native/evidence/<YYYY-MM-DD>-rerun/<tool>/``;
  4. aggregates everything into ``INDEX.md`` at the rerun root.

The driver is intentionally read-only with respect to source/build artefacts:
it never touches ``CMakeLists.txt``, the build tree, or any other scripts.
It only writes inside the rerun evidence directory.

Exit code: 0 when every per-tool report was written successfully (individual
verify-script failures are recorded in JSON, not propagated as driver errors).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import resource
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[2].resolve()
SCRIPTS_DIR = REPO / "fastgatk-native/scripts"
EVIDENCE_ROOT = REPO / "fastgatk-native/evidence"

# -----------------------------------------------------------------------------
# Tool → script mapping.  Mirrors the reasoning captured in
# /tmp/tool_mapping_final.txt for offline reproduction.
#
# Rules (longest composite first, simple prefix `bqsr` last):
#   apply_vqsr, apply_bqsr, analyze_covariates, gather_bqsr_reports, gather_bqsr,
#   gather_pileup, get_pileup, pileup_summaries, calculate_contamination,
#   annotate_intervals, count_bases_in_reference, compare_references,
#   check_reference_compatibility, fasta_alternate, fasta_reference, shift_fasta,
#   index_feature_file, count_reads, flag_stat, split_intervals, filter_intervals,
#   preprocess_intervals, depth_of_coverage, gather_tranches, variants_to_table,
#   variant_eval, validate_variants, variant_recalibrator, variant_filtration,
#   left_align, call_copy_ratio, model_segments, create_read_count_panel,
#   reblock, genotype_gvcf, combine_gvcfs, collect_f1r2, collect_allelic,
#   collect_read_counts, denoise_read_counts, learn_read, filter_mutect, mutect2,
#   genomicsdb, gather_vcfs, select_variants, mark_duplicates, sort_sam,
#   somatic, pairhmm, fragment_aggregation, malformed_input_fail_loud,
#   pon_gatk_bundled_byte_structure, cli_parity, gvcf_stream_overlapping_indels,
#   resource_limits, overlapping_quality_correction, activity_region,
#   real_assembly_graph, vqsr_scatter, cnv_somatic_e2e,
#   hdf5_simple_count_collection, bam_intervals.
# Simple prefixes:  hc_ → HaplotypeCaller, bqsr → BaseRecalibrator (last).
# -----------------------------------------------------------------------------

COMPOSITE_MATCHES: list[tuple[str, str]] = [
    ("apply_vqsr", "ApplyVQSR"),
    ("apply_bqsr", "ApplyBQSR"),
    ("analyze_covariates", "AnalyzeCovariates"),
    ("gather_bqsr_reports", "GatherBQSRReports"),
    ("gather_bqsr", "GatherBQSRReports"),
    ("gather_pileup", "GatherPileupSummaries"),
    ("get_pileup", "GetPileupSummaries"),
    ("pileup_summaries", "GetPileupSummaries"),
    ("calculate_contamination", "CalculateContamination"),
    ("annotate_intervals", "AnnotateIntervals"),
    ("count_bases_in_reference", "CountBasesInReference"),
    ("compare_references", "CompareReferences"),
    ("check_reference_compatibility", "CheckReferenceCompatibility"),
    ("fasta_alternate", "FastaAlternateReferenceMaker"),
    ("fasta_reference", "FastaReferenceMaker"),
    ("shift_fasta", "ShiftFasta"),
    ("index_feature_file", "IndexFeatureFile"),
    ("count_reads", "CountReads"),
    ("flag_stat", "FlagStat"),
    ("split_intervals", "SplitIntervals"),
    ("filter_intervals", "FilterIntervals"),
    ("preprocess_intervals", "PreprocessIntervals"),
    ("depth_of_coverage", "DepthOfCoverage"),
    ("gather_tranches", "GatherTranches"),
    ("variants_to_table", "VariantsToTable"),
    ("variant_eval", "VariantEval"),
    ("validate_variants", "ValidateVariants"),
    ("funcotate_segments", "FuncotateSegments"),
    ("filter_funcotations", "FilterFuncotations"),
    ("funcotator", "Funcotator"),
    ("variant_annotator", "VariantAnnotator"),
    ("variant_recalibrator", "VariantRecalibrator"),
    ("variant_filtration", "VariantFiltration"),
    ("left_align", "LeftAlignAndTrimVariants"),
    ("call_copy_ratio", "CallCopyRatioSegments"),
    ("model_segments", "ModelSegments"),
    ("create_read_count_panel", "CreateReadCountPanelOfNormals"),
    ("reblock", "ReblockGVCF"),
    ("genotype_gvcf", "GenotypeGVCFs"),
    ("combine_gvcfs", "CombineGVCFs"),
    ("collect_f1r2", "CollectF1R2Counts"),
    ("collect_allelic", "CollectAllelicCounts"),
    ("collect_read_counts", "CollectReadCounts"),
    ("denoise_read_counts", "DenoiseReadCounts"),
    ("learn_read", "LearnReadOrientationModel"),
    ("filter_mutect", "FilterMutectCalls"),
    ("mutect2", "Mutect2"),
    ("genomicsdb", "GenomicsDBImport"),
    ("gather_vcfs", "GatherVcfs"),
    ("select_variants", "SelectVariants"),
    ("mark_duplicates", "MarkDuplicates"),
    ("sort_sam", "SortSam"),
    ("somatic", "Mutect2"),
    ("pairhmm", "Mutect2"),
    ("fragment_aggregation", "Mutect2"),
    ("malformed_input_fail_loud", "VariantsToTable"),
    ("pon_gatk_bundled_byte_structure", "CreateReadCountPanelOfNormals"),
    ("cli_parity", "CountReads"),
    ("gvcf_stream_overlapping_indels", "HaplotypeCaller"),
    ("resource_limits", "HaplotypeCaller"),
    ("overlapping_quality_correction", "HaplotypeCaller"),
    ("activity_region", "HaplotypeCaller"),
    ("real_assembly_graph", "HaplotypeCaller"),
    ("vqsr_scatter", "VariantRecalibrator"),
    ("cnv_somatic_e2e", "Mutect2"),
    ("hdf5_simple_count_collection", "CollectReadCounts"),
    ("bam_intervals", "HaplotypeCaller"),
]

SIMPLE_MATCHES: list[tuple[str, str]] = [
    ("hc_", "HaplotypeCaller"),
    ("bqsr", "BaseRecalibrator"),
]

CANONICAL_TO_NATIVE_BINARY = {
    "HaplotypeCaller": "fastgatk-hc-call",
    "Mutect2": "fastgatk-mutect2",
    "GenomicsDBImport": "fastgatk-genomicsdb-import",
    "BaseRecalibrator": "fastgatk-bqsr",
    "ApplyBQSR": "fastgatk-apply-bqsr",
    "GatherBQSRReports": "fastgatk-gather-bqsr-reports",
    "AnalyzeCovariates": "fastgatk-analyze-covariates",
    "VariantsToTable": "fastgatk-variants-to-table",
    "VariantEval": "fastgatk-variant-eval",
    "ValidateVariants": "fastgatk-validate-variants",
    "GetPileupSummaries": "fastgatk-get-pileup-summaries",
    "CalculateContamination": "fastgatk-calculate-contamination",
    "GatherPileupSummaries": "fastgatk-gather-pileup-summaries",
    "LearnReadOrientationModel": "fastgatk-learn-read-orientation-model",
    "DenoiseReadCounts": "fastgatk-denoise-read-counts",
    "CreateReadCountPanelOfNormals": "fastgatk-create-read-count-panel-of-normals",
    "CallCopyRatioSegments": "fastgatk-call-copy-ratio-segments",
    "ModelSegments": "fastgatk-model-segments",
    "GatherTranches": "fastgatk-gather-tranches",
    "GenotypeGVCFs": "fastgatk-genotype-gvcf",
    "ReblockGVCF": "fastgatk-reblock-gvcf",
    "SelectVariants": "fastgatk-select-variants",
    "GatherVcfs": "fastgatk-gather-vcfs",
    "CombineGVCFs": "fastgatk-combine-gvcfs",
    "FilterMutectCalls": "fastgatk-filter-mutect-calls",
    "LeftAlignAndTrimVariants": "fastgatk-left-align-trim",
    "VariantFiltration": "fastgatk-variant-filtration",
    "ApplyVQSR": "fastgatk-apply-vqsr",
    "Funcotator": "fastgatk-funcotator",
    "VariantAnnotator": "fastgatk-variant-annotator",
    "VariantRecalibrator": "fastgatk-variant-recalibrator",
    "SortSam": "fastgatk-sort-sam",
    "MarkDuplicates": "fastgatk-mark-duplicates",
    "AnnotateIntervals": "fastgatk-annotate-intervals",
    "CountBasesInReference": "fastgatk-count-bases-in-reference",
    "CompareReferences": "fastgatk-compare-references",
    "CheckReferenceCompatibility": "fastgatk-check-reference-compatibility",
    "FastaReferenceMaker": "fastgatk-fasta-reference-maker",
    "FastaAlternateReferenceMaker": "fastgatk-fasta-alternate-reference-maker",
    "ShiftFasta": "fastgatk-shift-fasta",
    "IndexFeatureFile": "fastgatk-index-feature-file",
    "CountReads": "fastgatk-count-reads",
    "FlagStat": "fastgatk-flag-stat",
    "SplitIntervals": "fastgatk-split-intervals",
    "FilterIntervals": "fastgatk-filter-intervals",
    "PreprocessIntervals": "fastgatk-preprocess-intervals",
    "DepthOfCoverage": "fastgatk-depth-of-coverage",
    "CollectReadCounts": "fastgatk-collect-read-counts",
    "CollectF1R2Counts": "fastgatk-collect-f1r2-counts",
    "CollectAllelicCounts": "fastgatk-collect-allelic-counts",
}

# Script names that this driver EXCLUDES from per-tool reports (infra / meta).
# They are tracked in /tmp/tool_mapping_final.txt but not assigned to any tool.
EXCLUDED_INFRA_SCRIPTS = {
    "verify_kokkos_api_boundary.py",
    "verify_kokkos_backend_config.py",
    "verify_kokkos_backend_matrix.py",
    "verify_kernel_benchmark.py",
    "verify_benchmark_timing.py",
    "verify_giab_assets.py",
    "verify_giab_evaluation.py",
    "verify_giab_preflight.py",
    "verify_cli_alignment_doc.py",
    "verify_scope_contract.py",
    "verify_scheduler_retry_contract.py",
    "verify_slurm_wrapper.py",
    "verify_fixture_digests.py",
    "verify_launcher_contract.py",
    "verify_remote_staging.py",
    "verify_optional_boolean_contract.py",
    "verify_native.py",
    "verify_flow_hmer_option.py",
    "verify_gatk_oracle.py",
    "verify_wgs_benchmark_a6_oracle.py",
    "verify_flow_pairhmm_real_data_oracle.py",
    "verify_indel.py",
}

TOOL_NOTE = "(native-only bridge, no verify script)"


# ----------------------------------------------------------------------------- helpers

def classify(stem: str) -> str | None:
    """Return the canonical GATK Java tool name for a verify-script stem, or None.

    The function is total: every stem in the corpus maps to a known tool.
    Excluded infra scripts return ``None`` and are skipped silently.
    """
    if f"verify_{stem}.py" in EXCLUDED_INFRA_SCRIPTS:
        return None
    for needle, canon in COMPOSITE_MATCHES:
        if needle in stem:
            return canon
    # Simple prefix matching (with `bqsr` last).
    if stem.startswith("hc_"):
        return "HaplotypeCaller"
    if "bqsr" in stem:
        return "BaseRecalibrator"
    return None


def tool_slug(canonical: str) -> str:
    """Evidence-directory slug = native-binary stem, e.g. 'hc-call'."""
    return CANONICAL_TO_NATIVE_BINARY[canonical].replace("fastgatk-", "")


def build_tool_to_scripts() -> dict[str, list[str]]:
    """Group every verify_*.py under its canonical tool (alphabetised)."""
    out: dict[str, list[str]] = {}
    for script in sorted(SCRIPTS_DIR.glob("verify_*.py")):
        stem = script.name[len("verify_"):-len(".py")]
        canon = classify(stem)
        if canon is None:
            continue
        out.setdefault(canon, []).append(script.name)
    return out


def ncpu() -> int:
    """Return a sensible worker count, falling back to 1."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, (os.cpu_count() or 1))


# ----------------------------------------------------------------------------- runner

TAIL_LEN = 400

# Polling interval (seconds) for /proc/<pid>/status VmRSS/VmHWM.  5ms is
# the sweet spot: a 50-100ms native-binary invocation (typical of fastgatk-*
# synthetic-fixture runs) is reliably caught (10-20 polls per lifetime),
# while sampling overhead stays < 0.3% even for 30-minute scripts.  50ms
# polling missed many short-lived native children; 10ms still missed some
# in 2-way parallel mode under contention from disk-backed kernel cache
# traffic.
RSS_POLL_SECONDS = 0.005

# Process-classification: by /proc/<pid>/comm we decide whether a descendant
# is a native binary (fastgatk-*) or the GATK Java VM (comm == "java").  We
# carry per-class peak RSS and observed wall-clock so the rerun report can
# compare native-vs-Java resource use on the same workload.
NATIVE_COMM_PREFIXES = ("fastgatk-",)
JAVA_COMM_NAMES = {"java"}


def classify_comm(comm: str) -> str:
    """Return ``'native'``, ``'java'`` or ``'other'`` for a /proc/<pid>/comm."""
    comm = comm.strip()
    if not comm:
        return "other"
    if comm in JAVA_COMM_NAMES:
        return "java"
    if any(comm.startswith(p) for p in NATIVE_COMM_PREFIXES):
        return "native"
    return "other"


def read_proc_status(pid: int) -> tuple[int | None, int | None]:
    """Return (VmRSS_kB, VmHWM_kB) for ``pid`` from /proc/<pid>/status.

    Both fields may be ``None`` if the process has already exited or the
    status file is unreadable.
    """
    try:
        with open(f"/proc/{pid}/status", "r", encoding="utf-8") as handle:
            rss = None
            hwm = None
            for line in handle:
                if line.startswith("VmRSS:"):
                    rss = int(line.split()[1])
                elif line.startswith("VmHWM:"):
                    hwm = int(line.split()[1])
            return rss, hwm
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return None, None


def read_proc_comm(pid: int) -> str:
    try:
        with open(f"/proc/{pid}/comm", "r", encoding="utf-8") as handle:
            return handle.read()
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return ""


def read_proc_cmdline(pid: int) -> str:
    try:
        with open(f"/proc/{pid}/cmdline", "r", encoding="utf-8") as handle:
            raw = handle.read().replace("\x00", " ").strip()
            return raw
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return ""


def read_proc_stat(pid: int) -> tuple[int | None, int | None, int | None]:
    """Return (utime_ticks, stime_ticks, starttime_ticks) from /proc/<pid>/stat.

    All three are None when the file is unreadable.
    """
    try:
        with open(f"/proc/{pid}/stat", "r", encoding="utf-8") as handle:
            content = handle.read()
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return None, None, None
    # The comm field (2nd token) is wrapped in parens and may contain spaces
    # or parens; find the last ')' to split safely.
    close = content.rfind(")")
    if close == -1:
        return None, None, None
    fields = content[close + 1:].split()
    # After ')', field indices start at 3 (relative to /proc/PID/stat manual):
    # 3=state, 4=ppid, ... 14=utime, 15=stime, 22=starttime.
    # Because we dropped the comm field, the indices shift by 2.
    # new_index 14 (utime) -> fields[11], 15 (stime) -> fields[12],
    # 22 (starttime) -> fields[19].
    try:
        utime = int(fields[11])
        stime = int(fields[12])
        starttime = int(fields[19])
    except (IndexError, ValueError):
        return None, None, None
    return utime, stime, starttime


def descendants_of(root_pid: int) -> set[int]:
    """Walk the live process tree of ``root_pid`` via ``/proc/<pid>/task/<tid>/children``.

    The ``children`` file lists the *current* immediate children of the named
    thread.  This is updated synchronously by the kernel on fork/clone/exit
    and is the authoritative source for "what does pid X own right now".
    Critically, it does not require reading every ``/proc/<pid>/stat`` entry
    on the system, which is slow on a many-core machine with 1000+ procs.

    BFS from ``root_pid`` following the children file at each level.  Pids
    that exist but are already dead (we hit them via /proc/<pid>/stat in
    a previous iteration) are NOT visible here -- they will be picked up
    in the next poll's `/proc/*/stat` fallback.

    Returns the set of all pids seen, excluding ``root_pid``.
    """
    seen: set[int] = set()
    frontier: list[int] = [root_pid]
    while frontier:
        current = frontier.pop()
        # Each process has a main thread (tid == pid).  For thread-group-leader
        # processes, /proc/<pid>/task/<pid>/children lists all immediate
        # children of the group; for non-leader threads the file does not
        # exist.  We only need the leader's children.
        children_path = f"/proc/{current}/task/{current}/children"
        try:
            with open(children_path, "r", encoding="utf-8") as handle:
                tokens = handle.read().split()
            for tok in tokens:
                try:
                    child = int(tok)
                except ValueError:
                    continue
                if child not in seen and child != root_pid:
                    seen.add(child)
                    frontier.append(child)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            # Process exited between check and read -- fall back to /proc/<pid>/stat
            try:
                with open(f"/proc/{current}/stat", "r", encoding="utf-8") as handle:
                    content = handle.read()
                close = content.rfind(")")
                if close != -1:
                    fields = content[close + 1:].split()
                    if len(fields) >= 2:
                        # ppid = fields[1]; we want children, but if process
                        # is gone, /proc entry also gone.
                        pass
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                pass
    return seen


def sample_process_tree_loop(root_pid: int, samples: dict,
                            stop_event: threading.Event,
                            clk_tck: int = os.sysconf("SC_CLK_TCK")) -> None:
    """Background sampler: walk the process tree of ``root_pid`` periodically
    and record per-class peak RSS + per-process wall-clock.

    The ``samples`` dict is mutated in place.  Final layout::

        {
            "max_rss_kb": <int>,           # max across all descendants (root excluded)
            "max_hwm_kb": <int>,           # max across all descendants
            "native": {
                "peak_rss_kb": <int>,
                "peak_hwm_kb": <int>,
                "wall_clock_seconds_sum": <float>,  # sum of per-process elapsed
                "process_count": <int>,            # distinct native pids encountered
                "cmdlines": [str, ...]             # up to 8 distinct cmdlines
            },
            "java": {...},
            "other": {...},
            "processes_seen": <int>,       # total distinct pids visited
        }
    """
    max_rss = 0
    max_hwm = 0
    per_class: dict[str, dict] = {
        cls: {
            "peak_rss_kb": 0,
            "peak_hwm_kb": 0,
            "wall_clock_seconds_sum": 0.0,
            "process_count": 0,
            "cmdlines": [],
        }
        for cls in ("native", "java", "other")
    }
    seen_pids: set[int] = set()
    seen_cmdlines: set[str] = set()
    per_pid_first_seen: dict[int, float] = {}
    per_pid_class: dict[int, str] = {}

    while not stop_event.is_set():
        try:
            descendants = descendants_of(root_pid)
        except Exception:  # noqa: BLE001
            descendants = set()
        now = time.monotonic()
        live_pids = set()
        for pid in descendants:
            comm = read_proc_comm(pid)
            cls = classify_comm(comm)
            rss, hwm = read_proc_status(pid)
            cmdline = read_proc_cmdline(pid)
            live_pids.add(pid)
            if pid not in seen_pids:
                seen_pids.add(pid)
                per_pid_first_seen[pid] = now
                per_pid_class[pid] = cls
                if cmdline and cmdline not in seen_cmdlines and len(seen_cmdlines) < 32:
                    seen_cmdlines.add(cmdline)
                    bucket = per_class[cls]["cmdlines"]
                    if len(bucket) < 8:
                        bucket.append(cmdline[:200])
            # Per-class RSS peaks (use the classification we cached when
            # the pid first appeared; comm can change after exec but the
            # process identity does not).
            cls_now = per_pid_class.get(pid, cls)
            if rss is not None:
                if rss > max_rss:
                    max_rss = rss
                if rss > per_class[cls_now]["peak_rss_kb"]:
                    per_class[cls_now]["peak_rss_kb"] = rss
            if hwm is not None:
                if hwm > max_hwm:
                    max_hwm = hwm
                if hwm > per_class[cls_now]["peak_hwm_kb"]:
                    per_class[cls_now]["peak_hwm_kb"] = hwm
        # Charge wall-clock to each pid's class when it disappears.
        for pid in list(per_pid_first_seen):
            if pid not in live_pids:
                first = per_pid_first_seen.pop(pid)
                cls = per_pid_class.pop(pid)
                elapsed = now - first
                per_class[cls]["wall_clock_seconds_sum"] += elapsed
                per_class[cls]["process_count"] += 1
        stop_event.wait(RSS_POLL_SECONDS)

    samples["max_rss_kb"] = max_rss
    samples["max_hwm_kb"] = max_hwm
    samples["native"] = per_class["native"]
    samples["java"] = per_class["java"]
    samples["other"] = per_class["other"]
    samples["processes_seen"] = len(seen_pids)


def run_one(script_path: Path, repo: Path) -> dict:
    """Run a single verify_*.py; return a per-script result dict.

    Tracks both wall-clock time and process resource use:

    * ``elapsed_seconds`` / ``wall_clock_seconds`` — driver-side wall time.
    * ``peak_rss_bytes`` — peak resident set size observed for the verify
      script's Python process (via background polling of
      ``/proc/<pid>/status`` ``VmRSS``).  Captures native + Java + Python
      interpreter footprint while the script ran.
    * ``peak_vm_hwm_bytes`` — peak high-water-mark virtual size (``VmHWM``).
    * ``cpu_time_seconds`` — user + system CPU time recorded by
      ``resource.getrusage(RUSAGE_CHILDREN)`` after the subprocess exits.
      Aggregates native + java subprocesses spawned by the script.
    * ``exit_code``, ``timed_out``, ``stdout_tail``, ``stderr_tail``,
      ``passed``, ``skipped`` — as before.
    """
    started = time.monotonic()
    perf_start = time.perf_counter()
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    env.setdefault("FASTGATK_NATIVE_BUILD", str(repo / "fastgatk-native/build"))
    env["FASTGATK_REQUIRE_GATK_ORACLE"] = "0"  # oracle-missing does NOT fail the batch
    name = script_path.name
    timed_out = False
    exit_code = -1
    stdout = ""
    stderr = ""
    tree_samples: dict = {
        "max_rss_kb": 0, "max_hwm_kb": 0,
        "native": {"peak_rss_kb": 0, "peak_hwm_kb": 0, "wall_clock_seconds_sum": 0.0,
                   "process_count": 0, "cmdlines": []},
        "java":   {"peak_rss_kb": 0, "peak_hwm_kb": 0, "wall_clock_seconds_sum": 0.0,
                   "process_count": 0, "cmdlines": []},
        "other":  {"peak_rss_kb": 0, "peak_hwm_kb": 0, "wall_clock_seconds_sum": 0.0,
                   "process_count": 0, "cmdlines": []},
        "processes_seen": 0,
    }
    proc: subprocess.Popen | None = None
    sampler_stop = threading.Event()
    sampler_thread: threading.Thread | None = None

    try:
        proc = subprocess.Popen(
            [sys.executable, str(script_path)],
            cwd=str(repo),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        # Walk the script's entire process tree (Python + native fastgatk-*
        # + GATK java + helpers) and classify by /proc/<pid>/comm.
        sampler_thread = threading.Thread(
            target=sample_process_tree_loop,
            args=(proc.pid, tree_samples, sampler_stop),
            daemon=True,
        )
        sampler_thread.start()
        try:
            stdout, stderr = proc.communicate(timeout=600)
            exit_code = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            proc.kill()
            try:
                stdout, stderr = proc.communicate(timeout=10)
            except Exception:  # noqa: BLE001
                stdout, stderr = "", ""
            exit_code = -1
    except Exception as exc:  # noqa: BLE001
        # Popen / sample setup failed; degrade gracefully.
        stderr = (stderr or "") + f"\n[driver error: {exc!r}]\n"
        exit_code = -1
    finally:
        sampler_stop.set()
        if sampler_thread is not None:
            sampler_thread.join(timeout=2.0)
    elapsed = time.monotonic() - started
    wall_clock = time.perf_counter() - perf_start
    script_passed = exit_code == 0
    blob = (stderr + stdout).lower()
    # A script is "skipped" only when it explicitly took an oracle-not-
    # verified path, signalled by either:
    #   (a) the "NOT VERIFIED AGAINST GATK" banner that oracle_guard
    #       prints to stderr,
    #   (b) a {"status": "skip", ...} JSON line on stdout (used by
    #       scripts that detect the missing oracle and report skip via
    #       structured JSON instead of stderr),
    #   (c) an exit code of 77.
    oracle_skip = (
        "not verified against gatk" in blob
        or '"status": "skip"' in blob
        or '"status":"skip"' in blob
        or exit_code == 77
    )

    # Aggregate CPU time across all child processes the verify script
    # spawned (including native and Java subprocesses).  RUSAGE_CHILDREN
    # survives the death of the immediate child.
    try:
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_time = float(ru.ru_utime) + float(ru.ru_stime)
        max_rss_kb_children = int(ru.ru_maxrss)  # kilobytes on Linux
    except Exception:  # noqa: BLE001
        cpu_time = 0.0
        max_rss_kb_children = 0

    # Process-tree samples: per-class peak RSS + observed wall-clock sums.
    native_cls = tree_samples.get("native", {})
    java_cls = tree_samples.get("java", {})
    other_cls = tree_samples.get("other", {})

    def cls_metrics(cls_dict: dict) -> dict:
        return {
            "peak_rss_bytes": int(cls_dict.get("peak_rss_kb", 0)) * 1024,
            "peak_vm_hwm_bytes": int(cls_dict.get("peak_hwm_kb", 0)) * 1024,
            "wall_clock_seconds_sum": round(cls_dict.get("wall_clock_seconds_sum", 0.0), 3),
            "process_count": int(cls_dict.get("process_count", 0)),
            "sample_cmdlines": list(cls_dict.get("cmdlines", [])),
        }

    return {
        "name": name,
        "exit_code": exit_code,
        "elapsed_seconds": round(elapsed, 3),
        "wall_clock_seconds": round(wall_clock, 3),
        "peak_rss_bytes": int(tree_samples.get("max_rss_kb", 0)) * 1024,
        "peak_vm_hwm_bytes": int(tree_samples.get("max_hwm_kb", 0)) * 1024,
        "cpu_time_seconds": round(cpu_time, 3),
        "max_rss_kb_children_rusage": max_rss_kb_children,
        "processes_seen": int(tree_samples.get("processes_seen", 0)),
        "native": cls_metrics(native_cls),
        "java": cls_metrics(java_cls),
        "other": cls_metrics(other_cls),
        "stdout_tail": stdout[-TAIL_LEN:],
        "stderr_tail": stderr[-TAIL_LEN:],
        "passed": script_passed,
        "timed_out": timed_out,
        "skipped": oracle_skip,
        "resource_sampled": True,  # set by driver; older backfilled sidecars may have None/false
    }


def run_tool(tool_slug_name: str, scripts: list[str], workers: int) -> list[dict]:
    """Run every script for one tool, in parallel, and return per-script records."""
    records: list[dict] = []
    with ThreadPoolExecutor(max_workers=max(1, min(workers, len(scripts)))) as executor:
        futures = {
            executor.submit(run_one, SCRIPTS_DIR / s, REPO): s for s in scripts
        }
        for fut in as_completed(futures):
            try:
                records.append(fut.result())
            except Exception as exc:  # noqa: BLE001
                records.append({
                    "name": futures[fut],
                    "exit_code": -1,
                    "elapsed_seconds": 0.0,
                    "stdout_tail": "",
                    "stderr_tail": f"driver error: {exc!r}",
                    "passed": False,
                    "timed_out": False,
                    "skipped": False,
                })
    # Keep insertion order deterministic by name.
    records.sort(key=lambda r: r["name"])
    return records


# ----------------------------------------------------------------------------- reporters

def write_reports(canonical_tool: str, scripts: list[str], records: list[dict],
                  rerun_root: Path, report_date: str, ymd: str) -> dict:
    """Write per-tool Markdown and JSON sidecars. Returns a summary dict."""
    slug = tool_slug(canonical_tool)
    native = CANONICAL_TO_NATIVE_BINARY[canonical_tool]
    tool_dir = rerun_root / slug
    tool_dir.mkdir(parents=True, exist_ok=True)
    md_path = tool_dir / f"{slug}-rerun-{ymd}.md"
    json_path = tool_dir / f"{slug}-rerun-{ymd}.json"

    total = len(records)
    passed = sum(1 for r in records if r["passed"])
    failed = sum(1 for r in records if not r["passed"])
    skipped = sum(1 for r in records if r.get("skipped"))
    elapsed_total = round(sum(r["elapsed_seconds"] for r in records), 3)

    json_payload = {
        "schema_version": 1,
        "tool": canonical_tool,
        "native_binary": native,
        "report_date": report_date,
        "report_kind": "verify-rerun",
        "scripts_total": total,
        "scripts_passed": passed,
        "scripts_failed": failed,
        "scripts_skipped": skipped,
        "total_elapsed_seconds": elapsed_total,
        "scripts": records,
    }
    json_path.write_text(json.dumps(json_payload, indent=2, sort_keys=False) + "\n",
                         encoding="utf-8")

    md_lines = [
        f"# {canonical_tool} rerun report — {report_date}",
        "",
        "## Summary",
        f"- Total scripts: {total}",
        f"- Passed: {passed}    Failed: {failed}    Skipped (exit 77 / oracle-guard skip): {skipped}",
        f"- Total elapsed: {elapsed_total}s",
        "",
        "## Per-script results",
        "",
        "| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in records:
        verdict = "✅" if r["passed"] else ("⏱" if r["timed_out"] else "❌")
        cpu_t = r.get("cpu_time_seconds", 0.0)
        native_metrics = r.get("native", {})
        java_metrics = r.get("java", {})
        native_rss = native_metrics.get("peak_rss_bytes", 0)
        java_rss = java_metrics.get("peak_rss_bytes", 0)
        java_wallclock = java_metrics.get("wall_clock_seconds_sum", 0.0)
        md_lines.append(
            f"| `{r['name']}` | {r['exit_code']} | {r['elapsed_seconds']:.3f}s | "
            f"{cpu_t:.2f}s | "
            f"{native_rss / (1024 * 1024):.1f} MiB | "
            f"{java_rss / (1024 * 1024):.1f} MiB | "
            f"{java_wallclock:.2f}s | {verdict} |"
        )
    md_lines.extend([
        "",
        "## Re-run command",
        "",
        "```bash",
        f"python3 fastgatk-native/scripts/rerun_all_verify.py --tool {slug}",
        "```",
        "",
        "## JSON sidecar",
        "",
        f"See `{json_path.name}` for full stdout/stderr tails.",
        "",
    ])
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    return {
        "tool": canonical_tool,
        "slug": slug,
        "native_binary": native,
        "scripts_total": total,
        "scripts_passed": passed,
        "scripts_failed": failed,
        "scripts_skipped": skipped,
        "elapsed_seconds": elapsed_total,
        "scripts": [r["name"] for r in records],
        "report_md": str(md_path.relative_to(REPO)),
        "report_json": str(json_path.relative_to(REPO)),
    }


def write_index(rerun_root: Path, summary_rows: list[dict], report_date: str, ymd: str) -> None:
    """Aggregate all per-tool summaries into ``INDEX.md``."""
    total_scripts = sum(r["scripts_total"] for r in summary_rows)
    total_passed = sum(r["scripts_passed"] for r in summary_rows)
    total_failed = sum(r["scripts_failed"] for r in summary_rows)
    total_skipped = sum(r["scripts_skipped"] for r in summary_rows)
    elapsed_total = round(sum(r["elapsed_seconds"] for r in summary_rows), 3)

    lines = [
        f"# Verify-script re-run index — {report_date}",
        "",
        f"- Tools covered: **{len(summary_rows)}** (excluding `genomicsdb-export` {TOOL_NOTE})",
        f"- Scripts total: **{total_scripts}**",
        f"- Passed: **{total_passed}**    Failed: **{total_failed}**    Skipped: **{total_skipped}**",
        f"- Total elapsed: **{elapsed_total}s**",
        "",
        "Per-tool summaries (cross-ref with `fastgatk-native/docs/cli-alignment.md` "
        "*Bit-identical / bounded-parity evidence* table):",
        "",
        "| Tool | Native binary | Scripts | Passed | Failed | Skipped | Elapsed (s) | Report |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in summary_rows:
        lines.append(
            f"| `{r['tool']}` | `{r['native_binary']}` | {r['scripts_total']} | "
            f"{r['scripts_passed']} | {r['scripts_failed']} | {r['scripts_skipped']} | "
            f"{r['elapsed_seconds']} | "
            f"[{Path(r['report_md']).name}]({r['report_md']}) |"
        )
    lines.extend([
        "",
        "## Bridge binary with no verify script",
        "",
        f"- `fastgatk-genomicsdb-export` — {TOOL_NOTE}; cross-reference the docs.",
        "",
        "## GATK-alignment status",
        "",
        "Run `aggregate_alignment_status.py` after a re-run to cross-reference",
        "this index with the *Bit-identical / bounded-parity evidence* table",
        "in `fastgatk-native/docs/cli-alignment.md`.  It writes",
        "`ALIGNMENT_STATUS.md` next to this file.",
        "",
        "## Native-vs-Java resource comparison",
        "",
        "Run `compare_native_vs_java.py` to classify every descendant process",
        "of the rerun into `native` (comm starts with `fastgatk-`) vs `java`",
        "(comm == `java`) buckets, and emit `NATIVE_VS_JAVA.md` next to this",
        "file.  Anomalies (native slower than Java, or native RSS exceeding",
        "Java's) are auto-flagged in that report.",
        "",
        "## Workflow documentation",
        "",
        "Full workflow, schema details, and an anomaly follow-up playbook",
        "live in `fastgatk-native/docs/regression-evidence.md`.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "python3 fastgatk-native/scripts/rerun_all_verify.py --repo .",
        "```",
        "",
        "```bash",
        "python3 fastgatk-native/scripts/verify_rerun_report.py --repo .",
        "```",
        "",
        "```bash",
        "python3 fastgatk-native/scripts/aggregate_alignment_status.py --repo .",
        "```",
        "",
        "```bash",
        "python3 fastgatk-native/scripts/compare_native_vs_java.py --repo .",
        "```",
        "",
    ])
    (rerun_root / "INDEX.md").write_text("\n".join(lines), encoding="utf-8")


# ----------------------------------------------------------------------------- main

def main() -> int:
    today = dt.date.today()
    ymd = today.strftime("%Y%m%d")
    report_date = today.strftime("%Y-%m-%d")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=str(REPO),
                        help="Path to repository root (default: %(default)s)")
    parser.add_argument("--tool", default=None,
                        help="Run only the named tool slug (default: all 48)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Per-tool parallel workers (default: min(8, ncpu))")
    parser.add_argument("--rerun-root", default=None,
                        help=f"Override the evidence/<date>-rerun path (default: <date>=today)")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="Skip running scripts; rebuild INDEX.md from existing JSON sidecars")
    args = parser.parse_args()

    repo = Path(args.repo).resolve()
    if not (repo / "fastgatk-native/scripts").is_dir():
        raise SystemExit(f"no scripts dir at {repo}/fastgatk-native/scripts")
    if rerun_root_str := args.rerun_root:
        rerun_root = Path(rerun_root_str).resolve()
    else:
        rerun_root = EVIDENCE_ROOT / f"{report_date}-rerun"
    rerun_root.mkdir(parents=True, exist_ok=True)

    workers = args.workers if args.workers else min(8, ncpu())

    tool_map = build_tool_to_scripts()
    if not tool_map:
        raise SystemExit("no verify_*.py scripts found; refusing to write empty report")

    selected: dict[str, list[str]] = {}
    for canonical in sorted(tool_map):
        slug = tool_slug(canonical)
        if args.tool and args.tool != slug:
            continue
        selected[canonical] = tool_map[canonical]

    if args.tool and not selected:
        raise SystemExit(f"unknown --tool {args.tool!r}; no matching scripts in mapping")

    summary_rows: list[dict] = []
    print(f"# rerun_all_verify.py — {report_date}")
    print(f"# tool scope: {len(selected)} tools, "
          f"workers per tool: {workers}, rerun root: {rerun_root}")
    overall_started = time.monotonic()
    if args.aggregate_only:
        # Re-aggregate from existing JSON sidecars; do NOT re-run scripts.
        for canonical in sorted(selected):
            slug = tool_slug(canonical)
            json_path = rerun_root / slug / f"{slug}-rerun-{ymd}.json"
            if not json_path.is_file():
                raise SystemExit(f"--aggregate-only: missing {json_path}")
            data = json.loads(json_path.read_text(encoding="utf-8"))
            summary_rows.append({
                "tool": data.get("tool", canonical),
                "slug": slug,
                "native_binary": data.get("native_binary", CANONICAL_TO_NATIVE_BINARY[canonical]),
                "scripts_total": data.get("scripts_total", 0),
                "scripts_passed": data.get("scripts_passed", 0),
                "scripts_failed": data.get("scripts_failed", 0),
                "scripts_skipped": data.get("scripts_skipped", 0),
                "elapsed_seconds": data.get("total_elapsed_seconds", 0.0),
                "scripts": [s["name"] for s in data.get("scripts", [])],
                "report_md": str((rerun_root / slug / f"{slug}-rerun-{ymd}.md").relative_to(REPO)),
                "report_json": str(json_path.relative_to(REPO)),
            })
            print(f"   {slug}: {summary_rows[-1]['scripts_passed']}/{summary_rows[-1]['scripts_total']} passed")
    else:
        for canonical in sorted(selected):
            scripts = selected[canonical]
            print(f"\n## {canonical} ({tool_slug(canonical)}) — {len(scripts)} scripts")
            records = run_tool(tool_slug(canonical), scripts, workers)
            row = write_reports(canonical, scripts, records, rerun_root,
                                report_date, ymd)
            summary_rows.append(row)
            print(f"   passed={row['scripts_passed']}  failed={row['scripts_failed']}  "
                  f"skipped={row['scripts_skipped']}  elapsed={row['elapsed_seconds']}s")

    write_index(rerun_root, summary_rows, report_date, ymd)
    overall_elapsed = round(time.monotonic() - overall_started, 3)
    print(f"\n# DONE in {overall_elapsed}s. INDEX: {rerun_root / 'INDEX.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
