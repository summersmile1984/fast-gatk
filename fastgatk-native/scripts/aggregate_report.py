#!/usr/bin/env python3
"""Aggregate accuracy + efficiency + resource into a single comparison report.

This is the "run everything, produce one table" layer that the per-tool
oracles (accuracy) and per-tool ``benchmark_*.py`` (efficiency) never had.
It does NOT replace ``run_regression.sh``; it consumes its log for the
accuracy half and re-runs the cheap benchmark scripts for the efficiency half.

Sections
--------
* accuracy   -- parsed from the latest ``.diag/regression/<ts>/omp.log`` /
                 ``serial.log`` (per-test Passed/Failed -> per-tool counts).
* efficiency -- runs each tool's ``benchmark_*.py`` and reads its JSON
                 (p50/p95/warmup wall seconds, plus ``gatk_java`` when present).
* resource   -- ``benchmark_end_to_end.py`` (max RSS + FS read/write, HC/Mutect2
                 end-to-end).  Uniform per-tool RSS is not captured by the other
                 benchmark scripts (they time wall-clock only); see the note.

Output
------
``report/REPORT.md``      -- one markdown table (tool x accuracy/efficiency/resource)
``report/accuracy.json``  -- raw per-tool accuracy counts
``report/efficiency.json``-- raw per-tool benchmark JSON
``report/resources.json`` -- raw resource numbers

Usage
-----
    python3 fastgatk-native/scripts/aggregate_report.py                 # full
    python3 fastgatk-native/scripts/aggregate_report.py --skip-benchmarks
    python3 fastgatk-native/scripts/aggregate_report.py --regression-dir .diag/regression/20260914-115810
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "fastgatk-native" / "scripts"
BUILD = ROOT / "fastgatk-native" / "build"

# --------------------------------------------------------------------------
# Tool -> benchmark script (explicit, like compute_progress_score.TOOL_EVIDENCE)
# Several tools share one benchmark (BQSR/ApplyBQSR, the two FASTA makers,
# the read-metrics pair) and a few have none.
# --------------------------------------------------------------------------
TOOL_BENCHMARK: dict[str, str] = {
    "AnalyzeCovariates": "benchmark_analyze_covariates.py",
    "AnnotateIntervals": "benchmark_annotate_intervals.py",
    "ApplyBQSR": "benchmark_bqsr.py",
    "ApplyVQSR": "benchmark_apply_vqsr.py",
    "BaseRecalibrator": "benchmark_bqsr.py",
    "CalculateContamination": "benchmark_calculate_contamination.py",
    "CallCopyRatioSegments": "benchmark_call_copy_ratio_segments.py",
    "CheckReferenceCompatibility": "benchmark_check_reference_compatibility.py",
    "CollectAllelicCounts": "benchmark_collect_allelic_counts.py",
    "CollectF1R2Counts": "benchmark_collect_f1r2_counts.py",
    "CollectReadCounts": "benchmark_collect_read_counts.py",
    "CombineGVCFs": "benchmark_combine_gvcfs.py",
    "CompareReferences": "benchmark_compare_references.py",
    "CountBasesInReference": "benchmark_count_bases_in_reference.py",
    "CountReads": "benchmark_read_metrics.py",
    "CreateReadCountPanelOfNormals": "benchmark_create_read_count_panel_of_normals.py",
    "DenoiseReadCounts": "benchmark_denoise_read_counts.py",
    "DepthOfCoverage": "benchmark_depth_of_coverage.py",
    "FastaAlternateReferenceMaker": "benchmark_fasta_reference_tools.py",
    "FastaReferenceMaker": "benchmark_fasta_reference_tools.py",
    "FilterIntervals": "benchmark_filter_intervals.py",
    "FilterMutectCalls": "benchmark_filter_mutect_calls.py",
    "FlagStat": "benchmark_read_metrics.py",
    "GatherBQSRReports": "benchmark_gather_bqsr_reports.py",
    "GatherPileupSummaries": "benchmark_gather_pileup_summaries.py",
    "GatherTranches": "benchmark_gather_tranches.py",
    "GatherVcfs": "benchmark_gather_vcfs.py",
    "GenomicsDBImport": "benchmark_genomicsdb_import.py",
    "GenotypeGVCFs": "benchmark_genotype_gvcf.py",
    "GetPileupSummaries": "benchmark_get_pileup_summaries.py",
    "HaplotypeCaller": "benchmark_hc_broad.py",
    "IndexFeatureFile": "benchmark_index_feature_file.py",
    "LearnReadOrientationModel": "benchmark_learn_read_orientation_model.py",
    "LeftAlignAndTrimVariants": "benchmark_left_align.py",
    "MarkDuplicates": "benchmark_mark_duplicates.py",
    "ModelSegments": "benchmark_model_segments.py",
    "Mutect2": "benchmark_mutect2.py",
    "PreprocessIntervals": "benchmark_preprocess_intervals.py",
    "ReblockGVCF": "benchmark_reblock_gvcf.py",
    "SelectVariants": "benchmark_select_variants.py",
    "ShiftFasta": "benchmark_shift_fasta.py",
    "SortSam": "benchmark_sort_sam.py",
    "SplitIntervals": "benchmark_split_intervals.py",
    "ValidateVariants": "benchmark_validate_variants.py",
    "VariantEval": "benchmark_variant_eval.py",
    "VariantFiltration": "benchmark_variant_filtration.py",
    "VariantRecalibrator": "benchmark_variant_recalibrator.py",
    "VariantsToTable": "benchmark_variants_to_table.py",
}

# --------------------------------------------------------------------------
# Test-name prefix -> tool, for the accuracy half.  CTest names look like
# "fastgatk-<prefix>-<...>".  These prefixes are the kebab aliases used in
# CMakeLists.txt; the map is explicit so a new tool cannot silently attach
# itself to the wrong row.
# --------------------------------------------------------------------------
PREFIX_TOOL: dict[str, str] = {
    "hc": "HaplotypeCaller",
    "mutect2": "Mutect2",
    "bqsr": "BaseRecalibrator",
    "apply-bqsr": "ApplyBQSR",
    "apply-vqsr": "ApplyVQSR",
    "genotype-gvcf": "GenotypeGVCFs",
    "combine-gvcfs": "CombineGVCFs",
    "select-variants": "SelectVariants",
    "variant-filtration": "VariantFiltration",
    "variant-recalibrator": "VariantRecalibrator",
    "variant-eval": "VariantEval",
    "variants-to-table": "VariantsToTable",
    "mark-duplicates": "MarkDuplicates",
    "sort-sam": "SortSam",
    "reblock-gvcf": "ReblockGVCF",
    "reblock": "ReblockGVCF",
    "depth-of-coverage": "DepthOfCoverage",
    "filter-intervals": "FilterIntervals",
    "filter-mutect-calls": "FilterMutectCalls",
    "annotate-intervals": "AnnotateIntervals",
    "count-bases-in-reference": "CountBasesInReference",
    "count-reads": "CountReads",
    "compare-references": "CompareReferences",
    "check-reference-compatibility": "CheckReferenceCompatibility",
    "fasta-reference-maker": "FastaReferenceMaker",
    "fasta-alternate": "FastaAlternateReferenceMaker",
    "shift-fasta": "ShiftFasta",
    "index-feature-file": "IndexFeatureFile",
    "flag-stat": "FlagStat",
    "split-intervals": "SplitIntervals",
    "preprocess-intervals": "PreprocessIntervals",
    "analyze-covariates": "AnalyzeCovariates",
    "calculate-contamination": "CalculateContamination",
    "call-copy-ratio-segments": "CallCopyRatioSegments",
    "create-read-count-panel-of-normals": "CreateReadCountPanelOfNormals",
    "denoise-read-counts": "DenoiseReadCounts",
    "gather-bqsr": "GatherBQSRReports",
    "gather-tranches": "GatherTranches",
    "gather-vcfs": "GatherVcfs",
    "get-pileup": "GetPileupSummaries",
    "gather-pileup": "GatherPileupSummaries",
    "genomicsdb": "GenomicsDBImport",
    "learn-read-orientation-model": "LearnReadOrientationModel",
    "left-align": "LeftAlignAndTrimVariants",
    "model-segments": "ModelSegments",
    "validate-variants": "ValidateVariants",
    "collect-allelic-counts": "CollectAllelicCounts",
    "collect-f1r2-counts": "CollectF1R2Counts",
    "collect-read-counts": "CollectReadCounts",
}


def latest_regression_dir() -> Path | None:
    base = ROOT / ".diag" / "regression"
    if not base.is_dir():
        return None
    dirs = [p for p in base.iterdir() if p.is_dir()]
    return max(dirs, key=lambda p: p.stat().st_mtime) if dirs else None


def tool_for_test(name: str) -> str | None:
    # name like "fastgatk-hc-alleles-gatk-oracle"
    if not name.startswith("fastgatk-"):
        return None
    rest = name[len("fastgatk-"):]
    # longest prefix match
    for prefix in sorted(PREFIX_TOOL, key=len, reverse=True):
        if rest == prefix or rest.startswith(prefix + "-"):
            return PREFIX_TOOL[prefix]
    return None


def parse_accuracy(regression_dir: Path) -> dict[str, dict]:
    """Read omp.log + serial.log and return per-tool {pass, fail, total}."""
    acc: dict[str, dict] = {}
    pattern = re.compile(
        r"Test\s+#\d+:\s+(\S+)\s+\.+\s+(Passed|Failed|Skipped|Timeout)")
    for backend in ("omp", "serial"):
        log = regression_dir / f"{backend}.log"
        if not log.is_file():
            continue
        for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
            m = pattern.search(line)
            if not m:
                continue
            test_name, status = m.group(1), m.group(2)
            tool = tool_for_test(test_name)
            if tool is None:
                continue
            slot = acc.setdefault(tool, {"pass": 0, "fail": 0, "skip": 0, "total": 0})
            if status == "Passed":
                slot["pass"] += 1
            elif status == "Failed":
                slot["fail"] += 1
            else:
                slot["skip"] += 1
            slot["total"] += 1
    # sort descending by total then name
    return dict(sorted(acc.items(), key=lambda kv: (-kv[1]["total"], kv[0])))


def run_benchmark(script: str, timeout: int) -> dict:
    """Run one benchmark_*.py and parse its JSON (last stdout line)."""
    path = SCRIPTS / script
    env = os.environ.copy()
    env.setdefault("FASTGATK_NATIVE_BUILD", str(BUILD))
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    # keep benchmark fast and deterministic; java comparison is opt-out where
    # the benchmark supports it (see FASTGATK_SKIP_GATK_BENCHMARK).
    try:
        completed = subprocess.run(
            [sys.executable, str(path)], text=True, capture_output=True,
            timeout=timeout, env=env, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "script": script}
    except Exception as exc:  # noqa: BLE001
        return {"status": "error", "script": script, "error": str(exc)}
    if completed.returncode != 0:
        return {"status": "fail", "script": script,
                "stderr": completed.stderr[-400:]}
    # Benchmarks print JSON to stdout; some emit compact single-line JSON,
    # others pretty-print (indent=2).  Try the whole stdout first, then the
    # last line as a fallback.
    stdout = completed.stdout.strip()
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        pass
    for line in reversed(completed.stdout.splitlines()):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            return json.loads(line)
        except json.JSONDecodeError:
            continue
    return {"status": "no-json", "script": script}


def run_resource(end_to_end_timeout: int) -> dict:
    """benchmark_end_to_end.py is the only script that reports max RSS + FS IO."""
    return run_benchmark("benchmark_end_to_end.py", end_to_end_timeout)


def collect_efficiency(timeout: int) -> dict[str, dict]:
    eff: dict[str, dict] = {}
    for tool, script in sorted(TOOL_BENCHMARK.items()):
        result = run_benchmark(script, timeout)
        if "tool" not in result or result.get("status") not in (None, "pass"):
            result["_aggregate_status"] = result.get("status", "pass")
        eff[tool] = result
    return eff


def _sec(result: dict, *keys: str):
    for key in keys:
        if isinstance(result.get(key), (int, float)):
            return float(result[key])
    # some benchmarks nest the native timing under "native" (e.g. HC broad)
    nested = result.get("native")
    if isinstance(nested, dict):
        for key in keys:
            if isinstance(nested.get(key), (int, float)):
                return float(nested[key])
    return None


def render_markdown(accuracy: dict, efficiency: dict, resource: dict) -> str:
    lines = ["# fast-gatk 总对比报表", "",
             f"生成方式：`aggregate_report.py`（准确率解析自最新回归日志；",
             f"效率/资源来自 `benchmark_*.py`）", "",
             "| 工具 | 准确率(通过/总数) | native wall(s) | java wall(s) | 资源占用 |",
             "| --- | ---: | ---: | ---: | --- |"]
    tools = sorted(set(accuracy) | set(efficiency))
    for tool in tools:
        acc = accuracy.get(tool, {})
        acc_txt = f"{acc.get('pass', 0)}/{acc.get('total', 0)}" if acc else "—"
        if acc.get("fail"):
            acc_txt += f" (FAIL {acc['fail']})"
        eff = efficiency.get(tool, {})
        native_wall = _sec(eff, "p50_seconds", "wall_seconds", "warmup_seconds")
        java_wall = None
        if isinstance(eff.get("gatk_java"), dict):
            java_wall = _sec(eff["gatk_java"], "p50_seconds", "wall_seconds", "seconds")
        elif isinstance(eff.get("gatk_java"), list) and eff["gatk_java"]:
            java_wall = _sec(eff["gatk_java"][0], "seconds", "wall_seconds")
        eff_status = eff.get("_aggregate_status")
        if native_wall is None and eff_status and eff_status != "pass":
            native_txt = eff_status
        elif native_wall is None:
            native_txt = "—"
        else:
            native_txt = f"{native_wall:.4f}"
        java_txt = f"{java_wall:.4f}" if java_wall is not None else "—"
        # resource: end-to-end reports per-case max RSS; attach the case whose
        # name/summary tool matches this tool.
        res_txt = "—"
        cases = resource.get("cases") if isinstance(resource, dict) else None
        if isinstance(cases, list):
            for case in cases:
                if not isinstance(case, dict):
                    continue
                case_tool = case.get("name") or case.get("summary", {}).get("tool")
                if case_tool != tool:
                    continue
                rss = case.get("resource", {}).get("max_rss_kb")
                if rss:
                    res_txt = f"{rss} KB RSS"
                break
        lines.append(f"| {tool} | {acc_txt} | {native_txt} | {java_txt} | {res_txt} |")
    lines += ["",
              "> 注：",
              "> 1. 准确率 = 该工具所有 oracle/契约测试在双后端日志中的通过数/总数（失败标红）。",
              "> 2. native wall = benchmark 脚本报告的 p50/wall 墙钟；java wall 仅在 benchmark 同时跑",
              ">    Java 时才有值（多数 benchmark 只测 native，无 Java 对比）。",
              "> 3. 资源占用目前只有 `benchmark_end_to_end.py`（HC/Mutect2）报告 max RSS + FS 读写；",
              ">    其余工具的逐工具 RSS 未采集（benchmark 只测墙钟）。",
              "> 4. 若要「逐工具 native vs java 的效率 + RSS 对比」，需扩展各 `benchmark_*.py` 统一",
              ">    加 Java 计时与 `/usr/bin/time` 包裹。"]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regression-dir", default=None,
                        help="回归日志目录；默认取 .diag/regression/ 最新")
    parser.add_argument("--skip-benchmarks", action="store_true",
                        help="只做准确率，不跑 benchmark")
    parser.add_argument("--timeout", type=int, default=600,
                        help="每个 benchmark 的超时秒数")
    parser.add_argument("--outdir", default=str(ROOT / "report"))
    args = parser.parse_args()

    regression_dir = Path(args.regression_dir) if args.regression_dir else latest_regression_dir()
    if regression_dir is None:
        print("未找到回归日志目录；先跑 run_regression.sh 或用 --regression-dir 指定", file=sys.stderr)
        return 2
    accuracy = parse_accuracy(regression_dir)
    print(f"[accuracy] 从 {regression_dir} 解析到 {len(accuracy)} 个工具的测试结果")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.skip_benchmarks:
        efficiency: dict = {}
        resource: dict = {}
    else:
        print(f"[efficiency] 跑 {len(TOOL_BENCHMARK)} 个 benchmark（timeout={args.timeout}s）...")
        efficiency = collect_efficiency(args.timeout)
        print("[resource] 跑 benchmark_end_to_end.py ...")
        resource = run_resource(args.timeout)

    (outdir / "accuracy.json").write_text(
        json.dumps({"regression_dir": str(regression_dir), "tools": accuracy},
                   indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "efficiency.json").write_text(
        json.dumps(efficiency, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "resources.json").write_text(
        json.dumps(resource, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "REPORT.md").write_text(render_markdown(accuracy, efficiency, resource),
                                      encoding="utf-8")
    print(f"[done] 报表写入 {outdir}/REPORT.md（含 accuracy.json / efficiency.json / resources.json）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
