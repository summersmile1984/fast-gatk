#!/usr/bin/env python3
"""Describe sampled native/Java resources in the 2026-09-23 rerun sidecars.

This is not a paired benchmark: scripts use different invocation counts,
workloads and parameters. Historical resource data are unvalidated and include
manual estimates, cumulative CPU accounting and incomplete sampler finalization.
Re-aggregation cannot repair those source records or establish a speedup.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def load_tool(rerun_root: Path, slug: str) -> tuple[dict, list[dict]]:
    """Return ``(header, scripts)`` from the requested repository."""
    fp = rerun_root / slug / f"{slug}-rerun-20260923.json"
    if not fp.is_file():
        return {}, []
    data = json.loads(fp.read_text(encoding="utf-8"))
    return data, data.get("scripts", [])


def aggregate_metrics(scripts: list[dict]) -> dict:
    """Max recorded RSS; sum recorded lifetimes/counts, with explicit coverage.

    Unsampled records and absent/null fields are unavailable, not zero.
    RSS zero is not a measured process footprint. Partial aggregates describe
    only available records, never the complete workload.
    """
    result = {
        "sampled_scripts": sum(bool(s.get("resource_sampled")) for s in scripts),
        "total_scripts": len(scripts),
        "passed": sum(bool(s.get("passed")) and not s.get("skipped", False) for s in scripts),
        "skipped": sum(bool(s.get("skipped")) for s in scripts),
        "failed": sum(not s.get("passed") and not s.get("skipped") for s in scripts),
    }
    for cls_name in ("native", "java", "other"):
        records = [s.get(cls_name) or {} for s in scripts if s.get("resource_sampled")]
        bucket = {"measurement_scripts": {}, "cmdlines": []}
        for key in ("peak_rss_bytes", "wall_clock_seconds_sum", "process_count"):
            values = [record[key] for record in records if record.get(key) is not None
                      and (key != "peak_rss_bytes" or record[key] > 0)]
            bucket["measurement_scripts"][key] = len(values)
            bucket[key] = (max(values) if key == "peak_rss_bytes" else sum(values)) if values else None
        for record in records:
            bucket["cmdlines"].extend(record.get("sample_cmdlines") or [])
        result[cls_name] = bucket
    return result


def format_metric(bucket: dict, key: str, total: int) -> str:
    value = bucket[key]
    if value is None:
        return "unavailable"
    if key == "peak_rss_bytes":
        text = f"{value / (1024 * 1024):.1f} MiB"
    elif key == "wall_clock_seconds_sum":
        text = f"{value:.2f}s"
    else:
        text = str(value)
    measured = bucket["measurement_scripts"][key]
    return text if measured == total else f"{text} ({measured}/{total} scripts)"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=str(REPO))
    args = parser.parse_args()
    rerun_root = Path(args.repo).resolve() / "fastgatk-native/evidence/2026-09-23-rerun"
    tool_dirs = sorted(p for p in rerun_root.iterdir()
                       if p.is_dir() and (p / f"{p.name}-rerun-20260923.json").is_file())
    today = dt.date.today().isoformat()
    rows = []
    for tool_dir in tool_dirs:
        header, scripts = load_tool(rerun_root, tool_dir.name)
        if not header:
            continue
        rows.append((tool_dir.name, header, aggregate_metrics(scripts)))

    lines = [
        f"# Native and Java sampled resources — {today}",
        "",
        "## Historical data: unvalidated, not a benchmark",
        "",
        "The 2026-09-23 resource records are not valid performance evidence.",
        "They include manually estimated Mutect2 values, cumulative child CPU",
        "accounting, incomplete sampler finalization, and unlike workloads and",
        "invocation counts. Re-aggregation does not authenticate or repair them.",
        "No native optimization, speedup, memory reduction or leak is established.",
        "Fresh isolated Mutect2 evidence is pending in",
        "[work/mutect2-priority/](../../../work/mutect2-priority/).",
        "",
        "## Descriptive sidecar values",
        "",
        f"Tools represented: **{len(rows)}**; `genomicsdb-export` is a native-only",
        "bridge with no verify script (a registry fact, not a sampler inference).",
        "",
        "Exit-status counts are historical script outcomes, not output-parity",
        "proof. Consult each oracle's actual assertions and retained output.",
        "",
        "| Tool | Reported pass/fail/skip | Sampled scripts | Native recorded pids | Java recorded pids | Native lifetime sum | Java lifetime sum | Native max sampled RSS | Java max sampled RSS |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for slug, header, agg in rows:
        metrics = [format_metric(agg[cls], key, agg["total_scripts"])
                   for key in ("process_count", "wall_clock_seconds_sum", "peak_rss_bytes")
                   for cls in ("native", "java")]
        lines.append(
            f"| `{header.get('tool', slug)}` (`{slug}`) | "
            f"{header.get('scripts_passed', agg['passed'])}/"
            f"{header.get('scripts_failed', agg['failed'])}/"
            f"{header.get('scripts_skipped', agg['skipped'])} | "
            f"{agg['sampled_scripts']}/{agg['total_scripts']} | "
            + " | ".join(metrics) + " |"
        )
    lines.extend([
        "",
        "## Interpretation limits",
        "",
        "* RSS is the maximum of available per-script class peaks, never their",
        "  sum. The historical sampler takes individual process RSS maxima,",
        "  not simultaneous process-tree memory or an OS-certified peak.",
        "* Lifetime sums and pid counts are recorded sampler values. Historical",
        "  counts/lifetimes were charged when a pid disappeared; pids still",
        "  tracked at sampler shutdown were not finalized. Short-lived processes",
        "  may also be missed. Zero counts do not imply native-only execution.",
        "* Missing/null fields and unsampled scripts are unavailable, not zeros.",
        "  Zero RSS is also unavailable. Partial fields show their contributing",
        "  script count; otherwise all scripts supplied a recorded value.",
        "  Coverage is not evidence that a recorded measurement is valid.",
        "* Process classes follow `comm`: `java`, `fastgatk-*`, or other helpers.",
        "  This does not isolate equivalent commands or prove a GATK invocation.",
        "* Unequal invocation counts, helper tools, workload parameters and",
        "  concurrent execution prevent fair native-versus-Java comparisons.",
        "  There are deliberately no speedup ratios or anomaly rankings.",
        "",
    ])
    md_path = rerun_root / "NATIVE_VS_JAVA.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Descriptive historical resource report: {len(rows)} tools; not a benchmark.")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
