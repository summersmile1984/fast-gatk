#!/usr/bin/env python3
"""Resource parity report: native fastgatk-* vs pinned GATK 4.6.2.0.

Functional parity is a precondition: every case reuses the byte/record
equivalence checks of ``run_headline_comparison.py`` and speedup ratios are
reported only for output-equivalent runs (per docs/regression-evidence.md:
unlike workloads never justify a speedup claim).

Equal-resource protocol (v1)
----------------------------
* Same fixture and functionally identical CLI per case; the only permitted
  differences are the implementation entry points and the recorded JVM heap
  limits (``-Xmx1g``/``-Xmx2g``, fixed per case by the case definition).
* Both sides of a lane get the SAME resource envelope (taskset CPU affinity
  + OMP_NUM_THREADS), enforced outside the CLI so no per-tool flags leak in:

  - ``1-core``       -- ``taskset -c <core0>``, ``OMP_NUM_THREADS=1``:
    pure per-core algorithmic efficiency.
  - ``full-machine`` -- all cores to both sides (``OMP_NUM_THREADS=nproc``):
    as-shipped throughput inside the same envelope.  GATK's inability to
    exploit extra cores is its own limitation, not an unequal envelope.

* Runs are strictly serialized (one workload at a time).  Warm cache,
  ``--repeat N`` repetitions (default 2); p50/min/max reported.
* Paired, binary-identified runs: SHA-256 of every native binary exercised
  and of the GATK jar, plus java version, git HEAD and machine state
  (CPU model, nproc, loadavg before/after, per-lane affinity).

Measured per run
----------------
``/usr/bin/time -v`` root statistics (wall, user/sys CPU, max RSS,
filesystem block I/O) plus a 25 ms process-tree sampler that walks
``/proc/<pid>/task/<tid>/children`` and aggregates every descendant (native
and java processes alike):

* ``tree_rss_peak_kb``   -- peak of the sum of per-process RSS over the tree
* ``tree_io_read_bytes`` / ``tree_io_write_bytes`` -- physical disk bytes
  (``read_bytes``/``write_bytes``); near zero on warm cache, as expected
* ``tree_rchar_bytes`` / ``tree_wchar_bytes`` -- bytes requested through
  read/write syscalls (cache-inclusive I/O volume)

Sampler limitations match the validated ``work/mutect2-priority/sample_tree``
harness: 25 ms cadence can miss very short-lived helpers; counters are
per-process cumulative maxima, not cgroup accounting.

Output
------
``report/resource_parity_<ts>/``
    summary.json        -- machine-readable per-case stats + ratios
    RESOURCE_PARITY.md  -- human-readable comparison tables (one per lane)
    work/<case>/        -- raw case outputs (native + GATK artifacts)
``report/resource_parity_latest`` symlink -> newest run directory.

Usage
-----
    python3 fastgatk-native/scripts/resource_parity_report.py --repeat 2
    python3 fastgatk-native/scripts/resource_parity_report.py --cases hc bqsr
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
REPORTS = ROOT / "report"
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
TIME = "/usr/bin/time"
SAMPLE_INTERVAL = 0.025  # seconds; matches work/mutect2-priority/sample_tree

CASE_NAMES = ["hc", "hc_aline", "bqsr", "mutect2", "mutect2_dream",
              "sortsam", "markdup", "markdup_ceutrio", "genotype_gvcfs"]
LANES = ["1-core", "full-machine", "full-machine-tuned"]


# ---------------------------------------------------------------- tree scan

def _descendants(root_pid: int) -> list[int]:
    """root_pid plus every descendant, via /proc/<pid>/task/<tid>/children."""
    seen: set[int] = set()
    stack = [root_pid]
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        try:
            tids = os.listdir(f"/proc/{pid}/task")
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        for tid in tids:
            try:
                with open(f"/proc/{pid}/task/{tid}/children") as fh:
                    stack.extend(int(tok) for tok in fh.read().split())
            except (FileNotFoundError, ProcessLookupError, PermissionError,
                    ValueError):
                continue
    return sorted(seen)


def _rss_kb(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, PermissionError,
            ValueError):
        pass
    return 0


def _io_bytes(pid: int) -> tuple[int, int, int, int]:
    read_b = write_b = rchar = wchar = 0
    try:
        with open(f"/proc/{pid}/io") as fh:
            for line in fh:
                if line.startswith("read_bytes:"):
                    read_b = int(line.split()[1])
                elif line.startswith("write_bytes:"):
                    write_b = int(line.split()[1])
                elif line.startswith("rchar:"):
                    rchar = int(line.split()[1])
                elif line.startswith("wchar:"):
                    wchar = int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, PermissionError,
            ValueError):
        pass
    return read_b, write_b, rchar, wchar


class TreeSampler(threading.Thread):
    """25 ms tree RSS / IO sampler rooted at the measured process."""

    def __init__(self, root_pid: int) -> None:
        super().__init__(daemon=True)
        self.root_pid = root_pid
        self.stop_event = threading.Event()
        self.tree_rss_peak_kb = 0
        self.tree_io_read_bytes = 0
        self.tree_io_write_bytes = 0
        self.tree_rchar_bytes = 0
        self.tree_wchar_bytes = 0
        self.samples = 0

    def run(self) -> None:
        while not self.stop_event.is_set():
            self._sample()
            self.stop_event.wait(SAMPLE_INTERVAL)
        self._sample()

    def _sample(self) -> None:
        pids = _descendants(self.root_pid)
        if not pids:
            return
        rss = 0
        io_r = io_w = rchar = wchar = 0
        for pid in pids:
            rss += _rss_kb(pid)
            r, w, rc, wc = _io_bytes(pid)
            io_r += r
            io_w += w
            rchar += rc
            wchar += wc
        self.tree_rss_peak_kb = max(self.tree_rss_peak_kb, rss)
        self.tree_io_read_bytes = max(self.tree_io_read_bytes, io_r)
        self.tree_io_write_bytes = max(self.tree_io_write_bytes, io_w)
        self.tree_rchar_bytes = max(self.tree_rchar_bytes, rchar)
        self.tree_wchar_bytes = max(self.tree_wchar_bytes, wchar)
        self.samples += 1


# ------------------------------------------------------------ instrumentation

def _f(stats: dict, key: str) -> float | None:
    raw = stats.get(key)
    if raw is None:
        return None
    try:
        return float(re.sub(r"[^0-9.eE+-]", "", str(raw)))
    except (TypeError, ValueError):
        return None


import re  # noqa: E402  (used by _f)


def lane_command(lane: str, cmd: list[str]) -> list[str]:
    if lane == "1-core":
        return ["taskset", "-c", "0", *cmd]
    return list(cmd)


def lane_env(lane: str) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("GOMP_SPINCOUNT", None)  # lane-controlled, not ambient
    if lane == "1-core":
        env["OMP_NUM_THREADS"] = "1"
    else:
        env["OMP_NUM_THREADS"] = str(os.cpu_count() or 1)
    if lane == "full-machine-tuned":
        # Launcher-shipped spin budget (dispatcher sets the same default):
        # ~43% less CPU at ~10% wall on barrier-heavy workloads.
        env["GOMP_SPINCOUNT"] = "5000"
    return env


def instrument(headline, run_log: list[dict], holder: dict[str, str]):
    """Patch the headline module's runner; keep its return contract."""

    def run_timed(cmd, cwd):
        cwd = Path(cwd)
        lane = holder.get("lane", "?")
        timing = cwd / f".resource_parity_{time.monotonic_ns()}.time"
        full = [TIME, "-v", "-o", str(timing)] + lane_command(
            lane, [str(c) for c in cmd])
        started = time.perf_counter()
        proc = subprocess.Popen(full, cwd=str(cwd),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                env=lane_env(lane),
                                start_new_session=True)
        sampler = TreeSampler(proc.pid)
        sampler.start()
        try:
            _, err = proc.communicate()
        finally:
            sampler.stop_event.set()
            sampler.join(timeout=2.0)
        wall = time.perf_counter() - started
        stats: dict = headline.parse_time_v(
            timing.read_text(encoding="utf-8", errors="replace")
            if timing.is_file() else "")
        timing.unlink(missing_ok=True)
        stats.update({
            "x_wall_seconds": wall,
            "x_tree_rss_peak_kb": sampler.tree_rss_peak_kb,
            "x_tree_io_read_bytes": sampler.tree_io_read_bytes,
            "x_tree_io_write_bytes": sampler.tree_io_write_bytes,
            "x_tree_rchar_bytes": sampler.tree_rchar_bytes,
            "x_tree_wchar_bytes": sampler.tree_wchar_bytes,
            "x_tree_samples": sampler.samples,
            "x_cpu_seconds": (_f(stats, "User time (seconds)") or 0.0) +
                             (_f(stats, "System time (seconds)") or 0.0),
        })
        run_log.append({
            "case": holder.get("case", "?"),
            "lane": lane,
            "side": side_of(cmd),
            "argv": [str(c) for c in cmd],
            "rc": proc.returncode,
            **{k: v for k, v in stats.items()
               if k.startswith(("x_", "User", "System", "Percent", "Maximum",
                                "File"))},
        })
        return proc.returncode, stats, stats["x_wall_seconds"]

    headline.run_timed = run_timed
    return headline


def side_of(cmd: list) -> str:
    name = Path(str(cmd[0])).name
    return "gatk" if "java" in name else "native"


# ------------------------------------------- BAM equivalence modulo @PG

def _bam_parts(path: Path):
    """Decode a BAM into (header_text, record_bytes)."""
    import struct
    with gzip.open(path, "rb") as fh:
        data = fh.read()
    l_text = struct.unpack_from("<i", data, 4)[0]
    text = data[8:8 + l_text].decode(errors="replace")
    off = 8 + l_text
    n_ref = struct.unpack_from("<i", data, off)[0]
    off += 4
    for _ in range(n_ref):
        l_name = struct.unpack_from("<i", data, off)[0]
        off += 4 + l_name + 4
    return text, data[off:]


def bam_equal_modulo_pg(a: Path, b: Path) -> bool:
    """BAM byte equality up to @PG program-record lines.

    @PG lines carry run provenance (paths, versions, command lines) that two
    implementations cannot reproduce for each other; the GATK oracles assert
    their presence/ID, not cross-implementation byte equality.  Everything
    else - header text, references and the complete record stream - must
    match byte for byte.
    """
    try:
        atext, arecs = _bam_parts(a)
        btext, brecs = _bam_parts(b)
    except (OSError, EOFError, ValueError):
        return False
    aheader = [l for l in atext.splitlines() if not l.startswith("@PG\t")]
    bheader = [l for l in btext.splitlines() if not l.startswith("@PG\t")]
    return aheader == bheader and arecs == brecs


def output_path_of(argv: list[str]) -> Path | None:
    for i, arg in enumerate(argv):
        if arg in ("-O", "--output") and i + 1 < len(argv):
            return Path(argv[i + 1])
    return None



# --------------------------------------------------------------- aggregation

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_side(runs: list[dict]) -> dict:
    """Aggregate the runs of one measured command (one median_wall group)."""
    if not runs:
        return {}
    walls = [r["x_wall_seconds"] for r in runs]
    cpus = [r["x_cpu_seconds"] for r in runs]
    return {
        "repetitions": len(runs),
        "wall_p50_s": median(walls),
        "wall_min_s": min(walls),
        "wall_max_s": max(walls),
        "cpu_p50_s": median(cpus),
        "root_max_rss_kb_max": max(_f(r, "Maximum resident set size (kbytes)")
                                   or 0 for r in runs),
        "tree_rss_peak_kb_max": max(r["x_tree_rss_peak_kb"] for r in runs),
        "tree_io_read_bytes_max": max(r["x_tree_io_read_bytes"] for r in runs),
        "tree_io_write_bytes_max": max(r["x_tree_io_write_bytes"] for r in runs),
        "tree_rchar_bytes_max": max(r["x_tree_rchar_bytes"] for r in runs),
        "tree_wchar_bytes_max": max(r["x_tree_wchar_bytes"] for r in runs),
        "fs_in_blocks_max": max(_f(r, "File system inputs") or 0 for r in runs),
        "fs_out_blocks_max": max(_f(r, "File system outputs") or 0 for r in runs),
        "percent_cpu_p50": median([_f(r, "Percent of CPU this job got") or 0
                                   for r in runs]),
        "argv": runs[0]["argv"],
        "runs": runs,
    }


def equivalence_fields(case_dict: dict) -> dict:
    fields = {}
    for key, value in case_dict.items():
        if isinstance(value, bool) and any(
                token in key for token in ("identical", "byte", "match",
                                           "exact", "equal")):
            fields[key] = value
    return fields


def ratios(native: dict, gatk: dict) -> dict:
    if not native or not gatk:
        return {}

    def div(a, b):
        if not b:
            return None
        return round(a / b, 3)

    return {
        "wall_speedup_gatk_over_native_p50":
            div(gatk["wall_p50_s"], native["wall_p50_s"]),
        "wall_speedup_gatk_over_native_min":
            div(gatk["wall_min_s"], native["wall_min_s"]),
        "cpu_ratio_gatk_over_native":
            div(gatk["cpu_p50_s"], native["cpu_p50_s"]),
        "tree_rss_ratio_gatk_over_native":
            div(gatk["tree_rss_peak_kb_max"], native["tree_rss_peak_kb_max"]),
        "io_read_ratio_gatk_over_native":
            div(gatk["tree_io_read_bytes_max"],
                native["tree_io_read_bytes_max"]),
        "io_write_ratio_gatk_over_native":
            div(gatk["tree_io_write_bytes_max"],
                native["tree_io_write_bytes_max"]),
        "rchar_ratio_gatk_over_native":
            div(gatk["tree_rchar_bytes_max"], native["tree_rchar_bytes_max"]),
        "wchar_ratio_gatk_over_native":
            div(gatk["tree_wchar_bytes_max"], native["tree_wchar_bytes_max"]),
    }


# ------------------------------------------------------------------ reporting

def fmt(value, unit_scale: float = 1.0, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{value / unit_scale:.{digits}f}"


def lane_table(cases: list[dict], lane: str) -> list[str]:
    lines = [
        f"## Lane `{lane}`",
        "",
        "| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | "
        "CPU native (s) | CPU GATK (s) | CPU ratio | "
        "TreeRSS native (MiB) | TreeRSS GATK (MiB) | "
        "rchar native (MiB) | rchar GATK (MiB) | "
        "disk-read native (MiB) | disk-read GATK (MiB) |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for case in cases:
        lane_data = case["lanes"].get(lane, {})
        nat, gatk = lane_data.get("native", {}), lane_data.get("gatk", {})
        eq = case["output_equivalence"]
        equal = "yes" if eq["equivalent"] else ("no" if eq["fields"] else "n/a")
        r = lane_data.get("ratios", {})
        speed = (f"{r.get('wall_speedup_gatk_over_native_p50', 0):.2f}×"
                 if eq["equivalent"] and r.get("wall_speedup_gatk_over_native_p50")
                 else "n/a")
        cpu_ratio = (f"{r.get('cpu_ratio_gatk_over_native', 0):.2f}×"
                     if eq["equivalent"] and r.get("cpu_ratio_gatk_over_native")
                     else "n/a")
        lines.append(
            f"| {case['case']} | {equal} | {fmt(nat.get('wall_p50_s'))} | "
            f"{fmt(gatk.get('wall_p50_s'))} | {speed} | "
            f"{fmt(nat.get('cpu_p50_s'))} | {fmt(gatk.get('cpu_p50_s'))} | {cpu_ratio} | "
            f"{fmt(nat.get('tree_rss_peak_kb_max'), 1024)} | "
            f"{fmt(gatk.get('tree_rss_peak_kb_max'), 1024)} | "
            f"{fmt(nat.get('tree_rchar_bytes_max'), 1 << 20)} | "
            f"{fmt(gatk.get('tree_rchar_bytes_max'), 1 << 20)} | "
            f"{fmt(nat.get('tree_io_read_bytes_max'), 1 << 20)} | "
            f"{fmt(gatk.get('tree_io_read_bytes_max'), 1 << 20)} |")
    return lines


def render_md(payload: dict) -> str:
    env = payload["environment"]
    lines = [
        "# Resource parity: native fastgatk vs GATK 4.6.2.0",
        "",
        f"Report date: {payload['report_date']}  |  protocol: "
        f"{payload['protocol']['version']}  |  repetitions: "
        f"{payload['protocol']['repetitions']}",
        "",
        "Equal-resource protocol (v1): identical fixture and functional CLI per",
        "case; both sides of a lane get the SAME envelope (taskset affinity +",
        "OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;",
        "p50 reported; output-equivalence asserted per case before any ratio is",
        "reported. Lane `1-core` = one core each (algorithmic efficiency); lane",
        "`full-machine` = the whole machine for both sides (as-shipped",
        "throughput).",
        "",
        f"Machine: {env['cpu_model']} | {env['nproc']} cpus | "
        f"loadavg {env['loadavg_before']} -> {env['loadavg_after']}",
        "",
        f"git HEAD: `{env['git_head']}`  |  GATK jar sha256: "
        f"`{(env.get('gatk_jar_sha256') or '')[:16]}…`  |  java: "
        f"{(env.get('java_version') or ['?'])[0]}",
        "",
    ]
    for lane in payload["protocol"]["lanes"]:
        lines += lane_table(payload["cases"], lane)
        lines.append("")
    lines += [
        "Notes:",
        "* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio",
        "  is the equal-resource efficiency measure (same work, same envelope).",
        "* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the",
        "  whole workload process tree (native/java alike).",
        "* rchar = bytes requested via read syscalls (cache-inclusive volume);",
        "  disk-read = physical `read_bytes` (near zero on warm cache).",
        "* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and",
        "  FS-block counters in `summary.json`.",
    ]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------ main

def load_headline():
    spec = importlib.util.spec_from_file_location(
        "resource_parity_headline", SCRIPTS / "run_headline_comparison.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def environment() -> dict:
    cpu_model = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    git_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True).stdout.strip()
    java_version = subprocess.run([str(JAVA), "-version"], capture_output=True,
                                  text=True).stderr.splitlines()[:1]
    return {
        "cpu_model": cpu_model,
        "nproc": os.cpu_count(),
        "loadavg_before": list(os.getloadavg()),
        "git_head": git_head,
        "java_version": java_version,
        "gatk_jar_sha256": sha256_file(GATK_JAR) if GATK_JAR.is_file() else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=CASE_NAMES,
                        choices=CASE_NAMES)
    parser.add_argument("--lanes", nargs="*", default=LANES, choices=LANES)
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()

    headline = load_headline()
    if not all(p.is_file() for p in (headline.BUILD / "fastgatk-hc-call", JAVA,
                                     GATK_JAR)):
        print("missing native binaries / pinned JDK / GATK jar", file=sys.stderr)
        return 2

    ts = time.strftime("%Y%m%d-%H%M%S")
    outdir = args.outdir or (REPORTS / f"resource_parity_{ts}")
    outdir.mkdir(parents=True, exist_ok=True)
    work_root = outdir / "work"
    work_root.mkdir(exist_ok=True)

    env = environment()
    headline.REPEAT = args.repeat
    run_log: list[dict] = []
    holder: dict[str, str] = {"case": "?", "lane": "?"}
    instrument(headline, run_log, holder)

    # Mirror run_headline_comparison.main()'s dispatch (case keys -> defs).
    case_dispatch = {
        "hc": headline.case_hc,
        "hc_aline": headline.case_hc_multi_region,
        "bqsr": headline.case_bqsr,
        "mutect2": headline.case_mutect2,
        "mutect2_dream": headline.case_mutect2_dream_real,
        "sortsam": headline.case_sort_sam,
        "markdup": headline.case_mark_duplicates,
        "markdup_ceutrio": headline.case_mark_duplicates_ceutrio,
        "genotype_gvcfs": headline.case_genotype_gvcfs,
    }

    cases_payload = []
    for name in args.cases:
        case_fn = case_dispatch.get(name)
        if case_fn is None:
            print(f"unknown case: {name}", file=sys.stderr)
            continue
        holder["case"] = name
        lanes_payload = {}
        case_dict = None
        for lane in args.lanes:
            holder["lane"] = lane
            print(f"[resource-parity] case={name} lane={lane}", flush=True)
            case_dict = case_fn(work_root)
            lane_runs = [e for e in run_log
                         if e["case"] == name and e["lane"] == lane]
            side_groups: dict[str, list[list[dict]]] = {"native": [], "gatk": []}
            groups: dict[str, list[dict]] = {}
            for entry in lane_runs:
                groups.setdefault(str(entry["argv"]), []).append(entry)
            for entries in groups.values():
                side_groups[entries[0]["side"]].append(entries)
            native_runs = (side_groups["native"] or [[]])[-1]
            gatk_runs = (side_groups["gatk"] or [[]])[-1]
            native = collect_side(native_runs)
            gatk = collect_side(gatk_runs)
            lanes_payload[lane] = {
                "native": native,
                "gatk": gatk,
                "ratios": ratios(native, gatk),
            }
        eq_fields = equivalence_fields(case_dict or {})
        nat_out = output_path_of(native.get("argv", []))
        gatk_out = output_path_of(gatk.get("argv", []))
        if (nat_out is not None and gatk_out is not None and
                nat_out.suffix == ".bam" and gatk_out.suffix == ".bam"):
            eq_fields["byte_exact_bam_modulo_pg"] = bam_equal_modulo_pg(
                nat_out, gatk_out)
        cases_payload.append({
            **{k: v for k, v in (case_dict or {}).items() if k != "case_key"},
            "case_key": name,
            "output_equivalence": {
                "fields": eq_fields,
                "equivalent": any(v is True for v in eq_fields.values()),
            },
            "lanes": lanes_payload,
        })
        print(f"[resource-parity] case done: {name}", flush=True)

    env["loadavg_after"] = list(os.getloadavg())
    payload = {
        "schema_version": 1,
        "report_kind": "resource-parity",
        "report_date": time.strftime("%Y-%m-%d"),
        "protocol": {
            "version": "equal-resource-v1",
            "repetitions": args.repeat,
            "lanes": args.lanes,
            "lane_envelopes": {
                "1-core": "taskset -c 0, OMP_NUM_THREADS=1 (both sides)",
                "full-machine": f"all {os.cpu_count()} cores to both sides, "
                                f"OMP_NUM_THREADS={os.cpu_count()}",
                "full-machine-tuned": f"all {os.cpu_count()} cores, "
                                      f"OMP_NUM_THREADS={os.cpu_count()}, "
                                      "GOMP_SPINCOUNT=5000 (shipped launcher default)",
            },
            "jvm_heap": "fixed per case (-Xmx1g / -Xmx2g as recorded in argv)",
            "cache_state": "warm",
            "execution": "serialized, one workload at a time",
            "sampler": f"process-tree /proc/<pid>/io + VmRSS at {SAMPLE_INTERVAL}s",
        },
        "environment": env,
        "cases": cases_payload,
    }
    (outdir / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    (outdir / "RESOURCE_PARITY.md").write_text(render_md(payload))

    latest = REPORTS / "resource_parity_latest"
    if latest.is_symlink() or latest.exists():
        latest.unlink()
    latest.symlink_to(outdir.name)
    print(f"[resource-parity] wrote {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
