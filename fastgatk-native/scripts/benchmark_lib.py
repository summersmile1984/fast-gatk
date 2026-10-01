#!/usr/bin/env python3
"""Shared benchmark helpers: timed execution, disk logs, resource observations and a
matched Java/GATK baseline, so every ``benchmark_*.py`` reports the same
native-vs-Java comparison shape that ``aggregate_report.py`` consumes.

Conventions
-----------
* :func:`run_timed` wraps the workload in GNU time and records wall, max RSS,
  filesystem block counters, user/system CPU seconds, timestamps and log paths.
  Small callers retain captured stdout/stderr; large runs set ``capture_output=False``.
  An optional timeout cancels the entire process group, including surviving children.
* :func:`run_java_timed` runs the pinned GATK 4.6.2.0 jar with the same wrapper,
  so the Java side reports the same measurement fields.
* :func:`java_available` gates the optional Java baseline; benchmarks that add
  ``--include-java`` run the Java side only when the flag is set AND the pinned
  JDK/jar exist, keeping CTest runs (which never pass ``--include-java``) fast.
"""

from __future__ import annotations

import os
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def java_binary() -> Path:
    return Path(os.environ.get("JAVA", ROOT / "third_party/jdk17/bin/java"))


def gatk_jar() -> Path:
    return ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def java_available() -> bool:
    return java_binary().is_file() and gatk_jar().is_file()


def bgzip_and_index(plain: Path) -> Path:
    """bgzip + tabix a plain VCF so GATK can read it (requires an index).

    Uses ``-c`` (stdout) so the plain input is left intact, unlike ``bgzip``'s
    default which removes the source file.
    """
    bgzip = ROOT / "third_party/htslib-build/htslib-src/bgzip"
    tabix = ROOT / "third_party/htslib-build/htslib-src/tabix"
    gz = plain.with_suffix(plain.suffix + ".gz")
    with open(gz, "wb") as out:
        subprocess.run([str(bgzip), "-c", str(plain)], stdout=out, check=True)
    subprocess.run([str(tabix), "-f", "-p", "vcf", str(gz)], check=True,
                   capture_output=True)
    return gz


def write_simple_reference(reference: Path, contigs: list[tuple[str, int]]) -> Path:
    """Write a single-line-per-contig FASTA + ``.fai`` + ``.dict`` for GATK.

    Native reference-free benchmarks don't need this; it exists so a matched
    Java baseline can read the same synthetic reference without invoking the
    (slow, JVM) CreateSequenceDictionary tool.
    """
    text = "".join(f">{name}\n{'A' * length}\n" for name, length in contigs)
    reference.write_text(text, encoding="ascii")
    fai_lines = []
    offset = 0
    for name, length in contigs:
        header = len(f">{name}\n")
        fai_lines.append(f"{name}\t{length}\t{offset + header}\t{length}\t{length + 1}")
        offset += header + length + 1
    reference.with_name(reference.name + ".fai").write_text(
        "\n".join(fai_lines) + "\n", encoding="ascii")
    dict_lines = ["@HD\tVN:1.6"] + [
        f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs]
    reference.with_suffix(".dict").write_text(
        "\n".join(dict_lines) + "\n", encoding="ascii")
    return reference


def write_minimal_dict(reference: Path, contig: str, length: int) -> Path:
    """Write the minimal SAM-header ``.dict`` that GATK tools accept.

    GATK refuses to read a reference without a ``.dict``; native tools do not
    need one.  Several benchmarks therefore build only ``.fa`` + ``.fai`` for
    the native side and add this tiny ``.dict`` so the optional Java baseline
    can read the same reference.
    """
    dict_path = reference.with_suffix(".dict")
    dict_path.write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:{contig}\tLN:{length}\n", encoding="ascii")
    return dict_path


def run_timed(command: list[str], work: Path, name: str,
              check: bool = True, env: dict | None = None, *,
              capture_output: bool = True, timeout: float | None = None,
              cwd: Path | None = None, time_command: list[str] | None = None,
              terminate_grace_seconds: float = 5.0) -> dict:
    """Time the workload itself; stream logs to disk even for captured callers.

    Existing small benchmarks retain ``stdout``/``stderr`` strings. Large runs
    set ``capture_output=False`` and consume the log paths without loading them.
    FS counters are blocks, not bytes. ``time_command`` supports the pinned
    GNU time binary/loader inside the same resource envelope as the workload.
    """
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    timing_file = work / f"{name}.time"
    stdout_file, stderr_file = work / f"{name}.stdout.log", work / f"{name}.stderr.log"
    timing_file.unlink(missing_ok=True)
    actual = [*(time_command or ["/usr/bin/time"]), "-f",
              "FASTGATK_TIME\t%e\t%M\t%I\t%O\t%U\t%S", "-o", str(timing_file), *command]
    started_at = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    timed_out = False
    with stdout_file.open("wb") as stdout, stderr_file.open("wb") as stderr:
        process = subprocess.Popen(actual, stdout=stdout, stderr=stderr,
                                   env=env if env is not None else os.environ.copy(),
                                   cwd=cwd, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as failure:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=terminate_grace_seconds)
            except subprocess.TimeoutExpired:
                pass
            # Descendants may survive after their parent exits on SIGTERM.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            returncode = process.wait()
            if isinstance(failure, KeyboardInterrupt):
                raise
            timed_out = True
    wall = time.perf_counter() - start
    result = {
        "seconds": wall, "max_rss_kb": 0, "filesystem_inputs": 0,
        "filesystem_outputs": 0, "returncode": returncode,
        "stdout": stdout_file.read_text(errors="replace") if capture_output else "",
        "stderr": stderr_file.read_text(errors="replace") if capture_output else "",
        "stdout_path": str(stdout_file), "stderr_path": str(stderr_file),
        "timing_path": str(timing_file), "timing_measured": False,
        "user_seconds": None, "system_seconds": None, "cpu_seconds": None,
        "timeout": timed_out, "signal": -returncode if returncode < 0 else None,
        "started_at": started_at, "ended_at": datetime.now(timezone.utc).isoformat(),
        "argv": command, "timing_argv": actual,
    }
    if timing_file.is_file():
        for line in timing_file.read_text(errors="replace").splitlines():
            fields = line.split("\t")
            if len(fields) != 7 or fields[0] != "FASTGATK_TIME":
                continue
            try:
                result.update(time_wall_seconds=float(fields[1]), max_rss_kb=int(fields[2]),
                              filesystem_inputs=int(fields[3]), filesystem_outputs=int(fields[4]),
                              user_seconds=float(fields[5]), system_seconds=float(fields[6]))
            except ValueError:
                continue
            result["cpu_seconds"] = result["user_seconds"] + result["system_seconds"]
            result["timing_measured"] = True
    if check and (returncode != 0 or timed_out):
        raise RuntimeError(f"{name} exited {returncode}, timeout={timed_out}\n"
                           f"stdout: {stdout_file}\nstderr: {stderr_file}")
    return result


def run_java_timed(tool: str, args: list[str], work: Path, name: str,
                   check: bool = True, env: dict | None = None) -> dict:
    """Run the pinned GATK jar's *tool* under the same timed wrapper."""
    return run_timed(
        [str(java_binary()), "-Xmx1g", "-jar", str(gatk_jar()), tool, *args],
        work, name, check=check, env=env)


def summarize(samples: list[dict]) -> dict:
    """p50/warmup/p95 wall over :func:`run_timed` samples + peak max RSS."""
    if not samples:
        return {"repetitions": 0}
    timings = sorted(sample["seconds"] for sample in samples)
    return {
        "repetitions": len(samples),
        "warmup_seconds": timings[0],
        "p50_seconds": timings[len(timings) // 2],
        "p95_seconds": timings[-1],
        "max_rss_kb": max(sample["max_rss_kb"] for sample in samples),
        "filesystem_inputs": max(sample["filesystem_inputs"]
                                 for sample in samples),
        "filesystem_outputs": max(sample["filesystem_outputs"]
                                  for sample in samples),
    }
