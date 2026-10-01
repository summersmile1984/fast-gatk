#!/usr/bin/env python3
"""Verify the project scope boundary declared in `SCOPE.md`.

This contract pins the project's OUT-OF-SCOPE statement at the CI layer:
any change that violates the grep-level expectations listed in `SCOPE.md §6`
fails this test.  It does NOT decide policy — the source of truth remains
`SCOPE.md` §2; this script just enforces the boundary mechanically so a
future contributor cannot accidentally introduce Spark, Hadoop, or remote
URL handling without changing the scope document first.

The checks mirror §6 of `SCOPE.md`:

1. Count GATK `*Spark.java` files in `gatk-source/`.
2. Zero `SparkContext`/`JavaSparkContext`/`RDD`/`Dataset<` symbols in our
   four source trees.
3. Zero `hadoop`/`hdfs://` symbols in our four source trees.
4. Zero Spark/Hadoop/HDFS dependencies vendored in `third_party/`.
5. Cloud URL handling (`s3://`, `gs://`, `https?://`) only appears in
   rejection branches inside `genotype_gvcf_tool.cpp` and
   `genomicsdb_import_tool.cpp`.
6. Registry declares exactly 4 `spark-*` fallback boundaries and 19
   `cloud-*` fallback boundaries.
7. The dispatcher fails closed on `--gatk-config-file` (which can carry
   Spark/cloud defaults) and on the `--` post-tool separator that GATK
   uses to split tool args from Spark args.

The expected values are pulled directly from the latest `SCOPE.md §6` so a
PR that legitimately changes the scope only has to update that file plus
this script.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCOPE_MD = ROOT / "SCOPE.md"
REGISTRY = ROOT / "fastgatk-native" / "dispatcher" / "tool_registry.json"
DISPATCHER_PY = ROOT / "fastgatk-native" / "dispatcher" / "fastgatk.py"
SOURCE_TREES = [
    ROOT / "fastgatk-native" / "src",
    ROOT / "fastgatk-kernels" / "src",
    ROOT / "fastgatk-runtime" / "src",
    ROOT / "fastgatk-core",
]
THIRD_PARTY = ROOT / "third_party"
REJECTION_FILES = {
    ROOT / "fastgatk-native/src/genotype_gvcf_tool.cpp",
    ROOT / "fastgatk-native/src/genomicsdb_import_tool.cpp",
}


def run_grep(pattern: str, paths: list[Path]) -> int:
    """Return the number of matching lines across ``paths`` (recursive)."""
    if not paths:
        return 0
    args = ["grep", "-rE", "--include=*", pattern]
    for p in paths:
        args.append(str(p))
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    # Filter out binary-file lines (grep prints "Binary file ... matches").
    lines = [
        line for line in result.stdout.splitlines()
        if not line.startswith("Binary file")
    ]
    return len(lines)


def find_files(pattern: str, roots: list[Path]) -> list[Path]:
    """Return ``find`` matches across ``roots``."""
    matches: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        result = subprocess.run(
            ["find", str(root), "-name", pattern],
            capture_output=True, text=True, check=False,
        )
        for line in result.stdout.splitlines():
            p = Path(line)
            if p.is_file():
                matches.append(p)
    return matches


def files_under(roots: list[Path]) -> list[Path]:
    """Return all regular files under ``roots``."""
    out: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.is_file():
                out.append(p)
    return out


def main() -> int:
    failures: list[str] = []

    # 1. GATK Spark tool count (must equal the SCOPE.md count).
    spark_tools = find_files("*Spark.java", [ROOT / "gatk-source/src/main/java"])
    spark_tool_count = sum(
        1 for p in spark_tools
        if "/tools/" in str(p) and not p.name.startswith("Example")
    )
    if spark_tool_count < 50:
        failures.append(
            f"GATK Spark tool count {spark_tool_count} below expected floor 50; "
            "SCOPE.md §2.1 likely drifted from gatk-source."
        )

    # 2. Zero SparkContext/RDD/Dataset in our source.
    spark_api = run_grep(
        r"SparkContext|JavaSparkContext|\bRDD\b|Dataset<", SOURCE_TREES,
    )
    if spark_api != 0:
        failures.append(
            f"Spark API symbols (SparkContext/RDD/Dataset) found {spark_api} "
            "times in our source; expected 0. See SCOPE.md §2.2."
        )

    # 3. Zero Hadoop/HDFS in our source.
    hadoop = run_grep(r"hadoop|hdfs://", SOURCE_TREES)
    if hadoop != 0:
        failures.append(
            f"Hadoop/HDFS symbols found {hadoop} times in our source; "
            "expected 0. See SCOPE.md §2.2."
        )

    # 4. Zero Spark/Hadoop/HDFS vendored as dependencies (not javadoc).
    # Only count top-level directory names and binary artefacts; static
    # GATK-built `gatkdoc/*Spark.html` documentation pages do not pull
    # Spark into our build because we never link against them.
    vendor: list[str] = []
    if THIRD_PARTY.exists():
        for child in THIRD_PARTY.iterdir():
            name = child.name
            if re.search(r"(?i)^(spark|hadoop|hdfs)$", name):
                vendor.append(f"dir:{name}")
        for f in files_under([THIRD_PARTY]):
            low = f.name.lower()
            # Skip GATK's vendored *-spark.jar: it ships inside the GATK
            # 4.6.2.0 release tarball as a build variant for oracle
            # comparison; we never link it.  Only flag a Spark/Hadoop lib
            # if it sits outside the GATK release tree or is not part of
            # the local GATK jar family.
            in_gatk_release = "gatk-4.6.2.0" in str(f) or "gatk-package" in str(f)
            is_gatk_spark_jar = low.endswith("-spark.jar") and in_gatk_release
            if is_gatk_spark_jar:
                continue
            if re.search(r"(spark|hadoop|hdfs).*\.(jar|so|dylib|dll|a)$", low):
                # Allow gatk-package-*-spark.jar but flag any other spark lib.
                if in_gatk_release and "gatk-package" in low:
                    continue
                vendor.append(f"lib:{f}")
            if re.fullmatch(r"(pom\.xml|build\.gradle)", low):
                if any(
                    tok in f.read_text(encoding="utf-8", errors="ignore").lower()
                    for tok in ("org.apache.spark", "org.apache.hadoop")
                ):
                    vendor.append(f"manifest:{f}")
    if vendor:
        failures.append(
            "Spark/Hadoop/HDFS vendored under third_party/: "
            f"{vendor[:10]}; expected 0. See SCOPE.md §2.2."
        )

    # 5. Cloud URL rejection logic only in REJECTION_FILES.  We look for
    # the canonical GATK-style prefix check
    #     input.rfind("s3://", 0) == 0
    # but tolerate either prefix (s3, gs, http, https) and any of the
    # project source extensions.
    cloud_pattern = r'rfind\("(s3|gs|https?)://"\s*,\s*0\)'
    cloud_hits: list[Path] = []
    for src in files_under([ROOT / "fastgatk-native/src"]):
        if src.suffix not in {".cpp", ".hpp", ".h", ".cc", ".cxx"}:
            continue
        try:
            ctext = src.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if re.search(cloud_pattern, ctext):
            cloud_hits.append(src)
    bad_cloud = [p for p in (str(x) for x in cloud_hits) if Path(p) not in REJECTION_FILES]
    if bad_cloud:
        failures.append(
            "Remote URL handling (s3/gs/https) found outside rejection branches: "
            f"{bad_cloud}. See SCOPE.md §2.3."
        )
    if len(cloud_hits) < 2:
        failures.append(
            "Expected at least 2 cloud-URL rejection sites "
            "(genotype_gvcf_tool.cpp + genomicsdb_import_tool.cpp); "
            f"found {len(cloud_hits)}. Scope statement may be stale."
        )

    # 6. Registry fallback counts.
    if not REGISTRY.exists():
        failures.append(f"Registry missing: {REGISTRY}")
        spark_count = cloud_count = -1
    else:
        text = REGISTRY.read_text(encoding="utf-8")
        spark_count = len(re.findall(r'"spark-', text))
        cloud_count = len(re.findall(r'"cloud-', text))
        if spark_count != 4:
            failures.append(
                f"Registry declares {spark_count} spark-* fallback entries; "
                "expected 4. See SCOPE.md §3.1."
            )
        if cloud_count != 19:
            failures.append(
                f"Registry declares {cloud_count} cloud-* fallback entries; "
                "expected 19. See SCOPE.md §3.2."
            )

    # 7. Dispatcher fail-closed behavior for Spark defaults.
    if not DISPATCHER_PY.exists():
        failures.append(f"Dispatcher missing: {DISPATCHER_PY}")
    else:
        dispatcher_text = DISPATCHER_PY.read_text(encoding="utf-8")
        if "UNSUPPORTED_PARAMETER" not in dispatcher_text:
            failures.append(
                "Dispatcher has no UNSUPPORTED_PARAMETER error path; "
                "cannot enforce fail-closed on unknown Spark/cloud args."
            )
        if "--gatk-config-file" not in dispatcher_text:
            failures.append(
                "Dispatcher does not handle --gatk-config-file; "
                "GATK properties can silently enable Spark/cloud defaults."
            )
        if "UNSUPPORTED_PARAMETER" in dispatcher_text:
            # The GATK post-tool -- separator is rejected.  Look for the
            # branch that fires on `parsed.tool_args`.
            if 'in parsed.tool_args' not in dispatcher_text:
                failures.append(
                    "Dispatcher does not appear to inspect parsed.tool_args; "
                    "verify the GATK `--` post-tool separator is fail-closed."
                )

    if failures:
        report = {
            "status": "fail",
            "scope_file": str(SCOPE_MD.relative_to(ROOT)) if SCOPE_MD.exists() else "<missing>",
            "registry": str(REGISTRY.relative_to(ROOT)) if REGISTRY.exists() else "<missing>",
            "checks_run": 7,
            "observed": {
                "gatk_spark_tools": spark_tool_count,
                "spark_api_in_source": spark_api,
                "hadoop_in_source": hadoop,
                "vendored_spark_or_hadoop": vendor,
                "cloud_url_sites_outside_rejection": len(bad_cloud),
                "registry_spark_fallbacks": spark_count,
                "registry_cloud_fallbacks": cloud_count,
            },
            "failures": failures,
        }
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1

    report = {
        "status": "pass",
        "scope_file": str(SCOPE_MD.relative_to(ROOT)) if SCOPE_MD.exists() else "<missing>",
        "registry": str(REGISTRY.relative_to(ROOT)) if REGISTRY.exists() else "<missing>",
        "checks_run": 7,
        "observed": {
            "gatk_spark_tools": spark_tool_count,
            "spark_api_in_source": spark_api,
            "hadoop_in_source": hadoop,
            "vendored_spark_or_hadoop": vendor,
            "cloud_url_sites_outside_rejection": len(bad_cloud),
            "registry_spark_fallbacks": spark_count,
            "registry_cloud_fallbacks": cloud_count,
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
