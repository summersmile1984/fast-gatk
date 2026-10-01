#!/usr/bin/env python3
"""Build the 48-tool GATK 4.6.2.0 vs native CLI coverage inventory.

P0c of the no-Java-fallback 1:1 plan.  This is an inventory, not a claim that
every GATK option is already implemented.  P1 consumes ``must-implement``.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
REGISTRY = ROOT / "fastgatk-native/dispatcher/tool_registry.json"
GATK_JAVA = ROOT / "gatk-source/src/main/java"
CONSTANTS = ROOT / "gatk-source/src/main/java/org/broadinstitute/hellbender/cmdline/StandardArgumentDefinitions.java"
OUT = ROOT / "fastgatk-native/evidence/cli-coverage-48.json"

LAUNCHER_OPTIONS = {
    "--java-options", "-java-options", "--gatk-config-file",
    "--help", "-h", "--version", "--list", "--dry-run", "--dryrun",
    "--fallback", "--java-fallback",
}
DEFERRED_RE = re.compile(
    r"(spark|gcs|google|hadoop|cloud-prefetch|cloud-index-prefetch|"
    r"nio[_-]|requester-pays)",
    re.I,
)
DEFERRED_NAMES = {
    "-CIPB", "-CPB", "--cloud-prefetch-buffer", "--cloud-index-prefetch-buffer",
}
PICARD_VIA_GATK = {"SortSam", "MarkDuplicates", "GatherVcfs"}

CLASS_RE = re.compile(r"\b(?:public|protected|abstract|final|\s)+class\s+(\w+)")
EXTENDS_RE = re.compile(r"\bclass\s+\w+[^{]*\bextends\s+([\w.]+)")
PACKAGE_RE = re.compile(r"^package\s+([\w.]+);", re.M)
IMPORT_RE = re.compile(r"^import\s+([\w.]+);", re.M)
CONST_RE = re.compile(
    r"public\s+static\s+final\s+String\s+(\w+)\s*=\s*\"([^\"]*)\"\s*;"
)
ARG_BLOCK_RE = re.compile(
    r"@(?:Argument|ArgumentCollection)\s*\((.*?)\)",
    re.S,
)
ARG_TYPED_RE = re.compile(
    r"@Argument\s*\((.*?)\)\s*(?:public|protected|private)?\s*"
    r"([\w.<>,\s\?\[\]]+?)\s+(\w+)\s*[;=]",
    re.S,
)
FULL_RE = re.compile(r"fullName\s*=\s*([^,\s]+)")
SHORT_RE = re.compile(r"shortName\s*=\s*([^,\s]+)")
HIDDEN_RE = re.compile(r"\bhidden\s*=\s*true\b")
COLLECTION_TYPE_RE = re.compile(
    r"@(?:ArgumentCollection)[^{;]*?\s+([\w.]+)\s+\w+\s*[;=]",
    re.S,
)


def load_constants(paths: Iterable[Path]) -> dict[str, str]:
    table: dict[str, str] = {}
    for path in paths:
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for name, value in CONST_RE.findall(text):
            table[name] = value
    changed = True
    while changed:
        changed = False
        for name, value in list(table.items()):
            key = value.split(".")[-1]
            if key in table and table[key] != value:
                table[name] = table[key]
                changed = True
    return table


def index_java(root: Path) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for path in root.rglob("*.java"):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        pkg_match = PACKAGE_RE.search(text)
        pkg = pkg_match.group(1) if pkg_match else ""
        for name in CLASS_RE.findall(text):
            index[name] = path
            if pkg:
                index[f"{pkg}.{name}"] = path
    return index


def _as_cli(value: str, *, short: bool) -> str:
    if value.startswith("-"):
        return value
    if short or len(value) <= 3:
        return f"-{value}"
    return f"--{value}"


def resolve_name(raw: str, constants: dict[str, str], *, short: bool) -> str | None:
    raw = raw.strip().rstrip(",")
    if not raw:
        return None
    if raw.startswith('"') and raw.endswith('"'):
        return _as_cli(raw[1:-1], short=short)
    key = raw.split(".")[-1]
    if key in constants:
        return _as_cli(constants[key], short=short or key.endswith("SHORT_NAME"))
    return None


def _names_from_block(block: str, constants: dict[str, str]) -> list[str]:
    names: list[str] = []
    full = FULL_RE.search(block)
    short = SHORT_RE.search(block)
    if full:
        resolved = resolve_name(full.group(1), constants, short=False)
        if resolved:
            names.append(resolved)
    if short:
        resolved = resolve_name(short.group(1), constants, short=True)
        if resolved:
            names.append(resolved)
    return names


def parse_arguments(text: str, constants: dict[str, str]) -> tuple[set[str], set[str], set[str]]:
    options: set[str] = set()
    hidden: set[str] = set()
    flags: set[str] = set()
    for block, type_name, _field in ARG_TYPED_RE.findall(text):
        names = _names_from_block(block, constants)
        target = hidden if HIDDEN_RE.search(block) else options
        target.update(names)
        simple = type_name.strip().split("<", 1)[0].split(".", 1)[-1].strip()
        if simple in {"boolean", "Boolean"}:
            flags.update(names)
    # ArgumentCollection blocks and any @Argument the typed regex missed.
    for block in ARG_BLOCK_RE.findall(text):
        names = _names_from_block(block, constants)
        target = hidden if HIDDEN_RE.search(block) else options
        target.update(names)
    return options, hidden, flags


def superclasses(text: str, path: Path, index: dict[str, Path]) -> list[Path]:
    found: list[Path] = []
    match = EXTENDS_RE.search(text)
    if not match:
        return found
    name = match.group(1).split(".")[-1]
    parent = index.get(name)
    if parent and parent != path:
        found.append(parent)
    return found


def collection_types(text: str, index: dict[str, Path]) -> list[Path]:
    found: list[Path] = []
    for type_name in COLLECTION_TYPE_RE.findall(text):
        simple = type_name.split(".")[-1]
        path = index.get(simple) or index.get(type_name)
        if path:
            found.append(path)
    return found


def collect_tool_options(class_name: str, index: dict[str, Path],
                         constants: dict[str, str]) -> tuple[set[str], set[str], set[str], list[str]]:
    start = index.get(class_name)
    if start is None:
        return set(), set(), set(), []
    seen: set[Path] = set()
    queue = [start]
    options: set[str] = set()
    hidden: set[str] = set()
    flags: set[str] = set()
    files: list[str] = []
    while queue:
        path = queue.pop(0)
        if path in seen:
            continue
        seen.add(path)
        files.append(str(path.relative_to(ROOT)))
        text = path.read_text(encoding="utf-8", errors="replace")
        more, hid, more_flags = parse_arguments(text, constants)
        options.update(more)
        hidden.update(hid)
        flags.update(more_flags)
        queue.extend(superclasses(text, path, index))
        queue.extend(collection_types(text, index))
    return options, hidden, flags, files


def registry_options(entry: dict) -> set[str]:
    names: set[str] = set()
    for key in ("value_options", "flag_options", "repeatable_options",
                "optional_boolean_options"):
        for item in entry.get(key, []):
            names.add(str(item))
    return names


def classify(name: str) -> str:
    if name in LAUNCHER_OPTIONS:
        return "launcher"
    if name in DEFERRED_NAMES or DEFERRED_RE.search(name):
        return "deferred-cloud-or-spark"
    return "tool"


def main() -> int:
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    tools = registry["tools"]
    index = index_java(GATK_JAVA) if GATK_JAVA.is_dir() else {}
    constant_files = [CONSTANTS] if CONSTANTS.is_file() else []
    if GATK_JAVA.is_dir():
        constant_files.extend(GATK_JAVA.rglob("*.java"))
    constants = load_constants(constant_files)

    payload = {
        "schema_version": 1,
        "gatk_compatibility": registry.get("gatk_compatibility"),
        "runtime_java": registry.get("runtime_java", False),
        "tool_count": len(tools),
        "tools": {},
        "summary": {},
    }
    totals = defaultdict(int)
    for tool, entry in tools.items():
        native = registry_options(entry)
        gatk_opts, hidden, flags, sources = collect_tool_options(tool, index, constants)
        if not sources and tool in PICARD_VIA_GATK:
            bucket = {
                "status": entry.get("status"),
                "native_options": sorted(native),
                "gatk_options": [],
                "gatk_flags": [],
                "hidden_gatk_options": [],
                "implemented": sorted(native),
                "must-implement": [],
                "gatk-also-rejects": [],
                "belongs-to-other-tool": [],
                "deferred-cloud-or-spark": [],
                "launcher": [],
                "note": "Picard command wrapped by GATK; Java class is not in gatk-source.",
                "sources": [],
            }
            payload["tools"][tool] = bucket
            totals["picard-via-gatk"] += 1
            continue

        implemented = sorted(gatk_opts & native)
        remainder = gatk_opts - native
        must: list[str] = []
        deferred: list[str] = []
        launcher: list[str] = []
        other: list[str] = []
        for opt in sorted(remainder):
            kind = classify(opt)
            if kind == "launcher":
                launcher.append(opt)
            elif kind == "deferred-cloud-or-spark":
                deferred.append(opt)
            else:
                must.append(opt)
        extra = sorted(native - gatk_opts - LAUNCHER_OPTIONS)
        bucket = {
            "status": entry.get("status"),
            "native_options": sorted(native),
            "gatk_options": sorted(gatk_opts),
            "gatk_flags": sorted(flags & gatk_opts),
            "hidden_gatk_options": sorted(hidden),
            "implemented": implemented,
            "must-implement": must,
            "gatk-also-rejects": extra,
            "belongs-to-other-tool": other,
            "deferred-cloud-or-spark": deferred,
            "launcher": launcher,
            "sources": sources,
        }
        payload["tools"][tool] = bucket
        totals["must-implement"] += len(must)
        totals["implemented"] += len(implemented)
        totals["deferred-cloud-or-spark"] += len(deferred)
        totals["tools-with-gaps"] += 1 if must else 0

    payload["summary"] = dict(totals)
    payload["summary"]["tools"] = len(tools)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "pass",
        "path": str(OUT),
        "tools": len(tools),
        "must-implement": totals["must-implement"],
        "tools-with-gaps": totals["tools-with-gaps"],
        "deferred-cloud-or-spark": totals["deferred-cloud-or-spark"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
